// main.cpp  ---------------------------------------------------------------
// 主控程序：循环处理扫描式 GMTI 数据
// v 2025-08-07
// ------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include "config_structs.hpp" // 含 Config / GMTIOutput::Plane / Result 等声明
#include "GMTIProcessor.hpp"  // 你的处理类头文件
#include "dbs/DbsFusion.hpp"
#include "trackModule.hpp"
#include "TrackManager.hpp"
#include "pesudoTargetGen.hpp" // 伪目标生成函数声明

namespace {

bool deriveRuntimeConfig(Config &cfg)
{
    if (cfg.INFO_Type) {
        // 新协议：每个 PRT 256 字节头 + 4096 个采样点，每点双通道 complex float32。
        cfg.pkg_bytes = 256 + 4096 * 16;
    } else if (cfg.channel_mode == "separate") {
        // 旧协议：双文件单通道，每脉冲为包头 + (I,Q)*float32。
        cfg.pkg_bytes = cfg.info_len + cfg.pulse_len * 2 * 4;
    } else if (cfg.channel_mode == "interleaved") {
        // 旧协议：单文件双通道交织，每脉冲为包头 + (I1,Q1,I2,Q2)*float32。
        cfg.pkg_bytes = cfg.info_len + cfg.pulse_len * 4 * 4;
    } else {
        std::cerr << "[ERR] 未知 channel_mode: " << cfg.channel_mode << std::endl;
        return false;
    }

    cfg.lambda = C / cfg.fc;
    cfg.R_bin = C / (2.0 * cfg.fs);
    cfg.Rg.resize(cfg.rg_len);
    for (int i = 0; i < cfg.rg_len; ++i) {
        cfg.Rg[i] = cfg.R_min + i * cfg.R_bin;
    }
    const int procPulseNum = effectivePulseNum(cfg);
    cfg.fd_res = cfg.PRF / static_cast<double>(procPulseNum);

    // 初始占位；单 period 内会根据多普勒中心和动态支撑域重新计算。
    cfg.az_st = 1;
    cfg.az_ed = procPulseNum;
    cfg.az_center = (cfg.az_st + cfg.az_ed) / 2;
    cfg.Loc = true;

    return true;
}

} // namespace

int main(int argc, char **argv)
{
    TIMING_SCOPE(main_total);
    std::vector<double> MT_acc;
    std::vector<std::vector<double>> all_frames_MT;

#ifdef _OPENMP
    omp_set_num_threads(std::max(1, omp_get_max_threads()-1)); // 留一个给系统
    DBG("OpenMP is enabled. Using " << omp_get_max_threads() - 1 << " threads.");
#endif

    // ===== 1. 读取 XML 配置 ==============================================
    std::string xmlPath = "/home/shy/AIR/小长/GMTI程序/GMTI/GMTI_algorithm/temp_config.xml";
    if (argc >= 2)
        xmlPath = argv[1];

    GMTIProcessor proc;
    Config cfg; // 用于回传 cfg，方便 main 修改/查看
    if (!proc.readXmlParam(xmlPath, cfg)) { 
        // 内部也已写入私有 cfg_
        std::cerr << "[ERR] 读取 XML 失败: " << xmlPath << std::endl;
        return 1;
    }

    if (!deriveRuntimeConfig(cfg)) {
        return 1;
    }

    // ===== 2. 预读位姿信息 =================================================
    std::vector<std::vector<double>> posMatrix;
    int POS_num = 0;
    if (!cfg.INFO_Type) {
        // 老协议：从外部 POS 文件读取
        if (!proc.POS_dataread(cfg.Plane_POS_add, posMatrix, POS_num)) {
            std::cerr << "[ERR] 读取 POS 失败: " << cfg.Plane_POS_add << std::endl;
            return 1;
        }
    }
    
    std::vector<int> periodList;
    GMTIOutput::Plane plane{};

    if (cfg.wavepos_use_roi) {
        // ROI 选波位需要初始飞机位姿。新协议没有外部 POS 文件时，
        // extractPlanePV() 会自动从回波包头提取；默认范围选波位则不预读。
        if(!proc.extractPlanePV(posMatrix, cfg, plane)) {
            std::cerr << "[ERR] 飞机位置提取失败\n";
            return 1;
        }
    }

    if(!proc.makePeriodList(plane, cfg, periodList)) {
        std::cerr << "[ERR] 生成 periodList 失败\n";
        return 1;
    }
    MT_acc.reserve(periodList.size() * 40 * 8);

    // ===== 阶段 1: 检测与定位 ===============================================
    // 并行处理所有 period（在单 GPU 上通过多个 GMTIProcessor 实例并发运行）
    // FFT/cuFFT plan 在 processPeriodsParallel 内为每个 worker 单独初始化
    std::vector<GMTIOutput> periodResults;
    FusionGroupContext fusionCtx;
    const std::vector<std::vector<double>> posSource = cfg.INFO_Type ? std::vector<std::vector<double>>() : posMatrix;
    const bool processOk = cfg.enable_dbs_fusion
        ? proc.processPeriodsParallelFusion(periodList, cfg, posSource, fusionCtx, periodResults)
        : proc.processPeriodsParallel(periodList, cfg, posSource, periodResults);
    if (!processOk) {
        std::cerr << "[ERR] 并行处理部分 period 失败（见日志）。\n";
        // 继续尝试从成功的结果中收集 MT
    }
    if (cfg.enable_dbs_fusion && processOk) {
        if (!runDbsFusionImaging(fusionCtx, cfg, true)) {
            std::cerr << "[ERR] DBS 融合成像输出失败。\n";
        }
    }
    
    // 按输入 periodList 的顺序收集 MT 以保证确定性
    for (size_t i = 0; i < periodResults.size(); ++i) {
        auto &res = periodResults[i];
        if (!res.MT.empty()) {
            all_frames_MT.push_back(res.MT);
            MT_acc.insert(MT_acc.end(), std::make_move_iterator(res.MT.begin()),
                          std::make_move_iterator(res.MT.end()));
        }
        // std::cout << "[OK] Period " << periodList[i] << " 完成\n";
    }

    // 整个扫描周期只写一次检测结果：2 + N*36
    bool wroteCurrentDetections = false;
    if (!MT_acc.empty()) {
        if (!proc.writeResult(MT_acc, cfg)) {
            std::cerr << "[ERR] 整周期检测结果写盘失败\n";
        } else {
            wroteCurrentDetections = true;
        }
    } else {
        std::cerr << "[WARN] 本周期无检测目标，跳过检测结果写盘\n";
    }
    
    std::cout << "\n=== 检测定位阶段完成，开始航迹关联 ===\n";

    // ===== 阶段 2: 航迹关联 ===============================================
    // 调用 trackModule 进行航迹关联（从磁盘读取检测结果）
    // trackModule 参数：result_dir, idx_range, assoc_window
    if (!wroteCurrentDetections) {
        std::cout << "[TRACK][WARN] 本周期未写入新的检测结果文件，关联窗口中的当前编号文件可能是旧数据。" << std::endl;
    }
    TrackManager track_manager;
    trackModule(cfg, &track_manager);
    
    // ===== 后续处理（可选）===============================================
    // 注：trackModule 已经将航迹关联结果写到磁盘
    // 如需进一步处理轨迹，可在这里添加代码

    std::cout << "\n--- GMTI 扫描完成 ---\n";
    return 0;
}
