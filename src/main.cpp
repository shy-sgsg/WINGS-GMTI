// main.cpp  ---------------------------------------------------------------
// 主控程序：循环处理扫描式 GMTI 数据
// v 2025-08-07
// ------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <limits>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <cerrno>
#include <sys/stat.h>
#include "config_structs.hpp" // 含 Config / GMTIOutput::Plane / Result 等声明
#include "GMTIProcessor.hpp"  // 你的处理类头文件
#include "dbs/DbsFusion.hpp"
#include "trackModule.hpp"
#include "TrackManager.hpp"
#include "pesudoTargetGen.hpp" // 伪目标生成函数声明
#include "trig_lut.hpp"
#include "runtime_diagnostics.hpp"

namespace {

std::string joinPathMain(const std::string &a, const std::string &b)
{
    if (a.empty()) return b;
    if (a.back() == '/' || a.back() == '\\') return a + b;
    return a + "/" + b;
}

bool ensureDirMain(const std::string &dir)
{
    if (dir.empty()) return true;
    std::string cur;
    size_t i = 0;
    if (dir[0] == '/') {
        cur = "/";
        i = 1;
    }
    for (; i < dir.size(); ++i) {
        cur.push_back(dir[i]);
        if (dir[i] == '/' || i == dir.size() - 1) {
            if (cur == "/" || cur == "./" || cur == ".") continue;
            if (::mkdir(cur.c_str(), 0755) != 0 && errno != EEXIST) {
                return false;
            }
        }
    }
    return true;
}

std::vector<std::string> splitCsvMain(const std::string &line)
{
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
        if (ch == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (ch == ',' && !quoted) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    out.push_back(cur);
    return out;
}

double toDoubleMain(const std::string &s, double fallback = std::numeric_limits<double>::quiet_NaN())
{
    if (s.empty()) return fallback;
    char *end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return end && *end == '\0' ? v : fallback;
}

int toIntMain(const std::string &s, int fallback = -1)
{
    const double v = toDoubleMain(s);
    return std::isfinite(v) ? static_cast<int>(std::lround(v)) : fallback;
}

struct PcPeakRow {
    std::string case_id;
    std::string run_id;
    int period_id = -1;
    int beam_id = -1;
    double truth_range_m = std::numeric_limits<double>::quiet_NaN();
    double truth_range_sample_float = std::numeric_limits<double>::quiet_NaN();
    int truth_range_sample_int = -1;
    int pc_crop_start = 0;
    int expected_pc_bin_without_offset = -1;
    int pc_peak_bin = -1;
    double pc_peak_range_m = std::numeric_limits<double>::quiet_NaN();
    double pc_peak_power = std::numeric_limits<double>::quiet_NaN();
    int pc_peak_offset_bin = 0;
    double pc_peak_range_error_m = std::numeric_limits<double>::quiet_NaN();
    int global_peak_bin = -1;
    double global_peak_range_m = std::numeric_limits<double>::quiet_NaN();
    double global_peak_power = std::numeric_limits<double>::quiet_NaN();
    int valid = 0;
};

std::vector<PcPeakRow> readPcPeakRows(const std::string &path)
{
    std::ifstream in(path.c_str());
    std::vector<PcPeakRow> rows;
    if (!in) return rows;
    std::string line;
    if (!std::getline(in, line)) return rows;
    const auto header = splitCsvMain(line);
    std::map<std::string, size_t> col;
    for (size_t i = 0; i < header.size(); ++i) col[header[i]] = i;
    auto get = [&](const std::vector<std::string> &cells, const std::string &name) -> std::string {
        const auto it = col.find(name);
        return (it == col.end() || it->second >= cells.size()) ? std::string() : cells[it->second];
    };
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto cells = splitCsvMain(line);
        PcPeakRow r;
        r.case_id = get(cells, "case_id");
        r.run_id = get(cells, "run_id");
        r.period_id = toIntMain(get(cells, "period_id"), -1);
        r.beam_id = toIntMain(get(cells, "beam_id"), -1);
        r.truth_range_m = toDoubleMain(get(cells, "truth_range_m"));
        r.truth_range_sample_float = toDoubleMain(get(cells, "truth_range_sample_float"));
        r.truth_range_sample_int = toIntMain(get(cells, "truth_range_sample_int"), -1);
        r.pc_crop_start = toIntMain(get(cells, "pc_crop_start"), 0);
        r.expected_pc_bin_without_offset = toIntMain(get(cells, "expected_pc_bin_without_offset"), -1);
        r.pc_peak_bin = toIntMain(get(cells, "pc_peak_bin"), -1);
        r.pc_peak_range_m = toDoubleMain(get(cells, "pc_peak_range_m"));
        r.pc_peak_power = toDoubleMain(get(cells, "pc_peak_power"));
        r.pc_peak_offset_bin = toIntMain(get(cells, "pc_peak_offset_bin"), 0);
        r.pc_peak_range_error_m = toDoubleMain(get(cells, "pc_peak_range_error_m"));
        r.global_peak_bin = toIntMain(get(cells, "global_peak_bin"), -1);
        r.global_peak_range_m = toDoubleMain(get(cells, "global_peak_range_m"));
        r.global_peak_power = toDoubleMain(get(cells, "global_peak_power"));
        r.valid = toIntMain(get(cells, "valid"), 0);
        rows.push_back(r);
    }
    return rows;
}

void writeRangeAxisAudit(const Config &cfg)
{
    if (!cfg.debug_pc_peak) return;
    const std::string debug_dir = joinPathMain(cfg.result_add, "debug");
    if (!ensureDirMain(debug_dir)) return;
    const std::string path = joinPathMain(debug_dir, "range_axis_audit.csv");
    const double fs_hz = cfg.fs;
    const double sample_delay_us = cfg.has_sample_delay_us
        ? cfg.sample_delay_us
        : std::numeric_limits<double>::quiet_NaN();
    auto expectedRange = [&](int bin) -> double {
        if (!cfg.has_sample_delay_us || fs_hz <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const double raw_sample = static_cast<double>(cfg.range_crop_start + bin);
        return 0.5 * C * (sample_delay_us * 1.0e-6 + raw_sample / fs_hz);
    };
    const double front = cfg.Rg.empty() ? std::numeric_limits<double>::quiet_NaN() : cfg.Rg.front();
    const double back = cfg.Rg.empty() ? std::numeric_limits<double>::quiet_NaN() : cfg.Rg.back();
    const double expected_first = expectedRange(0);
    const double expected_last = expectedRange(std::max(0, cfg.rg_len - 1));
    const bool range_window_ok = std::isfinite(front) && std::isfinite(back) &&
                                 front >= 82000.0 && front <= 83000.0 &&
                                 back >= 92000.0 && back <= 93050.0;
    std::ofstream os(path.c_str());
    os << "fs_hz,sample_delay_us,pulse_len,rg_len,pc_crop_start,pc_crop_len,"
          "R_min,R_bin,cfg_Rg_size,cfg_Rg_front,cfg_Rg_back,"
          "expected_Rg_first,expected_Rg_last,Rg_first_error_m,Rg_last_error_m,"
          "range_axis_status\n";
    os << std::setprecision(15)
       << fs_hz << ","
       << sample_delay_us << ","
       << cfg.pulse_len << ","
       << cfg.rg_len << ","
       << cfg.range_crop_start << ","
       << cfg.range_compress_len << ","
       << cfg.R_min << ","
       << cfg.R_bin << ","
       << cfg.Rg.size() << ","
       << front << ","
       << back << ","
       << expected_first << ","
       << expected_last << ","
       << (front - expected_first) << ","
       << (back - expected_last) << ","
       << (range_window_ok ? "new_system_range_axis_consistent" : "range_axis_not_aligned_with_new_system_window")
       << "\n";
}

void writeFullPipelineCalibrationCheck(const Config &cfg,
                                       const std::vector<GMTIOutput::DetectionCsvRecord> &detections)
{
    if (!cfg.debug_pc_peak) return;
    const std::string debug_dir = joinPathMain(cfg.result_add, "debug");
    if (!ensureDirMain(debug_dir)) return;
    const std::vector<PcPeakRow> peaks = readPcPeakRows(joinPathMain(debug_dir, "pc_peak_check.csv"));
    if (peaks.empty()) return;

    const std::string csv_path = joinPathMain(debug_dir, "full_pipeline_calibration_check.csv");
    const std::string report_path = joinPathMain(debug_dir, "full_pipeline_calibration_report.md");
    std::ofstream os(csv_path.c_str());
    os << "case_id,run_id,period_id,beam_id,"
          "truth_range_m,truth_range_sample_float,truth_range_sample_int,"
          "pc_crop_start,expected_pc_bin_without_offset,"
          "pc_peak_bin,pc_peak_range_m,pc_peak_power,"
          "det_found,det_id,det_range_bin,det_range_m,det_row,det_col,det_power,det_amplitude,"
          "pc_to_det_bin_error,truth_to_det_bin_error,pc_to_det_range_error_m,"
          "match_status,failure_reason\n";

    int consistent = 0;
    int no_detection = 0;
    int range_axis_wrong = 0;
    int pc_invalid = 0;
    const int match_window = 20;
    for (const auto &p : peaks) {
        int best = -1;
        int best_bin_delta = 1000000000;
        double best_power = -1.0;
        if (p.valid && p.pc_peak_bin >= 0) {
            for (size_t i = 0; i < detections.size(); ++i) {
                const auto &d = detections[i];
                if (d.range_bin < 0) continue;
                const bool period_ok = (p.period_id <= 0 || d.period_id == p.period_id || d.period_id < 0);
                const bool beam_ok = (p.beam_id < 0 || d.beam_id == p.beam_id ||
                                      (d.beam_id >= 0 && std::abs(d.beam_id - p.beam_id) <= 1));
                if (!period_ok || !beam_ok) continue;
                const int delta = std::abs(d.range_bin - p.pc_peak_bin);
                if (delta > match_window) continue;
                const double power = std::isfinite(d.amplitude) ? d.amplitude * d.amplitude : 0.0;
                if (delta < best_bin_delta || (delta == best_bin_delta && power > best_power)) {
                    best = static_cast<int>(i);
                    best_bin_delta = delta;
                    best_power = power;
                }
            }
        }

        bool det_found = best >= 0;
        int det_id = best;
        int det_range_bin = -1;
        int det_row = -1;
        int det_col = -1;
        double det_range_m = std::numeric_limits<double>::quiet_NaN();
        double det_amp = std::numeric_limits<double>::quiet_NaN();
        double det_power = std::numeric_limits<double>::quiet_NaN();
        int pc_to_det_bin_error = 0;
        int truth_to_det_bin_error = 0;
        double pc_to_det_range_error_m = std::numeric_limits<double>::quiet_NaN();
        std::string match_status;
        std::string failure_reason;
        if (!p.valid || p.pc_peak_bin < 0) {
            match_status = "pc_peak_invalid";
            failure_reason = "no_valid_pc_peak";
            ++pc_invalid;
        } else if (!det_found) {
            match_status = "pc_peak_valid_but_no_detection";
            failure_reason = "cancelled_or_cfar_filtered_or_gated";
            ++no_detection;
        } else {
            const auto &d = detections[static_cast<size_t>(best)];
            det_range_bin = d.range_bin;
            det_range_m = d.range_m;
            det_row = d.row;
            det_col = d.col;
            det_amp = d.amplitude;
            det_power = std::isfinite(det_amp) ? det_amp * det_amp : std::numeric_limits<double>::quiet_NaN();
            pc_to_det_bin_error = det_range_bin - p.pc_peak_bin;
            truth_to_det_bin_error = det_range_bin - p.expected_pc_bin_without_offset;
            pc_to_det_range_error_m = det_range_m - p.pc_peak_range_m;
            if (std::abs(pc_to_det_bin_error) <= match_window &&
                std::isfinite(det_range_m) && std::isfinite(p.truth_range_m) &&
                std::abs(det_range_m - p.truth_range_m) > 1000.0) {
                match_status = "bin_consistent_range_axis_wrong";
                failure_reason = "detection_bin_matches_pc_peak_but_range_m_disagrees_with_truth";
                ++range_axis_wrong;
            } else if (std::abs(pc_to_det_bin_error) <= match_window) {
                match_status = "full_pipeline_consistent";
                failure_reason = "";
                ++consistent;
            } else {
                match_status = "pipeline_consistent_but_pc_mapping_wrong";
                failure_reason = "detected_bin_far_from_expected_truth_bin";
            }
        }

        os << p.case_id << ","
           << p.run_id << ","
           << p.period_id << ","
           << p.beam_id << ","
           << std::setprecision(15)
           << p.truth_range_m << ","
           << p.truth_range_sample_float << ","
           << p.truth_range_sample_int << ","
           << p.pc_crop_start << ","
           << p.expected_pc_bin_without_offset << ","
           << p.pc_peak_bin << ","
           << p.pc_peak_range_m << ","
           << p.pc_peak_power << ","
           << (det_found ? "true" : "false") << ","
           << (det_found ? det_id : -1) << ","
           << det_range_bin << ","
           << det_range_m << ","
           << det_row << ","
           << det_col << ","
           << det_power << ","
           << det_amp << ","
           << pc_to_det_bin_error << ","
           << truth_to_det_bin_error << ","
           << pc_to_det_range_error_m << ","
           << match_status << ","
           << failure_reason << "\n";
    }

    std::ofstream report(report_path.c_str());
    report << "# 全流程检测级标定与脉压峰一致性报告\n\n";
    report << "- pc_peak_check 行数：" << peaks.size() << "\n";
    report << "- detection 记录数：" << detections.size() << "\n";
    report << "- 匹配窗口：±" << match_window << " bin\n";
    report << "- full_pipeline_consistent：" << consistent << "\n";
    report << "- pc_peak_valid_but_no_detection：" << no_detection << "\n";
    report << "- bin_consistent_range_axis_wrong：" << range_axis_wrong << "\n";
    report << "- pc_peak_invalid：" << pc_invalid << "\n\n";
    const PcPeakRow *strongest = nullptr;
    for (const auto &p : peaks) {
        if (!strongest || p.global_peak_power > strongest->global_peak_power) {
            strongest = &p;
        }
    }
    if (strongest) {
        const auto &p = *strongest;
        report << "## 最强脉压峰行\n\n";
        report << "- period_id / beam_id：" << p.period_id << " / " << p.beam_id << "\n";
        report << "- expected_pc_bin_without_offset：" << p.expected_pc_bin_without_offset << "\n";
        report << "- pc_peak_bin：" << p.pc_peak_bin << "\n";
        report << "- pc_peak_power：" << p.pc_peak_power << "\n";
        report << "- pc_peak_offset_bin：" << p.pc_peak_offset_bin << "\n";
        report << "- pc_peak_range_error_m：" << p.pc_peak_range_error_m << "\n";
        report << "- global_peak_bin：" << p.global_peak_bin << "\n";
        report << "- global_peak_power：" << p.global_peak_power << "\n";
        report << "- global_peak_range_error_m：" << (p.global_peak_range_m - p.truth_range_m) << "\n";
        report << "- pc_peak_bin 是否在 expected 附近："
               << ((p.valid && std::abs(p.pc_peak_bin - p.expected_pc_bin_without_offset) <= 50) ? "是" : "否") << "\n";
        const bool global_in_expected_window =
            p.global_peak_bin >= 0 &&
            std::abs(p.global_peak_bin - p.expected_pc_bin_without_offset) <= 50;
        report << "- global_peak_bin 是否在 expected±50 内："
               << (global_in_expected_window ? "是" : "否") << "\n";
        if (!global_in_expected_window) {
            report << "- 判定：最强脉压峰不在 truth 预期窗口内；当前不能把 `expected_pc_bin_without_offset` "
                      "视为已标定通过。若 `pc_peak_bin` 落在搜索窗口边界，优先按全局峰偏移分析 "
                      "LFM 符号、匹配滤波对齐、`pc_crop_start` 或 `sample_delay`。\n";
        }
    }
    report << "\n## 结论说明\n\n";
    report << "- 本检查优先按 `range_bin` 匹配，不使用 `detection_results.csv` 的 `range_m` 作为第一匹配依据。\n";
    report << "- 若脉压峰存在但没有检测点落入窗口，说明目标可能在对消、CFAR、聚类或门控后未输出，不等同于脉压峰标定失败。\n";
    report << "- 若 bin 一致但 `range_m` 与 truth 明显不一致，应优先检查 `cfg.Rg` 距离轴。\n";
}

bool parseTrigModeArg(const std::string& arg, const char* next, bool& consumedNext)
{
    consumedNext = false;
    const std::string prefix = "--trig-mode=";
    if (arg.compare(0, prefix.size(), prefix) == 0) {
        return gmti::trig_lut::setModeFromString(arg.c_str() + prefix.size(), false);
    }
    if (arg == "--trig-mode") {
        if (!next) {
            std::cerr << "[ERR] --trig-mode 需要参数: lut|math|compare\n";
            return false;
        }
        consumedNext = true;
        return gmti::trig_lut::setModeFromString(next, false);
    }
    return true;
}

bool deriveRuntimeConfig(Config &cfg)
{
    if (cfg.INFO_Type) {
        // 新协议：每个 PRT 256 字节头，距离点数由 pulse_len 配置。
        cfg.pkg_bytes = cfg.info_len + cfg.pulse_len * 16;
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
    const bool use_new_range_axis = cfg.has_sample_delay_us && usesRangeCropWindow(cfg);
    for (int i = 0; i < cfg.rg_len; ++i) {
        if (use_new_range_axis) {
            const double raw_sample = static_cast<double>(cfg.range_crop_start + i);
            cfg.Rg[i] = 0.5 * C * (cfg.sample_delay_us * 1.0e-6 + raw_sample / cfg.fs);
        } else {
            cfg.Rg[i] = cfg.R_min + i * cfg.R_bin;
        }
    }
    if (use_new_range_axis && !cfg.Rg.empty()) {
        cfg.R_min = cfg.Rg.front();
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
    gmti::trig_lut::configureFromEnv();

    std::string xmlPath = "/home/shy/AIR/小长/GMTI程序/GMTI/GMTI_algorithm/temp_config.xml";
    std::string runtimeModeOverride;
    bool diagnosticsOverrideSet = false;
    bool diagnosticsOverride = true;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        bool consumedNext = false;
        if (arg == "--trig-mode" || arg.compare(0, std::string("--trig-mode=").size(), "--trig-mode=") == 0) {
            if (!parseTrigModeArg(arg, (i + 1 < argc) ? argv[i + 1] : nullptr, consumedNext)) {
                return 1;
            }
            if (consumedNext) ++i;
        } else if (arg.compare(0, std::string("--runtime-mode=").size(), "--runtime-mode=") == 0) {
            runtimeModeOverride = arg.substr(std::string("--runtime-mode=").size());
        } else if (arg == "--runtime-mode") {
            if (i + 1 >= argc) {
                std::cerr << "[ERR] --runtime-mode 需要参数: debug|release\n";
                return 1;
            }
            runtimeModeOverride = argv[++i];
        } else if (arg.compare(0, std::string("--runtime-diagnostics=").size(), "--runtime-diagnostics=") == 0) {
            const std::string value = arg.substr(std::string("--runtime-diagnostics=").size());
            diagnosticsOverrideSet = true;
            diagnosticsOverride = (value == "1" || value == "true" || value == "on");
        } else {
            xmlPath = arg;
        }
    }

    if (!gmti::trig_lut::initialize(true)) {
        std::cerr << "[ERR] 三角函数路径初始化失败\n";
        return 1;
    }
    // gmti::trig_lut::benchmark(1u << 22);

    std::vector<double> MT_acc;
    std::vector<std::vector<double>> all_frames_MT;

#ifdef _OPENMP
    omp_set_num_threads(std::max(1, omp_get_max_threads()-1)); // 留一个给系统
    DBG("OpenMP is enabled. Using " << omp_get_max_threads() - 1 << " threads.");
#endif

    // ===== 1. 读取 XML 配置 ==============================================
    GMTIProcessor proc;
    Config cfg; // 用于回传 cfg，方便 main 修改/查看
    const auto config_load_wall_start = std::chrono::system_clock::now();
    const auto config_load_steady_start = std::chrono::high_resolution_clock::now();
    if (!proc.readXmlParam(xmlPath, cfg)) { 
        // 内部也已写入私有 cfg_
        std::cerr << "[ERR] 读取 XML 失败: " << xmlPath << std::endl;
        return 1;
    }

    if (!deriveRuntimeConfig(cfg)) {
        return 1;
    }
    if (!runtimeModeOverride.empty()) {
        cfg.runtime_mode = runtimeModeOverride;
        std::transform(cfg.runtime_mode.begin(), cfg.runtime_mode.end(), cfg.runtime_mode.begin(),
                       [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        if (cfg.runtime_mode == "release" || cfg.runtime_mode == "formal" ||
            cfg.runtime_mode == "production") {
            cfg.runtime_diagnostics_enabled = false;
        } else if (cfg.runtime_mode == "debug" || cfg.runtime_mode == "trace") {
            cfg.runtime_diagnostics_enabled = true;
        }
    }
    if (diagnosticsOverrideSet) {
        cfg.runtime_diagnostics_enabled = diagnosticsOverride;
    }
    if (!cfg.runtime_diagnostics_enabled) {
        cfg.track_debug_level = 0;
        cfg.track_debug_dump = false;
        cfg.track_debug_dump_level = 0;
    }
    const auto config_load_wall_end = std::chrono::system_clock::now();
    const auto config_load_steady_end = std::chrono::high_resolution_clock::now();
    gmti::runtime::initializeRun(cfg, xmlPath, argv[0] ? argv[0] : "GMTI_core");
    writeRangeAxisAudit(cfg);
    gmti::runtime::recordTiming(
        "config_load",
        config_load_wall_start,
        config_load_wall_end,
        std::chrono::duration_cast<std::chrono::milliseconds>(
            config_load_steady_end - config_load_steady_start).count(),
        -1,
        "readXmlParam+deriveRuntimeConfig");

    // ===== 2. 预读位姿信息 =================================================
    std::vector<std::vector<double>> posMatrix;
    int POS_num = 0;
    if (!cfg.INFO_Type) {
        // 老协议：从外部 POS 文件读取
        if (!proc.POS_dataread(cfg.Plane_POS_add, posMatrix, POS_num)) {
            std::cerr << "[ERR] 读取 POS 失败: " << cfg.Plane_POS_add << std::endl;
            gmti::runtime::finishRun(cfg, false, 1, "POS_dataread failed");
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
            gmti::runtime::finishRun(cfg, false, 1, "extractPlanePV failed");
            return 1;
        }
    }

    if(!proc.makePeriodList(plane, cfg, periodList)) {
        std::cerr << "[ERR] 生成 periodList 失败\n";
        gmti::runtime::finishRun(cfg, false, 1, "makePeriodList failed");
        return 1;
    }
    MT_acc.reserve(periodList.size() * 40 * 8);

    // ===== 阶段 1: 检测与定位 ===============================================
    // 并行处理所有 period（在单 GPU 上通过多个 GMTIProcessor 实例并发运行）
    // FFT/cuFFT plan 在 processPeriodsParallel 内为每个 worker 单独初始化
    std::vector<GMTIOutput> periodResults;
    FusionGroupContext fusionCtx;
    const std::vector<std::vector<double>> posSource = cfg.INFO_Type ? std::vector<std::vector<double>>() : posMatrix;
    bool processOk = false;
    {
        TIMING_SCOPE(gmti_processing);
        processOk = cfg.enable_dbs_fusion
            ? proc.processPeriodsParallelFusion(periodList, cfg, posSource, fusionCtx, periodResults)
            : proc.processPeriodsParallel(periodList, cfg, posSource, periodResults);
    }
    if (!processOk) {
        std::cerr << "[ERR] 并行处理部分 period 失败（见日志）。\n";
        // 继续尝试从成功的结果中收集 MT
    }
    if (cfg.enable_dbs_fusion && processOk) {
        TIMING_SCOPE(dbs_processing);
        if (!runDbsFusionImaging(fusionCtx, cfg, true)) {
            std::cerr << "[WARN] DBS 融合成像输出失败，继续保留 GMTI 检测/CSV/标定流程。\n";
        }
    }
    
    // 按输入 periodList 的顺序收集 MT 以保证确定性
    std::vector<GMTIOutput::DetectionCsvRecord> detectionCsvRecords;
    for (size_t i = 0; i < periodResults.size(); ++i) {
        auto &res = periodResults[i];
        if (!res.MT.empty()) {
            all_frames_MT.push_back(res.MT);
            MT_acc.insert(MT_acc.end(), std::make_move_iterator(res.MT.begin()),
                          std::make_move_iterator(res.MT.end()));
        }
        detectionCsvRecords.insert(detectionCsvRecords.end(),
                                   res.detection_records.begin(),
                                   res.detection_records.end());
        // std::cout << "[OK] Period " << periodList[i] << " 完成\n";
    }

    // 整个扫描周期只写一次检测结果：2 + N*36
    bool wroteCurrentDetections = false;
    if (!MT_acc.empty()) {
        TIMING_SCOPE(result_write);
        if (!proc.writeResult(MT_acc, cfg)) {
            std::cerr << "[ERR] 整周期检测结果写盘失败\n";
        } else {
            wroteCurrentDetections = true;
        }
    } else {
        std::cerr << "[WARN] 本周期无检测目标，跳过检测结果写盘\n";
    }

    {
        TIMING_SCOPE(detection_result_csv_write);
        std::string sourceFile;
        if (cfg.result_file_id > 0 && cfg.result_file_id <= 99) {
            sourceFile = cfg.result_add;
            if (!sourceFile.empty() && sourceFile.back() != '/' && sourceFile.back() != '\\') {
                sourceFile.push_back('/');
            }
            char name[16];
            std::snprintf(name, sizeof(name), "GMTI%02d.bin", cfg.result_file_id);
            sourceFile += name;
        }
        if (!proc.writeDetectionCsv(detectionCsvRecords, cfg, sourceFile)) {
            std::cerr << "[ERR] detection_results.csv 写盘失败\n";
        }
    }
    writeFullPipelineCalibrationCheck(cfg, detectionCsvRecords);
    
    std::cout << "\n=== 检测定位阶段完成，开始航迹关联 ===\n";

    // ===== 阶段 2: 航迹关联 ===============================================
    // 调用 trackModule 进行航迹关联（从磁盘读取检测结果）
    // trackModule 参数：result_dir, idx_range, assoc_window
    if (!wroteCurrentDetections) {
        std::cout << "[TRACK][WARN] 本周期未写入新的检测结果文件，关联窗口中的当前编号文件可能是旧数据。" << std::endl;
    }
    TrackManager track_manager;
    {
        TIMING_SCOPE(tracking);
        trackModuleOnline(cfg, &track_manager);
    }
    
    // ===== 后续处理（可选）===============================================
    // 注：trackModule 已经将航迹关联结果写到磁盘
    // 如需进一步处理轨迹，可在这里添加代码

    std::cout << "\n--- GMTI 扫描完成 ---\n";
    gmti::runtime::finishRun(cfg, true, 0, processOk ? "normal exit" : "period processing reported failure");
    return 0;
}
