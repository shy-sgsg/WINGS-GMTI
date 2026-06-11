#ifndef CONFIG_STRUCTS_HPP
#define CONFIG_STRUCTS_HPP

#include <string>
#include <vector>
#include <complex>
#include <array>
#include <cstddef>

// 常量定义
const double EARTH_R = 6371.0;  // 地球半径，单位：km
const double C = 299792458.0;   // 光速，单位：m/s

// 配置结构体，存储从 XML 读取的所有配置参数
struct Config {
    // 文件路径相关
    std::string GMTI_Data_add;   // GMTI 数据路径1
    std::string GMTI_Data_add2;  // GMTI 数据路径2
    std::string GMTI_Data_new;   // 新协议 GMTI 数据路径
    std::string result_add;      // 结果文件夹路径
    std::string pipe_root_path = "/home/shy/pipe_test"; // POSIX 管道主路径
    std::string Plane_POS_add;   // 飞机位置文件路径
    std::string CAR_POS_add;     // 车辆位置文件路径
    std::string reffunc_add;     // 参考函数文件路径
    std::string channel_mode;    // 通道模式
    std::string iq_compose;      // IQ 组合方式
    std::vector<int> track_idx_range; // 航迹关联读取的周期索引
    int result_file_id = -1;     // 当前检测结果固定写入 GMTIxx.bin 的 xx，<=0 时才使用自动编号

    // 配置标志
    int INFO_Type;              // FPGA包头协议类型，1为新协议，0为旧协议
    int isPC;                   // 是否为PC，0表示不是，1表示是
    int hasRefFunc = 0;

    // 配置参数
    int info_len;               // 信息包长度
    int pulse_len;              // 脉冲长度
    int rg_len;                 // 距离采样长度
    int pulse_num;              // 脉冲数目
    int read_pulse_num = 0;     // 新协议实际读取脉冲数，<=0 表示读取 pulse_num
    int read_pulse_offset = -1; // 新协议读取起始偏移，<0 表示居中读取
    int process_pulse_num = 0;  // 当前处理矩阵脉冲数，<=0 表示使用 pulse_num
    int range_compress_len = 0; // 脉压/抽取后距离点数，>0 时覆盖 rg_len
    int pulse_dec;              // 脉冲压缩比例
    double fc;                  // 中心频率（GHz）
    double Br;                  // 带宽（MHz）
    double fs;                  // 采样频率（MHz）
    double Tr;                  // 脉冲宽度（秒）
    double PRF;                 // 脉冲重复频率（Hz）
    int az_count;               // 方位角计数
    double beamwidth_deg = 3.0; // 波束宽度，对应 XML boshu
    double loc_beam_gate_deg = -1.0; // 定位后波束方向门限半宽，<=0 时按 boshu 推导
    int week;                   // 周数
    double d_channel;           // 信道间隔
    double pf;                  // 固定参数
    double R_min;               // 最小距离
    double L0;                  // 参考经度
    double MT_nowz;             // 参考高度
    double calib_coef = 1.0;    // 校准系数
    int secBias;                // 偏差
    int skip_az_num;            // 跳过的方位脉冲数
    int Loc = 1;                // 是否使用飞机位置，1表示使用，0表示不使用

    // 新增参数
    int wavepos_st = 23;        // 波形位置起始索引
    int wavepos_ed = 27;        // 波形位置结束索引
    int wavepos_skip = 1;       // 波形位置跳过数
    int min_points = 11;        // 最小点数
    int min_len = 1;            // 最小轨迹长度

    // 跟踪输入帧顺序与调试控制
    bool mt_sort_by_utc = true; // true: 按 utc 重排; false: 保持原始 idx_range 顺序
    int track_debug_level = 0;  // 0: 关闭, 1: 输入帧摘要, 2: 关联细节
    int track_debug_frames = 5; // 仅打印前 N 帧调试信息, <=0 表示全部
    int track_debug_points = 3; // 每帧打印前 N 个点
    int track_idx_window = 4;   // 生成 track_idx_range 的窗长
    int track_truth_threshold = 3; // 判真阈值
    double track_gate_m = 500.0; // 航迹关联门限
    double track_v_max = 50.0;   // 航迹关联速度上限
    int track_confirm_window = 3; // 在线航迹 n屏选m确认窗口 N
    int track_confirm_hits = 2;   // 在线航迹 n屏选m确认命中数 M
    int track_max_missed = 2;     // Confirmed/Coasted 最大连续漏检保留帧数
    int track_tentative_max_missed = 1; // Tentative 最大连续漏检保留帧数
    double track_default_dt = 1.0; // UTC 无效时默认帧间隔，单位秒
    double track_chi2_gate = 9.21; // 二维量测 99% Mahalanobis 门限
    double track_tentative_gate_scale = 1.5; // Tentative 空间门限放宽倍率
    double track_tentative_chi2_scale = 1.5; // Tentative Mahalanobis 门限放宽倍率
    double track_dummy_cost = 1000.0; // Hungarian dummy 未匹配代价
    int track_linearity_window = 5; // 直线度评价窗口
    double track_min_linearity_confirm = 0.65; // Tentative 确认最小直线度
    double track_speed_smooth_weight = 0.5; // EN 差分速度平滑代价权重
    double track_heading_weight = 0.2; // EN 差分航向平滑代价权重
    double track_process_noise_pos = 25.0; // Kalman 位置过程噪声方差项
    double track_process_noise_vel = 10.0; // Kalman 速度过程噪声方差项
    double track_measurement_noise_pos = 50.0; // Kalman 位置量测噪声标准差
    bool track_debug_dump = true; // 在线 TrackManager 调试快照开关
    std::string track_debug_dir = ""; // 为空时使用 result_add/track_debug
    int track_debug_dump_level = 1; // 1: 周期摘要/CSV, 2: 逐航迹日志

    // 推导参数
    double lambda;              // 波长
    double R_bin;               // 距离采样分辨率
    std::vector<double> Rg;     // 距离向量
    double fd_res;              // 多普勒分辨率
    int rg_st, rg_ed, az_st, az_ed, az_center;  // 距离和方位支撑域
    int pkg_bytes;              // 每个脉冲数据包大小（字节）

    // ★ 新增：按 ROI 还是按 st..ed
    bool wavepos_use_roi = false;

    // ★ 新增：是否启用波位并行处理（true 表示并行，false 表示逐波位顺序处理）
    bool wavepos_parallel = true;
    bool enable_dbs_fusion = false;
    double dbs_out_res_m = 1.0; // DBS 成像输出分辨率，对应 XML raw_fenbianlv
    int dbs_beam_skip = 1;      // DBS 波位跳过数，对应 XML n_tiaoguo
    int dbs_range_skip = 1;     // DBS 距离抽样步长，对应 XML len_tiaoguo
    int dbs_interp_mode = 1;    // DBS 拼图插值模式，1=最近邻，2=双线性
    size_t dbs_max_mosaic_pixels = 200000000ULL; // DBS 拼图像素数上限，对应 XML dbs_max_mosaic_pixels

    // ★ 新增：是否估计误差角（false 时直接使用 XML 中的 squint_angle）
    bool estimate_error_angle = true;

    // ★ ROI 四角经纬度（度）：[lat1,lng1, lat2,lng2, lat3,lng3, lat4,lng4]
    std::array<double,4> roi_ll_deg;  // 调用前请填好

    // ★ 斜视角有效扫描范围（度），默认 [-25, 25]；离散波位由 wavepos_st/ed/skip 约束
    double scan_min_deg = -25.0;
    double scan_max_deg =  25.0;
    double lat_st;              
    double lat_ed;
    double lon_st;
    double lon_ed;

    int squint_side = 0;      // 右斜视侧，默认0
    double squint_angle = 0.0; // 斜视角度，单位：度
};

inline int effectivePulseNum(const Config& cfg)
{
    if (cfg.process_pulse_num > 0) {
        return cfg.process_pulse_num;
    }
    if (cfg.INFO_Type && cfg.read_pulse_num > 0) {
        return cfg.read_pulse_num;
    }
    return cfg.pulse_num;
}

// 定义结构体用于存储输出数据
struct GMTIOutput {
    // 原始数据
    std::vector<std::complex<double>> orgdata;

    // 时间戳
    double utcMid;

    // 飞机位置
    struct Plane {
        double V;    // 速度
        double E;    // 东向位置
        double N;    // 北向位置
        double H;    // 高度
        double V_angle;  // 飞机航向角
    } plane;

    // 数据处理结果
    struct Detect {
        std::vector<int> prow;
        std::vector<int> pcol;
        std::vector<double> row_af;
        std::vector<std::vector<double>> MTpos;
    } detect;

    // 多普勒频率
    std::vector<double> fa_c;

    // 相位
    std::vector<double> phi_diss;
    std::vector<double> phi_fit;

    // CSI 结果
    std::vector<std::complex<double>> GMTIclean;
    std::vector<std::complex<double>> echoF1;

    //
    std::vector<double> MT;
};

// 建议放到 config_structs.hpp 或 GMTIOutput 定义处
struct TargetSelection {
    std::vector<int> prow;                       // 行索引
    std::vector<int> pcol;                       // 列索引
    std::vector<double> A;                       // 幅度
    std::vector<std::complex<double>> ch1_Data;  // 通道1复数
    std::vector<std::complex<double>> ch2_Data;  // 通道2复数
};

struct Track {
    std::array<double,4> x{};
    std::array<double,16> P{};
    std::vector<std::array<double,2>> pos;
    std::vector<std::array<double,2>> kf;
    std::vector<std::array<double,4>> x_state;
    std::vector<double> time;
    double direction = 0.0;
    double range = 0.0;
    int last = 0;
    int missed = 0;
    int id = 0;
};

#endif // CONFIG_STRUCTS_HPP
