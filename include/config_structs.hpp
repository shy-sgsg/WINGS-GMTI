#ifndef CONFIG_STRUCTS_HPP
#define CONFIG_STRUCTS_HPP

#include <string>
#include <vector>
#include <complex>
#include <array>
#include <cstddef>
#include <limits>

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
    std::string iq_data_type = "float32"; // 原始 IQ 数据类型：float32 / int16
    int new_protocol_channel_count = 2;    // 新协议交织通道总数
    int new_protocol_read_channel_1 = 1;   // 新协议读取通道1（1-based）
    int new_protocol_read_channel_2 = 2;   // 新协议读取通道2（1-based）
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
    int range_fft_len = 0;      // 距离脉压 FFT 长度，<=0 表示使用 pulse_len
    int range_crop_start = 0;   // 距离脉压 IFFT 后截取起点，0-based
    int pulse_dec;              // 脉冲压缩比例
    double fc;                  // 中心频率（GHz）
    double Br;                  // 带宽（MHz）
    double fs;                  // 采样频率（MHz）
    double Tr;                  // 脉冲宽度（秒）
    double PRF;                 // 脉冲重复频率（Hz）
    bool has_sample_delay_us = false; // XML 显式提供采样/接收延迟时为 true
    double sample_delay_us = 0.0;      // 统一为 us；未提供时 dump 输出 null
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
    int wavepos_st = 0;        // 波形位置起始索引
    int wavepos_ed = 25;        // 波形位置结束索引
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
    // 在线 TrackManager 关联模式。旧滑窗 trackModule() 不使用这些模式开关。
    int track_distance_mode = 1;   // 0=Euclidean, 1=MahalanobisSquared
    int track_assignment_mode = 1; // 0=GreedyNearestNeighbor, 1=Hungarian
    bool track_use_distance_cost = true;
    bool track_use_speed_cost = true;
    bool track_use_heading_cost = true;
    bool track_use_detection_speed_cost = false;
    double track_distance_weight = 1.0;
    double track_detection_speed_weight = 0.0;
    bool track_use_euclidean_gate = true;
    bool track_use_mahalanobis_gate = true;
    bool track_use_max_speed_gate = true;
    bool track_use_heading_gate = false;
    bool track_use_detection_speed_gate = false;
    double track_gate_m = 300.0; // 欧氏物理距离门限，m
    double track_v_max = 100.0;  // 瞬时速度门限，m/s
    double track_heading_gate_deg = 90.0;
    double track_detection_speed_gate_mps = 30.0;
    double track_euclidean_cost_scale_m = 0.0;
    double track_mahalanobis_cost_scale = 0.0;
    double track_speed_cost_scale_mps = 0.0;
    double track_heading_cost_scale_deg = 180.0;
    double track_detection_speed_cost_scale_mps = 0.0;
    int track_confirm_window = 3; // 在线航迹 n屏选m确认窗口 N
    int track_confirm_hits = 2;   // 在线航迹 n屏选m确认命中数 M
    int track_max_missed = 2;     // Confirmed/Coasted 最大连续漏检保留帧数
    int track_tentative_max_missed = 1; // Tentative 最大连续漏检保留帧数
    double track_default_dt = 1.0; // UTC 无效时默认帧间隔，单位秒
    double track_chi2_gate = 9.21; // 二维量测 99% Mahalanobis 门限
    double track_tentative_gate_scale = 1.5; // Tentative 空间门限放宽倍率
    double track_tentative_chi2_scale = 1.5; // Tentative Mahalanobis 门限放宽倍率
    double track_dummy_cost = 1.5; // 真实匹配必须优于该漏检代价
    double track_invalid_cost = 1.0e12; // 有限大无效代价，避免 Hungarian 使用 inf
    bool track_allow_equal_dummy_cost = false;
    int track_linearity_window = 5; // 直线度评价窗口
    double track_min_linearity_confirm = 0.65; // Tentative 确认最小直线度
    double track_speed_smooth_weight = 0.0; // EN 差分速度平滑代价权重
    double track_heading_weight = 0.0; // EN 差分航向平滑代价权重
    double track_process_noise_pos = 25.0; // Kalman 位置过程噪声方差项
    double track_process_noise_vel = 10.0; // Kalman 速度过程噪声方差项
    double track_measurement_noise_pos = 50.0; // Kalman 位置量测噪声标准差
    bool track_debug_dump = true; // 在线 TrackManager 调试快照开关
    std::string track_debug_dir = ""; // 为空时使用 result_add/track_debug
    int track_debug_dump_level = 1; // 1: 周期摘要/CSV, 2: 逐航迹日志
    std::string runtime_mode = "debug"; // debug=输出可追溯诊断; release=实时运行优先
    bool runtime_diagnostics_enabled = true; // release 模式下默认由 loadXML 关闭
    bool debug_pc_peak = false; // P1.5: 输出脉压后、对消前单点峰值检查
    std::string pc_peak_scene_truth; // 可选 scene_truth.csv 路径；为空时按 result_add 推断
    bool motion_comp_enable = false; // 动目标定位时扣除目标自身运动多普勒
    bool motion_comp_analytic_enable = true; // true: 使用解析式运动补偿；false: 回退旧定位
    bool motion_comp_use_row_doppler = true; // true: af_total 使用检测所在 Doppler row
    std::string motion_comp_solver = "analytic"; // old / iterative / analytic / root1d / debug
    int motion_comp_iter = 8; // 运动补偿迭代次数
    double motion_comp_iter_tol_mps = 1.0e-4; // 迭代收敛阈值
    bool p38_refit_enable = true; // CFAR/聚类后是否执行 p38 二阶段重估
    int p38_refit_row_guard_bins = 2; // 检测点 row 方向排除半径
    int p38_refit_range_guard_bins = 2; // 检测点 range 方向排除半径
    double p38_refit_top_power_frac = 0.01; // 按 power_map 排除的顶部百分位
    int p38_refit_min_sample_count = 8; // p38_refit 最小样本行数
    double p38_refit_min_inlier_ratio = 0.60; // p38_refit 最小内点占比
    double p38_refit_max_rmse_rad = 0.60; // p38_refit 最大 RMSE
    double p38_refit_max_delta_k = 0.01; // p38_refit 与 p38_pre 的 k 限幅
    double p38_refit_max_delta_b_rad = 1.50; // p38_refit 与 p38_pre 的 b 限幅
    int ati_velocity_sign = 1; // ATI 相位到径向速度符号，必要时可设为 -1
    int ati_phase_to_velocity_sign = 1; // ATI 相位残差到径向速度的符号
    int motion_doppler_axis_sign = -1; // 径向速度到多普勒轴方向的符号
    double ati_phase_bias_rad = 0.0; // ATI 固定相位偏置
    double ati_vmax_mps = 1.6; // ATI 径向速度限幅，<=0 表示不限幅
    double motion_comp_denom_min = 1.0e-6; // 解析式分母过小则回退旧定位
    double motion_comp_root_grid_step_mps = 0.02; // root1d 粗搜索步长
    double motion_comp_root_cost_max = 0.25; // root1d 最大可接受 cost
    bool motion_comp_debug = false; // 输出运动多普勒补偿调试字段/日志

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
    bool enable_dbs_fusion = true;
    double dbs_out_res_m = 25.0; // DBS 成像输出分辨率，对应 XML raw_fenbianlv
    int dbs_beam_skip = 1;      // DBS 波位跳过数，对应 XML n_tiaoguo
    int dbs_range_skip = 1;     // DBS 距离抽样步长，对应 XML len_tiaoguo
    int dbs_interp_mode = 1;    // DBS 拼图插值模式，1=最近邻，2=双线性
    size_t dbs_max_mosaic_pixels = 200000000ULL; // DBS 拼图像素数上限，对应 XML dbs_max_mosaic_pixels
    double dbs_mosaic_margin_ratio = 0.05; // DBS 拼图范围比例裕量，对应 XML dbs_mosaic_margin_ratio
    double dbs_mosaic_margin_m = 200.0;    // DBS 拼图范围最小裕量（米），对应 XML dbs_mosaic_margin_m

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

inline int effectiveRangeFftLen(const Config& cfg)
{
    return cfg.range_fft_len > 0 ? cfg.range_fft_len : cfg.pulse_len;
}

inline int effectiveRangeCompressLen(const Config& cfg)
{
    return cfg.range_compress_len > 0 ? cfg.range_compress_len : cfg.rg_len;
}

inline bool usesRangeCropWindow(const Config& cfg)
{
    return effectiveRangeFftLen(cfg) != cfg.pulse_len ||
           cfg.range_crop_start != 0;
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

    struct DetectionCsvRecord {
        int period_id = -1;
        int beam_id = -1;
        double platform_e = std::numeric_limits<double>::quiet_NaN();
        double platform_n = std::numeric_limits<double>::quiet_NaN();
        double platform_h = std::numeric_limits<double>::quiet_NaN();
        double platform_v = std::numeric_limits<double>::quiet_NaN();
        double platform_v_angle_deg = std::numeric_limits<double>::quiet_NaN();
        double fd_ctr_wrapped = std::numeric_limits<double>::quiet_NaN();
        double fd_ctr_unwrapped = std::numeric_limits<double>::quiet_NaN();
        int range_bin = -1;
        int row = -1;
        int col = -1;
        double range_m = std::numeric_limits<double>::quiet_NaN();
        double theta_cmd_deg = std::numeric_limits<double>::quiet_NaN();
        double theta_true_deg = std::numeric_limits<double>::quiet_NaN();
        double e = std::numeric_limits<double>::quiet_NaN();
        double n = std::numeric_limits<double>::quiet_NaN();
        double lat = std::numeric_limits<double>::quiet_NaN();
        double lon = std::numeric_limits<double>::quiet_NaN();
        double utc = std::numeric_limits<double>::quiet_NaN();
        double amplitude = std::numeric_limits<double>::quiet_NaN();
        double radial_velocity_mps = std::numeric_limits<double>::quiet_NaN();
        double phase_rad = std::numeric_limits<double>::quiet_NaN();
        double p38_k = std::numeric_limits<double>::quiet_NaN();
        double p38_b = std::numeric_limits<double>::quiet_NaN();
        double phi_static_model_rad = std::numeric_limits<double>::quiet_NaN();
        std::string phi_static_model_name;
        double C_ati = std::numeric_limits<double>::quiet_NaN();
        double k_eff_static_phase_df = std::numeric_limits<double>::quiet_NaN();
        double phi_static_rad = std::numeric_limits<double>::quiet_NaN();
        double phi_static_total_rad = std::numeric_limits<double>::quiet_NaN();
        double phi_res_rad = std::numeric_limits<double>::quiet_NaN();
        double phi_static_at_zero = std::numeric_limits<double>::quiet_NaN();
        double phi_res_at_zero = std::numeric_limits<double>::quiet_NaN();
        double phi_static_geometry_rad = std::numeric_limits<double>::quiet_NaN();
        double af_phase = std::numeric_limits<double>::quiet_NaN();
        double af_total = std::numeric_limits<double>::quiet_NaN();
        double af_geometry = std::numeric_limits<double>::quiet_NaN();
        double af_motion = std::numeric_limits<double>::quiet_NaN();
        double phi_motion = std::numeric_limits<double>::quiet_NaN();
        double delta_t_s = std::numeric_limits<double>::quiet_NaN();
        double motion_comp_denom = std::numeric_limits<double>::quiet_NaN();
        double denom_without_k = std::numeric_limits<double>::quiet_NaN();
        double v_from_phase_raw = std::numeric_limits<double>::quiet_NaN();
        double v_from_phi_res = std::numeric_limits<double>::quiet_NaN();
        double v_iterative_mps = std::numeric_limits<double>::quiet_NaN();
        double v_analytic_mps = std::numeric_limits<double>::quiet_NaN();
        double v_root1d_mps = std::numeric_limits<double>::quiet_NaN();
        double v_old_mps = std::numeric_limits<double>::quiet_NaN();
        double af_geometry_old_hz = std::numeric_limits<double>::quiet_NaN();
        double af_geometry_iterative_hz = std::numeric_limits<double>::quiet_NaN();
        double af_geometry_analytic_hz = std::numeric_limits<double>::quiet_NaN();
        double af_geometry_root1d_hz = std::numeric_limits<double>::quiet_NaN();
        double root1d_cost = std::numeric_limits<double>::quiet_NaN();
        double p38_pre_k = std::numeric_limits<double>::quiet_NaN();
        double p38_pre_b = std::numeric_limits<double>::quiet_NaN();
        double p38_pre_rmse = std::numeric_limits<double>::quiet_NaN();
        double p38_refit_k = std::numeric_limits<double>::quiet_NaN();
        double p38_refit_b = std::numeric_limits<double>::quiet_NaN();
        double p38_refit_rmse = std::numeric_limits<double>::quiet_NaN();
        int p38_refit_sample_count = 0;
        double p38_refit_inlier_ratio = std::numeric_limits<double>::quiet_NaN();
        int p38_refit_valid = 0;
        double p38_used_k = std::numeric_limits<double>::quiet_NaN();
        double p38_used_b = std::numeric_limits<double>::quiet_NaN();
        std::string p38_used_source;
        double phi_static_pre_rad = std::numeric_limits<double>::quiet_NaN();
        double phi_res_pre_rad = std::numeric_limits<double>::quiet_NaN();
        double v_pre_mps = std::numeric_limits<double>::quiet_NaN();
        double phi_static_refit_rad = std::numeric_limits<double>::quiet_NaN();
        double phi_res_refit_rad = std::numeric_limits<double>::quiet_NaN();
        double v_refit_mps = std::numeric_limits<double>::quiet_NaN();
        double sinA_old = std::numeric_limits<double>::quiet_NaN();
        double sinA_comp = std::numeric_limits<double>::quiet_NaN();
        double sinA_used = std::numeric_limits<double>::quiet_NaN();
        double angle_from_sinA_deg = std::numeric_limits<double>::quiet_NaN();
        double theta_used_for_position_deg = std::numeric_limits<double>::quiet_NaN();
        double look_from_sinA_e = std::numeric_limits<double>::quiet_NaN();
        double look_from_sinA_n = std::numeric_limits<double>::quiet_NaN();
        double look_e_diff = std::numeric_limits<double>::quiet_NaN();
        double look_n_diff = std::numeric_limits<double>::quiet_NaN();
        double old_e = std::numeric_limits<double>::quiet_NaN();
        double old_n = std::numeric_limits<double>::quiet_NaN();
        double new_e = std::numeric_limits<double>::quiet_NaN();
        double new_n = std::numeric_limits<double>::quiet_NaN();
        int old_valid = 0;
        int comp_valid = 0;
        int old_invalid_comp_valid = 0;
        int motion_comp_valid = 0;
        int motion_comp_enable = 0;
        int motion_comp_used = 0;
        int motion_comp_fallback = 0;
        int p38_theory_sign = 0;
        int motion_doppler_axis_sign = 0;
        int ati_phase_to_velocity_sign = 0;
        std::string p38_mode;
        std::string geometry_calib_mode;
        std::string loc_used_mode;
        std::string motion_comp_status;
        std::string motion_comp_solver;
    };
    std::vector<DetectionCsvRecord> detection_records;
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
