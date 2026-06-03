#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <complex>
#include <cmath>
#include <chrono>
#include <algorithm> // lower_bound
#include <numeric>   // accumulate
#include <cstring>
#include <limits>
//#include <Eigen/Dense>
#include "GMTIProcessor.hpp"
#include "rangeCompress.hpp"
#include "geo/geoProj.hpp"
#include "rotation_xy.hpp"
#include "unwrap_fd.hpp"

using cudacd = cuFloatComplex;
using cd = std::complex<float>;

inline int sgn(double x, double eps = 1e-6) { return (x > eps) - (x < -eps); }
// 返回同样编号：1=SW,2=NE,3=NW,4=SE,0=未知
int flight_flag_by_sign(double Vx, double Vy, double eps = 1e-6)
{
    int sE = sgn(Vx, eps), sN = sgn(Vy, eps);
    if (sE == 0 && sN == 0)
        return 0;
    if (sN < 0 && sE < 0)
        return 1; // SW
    if (sN > 0 && sE > 0)
        return 2; // NE
    if (sN > 0 && sE < 0)
        return 3; // NW
    if (sN < 0 && sE > 0)
        return 4; // SE
    // 轴上用另一分量决定
    if (sE == 0)
        return (sN > 0) ? 3 : 4; // 正北→NW，正南→SE
    if (sN == 0)
        return (sE > 0) ? 2 : 1; // 正东→NE，正西→SW
    return 0;
}

static void compare_cpu_gpu_vector(const std::vector<std::complex<double>> &cpu,
                                   const std::vector<std::complex<double>> &gpu,
                                   const std::string &tag,
                                   double tol = 1e-8)
{
    if (cpu.size() != gpu.size()) {
        std::cerr << "[TEST] " << tag << " size mismatch: " << cpu.size() << " vs " << gpu.size() << "\n";
        return;
    }

    double maxAbs = 0.0;
    double maxRel = 0.0;
    size_t maxIdx = 0;
    for (size_t i = 0; i < cpu.size(); ++i) {
        std::complex<double> diff = cpu[i] - gpu[i];
        double absErr = std::abs(diff);
        double relErr = (std::abs(cpu[i]) > 1e-16) ? absErr / std::abs(cpu[i]) : absErr;
        if (absErr > maxAbs) {
            maxAbs = absErr;
            maxRel = relErr;
            maxIdx = i;
        }
    }

    std::cout << "[TEST] " << tag << " -> maxAbs=" << maxAbs << ", maxRel=" << maxRel
              << ", maxIdx=" << maxIdx << "\n";

    if (maxAbs > tol) {
        std::cout << "[TEST] " << tag << " FAILED at idx=" << maxIdx << ", cpu=" << cpu[maxIdx]
                  << ", gpu=" << gpu[maxIdx] << "\n";
    } else {
        std::cout << "[TEST] " << tag << " PASSED (tol=" << tol << ")\n";
    }
}

namespace {

static inline uint16_t load_u16_le(const uint8_t *p)
{
    return static_cast<uint16_t>(p[0] | (static_cast<uint16_t>(p[1]) << 8));
}

static inline int16_t load_i16_le(const uint8_t *p)
{
    return static_cast<int16_t>(load_u16_le(p));
}

static inline uint32_t load_u32_le(const uint8_t *p)
{
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

static inline uint64_t load_u64_le(const uint8_t *p)
{
    return static_cast<uint64_t>(p[0]) |
           (static_cast<uint64_t>(p[1]) << 8) |
           (static_cast<uint64_t>(p[2]) << 16) |
           (static_cast<uint64_t>(p[3]) << 24) |
           (static_cast<uint64_t>(p[4]) << 32) |
           (static_cast<uint64_t>(p[5]) << 40) |
           (static_cast<uint64_t>(p[6]) << 48) |
           (static_cast<uint64_t>(p[7]) << 56);
}

static inline float load_f32_le(const uint8_t *p)
{
    uint32_t raw = load_u32_le(p);
    float v = 0.0f;
    std::memcpy(&v, &raw, sizeof(float));
    return v;
}

static inline double load_f64_le(const uint8_t *p)
{
    uint64_t raw = load_u64_le(p);
    double v = 0.0;
    std::memcpy(&v, &raw, sizeof(double));
    return v;
}

static bool readPulseBlockNewProtocol(const Config &cfg,
                                      int beamskip,
                                      std::vector<std::complex<double>> &data1,
                                      std::vector<std::complex<double>> &data2,
                                      std::vector<double> &utc,
                                      double &theta_sq,
                                      std::vector<std::vector<double>> &posRaw)
{
    constexpr size_t kHeaderBytes = 256;
    constexpr size_t kSamplesPerPrt = 4096;
    constexpr size_t kBytesPerSample = 16; // ch1(I,Q) float32 + ch2(I,Q) float32
    constexpr size_t kPrtBytes = kHeaderBytes + kSamplesPerPrt * kBytesPerSample;

    if (cfg.pulse_len != static_cast<int>(kSamplesPerPrt)) {
        std::cerr << "[ERR] 新协议要求 pulse_len=4096，当前 pulse_len=" << cfg.pulse_len << std::endl;
        return false;
    }

    const size_t W = static_cast<size_t>(cfg.pulse_num);

    const std::string &echoPath = cfg.GMTI_Data_new.empty() ? cfg.GMTI_Data_add : cfg.GMTI_Data_new;
    std::ifstream fp(echoPath, std::ios::binary);
    if (!fp) {
        std::cerr << "[ERR] 无法打开新协议回波文件: " << echoPath << std::endl;
        return false;
    }
    
    fp.seekg(0, std::ios::end);
    const std::streamsize file_size = fp.tellg();
    if (file_size <= 0 || (static_cast<uint64_t>(file_size) % kPrtBytes) != 0U) {
        std::cerr << "[ERR] 新协议回波文件大小非法，不是完整 PRT 包整数倍: " << echoPath << std::endl;
        return false;
    }

    const size_t total_prt = static_cast<size_t>(static_cast<uint64_t>(file_size) / kPrtBytes);
    const size_t start_with_skip = (static_cast<size_t>(beamskip) - 1) * W + static_cast<size_t>(cfg.skip_az_num);
    const size_t start_without_skip = (static_cast<size_t>(beamskip) - 1) * W;

    size_t start_prt = start_with_skip;
    if (start_with_skip + W > total_prt) {
        if (start_without_skip + W <= total_prt) {
            std::cerr << "[WARN] 新协议读取带 skip_pulses 越界，自动按已裁剪文件读取（忽略 skip_pulses）: period="
                      << beamskip << " skip_pulses=" << cfg.skip_az_num << std::endl;
            start_prt = start_without_skip;
        } else {
            std::cerr << "[ERR] 新协议回波文件大小不足，无法读取 period=" << beamskip
                      << " (total_prt=" << total_prt
                      << ", need_end_prt=" << (start_with_skip + W)
                      << ", alt_need_end_prt=" << (start_without_skip + W) << ")" << std::endl;
            return false;
        }
    }

    const size_t start_byte = start_prt * kPrtBytes;
    if (static_cast<uint64_t>(file_size) < start_byte + W * kPrtBytes) {
        std::cerr << "[ERR] 新协议回波文件大小不足，无法读取 period=" << beamskip << std::endl;
        return false;
    }
    fp.seekg(static_cast<std::streamoff>(start_byte), std::ios::beg);

    data1.resize(W * kSamplesPerPrt);
    data2.resize(W * kSamplesPerPrt);
    utc.resize(W);
    posRaw.assign(W, std::vector<double>(7, 0.0));

    std::vector<uint8_t> packet(kPrtBytes);
    std::vector<double> fw_angle_deg(W, 0.0);

    for (size_t k = 0; k < W; ++k) {
        fp.read(reinterpret_cast<char *>(packet.data()), static_cast<std::streamsize>(kPrtBytes));
        if (fp.gcount() != static_cast<std::streamsize>(kPrtBytes)) {
            std::cerr << "[ERR] 读取新协议 PRT 失败, period=" << beamskip << " pulse=" << k << std::endl;
            return false;
        }

        const uint8_t *hdr = packet.data();
        utc[k] = static_cast<double>(load_f32_le(hdr + 16));
        posRaw[k][0] = utc[k];
        // New protocol stores latitude/longitude in degrees, while extractPlanePos expects radians.
        posRaw[k][1] = load_f64_le(hdr + 104) * M_PI / 180.0;
        posRaw[k][2] = load_f64_le(hdr + 112) * M_PI / 180.0;
        posRaw[k][3] = load_f64_le(hdr + 120);
        posRaw[k][4] = static_cast<double>(load_f32_le(hdr + 128));
        posRaw[k][5] = static_cast<double>(load_f32_le(hdr + 132));
        posRaw[k][6] = static_cast<double>(load_f32_le(hdr + 136));
        fw_angle_deg[k] = static_cast<double>(load_i16_le(hdr + 218)) / 100.0;

        const uint8_t *payload = hdr + kHeaderBytes;
        for (size_t n = 0; n < kSamplesPerPrt; ++n) {
            const size_t off = n * kBytesPerSample;
            const float ch1_i = load_f32_le(payload + off + 0);
            const float ch1_q = load_f32_le(payload + off + 4);
            const float ch2_i = load_f32_le(payload + off + 8);
            const float ch2_q = load_f32_le(payload + off + 12);
            data1[k * kSamplesPerPrt + n] = std::complex<double>(ch1_i, ch1_q);
            data2[k * kSamplesPerPrt + n] = std::complex<double>(ch2_i, ch2_q);
        }
    }

    theta_sq = fw_angle_deg[fw_angle_deg.size() / 2];
    
    return true;
}

static inline double normalize_azimuth_deg(double angle_deg)
{
    angle_deg = std::fmod(angle_deg, 360.0);
    if (angle_deg < 0.0)
        angle_deg += 360.0;
    return angle_deg;
}

static inline double wrap180_deg(double angle_deg)
{
    angle_deg = std::fmod(angle_deg + 180.0, 360.0);
    if (angle_deg < 0.0)
        angle_deg += 360.0;
    return angle_deg - 180.0;
}

} // namespace

double GMTIProcessor::estimateSquintAngleDeg(const GMTIOutput::Plane &plane, const Config &cfg, double fd_ctr)
{
    const double lambda = (cfg.lambda > 0.0) ? cfg.lambda : (C / (cfg.fc * 1e9));
    const double v = plane.V;

    if (std::isfinite(fd_ctr))
    {
        if (v <= 0.0 || lambda <= 0.0)
        {
            std::cout << "[SQUINT] fd_ctr=" << fd_ctr
                      << ", v=" << v
                      << ", lambda=" << lambda
                      << ", squint=0 (invalid v/lambda)" << std::endl;
            return 0.0;
        }

        const double ratio = fd_ctr * lambda / (2.0 * v);
        const double clipped = std::max(-1.0, std::min(1.0, ratio));
        const double squint_deg = std::asin(clipped) * 180.0 / M_PI;

        std::cout << "[SQUINT] fd_ctr=" << fd_ctr
                  << ", v=" << v
                  << ", lambda=" << lambda
                  << ", ratio=" << ratio
                  << ", squint=" << squint_deg
                  << std::endl;
        return squint_deg;
    }

    double refE = 0.0;
    double refN = 0.0;
    Gaussp3(cfg.roi_ll_deg[0], cfg.roi_ll_deg[1], cfg.L0, refE, refN);
    const double dE = refE - plane.E;
    const double dN = refN - plane.N;
    const double bearing_deg = std::atan2(dN, dE) * 180.0 / M_PI;
    const double squint_deg = wrap180_deg(bearing_deg - plane.V_angle + 90.0);
    
    std::cout << "[SQUINT] fd_ctr=nan"
              << ", v=" << v
              << ", lambda=" << lambda
              << ", refE=" << refE
              << ", refN=" << refN
              << ", planeE=" << plane.E
              << ", planeN=" << plane.N
              << ", bearing=" << bearing_deg
              << ", V_angle=" << plane.V_angle
              << ", squint=" << squint_deg
              << std::endl;
    
    return squint_deg;
}

double GMTIProcessor::estimateSquintAngleDeg(const std::vector<std::complex<double>> &data,
                                             const GMTIOutput::Plane &plane,
                                             const Config &cfg)
{
    double fd_ctr = 0.0;
    int start_pulse = 0;
    int window_pulses = 0;

    if (estimateCenterFdCtrFromData(data, cfg, fd_ctr, start_pulse, window_pulses))
    {
        std::cout << "[SQUINT] center-window start=" << start_pulse
                  << ", count=" << window_pulses
                  << ", fd_ctr=" << fd_ctr
                  << std::endl;
        return estimateSquintAngleDeg(plane, cfg, fd_ctr);
    }

    std::cout << "[SQUINT] center-window fd_ctr estimate failed, fallback to geometric estimate" << std::endl;
    return estimateSquintAngleDeg(plane, cfg);
}

// 辅助函数：进行脉压处理（支持GPU加速 + CPU回退）
bool GMTIProcessor::pulseCompression(std::vector<std::complex<double>> &data1,
                                     std::vector<std::complex<double>> &data2,
                                     const Config &cfg)
{
    TIMING_SCOPE(pulseCompression);
    // 对两个通道分别进行脉压
    if (data1.size() != (size_t)cfg.pulse_num * cfg.pulse_len ||
        data2.size() != (size_t)cfg.pulse_num * cfg.pulse_len)
    {
        std::cerr << "输入数据尺寸不匹配脉冲数和脉冲长度。" << std::endl;
        return false;
    }

    const int Lraw = cfg.pulse_len;
    const int M = cfg.rg_len;
    const int W = cfg.pulse_num;
    
    // #define FORCE_CPU_RANGE_COMPRESS 1 // 取消注释此行强制使用CPU脉压（调试用）
    bool try_gpu = true; // 默认尝试GPU
    #ifdef FORCE_CPU_RANGE_COMPRESS
    try_gpu = false;
    #endif

    // ========== 尝试GPU路径 ==========
    if (try_gpu && gpu_ptrs_.d1 != nullptr) {
        std::vector<std::complex<float>> data1_f(data1.size());
        std::vector<std::complex<float>> data2_f(data2.size());
        
        // 转换为float
        for (size_t i = 0; i < data1_f.size(); ++i) {
            data1_f[i] = std::complex<float>(static_cast<float>(data1[i].real()), static_cast<float>(data1[i].imag()));
            data2_f[i] = std::complex<float>(static_cast<float>(data2[i].real()), static_cast<float>(data2[i].imag()));
        }

        // 上传float数据到GPU
        size_t total_bytes = data1_f.size() * sizeof(std::complex<float>);
        if (cudaMemcpyAsync(gpu_ptrs_.d1, data1_f.data(), total_bytes, cudaMemcpyHostToDevice, stream_compute_) == cudaSuccess &&
            cudaMemcpyAsync(gpu_ptrs_.d2, data2_f.data(), total_bytes, cudaMemcpyHostToDevice, stream_compute_) == cudaSuccess)
        {
            // 调用GPU脉压设备路径：通道1
            if (rangeCompressCUFFT_device(Lraw, M, cfg)) {
                // GPU处理通道1成功，现在处理通道2
                // 首先交换d1和d2以处理通道2
                std::swap(gpu_ptrs_.d1, gpu_ptrs_.d2);
                if (rangeCompressCUFFT_device(Lraw, M, cfg)) {
                    // 都成功，交换回来，准备下载
                    std::swap(gpu_ptrs_.d1, gpu_ptrs_.d2);
                    
                    // 下载结果
                    std::vector<std::complex<float>> result1_f(W * M);
                    std::vector<std::complex<float>> result2_f(W * M);
                    size_t download_bytes = W * M * sizeof(std::complex<float>);
                    
                    if (cudaMemcpyAsync(result1_f.data(), gpu_ptrs_.d1, download_bytes, cudaMemcpyDeviceToHost, stream_compute_) == cudaSuccess &&
                        cudaMemcpyAsync(result2_f.data(), gpu_ptrs_.d2, download_bytes, cudaMemcpyDeviceToHost, stream_compute_) == cudaSuccess)
                    {
                        cudaStreamSynchronize(stream_compute_); // 等待下载完成
                        
                        // 转换回double并更新data1/data2
                        data1.resize(W * M);
                        data2.resize(W * M);
                        for (int i = 0; i < W * M; ++i) {
                            data1[i] = std::complex<double>(result1_f[i].real(), result1_f[i].imag());
                            data2[i] = std::complex<double>(result2_f[i].real(), result2_f[i].imag());
                        }
                        std::cout << "[GPU] 脉压处理成功 (GPU cuFFT)" << std::endl;
                        DBG("[GPU] 脉压处理成功 (GPU cuFFT)");
                        return true;
                    }
                }
                // 如果通道2失败，交换回来然后回退到CPU
                std::swap(gpu_ptrs_.d1, gpu_ptrs_.d2);
            }
        }
        std::cout << "[GPU] 脉压处理失败，回退到CPU" << std::endl;
        DBG("[GPU] 脉压处理失败，回退到CPU");
    }

    // ========== CPU回退路径 ==========
    std::vector<std::complex<double>> temp; // 一个通用中间缓存

    // --- 通道1 ---
    bool success1 = rangeCompressFFT(cfg, data1, temp);
    if (!success1)
        return false;
    data1.swap(temp); // 替换通道1结果

    // --- 通道2 ---
    success1 = rangeCompressFFT(cfg, data2, temp);
    if (!success1)
        return false;
    data2.swap(temp); // 替换通道2结果

    if (data1.size() != (size_t)cfg.pulse_num * cfg.rg_len ||
        data2.size() != (size_t)cfg.pulse_num * cfg.rg_len)
    {
        std::cerr << "输出数据尺寸不匹配脉冲数和距离长度。" << std::endl;
        return false;
    }

    DBG("[CPU] 脉压处理成功 (CPU FFTW)");
    return true;
}

// 辅助函数：计算多普勒频率轴
inline bool computeDoppler(const std::vector<std::complex<double>> &data, int k, const Config &cfg,
                           std::vector<double> &faAxis, double &fa2, const double &v, const double &theta_deg)
{
    // 获取数据的尺寸
    size_t Na = cfg.pulse_num; // 方位向点数
    size_t Nr = cfg.rg_len;    // 距离向点数

    // 相关矩阵初始化
    std::vector<std::complex<double>> R_m(Na, {0.0, 0.0});
    std::vector<std::complex<double>> R_tem(Nr, {0.0, 0.0});

    // 计算相关矩阵
    for (size_t mm = k; mm < Na; ++mm)
    {
        std::complex<double> R_m1 = {0.0, 0.0};
        for (size_t nn = 0; nn < Nr; ++nn)
        {
            // 计算相关系数（复数乘法）
            size_t index_mm = mm * Nr + nn;      // 计算一维数组的索引
            size_t index_k = (mm - k) * Nr + nn; // 偏移后的索引

            R_m1 += data[index_mm] * std::conj(data[index_k]);
            R_tem[nn] = R_m1;
        }
        R_m[mm] = R_m1 / static_cast<double>(Nr); // 归一化
    }

    // 计算最终的相关值
    std::complex<double> R_m2 = {0.0, 0.0};
    for (size_t kk = k; kk < Na; ++kk)
    {
        R_m2 += R_m[kk];
    }

    // 计算多普勒频率
    fa2 = (cfg.PRF / (2 * M_PI)) * std::arg(R_m2); // 使用复数的相位（`arg`）

    // 生成频率轴 faAxis
    faAxis.clear();
    for (double f = -cfg.PRF / 2 + cfg.fd_res; f <= cfg.PRF / 2; f += cfg.fd_res)
    {
        faAxis.push_back(f); // 填充频率轴
    }

    // fa2 = unwrap_prf_to_model(fa2 , cfg.PRF, theta_deg, v, cfg.fc); // 解除模糊

    // 将 fa2 加到 faAxis 上
    for (size_t i = 0; i < faAxis.size(); ++i)
    {
        faAxis[i] += fa2;
    }

    return true; // 返回 true 表示成功
}


bool estimateCenterFdCtrFromData(const std::vector<std::complex<double>> &data,
                                        const Config &cfg,
                                        double &fd_ctr,
                                        int &start_pulse,
                                        int &window_pulses)
{
    fd_ctr = 0.0;
    start_pulse = 0;
    window_pulses = 0;

    if (cfg.pulse_num < 2 || cfg.rg_len <= 0) {
        return false;
    }

    const int center_start = std::max(0, (cfg.pulse_num - 2) / 2);
    const int center_count = std::min(2, cfg.pulse_num);
    const size_t nr = static_cast<size_t>(cfg.rg_len);
    const size_t slice_size = static_cast<size_t>(center_count) * nr;
    if (data.size() < static_cast<size_t>(cfg.pulse_num) * nr || slice_size == 0) {
        return false;
    }

    std::vector<std::complex<double>> slice(slice_size);
    for (int p = 0; p < center_count; ++p) {
        const size_t src_base = static_cast<size_t>(center_start + p) * nr;
        const size_t dst_base = static_cast<size_t>(p) * nr;
        std::copy_n(data.begin() + static_cast<std::ptrdiff_t>(src_base), nr, slice.begin() + static_cast<std::ptrdiff_t>(dst_base));
    }

    Config local_cfg = cfg;
    local_cfg.pulse_num = center_count;
    local_cfg.fd_res = (center_count > 0) ? (cfg.PRF / double(center_count)) : cfg.fd_res;

    std::vector<double> faAxis;
    double fa_tmp = 0.0;
    if (!computeDoppler(slice, 1, local_cfg, faAxis, fa_tmp, 0.0, 0.0)) {
        return false;
    }

    fd_ctr = fa_tmp;
    start_pulse = center_start;
    window_pulses = center_count;
    return true;
}

// 处理一个周期的数据
bool GMTIProcessor::processOnePeriod(int periodIdx, const Config &cfg_, const std::vector<std::vector<double>> &posRaw, GMTIOutput &out)
{
    TIMING_SCOPE(processOnePeriod);
    Config cfg = cfg_; // 复制配置，局部修改
    // 0) 读取脉冲块（适配双文件/交织）
    std::vector<std::complex<double>> data1, data2;
    std::vector<double> utc;
    std::vector<std::uint8_t> headers;
    std::vector<std::vector<double>> echoPosRaw;

    // 读取脉冲块
    double theta_sq = 0.0; // 方位角平方
    bool readSuccess = false;
    if (cfg.INFO_Type) {
        readSuccess = readPulseBlockNewProtocol(cfg, periodIdx, data1, data2, utc, theta_sq, echoPosRaw);
    } else {
        readSuccess = readPulseBlock(cfg, periodIdx, data1, data2, utc, theta_sq);
    }
    if (!readSuccess)
    {
        std::cerr << "读取脉冲块失败。" << std::endl;
        return false;
    }
    DBG("读取脉冲块成功，脉冲数: " << utc.size());
    

    // 如果没有脉压，则需要脉压
    if (!cfg.isPC)
    {
        // 进行脉压处理（纯 CPU 路径）
        bool pc_ok = pulseCompression(data1, data2, cfg);
        if (!pc_ok) {
            std::cerr << "脉压处理失败。" << std::endl;
            return false;
        }
        DBG("脉压处理成功");
    }
    else
    {
        // 只进行抽取
        // 归一化 + 裁剪到 M，输出到 rc_out（W×M 行主）
        std::vector<std::complex<double>> rc_out;
        const int Lraw = cfg.pulse_len;
        const int W = cfg.pulse_num;
        const int M1 = cfg.rg_len;
        const int Lraw2M = Lraw / M1;
        rc_out.resize((size_t)W * M1);
        for (int k = 0; k < W; ++k)
        {
            for (int m = 0; m < M1; ++m)
            {
                rc_out[(size_t)k * M1 + m] = data1[(size_t)k * Lraw + m * Lraw2M];
            }
        }
        data1.swap(rc_out);

        rc_out.resize((size_t)W * M1);
        for (int k = 0; k < W; ++k)
        {
            for (int m = 0; m < M1; ++m)
            {
                rc_out[(size_t)k * M1 + m] = data2[(size_t)k * Lraw + m * Lraw2M];
            }
        }
        data2.swap(rc_out);
        DBG("跳过脉压，直接抽取数据成功");
    }

    // 脉压完成后再更新采样率，避免影响 rangeCompressFFT / cuFFT 的匹配滤波器构造
    cfg.fs = cfg.fs * cfg.rg_len / cfg.pulse_len;
    cfg.pulse_len = cfg.rg_len; // 更新脉冲长度为距离采样长度
    // 通用方位向抽取：cfg.pulse_dec 合 1 抽取
    if (data1.size() != (size_t)cfg.pulse_num * cfg.rg_len ||
        data2.size() != (size_t)cfg.pulse_num * cfg.rg_len)
    {
        ERR("数据尺寸不匹配脉冲数和距离长度，无法继续处理。");
        return false;
    }

    const int W_orig = cfg.pulse_num;
    const int M = cfg.rg_len;
    const int dec = cfg.pulse_dec;
    if (dec <= 0)
    {
        ERR("cfg.pulse_dec 必须大于 0（当前 = " << dec << "）。");
        return false;
    }
    if (W_orig % dec != 0)
    {
        ERR("cfg.pulse_num 不是" << dec << "的倍数，无法进行" << dec << "合 1 抽取（当前 pulse_num = " << W_orig << "）。");
        return false;
    }

    const int W_new = W_orig / dec;
    if (W_new == W_orig)
    {
        DBG("cfg.pulse_dec == 1，跳过抽取。");
    }
    else
    {
        // 临时缓冲，复用（data1 先处理，后处理 data2）
        std::vector<std::complex<double>> tmp((size_t)W_new * M);

        // 处理 data1
        for (int k_new = 0; k_new < W_new; ++k_new)
        {
            size_t out_base = (size_t)k_new * M;
            size_t in_base = (size_t)(dec * k_new) * M;
            for (int m = 0; m < M; ++m)
            {
                std::complex<double> acc = 0.0;
                for (int d = 0; d < dec; ++d)
                {
                    acc += data1[in_base + (size_t)d * M + m];
                }
                tmp[out_base + m] = acc;
            }
        }
        data1.swap(tmp);

        // 复用 tmp 处理 data2（重新分配确保大小正确）
        tmp.assign((size_t)W_new * M, std::complex<double>(0.0, 0.0));
        for (int k_new = 0; k_new < W_new; ++k_new)
        {
            size_t out_base = (size_t)k_new * M;
            size_t in_base = (size_t)(dec * k_new) * M;
            for (int m = 0; m < M; ++m)
            {
                std::complex<double> acc = 0.0;
                for (int d = 0; d < dec; ++d)
                {
                    acc += data2[in_base + (size_t)d * M + m];
                }
                tmp[out_base + m] = acc;
            }
        }
        data2.swap(tmp);
    }

    // ========== 抽取 utc 同步 ==========
    if ((int)utc.size() != W_orig)
    {
        ERR("utc 长度" << utc.size() << " 与 cfg.pulse_num(" << W_orig << ") 不匹配，无法同步抽取。");
        return false;
    }
    std::vector<double> utc_tmp(W_new);
    for (int k_new = 0; k_new < W_new; ++k_new)
    {
        // 平均法（推荐）
        double acc = 0.0;
        for (int d = 0; d < dec; ++d)
            acc += utc[dec * k_new + d];
        utc_tmp[k_new] = acc / (double)dec;

        // 若想使用取中值法（整数下标），可替换为：
        // utc_tmp[k_new] = utc[dec * k_new + dec/2];
    }
    utc.swap(utc_tmp);

    // 更新脉冲数
    cfg.pulse_num = W_new;
    cfg.PRF = cfg.PRF / dec; // 更新 PRF
    DBG("方位向" << dec << " 合 1 抽取完成：脉冲数 " << W_orig << " -> " << cfg.pulse_num << "，距离长度 M = " << M << "。");

    // 从帧头提取 UTC
    double utc_mean = 0.0;
    for (int i = 0; i < utc.size(); ++i)
    {
        utc_mean += utc[i];
    }
    utc_mean /= utc.size();

    for (int i = 0; i < utc.size(); ++i)
    {
        if (utc_mean - utc[i] > 0.6)
        {
            utc[i] += 1.0; // UTC 时间向前调整1秒
        }
        else if (utc_mean - utc[i] < -0.6)
        {
            utc[i] -= 1.0; // UTC 时间向后调整1秒
        }
    }

    out.utcMid = utc[utc.size() / 2];


    // 通道2 乘以系数，补偿幅度差异
    double coef = cfg.calib_coef; // 可根据实际情况调整
    const size_t N = data2.size();
    for (size_t i = 0; i < N; ++i)
        data2[i] *= coef;

    // 1) 飞机位姿/速度
    GMTIOutput::Plane plane;
    if (cfg.Loc)
    {
        const std::vector<std::vector<double>> &planePosSource = cfg.INFO_Type ? echoPosRaw : posRaw;
        bool planePosExtracted = extractPlanePos(utc, planePosSource, cfg, plane);
        if (!planePosExtracted)
        {
            std::cerr << "无法提取飞机位姿。" << std::endl;
            return false;
        }
        DBG("飞机位置提取成功: E=" << plane.E << " N=" << plane.N << " H=" << plane.H << " V=" << plane.V);
    }
    else
    {
        plane.V = 40 / 3.6; // 设置默认速度
    }

    // 2) 多普勒中心 & 动态支撑域
    double fa_ctr;
    std::vector<double> faAxis;
    double angle_deg = GMTIProcessor::estimateSquintAngleDeg(data1, plane, cfg);
    double theta_deg = angle_deg + theta_sq; // 斜视角（GMTI 内部估计）+ 帧内方位修正
    DBG("GMTI estimated angle_deg = " << angle_deg << " deg");

    computeDoppler(data1, 1, cfg, faAxis, fa_ctr, plane.V, theta_deg);
    DBG("多普勒中心计算成功: fa_ctr=" << fa_ctr);

    // 动态支撑域
    int az_center = 0, az_st = 0, az_ed = 0, rg_st = 0, rg_ed = 0;
    double fd_st = 0.0, fd_ed = 0.0, BW_az = 0.0;

    bool domainCalculated = computeDynamicSupportDomain(faAxis, fa_ctr, plane, cfg,
                                                        az_center, az_st, az_ed,
                                                        fd_st, fd_ed, BW_az,
                                                        rg_st, rg_ed);

    if (!domainCalculated)
    {
        std::cerr << "动态支撑域计算失败。" << std::endl;
        return false;
    }
    DBG("动态支撑域计算成功: az_st=" << az_st << " az_ed=" << az_ed << " rg_st=" << rg_st << " rg_ed=" << rg_ed);

    cfg.az_center = az_center;
    cfg.az_st = az_st;
    cfg.az_ed = az_ed;
    cfg.rg_st = rg_st;
    cfg.rg_ed = rg_ed;

    // 3) 通道对齐
    double baseline = cfg.d_channel; // 沿航迹基线长度，单位米
    double V = plane.V;     // 飞机速度，单位米/秒（日志里是50m/s）
    double PRF = cfg.PRF;   // 脉冲重复频率，单位Hz（之前算的是2000Hz）

    double delta_t = baseline / V; // 时间延迟，单位秒
    int skipInt_theory = static_cast<int>(std::round(delta_t * PRF));
    DBG("理论 skipInt = " << skipInt_theory);

    int skipInt = 1;
    std::vector<std::complex<double>> F1, F2, F1_r, F2_r;
    const size_t Na = static_cast<size_t>(cfg.pulse_num);
    const size_t Nr = static_cast<size_t>(cfg.rg_len);

    std::vector<std::complex<float>> data1_f(Na * Nr);
    std::vector<std::complex<float>> data2_f(Na * Nr);
    for (size_t i = 0; i < data1_f.size(); ++i) {
        data1_f[i] = std::complex<float>(static_cast<float>(data1[i].real()), static_cast<float>(data1[i].imag()));
        data2_f[i] = std::complex<float>(static_cast<float>(data2[i].real()), static_cast<float>(data2[i].imag()));
    }

    cuda_upload_async(data1_f, data2_f, Na, Nr);
    cuda_stage_align_async(skipInt, Na, Nr);
    cuda_stage_fft_async(Na, Nr);
    cuda_stage_dbs_async(static_cast<float>(fa_ctr), Na, Nr);

    double fa2 = unwrap_prf_to_model(fa_ctr, cfg.PRF, theta_deg, plane.V, cfg.fc); // 解除模糊
    double fa_shift = fa2 - fa_ctr;
    for (size_t i = 0; i < faAxis.size(); ++i) faAxis[i] += fa_shift;
    
    // --- 4) 距离向初相校正 ---
    double thre_rg = 0.1 * M_PI;
    std::vector<float> phi_fit;
    std::vector<double> cpu_phi_fit;
    std::vector<double> phi_diss_phase; // 存储相位值
    std::vector<int>    phi_diss_range; // 存储对应的距离向索引

    if (rg_correct_CUDA(cfg, thre_rg, phi_fit, phi_diss_phase, phi_diss_range)) {
        cuda_apply_rg_correction_async(phi_fit, Na, Nr);
    }
    else {
        std::cout << "相位校正失败。" << std::endl;
        return false;
    }
    DBG("最终对齐完成，使用移位脉冲数: " << skipInt);

    std::array<double, 2> p_38;
    if (!clutter_cancel_38_paper_1_p38_cuda(
        faAxis,
        cfg.az_st, cfg.rg_st, cfg.az_ed, cfg.rg_ed,
        cfg,
        p_38
    )) return false;

    // 计算用于 CFAR 的 phase_map（GPU 生成，避免回传大矩阵）
    std::vector<float> phase_map;
    if (!cuda_download_phase_map(phase_map, Na * Nr)) return false;

    // 对消
    cuda_stage_align_async(skipInt_theory, Na, Nr);
    cuda_stage_fft_async(Na, Nr);
    cuda_stage_dbs_async(static_cast<float>(fa_ctr), Na, Nr);
    if (rg_correct_CUDA(cfg, thre_rg, phi_fit, phi_diss_phase, phi_diss_range)) {
        cuda_apply_rg_correction_async(phi_fit, Na, Nr);
    }
    else {
        std::cout << "相位校正失败。" << std::endl;
        return false;
    }
     DBG("最终对齐完成，使用移位脉冲数: " << skipInt_theory);

    // --- 8) 最终 CSI 对消 ---
    std::array<double, 2> p_38_csi;
    std::vector<double> ph_trace;
    std::vector<double> fa_cut;
    if (!clutter_cancel_38_paper_1_cuda(
        faAxis,
        cfg.az_st, cfg.rg_st, cfg.az_ed, cfg.rg_ed,
        cfg,
        F1_r, p_38_csi, ph_trace, fa_cut
    )) return false;

    // --- 10) CFAR 处理 ---
    const int band_st = std::max(0, std::min(az_st, cfg.pulse_num - 1));
    const int band_ed = std::max(0, std::min(az_ed, cfg.pulse_num - 1));
    std::vector<int> prow, pcol;
    std::vector<float> mydata;
    int cfar_bnum = 16; // 背景单元数
    int c_num = 4;      // 保护单元数
    const bool use_gpu_cfar = true;
    bool cfarSuccess1 = false;
    bool cfarSuccess2 = false;
 
    std::vector<std::complex<float>> CSI_out;
    cfarSuccess1 = dpca_cfar2_fast_cuda(CSI_out, band_st, band_ed, cfg.pf, c_num, cfar_bnum, "GO", cfg, mydata);
    if (!cfarSuccess1)
    {
        std::cerr << "GPU CFAR 处理失败，尝试回退到 CPU 路径。" << std::endl;
        // CPU 回退路径：需要回传 F1/F2
        std::vector<std::complex<float>> F1_f, F2_f;
        cuda_download_sync(F1_f, F2_f, Na * Nr);
        F1.resize(F1_f.size());
        F2.resize(F2_f.size());
        for (size_t i = 0; i < F1_f.size(); ++i) {
            F1[i] = std::complex<double>(F1_f[i].real(), F1_f[i].imag());
            F2[i] = std::complex<double>(F2_f[i].real(), F2_f[i].imag());
        }
        if (!cuda_download_csi_sync(CSI_out, Na * Nr))
        {
            std::cerr << "CSI 回传失败。" << std::endl;
            return false;
        }
        std::vector<std::complex<double>> detect_data = F2;
        if (band_st <= band_ed)
        {
            for (int r = band_st; r <= band_ed; ++r)
            {
                const size_t off = static_cast<size_t>(r) * static_cast<size_t>(cfg.rg_len);
                for (size_t c = 0; c < static_cast<size_t>(cfg.rg_len); ++c) {
                    const auto v = CSI_out[off + c];
                    detect_data[off + c] = std::complex<double>(v.real(), v.imag());
                }
            }
        }
        std::vector<double> mydata_d;
        cfarSuccess2 = dpca_cfar2_fast(detect_data, cfg.pf, c_num, cfar_bnum, "GO", cfg, mydata_d, prow, pcol);
        if (cfarSuccess2) {
            mydata.assign(mydata_d.begin(), mydata_d.end());
        }
        if (cfarSuccess2) {
            DBG("合成检测数据完成(CPU)：杂波带内[" << band_st << ":" << band_ed << "]用CSI_out，带外用F2");
        }        
        else {
            std::cerr << "CPU CFAR 处理也失败。" << std::endl;
            return false;
        }
    }

    // --- 11) 聚类滤波 ---
    std::vector<int> prow_new, pcol_new;
    std::vector<float> refined_mydata;
    std::vector<float> phase_std_list;
    // cluster_filter(mydata, cfg.min_points, cfg, refined_mydata, prow_new, pcol_new);
    const bool use_gpu_cluster = true;
    bool cluster_ok = false;
    if (use_gpu_cluster)
    {
        cluster_ok = cluster_filter_gap_phase_cuda(mydata, phase_map, cfg.min_points, 2, 0.2f,
                               cfg, refined_mydata, prow_new, pcol_new, phase_std_list);
        if (!cluster_ok)
        {
            std::cerr << "GPU 聚类失败，回退到 CPU 路径。" << std::endl;
        }
    }
    if (!cluster_ok)
    {
        std::vector<double> mydata_d(mydata.begin(), mydata.end());
        std::vector<double> phase_map_d(phase_map.begin(), phase_map.end());
        std::vector<double> refined_mydata_d;
        std::vector<double> phase_std_list_d;
        cluster_ok = cluster_filter_gap_phase(mydata_d, phase_map_d, cfg.min_points, 2, 0.2,
                                              cfg, refined_mydata_d, prow_new, pcol_new, phase_std_list_d);
        if (cluster_ok) {
            refined_mydata.assign(refined_mydata_d.begin(), refined_mydata_d.end());
            phase_std_list.assign(phase_std_list_d.begin(), phase_std_list_d.end());
        }
    }

#ifdef DEBUG
    // 保存 CFAR 结果以便调试
    {
        std::ofstream fout("debug_CFI_mydata.bin", std::ios::binary);
        fout.write(reinterpret_cast<const char *>(mydata.data()), mydata.size() * sizeof(float));
        fout.close();
    }
#endif

    if (!cfarSuccess1 && !cfarSuccess2)
    {
        std::cerr << "CFAR 处理失败。" << std::endl;
        return false;
    }
    DBG("CFAR 处理成功，目标数: " << prow_new.size());


#ifdef DEBUG
    {
        std::ofstream fout1("debug_F1_final.bin", std::ios::binary);
        fout1.write(reinterpret_cast<const char *>(F1.data()), F1.size() * sizeof(std::complex<double>));
        fout1.close();

        std::ofstream fout2("debug_F2_final.bin", std::ios::binary);
        fout2.write(reinterpret_cast<const char *>(F2.data()), F2.size() * sizeof(std::complex<double>));
        fout2.close();
    }
#endif
    GMTIOutput::Detect targetSel;
    const bool use_gpu_target = true;
    bool targetDetection = false;
    if (use_gpu_target)
    {
        targetDetection = target_select_cuda(prow_new, pcol_new, cfg, targetSel);
        if (!targetDetection)
        {
            std::cerr << "GPU target_select 失败，回退到 CPU 路径。" << std::endl;
        }
    }
    if (!targetDetection)
    {
        std::vector<std::complex<float>> F1_f, F2_f;
        cuda_download_sync(F1_f, F2_f, Na * Nr);
        F1.resize(F1_f.size());
        F2.resize(F2_f.size());
        for (size_t i = 0; i < F1_f.size(); ++i) {
            F1[i] = std::complex<double>(F1_f[i].real(), F1_f[i].imag());
            F2[i] = std::complex<double>(F2_f[i].real(), F2_f[i].imag());
        }
        targetDetection = target_select(F1, F2, prow_new, pcol_new, cfg, targetSel);
    }

    double ref_phase = faAxis[az_center] * p_38[0] + p_38[1];

    auto deg2rad = [](double d)
    { return d * M_PI / 180.0; };

    // 结果容器
    const size_t L = targetSel.prow.size();
    std::vector<int> row_af(L);
    std::vector<double> af_ransac(L);
    std::vector<double> MT; // MT(i,:) = [lat,lng,cfg.MT_nowz,xP,yP]
    // MT.resize(L * 6, 0.0);

    // 参考相位：ref_phase 已在外部给定
    const double k = p_38[0]; // slope
    const double b = p_38[1]; // intercept
    DBG("多普勒斜率 k = " << k << ", 截距 b = " << b);
    const double thetaRot = plane.V_angle; // 度
    DBG("旋转角 thetaRot = " << thetaRot << " 度");
    const double cosT = std::abs(std::cos(deg2rad(thetaRot)));
    const double sinT = std::abs(std::sin(deg2rad(thetaRot)));
    const double VE = plane.V * std::cos(deg2rad(plane.V_angle));
    const double VN = plane.V * std::sin(deg2rad(plane.V_angle));

    int flag = flight_flag_by_sign(VE, VN); // 1/2/3/4
    flag += cfg_.squint_side * 4;           // 根据斜视侧调整方向标志

    for (size_t i = 0; i < L; ++i)
    {
        const int r = targetSel.prow[i];
        const int c = targetSel.pcol[i];
        const size_t off = static_cast<size_t>(r) * Nr + static_cast<size_t>(c);

        // dphi = angle(F1f * conj(F2f))
        // const std::complex<double> v = F1[off] * std::conj(F2[off]);
        // double dphi = std::arg(v);
        double dphi = static_cast<double>(phase_map[off]); 

        // dphi = phaseUnwrap(dphi, ref_phase);  —— 相对 ref_phase 的最短距离解缠
        {
            double diff = dphi - ref_phase;
            if (diff > M_PI)
                diff -= 2.0 * M_PI;
            if (diff < -M_PI)
                diff += 2.0 * M_PI;
            dphi = ref_phase + diff;
        }

        // af_ransac = (dphi - b) / k
        if (std::abs(k) < 1e-12)
        {
            af_ransac[i] = 0.0; // 防除零，可按需改成 continue/报错
        }
        else
        {
            af_ransac[i] = (dphi - b) / k;
        }

        // row_af(i) = round((af_ransac - faAxis(1)) / fd_res) + 1
        // —— 这里改为 0-based：不再 +1
        // row_af[i] = static_cast<int>(std::llround((af_ransac[i] - faAxis.front()) / cfg.fd_res));
        // 支撑域剔除
        // if (row_af[i] < az_st + 10 || row_af[i] > az_ed - 10) {
            // DBG("目标 " << i << " 因超出支撑域被剔除: row_af=" << row_af[i]
            //             << " 有效区间=[" << az_st + 10 << "," << az_ed - 10 << "]");
            // continue;
        // }

        // 几何定位（需要 plane.V, plane.H, plane.E, plane.N；cfg.Rg 为距离轴）
        const double sinA = (af_ransac[i]) * cfg.lambda / (2.0 * plane.V);
        const double Rg = cfg.Rg[static_cast<size_t>(c)];
        const double py = Rg * sinA;

        const double dz = (plane.H - cfg.MT_nowz);
        const double px2 = Rg * Rg - py * py - dz * dz;
        if (px2 < 0.0) {
            continue;
            DBG("目标 " << i << " 几何解算失败 (虚数根 px2 < 0): " << px2);
        }
        const double px = (px2 > 0.0) ? std::sqrt(px2) : 0.0; // 物理约束，负值置 0

        double xP, yP;

        rotation_xy(py, px, flag, cosT, sinT, xP, yP);

        xP += plane.E;
        yP += plane.N;

        double lat = 0.0, lng = 0.0;               // 由投影坐标反解经纬
        (void)Gaussp3RV(xP, yP, cfg.L0, lat, lng); // 假定返回 bool，忽略失败则保持 0

        const double dE = xP - plane.E;
        const double dN = yP - plane.N;
        const double target_azimuth_deg = normalize_azimuth_deg(std::atan2(dN, dE) * 180.0 / M_PI);  // 东为0°, 逆时针为正
        // 相对方向：顺时针为正 => 计算 plane - target，再映射到 [-180,180]
        const double direction = wrap180_deg(plane.V_angle - target_azimuth_deg);
        const double range = std::sqrt(dE * dE + dN * dN);
        
        // debug printing removed

        // MT(i,:) = [lat, lng, cfg.MT_nowz, xP, yP, utc, relative_dir, range]
        MT.push_back(lat);
        MT.push_back(lng);
        MT.push_back(cfg.MT_nowz);
        MT.push_back(xP);
        MT.push_back(yP);
        MT.push_back(out.utcMid);
        MT.push_back(direction);
        MT.push_back(range);
        DBG("目标 " << i << " 定位成功: Lat=" << lat << " Lng=" << lng);
    }

    out.detect = targetSel;
    out.MT = MT;

    if (!targetDetection)
    {
        std::cerr << "动目标定位失败。" << std::endl;
        return false;
    }
//    DBG("动目标定位成功");

    return true; // 所有处理成功，返回 true
}

bool GMTIProcessor::extractPlanePos(const std::vector<double> &t_utc,
                                    const std::vector<std::vector<double>> &POS, // [POS_num][17]
                                    const Config &cfg,
                                    GMTIOutput::Plane &plane)
{
    // 1) 退化条件
    const bool empty_or_zero = t_utc.empty() ||
                               std::all_of(t_utc.begin(), t_utc.end(), [](double x)
                                           { return x == 0.0; });

    if (empty_or_zero)
    {
        plane.E = plane.N = 0.0;
        plane.H = cfg.MT_nowz;
        plane.V = 40.0 / 3.6; // 40 km/h -> m/s
        plane.V_angle = 0.0;
        return true;
    }

    // 2) 基本检查
    const size_t pos_num = POS.size();
    if (pos_num == 0 || POS[0].size() < 4)
        return false;

    // 拿出 POS 时间轴（秒）
    std::vector<double> t_pos(pos_num);
    for (size_t i = 0; i < pos_num; ++i)
        t_pos[i] = POS[i][0];

    // 3) 逐个 t_utc 做线性插值 -> lat/lon/alt（弧度/弧度/米）
    const size_t M = t_utc.size();
    std::vector<double> lat_rad(M), lon_rad(M), alt_m(M);
    lat_rad.reserve(M);
    lon_rad.reserve(M);
    alt_m.reserve(M);

    auto interp_scalar = [&](double t, int col) -> double
    {
        // 边界：落在两端时取端点
        if (t <= t_pos.front())
            return POS.front()[col];
        if (t >= t_pos.back())
            return POS.back()[col];

        // lower_bound 找到第一个 >= t 的位置
        auto it = std::lower_bound(t_pos.begin(), t_pos.end(), t);
        size_t ir = size_t(it - t_pos.begin());
        size_t il = ir - 1;

        double t1 = t_pos[il], t2 = t_pos[ir];
        double y1 = POS[il][col], y2 = POS[ir][col];
        double r = (t - t1) / (t2 - t1);
        return y1 + r * (y2 - y1);
    };

    for (size_t i = 0; i < M; ++i)
    {
        double lat_val = interp_scalar(t_utc[i], 1);
        double lon_val = interp_scalar(t_utc[i], 2);

        // Compatibility: old POS usually stores radians; some new inputs may contain degrees.
        lat_rad[i] = (std::abs(lat_val) > M_PI) ? (lat_val * M_PI / 180.0) : lat_val;
        lon_rad[i] = (std::abs(lon_val) > M_PI) ? (lon_val * M_PI / 180.0) : lon_val;
        alt_m[i] = interp_scalar(t_utc[i], 3);   // 列4：高度(米)
    }

    // 4) 经纬 -> 投影 EN
    std::vector<double> E(M), N(M);
    for (size_t i = 0; i < M; ++i)
    {
        double lat_deg = lat_rad[i] * 180.0 / M_PI;
        double lon_deg = lon_rad[i] * 180.0 / M_PI;
        double e = 0.0, n = 0.0;
        Gaussp3(lat_deg, lon_deg, cfg.L0, e, n); // 若你的实现返回bool
        E[i] = e;
        N[i] = n;
    }

    // 5) 平均高度/位置
    auto mean = [](const std::vector<double> &v) -> double
    {
        if (v.empty())
            return 0.0;
        double s = std::accumulate(v.begin(), v.end(), 0.0);
        return s / double(v.size());
    };

    plane.H = mean(alt_m);
    plane.E = mean(E);
    plane.N = mean(N);

    // 6) 速度估计（与 MATLAB 一致：Δ/脉冲数 * PRF）
    //   注：若 t_utc 非等间隔，用下行替代：
    //   const double dt = (t_utc.back() - t_utc.front()); // 秒
    //   plane.Vx = (E.back() - E.front()) / dt; ...
    const double scale = (cfg.pulse_num > 0) ? (cfg.PRF / double(cfg.pulse_num)) : 0.0;
    double Vx = (E.back() - E.front()) * scale;
    double Vy = (N.back() - N.front()) * scale;
    double Vz = (alt_m.back() - alt_m.front()) * scale;

    plane.V = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
    plane.V_angle = std::atan2(Vy, Vx) * 180.0 / M_PI; // 东北平面速度方向角(°)

    DBG("飞机速度分量: Vx=" << Vx << " Vy=" << Vy << " Vz=" << Vz);

    return true;
}

bool GMTIProcessor::extractPlanePVFromEcho(const Config &cfg,
                                           GMTIOutput::Plane &plane)
{
    std::vector<std::complex<double>> data1, data2;
    std::vector<double> utc;
    std::vector<std::vector<double>> echoPosRaw;
    double theta_sq = 0.0;

    if (!readPulseBlockNewProtocol(cfg, 1, data1, data2, utc, theta_sq, echoPosRaw)) {
        std::cerr << "extractPlanePVFromEcho: failed to read new-protocol echo" << std::endl;
        return false;
    }

    if (!extractPlanePos(utc, echoPosRaw, cfg, plane)) {
        std::cerr << "extractPlanePVFromEcho: failed to extract plane from embedded echo pose" << std::endl;
        return false;
    }

    return true;
}

// Compute dataset-level squint by reading center period(s) and estimating fd_ctr
// from center-window data. This replicates DBS semantics: estimate fd_ctr from
// the central period(s), unwrap PRF ambiguity, convert to squint angle and
// average when two centers exist.
bool GMTIProcessor::computeDatasetSquintFromCenter(const std::vector<int> &periodList,
                                                   const Config &cfg,
                                                   const std::vector<std::vector<double>> &posRaw,
                                                   double &out_squint)
{
    out_squint = 0.0;
    if (periodList.empty()) return false;

    // determine center period(s)
    std::vector<int> centers;
    size_t n = periodList.size();
    if (n % 2 == 1) {
        centers.push_back(periodList[n/2]);
    } else {
        centers.push_back(periodList[n/2 - 1]);
        centers.push_back(periodList[n/2]);
    }

    std::vector<double> angles_deg;
    for (int per : centers) {
        std::vector<std::complex<double>> data1, data2;
        std::vector<double> utc;
        std::vector<std::vector<double>> echoPosRaw;
        double theta_sq_local = 0.0;

        bool readOk = false;
        if (cfg.INFO_Type) {
            readOk = readPulseBlockNewProtocol(cfg, per, data1, data2, utc, theta_sq_local, echoPosRaw);
        } else {
            readOk = readPulseBlock(cfg, per, data1, data2, utc, theta_sq_local);
        }
        if (!readOk) {
            std::cerr << "computeDatasetSquintFromCenter: failed to read period " << per << std::endl;
            return false;
        }

        // Range-compress if necessary
        Config local_cfg = cfg;
        if (!local_cfg.isPC) {
            if (!pulseCompression(data1, data2, local_cfg)) {
                std::cerr << "computeDatasetSquintFromCenter: pulseCompression failed for period " << per << std::endl;
                return false;
            }
        } else {
            // Extract/normalize path (same as processOnePeriod extraction)
            const int Lraw = local_cfg.pulse_len;
            const int W = local_cfg.pulse_num;
            const int M1 = local_cfg.rg_len;
            const int Lraw2M = Lraw / M1;
            std::vector<std::complex<double>> rc_out((size_t)W * M1);
            for (int k = 0; k < W; ++k) {
                for (int m = 0; m < M1; ++m) {
                    rc_out[(size_t)k * M1 + m] = data1[(size_t)k * Lraw + m * Lraw2M];
                }
            }
            data1.swap(rc_out);
            rc_out.assign((size_t)W * M1, std::complex<double>(0.0,0.0));
            for (int k = 0; k < W; ++k) {
                for (int m = 0; m < M1; ++m) {
                    rc_out[(size_t)k * M1 + m] = data2[(size_t)k * Lraw + m * Lraw2M];
                }
            }
            data2.swap(rc_out);
        }

        // After range compression/extraction, update sampling parameters like processOnePeriod does
        local_cfg.fs = local_cfg.fs * local_cfg.rg_len / local_cfg.pulse_len;
        local_cfg.pulse_len = local_cfg.rg_len;

        // Apply azimuth decimation (pulse_dec) same as processOnePeriod
        const int dec = local_cfg.pulse_dec;
        if (dec <= 0) {
            std::cerr << "computeDatasetSquintFromCenter: invalid pulse_dec" << std::endl;
            return false;
        }
        const int W_orig = local_cfg.pulse_num;
        if (W_orig % dec != 0) {
            std::cerr << "computeDatasetSquintFromCenter: pulse_num not divisible by pulse_dec" << std::endl;
            return false;
        }
        const int W_new = W_orig / dec;
        if (W_new != W_orig) {
            // decimate data1 in-place
            std::vector<std::complex<double>> tmp((size_t)W_new * local_cfg.rg_len);
            for (int k_new = 0; k_new < W_new; ++k_new) {
                size_t out_base = (size_t)k_new * local_cfg.rg_len;
                size_t in_base = (size_t)(dec * k_new) * local_cfg.rg_len;
                for (int m = 0; m < local_cfg.rg_len; ++m) {
                    std::complex<double> acc = 0.0;
                    for (int d = 0; d < dec; ++d) acc += data1[in_base + (size_t)d * local_cfg.rg_len + m];
                    tmp[out_base + m] = acc;
                }
            }
            data1.swap(tmp);
            tmp.assign((size_t)W_new * local_cfg.rg_len, std::complex<double>(0.0,0.0));
            for (int k_new = 0; k_new < W_new; ++k_new) {
                size_t out_base = (size_t)k_new * local_cfg.rg_len;
                size_t in_base = (size_t)(dec * k_new) * local_cfg.rg_len;
                for (int m = 0; m < local_cfg.rg_len; ++m) {
                    std::complex<double> acc = 0.0;
                    for (int d = 0; d < dec; ++d) acc += data2[in_base + (size_t)d * local_cfg.rg_len + m];
                    tmp[out_base + m] = acc;
                }
            }
            data2.swap(tmp);
            local_cfg.pulse_num = W_new;
            local_cfg.PRF = local_cfg.PRF / dec;
        }

        // Extract plane for this period (use echoPosRaw if available)
        GMTIOutput::Plane plane;
        bool plane_ok = false;
        if (!echoPosRaw.empty()) {
            plane_ok = extractPlanePos(utc, echoPosRaw, local_cfg, plane);
        } else if (!posRaw.empty()) {
            plane_ok = extractPlanePos(utc, posRaw, local_cfg, plane);
        } else {
            plane_ok = extractPlanePVFromEcho(local_cfg, plane);
        }
        if (!plane_ok) {
            std::cerr << "computeDatasetSquintFromCenter: failed to extract plane for period " << per << std::endl;
            return false;
        }

        // Estimate center fd_ctr from data (center window)
        double fd_ctr = 0.0;
        int start_pulse = 0, window_pulses = 0;
        if (!estimateCenterFdCtrFromData(data1, local_cfg, fd_ctr, start_pulse, window_pulses)) {
            std::cerr << "computeDatasetSquintFromCenter: center fd estimate failed for period " << per << std::endl;
            return false;
        }

        // Unwrap PRF ambiguity using model
        double fd_unwrapped = unwrap_prf_to_model(fd_ctr, local_cfg.PRF, theta_sq_local, plane.V, local_cfg.fc);

        // Convert to squint angle (deg)
        const double lambda = (local_cfg.lambda > 0.0) ? local_cfg.lambda : (3.0e8 / (local_cfg.fc * 1e9));
        if (plane.V <= 0.0 || lambda <= 0.0) {
            std::cerr << "computeDatasetSquintFromCenter: invalid V/lambda" << std::endl;
            return false;
        }
        double ratio = -fd_unwrapped * lambda / (2.0 * plane.V);
        if (ratio > 1.0) ratio = 1.0;
        if (ratio < -1.0) ratio = -1.0;
        double angle_deg = std::asin(ratio) * 180.0 / M_PI;

        // Use the raw local angle as the reference and take the difference.
        // This is the actual error angle that should be applied globally.
        double bias_deg = wrap180_deg(angle_deg - theta_sq_local);

        std::cout << "[SQUINT] period=" << per
                  << ", theta_sq_local=" << theta_sq_local
                  << ", estimated_angle=" << angle_deg
                  << ", bias_angle=" << bias_deg
                  << std::endl;

        angles_deg.push_back(bias_deg);
    }

    if (angles_deg.empty()) return false;
    double sum = 0.0; for (double a : angles_deg) sum += a;
    out_squint = sum / double(angles_deg.size());
    return true;
}

bool GMTIProcessor::cuda_upload_async(const std::vector<cd> &data1,
                                     const std::vector<cd> &data2,
                                     size_t Na, size_t Nr) {
    size_t total_bytes = Na * Nr * sizeof(cd);
    
    // 异步拷贝数据到 GPU 原始区 (d1, d2)
    cudaMemcpyAsync(gpu_ptrs_.d1, data1.data(), total_bytes, cudaMemcpyHostToDevice, stream_compute_);
    cudaMemcpyAsync(gpu_ptrs_.d2, data2.data(), total_bytes, cudaMemcpyHostToDevice, stream_compute_);
    
    return (cudaGetLastError() == cudaSuccess);
}

bool GMTIProcessor::cuda_download_sync(std::vector<cd> &out1, std::vector<cd> &out2, 
                                      size_t total) {
    out1.resize(total);
    out2.resize(total);

    // 1. 异步回传结果
    cudaMemcpyAsync(out1.data(), gpu_ptrs_.t1, total * sizeof(cd), cudaMemcpyDeviceToHost, stream_compute_);
    cudaMemcpyAsync(out2.data(), gpu_ptrs_.t2, total * sizeof(cd), cudaMemcpyDeviceToHost, stream_compute_);

    // 2. 显式流同步：确保 CPU 下一行代码拿到的 out1/out2 是完整的
    cudaStreamSynchronize(stream_compute_);

    return true;
}

bool GMTIProcessor::debug_compare_range(const Config &cfg, int periodIdx)
{
    Config local_cfg = cfg;
    std::vector<std::complex<double>> data1, data2;
    std::vector<double> utc;
    double theta_sq = 0.0;

    if (!readPulseBlock(local_cfg, periodIdx, data1, data2, utc, theta_sq)) {
        std::cerr << "debug_compare_range: readPulseBlock failed" << std::endl;
        return false;
    }

    if (!local_cfg.isPC) {
        // CPU reference path
        std::vector<std::complex<double>> cpu_out;
        if (!rangeCompressFFT(local_cfg, data1, cpu_out)) {
            std::cerr << "debug_compare_range: CPU rangeCompressFFT failed" << std::endl;
            return false;
        }

        // GPU cuFFT path (host-side convenience wrapper)
        std::vector<std::complex<float>> gpu_in(data1.size());
        for (size_t i = 0; i < data1.size(); ++i) {
            gpu_in[i] = std::complex<float>((float)data1[i].real(), (float)data1[i].imag());
        }
        std::vector<std::complex<float>> gpu_out;
        if (!rangeCompressCUFFT(gpu_in, gpu_out, local_cfg)) {
            std::cerr << "debug_compare_range: GPU rangeCompressCUFFT failed" << std::endl;
            return false;
        }

        const size_t total = std::min(cpu_out.size(), gpu_out.size());
        double max_abs = 0.0;
        double sum_abs = 0.0;
        size_t max_idx = 0;
        std::vector<std::pair<double, size_t>> top;
        top.reserve(10);

        for (size_t i = 0; i < total; ++i) {
            const double cre = cpu_out[i].real();
            const double cim = cpu_out[i].imag();
            const double gre = (double)gpu_out[i].real();
            const double gim = (double)gpu_out[i].imag();
            const double abs_err = std::hypot(cre - gre, cim - gim);
            sum_abs += abs_err;
            if (abs_err > max_abs) {
                max_abs = abs_err;
                max_idx = i;
            }
            if (top.size() < 10) {
                top.emplace_back(abs_err, i);
                std::sort(top.begin(), top.end(), [](const std::pair<double,size_t>& a, const std::pair<double,size_t>& b){ return a.first > b.first; });
            } else if (abs_err > top.back().first) {
                top.back() = std::make_pair(abs_err, i);
                std::sort(top.begin(), top.end(), [](const std::pair<double,size_t>& a, const std::pair<double,size_t>& b){ return a.first > b.first; });
            }
        }

        const double mean_abs = total ? sum_abs / double(total) : 0.0;
        std::cout << "[RANGE-COMPARE] total=" << total
                  << " max_abs=" << max_abs
                  << " mean_abs=" << mean_abs
                  << " max_idx=" << max_idx << std::endl;
        std::cout << "[RANGE-COMPARE] cpu[max]=(" << cpu_out[max_idx].real() << "," << cpu_out[max_idx].imag()
                  << ") gpu[max]=(" << gpu_out[max_idx].real() << "," << gpu_out[max_idx].imag() << ")" << std::endl;
        std::cout << "[RANGE-COMPARE] top differences:" << std::endl;
        for (const auto &e : top) {
            const size_t i = e.second;
            std::cout << "  idx=" << i << " abs=" << e.first
                      << " cpu=(" << cpu_out[i].real() << "," << cpu_out[i].imag() << ")"
                      << " gpu=(" << gpu_out[i].real() << "," << gpu_out[i].imag() << ")" << std::endl;
        }
        return true;
    }

    std::cerr << "debug_compare_range: cfg.isPC=true not supported for this compare" << std::endl;
    return false;
}

// Debug helper: download device buffer `d1` (assumed cuFloatComplex) into host float complex vector
bool GMTIProcessor::debug_download_d1(std::vector<std::complex<float>> &out, size_t total) {
    if (gpu_ptrs_.d1 == nullptr) return false;
    out.resize(total);
    CUDA_CHECK(cudaMemcpyAsync(out.data(), gpu_ptrs_.d1, total * sizeof(cuFloatComplex), cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));
    return true;
}
