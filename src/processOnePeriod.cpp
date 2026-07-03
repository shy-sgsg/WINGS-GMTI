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
#include <cstdlib>
#include <limits>
#include <sstream>
#include <iomanip>
#include <map>
#include <cerrno>
#include <mutex>
#include <sys/stat.h>
//#include <Eigen/Dense>
#include "GMTIProcessor.hpp"
#include "rangeCompress.hpp"
#include "geo/geoProj.hpp"
#include "rotation_xy.hpp"
#include "unwrap_fd.hpp"
#include "dbs/NewProtocolReader.hpp"
#include "trig_lut.hpp"

using cudacd = cuFloatComplex;
using cd = std::complex<float>;

inline int sgn(double x, double eps = 1e-6) { return (x > eps) - (x < -eps); }

namespace {

std::string joinPathLocal(const std::string &a, const std::string &b)
{
    if (a.empty()) return b;
    if (a.back() == '/' || a.back() == '\\') return a + b;
    return a + "/" + b;
}

bool ensureDirLocal(const std::string &dir)
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

std::string parentDirLocal(const std::string &path)
{
    const size_t pos = path.find_last_of("/\\");
    return pos == std::string::npos ? std::string() : path.substr(0, pos);
}

std::string inferSceneTruthPath(const Config &cfg)
{
    if (!cfg.pc_peak_scene_truth.empty()) {
        return cfg.pc_peak_scene_truth;
    }
    const std::string out_parent = parentDirLocal(cfg.result_add);
    return joinPathLocal(joinPathLocal(out_parent, "truth"), "scene_truth.csv");
}

std::vector<std::string> splitCsvLine(const std::string &line)
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

double toDoubleLocal(const std::string &s, double fallback = std::numeric_limits<double>::quiet_NaN())
{
    if (s.empty()) return fallback;
    char *end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return end && *end == '\0' ? v : fallback;
}

int toIntLocal(const std::string &s, int fallback = -1)
{
    const double v = toDoubleLocal(s);
    return std::isfinite(v) ? static_cast<int>(std::lround(v)) : fallback;
}

struct PcPeakTruth {
    bool ok = false;
    double range_m = std::numeric_limits<double>::quiet_NaN();
    double range_sample_float = std::numeric_limits<double>::quiet_NaN();
    int range_sample_int = -1;
    int beam_id = -1;
    double theta_cmd_deg = std::numeric_limits<double>::quiet_NaN();
};

PcPeakTruth loadPcPeakTruth(const Config &cfg)
{
    PcPeakTruth truth;
    const std::string path = inferSceneTruthPath(cfg);
    std::ifstream in(path.c_str());
    if (!in) {
        return truth;
    }
    std::string line;
    if (!std::getline(in, line)) {
        return truth;
    }
    const std::vector<std::string> header = splitCsvLine(line);
    std::map<std::string, size_t> col;
    for (size_t i = 0; i < header.size(); ++i) col[header[i]] = i;
    auto get = [&](const std::vector<std::string> &cells, const std::string &name) -> std::string {
        const auto it = col.find(name);
        return (it == col.end() || it->second >= cells.size()) ? std::string() : cells[it->second];
    };

    std::vector<std::string> best;
    double best_score = -1.0e300;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> cells = splitCsvLine(line);
        const std::string type = get(cells, "type");
        double score = 0.0;
        if (type == "single_point") score += 1.0e9;
        const double amp = toDoubleLocal(get(cells, "amplitude"), 0.0);
        const double rcs = toDoubleLocal(get(cells, "rcs_db"), -999.0);
        score += amp * 1.0e6 + rcs;
        if (score > best_score) {
            best_score = score;
            best = cells;
        }
    }
    if (best.empty()) {
        return truth;
    }
    truth.range_m = toDoubleLocal(get(best, "range_m"));
    if (!std::isfinite(truth.range_m)) {
        truth.range_m = toDoubleLocal(get(best, "initial_range_m"));
    }
    truth.range_sample_float = toDoubleLocal(get(best, "range_sample_float"));
    if (!std::isfinite(truth.range_sample_float) &&
        cfg.has_sample_delay_us && std::isfinite(truth.range_m) && cfg.fs > 0.0) {
        const double tau_abs_sec = 2.0 * truth.range_m / C;
        const double tau_rel_sec = tau_abs_sec - cfg.sample_delay_us * 1.0e-6;
        truth.range_sample_float = tau_rel_sec * cfg.fs;
    }
    truth.range_sample_int = toIntLocal(get(best, "range_sample_int"), -1);
    if (truth.range_sample_int < 0 && std::isfinite(truth.range_sample_float)) {
        truth.range_sample_int = static_cast<int>(std::lround(truth.range_sample_float));
    }
    truth.beam_id = toIntLocal(get(best, "beam_id_0based"), -1);
    if (truth.beam_id < 0) {
        truth.beam_id = toIntLocal(get(best, "beam_id"), -1);
    }
    truth.theta_cmd_deg = toDoubleLocal(get(best, "theta_cmd_deg"));
    if (!std::isfinite(truth.theta_cmd_deg)) {
        truth.theta_cmd_deg = toDoubleLocal(get(best, "initial_azimuth_deg"));
    }
    truth.ok = std::isfinite(truth.range_m) && std::isfinite(truth.range_sample_float);
    return truth;
}

void writePcPeakDebug(const Config &cfg,
                      int beamIdx,
                      const std::vector<std::complex<float>> &data1,
                      const std::vector<std::complex<float>> &data2)
{
    static std::mutex pc_peak_mutex;
    if (!cfg.debug_pc_peak || data1.empty() || data2.empty()) {
        return;
    }
    std::lock_guard<std::mutex> lock(pc_peak_mutex);
    const int Na = effectivePulseNum(cfg);
    const int Nr = cfg.rg_len;
    if (Na <= 0 || Nr <= 0 ||
        data1.size() < static_cast<size_t>(Na) * static_cast<size_t>(Nr) ||
        data2.size() < static_cast<size_t>(Na) * static_cast<size_t>(Nr)) {
        return;
    }
    const PcPeakTruth truth = loadPcPeakTruth(cfg);
    const int expected = truth.ok ? (truth.range_sample_int - cfg.range_crop_start) : -1;
    const int left = expected >= 0 ? std::max(0, expected - 50) : -1;
    const int right = expected >= 0 ? std::min(Nr - 1, expected + 50) : -1;

    std::vector<double> power(static_cast<size_t>(Nr), 0.0);
    for (int c = 0; c < Nr; ++c) {
        double acc = 0.0;
        for (int r = 0; r < Na; ++r) {
            const size_t off = static_cast<size_t>(r) * static_cast<size_t>(Nr) + static_cast<size_t>(c);
            acc += 0.5 * (std::norm(data1[off]) + std::norm(data2[off]));
        }
        power[static_cast<size_t>(c)] = acc / static_cast<double>(Na);
    }

    auto maxInRange = [&](int l, int r) {
        int idx = -1;
        double val = -1.0;
        if (l <= r && l >= 0 && r < Nr) {
            for (int c = l; c <= r; ++c) {
                if (power[static_cast<size_t>(c)] > val) {
                    val = power[static_cast<size_t>(c)];
                    idx = c;
                }
            }
        }
        return std::pair<int, double>(idx, val);
    };
    const auto local_peak = maxInRange(left, right);
    const auto global_peak = maxInRange(0, Nr - 1);
    const bool valid = truth.ok && expected >= 0 && expected < Nr &&
                       local_peak.first >= 0 && local_peak.second > 0.0;

    const std::string debug_dir = joinPathLocal(cfg.result_add, "debug");
    if (!ensureDirLocal(debug_dir)) {
        return;
    }

    {
        char name[128];
        std::snprintf(name, sizeof(name), "pc_range_profile_beam%02d.csv", beamIdx);
        std::ofstream prof(joinPathLocal(debug_dir, name).c_str());
        prof << "range_bin,range_m,power,in_expected_window\n";
        for (int c = 0; c < Nr; ++c) {
            const double range_m = (c >= 0 && c < static_cast<int>(cfg.Rg.size()))
                ? cfg.Rg[static_cast<size_t>(c)]
                : cfg.R_min + static_cast<double>(c) * cfg.R_bin;
            prof << c << "," << std::setprecision(15) << range_m << ","
                 << power[static_cast<size_t>(c)] << ","
                 << ((left >= 0 && c >= left && c <= right) ? 1 : 0) << "\n";
        }
    }

    const std::string check_path = joinPathLocal(debug_dir, "pc_peak_check.csv");
    const bool need_header = !static_cast<bool>(std::ifstream(check_path.c_str()));
    std::ofstream os(check_path.c_str(), std::ios::app);
    if (!os) return;
    if (need_header) {
        os << "case_id,run_id,period_id,beam_id,"
              "truth_range_m,truth_range_sample_float,truth_range_sample_int,"
              "pc_crop_start,expected_pc_bin_without_offset,"
              "search_win_left,search_win_right,"
              "pc_peak_bin,pc_peak_range_m,pc_peak_power,"
              "pc_peak_offset_bin,pc_peak_range_error_m,"
              "global_peak_bin,global_peak_range_m,global_peak_power,"
              "valid\n";
    }
    auto rangeAt = [&](int bin) -> double {
        if (bin >= 0 && bin < static_cast<int>(cfg.Rg.size())) return cfg.Rg[static_cast<size_t>(bin)];
        return std::numeric_limits<double>::quiet_NaN();
    };
    const double pc_range = rangeAt(local_peak.first);
    const double global_range = rangeAt(global_peak.first);
    os << gmti::runtime::caseId() << ","
       << gmti::runtime::runId() << ","
       << 0 << ","
       << beamIdx << ","
       << std::setprecision(15)
       << truth.range_m << ","
       << truth.range_sample_float << ","
       << truth.range_sample_int << ","
       << cfg.range_crop_start << ","
       << expected << ","
       << left << ","
       << right << ","
       << local_peak.first << ","
       << pc_range << ","
       << local_peak.second << ","
       << (local_peak.first >= 0 && expected >= 0 ? local_peak.first - expected : 0) << ","
       << (std::isfinite(pc_range) && std::isfinite(truth.range_m) ? pc_range - truth.range_m : std::numeric_limits<double>::quiet_NaN()) << ","
       << global_peak.first << ","
       << global_range << ","
       << global_peak.second << ","
       << (valid ? 1 : 0) << "\n";
}

bool pcProfileEnvEnabled()
{
    const char *v = std::getenv("GMTI_DEBUG_PC_PROFILE");
    return v && std::string(v) == "1";
}

bool pcProfile2dEnvEnabled()
{
    const char *v = std::getenv("GMTI_PC_DUMP_2D");
    return v && std::string(v) == "1";
}

bool pcProfileBeamEnabled(int beamIdx)
{
    const char *v = std::getenv("GMTI_PC_DUMP_BEAMS");
    if (!v || !*v) return true;
    std::stringstream ss(v);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (toIntLocal(item, -1) == beamIdx) return true;
    }
    return false;
}

void writeManualPcProfile(const Config &cfg,
                          int beamIdx,
                          const std::vector<std::complex<float>> &data1,
                          const std::vector<std::complex<float>> &data2)
{
    if (!pcProfileEnvEnabled() || !pcProfileBeamEnabled(beamIdx) || data1.empty()) return;
    const int Na = effectivePulseNum(cfg);
    const int Nr = cfg.rg_len;
    const bool has_ch2 = data2.size() >= static_cast<size_t>(Na) * static_cast<size_t>(Nr);
    if (Na <= 0 || Nr <= 0 ||
        data1.size() < static_cast<size_t>(Na) * static_cast<size_t>(Nr)) {
        return;
    }
    const std::string debug_dir = joinPathLocal(cfg.result_add, "debug");
    if (!ensureDirLocal(debug_dir)) return;

    auto rangeAt = [&](int bin) {
        if (bin >= 0 && bin < static_cast<int>(cfg.Rg.size())) return cfg.Rg[static_cast<size_t>(bin)];
        return cfg.R_min + static_cast<double>(bin) * cfg.R_bin;
    };

    std::vector<double> power(static_cast<size_t>(Nr), 0.0);
    int peak_bin = 0;
    double peak_power = -1.0;
    for (int c = 0; c < Nr; ++c) {
        double acc = 0.0;
        for (int p = 0; p < Na; ++p) {
            const size_t off = static_cast<size_t>(p) * static_cast<size_t>(Nr) + static_cast<size_t>(c);
            const double p1 = std::norm(data1[off]);
            if (has_ch2) {
                acc += 0.5 * (p1 + std::norm(data2[off]));
            } else {
                acc += p1;
            }
        }
        power[static_cast<size_t>(c)] = acc / static_cast<double>(Na);
        if (power[static_cast<size_t>(c)] > peak_power) {
            peak_power = power[static_cast<size_t>(c)];
            peak_bin = c;
        }
    }

    char name[128];
    std::snprintf(name, sizeof(name), "pc_range_profile_beam%02d.csv", beamIdx);
    const std::string profile_path = joinPathLocal(debug_dir, name);
    std::ofstream prof(profile_path.c_str());
    prof << "bin,range_m,power,power_db\n";
    for (int c = 0; c < Nr; ++c) {
        const double pwr = power[static_cast<size_t>(c)];
        prof << c << "," << std::setprecision(15) << rangeAt(c) << ","
             << pwr << "," << 10.0 * std::log10(pwr + 1.0e-30) << "\n";
    }

    if (pcProfile2dEnvEnabled()) {
        std::snprintf(name, sizeof(name), "pc_amplitude_2d_beam%02d.csv", beamIdx);
        std::ofstream out2d(joinPathLocal(debug_dir, name).c_str());
        out2d << "pulse_idx,bin,range_m,amp_ch1,amp_ch2,amp_mean,power_mean,power_db\n";
        for (int p = 0; p < Na; ++p) {
            for (int c = 0; c < Nr; ++c) {
                const size_t off = static_cast<size_t>(p) * static_cast<size_t>(Nr) + static_cast<size_t>(c);
                const double a1 = std::abs(data1[off]);
                if (has_ch2) {
                    const double a2 = std::abs(data2[off]);
                    const double pm = 0.5 * (a1 * a1 + a2 * a2);
                    out2d << p << "," << c << "," << std::setprecision(15) << rangeAt(c) << ","
                          << a1 << "," << a2 << "," << 0.5 * (a1 + a2) << ","
                          << pm << "," << 10.0 * std::log10(pm + 1.0e-30) << "\n";
                } else {
                    const double pm = a1 * a1;
                    out2d << p << "," << c << "," << std::setprecision(15) << rangeAt(c) << ","
                          << a1 << ",," << a1 << "," << pm << ","
                          << 10.0 * std::log10(pm + 1.0e-30) << "\n";
                }
            }
        }
    }

    const double peak_db = 10.0 * std::log10(peak_power + 1.0e-30);
    std::cout << "[PC_PROFILE] wrote " << profile_path
              << ", beam=" << beamIdx
              << ", Na=" << Na
              << ", Nr=" << Nr
              << ", peak_bin=" << peak_bin
              << ", peak_range_m=" << rangeAt(peak_bin)
              << ", peak_power_db=" << peak_db << std::endl;
}

} // namespace
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

static inline double beam_center_relative_dir_deg(int squint_side, double theta_deg)
{
    const double side_dir = (squint_side == 1) ? -90.0 : 90.0;
    return wrap180_deg(side_dir - theta_deg);
}

static inline double location_beam_gate_deg(const Config &cfg)
{
    if (cfg.loc_beam_gate_deg > 0.0) {
        return cfg.loc_beam_gate_deg;
    }
    return std::max(0.0, cfg.beamwidth_deg * 0.5) + 0.5;
}

} // namespace

double GMTIProcessor::estimateSquintAngleDeg(const GMTIOutput::Plane &plane, const Config &cfg, double fd_ctr)
{
    const double lambda = (cfg.lambda > 0.0) ? cfg.lambda : (C / cfg.fc);
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

        const double ratio = -fd_ctr * lambda / (2.0 * v);
        const double clipped = std::max(-1.0, std::min(1.0, ratio));
        const double squint_deg = gmti::trig_lut::asin(clipped) * 180.0 / M_PI;

        // std::cout << "[SQUINT] fd_ctr=" << fd_ctr
        //           << ", v=" << v
        //           << ", lambda=" << lambda
        //           << ", ratio=" << ratio
        //           << ", squint=" << squint_deg
        //           << std::endl;
        return squint_deg;
    }

    double refE = 0.0;
    double refN = 0.0;
    Gaussp3(cfg.roi_ll_deg[0], cfg.roi_ll_deg[1], cfg.L0, refE, refN);
    const double dE = refE - plane.E;
    const double dN = refN - plane.N;
    const double bearing_deg = gmti::trig_lut::atan2(dN, dE) * 180.0 / M_PI;
    const double squint_deg = wrap180_deg(bearing_deg - plane.V_angle + 90.0);
    
    // std::cout << "[SQUINT] fd_ctr=nan"
    //           << ", v=" << v
    //           << ", lambda=" << lambda
    //           << ", refE=" << refE
    //           << ", refN=" << refN
    //           << ", planeE=" << plane.E
    //           << ", planeN=" << plane.N
    //           << ", bearing=" << bearing_deg
    //           << ", V_angle=" << plane.V_angle
    //           << ", squint=" << squint_deg
    //           << std::endl;
    
    return squint_deg;
}

double GMTIProcessor::estimateSquintAngleDeg(const std::vector<std::complex<float>> &data,
                                             const GMTIOutput::Plane &plane,
                                             const Config &cfg)
{
    double fd_ctr = 0.0;
    int start_pulse = 0;
    int window_pulses = 0;

    if (estimateCenterFdCtrFromData(data, cfg, fd_ctr, start_pulse, window_pulses))
    {
        // std::cout << "[SQUINT] center-window start=" << start_pulse
        //           << ", count=" << window_pulses
        //           << ", fd_ctr=" << fd_ctr
        //           << std::endl;
        return estimateSquintAngleDeg(plane, cfg, fd_ctr);
    }

    // std::cout << "[SQUINT] center-window fd_ctr estimate failed, fallback to geometric estimate" << std::endl;
    return estimateSquintAngleDeg(plane, cfg);
}

bool GMTIProcessor::pulseCompressionGpuResident(const std::vector<std::complex<float>> &data1,
                                                const std::vector<std::complex<float>> &data2,
                                                const Config &cfg)
{
    const int W = effectivePulseNum(cfg);
    if (W <= 0 || cfg.pulse_len <= 0 || cfg.rg_len <= 0) {
        return false;
    }
    if (data1.size() != static_cast<size_t>(W) * static_cast<size_t>(cfg.pulse_len) ||
        data2.size() != static_cast<size_t>(W) * static_cast<size_t>(cfg.pulse_len) ||
        gpu_ptrs_.d1 == nullptr || gpu_ptrs_.d2 == nullptr) {
        return false;
    }

    const size_t total_bytes = data1.size() * sizeof(std::complex<float>);
    if (cudaMemcpyAsync(gpu_ptrs_.d1, data1.data(), total_bytes, cudaMemcpyHostToDevice, stream_compute_) != cudaSuccess ||
        cudaMemcpyAsync(gpu_ptrs_.d2, data2.data(), total_bytes, cudaMemcpyHostToDevice, stream_compute_) != cudaSuccess) {
        return false;
    }

    return rangeCompressCUFFT_device_pair(cfg.pulse_len, cfg.rg_len, cfg);
}

// 辅助函数：进行脉压处理（支持GPU加速 + CPU回退）
bool GMTIProcessor::pulseCompression(std::vector<std::complex<float>> &data1,
                                     std::vector<std::complex<float>> &data2,
                                     const Config &cfg)
{
    // TIMING_SCOPE(pulseCompression);
    // 对两个通道分别进行脉压
    const int W = effectivePulseNum(cfg);
    if (data1.size() != (size_t)W * cfg.pulse_len ||
        data2.size() != (size_t)W * cfg.pulse_len)
    {
        std::cerr << "输入数据尺寸不匹配脉冲数和脉冲长度。" << std::endl;
        return false;
    }

    // GPU 脉压使用 pulseCompressionGpuResident()，成功后不下载整幅矩阵。
    // 本函数只作为 CPU 回退路径使用。
    auto compress_one = [&](std::vector<std::complex<float>>& data) -> bool {
        std::vector<std::complex<double>> in_d(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            in_d[i] = std::complex<double>(data[i].real(), data[i].imag());
        }
        std::vector<std::complex<double>> out_d;
        if (!rangeCompressFFT(cfg, in_d, out_d)) {
            return false;
        }
        data.resize(out_d.size());
        for (size_t i = 0; i < out_d.size(); ++i) {
            data[i] = std::complex<float>(static_cast<float>(out_d[i].real()),
                                          static_cast<float>(out_d[i].imag()));
        }
        return true;
    };

    if (!compress_one(data1))
        return false;
    if (!compress_one(data2))
        return false;

    if (data1.size() != (size_t)W * cfg.rg_len ||
        data2.size() != (size_t)W * cfg.rg_len)
    {
        std::cerr << "输出数据尺寸不匹配脉冲数和距离长度。" << std::endl;
        return false;
    }

    DBG("[CPU] 脉压处理成功 (CPU FFTW)");
    return true;
}

// 辅助函数：计算多普勒频率轴
inline bool computeDoppler(const std::vector<std::complex<float>> &data, int k, const Config &cfg,
                           std::vector<double> &faAxis, double &fa2, const double &v, const double &theta_deg)
{
    // 获取数据的尺寸
    size_t Na = effectivePulseNum(cfg); // 方位向点数
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
    faAxis.resize(Na);
    const double df = (Na > 0) ? (cfg.PRF / static_cast<double>(Na)) : 0.0;
    for (size_t i = 0; i < Na; ++i)
    {
        faAxis[i] = -0.5 * cfg.PRF + static_cast<double>(i) * df;
    }

    // fa2 = unwrap_prf_to_model(fa2 , cfg.PRF, theta_deg, v, cfg.fc); // 解除模糊

    // 将 fa2 加到 faAxis 上
    for (size_t i = 0; i < faAxis.size(); ++i)
    {
        faAxis[i] += fa2;
    }

    return true; // 返回 true 表示成功
}


bool estimateCenterFdCtrFromData(const std::vector<std::complex<float>> &data,
                                        const Config &cfg,
                                        double &fd_ctr,
                                        int &start_pulse,
                                        int &window_pulses)
{
    fd_ctr = 0.0;
    start_pulse = 0;
    window_pulses = 0;

    const int W = effectivePulseNum(cfg);
    if (W < 2 || cfg.rg_len <= 0) {
        return false;
    }

    const size_t nr = static_cast<size_t>(cfg.rg_len);
    if (data.size() < static_cast<size_t>(W) * nr) {
        return false;
    }

    std::vector<double> faAxis;
    double fa_tmp = 0.0;
    if (!computeDoppler(data, 1, cfg, faAxis, fa_tmp, 0.0, 0.0)) {
        return false;
    }

    fd_ctr = fa_tmp;
    start_pulse = 0;
    window_pulses = W;
    return true;
}

// 处理一个周期的数据
bool GMTIProcessor::processOnePeriod(int periodIdx, const Config &cfg_, const std::vector<std::vector<double>> &posRaw, GMTIOutput &out)
{
    TIMING_SCOPE(processOnePeriod);
    Config cfg = cfg_; // 复制配置，局部修改
    const std::string timing_beam_extra = "beam_id=" + std::to_string(periodIdx);
    // 0) 读取脉冲块（适配双文件/交织）
    std::vector<std::complex<float>> data1, data2;
    std::vector<double> utc;
    std::vector<std::uint8_t> headers;
    std::vector<std::vector<double>> echoPosRaw;

    // 读取脉冲块
    double theta_sq = 0.0; // 方位角平方
    bool readSuccess = false;
    {
        gmti::runtime::TimingScope timing_data_read("data_read", 0, timing_beam_extra);
        if (cfg.INFO_Type) {
            readSuccess = readPulseBlockNewProtocol(cfg, periodIdx, data1, data2, utc, theta_sq, echoPosRaw);
        } else {
            readSuccess = readPulseBlock(cfg, periodIdx, data1, data2, utc, theta_sq);
        }
    }
    if (!readSuccess)
    {
        std::cerr << "读取脉冲块失败。" << std::endl;
        return false;
    }
    if (cfg.INFO_Type) {
        cfg.process_pulse_num = static_cast<int>(utc.size());
    }
    DBG("读取脉冲块成功，脉冲数: " << utc.size());
    
    bool data_on_gpu = false;

    {
        gmti::runtime::TimingScope timing_pulse_compression("pulse_compression", 0, timing_beam_extra);
        // 如果没有脉压，则需要脉压
        if (!cfg.isPC)
        {
            // P1.5 debug needs a stable host-side matrix immediately after pulse compression.
            // Keep the normal GPU-resident path unchanged when debug_pc_peak is disabled.
            data_on_gpu = pulseCompressionGpuResident(data1, data2, cfg);
            bool pc_ok = data_on_gpu;
            if (!pc_ok) {
                pc_ok = pulseCompression(data1, data2, cfg);
            }
            if (!pc_ok) {
                std::cerr << "脉压处理失败。" << std::endl;
                return false;
            }
            DBG(data_on_gpu ? "脉压处理成功，结果驻留GPU" : "脉压处理成功，结果位于CPU内存");
        }
        else
        {
            // 只进行抽取
            // 归一化 + 裁剪到 M，输出到 rc_out（W×M 行主）
            std::vector<std::complex<float>> rc_out;
            const int Lraw = cfg.pulse_len;
            const int W = effectivePulseNum(cfg);
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
    }

    // 脉压完成后再更新采样率，避免影响 rangeCompressFFT / cuFFT 的匹配滤波器构造
    if (!usesRangeCropWindow(cfg)) {
        cfg.fs = cfg.fs * cfg.rg_len / cfg.pulse_len;
    }
    cfg.pulse_len = cfg.rg_len; // 更新脉冲长度为距离采样长度
    // 通用方位向抽取：cfg.pulse_dec 合 1 抽取
    if (!data_on_gpu &&
        (data1.size() != (size_t)effectivePulseNum(cfg) * cfg.rg_len ||
         data2.size() != (size_t)effectivePulseNum(cfg) * cfg.rg_len))
    {
        ERR("数据尺寸不匹配脉冲数和距离长度，无法继续处理。");
        return false;
    }

    const int W_orig = effectivePulseNum(cfg);
    const int M = cfg.rg_len;
    const int dec = cfg.pulse_dec;
    if (dec <= 0)
    {
        ERR("cfg.pulse_dec 必须大于 0（当前 = " << dec << "）。");
        return false;
    }
    if (W_orig % dec != 0)
    {
        ERR("处理脉冲数不是" << dec << "的倍数，无法进行" << dec << "合 1 抽取（当前 process_pulse_num = " << W_orig << "）。");
        return false;
    }

    const int W_new = W_orig / dec;
    if (data_on_gpu)
    {
        if (!cuda_stage_az_decimate_async(W_orig, M, dec)) {
            ERR("GPU 方位向抽取失败。");
            return false;
        }
        data1.clear();
        data2.clear();
        data1.shrink_to_fit();
        data2.shrink_to_fit();
        if (W_new == W_orig) {
            DBG("cfg.pulse_dec == 1，GPU数据保持原方位长度。");
        } else {
            DBG("GPU 方位向" << dec << " 合 1 抽取完成：处理脉冲数 " << W_orig << " -> " << W_new);
        }
    }
    else if (W_new == W_orig)
    {
        DBG("cfg.pulse_dec == 1，跳过抽取。");
    }
    else
    {
        // 临时缓冲，复用（data1 先处理，后处理 data2）
        std::vector<std::complex<float>> tmp((size_t)W_new * M);

        // 处理 data1
        for (int k_new = 0; k_new < W_new; ++k_new)
        {
            size_t out_base = (size_t)k_new * M;
            size_t in_base = (size_t)(dec * k_new) * M;
            for (int m = 0; m < M; ++m)
            {
                std::complex<float> acc = 0.0f;
                for (int d = 0; d < dec; ++d)
                {
                    acc += data1[in_base + (size_t)d * M + m];
                }
                tmp[out_base + m] = acc;
            }
        }
        data1.swap(tmp);

        // 复用 tmp 处理 data2（重新分配确保大小正确）
        tmp.assign((size_t)W_new * M, std::complex<float>(0.0f, 0.0f));
        for (int k_new = 0; k_new < W_new; ++k_new)
        {
            size_t out_base = (size_t)k_new * M;
            size_t in_base = (size_t)(dec * k_new) * M;
            for (int m = 0; m < M; ++m)
            {
                std::complex<float> acc = 0.0f;
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
        ERR("utc 长度" << utc.size() << " 与处理脉冲数(" << W_orig << ") 不匹配，无法同步抽取。");
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
    cfg.process_pulse_num = W_new;
    cfg.PRF = cfg.PRF / dec; // 更新 PRF
    DBG("方位向" << dec << " 合 1 抽取完成：处理脉冲数 " << W_orig << " -> " << cfg.process_pulse_num << "，距离长度 M = " << M << "。");

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
    if (data_on_gpu) {
        if (!cuda_scale_channel2_async(static_cast<float>(coef),
                                       static_cast<size_t>(effectivePulseNum(cfg)) * static_cast<size_t>(cfg.rg_len))) {
            ERR("GPU 通道2幅度校准失败。");
            return false;
        }
        if (pcProfileEnvEnabled() && pcProfileBeamEnabled(periodIdx)) {
            std::vector<std::complex<float>> pc1, pc2;
            const size_t total = static_cast<size_t>(effectivePulseNum(cfg)) * static_cast<size_t>(cfg.rg_len);
            pc1.resize(total);
            pc2.resize(total);
            cudaMemcpyAsync(pc1.data(), gpu_ptrs_.d1, total * sizeof(cd), cudaMemcpyDeviceToHost, stream_compute_);
            cudaMemcpyAsync(pc2.data(), gpu_ptrs_.d2, total * sizeof(cd), cudaMemcpyDeviceToHost, stream_compute_);
            cudaStreamSynchronize(stream_compute_);
            writeManualPcProfile(cfg, periodIdx, pc1, pc2);
        }
    } else {
        const size_t N = data2.size();
        for (size_t i = 0; i < N; ++i)
            data2[i] *= coef;
        writeManualPcProfile(cfg, periodIdx, data1, data2);
    }

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
    double fa_ctr = 0.0;
    std::vector<double> faAxis;
    double angle_deg = 0.0;
    if (data_on_gpu) {
        if (!cuda_compute_fd_ctr_from_d1(1, cfg, fa_ctr)) {
            ERR("GPU 多普勒中心估计失败。");
            return false;
        }
        angle_deg = GMTIProcessor::estimateSquintAngleDeg(plane, cfg, fa_ctr);
        const size_t axisN = static_cast<size_t>(effectivePulseNum(cfg));
        faAxis.resize(axisN);
        const double df = (axisN > 0) ? (cfg.PRF / static_cast<double>(axisN)) : 0.0;
        for (size_t i = 0; i < axisN; ++i) {
            faAxis[i] = -0.5 * cfg.PRF + static_cast<double>(i) * df + fa_ctr;
        }
    } else {
        angle_deg = GMTIProcessor::estimateSquintAngleDeg(data1, plane, cfg);
    }
    double theta_deg = angle_deg + theta_sq; // 斜视角（GMTI 内部估计）+ 帧内方位修正
    DBG("GMTI estimated angle_deg = " << angle_deg << " deg");

    if (!data_on_gpu) {
        computeDoppler(data1, 1, cfg, faAxis, fa_ctr, plane.V, theta_deg);
    }
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

    int skipInt = 0;
    std::vector<std::complex<double>> F1, F2, F1_r, F2_r;
    const size_t Na = static_cast<size_t>(effectivePulseNum(cfg));
    const size_t Nr = static_cast<size_t>(cfg.rg_len);

    std::vector<float> phi_fit;
    std::vector<double> cpu_phi_fit;
    std::vector<double> phi_diss_phase; // 存储相位值
    std::vector<int>    phi_diss_range; // 存储对应的距离向索引

    std::array<double, 2> p_38;
    std::vector<float> phase_map;
    std::array<double, 2> p_38_csi;
    std::vector<double> ph_trace;
    std::vector<double> fa_cut;
    {
        gmti::runtime::TimingScope timing_channel_cancellation("channel_cancellation", 0, timing_beam_extra);
        if (!data_on_gpu) {
            cuda_upload_async(data1, data2, Na, Nr);
        }

        if (!cuda_stage_align_async(skipInt, Na, Nr)) {
            return false;
        }
        if (!cuda_stage_fft_async(Na, Nr)) {
            return false;
        }
        if (!cuda_stage_dbs_async(static_cast<float>(fa_ctr), static_cast<float>(cfg.PRF), Na, Nr)) {
            return false;
        }

        double fa2 = unwrap_prf_to_model(fa_ctr, cfg.PRF, theta_deg, plane.V, cfg.fc); // 解除模糊
        double fa_shift = fa2 - fa_ctr;
        for (size_t i = 0; i < faAxis.size(); ++i) faAxis[i] += fa_shift;

        // --- 4) 距离向初相校正 ---
        double thre_rg = 0.1 * M_PI;

        if (rg_correct_CUDA(cfg, thre_rg, phi_fit, phi_diss_phase, phi_diss_range)) {
            cuda_apply_rg_correction_async(phi_fit, Na, Nr);
        }
        else {
            std::cout << "相位校正失败。" << std::endl;
            return false;
        }
        DBG("最终对齐完成，使用移位脉冲数: " << skipInt);

        if (!clutter_cancel_38_paper_1_p38_cuda(
            faAxis,
            cfg.az_st, cfg.rg_st, cfg.az_ed, cfg.rg_ed,
            cfg,
            p_38
        )) return false;

        // 计算用于 CFAR 的 phase_map（GPU 生成，避免回传大矩阵）
        if (!cuda_download_phase_map(phase_map, Na * Nr)) return false;

        // 对消
        cuda_stage_align_async(skipInt_theory, Na, Nr);
        cuda_stage_fft_async(Na, Nr);
        cuda_stage_dbs_async(static_cast<float>(fa_ctr), static_cast<float>(cfg.PRF), Na, Nr);
        if (rg_correct_CUDA(cfg, thre_rg, phi_fit, phi_diss_phase, phi_diss_range)) {
            cuda_apply_rg_correction_async(phi_fit, Na, Nr);
        }
        else {
            std::cout << "相位校正失败。" << std::endl;
            return false;
        }
         DBG("最终对齐完成，使用移位脉冲数: " << skipInt_theory);

        // --- 8) 最终 CSI 对消 ---
        if (!clutter_cancel_38_paper_1_cuda(
            faAxis,
            cfg.az_st, cfg.rg_st, cfg.az_ed, cfg.rg_ed,
            cfg,
            F1_r, p_38_csi, ph_trace, fa_cut
        )) return false;
    }

    // --- 10) CFAR 处理 ---
    GMTIOutput::Detect targetSel;
    bool targetDetection = false;
    const int band_st = std::max(0, std::min(az_st, effectivePulseNum(cfg) - 1));
    const int band_ed = std::max(0, std::min(az_ed, effectivePulseNum(cfg) - 1));
    {
        gmti::runtime::TimingScope timing_detection("detection", 0, timing_beam_extra);
        std::vector<int> prow, pcol;
        std::vector<float> mydata;
        int cfar_bnum = 16; // 背景单元数
        int c_num = 4;      // 保护单元数
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
        const bool use_gpu_target = true;
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
    const double cosT = std::abs(gmti::trig_lut::cos(deg2rad(thetaRot)));
    const double sinT = std::abs(gmti::trig_lut::sin(deg2rad(thetaRot)));
    const double VE = plane.V * gmti::trig_lut::cos(deg2rad(plane.V_angle));
    const double VN = plane.V * gmti::trig_lut::sin(deg2rad(plane.V_angle));

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
        const double target_azimuth_deg = normalize_azimuth_deg(gmti::trig_lut::atan2(dN, dE) * 180.0 / M_PI);  // 东为0°, 逆时针为正
        // 相对方向：顺时针为正 => 计算 plane - target，再映射到 [-180,180]
        const double direction = wrap180_deg(plane.V_angle - target_azimuth_deg);
        const double beam_center_dir = beam_center_relative_dir_deg(cfg.squint_side, theta_deg);
        const double beam_half_width = location_beam_gate_deg(cfg);
        const double beam_dir_err = wrap180_deg(direction - beam_center_dir);
        if (std::abs(beam_dir_err) > beam_half_width) {
            DBG("目标 " << i << " 因方向超出波束被剔除: direction=" << direction
                        << " beam_center_dir=" << beam_center_dir
                        << " err=" << beam_dir_err
                        << " gate=" << beam_half_width);
            continue;
        }
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
        GMTIOutput::DetectionCsvRecord rec;
        rec.period_id = 0;
        rec.beam_id = periodIdx;
        rec.range_bin = c;
        rec.row = r;
        rec.col = c;
        rec.range_m = Rg;
        rec.theta_cmd_deg = theta_sq;
        rec.theta_true_deg = theta_deg;
        rec.e = xP;
        rec.n = yP;
        rec.lat = lat;
        rec.lon = lng;
        rec.utc = out.utcMid;
        out.detection_records.push_back(rec);
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

bool GMTIProcessor::processOnePeriodFusionCache(int periodIdx,
                                                const Config &cfg_,
                                                const std::vector<std::vector<double>> &posRaw,
                                                size_t slot,
                                                FusionGroupContext &ctx)
{
    Config cfg = cfg_;
    const std::string timing_beam_extra = "beam_id=" + std::to_string(periodIdx);
    std::vector<std::complex<float>> data1, data2;
    std::vector<double> utc;
    std::vector<std::vector<double>> echoPosRaw;
    double theta_sq = 0.0;

    bool readSuccess = false;
    {
        gmti::runtime::TimingScope timing_data_read("data_read", 0, timing_beam_extra);
        if (cfg.INFO_Type) {
            readSuccess = readPulseBlockNewProtocol(cfg, periodIdx, data1, data2, utc, theta_sq, echoPosRaw);
        } else {
            readSuccess = readPulseBlock(cfg, periodIdx, data1, data2, utc, theta_sq);
        }
    }
    if (!readSuccess) {
        return false;
    }
    if (cfg.INFO_Type) {
        cfg.process_pulse_num = static_cast<int>(utc.size());
    }

    bool data_on_gpu = false;
    {
        gmti::runtime::TimingScope timing_pulse_compression("pulse_compression", 0, timing_beam_extra);
        if (!cfg.isPC) {
            // P1.5 debug needs a stable host-side matrix immediately after pulse compression.
            // Keep the normal GPU-resident path unchanged when debug_pc_peak is disabled.
            data_on_gpu = pulseCompressionGpuResident(data1, data2, cfg);
            bool pc_ok = data_on_gpu;
            if (!pc_ok) {
                pc_ok = pulseCompression(data1, data2, cfg);
            }
            if (!pc_ok) {
                return false;
            }
        } else {
            const int Lraw = cfg.pulse_len;
            const int W = effectivePulseNum(cfg);
            const int M1 = cfg.rg_len;
            const int Lraw2M = Lraw / M1;
            std::vector<std::complex<float>> rc_out((size_t)W * M1);
            for (int k = 0; k < W; ++k) {
                for (int m = 0; m < M1; ++m) {
                    rc_out[(size_t)k * M1 + m] = data1[(size_t)k * Lraw + m * Lraw2M];
                }
            }
            data1.swap(rc_out);

            rc_out.assign((size_t)W * M1, std::complex<float>(0.0f, 0.0f));
            for (int k = 0; k < W; ++k) {
                for (int m = 0; m < M1; ++m) {
                    rc_out[(size_t)k * M1 + m] = data2[(size_t)k * Lraw + m * Lraw2M];
                }
            }
            data2.swap(rc_out);
        }
    }

    if (!usesRangeCropWindow(cfg)) {
        cfg.fs = cfg.fs * cfg.rg_len / cfg.pulse_len;
    }
    cfg.pulse_len = cfg.rg_len;

    const int W_orig = effectivePulseNum(cfg);
    const int M = cfg.rg_len;
    const int dec = cfg.pulse_dec;
    if (dec <= 0 || W_orig % dec != 0) {
        return false;
    }
    const int W_new = W_orig / dec;
    if (data_on_gpu) {
        if (!cuda_stage_az_decimate_async(W_orig, M, dec)) {
            return false;
        }
        data1.clear();
        data2.clear();
        data1.shrink_to_fit();
        data2.shrink_to_fit();
    } else if (W_new != W_orig) {
        std::vector<std::complex<float>> tmp((size_t)W_new * M);
        for (int k_new = 0; k_new < W_new; ++k_new) {
            const size_t out_base = (size_t)k_new * M;
            const size_t in_base = (size_t)(dec * k_new) * M;
            for (int m = 0; m < M; ++m) {
                std::complex<float> acc = 0.0f;
                for (int d = 0; d < dec; ++d) {
                    acc += data1[in_base + (size_t)d * M + m];
                }
                tmp[out_base + m] = acc;
            }
        }
        data1.swap(tmp);

        tmp.assign((size_t)W_new * M, std::complex<float>(0.0f, 0.0f));
        for (int k_new = 0; k_new < W_new; ++k_new) {
            const size_t out_base = (size_t)k_new * M;
            const size_t in_base = (size_t)(dec * k_new) * M;
            for (int m = 0; m < M; ++m) {
                std::complex<float> acc = 0.0f;
                for (int d = 0; d < dec; ++d) {
                    acc += data2[in_base + (size_t)d * M + m];
                }
                tmp[out_base + m] = acc;
            }
        }
        data2.swap(tmp);
    }

    if (W_new != W_orig) {
        std::vector<double> utc_tmp(W_new);
        for (int k_new = 0; k_new < W_new; ++k_new) {
            double acc = 0.0;
            for (int d = 0; d < dec; ++d) {
                acc += utc[dec * k_new + d];
            }
            utc_tmp[k_new] = acc / static_cast<double>(dec);
        }
        utc.swap(utc_tmp);
        cfg.process_pulse_num = W_new;
        cfg.PRF = cfg.PRF / dec;
    }

    double utc_mean = 0.0;
    for (double t : utc) {
        utc_mean += t;
    }
    if (!utc.empty()) {
        utc_mean /= static_cast<double>(utc.size());
    }
    for (double &t : utc) {
        if (utc_mean - t > 0.6) {
            t += 1.0;
        } else if (utc_mean - t < -0.6) {
            t -= 1.0;
        }
    }

    const double coef = cfg.calib_coef;
    if (data_on_gpu) {
        if (!cuda_scale_channel2_async(static_cast<float>(coef),
                                       static_cast<size_t>(effectivePulseNum(cfg)) * static_cast<size_t>(cfg.rg_len))) {
            return false;
        }
        if (pcProfileEnvEnabled() && pcProfileBeamEnabled(periodIdx)) {
            std::vector<std::complex<float>> pc1, pc2;
            const size_t total = static_cast<size_t>(effectivePulseNum(cfg)) * static_cast<size_t>(cfg.rg_len);
            pc1.resize(total);
            pc2.resize(total);
            cudaMemcpyAsync(pc1.data(), gpu_ptrs_.d1, total * sizeof(cd), cudaMemcpyDeviceToHost, stream_compute_);
            cudaMemcpyAsync(pc2.data(), gpu_ptrs_.d2, total * sizeof(cd), cudaMemcpyDeviceToHost, stream_compute_);
            cudaStreamSynchronize(stream_compute_);
            writeManualPcProfile(cfg, periodIdx, pc1, pc2);
        }
    } else {
        for (auto &v : data2) {
            v *= coef;
        }
        writeManualPcProfile(cfg, periodIdx, data1, data2);
    }

    GMTIOutput::Plane plane;
    if (cfg.Loc) {
        const std::vector<std::vector<double>> &planePosSource = cfg.INFO_Type ? echoPosRaw : posRaw;
        if (!extractPlanePos(utc, planePosSource, cfg, plane)) {
            return false;
        }
    } else {
        plane.V = 40.0 / 3.6;
        plane.E = plane.N = 0.0;
        plane.H = cfg.MT_nowz;
        plane.V_angle = 0.0;
    }

    std::vector<double> faAxis;
    double fa_ctr = 0.0;
    double angle_deg = 0.0;
    if (data_on_gpu) {
        if (!cuda_compute_fd_ctr_from_d1(1, cfg, fa_ctr)) {
            return false;
        }
        angle_deg = GMTIProcessor::estimateSquintAngleDeg(plane, cfg, fa_ctr);
        const size_t axisN = static_cast<size_t>(effectivePulseNum(cfg));
        faAxis.resize(axisN);
        const double df = (axisN > 0) ? (cfg.PRF / static_cast<double>(axisN)) : 0.0;
        for (size_t i = 0; i < axisN; ++i) {
            faAxis[i] = -0.5 * cfg.PRF + static_cast<double>(i) * df + fa_ctr;
        }
    } else {
        angle_deg = GMTIProcessor::estimateSquintAngleDeg(data1, plane, cfg);
    }
    const double theta_deg = angle_deg + theta_sq;
    if (!data_on_gpu && !computeDoppler(data1, 1, cfg, faAxis, fa_ctr, plane.V, theta_deg)) {
        return false;
    }
    // std::cout << "[fusion][fd] beam=" << periodIdx
    //           << " theta_sq=" << theta_sq
    //           << " angle_from_fd=" << angle_deg
    //           << " theta_for_support=" << theta_deg
    //           << " fd_wrapped=" << fa_ctr
    //           << " PRF=" << cfg.PRF << std::endl;

    int az_center = 0, az_st = 0, az_ed = 0, rg_st = 0, rg_ed = 0;
    double fd_st = 0.0, fd_ed = 0.0, BW_az = 0.0;
    if (!computeDynamicSupportDomain(faAxis, fa_ctr, plane, cfg,
                                     az_center, az_st, az_ed,
                                     fd_st, fd_ed, BW_az,
                                     rg_st, rg_ed)) {
        return false;
    }
    cfg.az_center = az_center;
    cfg.az_st = az_st;
    cfg.az_ed = az_ed;
    cfg.rg_st = rg_st;
    cfg.rg_ed = rg_ed;

    const double baseline = cfg.d_channel;
    const double delta_t = (plane.V > 0.0) ? (baseline / plane.V) : 0.0;
    const int skipInt_theory = static_cast<int>(std::round(delta_t * cfg.PRF));
    const int skipInt = 0;

    const size_t Na = static_cast<size_t>(effectivePulseNum(cfg));
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    FusionBeamMeta beamMeta;
    std::array<double, 2> p_38;
    std::vector<float> phase_map;
    {
        gmti::runtime::TimingScope timing_channel_cancellation("channel_cancellation", 0, timing_beam_extra);
        if (!data_on_gpu) {
            if (!cuda_upload_async(data1, data2, Na, Nr)) {
                return false;
            }
        }
        if (!cuda_stage_align_async(skipInt, Na, Nr)) {
            return false;
        }
        if (!cuda_stage_fft_async(Na, Nr)) {
            return false;
        }
        if (!cuda_stage_dbs_async(static_cast<float>(fa_ctr), static_cast<float>(cfg.PRF), Na, Nr)) {
            return false;
        }

        beamMeta.beam_index = periodIdx;
        beamMeta.slot = static_cast<int>(slot);
        beamMeta.theta_sq = theta_sq;
        beamMeta.theta_true = theta_sq;
        beamMeta.fd_ctr_wrapped = fa_ctr;
        beamMeta.fd_ctr_unwrapped = fa_ctr;
        beamMeta.utc_mid = utc.empty() ? 0.0 : utc[utc.size() / 2];
        beamMeta.plane = plane;
        beamMeta.PRF = cfg.PRF;
        beamMeta.fc_hz = cfg.fc;
        beamMeta.lambda = (cfg.lambda > 0.0) ? cfg.lambda : (C / cfg.fc); 

        if (!exportDbsCacheAfterRecenter(cfg, beamMeta, slot, Na, Nr, ctx.rd, ctx.meta)) {
            return false;
        }

        const double thre_rg = 0.1 * M_PI;
        std::vector<float> phi_fit;
        std::vector<double> phi_diss_phase;
        std::vector<int> phi_diss_range;
        if (!rg_correct_CUDA(cfg, thre_rg, phi_fit, phi_diss_phase, phi_diss_range)) {
            return false;
        }
        cuda_apply_rg_correction_async(phi_fit, Na, Nr);

        if (!clutter_cancel_38_paper_1_p38_cuda(
                faAxis,
                cfg.az_st, cfg.rg_st, cfg.az_ed, cfg.rg_ed,
                cfg,
                p_38)) {
            return false;
        }

        if (!cuda_download_phase_map(phase_map, Na * Nr)) {
            return false;
        }

        if (!cuda_stage_align_async(skipInt_theory, Na, Nr)) {
            return false;
        }
        if (!cuda_stage_fft_async(Na, Nr)) {
            return false;
        }
        if (!cuda_stage_dbs_async(static_cast<float>(fa_ctr), static_cast<float>(cfg.PRF), Na, Nr)) {
            return false;
        }
        if (!rg_correct_CUDA(cfg, thre_rg, phi_fit, phi_diss_phase, phi_diss_range)) {
            return false;
        }
        cuda_apply_rg_correction_async(phi_fit, Na, Nr);

        std::array<double, 2> p_38_csi;
        std::vector<double> ph_trace;
        std::vector<double> fa_cut;
        std::vector<std::complex<double>> csi_trace;
        if (!clutter_cancel_38_paper_1_cuda(
                faAxis,
                cfg.az_st, cfg.rg_st, cfg.az_ed, cfg.rg_ed,
                cfg,
                csi_trace, p_38_csi, ph_trace, fa_cut)) {
            return false;
        }
    }

    const int band_st = std::max(0, std::min(az_st, effectivePulseNum(cfg) - 1));
    const int band_ed = std::max(0, std::min(az_ed, effectivePulseNum(cfg) - 1));
    size_t cfar_hits = 0;
    std::vector<int> prow_new, pcol_new;
    std::vector<float> refined_mydata;
    GMTIOutput::Detect targetSel;
    {
        gmti::runtime::TimingScope timing_detection("detection", 0, timing_beam_extra);
        std::vector<float> mydata;
        std::vector<std::complex<float>> CSI_out;
        bool cfar_ok = dpca_cfar2_fast_cuda(CSI_out, band_st, band_ed,
                                            cfg.pf, 4, 16, "GO", cfg, mydata);
        if (!cfar_ok) {
            std::vector<std::complex<float>> F1_f, F2_f;
            if (!cuda_download_sync(F1_f, F2_f, Na * Nr)) {
                return false;
            }
            if (!cuda_download_csi_sync(CSI_out, Na * Nr)) {
                return false;
            }

            std::vector<std::complex<double>> detect_data(F2_f.size());
            for (size_t i = 0; i < F2_f.size(); ++i) {
                detect_data[i] = std::complex<double>(F2_f[i].real(), F2_f[i].imag());
            }
            if (band_st <= band_ed) {
                for (int r = band_st; r <= band_ed; ++r) {
                    const size_t off = static_cast<size_t>(r) * Nr;
                    for (size_t c = 0; c < Nr; ++c) {
                        const auto v = CSI_out[off + c];
                        detect_data[off + c] = std::complex<double>(v.real(), v.imag());
                    }
                }
            }

            std::vector<double> mydata_d;
            std::vector<int> prow_cpu, pcol_cpu;
            cfar_ok = dpca_cfar2_fast(detect_data, cfg.pf, 4, 16, "GO",
                                      cfg, mydata_d, prow_cpu, pcol_cpu);
            if (cfar_ok) {
                mydata.assign(mydata_d.begin(), mydata_d.end());
            }
        }
        if (!cfar_ok) {
            return false;
        }
        cfar_hits = static_cast<size_t>(
            std::count_if(mydata.begin(), mydata.end(), [](float v) { return v > 0.0f; }));

        std::vector<float> phase_std_list;
        bool cluster_ok = cluster_filter_gap_phase_cuda(mydata, phase_map, cfg.min_points, 2, 0.2f,
                                                        cfg, refined_mydata, prow_new, pcol_new,
                                                        phase_std_list);
        if (!cluster_ok) {
            std::vector<double> mydata_d(mydata.begin(), mydata.end());
            std::vector<double> phase_map_d(phase_map.begin(), phase_map.end());
            std::vector<double> refined_mydata_d;
            std::vector<double> phase_std_list_d;
            cluster_ok = cluster_filter_gap_phase(mydata_d, phase_map_d, cfg.min_points, 2, 0.2,
                                                  cfg, refined_mydata_d, prow_new, pcol_new,
                                                  phase_std_list_d);
            if (cluster_ok) {
                refined_mydata.assign(refined_mydata_d.begin(), refined_mydata_d.end());
                phase_std_list.assign(phase_std_list_d.begin(), phase_std_list_d.end());
            }
        }
        if (!cluster_ok) {
            return false;
        }

        bool targetDetection = target_select_cuda(prow_new, pcol_new, cfg, targetSel);
        if (!targetDetection) {
            std::vector<std::complex<float>> F1_f, F2_f;
            if (!cuda_download_sync(F1_f, F2_f, Na * Nr)) {
                return false;
            }
            std::vector<std::complex<double>> F1(F1_f.size()), F2(F2_f.size());
            for (size_t i = 0; i < F1_f.size(); ++i) {
                F1[i] = std::complex<double>(F1_f[i].real(), F1_f[i].imag());
                F2[i] = std::complex<double>(F2_f[i].real(), F2_f[i].imag());
            }
            targetDetection = target_select(F1, F2, prow_new, pcol_new, cfg, targetSel);
        }
        if (!targetDetection) {
            return false;
        }
    }

    // std::cout << "[fusion][detect] beam=" << periodIdx
    //           << " slot=" << slot
    //           << " cfar_hits=" << cfar_hits
    //           << " clusters=" << prow_new.size()
    //           << " selected=" << targetSel.prow.size()
    //           << " min_points=" << cfg.min_points
    //           << " pf=" << cfg.pf << std::endl;
    (void)cfar_hits;

    beamMeta.phase_slope = p_38[0];
    beamMeta.phase_intercept = p_38[1];
    beamMeta.az_center = az_center;

    if (slot < ctx.detections.size()) {
        auto &raw = ctx.detections[slot];
        raw.clear();
        raw.reserve(targetSel.prow.size());
        const double ref_phase = faAxis[static_cast<size_t>(az_center)] * p_38[0] + p_38[1];
        for (size_t i = 0; i < targetSel.prow.size(); ++i) {
            const int r = targetSel.prow[i];
            const int c = targetSel.pcol[i];
            if (r < 0 || c < 0 || r >= static_cast<int>(Na) || c >= static_cast<int>(Nr)) {
                continue;
            }
            const size_t off = static_cast<size_t>(r) * Nr + static_cast<size_t>(c);
            double dphi = static_cast<double>(phase_map[off]);
            double diff = dphi - ref_phase;
            if (diff > M_PI) {
                diff -= 2.0 * M_PI;
            }
            if (diff < -M_PI) {
                diff += 2.0 * M_PI;
            }
            dphi = ref_phase + diff;

            const double af_wrapped = (std::abs(p_38[0]) < 1e-12)
                ? ((r < static_cast<int>(faAxis.size())) ? faAxis[static_cast<size_t>(r)] : fa_ctr)
                : ((dphi - p_38[1]) / p_38[0]);

            DetectionRaw d;
            d.beam_index = periodIdx;
            d.slot = static_cast<int>(slot);
            d.prow = r;
            d.pcol = c;
            d.range_m = cfg.Rg.empty() ? (cfg.R_min + static_cast<double>(c) * cfg.R_bin)
                                       : cfg.Rg[static_cast<size_t>(c)];
            d.af_wrapped = af_wrapped;
            d.phase = dphi;
            d.amplitude = (i < refined_mydata.size()) ? static_cast<double>(refined_mydata[i]) : 0.0;
            d.utc_mid = beamMeta.utc_mid;
            raw.push_back(d);
        }
    }

    if (slot < ctx.beam_meta.size()) {
        ctx.beam_meta[slot] = beamMeta;
    }
    if (slot < ctx.done.size()) {
        ctx.done[slot] = 1;
    }

    return true;
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

    // 6) 速度估计
    // 新协议包头已经携带 vn/ve/vd，优先直接使用，避免 float32 UTC
    // 分辨率或 POS 时间范围钳制造成首末位置相同、反推速度为 0。
    // POS 列定义为 [utc, lat, lon, alt, vn, ve, vd]，而投影坐标
    // Vx/Vy 分别对应 east/north，因此 Vx=ve, Vy=vn, Vz=-vd。
    double Vx = 0.0;
    double Vy = 0.0;
    double Vz = 0.0;
    if (cfg.INFO_Type &&
        std::all_of(POS.begin(), POS.end(),
                    [](const std::vector<double> &row) {
                        return row.size() >= 7 &&
                               std::isfinite(row[4]) &&
                               std::isfinite(row[5]) &&
                               std::isfinite(row[6]);
                    })) {
        double sum_vn = 0.0;
        double sum_ve = 0.0;
        double sum_vd = 0.0;
        for (const auto &row : POS) {
            sum_vn += row[4];
            sum_ve += row[5];
            sum_vd += row[6];
        }
        const double inv_count = 1.0 / static_cast<double>(POS.size());
        Vx = sum_ve * inv_count;
        Vy = sum_vn * inv_count;
        Vz = -sum_vd * inv_count;
    } else {
        // 旧协议保持原算法：首末位置差 / 波位时长。
        const int W = effectivePulseNum(cfg);
        const double scale = (W > 0) ? (cfg.PRF / double(W)) : 0.0;
        Vx = (E.back() - E.front()) * scale;
        Vy = (N.back() - N.front()) * scale;
        Vz = (alt_m.back() - alt_m.front()) * scale;
    }

    plane.V = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
    plane.V_angle = gmti::trig_lut::atan2(Vy, Vx) * 180.0 / M_PI; // 东北平面速度方向角(°)

    DBG("飞机速度分量: Vx=" << Vx << " Vy=" << Vy << " Vz=" << Vz);

    return true;
}

bool GMTIProcessor::extractPlanePVFromEcho(const Config &cfg,
                                           GMTIOutput::Plane &plane)
{
    Config local_cfg = cfg;
    std::vector<std::complex<float>> data1, data2;
    std::vector<double> utc;
    std::vector<std::vector<double>> echoPosRaw;
    double theta_sq = 0.0;

    if (!readPulseBlockNewProtocol(local_cfg, 1, data1, data2, utc, theta_sq, echoPosRaw)) {
        std::cerr << "extractPlanePVFromEcho: failed to read new-protocol echo" << std::endl;
        return false;
    }
    local_cfg.process_pulse_num = static_cast<int>(utc.size());

    if (!extractPlanePos(utc, echoPosRaw, local_cfg, plane)) {
        std::cerr << "extractPlanePVFromEcho: failed to extract plane from embedded echo pose" << std::endl;
        return false;
    }

    return true;
}

// Compute dataset-level squint by reading center period(s) and estimating fd_ctr.
// This replicates DBS semantics: estimate fd_ctr from the central period(s),
// unwrap PRF ambiguity, convert to squint angle and average when two centers exist.
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
        std::vector<std::complex<float>> data1, data2;
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
        if (local_cfg.INFO_Type) {
            local_cfg.process_pulse_num = static_cast<int>(utc.size());
        }
        if (!local_cfg.isPC) {
            if (!pulseCompression(data1, data2, local_cfg)) {
                std::cerr << "computeDatasetSquintFromCenter: pulseCompression failed for period " << per << std::endl;
                return false;
            }
        } else {
            // Extract/normalize path (same as processOnePeriod extraction)
            const int Lraw = local_cfg.pulse_len;
            const int W = effectivePulseNum(local_cfg);
            const int M1 = local_cfg.rg_len;
            const int Lraw2M = Lraw / M1;
            std::vector<std::complex<float>> rc_out((size_t)W * M1);
            for (int k = 0; k < W; ++k) {
                for (int m = 0; m < M1; ++m) {
                    rc_out[(size_t)k * M1 + m] = data1[(size_t)k * Lraw + m * Lraw2M];
                }
            }
            data1.swap(rc_out);
            rc_out.assign((size_t)W * M1, std::complex<float>(0.0f,0.0f));
            for (int k = 0; k < W; ++k) {
                for (int m = 0; m < M1; ++m) {
                    rc_out[(size_t)k * M1 + m] = data2[(size_t)k * Lraw + m * Lraw2M];
                }
            }
            data2.swap(rc_out);
        }

        // After range compression/extraction, update sampling parameters like processOnePeriod does
        if (!usesRangeCropWindow(local_cfg)) {
            local_cfg.fs =
                local_cfg.fs * local_cfg.rg_len / local_cfg.pulse_len;
        }
        local_cfg.pulse_len = local_cfg.rg_len;

        // Apply azimuth decimation (pulse_dec) same as processOnePeriod
        const int dec = local_cfg.pulse_dec;
        if (dec <= 0) {
            std::cerr << "computeDatasetSquintFromCenter: invalid pulse_dec" << std::endl;
            return false;
        }
        const int W_orig = effectivePulseNum(local_cfg);
        if (W_orig % dec != 0) {
            std::cerr << "computeDatasetSquintFromCenter: pulse_num not divisible by pulse_dec" << std::endl;
            return false;
        }
        const int W_new = W_orig / dec;
        if (W_new != W_orig) {
            // decimate data1 in-place
            std::vector<std::complex<float>> tmp((size_t)W_new * local_cfg.rg_len);
            for (int k_new = 0; k_new < W_new; ++k_new) {
                size_t out_base = (size_t)k_new * local_cfg.rg_len;
                size_t in_base = (size_t)(dec * k_new) * local_cfg.rg_len;
                for (int m = 0; m < local_cfg.rg_len; ++m) {
                    std::complex<float> acc = 0.0f;
                    for (int d = 0; d < dec; ++d) acc += data1[in_base + (size_t)d * local_cfg.rg_len + m];
                    tmp[out_base + m] = acc;
                }
            }
            data1.swap(tmp);
            tmp.assign((size_t)W_new * local_cfg.rg_len, std::complex<float>(0.0f,0.0f));
            for (int k_new = 0; k_new < W_new; ++k_new) {
                size_t out_base = (size_t)k_new * local_cfg.rg_len;
                size_t in_base = (size_t)(dec * k_new) * local_cfg.rg_len;
                for (int m = 0; m < local_cfg.rg_len; ++m) {
                    std::complex<float> acc = 0.0f;
                    for (int d = 0; d < dec; ++d) acc += data2[in_base + (size_t)d * local_cfg.rg_len + m];
                    tmp[out_base + m] = acc;
                }
            }
            data2.swap(tmp);
            local_cfg.process_pulse_num = W_new;
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

        // Estimate wrapped fd_ctr from data.
        double fd_ctr = 0.0;
        int start_pulse = 0, window_pulses = 0;
        if (!estimateCenterFdCtrFromData(data1, local_cfg, fd_ctr, start_pulse, window_pulses)) {
            std::cerr << "computeDatasetSquintFromCenter: center fd estimate failed for period " << per << std::endl;
            return false;
        }

        // Unwrap PRF ambiguity using model
        double fd_unwrapped = unwrap_prf_to_model(fd_ctr, local_cfg.PRF, theta_sq_local, plane.V, local_cfg.fc);

        // Convert to squint angle (deg)
        const double lambda = (local_cfg.lambda > 0.0) ? local_cfg.lambda : (C / local_cfg.fc);
        if (plane.V <= 0.0 || lambda <= 0.0) {
            std::cerr << "computeDatasetSquintFromCenter: invalid V/lambda" << std::endl;
            return false;
        }
        double ratio = -fd_unwrapped * lambda / (2.0 * plane.V);
        if (ratio > 1.0) ratio = 1.0;
        if (ratio < -1.0) ratio = -1.0;
        double angle_deg = gmti::trig_lut::asin(ratio) * 180.0 / M_PI;

        // Use the raw local angle as the reference and take the difference.
        // This is the actual error angle that should be applied globally.
        double bias_deg = wrap180_deg(angle_deg - theta_sq_local);

        // std::cout << "[SQUINT] beam=" << per
        //           << ", theta_sq_local=" << theta_sq_local
        //           << ", estimated_angle=" << angle_deg
        //           << ", bias_angle=" << bias_deg
        //           << std::endl;

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

bool GMTIProcessor::exportDbsCacheAfterRecenter(const Config& cfg,
                                                const FusionBeamMeta& beamMeta,
                                                size_t slot,
                                                size_t Na,
                                                size_t Nr,
                                                RDData& rd,
                                                MetaPack& meta)
{
    if (Na == 0 || Nr == 0 || gpu_ptrs_.t1 == nullptr) {
        return false;
    }
    if (slot >= rd.amp.size() || slot >= rd.fd_axis.size() ||
        slot >= rd.rg_axis.size() || slot >= meta.beams.size()) {
        return false;
    }

    const size_t total = Na * Nr;
    std::vector<std::complex<float>> h_az(total);
    cudaMemcpyAsync(h_az.data(), gpu_ptrs_.t1, total * sizeof(cd),
                    cudaMemcpyDeviceToHost, stream_compute_);
    cudaStreamSynchronize(stream_compute_);

    Image2D<float> amp(static_cast<int>(Na), static_cast<int>(Nr));
    for (size_t r = 0; r < Na; ++r) {
        for (size_t c = 0; c < Nr; ++c) {
            amp.at(static_cast<int>(r), static_cast<int>(c)) =
                std::abs(h_az[r * Nr + c]);
        }
    }

    std::vector<float> fd_axis(Na, 0.0f);
    const double df = (Na > 0) ? (cfg.PRF / static_cast<double>(Na)) : 0.0;
    for (size_t r = 0; r < Na; ++r) {
        fd_axis[r] = static_cast<float>(-0.5 * cfg.PRF + static_cast<double>(r) * df);
    }

    std::vector<float> rg_axis(Nr, 0.0f);
    for (size_t c = 0; c < Nr; ++c) {
        if (c < cfg.Rg.size()) {
            rg_axis[c] = static_cast<float>(cfg.Rg[c]);
        } else {
            rg_axis[c] = static_cast<float>(cfg.R_min + static_cast<double>(c) * cfg.R_bin);
        }
    }

    rd.nEff = static_cast<int>(Na);
    rd.amp[slot] = amp;
    rd.fd_axis[slot] = fd_axis;
    rd.rg_axis[slot] = rg_axis;

    MetaPerBeam m;
    m.vN = static_cast<float>(beamMeta.plane.V * gmti::trig_lut::sin(beamMeta.plane.V_angle * M_PI / 180.0));
    m.vE = static_cast<float>(beamMeta.plane.V * gmti::trig_lut::cos(beamMeta.plane.V_angle * M_PI / 180.0));
    m.vU = 0.0f;
    m.x = static_cast<float>(beamMeta.plane.E);
    m.y = static_cast<float>(beamMeta.plane.N);
    m.z = static_cast<float>(beamMeta.plane.H);
    m.fd_ctr = static_cast<float>(beamMeta.fd_ctr_wrapped);
    m.angle_deg = static_cast<float>(beamMeta.theta_sq);
    meta.beams[slot] = m;

    return true;
}

bool GMTIProcessor::debug_compare_range(const Config &cfg, int periodIdx)
{
    Config local_cfg = cfg;
    std::vector<std::complex<float>> data1, data2;
    std::vector<double> utc;
    double theta_sq = 0.0;

    if (!readPulseBlock(local_cfg, periodIdx, data1, data2, utc, theta_sq)) {
        std::cerr << "debug_compare_range: readPulseBlock failed" << std::endl;
        return false;
    }

    if (!local_cfg.isPC) {
        // CPU reference path
        std::vector<std::complex<double>> data1_d(data1.size());
        for (size_t i = 0; i < data1.size(); ++i) {
            data1_d[i] = std::complex<double>(data1[i].real(), data1[i].imag());
        }
        std::vector<std::complex<double>> cpu_out;
        if (!rangeCompressFFT(local_cfg, data1_d, cpu_out)) {
            std::cerr << "debug_compare_range: CPU rangeCompressFFT failed" << std::endl;
            return false;
        }

        // GPU cuFFT path (host-side convenience wrapper)
        std::vector<std::complex<float>> gpu_out;
        if (!rangeCompressCUFFT(data1, gpu_out, local_cfg)) {
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
