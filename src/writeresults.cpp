#include "GMTIProcessor.hpp"
#include "geo/geoProj.hpp" // 包含 Gaussp3RV
#include "runtime_diagnostics.hpp"

#include <vector>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h> // opendir, readdir, closedir
#include <cstdio>   // std::snprintf

static inline bool starts_with(const char *s, const char *pfx)
{
    while (*pfx)
    {
        if (*s++ != *pfx++)
            return false;
    }
    return true;
}
static inline bool ends_with(const char *s, const char *sfx)
{
    const size_t ls = std::strlen(s), lr = std::strlen(sfx);
    return (ls >= lr) && (std::strcmp(s + (ls - lr), sfx) == 0);
}

static inline double wrap180_deg(double angle_deg)
{
    angle_deg = std::fmod(angle_deg + 180.0, 360.0);
    if (angle_deg < 0.0)
        angle_deg += 360.0;
    return angle_deg - 180.0;
}

static std::string csv_escape_local(const std::string &s)
{
    if (s.find_first_of(",\"\n\r") == std::string::npos) {
        return s;
    }
    std::string out = "\"";
    for (char ch : s) {
        out += (ch == '"') ? "\"\"" : std::string(1, ch);
    }
    out += "\"";
    return out;
}

static std::string csv_num(double v)
{
    if (!std::isfinite(v)) {
        return "";
    }
    std::ostringstream os;
    os << std::setprecision(15) << v;
    return os.str();
}

static std::string csv_int_or_blank(int v)
{
    return v >= 0 ? std::to_string(v) : "";
}

static double wrap_pi_csv(double x)
{
    x = std::fmod(x + M_PI, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x - M_PI;
}

static double clamp_unit_local(double v)
{
    if (!std::isfinite(v)) {
        return v;
    }
    return std::max(-1.0, std::min(1.0, v));
}

static double deg_from_rad_local(double v)
{
    return v * 180.0 / M_PI;
}

static double rad_from_deg_local(double v)
{
    return v * M_PI / 180.0;
}

static double normalize_theta_deg_local(double theta_deg)
{
    theta_deg = std::fmod(theta_deg + 180.0, 360.0);
    if (theta_deg < 0.0) {
        theta_deg += 360.0;
    }
    return theta_deg - 180.0;
}

static double theta_from_position_deg_local(const Config &cfg,
                                            double platform_v_angle_deg,
                                            double platform_e,
                                            double platform_n,
                                            double target_e,
                                            double target_n)
{
    if (!std::isfinite(platform_v_angle_deg) ||
        !std::isfinite(platform_e) ||
        !std::isfinite(platform_n) ||
        !std::isfinite(target_e) ||
        !std::isfinite(target_n)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double dE = target_e - platform_e;
    const double dN = target_n - platform_n;
    const double target_azimuth_deg = deg_from_rad_local(std::atan2(dN, dE));
    const double side_dir = (cfg.squint_side == 1) ? -90.0 : 90.0;
    return normalize_theta_deg_local(side_dir - (platform_v_angle_deg - target_azimuth_deg));
}

static double theta_from_af_geometry_deg_local(const Config &cfg,
                                               double platform_v,
                                               double af_geometry_hz,
                                               double lambda)
{
    if (!(platform_v > 0.0) || !(lambda > 0.0) || !std::isfinite(af_geometry_hz)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double sin_a = clamp_unit_local(af_geometry_hz * lambda / (2.0 * platform_v));
    if (!std::isfinite(sin_a)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double theta_deg = deg_from_rad_local(std::asin(sin_a));
    (void)cfg;
    return theta_deg;
}

static double horizontal_range_from_slant(double slant_range_m, double dz_m)
{
    if (!std::isfinite(slant_range_m) || !std::isfinite(dz_m)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double v = slant_range_m * slant_range_m - dz_m * dz_m;
    return (v > 0.0) ? std::sqrt(v) : 0.0;
}

static double slant_range_from_horizontal(double ground_range_m, double dz_m)
{
    if (!std::isfinite(ground_range_m) || !std::isfinite(dz_m)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::sqrt(std::max(0.0, ground_range_m * ground_range_m + dz_m * dz_m));
}

static std::string dirname_local(std::string path)
{
    while (!path.empty() && (path.back() == '/' || path.back() == '\\')) {
        path.pop_back();
    }
    const size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return ".";
    if (pos == 0) return "/";
    return path.substr(0, pos);
}

static bool file_exists_local(const std::string &path)
{
    std::ifstream in(path.c_str());
    return static_cast<bool>(in);
}

static std::vector<std::string> split_csv_simple(const std::string &line)
{
    std::vector<std::string> out;
    std::string cur;
    bool in_quote = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
        if (ch == '"') {
            if (in_quote && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quote = !in_quote;
            }
        } else if (ch == ',' && !in_quote) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    out.push_back(cur);
    if (!out.empty() && !out.back().empty() && out.back().back() == '\r') {
        out.back().pop_back();
    }
    return out;
}

static int csv_col_index(const std::vector<std::string> &header, const std::string &name)
{
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) return static_cast<int>(i);
    }
    return -1;
}

static double csv_field_double(const std::vector<std::string> &row, int idx)
{
    if (idx < 0 || idx >= static_cast<int>(row.size()) || row[static_cast<size_t>(idx)].empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::strtod(row[static_cast<size_t>(idx)].c_str(), nullptr);
}

static int csv_field_int(const std::vector<std::string> &row, int idx, int fallback = -1)
{
    if (idx < 0 || idx >= static_cast<int>(row.size()) || row[static_cast<size_t>(idx)].empty()) {
        return fallback;
    }
    return static_cast<int>(std::strtol(row[static_cast<size_t>(idx)].c_str(), nullptr, 10));
}

struct MovingTruthCsvRow {
    bool valid = false;
    int period_id = -1;
    int beam_id = -1;
    int range_bin = -1;
    int row_truth = -1;
    double phi_total_truth_rad = std::numeric_limits<double>::quiet_NaN();
    double phi_static_truth_rad = std::numeric_limits<double>::quiet_NaN();
    double phi_motion_truth_rad = std::numeric_limits<double>::quiet_NaN();
    double v_truth_mps = std::numeric_limits<double>::quiet_NaN();
    double af_motion_truth_hz = std::numeric_limits<double>::quiet_NaN();
    double af_geometry_truth_hz = std::numeric_limits<double>::quiet_NaN();
    double af_total_truth_hz = std::numeric_limits<double>::quiet_NaN();
    double target_e_truth = std::numeric_limits<double>::quiet_NaN();
    double target_n_truth = std::numeric_limits<double>::quiet_NaN();
};

static double row_truth_float_on_fa_axis(const MovingTruthCsvRow &truth, const Config &cfg)
{
    if (!(cfg.fd_res > 0.0) ||
        !std::isfinite(truth.af_total_truth_hz) ||
        !std::isfinite(truth.af_geometry_truth_hz) ||
        !(cfg.PRF > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double first = truth.af_geometry_truth_hz - 0.5 * cfg.PRF;
    return (truth.af_total_truth_hz - first) / cfg.fd_res;
}

static int row_truth_raw_index(const MovingTruthCsvRow &truth, const Config &cfg)
{
    const double v = row_truth_float_on_fa_axis(truth, cfg);
    return std::isfinite(v) ? static_cast<int>(std::llround(v)) : -1;
}

static int row_truth_after_recenter_index(const MovingTruthCsvRow &truth,
                                         const GMTIOutput::DetectionCsvRecord &r,
                                         const Config &cfg)
{
    if (!(cfg.fd_res > 0.0) ||
        !std::isfinite(truth.af_total_truth_hz) ||
        !std::isfinite(r.fd_ctr_unwrapped) ||
        !(cfg.PRF > 0.0)) {
        return -1;
    }
    const double first = r.fd_ctr_unwrapped - 0.5 * cfg.PRF;
    const double row_float = (truth.af_total_truth_hz - first) / cfg.fd_res;
    return static_cast<int>(std::llround(row_float));
}

static std::string moving_truth_path_from_cfg(const Config &cfg)
{
    const std::string result_dir = cfg.result_add.empty() ? "." : cfg.result_add;
    const std::string root = dirname_local(result_dir);
    const std::string p1 = root + "/truth/moving_target_truth.csv";
    if (file_exists_local(p1)) return p1;
    const std::string p2 = result_dir + "/../truth/moving_target_truth.csv";
    if (file_exists_local(p2)) return p2;
    return p1;
}

static std::vector<MovingTruthCsvRow> load_moving_truth_rows(const Config &cfg)
{
    std::vector<MovingTruthCsvRow> rows;
    std::ifstream in(moving_truth_path_from_cfg(cfg).c_str());
    if (!in) return rows;

    std::string line;
    if (!std::getline(in, line)) return rows;
    const std::vector<std::string> header = split_csv_simple(line);
    const int c_period = csv_col_index(header, "period_id");
    const int c_beam = csv_col_index(header, "beam_id");
    const int c_range = csv_col_index(header, "range_bin");
    const int c_row = csv_col_index(header, "row_truth");
    const int c_phi_total = csv_col_index(header, "phi_total_truth_rad");
    const int c_phi_static = csv_col_index(header, "phi_static_truth_rad");
    const int c_phi_motion = csv_col_index(header, "phi_motion_truth_rad");
    const int c_v = csv_col_index(header, "v_truth_mps");
    const int c_af_motion = csv_col_index(header, "af_motion_truth_hz");
    const int c_af_geometry = csv_col_index(header, "af_geometry_truth_hz");
    const int c_af_total = csv_col_index(header, "af_total_truth_hz");
    const int c_e = csv_col_index(header, "target_e_truth");
    const int c_n = csv_col_index(header, "target_n_truth");

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> f = split_csv_simple(line);
        MovingTruthCsvRow r;
        r.valid = true;
        r.period_id = csv_field_int(f, c_period);
        r.beam_id = csv_field_int(f, c_beam);
        r.range_bin = csv_field_int(f, c_range);
        r.row_truth = csv_field_int(f, c_row);
        r.phi_total_truth_rad = csv_field_double(f, c_phi_total);
        r.phi_static_truth_rad = csv_field_double(f, c_phi_static);
        r.phi_motion_truth_rad = csv_field_double(f, c_phi_motion);
        r.v_truth_mps = csv_field_double(f, c_v);
        r.af_motion_truth_hz = csv_field_double(f, c_af_motion);
        r.af_geometry_truth_hz = csv_field_double(f, c_af_geometry);
        r.af_total_truth_hz = csv_field_double(f, c_af_total);
        r.target_e_truth = csv_field_double(f, c_e);
        r.target_n_truth = csv_field_double(f, c_n);
        rows.push_back(r);
    }
    return rows;
}

static const MovingTruthCsvRow *match_moving_truth(
    const GMTIOutput::DetectionCsvRecord &r,
    const std::vector<MovingTruthCsvRow> &truths)
{
    const MovingTruthCsvRow *best = nullptr;
    int best_score = 1 << 30;
    for (size_t i = 0; i < truths.size(); ++i) {
        const MovingTruthCsvRow &t = truths[i];
        if (!t.valid) continue;
        if (r.beam_id >= 0 && t.beam_id >= 0 && r.beam_id != t.beam_id) continue;
        if (r.period_id >= 0 && t.period_id >= 0 && r.period_id != t.period_id) continue;
        const int dr = (r.range_bin >= 0 && t.range_bin >= 0) ? std::abs(r.range_bin - t.range_bin) : 9999;
        if (dr > 2) continue;
        const int da = (r.row >= 0 && t.row_truth >= 0) ? std::abs(r.row - t.row_truth) : 0;
        const int score = dr * 1000 + da;
        if (score < best_score) {
            best_score = score;
            best = &t;
        }
    }
    return best;
}

static std::string result_file_from_cfg(const Config &cfg)
{
    if (cfg.result_file_id <= 0 || cfg.result_file_id > 99) {
        return "";
    }
    std::string dir = cfg.result_add;
    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') {
        dir.push_back('/');
    }
    char name[16];
    std::snprintf(name, sizeof(name), "GMTI%02d.bin", cfg.result_file_id);
    return dir + name;
}

bool GMTIProcessor::nextGMTIFileName(const std::string &dir,
                                     std::string &out_path,
                                     int &out_index) const
{
    out_path.clear();
    out_index = 0;

    // 1) 确保目录存在
    std::string d = dir;
    if (!d.empty() && d.back() != '/' && d.back() != '\\')
        d.push_back('/');
    if (!mkdir_p(d))
    {
        std::cerr << "nextGMTIFileName: 创建目录失败: " << d << "\n";
        return false;
    }

    // 2) 打开目录并扫描现有文件，匹配 "GMTIxx.bin"
    int max_idx = 0;

    DIR *dp = ::opendir(d.c_str());
    if (dp)
    {
        while (dirent *ent = ::readdir(dp))
        {
            const char *name = ent->d_name;
            // 过滤掉 "." ".."
            if (name[0] == '.')
                continue;

            // 格式严格：长度==10，前缀GMTI，后缀.bin，中间两位数字
            // "GMTI" + 2 + ".bin" = 4 + 2 + 4 = 10
            const size_t len = std::strlen(name);
            if (len != 10)
                continue;
            if (!starts_with(name, "GMTI"))
                continue;
            if (!ends_with(name, ".bin"))
                continue;
            if (name[4] < '0' || name[4] > '9')
                continue;
            if (name[5] < '0' || name[5] > '9')
                continue;

            int idx = (name[4] - '0') * 10 + (name[5] - '0');
            if (idx >= 1 && idx <= 99)
            {
                if (idx > max_idx)
                    max_idx = idx;
            }
        }
        ::closedir(dp);
    }
    else
    {
        // 目录不可读也不致命（可能刚创建），从 01 开始
    }

    // 3) 下一个序号
    const int next_idx = max_idx + 1;
    if (next_idx > 99)
    {
        std::cerr << "nextGMTIFileName: 已达到 99 个文件，无法继续编号。\n";
        return false;
    }

    // 4) 组装路径
    char fname[16];
    std::snprintf(fname, sizeof(fname), "GMTI%02d.bin", next_idx);
    out_path = d + fname;
    out_index = next_idx;
    return true;
}

static inline void put_u8(std::vector<uint8_t> &buf, size_t pos, uint8_t v)
{
    buf[pos] = v;
}

static inline void put_u16_le(std::vector<uint8_t> &buf, size_t pos, uint16_t v)
{
    buf[pos + 0] = static_cast<uint8_t>(v & 0xFF);
    buf[pos + 1] = static_cast<uint8_t>((v >> 8) & 0xFF);
}

static inline void put_f32_le(std::vector<uint8_t> &buf, size_t pos, float v)
{
    static_assert(sizeof(float) == 4, "float must be 4 bytes");
    uint32_t u;
    std::memcpy(&u, &v, sizeof(u));
    // 写为 little-endian
    buf[pos + 0] = static_cast<uint8_t>(u & 0xFF);
    buf[pos + 1] = static_cast<uint8_t>((u >> 8) & 0xFF);
    buf[pos + 2] = static_cast<uint8_t>((u >> 16) & 0xFF);
    buf[pos + 3] = static_cast<uint8_t>((u >> 24) & 0xFF);
}

static inline void put_f64_le(std::vector<uint8_t> &buf, size_t pos, double v)
{
    static_assert(sizeof(double) == sizeof(uint64_t), "double must be 8 bytes");
    uint64_t u = 0;
    std::memcpy(&u, &v, sizeof(u));
    for (int i = 0; i < 8; ++i)
    {
        buf[pos + static_cast<size_t>(i)] = static_cast<uint8_t>((u >> (8 * i)) & 0xFF);
    }
}

static inline void put_i32_le(std::vector<uint8_t> &buf, size_t pos, int32_t v);
static inline int32_t quant_deg_to_i32(double deg);

static inline uint16_t quant_speed_to_u16(double speed)
{
    const double q = std::round(std::max(0.0, speed) / 0.01);
    if (q > 65535.0)
        return 65535;
    return static_cast<uint16_t>(q);
}

bool GMTIProcessor::writeResult(const std::vector<double> &res, const Config &cfg)
{
    // 布局：2 + N*36
    // 单记录：id(u16) + lon(i32) + lat(i32) + speed(u16) + direction(f64) + range(f64) + utc(f64)
    constexpr size_t REC_BYTES = 36;
    constexpr size_t HEADER_BYTES = 2;
    const size_t n_targets = res.size() / 8;
    const size_t FILE_BYTES = HEADER_BYTES + n_targets * REC_BYTES;

    std::vector<uint8_t> buf(FILE_BYTES, 0u);

    // [0..1] 写总数（u16, LE）
    put_u16_le(buf, 0, static_cast<uint16_t>(n_targets));

    // 逐目标写 36 字节记录（含每目标 utc）
    for (size_t i = 0; i < n_targets; ++i)
    {
        const size_t base = HEADER_BYTES + i * REC_BYTES;

        const double lat = res[i * 8 + 0];
        const double lng = res[i * 8 + 1];
        const double target_utc = res[i * 8 + 5];
        const double direction = wrap180_deg(res[i * 8 + 6]);
        const double range = res[i * 8 + 7];

        // 0..1: id (u16, LE)
        put_u16_le(buf, base + 0, static_cast<uint16_t>(std::max<size_t>(1, std::min<size_t>(65535, i + 1))));

        // 2..5: 经度（int32）
        put_i32_le(buf, base + 2, quant_deg_to_i32(lng));
        // 6..9: 纬度（int32）
        put_i32_le(buf, base + 6, quant_deg_to_i32(lat));
        // 10..11: 速度（u16），当前周期检测结果暂写 0
        put_u16_le(buf, base + 10, 0);
        // 12..19: 方向（double）
        put_f64_le(buf, base + 12, direction);
        // 20..27: 距离（double）
        put_f64_le(buf, base + 20, range);
        // 28..35: 目标 utc（double）
        put_f64_le(buf, base + 28, target_utc);
    }

    // 目标路径：优先使用当前回波文件编号，避免扫描目录自动加一导致关联窗口错位。
    std::string outpath;
    int idx = 0;
    if (cfg.result_file_id > 0 && cfg.result_file_id <= 99) {
        std::string d = cfg.result_add;
        if (!d.empty() && d.back() != '/' && d.back() != '\\') {
            d.push_back('/');
        }
        if (!mkdir_p(d)) {
            std::cerr << "writeResult: 创建目录失败: " << d << "\n";
            return false;
        }
        char fname[16];
        std::snprintf(fname, sizeof(fname), "GMTI%02d.bin", cfg.result_file_id);
        outpath = d + fname;
        idx = cfg.result_file_id;
    } else if (!nextGMTIFileName(cfg.result_add, outpath, idx)) {
        return false;
    }

    std::ofstream ofs(outpath.c_str(), std::ios::binary);
    if (!ofs)
    {
        std::cerr << "writeResult: 无法打开输出文件: " << outpath << "\n";
        return false;
    }
    ofs.write(reinterpret_cast<const char *>(buf.data()),
              static_cast<std::streamsize>(buf.size()));
    if (!ofs)
    {
        std::cerr << "writeResult: 写文件失败: " << outpath << "\n";
        return false;
    }
    ofs.close();
    std::cout << "[GMTI] writeResult: 写入当前检测结果 "
              << outpath << "，目标数: " << n_targets << std::endl;
    return true;
}

bool GMTIProcessor::writeDetectionCsv(
    const std::vector<GMTIOutput::DetectionCsvRecord> &records,
    const Config &cfg,
    const std::string &source_file) const
{
    if (!mkdir_p(cfg.result_add)) {
        std::cerr << "writeDetectionCsv: 创建目录失败: " << cfg.result_add << "\n";
        return false;
    }

    std::string outpath = cfg.result_add;
    if (!outpath.empty() && outpath.back() != '/' && outpath.back() != '\\') {
        outpath.push_back('/');
    }
    outpath += "detection_results.csv";

    std::ofstream os(outpath.c_str());
    if (!os) {
        std::cerr << "writeDetectionCsv: 无法打开输出文件: " << outpath << "\n";
        return false;
    }

    const std::string case_id = gmti::runtime::caseId();
    const std::string run_id = gmti::runtime::runId();
    std::string result_id = gmti::runtime::resultId();
    if (result_id.empty() && cfg.result_file_id > 0) {
        char name[16];
        std::snprintf(name, sizeof(name), "GMTI%02d", cfg.result_file_id);
        result_id = name;
    }
    const std::string src = source_file.empty() ? result_file_from_cfg(cfg) : source_file;
    const std::vector<MovingTruthCsvRow> moving_truth_rows = load_moving_truth_rows(cfg);
    struct CsvWriteEntry {
        const GMTIOutput::DetectionCsvRecord *rec = nullptr;
        const MovingTruthCsvRow *truth = nullptr;
        double power = std::numeric_limits<double>::quiet_NaN();
        double position_error_m = std::numeric_limits<double>::quiet_NaN();
    };
    std::vector<CsvWriteEntry> write_entries;
    write_entries.reserve(records.size());
    auto same_nms_cluster = [](const GMTIOutput::DetectionCsvRecord &a,
                               const GMTIOutput::DetectionCsvRecord &b) {
        if (a.beam_id != b.beam_id) return false;
        if (a.range_bin >= 0 && b.range_bin >= 0 && std::abs(a.range_bin - b.range_bin) > 1) return false;
        if (a.row >= 0 && b.row >= 0 && std::abs(a.row - b.row) > 1) return false;
        if (std::isfinite(a.radial_velocity_mps) && std::isfinite(b.radial_velocity_mps) &&
            std::abs(a.radial_velocity_mps - b.radial_velocity_mps) > 0.35) return false;
        return true;
    };
    auto prefer_nms_entry = [](const CsvWriteEntry &lhs, const CsvWriteEntry &rhs) {
        const bool lhs_power_ok = std::isfinite(lhs.power) && lhs.power > 0.0;
        const bool rhs_power_ok = std::isfinite(rhs.power) && rhs.power > 0.0;
        if (lhs_power_ok != rhs_power_ok) {
            return rhs_power_ok;
        }
        if (lhs_power_ok && rhs_power_ok) {
            if (std::abs(lhs.power - rhs.power) > 1.0e-12) {
                return rhs.power > lhs.power;
            }
        }
        const bool lhs_pos_ok = std::isfinite(lhs.position_error_m);
        const bool rhs_pos_ok = std::isfinite(rhs.position_error_m);
        if (lhs_pos_ok != rhs_pos_ok) {
            return rhs_pos_ok;
        }
        if (lhs_pos_ok && rhs_pos_ok && std::abs(lhs.position_error_m - rhs.position_error_m) > 1.0e-9) {
            return rhs.position_error_m < lhs.position_error_m;
        }
        if (lhs.truth && rhs.truth) {
            const int lhs_range_gap = std::abs(lhs.rec->range_bin - lhs.truth->range_bin);
            const int rhs_range_gap = std::abs(rhs.rec->range_bin - rhs.truth->range_bin);
            if (lhs_range_gap != rhs_range_gap) {
                return rhs_range_gap < lhs_range_gap;
            }
            const int lhs_row_gap = std::abs(lhs.rec->row - lhs.truth->row_truth);
            const int rhs_row_gap = std::abs(rhs.rec->row - rhs.truth->row_truth);
            if (lhs_row_gap != rhs_row_gap) {
                return rhs_row_gap < lhs_row_gap;
            }
        }
        return false;
    };
    for (const auto &r : records) {
        CsvWriteEntry cur;
        cur.rec = &r;
        cur.truth = match_moving_truth(r, moving_truth_rows);
        cur.power = std::isfinite(r.amplitude) ? (r.amplitude * r.amplitude)
                                               : std::numeric_limits<double>::quiet_NaN();
        if (cur.truth) {
            cur.position_error_m = std::hypot(r.new_e - cur.truth->target_e_truth,
                                              r.new_n - cur.truth->target_n_truth);
        }
        bool merged = false;
        for (auto &kept : write_entries) {
            if (!same_nms_cluster(*kept.rec, r)) {
                continue;
            }
            if (prefer_nms_entry(kept, cur)) {
                kept = cur;
            }
            merged = true;
            break;
        }
        if (!merged) {
            write_entries.push_back(cur);
        }
    }

    os << "case_id,run_id,result_id,period_id,beam_id,det_id,range_bin,range_m,"
          "theta_cmd_deg,theta_true_deg,row,col,e,n,lat,lon,utc,amplitude,power,"
          "v_radial_mps,k,b,phi_static_model_name,phi_static_model_rad,C_ati,k_eff_static_phase_df,"
          "p38_mode,geometry_calib_mode,p38_theory_sign,"
          "motion_doppler_axis_sign,ati_phase_to_velocity_sign,"
          "phase_rad,phi_static_total_rad,phi_static_at_zero,phi_res_at_zero,phi_res_rad,"
          "phi_static_geometry_rad,phi_motion_rad,delta_t_s,denom,"
          "denom_without_k,v_from_phase_raw,v_from_phi_res,"
          "v_old_mps,v_iterative_mps,v_analytic_mps,v_root1d_mps,"
          "af_geometry_old_hz,af_geometry_iterative_hz,af_geometry_analytic_hz,af_geometry_root1d_hz,root1d_cost,"
          "p38_pre_k,p38_pre_b,p38_pre_rmse,p38_refit_k,p38_refit_b,p38_refit_rmse,p38_refit_sample_count,p38_refit_inlier_ratio,p38_refit_valid,p38_used_k,p38_used_b,p38_used_source,"
          "phi_static_pre_rad,phi_res_pre_rad,v_pre_mps,phi_static_refit_rad,phi_res_refit_rad,v_refit_mps,"
          "af_phase_hz,af_total_hz,af_geometry_hz,af_motion_hz,"
          "motion_comp_valid,motion_comp_enable,motion_comp_used,motion_comp_fallback,"
          "sinA_old,sinA_comp,sinA_used,old_valid,comp_valid,old_invalid_comp_valid,"
          "loc_used_mode,motion_comp_status,motion_comp_solver,angle_from_sinA_deg,theta_used_for_position_deg,"
          "look_from_sinA_e,look_from_sinA_n,look_e_diff,look_n_diff,old_e,old_n,new_e,new_n,"
          "truth_period_id,truth_beam_id,truth_range_bin,row_truth,row_truth_raw,row_truth_after_fftshift,row_truth_after_recenter,row_truth_float_on_faAxis,"
          "phi_total_truth_rad,phi_static_truth_rad,phi_motion_truth_rad,"
          "v_truth_mps,af_motion_truth_hz,af_geometry_truth_hz,af_total_truth_hz,af_total_detected_hz,"
          "target_e_truth,target_n_truth,"
          "range_m_detected,range_m_truth,ground_range_detected,ground_range_truth,"
          "theta_from_af_geometry_deg,theta_from_truth_position_deg,theta_from_new_position_deg,"
          "look_e_truth,look_n_truth,look_e_used,look_n_used,dx,dy,position_error_m,"
          "phase_truth_error_rad,phi_static_pre_error_rad,phi_static_refit_error_rad,"
          "phi_res_pre_error_rad,phi_res_refit_error_rad,v_error_mps,"
          "af_total_error_hz,af_motion_error_hz,af_geometry_error_hz,"
          "source_file\n";

    for (size_t i = 0; i < write_entries.size(); ++i) {
        const auto &entry = write_entries[i];
        const auto &r = *entry.rec;
        const MovingTruthCsvRow *truth = entry.truth;
        const double power = entry.power;
        os << csv_escape_local(case_id) << ","
           << csv_escape_local(run_id) << ","
           << csv_escape_local(result_id) << ","
           << csv_int_or_blank(r.period_id) << ","
           << csv_int_or_blank(r.beam_id) << ","
           << i << ","
           << csv_int_or_blank(r.range_bin) << ","
           << csv_num(r.range_m) << ","
           << csv_num(r.theta_cmd_deg) << ","
           << csv_num(r.theta_true_deg) << ","
           << csv_int_or_blank(r.row) << ","
           << csv_int_or_blank(r.col) << ","
           << csv_num(r.e) << ","
           << csv_num(r.n) << ","
           << csv_num(r.lat) << ","
           << csv_num(r.lon) << ","
           << csv_num(r.utc) << ","
           << csv_num(r.amplitude) << ","
           << csv_num(power) << ","
           << csv_num(r.radial_velocity_mps) << ","
           << csv_num(r.p38_k) << ","
           << csv_num(r.p38_b) << ","
           << csv_escape_local(r.phi_static_model_name) << ","
           << csv_num(r.phi_static_model_rad) << ","
           << csv_num(r.C_ati) << ","
           << csv_num(r.k_eff_static_phase_df) << ","
           << csv_escape_local(r.p38_mode) << ","
           << csv_escape_local(r.geometry_calib_mode) << ","
           << r.p38_theory_sign << ","
           << r.motion_doppler_axis_sign << ","
           << r.ati_phase_to_velocity_sign << ","
           << csv_num(r.phase_rad) << ","
           << csv_num(r.phi_static_total_rad) << ","
           << csv_num(r.phi_static_at_zero) << ","
           << csv_num(r.phi_res_at_zero) << ","
           << csv_num(r.phi_res_rad) << ","
           << csv_num(r.phi_static_geometry_rad) << ","
           << csv_num(r.phi_motion) << ","
           << csv_num(r.delta_t_s) << ","
           << csv_num(r.motion_comp_denom) << ","
           << csv_num(r.denom_without_k) << ","
           << csv_num(r.v_from_phase_raw) << ","
           << csv_num(r.v_from_phi_res) << ","
           << csv_num(r.v_old_mps) << ","
           << csv_num(r.v_iterative_mps) << ","
           << csv_num(r.v_analytic_mps) << ","
           << csv_num(r.v_root1d_mps) << ","
           << csv_num(r.af_geometry_old_hz) << ","
           << csv_num(r.af_geometry_iterative_hz) << ","
           << csv_num(r.af_geometry_analytic_hz) << ","
           << csv_num(r.af_geometry_root1d_hz) << ","
           << csv_num(r.root1d_cost) << ","
           << csv_num(r.p38_pre_k) << ","
           << csv_num(r.p38_pre_b) << ","
           << csv_num(r.p38_pre_rmse) << ","
           << csv_num(r.p38_refit_k) << ","
           << csv_num(r.p38_refit_b) << ","
           << csv_num(r.p38_refit_rmse) << ","
           << r.p38_refit_sample_count << ","
           << csv_num(r.p38_refit_inlier_ratio) << ","
           << r.p38_refit_valid << ","
           << csv_num(r.p38_used_k) << ","
           << csv_num(r.p38_used_b) << ","
           << csv_escape_local(r.p38_used_source) << ","
           << csv_num(r.phi_static_pre_rad) << ","
           << csv_num(r.phi_res_pre_rad) << ","
           << csv_num(r.v_pre_mps) << ","
           << csv_num(r.phi_static_refit_rad) << ","
           << csv_num(r.phi_res_refit_rad) << ","
           << csv_num(r.v_refit_mps) << ","
           << csv_num(r.af_phase) << ","
           << csv_num(r.af_total) << ","
           << csv_num(r.af_geometry) << ","
           << csv_num(r.af_motion) << ","
           << r.motion_comp_valid << ","
           << r.motion_comp_enable << ","
           << r.motion_comp_used << ","
           << r.motion_comp_fallback << ","
           << csv_num(r.sinA_old) << ","
           << csv_num(r.sinA_comp) << ","
           << csv_num(r.sinA_used) << ","
           << r.old_valid << ","
           << r.comp_valid << ","
           << r.old_invalid_comp_valid << ","
           << csv_escape_local(r.loc_used_mode) << ","
           << csv_escape_local(r.motion_comp_status) << ","
           << csv_escape_local(r.motion_comp_solver) << ","
           << csv_num(r.angle_from_sinA_deg) << ","
           << csv_num(r.theta_used_for_position_deg) << ","
           << csv_num(r.look_from_sinA_e) << ","
           << csv_num(r.look_from_sinA_n) << ","
           << csv_num(r.look_e_diff) << ","
           << csv_num(r.look_n_diff) << ","
           << csv_num(r.old_e) << ","
           << csv_num(r.old_n) << ","
           << csv_num(r.new_e) << ","
           << csv_num(r.new_n) << ",";
        if (truth) {
            const double dx = r.new_e - truth->target_e_truth;
            const double dy = r.new_n - truth->target_n_truth;
            const double position_error_m = std::hypot(dx, dy);
            const double dz = std::isfinite(r.platform_h) ? (r.platform_h - cfg.MT_nowz) : std::numeric_limits<double>::quiet_NaN();
            const double ground_range_detected = horizontal_range_from_slant(r.range_m, dz);
            const double range_m_detected = r.range_m;
            const double ground_range_truth =
                std::hypot(truth->target_e_truth - r.platform_e, truth->target_n_truth - r.platform_n);
            const double range_m_truth = slant_range_from_horizontal(ground_range_truth, dz);
            const double lambda = (cfg.lambda > 0.0) ? cfg.lambda : ((cfg.fc > 0.0) ? (C / (cfg.fc * 1.0e9)) : std::numeric_limits<double>::quiet_NaN());
            const double theta_from_af_geometry_deg =
                theta_from_af_geometry_deg_local(cfg, r.platform_v, r.af_geometry, lambda);
            const double theta_from_truth_position_deg =
                theta_from_position_deg_local(cfg, r.platform_v_angle_deg, r.platform_e, r.platform_n,
                                              truth->target_e_truth, truth->target_n_truth);
            const double theta_from_new_position_deg =
                theta_from_position_deg_local(cfg, r.platform_v_angle_deg, r.platform_e, r.platform_n,
                                              r.new_e, r.new_n);
            const double look_e_truth = (ground_range_truth > 1.0e-9)
                ? ((truth->target_e_truth - r.platform_e) / ground_range_truth)
                : 0.0;
            const double look_n_truth = (ground_range_truth > 1.0e-9)
                ? ((truth->target_n_truth - r.platform_n) / ground_range_truth)
                : 0.0;
            const double look_e_used = (position_error_m >= 0.0 && std::isfinite(r.new_e) && std::isfinite(r.new_n))
                ? ((r.new_e - r.platform_e) / std::max(1.0e-9, std::hypot(r.new_e - r.platform_e, r.new_n - r.platform_n)))
                : 0.0;
            const double look_n_used = (position_error_m >= 0.0 && std::isfinite(r.new_e) && std::isfinite(r.new_n))
                ? ((r.new_n - r.platform_n) / std::max(1.0e-9, std::hypot(r.new_e - r.platform_e, r.new_n - r.platform_n)))
                : 0.0;
            const int row_truth_raw = row_truth_raw_index(*truth, cfg);
            const int row_truth_after_fftshift = row_truth_raw;
            const int row_truth_after_recenter = row_truth_after_recenter_index(*truth, r, cfg);
            const double row_truth_float_on_faAxis = row_truth_float_on_fa_axis(*truth, cfg);
            const double phase_truth_error_rad =
                wrap_pi_csv(r.phase_rad - truth->phi_total_truth_rad);
            const double phi_static_pre_error_rad =
                wrap_pi_csv(r.phi_static_pre_rad - truth->phi_static_truth_rad);
            const double phi_static_refit_error_rad =
                wrap_pi_csv(r.phi_static_refit_rad - truth->phi_static_truth_rad);
            const double phi_res_pre_error_rad =
                wrap_pi_csv(r.phi_res_pre_rad - truth->phi_motion_truth_rad);
            const double phi_res_refit_error_rad =
                wrap_pi_csv(r.phi_res_refit_rad - truth->phi_motion_truth_rad);
            const double v_error_mps = r.radial_velocity_mps - truth->v_truth_mps;
            const double af_total_error_hz = r.af_total - truth->af_total_truth_hz;
            const double af_motion_error_hz = r.af_motion - truth->af_motion_truth_hz;
            const double af_geometry_error_hz = r.af_geometry - truth->af_geometry_truth_hz;
            os << csv_int_or_blank(truth->period_id) << ","
               << csv_int_or_blank(truth->beam_id) << ","
               << csv_int_or_blank(truth->range_bin) << ","
               << csv_int_or_blank(truth->row_truth) << ","
               << csv_int_or_blank(row_truth_raw) << ","
               << csv_int_or_blank(row_truth_after_fftshift) << ","
               << csv_int_or_blank(row_truth_after_recenter) << ","
               << csv_num(row_truth_float_on_faAxis) << ","
               << csv_num(truth->phi_total_truth_rad) << ","
               << csv_num(truth->phi_static_truth_rad) << ","
               << csv_num(truth->phi_motion_truth_rad) << ","
               << csv_num(truth->v_truth_mps) << ","
               << csv_num(truth->af_motion_truth_hz) << ","
               << csv_num(truth->af_geometry_truth_hz) << ","
               << csv_num(truth->af_total_truth_hz) << ","
               << csv_num(r.af_total) << ","
               << csv_num(truth->target_e_truth) << ","
               << csv_num(truth->target_n_truth) << ","
               << csv_num(range_m_detected) << ","
               << csv_num(range_m_truth) << ","
               << csv_num(ground_range_detected) << ","
               << csv_num(ground_range_truth) << ","
               << csv_num(theta_from_af_geometry_deg) << ","
               << csv_num(theta_from_truth_position_deg) << ","
               << csv_num(theta_from_new_position_deg) << ","
               << csv_num(look_e_truth) << ","
               << csv_num(look_n_truth) << ","
               << csv_num(look_e_used) << ","
               << csv_num(look_n_used) << ","
               << csv_num(dx) << ","
               << csv_num(dy) << ","
               << csv_num(position_error_m) << ","
               << csv_num(phase_truth_error_rad) << ","
               << csv_num(phi_static_pre_error_rad) << ","
               << csv_num(phi_static_refit_error_rad) << ","
               << csv_num(phi_res_pre_error_rad) << ","
               << csv_num(phi_res_refit_error_rad) << ","
               << csv_num(v_error_mps) << ","
               << csv_num(af_total_error_hz) << ","
               << csv_num(af_motion_error_hz) << ","
               << csv_num(af_geometry_error_hz) << ",";
        } else {
            for (int blank_i = 0; blank_i < 41; ++blank_i) {
                if (blank_i > 0) {
                    os << ",";
                }
            }
        }
        os
           << csv_escape_local(src) << "\n";
    }

    if (!os) {
        std::cerr << "writeDetectionCsv: 写文件失败: " << outpath << "\n";
        return false;
    }
    std::cout << "[GMTI] detection_results_csv = " << outpath
              << "，目标数: " << records.size() << std::endl;
    return true;
}

// ---------- 小工具：递归创建目录 ----------
bool GMTIProcessor::mkdir_p(const std::string &dir)
{
    if (dir.empty())
        return true;
    std::string cur;
    size_t i = 0;
    if (dir[0] == '/')
    {
        cur = "/";
        i = 1;
    }
    for (; i < dir.size(); ++i)
    {
        cur.push_back(dir[i]);
        if (dir[i] == '/' || i == dir.size() - 1)
        {
            if (cur == "/" || cur == "./" || cur == ".")
                continue;
            if (::mkdir(cur.c_str(), 0755) != 0)
            {
                if (errno == EEXIST)
                    continue;
                if (errno == ENOENT)
                    continue; // 父目录尚未就绪，继续循环
                std::cerr << "mkdir_p failed: " << cur << " errno=" << errno << "\n";
                return false;
            }
        }
    }
    return true;
}

// ---------- 占位：xy(m) -> (lat,lng) 度 ----------

// ---------- 写入帮助：LE 编码 ----------
static inline void put_i32_le(std::vector<uint8_t> &buf, size_t pos, int32_t v)
{
    uint32_t u = static_cast<uint32_t>(v);
    buf[pos + 0] = static_cast<uint8_t>(u & 0xFF);
    buf[pos + 1] = static_cast<uint8_t>((u >> 8) & 0xFF);
    buf[pos + 2] = static_cast<uint8_t>((u >> 16) & 0xFF);
    buf[pos + 3] = static_cast<uint8_t>((u >> 24) & 0xFF);
}

// 量化：度 -> int32（LSB=8.38191e-8 度），四舍五入并饱和
static inline int32_t quant_deg_to_i32(double deg)
{
    const double LSB = 8.38191e-8; // deg / LSB
    double q = deg / LSB;
    // 四舍五入
    if (q >= 0.0)
        q = std::floor(q + 0.5);
    else
        q = std::ceil(q - 0.5);
    // 饱和到 int32
    if (q > 2147483647.0)
        q = 2147483647.0;
    if (q < -2147483648.0)
        q = -2147483648.0;
    return static_cast<int32_t>(q);
}

// ---------- 核心：写入所有 tracks 的所有 pos ----------
bool GMTIProcessor::writeTracksBinary(const std::vector<Track> &tracks,
                                      double utcMid,
                                      const GMTIOutput::Plane &plane,
                                      const Config &cfg)
{
    (void)utcMid;
    (void)plane;

    // 逐目标写 28 字节记录，文件尾追加 8 字节 utc
    const size_t REC_BYTES = 28;
    const size_t HEADER_BYTES = 2;
    std::vector<uint8_t> buf;

    // 1) 将所有轨迹的滤波点 kf 展开为点序列（按 track.id 升序；每条内按时间顺序）
    struct Rec
    {
        uint16_t id;
        double x;
        double y;
        double speed;
        double direction;
        double range;
    };
    std::vector<Rec> points;
    points.reserve(64);

    // 拷贝并排序轨迹（确保 id 升序；id<=0 的放后面）
    std::vector<const Track *> order;
    order.reserve(tracks.size());
    for (size_t i = 0; i < tracks.size(); ++i)
        order.push_back(&tracks[i]);
    std::stable_sort(order.begin(), order.end(),
                     [](const Track *a, const Track *b)
                     {
                         const int ia = (a->id > 0) ? a->id : 0x7fffffff;
                         const int ib = (b->id > 0) ? b->id : 0x7fffffff;
                         if (ia != ib)
                             return ia < ib;
                         return a < b;
                     });

    for (size_t ti = 0; ti < order.size(); ++ti)
    {
        const Track &tr = *order[ti];
        const uint16_t tid = (tr.id > 0 && tr.id <= 65535) ? static_cast<uint16_t>(tr.id)
                                                           : static_cast<uint16_t>(std::min<size_t>(65535, ti + 1));
        // 只写入每条轨迹的第一个滤波点
        if (tr.kf.empty())
            continue;
        double speed = 0.0;
        if (tr.pos.size() >= 2 && tr.time.size() >= 2)
        {
            const size_t n = tr.pos.size();
            const double dx = tr.pos[n - 1][0] - tr.pos[n - 2][0];
            const double dy = tr.pos[n - 1][1] - tr.pos[n - 2][1];
            const double dt = tr.time[n - 1] - tr.time[n - 2];
            speed = std::sqrt(dx * dx + dy * dy) / (dt > 1e-6 ? dt : 1.0);
        }
        Rec r;
        r.id = tid;
        r.x = tr.kf.back()[0];
        r.y = tr.kf.back()[1];
        r.speed = speed;
        r.direction = tr.direction;
        r.range = tr.range;
        points.push_back(r);
    }

    std::cout << "[GMTI] writeTracksBinary: 共 " << tracks.size() << " 个目标，"
              << points.size() << " 个点。\n";

    // 2 字节目标数 + N*28 字节记录 + 8 字节 utc
    buf.reserve(HEADER_BYTES + points.size() * REC_BYTES + 8);
    buf.resize(HEADER_BYTES, 0u);
    put_u16_le(buf, 0, static_cast<uint16_t>(std::min<size_t>(65535, points.size())));

    // 逐点写 28 字节记录
    for (size_t i = 0; i < points.size(); ++i)
    {
        const size_t base = HEADER_BYTES + i * REC_BYTES;
        buf.resize(base + REC_BYTES, 0u);

        // 2 byte id
        put_u16_le(buf, base + 0, points[i].id);

        // 坐标转换：x,y(m) → lat,lng(deg)
        double lat = 0.0, lng = 0.0;
        Gaussp3RV(points[i].x, points[i].y, cfg.L0, lat, lng);

        // 2..5 经度，6..9 纬度
        put_i32_le(buf, base + 2, quant_deg_to_i32(lng));
        put_i32_le(buf, base + 6, quant_deg_to_i32(lat));
        // 10..11 速度：轨迹末端速度估计
        put_u16_le(buf, base + 10, quant_speed_to_u16(points[i].speed));
        // 12..19 方向，20..27 距离
        put_f64_le(buf, base + 12, wrap180_deg(points[i].direction));
        put_f64_le(buf, base + 20, points[i].range);
    }

    put_f64_le(buf, HEADER_BYTES + points.size() * REC_BYTES, utcMid);

    // 输出到 result_add 目录
    std::string outpath;
    int idx = 0;
    if (!nextGMTIFileName(cfg.result_add, outpath, idx))
    {
        return false; // 或者回退到固定文件名
    }

    std::ofstream ofs(outpath.c_str(), std::ios::binary);
    if (!ofs)
    {
        std::cerr << "writeTracksBinary: 无法打开输出文件: " << outpath << "\n";
        return false;
    }
    ofs.write(reinterpret_cast<const char *>(buf.data()),
              static_cast<std::streamsize>(buf.size()));
    if (!ofs)
    {
        std::cerr << "writeTracksBinary: 写文件失败: " << outpath << "\n";
        return false;
    }
    std::cout << "writeTracksBinary: 写文件成功: " << outpath << "\n";
    return true;
}
