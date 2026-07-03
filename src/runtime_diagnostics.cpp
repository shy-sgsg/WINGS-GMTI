#include "runtime_diagnostics.hpp"

#include "config_structs.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#ifndef GMTI_BUILD_TYPE
#ifdef NDEBUG
#define GMTI_BUILD_TYPE "Release"
#else
#define GMTI_BUILD_TYPE "Debug"
#endif
#endif

namespace gmti {
namespace runtime {
namespace {

struct TimingMetric {
    std::string case_id;
    std::string run_id;
    std::string result_id;
    int period_id;
    std::string scope_name;
    std::string start_time;
    std::string end_time;
    long long elapsed_ms;
    std::string extra;
};

struct RunState {
    bool initialized = false;
    bool diagnostics_enabled = false;
    std::string case_id;
    std::string run_id;
    std::string result_id;
    std::string start_time;
    std::string end_time;
    std::string git_commit = "unknown";
    std::string build_type = GMTI_BUILD_TYPE;
    std::string executable;
    std::string config_path;
    std::string working_directory;
    std::string input_data_path;
    std::string output_dir;
    RunPaths paths;
    std::vector<TimingMetric> timings;
};

std::mutex& stateMutex()
{
    static std::mutex m;
    return m;
}

RunState& state()
{
    static RunState s;
    return s;
}

std::string jsonEscape(const std::string& s)
{
    std::ostringstream os;
    for (char ch : s) {
        switch (ch) {
        case '\\': os << "\\\\"; break;
        case '"': os << "\\\""; break;
        case '\n': os << "\\n"; break;
        case '\r': os << "\\r"; break;
        case '\t': os << "\\t"; break;
        default:
            if (static_cast<unsigned char>(ch) < 0x20) {
                os << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                   << static_cast<int>(static_cast<unsigned char>(ch))
                   << std::dec << std::setfill(' ');
            } else {
                os << ch;
            }
        }
    }
    return os.str();
}

std::string q(const std::string& s)
{
    return "\"" + jsonEscape(s) + "\"";
}

std::string jsonNullableDouble(bool available, double value)
{
    if (!available) {
        return "null";
    }
    std::ostringstream os;
    os << std::setprecision(15) << value;
    return os.str();
}

std::string csvEscape(const std::string& s)
{
    if (s.find_first_of(",\"\n\r") == std::string::npos) {
        return s;
    }
    std::string out = "\"";
    for (char ch : s) {
        if (ch == '"') out += "\"\"";
        else out.push_back(ch);
    }
    out += "\"";
    return out;
}

std::string pathJoin(std::string dir, const std::string& file)
{
    if (dir.empty()) {
        return file;
    }
    if (dir.back() != '/' && dir.back() != '\\') {
        dir.push_back('/');
    }
    return dir + file;
}

bool mkdirP(const std::string& dir)
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
            if (cur == "/" || cur == "./" || cur == ".") {
                continue;
            }
            if (::mkdir(cur.c_str(), 0755) != 0 && errno != EEXIST) {
                return false;
            }
        }
    }
    return true;
}

bool fileExists(const std::string& path)
{
    struct stat st {};
    return !path.empty() && ::stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

std::string basenameNoExt(const std::string& path)
{
    const size_t slash = path.find_last_of("/\\");
    std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
    const size_t dot = name.find_last_of('.');
    if (dot != std::string::npos && dot > 0) {
        name = name.substr(0, dot);
    }
    return name.empty() ? "unknown_case" : name;
}

std::string currentWorkingDirectory()
{
    char buf[4096];
    if (::getcwd(buf, sizeof(buf))) {
        return std::string(buf);
    }
    return "unknown";
}

std::string makeRunIdFromTime()
{
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
    localtime_r(&t, &tm);
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(
        now.time_since_epoch()).count() % 1000000;
    std::ostringstream os;
    os << std::put_time(&tm, "%Y%m%d_%H%M%S") << "_"
       << std::setw(6) << std::setfill('0') << micros;
    return os.str();
}

std::string runCommandOneLine(const char* cmd)
{
    FILE* pipe = ::popen(cmd, "r");
    if (!pipe) return "unknown";
    char buf[256];
    std::string out;
    if (std::fgets(buf, sizeof(buf), pipe)) {
        out = buf;
    }
    const int rc = ::pclose(pipe);
    if (rc != 0 || out.empty()) {
        return "unknown";
    }
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r' || out.back() == ' ')) {
        out.pop_back();
    }
    return out.empty() ? "unknown" : out;
}

std::string protocolType(const Config& cfg)
{
    return cfg.INFO_Type ? "new_protocol" : "legacy_protocol";
}

std::string payloadFormat(const Config& cfg)
{
    if (cfg.INFO_Type) return "float32_ch1_iq_ch2_iq";
    if (cfg.channel_mode == "separate") return "float32_single_channel_iq";
    if (cfg.channel_mode == "interleaved") return "float32_ch1_iq_ch2_iq";
    return "not_available";
}

std::string assignmentModeName(int v)
{
    if (v == 0) return "greedy_nearest_neighbor";
    if (v == 1) return "hungarian";
    return "unknown";
}

std::string distanceModeName(int v)
{
    if (v == 0) return "euclidean";
    if (v == 1) return "mahalanobis_squared";
    return "unknown";
}

int beamCount(const Config& cfg)
{
    if (cfg.wavepos_skip <= 0 || cfg.wavepos_ed < cfg.wavepos_st) {
        return cfg.az_count;
    }
    return (cfg.wavepos_ed - cfg.wavepos_st) / cfg.wavepos_skip + 1;
}

double scanStepDeg(const Config& cfg)
{
    const int count = beamCount(cfg);
    if (count <= 1) return 0.0;
    return (cfg.scan_max_deg - cfg.scan_min_deg) / static_cast<double>(count - 1);
}

long long bytesPerPeriod(const Config& cfg)
{
    return static_cast<long long>(beamCount(cfg)) *
           static_cast<long long>(effectivePulseNum(cfg)) *
           static_cast<long long>(cfg.pkg_bytes);
}

std::vector<std::string> listFilesWithPrefixSuffix(const std::string& dir,
                                                   const std::string& prefix,
                                                   const std::string& suffix,
                                                   const std::string& exclude_suffix = "")
{
    std::vector<std::string> out;
    DIR* dp = ::opendir(dir.c_str());
    if (!dp) return out;
    while (dirent* ent = ::readdir(dp)) {
        const std::string name = ent->d_name;
        if (name == "." || name == "..") continue;
        const bool prefix_ok = prefix.empty() || name.compare(0, prefix.size(), prefix) == 0;
        const bool suffix_ok = suffix.empty() ||
            (name.size() >= suffix.size() &&
             name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0);
        const bool exclude_ok = exclude_suffix.empty() ||
            name.size() < exclude_suffix.size() ||
            name.compare(name.size() - exclude_suffix.size(),
                         exclude_suffix.size(),
                         exclude_suffix) != 0;
        if (prefix_ok && suffix_ok && exclude_ok) {
            out.push_back(pathJoin(dir, name));
        }
    }
    ::closedir(dp);
    std::sort(out.begin(), out.end());
    return out;
}

std::string jsonStringArray(const std::vector<std::string>& values, int indent)
{
    if (values.empty()) return "[]";
    std::ostringstream os;
    os << "[\n";
    const std::string pad(static_cast<size_t>(indent), ' ');
    for (size_t i = 0; i < values.size(); ++i) {
        os << pad << q(values[i]);
        if (i + 1 < values.size()) os << ",";
        os << "\n";
    }
    os << std::string(static_cast<size_t>(indent - 2), ' ') << "]";
    return os.str();
}

std::string latestFile(const std::vector<std::string>& files)
{
    return files.empty() ? "" : files.back();
}

void writeTimingCsvUnlocked()
{
    const RunState& s = state();
    if (s.paths.timing_metrics_csv.empty()) return;
    std::ofstream ofs(s.paths.timing_metrics_csv.c_str());
    if (!ofs) {
        std::cerr << "[TIMING][WARN] cannot write " << s.paths.timing_metrics_csv << std::endl;
        return;
    }
    ofs << "case_id,run_id,result_id,period_id,scope_name,start_time,end_time,elapsed_ms,extra\n";
    for (const auto& m : s.timings) {
        ofs << csvEscape(m.case_id) << ","
            << csvEscape(m.run_id) << ","
            << csvEscape(m.result_id) << ","
            << m.period_id << ","
            << csvEscape(m.scope_name) << ","
            << csvEscape(m.start_time) << ","
            << csvEscape(m.end_time) << ","
            << m.elapsed_ms << ","
            << csvEscape(m.extra) << "\n";
    }
}

void writeManifest(const Config& cfg, bool normal_exit, int exit_code,
                   const std::string& notes)
{
    RunState& s = state();
    if (s.paths.run_manifest_json.empty()) return;

    const std::vector<std::string> gmti_bins = listFilesWithPrefixSuffix(cfg.result_add, "GMTI", ".bin", "_track.bin");
    const std::vector<std::string> gmti_tracks = listFilesWithPrefixSuffix(cfg.result_add, "GMTI", "_track.bin");
    const std::vector<std::string> pngs = listFilesWithPrefixSuffix(cfg.result_add, "", ".png");
    const std::vector<std::string> txts = listFilesWithPrefixSuffix(cfg.result_add, "", ".txt");
    const std::string detection_csv = pathJoin(cfg.result_add, "detection_results.csv");
    const std::string calibration_dir = pathJoin(cfg.result_add, "calibration");
    const std::string single_point_calibration_json =
        pathJoin(calibration_dir, "single_point_bin_calibration.json");
    const std::string single_point_calibration_csv =
        pathJoin(calibration_dir, "single_point_bin_calibration.csv");
    const std::string single_point_calibration_report =
        pathJoin(calibration_dir, "single_point_calibration_report.md");

    std::ofstream os(s.paths.run_manifest_json.c_str());
    if (!os) {
        std::cerr << "[MANIFEST][WARN] cannot write " << s.paths.run_manifest_json << std::endl;
        return;
    }
    os << "{\n";
    os << "  \"case_id\": " << q(s.case_id) << ",\n";
    os << "  \"run_id\": " << q(s.run_id) << ",\n";
    os << "  \"result_id\": " << q(s.result_id) << ",\n";
    os << "  \"git_commit\": " << q(s.git_commit) << ",\n";
    os << "  \"build_type\": " << q(s.build_type) << ",\n";
    os << "  \"executable\": " << q(s.executable) << ",\n";
    os << "  \"config_path\": " << q(s.config_path) << ",\n";
    os << "  \"input_data_path\": " << q(s.input_data_path) << ",\n";
    os << "  \"output_dir\": " << q(s.output_dir) << ",\n";
    os << "  \"start_time\": " << q(s.start_time) << ",\n";
    os << "  \"end_time\": " << q(s.end_time) << ",\n";
    os << "  \"normal_exit\": " << (normal_exit ? "true" : "false") << ",\n";
    os << "  \"exit_code\": " << exit_code << ",\n";
    os << "  \"generated_files\": {\n";
    os << "    \"runtime_config_dump_json\": " << q(s.paths.runtime_config_json) << ",\n";
    os << "    \"runtime_config_dump_txt\": " << q(s.paths.runtime_config_txt) << ",\n";
    os << "    \"timing_metrics_csv\": " << q(s.paths.timing_metrics_csv) << ",\n";
    os << "    \"algorithm_log\": " << q("not_available") << ",\n";
    os << "    \"gmti_bin\": " << q(latestFile(gmti_bins)) << ",\n";
    os << "    \"gmti_track_bin\": " << q(latestFile(gmti_tracks)) << ",\n";
    os << "    \"detection_results_csv\": " << q(fileExists(detection_csv) ? detection_csv : "") << ",\n";
    os << "    \"single_point_bin_calibration_json\": " << q(fileExists(single_point_calibration_json) ? single_point_calibration_json : "") << ",\n";
    os << "    \"single_point_bin_calibration_csv\": " << q(fileExists(single_point_calibration_csv) ? single_point_calibration_csv : "") << ",\n";
    os << "    \"single_point_calibration_report\": " << q(fileExists(single_point_calibration_report) ? single_point_calibration_report : "") << ",\n";
    os << "    \"png_outputs\": " << jsonStringArray(pngs, 6) << ",\n";
    os << "    \"txt_outputs\": " << jsonStringArray(txts, 6) << "\n";
    os << "  },\n";
    os << "  \"notes\": " << q(notes) << "\n";
    os << "}\n";
}

void writeRuntimeConfigJson(const Config& cfg)
{
    RunState& s = state();
    std::ofstream os(s.paths.runtime_config_json.c_str());
    if (!os) {
        std::cerr << "[PARAM_DUMP][WARN] cannot write "
                  << s.paths.runtime_config_json << std::endl;
        return;
    }
    const int proc_pulses = effectivePulseNum(cfg);
    const double period_time_sec =
        cfg.PRF > 0.0 ? static_cast<double>(beamCount(cfg)) * proc_pulses / cfg.PRF : 0.0;
    const double doppler_bin_hz = proc_pulses > 0 ? cfg.PRF / proc_pulses : 0.0;
    const double velocity_res = cfg.lambda * doppler_bin_hz / 2.0;
    const double first_blind = cfg.lambda * cfg.PRF / 2.0;

    os << "{\n";
    os << "  \"run_info\": {\n";
    os << "    \"case_id\": " << q(s.case_id) << ",\n";
    os << "    \"run_id\": " << q(s.run_id) << ",\n";
    os << "    \"result_id\": " << q(s.result_id) << ",\n";
    os << "    \"start_time\": " << q(s.start_time) << ",\n";
    os << "    \"git_commit\": " << q(s.git_commit) << ",\n";
    os << "    \"build_type\": " << q(s.build_type) << ",\n";
    os << "    \"runtime_mode\": " << q(cfg.runtime_mode) << ",\n";
    os << "    \"runtime_diagnostics_enabled\": " << (cfg.runtime_diagnostics_enabled ? "true" : "false") << ",\n";
    os << "    \"executable\": " << q(s.executable) << ",\n";
    os << "    \"config_path\": " << q(s.config_path) << ",\n";
    os << "    \"working_directory\": " << q(s.working_directory) << "\n";
    os << "  },\n";
    os << "  \"file_layout\": {\n";
    os << "    \"input_data_path\": " << q(s.input_data_path) << ",\n";
    os << "    \"output_dir\": " << q(cfg.result_add) << ",\n";
    os << "    \"protocol_type\": " << q(protocolType(cfg)) << ",\n";
    os << "    \"header_size_bytes\": " << cfg.info_len << ",\n";
    os << "    \"payload_format\": " << q(payloadFormat(cfg)) << ",\n";
    os << "    \"pulse_len\": " << cfg.pulse_len << ",\n";
    os << "    \"read_pulse_num\": " << cfg.read_pulse_num << ",\n";
    os << "    \"pulse_num\": " << cfg.pulse_num << ",\n";
    os << "    \"beam_count\": " << beamCount(cfg) << ",\n";
    os << "    \"wavepos_st\": " << cfg.wavepos_st << ",\n";
    os << "    \"wavepos_ed\": " << cfg.wavepos_ed << ",\n";
    os << "    \"wavepos_skip\": " << cfg.wavepos_skip << ",\n";
    os << "    \"bytes_per_prt\": " << cfg.pkg_bytes << ",\n";
    os << "    \"prt_per_period\": " << beamCount(cfg) * proc_pulses << "\n";
    os << "  },\n";
    os << "  \"waveform\": {\n";
    os << "    \"fc_hz\": " << std::setprecision(15) << cfg.fc << ",\n";
    os << "    \"Br_hz\": " << cfg.Br << ",\n";
    os << "    \"fs_hz\": " << cfg.fs << ",\n";
    os << "    \"Tr_sec\": " << cfg.Tr << ",\n";
    os << "    \"PRF_hz\": " << cfg.PRF << ",\n";
    os << "    \"lambda_m\": " << cfg.lambda << ",\n";
    os << "    \"sample_delay_us\": "
       << jsonNullableDouble(cfg.has_sample_delay_us, cfg.sample_delay_us) << "\n";
    os << "  },\n";
    os << "  \"pulse_compression\": {\n";
    os << "    \"range_fft_len\": " << effectiveRangeFftLen(cfg) << ",\n";
    os << "    \"fft_pulse_len\": " << effectiveRangeFftLen(cfg) << ",\n";
    os << "    \"range_crop_start\": " << cfg.range_crop_start << ",\n";
    os << "    \"range_compress_len\": " << cfg.range_compress_len << ",\n";
    os << "    \"pc_crop_start\": " << cfg.range_crop_start << ",\n";
    os << "    \"pc_crop_len\": " << cfg.range_compress_len << "\n";
    os << "  },\n";
    os << "  \"geometry\": {\n";
    os << "    \"scan_min_deg\": " << cfg.scan_min_deg << ",\n";
    os << "    \"scan_max_deg\": " << cfg.scan_max_deg << ",\n";
    os << "    \"scan_step_deg\": " << scanStepDeg(cfg) << ",\n";
    os << "    \"az_count\": " << cfg.az_count << ",\n";
    os << "    \"beam_width_deg\": " << cfg.beamwidth_deg << ",\n";
    os << "    \"squint_angle\": " << cfg.squint_angle << ",\n";
    os << "    \"estimate_error_angle\": " << (cfg.estimate_error_angle ? "true" : "false") << ",\n";
    os << "    \"beam_pointing_bias\": null,\n";
    os << "    \"beam_pointing_bias_source\": " << q("no_explicit_algorithm_config_field") << ",\n";
    os << "    \"loc_beam_gate_deg\": " << cfg.loc_beam_gate_deg << "\n";
    os << "  },\n";
    os << "  \"channel_and_gmti\": {\n";
    os << "    \"d_chan\": " << cfg.d_channel << ",\n";
    os << "    \"calib_coef\": " << cfg.calib_coef << ",\n";
    os << "    \"channel_phase_enabled\": " << q("called_in_processOnePeriod") << ",\n";
    os << "    \"channel_phase_source\": " << q("rg_correct_CUDA applied in processOnePeriod/processOnePeriodFusionCache") << ",\n";
    os << "    \"cancellation_enabled\": " << q("called_in_processOnePeriod") << ",\n";
    os << "    \"cancellation_source\": " << q("clutter_cancel_38_paper_1_p38_cuda and clutter_cancel_38_paper_1_cuda called in processOnePeriod/processOnePeriodFusionCache") << "\n";
    os << "  },\n";
    os << "  \"detection\": {\n";
    os << "    \"pf\": " << cfg.pf << ",\n";
    os << "    \"threshold\": null,\n";
    os << "    \"threshold_mode\": " << q("runtime_adaptive") << ",\n";
    os << "    \"min_points\": " << cfg.min_points << ",\n";
    os << "    \"min_len\": " << cfg.min_len << ",\n";
    os << "    \"cfar_type\": " << q("GO") << ",\n";
    os << "    \"cfar_guard_cells\": 4,\n";
    os << "    \"cfar_background_cells\": 16,\n";
    os << "    \"cluster_max_gap\": 2,\n";
    os << "    \"cluster_max_phase_std\": 0.2,\n";
    os << "    \"cfar_params\": " << q("not_structured_yet") << "\n";
    os << "  },\n";
    os << "  \"tracking\": {\n";
    os << "    \"enabled\": true,\n";
    os << "    \"assignment_mode\": " << q(assignmentModeName(cfg.track_assignment_mode)) << ",\n";
    os << "    \"distance_mode\": " << q(distanceModeName(cfg.track_distance_mode)) << ",\n";
    os << "    \"confirm_hits\": " << cfg.track_confirm_hits << ",\n";
    os << "    \"delete_misses\": " << cfg.track_max_missed << ",\n";
    os << "    \"gate_threshold\": " << cfg.track_gate_m << ",\n";
    os << "    \"track_manager_params\": {\n";
    os << "      \"track_idx_window\": " << cfg.track_idx_window << ",\n";
    os << "      \"track_truth_threshold\": " << cfg.track_truth_threshold << ",\n";
    os << "      \"track_confirm_window\": " << cfg.track_confirm_window << ",\n";
    os << "      \"track_tentative_max_missed\": " << cfg.track_tentative_max_missed << ",\n";
    os << "      \"track_v_max\": " << cfg.track_v_max << ",\n";
    os << "      \"track_chi2_gate\": " << cfg.track_chi2_gate << ",\n";
    os << "      \"track_default_dt\": " << cfg.track_default_dt << ",\n";
    os << "      \"track_process_noise_pos\": " << cfg.track_process_noise_pos << ",\n";
    os << "      \"track_process_noise_vel\": " << cfg.track_process_noise_vel << ",\n";
    os << "      \"track_measurement_noise_pos\": " << cfg.track_measurement_noise_pos << "\n";
    os << "    }\n";
    os << "  },\n";
    os << "  \"derived\": {\n";
    os << "    \"period_time_sec\": " << period_time_sec << ",\n";
    os << "    \"realtime_threshold_80pct_sec\": " << period_time_sec * 0.8 << ",\n";
    os << "    \"prt_per_period\": " << beamCount(cfg) * proc_pulses << ",\n";
    os << "    \"bytes_per_period\": " << bytesPerPeriod(cfg) << ",\n";
    os << "    \"range_resolution_m\": " << (cfg.Br > 0.0 ? C / (2.0 * cfg.Br) : 0.0) << ",\n";
    os << "    \"range_bin_m\": " << cfg.R_bin << ",\n";
    os << "    \"doppler_bin_hz\": " << doppler_bin_hz << ",\n";
    os << "    \"velocity_resolution_mps\": " << velocity_res << ",\n";
    os << "    \"first_blind_speed_mps\": " << first_blind << "\n";
    os << "  }\n";
    os << "}\n";
}

void writeRuntimeConfigTxt(const Config& cfg)
{
    RunState& s = state();
    std::ofstream os(s.paths.runtime_config_txt.c_str());
    if (!os) {
        std::cerr << "[PARAM_DUMP][WARN] cannot write "
                  << s.paths.runtime_config_txt << std::endl;
        return;
    }
    const int proc_pulses = effectivePulseNum(cfg);
    const double period_time_sec =
        cfg.PRF > 0.0 ? static_cast<double>(beamCount(cfg)) * proc_pulses / cfg.PRF : 0.0;

    os << "[RUN]\n";
    os << "case_id = " << s.case_id << "\n";
    os << "run_id = " << s.run_id << "\n";
    os << "config_path = " << s.config_path << "\n";
    os << "output_dir = " << cfg.result_add << "\n";
    os << "git_commit = " << s.git_commit << "\n\n";
    os << "runtime_mode = " << cfg.runtime_mode << "\n";
    os << "runtime_diagnostics_enabled = "
       << (cfg.runtime_diagnostics_enabled ? "true" : "false") << "\n\n";

    os << "[FILE_LAYOUT]\n";
    os << "input_data_path = " << s.input_data_path << "\n";
    os << "pulse_len = " << cfg.pulse_len << "\n";
    os << "read_pulse_num = " << cfg.read_pulse_num << "\n";
    os << "pulse_num = " << cfg.pulse_num << "\n";
    os << "beam_count = " << beamCount(cfg) << "\n";
    os << "bytes_per_prt = " << cfg.pkg_bytes << "\n\n";

    os << "[WAVEFORM]\n";
    os << "fc_hz = " << std::setprecision(15) << cfg.fc << "\n";
    os << "Br_hz = " << cfg.Br << "\n";
    os << "fs_hz = " << cfg.fs << "\n";
    os << "Tr_sec = " << cfg.Tr << "\n";
    os << "PRF_hz = " << cfg.PRF << "\n";
    os << "lambda_m = " << cfg.lambda << "\n\n";
    os << "sample_delay_us = "
       << jsonNullableDouble(cfg.has_sample_delay_us, cfg.sample_delay_us) << "\n\n";

    os << "[PULSE_COMPRESSION]\n";
    os << "range_fft_len = " << effectiveRangeFftLen(cfg) << "\n";
    os << "pc_crop_start = " << cfg.range_crop_start << "\n";
    os << "pc_crop_len = " << cfg.range_compress_len << "\n\n";

    os << "[GEOMETRY]\n";
    os << "scan_min_deg = " << cfg.scan_min_deg << "\n";
    os << "scan_max_deg = " << cfg.scan_max_deg << "\n";
    os << "scan_step_deg = " << scanStepDeg(cfg) << "\n";
    os << "az_count = " << cfg.az_count << "\n";
    os << "squint_angle = " << cfg.squint_angle << "\n";
    os << "estimate_error_angle = " << (cfg.estimate_error_angle ? "true" : "false") << "\n\n";
    os << "beam_pointing_bias = null\n";
    os << "beam_pointing_bias_source = no_explicit_algorithm_config_field\n\n";

    os << "[CHANNEL_AND_GMTI]\n";
    os << "d_chan = " << cfg.d_channel << "\n";
    os << "calib_coef = " << cfg.calib_coef << "\n";
    os << "channel_phase_enabled = called_in_processOnePeriod\n";
    os << "channel_phase_source = rg_correct_CUDA applied in processOnePeriod/processOnePeriodFusionCache\n";
    os << "cancellation_enabled = called_in_processOnePeriod\n";
    os << "cancellation_source = clutter_cancel_38_paper_1_p38_cuda and clutter_cancel_38_paper_1_cuda called in processOnePeriod/processOnePeriodFusionCache\n\n";

    os << "[DETECTION]\n";
    os << "pf = " << cfg.pf << "\n";
    os << "min_points = " << cfg.min_points << "\n";
    os << "min_len = " << cfg.min_len << "\n";
    os << "threshold = null\n";
    os << "threshold_mode = runtime_adaptive\n";
    os << "cfar_type = GO\n";
    os << "cfar_guard_cells = 4\n";
    os << "cfar_background_cells = 16\n";
    os << "cluster_max_gap = 2\n";
    os << "cluster_max_phase_std = 0.2\n";
    os << "cfar_params = not_structured_yet\n\n";

    os << "[TRACKING]\n";
    os << "assignment_mode = " << assignmentModeName(cfg.track_assignment_mode) << "\n";
    os << "distance_mode = " << distanceModeName(cfg.track_distance_mode) << "\n";
    os << "confirm_hits = " << cfg.track_confirm_hits << "\n";
    os << "delete_misses = " << cfg.track_max_missed << "\n";
    os << "track_gate_m = " << cfg.track_gate_m << "\n\n";

    os << "[DERIVED]\n";
    os << "period_time_sec = " << period_time_sec << "\n";
    os << "realtime_threshold_80pct_sec = " << period_time_sec * 0.8 << "\n";
    os << "prt_per_period = " << beamCount(cfg) * proc_pulses << "\n";
}

} // namespace

std::string currentIsoTime()
{
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
    localtime_r(&t, &tm);
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(
        now.time_since_epoch()).count() % 1000000;
    std::ostringstream os;
    os << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S") << "."
       << std::setw(6) << std::setfill('0') << micros;
    return os.str();
}

void initializeRun(const Config& cfg,
                   const std::string& config_path,
                   const std::string& executable_path)
{
    RunState& s = state();
    s.initialized = true;
    s.diagnostics_enabled = cfg.runtime_diagnostics_enabled;
    s.case_id = basenameNoExt(config_path);
    s.run_id = makeRunIdFromTime();
    std::ostringstream rid;
    rid << "GMTI";
    if (cfg.result_file_id > 0) {
        rid << std::setw(2) << std::setfill('0') << cfg.result_file_id;
    } else {
        rid << "auto";
    }
    s.result_id = rid.str();
    s.start_time = currentIsoTime();
    s.end_time.clear();
    s.git_commit = runCommandOneLine("git rev-parse --short HEAD 2>/dev/null");
    s.executable = executable_path;
    s.config_path = config_path;
    s.working_directory = currentWorkingDirectory();
    s.input_data_path = cfg.INFO_Type ? cfg.GMTI_Data_new : cfg.GMTI_Data_add;
    s.output_dir = cfg.result_add;
    s.paths.runtime_config_json = pathJoin(cfg.result_add, "runtime_config_dump.json");
    s.paths.runtime_config_txt = pathJoin(cfg.result_add, "runtime_config_dump.txt");
    s.paths.timing_metrics_csv = pathJoin(cfg.result_add, "timing_metrics.csv");
    s.paths.run_manifest_json = pathJoin(cfg.result_add, "run_manifest.json");

    if (!s.diagnostics_enabled) {
        std::cout << "[RUNTIME] runtime_mode=" << cfg.runtime_mode
                  << ", diagnostics disabled for realtime run" << std::endl;
        return;
    }

    if (!mkdirP(cfg.result_add)) {
        std::cerr << "[PARAM_DUMP][WARN] cannot create output dir: "
                  << cfg.result_add << std::endl;
    }
    writeRuntimeConfigDump(cfg, config_path, executable_path);
    writeTimingCsvUnlocked();
    writeManifest(cfg, false, -1, "run started");

    std::cout << "[PARAM_DUMP] runtime_config_dump = "
              << s.paths.runtime_config_json << std::endl;
    std::cout << "[PARAM_DUMP] pulse_len=" << cfg.pulse_len
              << ", read_pulse_num=" << cfg.read_pulse_num
              << ", fc=" << cfg.fc
              << ", Br=" << cfg.Br
              << ", fs=" << cfg.fs
              << ", Tr=" << cfg.Tr
              << ", PRF=" << cfg.PRF << std::endl;
    std::cout << "[PARAM_DUMP] scan_min=" << cfg.scan_min_deg
              << ", scan_max=" << cfg.scan_max_deg
              << ", az_count=" << cfg.az_count
              << ", pc_crop_start=" << cfg.range_crop_start
              << ", pc_crop_len=" << cfg.range_compress_len << std::endl;
    std::cout << "[MANIFEST] run_manifest = "
              << s.paths.run_manifest_json << std::endl;
    std::cout << "[TIMING] timing_metrics = "
              << s.paths.timing_metrics_csv << std::endl;
}

void finishRun(const Config& cfg, bool normal_exit, int exit_code,
               const std::string& notes)
{
    RunState& s = state();
    if (!s.initialized || !s.diagnostics_enabled) return;
    s.end_time = currentIsoTime();
    writeTimingCsvUnlocked();
    writeManifest(cfg, normal_exit, exit_code, notes);
}

void writeRuntimeConfigDump(const Config& cfg,
                            const std::string&,
                            const std::string&)
{
    if (!diagnosticsEnabled()) return;
    writeRuntimeConfigJson(cfg);
    writeRuntimeConfigTxt(cfg);
}

void recordTiming(const char* scope_name,
                  std::chrono::system_clock::time_point start_time,
                  std::chrono::system_clock::time_point end_time,
                  long long elapsed_ms,
                  int period_id,
                  const std::string& extra)
{
    std::lock_guard<std::mutex> lock(stateMutex());
    RunState& s = state();
    if (!s.initialized || !s.diagnostics_enabled) return;

    auto toIso = [](std::chrono::system_clock::time_point tp) {
        const std::time_t t = std::chrono::system_clock::to_time_t(tp);
        std::tm tm{};
        localtime_r(&t, &tm);
        const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(
            tp.time_since_epoch()).count() % 1000000;
        std::ostringstream os;
        os << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S") << "."
           << std::setw(6) << std::setfill('0') << micros;
        return os.str();
    };

    TimingMetric m;
    m.case_id = s.case_id;
    m.run_id = s.run_id;
    m.result_id = s.result_id;
    m.period_id = period_id;
    m.scope_name = scope_name ? scope_name : "";
    m.start_time = toIso(start_time);
    m.end_time = toIso(end_time);
    m.elapsed_ms = elapsed_ms;
    m.extra = extra;
    s.timings.push_back(m);
    writeTimingCsvUnlocked();
}

void flushTimingMetrics()
{
    std::lock_guard<std::mutex> lock(stateMutex());
    if (!state().diagnostics_enabled) return;
    writeTimingCsvUnlocked();
}

bool diagnosticsEnabled()
{
    std::lock_guard<std::mutex> lock(stateMutex());
    return state().initialized && state().diagnostics_enabled;
}

const std::string& runId()
{
    return state().run_id;
}

const std::string& caseId()
{
    return state().case_id;
}

const std::string& resultId()
{
    return state().result_id;
}

RunPaths paths()
{
    return state().paths;
}

TimingScope::TimingScope(const char* name, int period_id, const std::string& extra)
    : name_(name),
      period_id_(period_id),
      extra_(extra),
      steady_start_(std::chrono::high_resolution_clock::now()),
      wall_start_(std::chrono::system_clock::now())
{
}

TimingScope::~TimingScope()
{
    const auto steady_end = std::chrono::high_resolution_clock::now();
    const auto wall_end = std::chrono::system_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        steady_end - steady_start_);
    if (diagnosticsEnabled()) {
        std::cout << "[TIMING-SCOPE] " << name_ << ": "
                  << duration.count() << " ms" << std::endl;
        recordTiming(name_, wall_start_, wall_end, duration.count(), period_id_, extra_);
    }
}

} // namespace runtime
} // namespace gmti
