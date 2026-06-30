#include "stage3_config.h"

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

namespace gmti {
namespace stage3 {

namespace {

bool mkdirOne(const std::string &p, std::string &err)
{
    if (p.empty()) return true;
    if (::mkdir(p.c_str(), 0775) == 0 || errno == EEXIST) return true;
    err = "mkdir failed for " + p + ": " + std::strerror(errno);
    return false;
}

bool mkdirRecursive(const std::string &p, std::string &err)
{
    std::string cur;
    for (size_t i = 0; i < p.size(); ++i) {
        cur.push_back(p[i]);
        if (p[i] == '/' || i + 1 == p.size()) {
            if (cur == "/" || cur.empty()) continue;
            if (cur.size() > 1 && cur[cur.size() - 1] == '/') cur.resize(cur.size() - 1);
            if (!cur.empty() && !mkdirOne(cur, err)) return false;
            if (i + 1 < p.size() && (cur.empty() || cur[cur.size() - 1] != '/')) cur.push_back('/');
        }
    }
    return true;
}

bool readBool(const std::string &s)
{
    return !(s == "false" || s == "0" || s == "no");
}

} // namespace

bool parseStage3CommandLine(int argc, char **argv, Stage3RunOptions &opt, std::string &err)
{
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        const char *v = (i + 1 < argc) ? argv[i + 1] : nullptr;
        if (k == "--stage3-config" && v) opt.stage3_config = argv[++i];
        else if (k == "--sar-image" && v) opt.sar_image = argv[++i];
        else if (k == "--georef" && v) opt.georef_path = argv[++i];
        else if (k == "--scatterer-csv" && v) opt.scatterer_csv = argv[++i];
        else if (k == "--output-dir" && v) opt.output_dir = argv[++i];
        else if (k == "--target-config" && v) opt.target_config = argv[++i];
        else if (k == "--target-enabled" && v) opt.target_enabled = readBool(argv[++i]);
        else if (k == "--validate" && v) opt.validate = readBool(argv[++i]);
        else if (k == "--extract-only" && v) opt.extract_only = readBool(argv[++i]);
        else if (k == "--forward-only" && v) opt.forward_only = readBool(argv[++i]);
        else if (k == "--scene-mode" && v) opt.scene_mode = argv[++i];
        else if (k == "--period-start" && v) opt.period_start = std::atoi(argv[++i]);
        else if (k == "--period-count" && v) opt.period_count = std::atoi(argv[++i]);
        else if (k == "--single-scatterer-range" && v) opt.single_scatterer_range_m = std::atof(argv[++i]);
        else if (k == "--single-scatterer-azimuth" && v) opt.single_scatterer_azimuth_deg = std::atof(argv[++i]);
        else if (k == "--single-scatterer-amplitude" && v) opt.single_scatterer_amplitude = std::atof(argv[++i]);
        else {
            err = "unknown or incomplete option: " + k;
            return false;
        }
    }
    return true;
}

bool ensureStage3Dirs(const std::string &output_dir, std::string &err)
{
    const char *names[] = {"config", "data", "truth", "debug", "figures", "logs", "reports"};
    if (!mkdirRecursive(output_dir, err)) return false;
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i) {
        if (!mkdirRecursive(pathJoin(output_dir, names[i]), err)) return false;
    }
    return true;
}

bool writeDefaultStage3Config(const std::string &path, std::string &err)
{
    std::ofstream os(path.c_str());
    if (!os) {
        err = "failed to open " + path;
        return false;
    }
    os << "{\n"
       << "  \"system\": {\"fc_ghz\": 16.0, \"bandwidth_mhz\": 50.0, \"fs_mhz\": 60.0,\n"
       << "    \"pulse_width_us\": 130.0, \"prf_hz\": 1300.0, \"sample_delay_us\": 488.0,\n"
       << "    \"sample_window_us\": 197.0, \"ddc_len\": 11820, \"fft_len\": 12288,\n"
       << "    \"pc_crop_start\": 3864, \"pc_crop_len\": 4096, \"scan_min_deg\": -60.0,\n"
       << "    \"scan_step_deg\": 2.0, \"beam_count\": 61, \"beam_width_deg\": 2.28,\n"
       << "    \"pulse_num\": 130, \"platform_height_m\": 6000.0,\n"
       << "    \"platform_speed_mps\": 60.0, \"d_chan_m\": 0.17},\n"
       << "  \"sar_input\": {\"image_path\": \"data/sar_highres.tif\", \"input_type\": \"geotiff\",\n"
       << "    \"georef_path\": \"\", \"image_value_type\": \"intensity\", \"nodata_value\": 0,\n"
       << "    \"crop_enabled\": true, \"crop_range_min_m\": 82800.0, \"crop_range_max_m\": 93000.0,\n"
       << "    \"crop_azimuth_min_deg\": -60.0, \"crop_azimuth_max_deg\": 60.0},\n"
       << "  \"scatterer_extraction\": {\"method\": \"grid_adaptive_sampling\", \"max_scatterers\": 50000,\n"
       << "    \"min_scatterers\": 2000, \"intensity_threshold_percentile\": 60.0,\n"
       << "    \"strong_threshold_percentile\": 99.5, \"grid_cell_m\": 10.0,\n"
       << "    \"max_scatterers_per_cell\": 3, \"amplitude_scale\": 1.0,\n"
       << "    \"amplitude_gamma\": 0.5, \"phase_mode\": \"random_uniform\", \"random_seed\": 202606},\n"
       << "  \"forward\": {\"platform_mode\": \"ideal_straight\", \"beam_pattern_mode\": \"gaussian\",\n"
       << "    \"beam_gain_threshold\": 0.01, \"channel_phase_mode\": \"baseline_approx\",\n"
       << "    \"chirp_phase_sign\": 1, \"carrier_phase_sign\": -1, \"period_start\": 0,\n"
       << "    \"period_count\": 1, \"enable_openmp\": true},\n"
       << "  \"targets\": {\"enabled\": true, \"target_config_path\": \"targets.json\",\n"
       << "    \"target_snr_db\": 25.0, \"amplitude_mode\": \"snr_db\"}\n"
       << "}\n";
    return true;
}

std::string pathJoin(const std::string &a, const std::string &b)
{
    if (a.empty()) return b;
    if (a[a.size() - 1] == '/') return a + b;
    return a + "/" + b;
}

} // namespace stage3
} // namespace gmti
