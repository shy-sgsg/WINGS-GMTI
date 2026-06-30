#pragma once

#include <string>

namespace gmti {
namespace stage3 {

struct Stage3SystemConfig {
    double fc_ghz = 16.0;
    double bandwidth_mhz = 50.0;
    double fs_mhz = 60.0;
    double pulse_width_us = 130.0;
    double prf_hz = 1300.0;
    double sample_delay_us = 488.0;
    double sample_window_us = 197.0;
    int ddc_len = 11820;
    int fft_len = 12288;
    int pc_crop_start = 3864;
    int pc_crop_len = 4096;
    double scan_min_deg = -60.0;
    double scan_step_deg = 2.0;
    int beam_count = 61;
    double beam_width_deg = 2.28;
    int pulse_num = 130;
    double platform_height_m = 6000.0;
    double platform_speed_mps = 60.0;
    double d_chan_m = 0.17;
};

struct Stage3SarInputConfig {
    std::string image_path = "data/sar_highres.tif";
    std::string input_type = "geotiff";
    std::string georef_path;
    std::string image_value_type = "intensity";
    double nodata_value = 0.0;
    bool crop_enabled = true;
    double crop_range_min_m = 82800.0;
    double crop_range_max_m = 93000.0;
    double crop_azimuth_min_deg = -60.0;
    double crop_azimuth_max_deg = 60.0;
};

struct Stage3ScattererExtractionConfig {
    std::string method = "grid_adaptive_sampling";
    int max_scatterers = 50000;
    int min_scatterers = 2000;
    double intensity_threshold_percentile = 60.0;
    double strong_threshold_percentile = 99.5;
    double grid_cell_m = 10.0;
    int max_scatterers_per_cell = 3;
    double amplitude_scale = 1.0;
    double amplitude_gamma = 0.5;
    std::string phase_mode = "random_uniform";
    unsigned int random_seed = 202606;
};

struct Stage3ForwardConfig {
    std::string platform_mode = "ideal_straight";
    std::string beam_pattern_mode = "gaussian";
    double beam_gain_threshold = 0.01;
    std::string channel_phase_mode = "baseline_approx";
    int chirp_phase_sign = 1;
    int carrier_phase_sign = -1;
    int period_start = 0;
    int period_count = 1;
    bool enable_openmp = true;
};

struct Stage3TargetsConfig {
    bool enabled = true;
    std::string target_config_path = "targets.json";
    double target_snr_db = 25.0;
    std::string amplitude_mode = "snr_db";
};

struct Stage3Config {
    Stage3SystemConfig system;
    Stage3SarInputConfig sar_input;
    Stage3ScattererExtractionConfig extraction;
    Stage3ForwardConfig forward;
    Stage3TargetsConfig targets;
};

struct Stage3RunOptions {
    std::string stage3_config = "stage3_config.json";
    std::string sar_image;
    std::string georef_path;
    std::string scatterer_csv;
    std::string output_dir = "outputs/stage3";
    std::string target_config;
    std::string scene_mode = "sar";
    bool target_enabled = true;
    bool validate = true;
    bool extract_only = false;
    bool forward_only = false;
    int period_start = -1;
    int period_count = -1;
    double single_scatterer_range_m = 85000.0;
    double single_scatterer_azimuth_deg = 0.0;
    double single_scatterer_amplitude = 1.0;
};

bool parseStage3CommandLine(int argc, char **argv, Stage3RunOptions &opt, std::string &err);
bool ensureStage3Dirs(const std::string &output_dir, std::string &err);
bool writeDefaultStage3Config(const std::string &path, std::string &err);
std::string pathJoin(const std::string &a, const std::string &b);

} // namespace stage3
} // namespace gmti
