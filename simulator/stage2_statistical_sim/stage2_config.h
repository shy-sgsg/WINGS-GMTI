#pragma once

#include "../target_injection/target_config.h"

#include <string>

namespace gmti {
namespace stage2 {

struct AreaClutterConfig {
    bool enabled = true;
    std::string model = "rayleigh_lognormal_texture";
    int scatterer_count = 1000;
    double mean_power = 1.0;
    double texture_sigma = 0.4;
    double spatial_cell_m = 30.0;
};

struct StrongScattererConfig {
    bool enabled = true;
    int count = 20;
    double rcs_db_min = 10.0;
    double rcs_db_max = 30.0;
};

struct LineScattererConfig {
    bool enabled = true;
    int line_count = 1;
    int points_per_line = 50;
    double rcs_db = 12.0;
};

struct ThermalNoiseConfig {
    bool enabled = true;
    double noise_power = 0.01;
};

struct SceneConfig {
    double range_min_m = 82800.0;
    double range_max_m = 93000.0;
    double azimuth_min_deg = -60.0;
    double azimuth_max_deg = 60.0;
    double ground_z_m = 0.0;
    AreaClutterConfig area;
    StrongScattererConfig strong;
    LineScattererConfig line;
    ThermalNoiseConfig noise;
};

struct SimulationConfig {
    int period_start = 0;
    int period_count = 1;
    uint32_t random_seed = 202606;
    std::string platform_mode = "ideal_straight";
    std::string beam_pattern_mode = "gaussian";
    std::string channel_phase_mode = "baseline_approx";
    int chirp_phase_sign = 1;
    int carrier_phase_sign = -1;
    std::string output_format = "compatible_with_algorithm";
};

struct TargetStage2Config {
    bool enabled = true;
    std::string target_config_path = "targets.json";
    double target_snr_db = 25.0;
    std::string amplitude_mode = "snr_db";
};

struct Stage2Config {
    gmti::target_injection::RadarConfig radar;
    double platform_height_m = 6000.0;
    double platform_speed_mps = 60.0;
    SimulationConfig sim;
    SceneConfig scene;
    TargetStage2Config target;
};

struct Stage2RunOptions {
    std::string stage2_config = "stage2_config.json";
    std::string target_config = "";
    std::string output_dir = "outputs/stage2";
    std::string scene_mode = "full";
    bool target_enabled = true;
    bool validate = true;
    int period_start = -1;
    int period_count = -1;
    double single_scatterer_range_m = 85000.0;
    double single_scatterer_azimuth_deg = 0.0;
    double single_scatterer_amplitude = 1.0;
};

bool parseStage2CommandLine(int argc, char **argv, Stage2RunOptions &opt);
bool loadStage2Config(const std::string &path, Stage2Config &cfg, std::string &err);
bool writeDefaultStage2Config(const std::string &path, std::string &err);
bool writeStage2OutputConfig(const Stage2Config &cfg,
                             const std::string &out_xml,
                             const std::string &data_file,
                             std::string &err);
gmti::target_injection::TargetGlobalConfig makeTargetGlobal(const Stage2Config &cfg);

} // namespace stage2
} // namespace gmti
