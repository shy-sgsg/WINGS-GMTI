#pragma once

#include "../common/SimulationGeometry.h"
#include "target_common.h"

#include <limits>
#include <string>

namespace gmti {
namespace target_injection {

struct RadarConfig {
    std::string input_config;
    std::string input_data_file;
    std::string output_data_file;
    int pulse_len = 11820;
    int pulse_num = 130;
    int beam_count = 61;
    int range_fft_len = 12288;
    int range_crop_start = 3864;
    int range_crop_len = 4096;
    double scan_min_deg = -60.0;
    double scan_step_deg = 2.0;
    double beam_width_deg = 2.28;
    double fc_hz = 16.0e9;
    double br_hz = 50.0e6;
    double fs_hz = 60.0e6;
    double tr_sec = 130.0e-6;
    double prf_hz = 1300.0;
    double sample_delay_sec = 488.0e-6;
    double d_chan_m = 0.17;
};

struct TargetGlobalConfig {
    std::string coordinate_mode = "project_local";
    std::string visibility_mode = "hard_gate";
    std::string amplitude_mode = "snr_db";
    std::string channel_phase_mode = "baseline_approx";
    std::string platform_mode = "ideal_platform";
    double target_snr_db = 30.0;
    double direct_amplitude = 1.0;
    double beam_gain_threshold = 0.05;
    double rms_floor = 1.0e-6;
    int chirp_phase_sign = 1;
    int carrier_phase_sign = -1;
    double platform_speed_mps = 60.0;
    double platform_height_m = 6000.0;
    double platform_origin_lat_deg = 40.45121057;
    double platform_origin_lon_deg = 116.98377429;
    double platform_origin_alt_m = 0.0;
    double projection_ref_lon_deg = 117.0;
    gmti::sim_geometry::Stage2GeometryConfig geometry;
};

struct TargetConfig {
    int id = 1;
    std::string name = "strong_cooperative_target";
    std::string motion_model = "constant_velocity";
    std::string init_mode = "range_azimuth";
    bool enabled = true;
    int start_period = 0;
    int end_period = 10;
    Vec3 p0;
    Vec3 v;
    double slant_range_m = 85000.0;
    double azimuth_deg = 0.0;
    double height_m = 0.0;
    double radial_velocity_mps = 10.0;
    double cross_velocity_mps = 0.0;
    bool has_ref_geometry = false;
    int ref_beam_id = -1;
    int ref_pulse_idx = -1;
    double ref_time_s = 0.0;
    Vec3 ref_platform;
    Vec3 ref_target;
    double ref_range_m = 0.0;
    double ref_range_sample_float = 0.0;
    int ref_range_sample_int = 0;
    double echo_delay_sample_center_used = 0.0;
    double override_speed_mps = std::numeric_limits<double>::quiet_NaN();
    double override_rcs_db = std::numeric_limits<double>::quiet_NaN();
    std::string velocity_mode;
    double target_ve_mps = std::numeric_limits<double>::quiet_NaN();
    double target_vn_mps = std::numeric_limits<double>::quiet_NaN();
    double target_vr_self_mps = std::numeric_limits<double>::quiet_NaN();
    double target_vt_self_mps = std::numeric_limits<double>::quiet_NaN();
    double af_motion_truth_hz = std::numeric_limits<double>::quiet_NaN();
};

struct RunOptions {
    std::string input_config = "outputs/stage1/config/temp_config_stage1_newsystem.xml";
    std::string input_data_dir = "outputs/stage1/data";
    std::string target_config = "targets.json";
    std::string output_dir = "outputs/stage1_target";
    std::string background_mode = "copy";
    std::string visibility_mode;
    std::string amplitude_mode;
    int period_start = 0;
    int period_count = 1;
    int beam_start = 0;
    int beam_count = -1;
    int pulse_start = 0;
    int pulse_count = -1;
    double target_snr_db = -999.0;
    double direct_amplitude = -1.0;
    bool write_truth = true;
    bool validate = true;
    bool write_config = true;
};

struct InjectionConfig {
    RadarConfig radar;
    TargetGlobalConfig global;
    TargetConfig target;
    RunOptions run;
};

bool parseCommandLine(int argc, char **argv, RunOptions &opt);
bool loadRadarConfig(const std::string &xml_path, RadarConfig &cfg, std::string &err);
bool loadTargetConfig(const std::string &json_path,
                      const RadarConfig &radar,
                      TargetGlobalConfig &global,
                      TargetConfig &target,
                      std::string &err);
void applyRunOverrides(const RunOptions &run, TargetGlobalConfig &global);
bool writeOutputConfig(const std::string &input_xml,
                       const std::string &out_xml,
                       const std::string &out_data_file,
                       std::string &err);
bool writeDefaultTargetsJson(const std::string &path, std::string &err);

} // namespace target_injection
} // namespace gmti
