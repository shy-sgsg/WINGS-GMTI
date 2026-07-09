#pragma once

#include "lfm_echo_generator.h"

#include <fstream>
#include <limits>
#include <map>

namespace gmti {
namespace target_injection {

struct BeamSummaryAccumulator {
    int period_id = 0;
    int beam_id = 0;
    int target_id = 0;
    std::string target_name;
    int rows = 0;
    int visible_pulse_count = 0;
    double sum_range_m = 0.0;
    double sum_range_sample = 0.0;
    double sum_theta_cmd_deg = 0.0;
    double sum_target_azimuth_deg = 0.0;
    double sum_beam_gain = 0.0;
    double sum_radial_velocity_mps = 0.0;
    double sum_target_amplitude = 0.0;
    int injected_sample_count = 0;
    bool has_ref_geometry = false;
    int ref_pulse_idx = -1;
    double ref_time_s = 0.0;
    Vec3 ref_platform;
    Vec3 ref_target;
    double ref_platform_e = 0.0;
    double ref_platform_n = 0.0;
    double ref_platform_lat = 0.0;
    double ref_platform_lon = 0.0;
    double ref_platform_ve = 0.0;
    double ref_platform_vn = 0.0;
    double ref_target_e = 0.0;
    double ref_target_n = 0.0;
    double ref_target_lat = 0.0;
    double ref_target_lon = 0.0;
    double ref_range_m = 0.0;
    double ref_range_sample_float = 0.0;
    int expected_range_bin = -1;
    double echo_delay_sample_center_used = 0.0;
    double look_e = 0.0;
    double look_n = 0.0;
    double ground_range_m = 0.0;
    double slant_range_m = 0.0;
    std::string geometry_config_name;
    double moving_target_speed_mps = std::numeric_limits<double>::quiet_NaN();
    double rcs_db = std::numeric_limits<double>::quiet_NaN();
    double snr_db = std::numeric_limits<double>::quiet_NaN();
    double target_ve_mps = std::numeric_limits<double>::quiet_NaN();
    double target_vn_mps = std::numeric_limits<double>::quiet_NaN();
    double target_vr_self_mps = std::numeric_limits<double>::quiet_NaN();
    double target_vt_self_mps = std::numeric_limits<double>::quiet_NaN();
    double af_motion_truth_hz = std::numeric_limits<double>::quiet_NaN();
    double af_geometry_truth_hz = std::numeric_limits<double>::quiet_NaN();
    double af_total_truth_hz = std::numeric_limits<double>::quiet_NaN();
    double phi_total_truth_rad = std::numeric_limits<double>::quiet_NaN();
    double phi_static_truth_rad = std::numeric_limits<double>::quiet_NaN();
    double phi_motion_truth_rad = std::numeric_limits<double>::quiet_NaN();
    int row_truth = -1;
};

class TruthWriter {
public:
    void setCaseId(const std::string &case_id) { case_id_ = case_id; }
    bool open(const std::string &truth_dir, std::string &err);
    void writePulse(const PulseTruth &t);
    void writeSummary();
    void close();

private:
    std::ofstream pulse_;
    std::ofstream summary_;
    std::ofstream moving_;
    std::string case_id_ = "stage2";
    std::map<std::string, BeamSummaryAccumulator> acc_;
};

bool writeTargetInjectionReport(const std::string &path,
                                const InjectionConfig &cfg,
                                const InjectionStats &stats,
                                const std::string &notes,
                                std::string &err);

} // namespace target_injection
} // namespace gmti
