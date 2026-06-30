#pragma once

#include "target_config.h"
#include "target_motion_model.h"

namespace gmti {
namespace target_injection {

struct PlatformState {
    Vec3 position;
    Vec3 velocity;
};

struct GeometrySample {
    double time_sec = 0.0;
    double theta_cmd_deg = 0.0;
    double theta_true_deg = 0.0;
    PlatformState platform;
    TargetState target;
    Vec3 los;
    Vec3 los_unit;
    double range_m = 0.0;
    double tau_abs_sec = 0.0;
    double tau_rel_sec = 0.0;
    double range_sample_float = 0.0;
    int range_sample_int = 0;
    bool in_range_window = false;
    double target_azimuth_deg = 0.0;
    double angle_error_deg = 0.0;
    double radial_velocity_mps = 0.0;
};

double pulseTimeSec(const RadarConfig &radar, int period_id, int beam_id, int pulse_id);
PlatformState evaluatePlatformState(const TargetGlobalConfig &global, double time_sec);
GeometrySample evaluateGeometry(const RadarConfig &radar,
                                const TargetGlobalConfig &global,
                                const TargetConfig &target,
                                int period_id,
                                int beam_id,
                                int pulse_id);

} // namespace target_injection
} // namespace gmti

