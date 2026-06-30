#include "radar_geometry.h"

#include <algorithm>

namespace gmti {
namespace target_injection {

double pulseTimeSec(const RadarConfig &radar, int period_id, int beam_id, int pulse_id)
{
    const double t_beam = static_cast<double>(radar.pulse_num) / radar.prf_hz;
    const double t_cycle = static_cast<double>(radar.beam_count) * t_beam;
    return static_cast<double>(period_id) * t_cycle +
           static_cast<double>(beam_id) * t_beam +
           static_cast<double>(pulse_id) / radar.prf_hz;
}

PlatformState evaluatePlatformState(const TargetGlobalConfig &global, double time_sec)
{
    PlatformState p;
    p.position = Vec3{global.platform_speed_mps * time_sec, 0.0, global.platform_height_m};
    p.velocity = Vec3{global.platform_speed_mps, 0.0, 0.0};
    return p;
}

GeometrySample evaluateGeometry(const RadarConfig &radar,
                                const TargetGlobalConfig &global,
                                const TargetConfig &target,
                                int period_id,
                                int beam_id,
                                int pulse_id)
{
    GeometrySample g;
    g.time_sec = pulseTimeSec(radar, period_id, beam_id, pulse_id);
    g.theta_cmd_deg = radar.scan_min_deg + radar.scan_step_deg * static_cast<double>(beam_id);
    g.theta_true_deg = g.theta_cmd_deg;
    g.platform = evaluatePlatformState(global, g.time_sec);
    g.target = evaluateTargetState(target, g.time_sec);
    g.los = g.target.position - g.platform.position;
    g.range_m = norm(g.los);
    if (g.range_m > 0.0) {
        g.los_unit = g.los * (1.0 / g.range_m);
    }
    g.tau_abs_sec = 2.0 * g.range_m / kC;
    g.tau_rel_sec = g.tau_abs_sec - radar.sample_delay_sec;
    g.range_sample_float = g.tau_rel_sec * radar.fs_hz;
    g.range_sample_int = static_cast<int>(std::floor(g.range_sample_float + 0.5));
    g.in_range_window = (g.range_sample_float >= 0.0 &&
                         g.range_sample_float < static_cast<double>(radar.pulse_len));
    g.target_azimuth_deg = rad2deg(std::atan2(g.los.y, g.los.x));
    g.angle_error_deg = wrapTo180(g.target_azimuth_deg - g.theta_true_deg);
    const Vec3 rel_v = g.target.velocity - g.platform.velocity;
    g.radial_velocity_mps = dot(rel_v, g.los_unit);
    return g;
}

} // namespace target_injection
} // namespace gmti

