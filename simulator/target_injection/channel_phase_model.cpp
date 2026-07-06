#include "channel_phase_model.h"

namespace gmti {
namespace target_injection {

double channelPhaseRad(const RadarConfig &radar,
                       const TargetGlobalConfig &global,
                       const GeometrySample &geom)
{
    (void)global;
    const double lambda = kC / radar.fc_hz;
    const double look_angle_rad = deg2rad(geom.theta_true_deg);
    const double phi_static =
        2.0 * kPi * radar.d_chan_m * std::sin(look_angle_rad) / lambda;
    const double platform_speed = norm(geom.platform.velocity);
    const double delta_t = (platform_speed > 1.0e-9)
        ? (radar.d_chan_m / platform_speed)
        : 0.0;
    const double target_radial_velocity = dot(geom.target.velocity, geom.los_unit);
    const double phi_motion =
        4.0 * kPi * target_radial_velocity * delta_t / lambda;
    return phi_static + phi_motion;
}

} // namespace target_injection
} // namespace gmti
