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
    return 2.0 * kPi * radar.d_chan_m * std::sin(look_angle_rad) / lambda;
}

} // namespace target_injection
} // namespace gmti

