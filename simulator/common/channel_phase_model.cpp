#include "channel_phase_model.h"
#include "radar_geometry.h"

#include <cmath>

namespace gmti {
namespace sim_common {

namespace {
const double kPi = 3.14159265358979323846;
}

double baselineApproxPhaseRad(double d_chan_m, double theta_true_deg, double wavelength_m)
{
    if (wavelength_m == 0.0) return 0.0;
    return 2.0 * kPi * d_chan_m * std::sin(deg2rad(theta_true_deg)) / wavelength_m;
}

} // namespace sim_common
} // namespace gmti
