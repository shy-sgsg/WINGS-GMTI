#include "beam_pattern.h"

#include <cmath>

namespace gmti {
namespace sim_common {

double gaussianBeamGain(double angle_error_deg, double beam_width_deg)
{
    const double sigma = beam_width_deg / 2.355;
    if (sigma <= 0.0) return 0.0;
    return std::exp(-0.5 * (angle_error_deg / sigma) * (angle_error_deg / sigma));
}

} // namespace sim_common
} // namespace gmti
