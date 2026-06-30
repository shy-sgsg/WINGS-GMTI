#include "beam_visibility.h"

#include <cmath>

namespace gmti {
namespace target_injection {

VisibilityResult evaluateVisibility(const RadarConfig &radar,
                                    const TargetGlobalConfig &global,
                                    double angle_error_deg)
{
    VisibilityResult r;
    const double ae = std::fabs(angle_error_deg);
    if (global.visibility_mode == "gaussian") {
        const double sigma = radar.beam_width_deg / 2.355;
        r.beam_gain = std::exp(-0.5 * (angle_error_deg / sigma) * (angle_error_deg / sigma));
        r.visible = r.beam_gain >= global.beam_gain_threshold;
    } else {
        r.visible = ae <= radar.beam_width_deg * 0.5;
        r.beam_gain = r.visible ? 1.0 : 0.0;
    }
    return r;
}

} // namespace target_injection
} // namespace gmti

