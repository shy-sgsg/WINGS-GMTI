#pragma once

#include "radar_geometry.h"

namespace gmti {
namespace target_injection {

struct VisibilityResult {
    double beam_gain = 0.0;
    bool visible = false;
};

VisibilityResult evaluateVisibility(const RadarConfig &radar,
                                    const TargetGlobalConfig &global,
                                    double angle_error_deg);

} // namespace target_injection
} // namespace gmti

