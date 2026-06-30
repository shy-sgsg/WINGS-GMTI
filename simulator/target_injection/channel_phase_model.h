#pragma once

#include "radar_geometry.h"

namespace gmti {
namespace target_injection {

double channelPhaseRad(const RadarConfig &radar,
                       const TargetGlobalConfig &global,
                       const GeometrySample &geom);

} // namespace target_injection
} // namespace gmti

