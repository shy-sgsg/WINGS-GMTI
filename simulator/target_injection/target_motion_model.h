#pragma once

#include "target_config.h"

namespace gmti {
namespace target_injection {

struct TargetState {
    Vec3 position;
    Vec3 velocity;
};

TargetState evaluateTargetState(const TargetConfig &target, double time_sec);

} // namespace target_injection
} // namespace gmti

