#include "target_motion_model.h"

namespace gmti {
namespace target_injection {

TargetState evaluateTargetState(const TargetConfig &target, double time_sec)
{
    TargetState s;
    s.position = target.p0 + target.v * time_sec;
    s.velocity = target.v;
    return s;
}

} // namespace target_injection
} // namespace gmti

