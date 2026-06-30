#pragma once

#include "../target_injection/target_common.h"

#include <string>
#include <vector>

namespace gmti {
namespace stage2 {

struct Scatterer {
    int id = 0;
    std::string type = "area";
    gmti::target_injection::Vec3 position;
    double amplitude = 1.0;
    double phase_rad = 0.0;
    double rcs_db = 0.0;
    double initial_range_m = 0.0;
    double initial_azimuth_deg = 0.0;
};

typedef std::vector<Scatterer> ScattererList;

} // namespace stage2
} // namespace gmti

