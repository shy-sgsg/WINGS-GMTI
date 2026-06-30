#pragma once

#include "sar_georef.h"
#include "sar_image_reader.h"
#include "stage3_config.h"

#include <string>
#include <vector>

namespace gmti {
namespace stage3 {

struct SarScatterer {
    int scatterer_id = 0;
    int source_row = 0;
    int source_col = 0;
    GroundPoint pos;
    double sar_intensity_raw = 0.0;
    double sar_intensity_norm = 0.0;
    double amplitude = 0.0;
    double phase_rad = 0.0;
    std::string scatterer_type = "area";
    double initial_slant_range_m = 0.0;
    double initial_azimuth_deg = 0.0;
    std::string selected_by = "grid_adaptive_sampling";
    int grid_id = 0;
};

bool extractSarScatterers(const SarImage &normalized,
                          const ImageCornerGeoref &georef,
                          const Stage3Config &cfg,
                          std::vector<SarScatterer> &scatterers,
                          std::string &err);

} // namespace stage3
} // namespace gmti
