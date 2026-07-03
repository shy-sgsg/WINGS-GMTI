#pragma once

#include "stage2_config.h"
#include "stage2_scatterer.h"

#include <random>

namespace gmti {
namespace stage2 {

Scatterer makeScattererFromRangeAzimuth(const Stage2Config &cfg,
                                        int id,
                                        const std::string &type,
                                        double slant_range_m,
                                        double azimuth_deg,
                                        double z_m,
                                        double amplitude,
                                        double phase_rad,
                                        double rcs_db);

ScattererList generateScene(const Stage2Config &cfg,
                            const Stage2RunOptions &opt,
                            std::mt19937 &rng);

bool writeSceneTruth(const std::string &truth_dir,
                     const Stage2Config &cfg,
                     const ScattererList &scatterers,
                     std::string &err);

} // namespace stage2
} // namespace gmti
