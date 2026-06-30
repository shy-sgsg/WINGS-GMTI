#pragma once

#include "stage2_config.h"
#include "stage2_lfm_forward.h"
#include "stage2_scatterer.h"

namespace gmti {
namespace stage2 {

bool writeStage2Report(const std::string &path,
                       const Stage2Config &cfg,
                       const Stage2RunOptions &opt,
                       const ScattererList &scatterers,
                       const Stage2Stats &stats,
                       double elapsed_sec,
                       std::string &err);

} // namespace stage2
} // namespace gmti

