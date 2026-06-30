#pragma once

#include "sar_scatterer_extractor.h"
#include "stage3_config.h"

#include <string>
#include <vector>

namespace gmti {
namespace stage3 {

bool forwardSarSceneToDdc(const Stage3Config &cfg,
                          const std::vector<SarScatterer> &scatterers,
                          const std::string &output_path,
                          std::string &err);

} // namespace stage3
} // namespace gmti
