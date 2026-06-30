#pragma once

#include "sar_scatterer_extractor.h"

#include <string>
#include <vector>

namespace gmti {
namespace stage3 {

bool limitScattererCount(int max_scatterers,
                         std::vector<SarScatterer> &scatterers,
                         std::string &err);

} // namespace stage3
} // namespace gmti
