#pragma once

#include "stage3_config.h"

#include <string>

namespace gmti {
namespace stage3 {

bool writeStage3PlanningReport(const std::string &path,
                               const Stage3RunOptions &opt,
                               std::string &err);

} // namespace stage3
} // namespace gmti
