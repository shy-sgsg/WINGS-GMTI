#pragma once

#include <string>

namespace gmti {
namespace sim_common {

bool writeTextFile(const std::string &path, const std::string &content, std::string &err);

} // namespace sim_common
} // namespace gmti
