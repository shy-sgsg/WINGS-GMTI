#include "truth_writer.h"

#include <fstream>

namespace gmti {
namespace sim_common {

bool writeTextFile(const std::string &path, const std::string &content, std::string &err)
{
    std::ofstream os(path.c_str());
    if (!os) {
        err = "failed to open " + path;
        return false;
    }
    os << content;
    return true;
}

} // namespace sim_common
} // namespace gmti
