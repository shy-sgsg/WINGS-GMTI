#include "target_common.h"

#include <cerrno>
#include <sys/stat.h>
#include <sys/types.h>

namespace gmti {
namespace target_injection {

double wrapTo180(double deg)
{
    while (deg > 180.0) deg -= 360.0;
    while (deg <= -180.0) deg += 360.0;
    return deg;
}

bool ensureDir(const std::string &path)
{
    if (path.empty()) return false;
    if (::mkdir(path.c_str(), 0755) == 0) return true;
    return errno == EEXIST;
}

std::string joinPath(const std::string &a, const std::string &b)
{
    if (a.empty()) return b;
    if (a[a.size() - 1] == '/') return a + b;
    return a + "/" + b;
}

} // namespace target_injection
} // namespace gmti

