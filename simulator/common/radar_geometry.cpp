#include "radar_geometry.h"

#include <cmath>

namespace gmti {
namespace sim_common {

namespace {
const double kPi = 3.14159265358979323846;
}

double norm(const Vec3 &v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double rad2deg(double x)
{
    return x * 180.0 / kPi;
}

double deg2rad(double x)
{
    return x * kPi / 180.0;
}

double wrapTo180(double deg)
{
    while (deg > 180.0) deg -= 360.0;
    while (deg <= -180.0) deg += 360.0;
    return deg;
}

double azimuthDeg(const Vec3 &from, const Vec3 &to)
{
    return rad2deg(std::atan2(to.y - from.y, to.x - from.x));
}

} // namespace sim_common
} // namespace gmti
