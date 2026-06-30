#pragma once

namespace gmti {
namespace sim_common {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

double norm(const Vec3 &v);
double rad2deg(double x);
double deg2rad(double x);
double wrapTo180(double deg);
double azimuthDeg(const Vec3 &from, const Vec3 &to);

} // namespace sim_common
} // namespace gmti
