#pragma once

#include <cmath>
#include <complex>
#include <cstdint>
#include <string>

namespace gmti {
namespace target_injection {

const double kPi = 3.14159265358979323846;
const double kC = 299792458.0;

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() {}
    Vec3(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in) {}
};

inline Vec3 operator+(const Vec3 &a, const Vec3 &b)
{
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vec3 operator-(const Vec3 &a, const Vec3 &b)
{
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vec3 operator*(const Vec3 &a, double s)
{
    return Vec3{a.x * s, a.y * s, a.z * s};
}

inline double dot(const Vec3 &a, const Vec3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(const Vec3 &a)
{
    return std::sqrt(dot(a, a));
}

inline double deg2rad(double x)
{
    return x * kPi / 180.0;
}

inline double rad2deg(double x)
{
    return x * 180.0 / kPi;
}

double wrapTo180(double deg);
bool ensureDir(const std::string &path);
std::string joinPath(const std::string &a, const std::string &b);

} // namespace target_injection
} // namespace gmti
