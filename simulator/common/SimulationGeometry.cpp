#include "SimulationGeometry.h"

#include "geo/geoProj.hpp"

#include <algorithm>
#include <cmath>

namespace gmti {
namespace sim_geometry {

namespace {

const double kPi = 3.14159265358979323846;

double deg2rad(double v)
{
    return v * kPi / 180.0;
}

void addLocalAxis(double local_value,
                  const std::string &axis,
                  double &east,
                  double &north)
{
    if (axis == "east" || axis == "track_right") {
        east += local_value;
    } else if (axis == "west" || axis == "track_left") {
        east -= local_value;
    } else if (axis == "north" || axis == "track_forward") {
        north += local_value;
    } else if (axis == "south" || axis == "track_backward") {
        north -= local_value;
    }
}

void axisCoefficients(const Stage2GeometryConfig &cfg,
                      double &xe,
                      double &xn,
                      double &ye,
                      double &yn)
{
    xe = xn = ye = yn = 0.0;
    addLocalAxis(1.0, cfg.local_x_axis, xe, xn);
    addLocalAxis(1.0, cfg.local_y_axis, ye, yn);
}

} // namespace

ENUPoint localToEnu(const LocalPoint &p, const Stage2GeometryConfig &cfg)
{
    double xe = 0.0, xn = 0.0, ye = 0.0, yn = 0.0;
    axisCoefficients(cfg, xe, xn, ye, yn);
    ENUPoint out;
    out.e = xe * p.x + ye * p.y;
    out.n = xn * p.x + yn * p.y;
    out.u = p.z;
    return out;
}

LocalPoint enuToLocal(const ENUPoint &p, const Stage2GeometryConfig &cfg)
{
    double xe = 0.0, xn = 0.0, ye = 0.0, yn = 0.0;
    axisCoefficients(cfg, xe, xn, ye, yn);
    const double det = xe * yn - ye * xn;
    LocalPoint out;
    if (std::abs(det) > 1.0e-12) {
        out.x = (p.e * yn - ye * p.n) / det;
        out.y = (xe * p.n - p.e * xn) / det;
    }
    out.z = p.u;
    return out;
}

ENUVelocity localVelocityToEnu(const LocalVelocity &v, const Stage2GeometryConfig &cfg)
{
    const ENUPoint p = localToEnu(LocalPoint{v.vx, v.vy, v.vz}, cfg);
    ENUVelocity out;
    out.ve = p.e;
    out.vn = p.n;
    out.vu = p.u;
    return out;
}

LocalVelocity enuVelocityToLocal(const ENUVelocity &v, const Stage2GeometryConfig &cfg)
{
    const LocalPoint p = enuToLocal(ENUPoint{v.ve, v.vn, v.vu}, cfg);
    LocalVelocity out;
    out.vx = p.x;
    out.vy = p.y;
    out.vz = p.z;
    return out;
}

LookVectorEN makeAlgorithmLookVectorEN(double platform_ve,
                                       double platform_vn,
                                       double theta_deg,
                                       const Stage2GeometryConfig &cfg)
{
    const double heading_deg =
        (cfg.platform_heading_source == "fixed_angle")
            ? cfg.platform_heading_deg
            : std::atan2(platform_vn, platform_ve) * 180.0 / kPi;
    const double side_dir = (cfg.squint_side == 1) ? -90.0 : 90.0;
    const double beam_center_dir = side_dir - (theta_deg + cfg.beam_theta_offset_deg);
    const double target_azimuth = heading_deg - beam_center_dir;
    double east = std::cos(deg2rad(target_azimuth));
    double north = std::sin(deg2rad(target_azimuth));

    const double norm = std::sqrt(east * east + north * north);
    LookVectorEN out;
    if (norm > 1.0e-12) {
        out.east = east / norm;
        out.north = north / norm;
    }
    return out;
}

LookVectorEN computeLookFromSinA(double sinA,
                                 double platform_ve,
                                 double platform_vn,
                                 int look_side)
{
    LookVectorEN out;
    if (!std::isfinite(sinA) ||
        !std::isfinite(platform_ve) ||
        !std::isfinite(platform_vn)) {
        out.east = std::numeric_limits<double>::quiet_NaN();
        out.north = std::numeric_limits<double>::quiet_NaN();
        return out;
    }

    const double vnorm = std::sqrt(platform_ve * platform_ve + platform_vn * platform_vn);
    if (!(vnorm > 1.0e-12)) {
        out.east = std::numeric_limits<double>::quiet_NaN();
        out.north = std::numeric_limits<double>::quiet_NaN();
        return out;
    }

    const double along_e = platform_ve / vnorm;
    const double along_n = platform_vn / vnorm;
    const double left_e = -along_n;
    const double left_n = along_e;
    const double right_e = along_n;
    const double right_n = -along_e;
    const double cross = std::sqrt(std::max(0.0, 1.0 - sinA * sinA));
    // Keep the mapping consistent with makeAlgorithmLookVectorEN():
    // squint_side/look_side == 1 is left-looking, 0 is right-looking.
    const bool use_left = (look_side == 1);
    const double side_e = use_left ? left_e : right_e;
    const double side_n = use_left ? left_n : right_n;

    out.east = cross * side_e + sinA * along_e;
    out.north = cross * side_n + sinA * along_n;

    const double norm = std::sqrt(out.east * out.east + out.north * out.north);
    if (norm > 1.0e-12) {
        out.east /= norm;
        out.north /= norm;
    } else {
        out.east = std::numeric_limits<double>::quiet_NaN();
        out.north = std::numeric_limits<double>::quiet_NaN();
    }
    return out;
}

double slantRangeToGroundRange(double slant_range_m,
                               double platform_height_m,
                               double target_height_m,
                               const Stage2GeometryConfig &cfg)
{
    if (!cfg.use_ground_range_for_position || cfg.range_geometry == "slant") {
        return slant_range_m;
    }
    const double dz = platform_height_m - target_height_m;
    return std::sqrt(std::max(0.0, slant_range_m * slant_range_m - dz * dz));
}

GeoPoint enuToLatLon(double e,
                     double n,
                     double h,
                     const Stage2GeometryConfig &cfg)
{
    double e0 = 0.0;
    double n0 = 0.0;
    Gaussp3(cfg.platform_origin_lat_deg,
            cfg.platform_origin_lon_deg,
            cfg.projection_ref_lon_deg,
            e0,
            n0);
    GeoPoint out;
    out.e = e0 + e;
    out.n = n0 + n;
    out.h = h;
    (void)Gaussp3RV(out.e, out.n, cfg.projection_ref_lon_deg, out.lat, out.lon);
    return out;
}

ENUPoint latLonToEnu(double lat_deg,
                     double lon_deg,
                     double h,
                     const Stage2GeometryConfig &cfg)
{
    double e0 = 0.0;
    double n0 = 0.0;
    double e = 0.0;
    double n = 0.0;
    Gaussp3(cfg.platform_origin_lat_deg,
            cfg.platform_origin_lon_deg,
            cfg.projection_ref_lon_deg,
            e0,
            n0);
    Gaussp3(lat_deg, lon_deg, cfg.projection_ref_lon_deg, e, n);
    ENUPoint out;
    out.e = e - e0;
    out.n = n - n0;
    out.u = h;
    return out;
}

PosSample makeProtocolPosSample(const LocalPoint &position_local,
                                const LocalVelocity &velocity_local,
                                const Stage2GeometryConfig &cfg)
{
    const ENUPoint enu = localToEnu(position_local, cfg);
    const GeoPoint geo = enuToLatLon(enu.e, enu.n, position_local.z, cfg);
    const ENUVelocity vel = localVelocityToEnu(velocity_local, cfg);
    PosSample out;
    out.lat_deg = geo.lat;
    out.lon_deg = geo.lon;
    out.height_m = geo.h;
    out.vn_mps = vel.vn;
    out.ve_mps = vel.ve;
    out.vd_mps = -vel.vu;
    return out;
}

LocalPoint makePointFromRangeAzimuth(const LocalPoint &ref_platform_local,
                                     const LocalVelocity &ref_velocity_local,
                                     double slant_or_ground_range_m,
                                     double theta_deg,
                                     double platform_height_m,
                                     double target_height_m,
                                     const Stage2GeometryConfig &cfg)
{
    const ENUPoint ref_enu = localToEnu(ref_platform_local, cfg);
    const ENUVelocity ref_vel = localVelocityToEnu(ref_velocity_local, cfg);
    const double ground_range =
        slantRangeToGroundRange(slant_or_ground_range_m, platform_height_m, target_height_m, cfg);
    const LookVectorEN look = makeAlgorithmLookVectorEN(ref_vel.ve, ref_vel.vn, theta_deg, cfg);
    ENUPoint target_enu;
    target_enu.e = ref_enu.e + ground_range * look.east;
    target_enu.n = ref_enu.n + ground_range * look.north;
    target_enu.u = target_height_m;
    return enuToLocal(target_enu, cfg);
}

} // namespace sim_geometry
} // namespace gmti
