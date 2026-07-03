#pragma once

#include <string>

namespace gmti {
namespace sim_geometry {

struct LocalPoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    LocalPoint() {}
    LocalPoint(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

struct ENUPoint {
    double e = 0.0;
    double n = 0.0;
    double u = 0.0;
    ENUPoint() {}
    ENUPoint(double e_, double n_, double u_) : e(e_), n(n_), u(u_) {}
};

struct GeoPoint {
    double e = 0.0;
    double n = 0.0;
    double h = 0.0;
    double lat = 0.0;
    double lon = 0.0;
};

struct LookVectorEN {
    double east = 0.0;
    double north = 0.0;
};

struct LocalVelocity {
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    LocalVelocity() {}
    LocalVelocity(double vx_, double vy_, double vz_) : vx(vx_), vy(vy_), vz(vz_) {}
};

struct ENUVelocity {
    double ve = 0.0;
    double vn = 0.0;
    double vu = 0.0;
    ENUVelocity() {}
    ENUVelocity(double ve_, double vn_, double vu_) : ve(ve_), vn(vn_), vu(vu_) {}
};

struct PosSample {
    double lat_deg = 0.0;
    double lon_deg = 0.0;
    double height_m = 0.0;
    double vn_mps = 0.0;
    double ve_mps = 0.0;
    double vd_mps = 0.0;
};

struct Stage2GeometryConfig {
    std::string geometry_config_name = "algorithm_axis_x_north";
    std::string local_x_axis = "north";
    std::string local_y_axis = "east";
    std::string platform_heading_source = "velocity";
    double platform_heading_deg = 0.0;
    std::string beam_angle_reference = "algorithm";
    std::string beam_zero_direction = "algorithm";
    std::string beam_positive_direction = "algorithm";
    double beam_theta_offset_deg = 0.0;
    std::string range_geometry = "algorithm";
    bool use_ground_range_for_position = true;
    double platform_origin_lat_deg = 40.4512107203;
    double platform_origin_lon_deg = 116.985931582;
    double platform_origin_alt_m = 6000.0;
    double projection_ref_lon_deg = 117.0;
    int squint_side = 1;
};

ENUPoint localToEnu(const LocalPoint &p, const Stage2GeometryConfig &cfg);
LocalPoint enuToLocal(const ENUPoint &p, const Stage2GeometryConfig &cfg);
ENUVelocity localVelocityToEnu(const LocalVelocity &v, const Stage2GeometryConfig &cfg);
LocalVelocity enuVelocityToLocal(const ENUVelocity &v, const Stage2GeometryConfig &cfg);

LookVectorEN makeAlgorithmLookVectorEN(double platform_ve,
                                       double platform_vn,
                                       double theta_deg,
                                       const Stage2GeometryConfig &cfg);

double slantRangeToGroundRange(double slant_range_m,
                               double platform_height_m,
                               double target_height_m,
                               const Stage2GeometryConfig &cfg);

GeoPoint enuToLatLon(double e,
                     double n,
                     double h,
                     const Stage2GeometryConfig &cfg);

ENUPoint latLonToEnu(double lat_deg,
                     double lon_deg,
                     double h,
                     const Stage2GeometryConfig &cfg);

PosSample makeProtocolPosSample(const LocalPoint &position_local,
                                const LocalVelocity &velocity_local,
                                const Stage2GeometryConfig &cfg);

LocalPoint makePointFromRangeAzimuth(const LocalPoint &ref_platform_local,
                                     const LocalVelocity &ref_velocity_local,
                                     double slant_or_ground_range_m,
                                     double theta_deg,
                                     double platform_height_m,
                                     double target_height_m,
                                     const Stage2GeometryConfig &cfg);

} // namespace sim_geometry
} // namespace gmti
