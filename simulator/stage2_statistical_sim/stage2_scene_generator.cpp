#include "stage2_scene_generator.h"

#include "../common/SimulationGeometry.h"
#include "../target_injection/target_common.h"
#include "../target_injection/radar_geometry.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <limits>

namespace gmti {
namespace stage2 {

using gmti::target_injection::Vec3;
using gmti::target_injection::ensureDir;
using gmti::target_injection::evaluatePlatformState;
using gmti::target_injection::joinPath;
using gmti::target_injection::kC;
using gmti::target_injection::kPi;
using gmti::target_injection::norm;
using gmti::target_injection::pulseTimeSec;

namespace {

bool isTargetOnlyMode(const std::string &scene_mode)
{
    return scene_mode == "target_only" ||
           scene_mode == "cooperative_target_only" ||
           scene_mode == "empty_target";
}

double uni(std::mt19937 &rng, double a, double b)
{
    std::uniform_real_distribution<double> d(a, b);
    return d(rng);
}

double rayleigh(std::mt19937 &rng, double sigma)
{
    std::uniform_real_distribution<double> d(1.0e-12, 1.0);
    return sigma * std::sqrt(-2.0 * std::log(d(rng)));
}

gmti::sim_geometry::GeoPoint projectLocalPoint(const Stage2Config &cfg, const Vec3 &local)
{
    const gmti::sim_geometry::ENUPoint enu =
        gmti::sim_geometry::localToEnu(
            gmti::sim_geometry::LocalPoint(local.x, local.y, local.z),
            cfg.geometry);
    return gmti::sim_geometry::enuToLatLon(enu.e, enu.n, local.z, cfg.geometry);
}

void appendArea(const Stage2Config &cfg, std::mt19937 &rng, ScattererList &out, int &next_id)
{
    if (!cfg.scene.area.enabled) return;
    const double sigma = std::sqrt(std::max(1.0e-12, cfg.scene.area.mean_power) / 2.0);
    std::lognormal_distribution<double> texture(0.0, cfg.scene.area.texture_sigma);
    for (int i = 0; i < cfg.scene.area.scatterer_count; ++i) {
        const double r = uni(rng, cfg.scene.range_min_m, cfg.scene.range_max_m);
        const double az = uni(rng, cfg.scene.azimuth_min_deg, cfg.scene.azimuth_max_deg);
        double amp = rayleigh(rng, sigma);
        if (cfg.scene.area.model == "rayleigh_lognormal_texture") {
            amp *= texture(rng);
        }
        out.push_back(makeScattererFromRangeAzimuth(cfg, next_id++, "area", r, az, cfg.scene.ground_z_m,
                                                    amp, uni(rng, 0.0, 2.0 * kPi), 20.0 * std::log10(std::max(amp, 1e-12))));
    }
}

void appendStrong(const Stage2Config &cfg, std::mt19937 &rng, ScattererList &out, int &next_id)
{
    if (!cfg.scene.strong.enabled) return;
    for (int i = 0; i < cfg.scene.strong.count; ++i) {
        const double r = uni(rng, cfg.scene.range_min_m, cfg.scene.range_max_m);
        const double az = uni(rng, cfg.scene.azimuth_min_deg, cfg.scene.azimuth_max_deg);
        const double rcs = uni(rng, cfg.scene.strong.rcs_db_min, cfg.scene.strong.rcs_db_max);
        const double amp = std::pow(10.0, rcs / 20.0);
        out.push_back(makeScattererFromRangeAzimuth(cfg, next_id++, "strong", r, az, cfg.scene.ground_z_m,
                                                    amp, uni(rng, 0.0, 2.0 * kPi), rcs));
    }
}

void appendLines(const Stage2Config &cfg, std::mt19937 &rng, ScattererList &out, int &next_id)
{
    if (!cfg.scene.line.enabled) return;
    std::normal_distribution<double> jitter(0.0, 1.5);
    for (int l = 0; l < cfg.scene.line.line_count; ++l) {
        const double r0 = uni(rng, cfg.scene.range_min_m, cfg.scene.range_max_m);
        const double a0 = uni(rng, cfg.scene.azimuth_min_deg, cfg.scene.azimuth_max_deg);
        const double r1 = uni(rng, cfg.scene.range_min_m, cfg.scene.range_max_m);
        const double a1 = uni(rng, cfg.scene.azimuth_min_deg, cfg.scene.azimuth_max_deg);
        const int npts = std::max(1, cfg.scene.line.points_per_line);
        for (int i = 0; i < npts; ++i) {
            const double w = npts == 1 ? 0.0 : static_cast<double>(i) / static_cast<double>(npts - 1);
            const double r = r0 + (r1 - r0) * w;
            const double az = a0 + (a1 - a0) * w;
            const double rcs = cfg.scene.line.rcs_db + jitter(rng);
            const double amp = std::pow(10.0, rcs / 20.0);
            out.push_back(makeScattererFromRangeAzimuth(cfg, next_id++, "line", r, az, cfg.scene.ground_z_m,
                                                        amp, uni(rng, 0.0, 2.0 * kPi), rcs));
        }
    }
}

Scatterer makeSinglePointScatterer(const Stage2Config &cfg,
                                   const Stage2RunOptions &opt,
                                   int id)
{
    const int beam_from_az = (cfg.radar.scan_step_deg != 0.0)
        ? static_cast<int>(std::floor(
              (opt.single_scatterer_azimuth_deg - cfg.radar.scan_min_deg) /
              cfg.radar.scan_step_deg + 0.5))
        : 0;
    const int ref_beam_id = std::max(
        0,
        std::min(cfg.radar.beam_count - 1,
                 opt.single_point_beam_id_1based >= 1 ? opt.single_point_beam_id_1based - 1 : beam_from_az));
    const int ref_pulse_idx = std::max(0, cfg.radar.pulse_num / 2);
    const double theta_cmd_deg = cfg.radar.scan_min_deg +
                                 static_cast<double>(ref_beam_id) * cfg.radar.scan_step_deg;
    const double range_sample_float = opt.single_point_expected_bin >= 0
        ? static_cast<double>(cfg.radar.range_crop_start + opt.single_point_expected_bin)
        : (2.0 * opt.single_scatterer_range_m / kC - cfg.radar.sample_delay_sec) * cfg.radar.fs_hz;
    const int range_sample_int = static_cast<int>(std::floor(range_sample_float + 0.5));
    const double range_m = opt.single_point_expected_bin >= 0
        ? 0.5 * kC * (cfg.radar.sample_delay_sec + range_sample_float / cfg.radar.fs_hz)
        : opt.single_scatterer_range_m;

    gmti::target_injection::TargetGlobalConfig global = makeTargetGlobal(cfg);
    const double ref_time_s = pulseTimeSec(cfg.radar, cfg.sim.period_start, ref_beam_id, ref_pulse_idx);
    const Vec3 p_ref = evaluatePlatformState(global, ref_time_s).position;
    const Vec3 v_ref = evaluatePlatformState(global, ref_time_s).velocity;
    const gmti::sim_geometry::LocalPoint local =
        gmti::sim_geometry::makePointFromRangeAzimuth(
            gmti::sim_geometry::LocalPoint(p_ref.x, p_ref.y, p_ref.z),
            gmti::sim_geometry::LocalVelocity(v_ref.x, v_ref.y, v_ref.z),
            range_m,
            theta_cmd_deg,
            p_ref.z,
            cfg.scene.ground_z_m,
            cfg.geometry);

    Scatterer s;
    s.id = id;
    s.type = "single_point";
    s.initial_range_m = range_m;
    s.initial_azimuth_deg = theta_cmd_deg;
    s.amplitude = opt.single_scatterer_amplitude;
    s.phase_rad = 0.0;
    s.rcs_db = 20.0 * std::log10(std::max(opt.single_scatterer_amplitude, 1e-12));
    s.position = Vec3(local.x, local.y, cfg.scene.ground_z_m);
    s.has_ref_geometry = true;
    s.ref_beam_id = ref_beam_id;
    s.ref_pulse_idx = ref_pulse_idx;
    s.ref_time_s = ref_time_s;
    s.ref_platform = p_ref;
    s.range_m_ref = range_m;
    s.range_sample_float_ref = range_sample_float;
    s.range_sample_int_ref = range_sample_int;

    const double sanity = std::fabs(norm(s.position - p_ref) - range_m);
    if (sanity > 1.0e-3) {
        std::cerr << "[stage2][WARN] single_point geometry sanity mismatch: "
                  << sanity << " m\n";
    }
    return s;
}

} // namespace

Scatterer makeScattererFromRangeAzimuth(const Stage2Config &cfg,
                                        int id,
                                        const std::string &type,
                                        double slant_range_m,
                                        double azimuth_deg,
                                        double z_m,
                                        double amplitude,
                                        double phase_rad,
                                        double rcs_db)
{
    Scatterer s;
    s.id = id;
    s.type = type;
    s.initial_range_m = slant_range_m;
    s.initial_azimuth_deg = azimuth_deg;
    s.amplitude = amplitude;
    s.phase_rad = phase_rad;
    s.rcs_db = rcs_db;
    const gmti::sim_geometry::LocalPoint local =
        gmti::sim_geometry::makePointFromRangeAzimuth(
            gmti::sim_geometry::LocalPoint(0.0, 0.0, cfg.platform_height_m),
            gmti::sim_geometry::LocalVelocity(cfg.platform_speed_mps, 0.0, 0.0),
            slant_range_m,
            azimuth_deg,
            cfg.platform_height_m,
            z_m,
            cfg.geometry);
    s.position = Vec3(local.x, local.y, z_m);
    return s;
}

ScattererList generateScene(const Stage2Config &cfg,
                            const Stage2RunOptions &opt,
                            std::mt19937 &rng)
{
    ScattererList out;
    int next_id = 1;
    if (isTargetOnlyMode(opt.scene_mode)) {
        return out;
    }
    if (opt.scene_mode == "point_target_only") {
        if (opt.single_point_beam_id_1based < 0 && opt.single_point_expected_bin < 0) {
            out.push_back(makeScattererFromRangeAzimuth(cfg, next_id++, "single_point",
                                                        opt.single_scatterer_range_m,
                                                        opt.single_scatterer_azimuth_deg,
                                                        cfg.scene.ground_z_m,
                                                        opt.single_scatterer_amplitude,
                                                        0.0,
                                                        20.0 * std::log10(std::max(opt.single_scatterer_amplitude, 1e-12))));
            return out;
        }
        out.push_back(makeSinglePointScatterer(cfg, opt, next_id++));
        return out;
    }
    if (opt.scene_mode == "noise_only") {
        return out;
    }
    if (opt.scene_mode == "area_clutter_only") {
        appendArea(cfg, rng, out, next_id);
        return out;
    }
    if (opt.scene_mode == "clutter_only" || opt.scene_mode == "full") {
        appendArea(cfg, rng, out, next_id);
        appendStrong(cfg, rng, out, next_id);
        appendLines(cfg, rng, out, next_id);
    }
    return out;
}

bool writeSceneTruth(const std::string &truth_dir,
                     const Stage2Config &cfg,
                     const ScattererList &scatterers,
                     std::string &err)
{
    if (!ensureDir(truth_dir)) {
        err = "failed to create truth dir";
        return false;
    }
    std::ofstream scene(joinPath(truth_dir, "scene_truth.csv").c_str());
    std::ofstream area(joinPath(truth_dir, "area_clutter_scatterers.csv").c_str());
    std::ofstream strong(joinPath(truth_dir, "strong_scatterers.csv").c_str());
    std::ofstream line(joinPath(truth_dir, "line_scatterers.csv").c_str());
    if (!scene || !area || !strong || !line) {
        err = "failed to open scene truth csv";
        return false;
    }
    const char *header =
        "scatterer_id,type,x_m,y_m,z_m,amplitude,phase_rad,rcs_db,"
        "initial_range_m,initial_azimuth_deg,range_m,tau_abs_us,tau_rel_us,"
        "range_sample_float,range_sample_int,beam_id,beam_id_0based,beam_id_1based,"
        "theta_cmd_deg,theta_true_deg,target_azimuth_deg,angle_error_deg,visible_by_beam,"
        "ref_beam_id,ref_pulse_idx,ref_time_s,ref_platform_x,ref_platform_y,ref_platform_z,"
        "target_x,target_y,target_z,range_m_ref,range_sample_float_ref,range_sample_int_ref,"
        "ref_platform_e,ref_platform_n,ref_platform_lat,ref_platform_lon,"
        "target_e,target_n,target_lat,target_lon,"
        "echo_delay_sample_center_used,echo_delay_sample_min_used,echo_delay_sample_max_used\n";
    scene << header;
    area << header;
    strong << header;
    line << header;
    for (size_t i = 0; i < scatterers.size(); ++i) {
        const Scatterer &s = scatterers[i];
        const double range_m = s.initial_range_m;
        const double tau_abs_sec = 2.0 * range_m / 299792458.0;
        const double tau_rel_sec = tau_abs_sec - cfg.radar.sample_delay_sec;
        const double range_sample_float = tau_rel_sec * cfg.radar.fs_hz;
        const int range_sample_int = static_cast<int>(std::floor(range_sample_float + 0.5));
        int beam_id_0based = 0;
        if (cfg.radar.scan_step_deg != 0.0) {
            beam_id_0based = static_cast<int>(std::floor(
                (s.initial_azimuth_deg - cfg.radar.scan_min_deg) / cfg.radar.scan_step_deg + 0.5));
        }
        beam_id_0based = std::max(0, std::min(cfg.radar.beam_count - 1, beam_id_0based));
        const int beam_id_1based = beam_id_0based + 1;
        const double theta_cmd_deg = cfg.radar.scan_min_deg +
                                     cfg.radar.scan_step_deg * static_cast<double>(beam_id_0based);
        const double theta_true_deg = theta_cmd_deg;
        double angle_error_deg = s.initial_azimuth_deg - theta_true_deg;
        while (angle_error_deg > 180.0) angle_error_deg -= 360.0;
        while (angle_error_deg < -180.0) angle_error_deg += 360.0;
        const bool visible_by_beam = std::fabs(angle_error_deg) <= cfg.radar.beam_width_deg * 0.5;
        double echo_center = std::numeric_limits<double>::quiet_NaN();
        double echo_min = std::numeric_limits<double>::quiet_NaN();
        double echo_max = std::numeric_limits<double>::quiet_NaN();
        if (s.has_ref_geometry) {
            gmti::target_injection::TargetGlobalConfig global = makeTargetGlobal(cfg);
            const double r_center = norm(s.position - s.ref_platform);
            echo_center = (2.0 * r_center / kC - cfg.radar.sample_delay_sec) * cfg.radar.fs_hz;
            echo_min = std::numeric_limits<double>::infinity();
            echo_max = -std::numeric_limits<double>::infinity();
            for (int p = 0; p < cfg.radar.pulse_num; ++p) {
                const double t = pulseTimeSec(cfg.radar, cfg.sim.period_start, s.ref_beam_id, p);
                const Vec3 pp = evaluatePlatformState(global, t).position;
                const double rr = norm(s.position - pp);
                const double sample = (2.0 * rr / kC - cfg.radar.sample_delay_sec) * cfg.radar.fs_hz;
                echo_min = std::min(echo_min, sample);
                echo_max = std::max(echo_max, sample);
            }
        }
        std::ostringstream row;
        row << std::setprecision(12)
            << s.id << "," << s.type << ","
            << s.position.x << "," << s.position.y << "," << s.position.z << ","
            << s.amplitude << "," << s.phase_rad << "," << s.rcs_db << ","
            << s.initial_range_m << "," << s.initial_azimuth_deg << ","
            << range_m << ","
            << tau_abs_sec * 1.0e6 << ","
            << tau_rel_sec * 1.0e6 << ","
            << range_sample_float << ","
            << range_sample_int << ","
            << beam_id_1based << ","
            << beam_id_0based << ","
            << beam_id_1based << ","
            << theta_cmd_deg << ","
            << theta_true_deg << ","
            << s.initial_azimuth_deg << ","
            << angle_error_deg << ","
            << (visible_by_beam ? 1 : 0) << ",";
        if (s.has_ref_geometry) {
            row << s.ref_beam_id << ","
                << s.ref_pulse_idx << ","
                << s.ref_time_s << ","
                << s.ref_platform.x << ","
                << s.ref_platform.y << ","
                << s.ref_platform.z << ","
                << s.position.x << ","
                << s.position.y << ","
                << s.position.z << ","
                << s.range_m_ref << ","
                << s.range_sample_float_ref << ","
                << s.range_sample_int_ref << ",";
            const gmti::sim_geometry::GeoPoint ref_platform_geo = projectLocalPoint(cfg, s.ref_platform);
            const gmti::sim_geometry::GeoPoint target_geo = projectLocalPoint(cfg, s.position);
            row << ref_platform_geo.e << ","
                << ref_platform_geo.n << ","
                << ref_platform_geo.lat << ","
                << ref_platform_geo.lon << ","
                << target_geo.e << ","
                << target_geo.n << ","
                << target_geo.lat << ","
                << target_geo.lon << ","
                << echo_center << ","
                << echo_min << ","
                << echo_max << "\n";
        } else {
            for (int empty_col = 0; empty_col < 23; ++empty_col) {
                if (empty_col > 0) row << ",";
            }
            row << "\n";
        }
        scene << row.str();
        if (s.type == "area") area << row.str();
        else if (s.type == "strong") strong << row.str();
        else if (s.type == "line") line << row.str();
    }
    return true;
}

} // namespace stage2
} // namespace gmti
