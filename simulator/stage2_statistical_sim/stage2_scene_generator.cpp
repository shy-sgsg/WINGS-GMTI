#include "stage2_scene_generator.h"

#include "../target_injection/target_common.h"

#include <cmath>
#include <fstream>
#include <iomanip>

namespace gmti {
namespace stage2 {

using gmti::target_injection::Vec3;
using gmti::target_injection::deg2rad;
using gmti::target_injection::ensureDir;
using gmti::target_injection::joinPath;
using gmti::target_injection::kPi;

namespace {

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
        out.push_back(makeScattererFromRangeAzimuth(next_id++, "area", r, az, cfg.scene.ground_z_m,
                                                    cfg.platform_height_m, amp, uni(rng, 0.0, 2.0 * kPi), 20.0 * std::log10(std::max(amp, 1e-12))));
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
        out.push_back(makeScattererFromRangeAzimuth(next_id++, "strong", r, az, cfg.scene.ground_z_m,
                                                    cfg.platform_height_m, amp, uni(rng, 0.0, 2.0 * kPi), rcs));
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
            out.push_back(makeScattererFromRangeAzimuth(next_id++, "line", r, az, cfg.scene.ground_z_m,
                                                        cfg.platform_height_m, amp, uni(rng, 0.0, 2.0 * kPi), rcs));
        }
    }
}

} // namespace

Scatterer makeScattererFromRangeAzimuth(int id,
                                        const std::string &type,
                                        double slant_range_m,
                                        double azimuth_deg,
                                        double z_m,
                                        double platform_height_m,
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
    const double dz = z_m - platform_height_m;
    const double horiz = std::sqrt(std::max(0.0, slant_range_m * slant_range_m - dz * dz));
    const double az = deg2rad(azimuth_deg);
    s.position = Vec3(horiz * std::cos(az), horiz * std::sin(az), z_m);
    return s;
}

ScattererList generateScene(const Stage2Config &cfg,
                            const Stage2RunOptions &opt,
                            std::mt19937 &rng)
{
    ScattererList out;
    int next_id = 1;
    if (opt.scene_mode == "point_target_only") {
        out.push_back(makeScattererFromRangeAzimuth(next_id++, "single_point",
                                                    opt.single_scatterer_range_m,
                                                    opt.single_scatterer_azimuth_deg,
                                                    cfg.scene.ground_z_m,
                                                    cfg.platform_height_m,
                                                    opt.single_scatterer_amplitude,
                                                    0.0,
                                                    20.0 * std::log10(std::max(opt.single_scatterer_amplitude, 1e-12))));
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
    const char *header = "scatterer_id,type,x_m,y_m,z_m,amplitude,phase_rad,rcs_db,initial_range_m,initial_azimuth_deg\n";
    scene << header;
    area << header;
    strong << header;
    line << header;
    for (size_t i = 0; i < scatterers.size(); ++i) {
        const Scatterer &s = scatterers[i];
        std::ostringstream row;
        row << std::setprecision(12)
            << s.id << "," << s.type << ","
            << s.position.x << "," << s.position.y << "," << s.position.z << ","
            << s.amplitude << "," << s.phase_rad << "," << s.rcs_db << ","
            << s.initial_range_m << "," << s.initial_azimuth_deg << "\n";
        scene << row.str();
        if (s.type == "area") area << row.str();
        else if (s.type == "strong") strong << row.str();
        else if (s.type == "line") line << row.str();
    }
    return true;
}

} // namespace stage2
} // namespace gmti
