#include "stage2_config.h"
#include "stage2_lfm_forward.h"
#include "stage2_scene_generator.h"
#include "stage2_validator.h"

#include "../common/SimulationGeometry.h"
#include "../target_injection/lfm_echo_generator.h"
#include "../target_injection/radar_geometry.h"
#include "../target_injection/truth_writer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using namespace gmti::stage2;
using namespace gmti::target_injection;

namespace {

void makeDirs(const std::string &root)
{
    ensureDir(root);
    ensureDir(joinPath(root, "config"));
    ensureDir(joinPath(root, "data"));
    ensureDir(joinPath(root, "truth"));
    ensureDir(joinPath(root, "debug"));
    ensureDir(joinPath(root, "logs"));
    ensureDir(joinPath(root, "reports"));
    ensureDir(joinPath(root, "figures"));
}

void writePlaceholderEval(const std::string &out_dir)
{
    std::ofstream det(joinPath(joinPath(out_dir, "reports"), "detection_eval.csv").c_str());
    det << "matched_count,missed_count,false_alarm_count,detection_rate,mean_range_error_m,max_range_error_m\n";
    det << "0,0,0,0,0,0\n";
    std::ofstream trk(joinPath(joinPath(out_dir, "reports"), "tracking_eval.csv").c_str());
    trk << "confirmed_track_count,track_start_delay_periods,id_switch_count,track_break_count,track_continuity_rate\n";
    trk << "0,0,0,0,0\n";
}

void copyTextFile(const std::string &src, const std::string &dst)
{
    std::ifstream in(src.c_str());
    std::ofstream out(dst.c_str());
    out << in.rdbuf();
}

bool hasMovingTargetOverride(const Stage2RunOptions &opt)
{
    return opt.moving_target_beam_id_1based >= 0 ||
           opt.moving_target_expected_bin >= 0 ||
           opt.moving_target_speed_mps != 0.0 ||
           !opt.moving_target_velocity_mode.empty() ||
           std::isfinite(opt.moving_target_ve_mps) ||
           std::isfinite(opt.moving_target_vn_mps) ||
           std::isfinite(opt.moving_target_radial_speed_mps) ||
           std::isfinite(opt.moving_target_tangential_speed_mps) ||
           opt.moving_target_rcs_db > -900.0;
}

bool isTargetOnlyMode(const std::string &scene_mode)
{
    return scene_mode == "target_only" ||
           scene_mode == "cooperative_target_only" ||
           scene_mode == "empty_target";
}

void applyMovingTargetOverride(const Stage2Config &cfg,
                               const Stage2RunOptions &opt,
                               TargetGlobalConfig &global,
                               TargetConfig &target)
{
    if (!hasMovingTargetOverride(opt)) return;

    const int beam_id_0based = std::max(
        0,
        std::min(cfg.radar.beam_count - 1,
                 opt.moving_target_beam_id_1based >= 1 ? opt.moving_target_beam_id_1based - 1 : 0));
    const int ref_pulse_idx = std::max(0, cfg.radar.pulse_num / 2);
    const double theta_cmd_deg =
        cfg.radar.scan_min_deg + cfg.radar.scan_step_deg * static_cast<double>(beam_id_0based);
    const double range_sample_float = opt.moving_target_expected_bin >= 0
        ? static_cast<double>(cfg.radar.range_crop_start + opt.moving_target_expected_bin)
        : (2.0 * target.slant_range_m / kC - cfg.radar.sample_delay_sec) * cfg.radar.fs_hz;
    const int range_sample_int = static_cast<int>(std::floor(range_sample_float + 0.5));
    const double range_m = 0.5 * kC *
        (cfg.radar.sample_delay_sec + range_sample_float / cfg.radar.fs_hz);

    const double ref_time_s = pulseTimeSec(cfg.radar, cfg.sim.period_start, beam_id_0based, ref_pulse_idx);
    const PlatformState ref_platform = evaluatePlatformState(global, ref_time_s);
    const gmti::sim_geometry::LocalPoint target_ref_local =
        gmti::sim_geometry::makePointFromRangeAzimuth(
            gmti::sim_geometry::LocalPoint(ref_platform.position.x,
                                           ref_platform.position.y,
                                           ref_platform.position.z),
            gmti::sim_geometry::LocalVelocity(ref_platform.velocity.x,
                                              ref_platform.velocity.y,
                                              ref_platform.velocity.z),
            range_m,
            theta_cmd_deg,
            ref_platform.position.z,
            target.height_m,
            cfg.geometry);
    const Vec3 target_ref(target_ref_local.x, target_ref_local.y, target.height_m);
    const gmti::sim_geometry::ENUVelocity ref_vel_en =
        gmti::sim_geometry::localVelocityToEnu(
            gmti::sim_geometry::LocalVelocity(ref_platform.velocity.x,
                                              ref_platform.velocity.y,
                                              ref_platform.velocity.z),
            cfg.geometry);
    const gmti::sim_geometry::LookVectorEN look =
        gmti::sim_geometry::makeAlgorithmLookVectorEN(ref_vel_en.ve,
                                                      ref_vel_en.vn,
                                                      theta_cmd_deg,
                                                      cfg.geometry);
    const gmti::sim_geometry::ENUPoint ref_enu =
        gmti::sim_geometry::localToEnu(
            gmti::sim_geometry::LocalPoint(ref_platform.position.x,
                                           ref_platform.position.y,
                                           ref_platform.position.z),
            cfg.geometry);
    const gmti::sim_geometry::ENUPoint target_enu =
        gmti::sim_geometry::localToEnu(target_ref_local, cfg.geometry);
    const double de = target_enu.e - ref_enu.e;
    const double dn = target_enu.n - ref_enu.n;
    const double ground_range = std::sqrt(de * de + dn * dn);
    const double ur_e = ground_range > 1.0e-9 ? de / ground_range : look.east;
    const double ur_n = ground_range > 1.0e-9 ? dn / ground_range : look.north;
    const double ut_e = -ur_n;
    const double ut_n = ur_e;
    std::string velocity_mode = opt.moving_target_velocity_mode.empty()
        ? "radial"
        : opt.moving_target_velocity_mode;
    if (velocity_mode != "enu" && velocity_mode != "radial" && velocity_mode != "tangential") {
        std::cerr << "[stage2][WARN] unknown moving target velocity mode: "
                  << velocity_mode << ", fallback to radial\n";
        velocity_mode = "radial";
    }
    const double legacy_speed = opt.moving_target_speed_mps != 0.0
        ? opt.moving_target_speed_mps
        : target.radial_velocity_mps;
    double target_ve = 0.0;
    double target_vn = 0.0;
    if (velocity_mode == "enu") {
        target_ve = std::isfinite(opt.moving_target_ve_mps) ? opt.moving_target_ve_mps : 0.0;
        target_vn = std::isfinite(opt.moving_target_vn_mps) ? opt.moving_target_vn_mps : 0.0;
    } else {
        const double radial_speed = std::isfinite(opt.moving_target_radial_speed_mps)
            ? opt.moving_target_radial_speed_mps
            : (velocity_mode == "tangential" ? 0.0 : legacy_speed);
        const double tangential_speed = std::isfinite(opt.moving_target_tangential_speed_mps)
            ? opt.moving_target_tangential_speed_mps
            : (velocity_mode == "tangential" ? legacy_speed : 0.0);
        target_ve = radial_speed * ur_e + tangential_speed * ut_e;
        target_vn = radial_speed * ur_n + tangential_speed * ut_n;
    }
    const double target_vr_self = target_ve * ur_e + target_vn * ur_n;
    const double target_vt_self = target_ve * ut_e + target_vn * ut_n;
    const double lambda = kC / cfg.radar.fc_hz;
    const double af_motion_truth = lambda > 0.0 ? 2.0 * target_vr_self / lambda : 0.0;
    const gmti::sim_geometry::LocalVelocity target_vel_local =
        gmti::sim_geometry::enuVelocityToLocal(
            gmti::sim_geometry::ENUVelocity(target_ve, target_vn, 0.0),
            cfg.geometry);
    target.v = Vec3(target_vel_local.vx, target_vel_local.vy, target_vel_local.vz);
    target.p0 = target_ref - target.v * ref_time_s;
    target.slant_range_m = range_m;
    target.azimuth_deg = theta_cmd_deg;
    target.init_mode = "project_local";

    target.has_ref_geometry = true;
    target.ref_beam_id = beam_id_0based;
    target.ref_pulse_idx = ref_pulse_idx;
    target.ref_time_s = ref_time_s;
    target.ref_platform = ref_platform.position;
    target.ref_target = target_ref;
    target.ref_range_m = range_m;
    target.ref_range_sample_float = range_sample_float;
    target.ref_range_sample_int = range_sample_int;
    target.echo_delay_sample_center_used =
        (2.0 * norm(target_ref - ref_platform.position) / kC - cfg.radar.sample_delay_sec) *
        cfg.radar.fs_hz;
    target.override_speed_mps = std::sqrt(target_ve * target_ve + target_vn * target_vn);
    target.velocity_mode = velocity_mode;
    target.target_ve_mps = target_ve;
    target.target_vn_mps = target_vn;
    target.target_vr_self_mps = target_vr_self;
    target.target_vt_self_mps = target_vt_self;
    target.af_motion_truth_hz = af_motion_truth;

    if (opt.moving_target_rcs_db > -900.0) {
        global.amplitude_mode = "direct_amplitude";
        global.direct_amplitude = std::pow(10.0, opt.moving_target_rcs_db / 20.0);
        target.override_rcs_db = opt.moving_target_rcs_db;
    }
}

} // namespace

int main(int argc, char **argv)
{
    Stage2RunOptions opt;
    parseStage2CommandLine(argc, argv, opt);

    Stage2Config cfg;
    std::string err;
    if (!loadStage2Config(opt.stage2_config, cfg, err)) {
        std::cerr << "[stage2][ERR] " << err << "\n";
        return 1;
    }
    if (opt.period_start >= 0) cfg.sim.period_start = opt.period_start;
    if (opt.period_count > 0) cfg.sim.period_count = opt.period_count;
    if (!opt.target_config.empty()) cfg.target.target_config_path = opt.target_config;

    makeDirs(opt.output_dir);
    copyTextFile(opt.stage2_config, joinPath(joinPath(opt.output_dir, "config"), "stage2_config.json"));
    std::mt19937 rng(cfg.sim.random_seed);
    ScattererList scatterers = generateScene(cfg, opt, rng);
    if (!writeSceneTruth(joinPath(opt.output_dir, "truth"), cfg, scatterers, err)) {
        std::cerr << "[stage2][ERR] " << err << "\n";
        return 1;
    }

    TargetGlobalConfig global = makeTargetGlobal(cfg);
    TargetConfig target;
    bool target_loaded = false;
    if (opt.target_enabled && cfg.target.enabled) {
        if (!loadTargetConfig(cfg.target.target_config_path, cfg.radar, global, target, err)) {
            std::cerr << "[stage2][ERR] " << err << "\n";
            return 1;
        }
        global.amplitude_mode = cfg.target.amplitude_mode;
        global.target_snr_db = cfg.target.target_snr_db;
        applyMovingTargetOverride(cfg, opt, global, target);
        target_loaded = true;
    }

    const size_t packet_bytes = 256U + static_cast<size_t>(cfg.radar.pulse_len) * 16U;
    std::vector<uint8_t> packet(packet_bytes, 0U);
    const std::string data_file = joinPath(joinPath(opt.output_dir, "data"), "stage2_statistical_newprotocol.bin");
    std::ofstream out(data_file.c_str(), std::ios::binary);
    if (!out) {
        std::cerr << "[stage2][ERR] failed to open output data\n";
        return 1;
    }
    TruthWriter target_truth;
    if (target_loaded && !target_truth.open(joinPath(opt.output_dir, "truth"), err)) {
        std::cerr << "[stage2][ERR] " << err << "\n";
        return 1;
    }

    Stage2Stats stats;
    const auto t0 = std::chrono::steady_clock::now();
    uint32_t prt_counter = 0;
    for (int pp = 0; pp < cfg.sim.period_count; ++pp) {
        const int period_id = cfg.sim.period_start + pp;
        for (int b = 0; b < cfg.radar.beam_count; ++b) {
            for (int m = 0; m < cfg.radar.pulse_num; ++m) {
                const double t = pulseTimeSec(cfg.radar, period_id, b, m);
                const double theta = cfg.radar.scan_min_deg + cfg.radar.scan_step_deg * static_cast<double>(b);
                fillZeroPacketHeader(packet, cfg.radar, global, prt_counter, t, theta);
                if (opt.scene_mode != "noise_only") {
                    addScatterersToPacket(packet, cfg.radar, global, scatterers, period_id, b, m,
                                          global.beam_gain_threshold, stats);
                }
                if (cfg.scene.noise.enabled &&
                    opt.scene_mode != "point_target_only" &&
                    !isTargetOnlyMode(opt.scene_mode)) {
                    addThermalNoise(packet, cfg.radar, cfg.scene.noise.noise_power, rng, stats);
                }
                if (target_loaded) {
                    PulseTruth pt = injectOnePulse(packet, cfg.radar, global, target, period_id, b, m);
                    if (pt.injection_enabled) {
                        ++stats.target_pulses_injected;
                        stats.target_samples_injected += static_cast<uint64_t>(pt.injected_sample_count);
                    }
                    target_truth.writePulse(pt);
                }
                scanPacketStats(packet, cfg.radar, stats);
                out.write(reinterpret_cast<const char *>(&packet[0]), static_cast<std::streamsize>(packet.size()));
                if (!out) {
                    std::cerr << "[stage2][ERR] write failed\n";
                    return 1;
                }
                ++stats.packets_written;
                ++prt_counter;
            }
        }
    }
    if (target_loaded) {
        target_truth.writeSummary();
        target_truth.close();
        std::rename(joinPath(joinPath(opt.output_dir, "truth"), "truth_pulse.csv").c_str(),
                    joinPath(joinPath(opt.output_dir, "truth"), "target_truth_pulse.csv").c_str());
        std::rename(joinPath(joinPath(opt.output_dir, "truth"), "truth_beam_summary.csv").c_str(),
                    joinPath(joinPath(opt.output_dir, "truth"), "target_truth_beam_summary.csv").c_str());
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(t1 - t0).count();
    writeStage2OutputConfig(cfg, joinPath(joinPath(opt.output_dir, "config"), "temp_config_stage2_newsystem.xml"), data_file, err);
    writePlaceholderEval(opt.output_dir);
    {
        std::ofstream log(joinPath(joinPath(opt.output_dir, "logs"), "simulate_stage2.log").c_str());
        log << "scene_mode=" << opt.scene_mode << "\n";
        log << "period_start=" << cfg.sim.period_start << "\n";
        log << "period_count=" << cfg.sim.period_count << "\n";
        log << "scatterers=" << scatterers.size() << "\n";
        log << "packets_written=" << stats.packets_written << "\n";
    }
    if (!writeStage2Report(joinPath(joinPath(opt.output_dir, "reports"), "stage2_simulation_report.md"),
                           cfg, opt, scatterers, stats, elapsed, err)) {
        std::cerr << "[stage2][ERR] " << err << "\n";
        return 1;
    }

    std::cout << "[stage2] output=" << data_file
              << " packets=" << stats.packets_written
              << " scatterers=" << scatterers.size()
              << " scatterer_echoes=" << stats.scatterer_echoes
              << " target_pulses=" << stats.target_pulses_injected
              << " elapsed_ms=" << elapsed * 1000.0 << "\n";
    return (stats.has_nan || stats.has_inf) ? 2 : 0;
}
