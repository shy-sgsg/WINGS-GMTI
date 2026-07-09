#include "stage2_config.h"
#include "stage2_lfm_forward.h"
#include "stage2_scene_generator.h"
#include "stage2_validator.h"

#include "../common/SimulationGeometry.h"
#include "../target_injection/lfm_echo_generator.h"
#include "../target_injection/radar_geometry.h"
#include "../target_injection/truth_writer.h"
#include "dbs/NewProtocolLayout.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
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

bool isTargetOnlyMode(const std::string &scene_mode)
{
    return scene_mode == "target_only" ||
           scene_mode == "cooperative_target_only" ||
           scene_mode == "empty_target";
}

bool isContinuousAreaModel(const std::string &model)
{
    return model == "continuous_texture" ||
           model == "continuous_surface" ||
           model == "continuous_grid" ||
           model == "grid_texture";
}

void applyClutterOverrides(Stage2Config &cfg, const Stage2RunOptions &opt)
{
    if (std::isfinite(opt.clutter_amplitude_scale)) {
        cfg.scene.clutter_amplitude_scale = opt.clutter_amplitude_scale;
    }
    if (opt.area_clutter_scatterer_count >= 0) {
        cfg.scene.area.scatterer_count = opt.area_clutter_scatterer_count;
    }
    if (std::isfinite(opt.area_clutter_mean_power)) {
        cfg.scene.area.mean_power = opt.area_clutter_mean_power;
    }
    if (std::isfinite(opt.area_clutter_texture_sigma)) {
        cfg.scene.area.texture_sigma = opt.area_clutter_texture_sigma;
    }
    if (opt.strong_scatterer_count >= 0) {
        cfg.scene.strong.count = opt.strong_scatterer_count;
    }
    if (std::isfinite(opt.strong_rcs_db_min)) {
        cfg.scene.strong.rcs_db_min = opt.strong_rcs_db_min;
    }
    if (std::isfinite(opt.strong_rcs_db_max)) {
        cfg.scene.strong.rcs_db_max = opt.strong_rcs_db_max;
    }
    if (opt.line_scatterer_count >= 0) {
        cfg.scene.line.line_count = opt.line_scatterer_count;
    }
    if (opt.line_points_per_line >= 0) {
        cfg.scene.line.points_per_line = opt.line_points_per_line;
    }
    if (std::isfinite(opt.line_rcs_db)) {
        cfg.scene.line.rcs_db = opt.line_rcs_db;
    }
    if (std::isfinite(opt.noise_power)) {
        cfg.scene.noise.noise_power = opt.noise_power;
    }
}

void writeScenarioResolved(const Stage2RunConfig &run)
{
    std::ofstream out(joinPath(run.output_dir, "scenario_resolved.json").c_str());
    out << std::setprecision(12);
    out << "{\n";
    out << "  \"case_id\": \"" << run.case_id << "\",\n";
    out << "  \"output_dir\": \"" << run.output_dir << "\",\n";
    out << "  \"scene_mode\": \"" << run.scene_mode << "\",\n";
    out << "  \"beam_index_base\": " << run.beam_index_base << ",\n";
    out << "  \"scene\": {\n";
    out << "    \"mode\": \"" << run.scene_mode << "\",\n";
    out << "    \"range_min_m\": " << run.cfg.scene.range_min_m << ",\n";
    out << "    \"range_max_m\": " << run.cfg.scene.range_max_m << ",\n";
    out << "    \"azimuth_min_deg\": " << run.cfg.scene.azimuth_min_deg << ",\n";
    out << "    \"azimuth_max_deg\": " << run.cfg.scene.azimuth_max_deg << ",\n";
    out << "    \"ground_z_m\": " << run.cfg.scene.ground_z_m << ",\n";
    out << "    \"clutter_amplitude_scale\": " << run.cfg.scene.clutter_amplitude_scale << ",\n";
    out << "    \"single_point\": {\n";
    out << "      \"range_m\": " << run.legacy_scene_options.single_scatterer_range_m << ",\n";
    out << "      \"beam_id\": " << run.legacy_scene_options.single_point_beam_id_1based << ",\n";
    out << "      \"expected_bin\": " << run.legacy_scene_options.single_point_expected_bin << ",\n";
    out << "      \"azimuth_deg\": " << run.legacy_scene_options.single_scatterer_azimuth_deg << ",\n";
    out << "      \"amplitude\": " << run.legacy_scene_options.single_scatterer_amplitude << "\n";
    out << "    },\n";
    out << "    \"area_clutter\": {\n";
    out << "      \"enabled\": " << (run.cfg.scene.area.enabled ? "true" : "false") << ",\n";
    out << "      \"model\": \"" << run.cfg.scene.area.model << "\",\n";
    out << "      \"scatterer_count\": " << run.cfg.scene.area.scatterer_count << ",\n";
    out << "      \"mean_power\": " << run.cfg.scene.area.mean_power << ",\n";
    out << "      \"texture_sigma\": " << run.cfg.scene.area.texture_sigma << ",\n";
    out << "      \"spatial_cell_m\": " << run.cfg.scene.area.spatial_cell_m << "\n";
    out << "    },\n";
    out << "    \"strong_scatterers\": {\n";
    out << "      \"enabled\": " << (run.cfg.scene.strong.enabled ? "true" : "false") << ",\n";
    out << "      \"count\": " << run.cfg.scene.strong.count << ",\n";
    out << "      \"rcs_db_min\": " << run.cfg.scene.strong.rcs_db_min << ",\n";
    out << "      \"rcs_db_max\": " << run.cfg.scene.strong.rcs_db_max << "\n";
    out << "    },\n";
    out << "    \"line_scatterers\": {\n";
    out << "      \"enabled\": " << (run.cfg.scene.line.enabled ? "true" : "false") << ",\n";
    out << "      \"line_count\": " << run.cfg.scene.line.line_count << ",\n";
    out << "      \"points_per_line\": " << run.cfg.scene.line.points_per_line << ",\n";
    out << "      \"rcs_db\": " << run.cfg.scene.line.rcs_db << "\n";
    out << "    },\n";
    out << "    \"thermal_noise\": {\n";
    out << "      \"enabled\": " << (run.cfg.scene.noise.enabled ? "true" : "false") << ",\n";
    out << "      \"noise_power\": " << run.cfg.scene.noise.noise_power << "\n";
    out << "    }\n";
    out << "  },\n";
    out << "  \"scan\": {\n";
    out << "    \"scan_min_deg\": " << run.cfg.radar.scan_min_deg << ",\n";
    out << "    \"scan_step_deg\": " << run.cfg.radar.scan_step_deg << ",\n";
    out << "    \"beam_count\": " << run.cfg.radar.beam_count << "\n";
    out << "  },\n";
    out << "  \"targets\": [\n";
    for (size_t i = 0; i < run.targets.size(); ++i) {
        const Stage2RunTarget &t = run.targets[i];
        out << "    {\n";
        out << "      \"target_id\": \"" << t.target_id << "\",\n";
        out << "      \"enabled\": " << (t.enabled ? "true" : "false") << ",\n";
        out << "      \"init_type\": \"" << t.init_type << "\",\n";
        out << "      \"beam_id\": " << t.beam_id << ",\n";
        out << "      \"expected_bin\": " << t.expected_bin << ",\n";
        out << "      \"theta_cmd_deg\": " << t.theta_cmd_deg << ",\n";
        out << "      \"azimuth_deg\": " << t.azimuth_deg << ",\n";
        out << "      \"azimuth_offset_deg\": " << t.azimuth_offset_deg << ",\n";
        out << "      \"motion_type\": \"" << t.motion_type << "\",\n";
        out << "      \"ve_mps\": " << t.ve_mps << ",\n";
        out << "      \"vn_mps\": " << t.vn_mps << ",\n";
        out << "      \"amplitude_type\": \"" << t.amplitude_type << "\",\n";
        out << "      \"snr_db\": " << t.snr_db << ",\n";
        out << "      \"visibility_type\": \"" << t.visibility_type << "\"\n";
        out << "    }" << (i + 1 < run.targets.size() ? "," : "") << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

int generateStage2Data(const Stage2RunConfig &run)
{
    Stage2Config cfg = run.cfg;
    Stage2RunOptions scene_opt = run.use_legacy_scene_options
        ? run.legacy_scene_options
        : Stage2RunOptions();
    scene_opt.output_dir = run.output_dir;
    scene_opt.scene_mode = run.scene_mode;
    scene_opt.target_enabled = !run.targets.empty();

    std::string err;
    makeDirs(run.output_dir);
    writeScenarioResolved(run);

    std::mt19937 rng(cfg.sim.random_seed);
    ScattererList scatterers = generateScene(cfg, scene_opt, rng);
    if (!writeSceneTruth(joinPath(run.output_dir, "truth"), cfg, scatterers, err)) {
        std::cerr << "[stage2][ERR] " << err << "\n";
        return 1;
    }

    TargetGlobalConfig global = run.global;
    global.amplitude_mode = cfg.target.amplitude_mode;
    global.target_snr_db = cfg.target.target_snr_db;

    std::vector<TargetConfig> targets;
    for (size_t i = 0; i < run.targets.size(); ++i) {
        if (run.targets[i].enabled) {
            targets.push_back(run.targets[i].target);
        }
    }

    const size_t packet_bytes = gmti::new_protocol::packetBytes(
        static_cast<size_t>(cfg.radar.pulse_len),
        static_cast<size_t>(std::max(2, cfg.radar.new_protocol_channel_count)),
        cfg.radar.iq_data_type.empty() ? "float32" : cfg.radar.iq_data_type);
    std::vector<uint8_t> packet(packet_bytes, 0U);
    const std::string data_file =
        joinPath(joinPath(run.output_dir, "data"), "stage2_statistical_newprotocol.bin");
    std::ofstream out(data_file.c_str(), std::ios::binary);
    if (!out) {
        std::cerr << "[stage2][ERR] failed to open output data\n";
        return 1;
    }

    TruthWriter target_truth;
    const bool write_target_truth = run.truth_output && !targets.empty();
    if (write_target_truth) {
        target_truth.setCaseId(run.case_id);
        if (!target_truth.open(joinPath(run.output_dir, "truth"), err)) {
            std::cerr << "[stage2][ERR] " << err << "\n";
            return 1;
        }
    }

    Stage2Stats stats;
    const auto t0 = std::chrono::steady_clock::now();
    uint32_t prt_counter = 0;
    for (int pp = 0; pp < cfg.sim.period_count; ++pp) {
        const int period_id = cfg.sim.period_start + pp;
        for (int b = 0; b < cfg.radar.beam_count; ++b) {
            for (int m = 0; m < cfg.radar.pulse_num; ++m) {
                const double t = pulseTimeSec(cfg.radar, period_id, b, m);
                const double theta =
                    cfg.radar.scan_min_deg + cfg.radar.scan_step_deg * static_cast<double>(b);
                fillZeroPacketHeader(packet, cfg.radar, global, prt_counter, t, theta);
                if (run.scene_mode != "noise_only") {
                    if (cfg.scene.area.enabled && isContinuousAreaModel(cfg.scene.area.model)) {
                        addContinuousAreaClutter(packet, cfg.radar, global, cfg.scene,
                                                 cfg.sim.random_seed, period_id, b, m, stats);
                    } else {
                        addScatterersToPacket(packet, cfg.radar, global, scatterers, period_id, b, m,
                                              global.beam_gain_threshold, stats);
                    }
                }
                if (cfg.scene.noise.enabled &&
                    run.scene_mode != "point_target_only" &&
                    !isTargetOnlyMode(run.scene_mode)) {
                    addThermalNoise(packet, cfg.radar, cfg.scene.noise.noise_power, rng, stats);
                }
                for (size_t ti = 0; ti < targets.size(); ++ti) {
                    PulseTruth pt =
                        injectOnePulse(packet, cfg.radar, global, targets[ti], period_id, b, m);
                    if (pt.injection_enabled) {
                        ++stats.target_pulses_injected;
                        stats.target_samples_injected +=
                            static_cast<uint64_t>(pt.injected_sample_count);
                    }
                    if (write_target_truth) {
                        target_truth.writePulse(pt);
                    }
                }
                scanPacketStats(packet, cfg.radar, stats);
                out.write(reinterpret_cast<const char *>(&packet[0]),
                          static_cast<std::streamsize>(packet.size()));
                if (!out) {
                    std::cerr << "[stage2][ERR] write failed\n";
                    return 1;
                }
                ++stats.packets_written;
                ++prt_counter;
            }
        }
    }
    if (write_target_truth) {
        target_truth.writeSummary();
        target_truth.close();
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(t1 - t0).count();
    if (!writeStage2OutputConfig(cfg,
                                 joinPath(joinPath(run.output_dir, "config"),
                                          "temp_config_stage2_newsystem.xml"),
                                 data_file,
                                 err)) {
        std::cerr << "[stage2][ERR] " << err << "\n";
        return 1;
    }
    writePlaceholderEval(run.output_dir);
    {
        std::ofstream log(joinPath(joinPath(run.output_dir, "logs"), "simulate_stage2.log").c_str());
        log << "case_id=" << run.case_id << "\n";
        log << "scene_mode=" << run.scene_mode << "\n";
        log << "period_start=" << cfg.sim.period_start << "\n";
        log << "period_count=" << cfg.sim.period_count << "\n";
        log << "scatterers=" << scatterers.size() << "\n";
        log << "targets=" << targets.size() << "\n";
        log << "clutter_amplitude_scale=" << cfg.scene.clutter_amplitude_scale << "\n";
        log << "area_clutter_model=" << cfg.scene.area.model << "\n";
        log << "area_clutter_scatterer_count=" << cfg.scene.area.scatterer_count << "\n";
        log << "continuous_area_packets=" << stats.continuous_area_packets << "\n";
        log << "continuous_area_samples=" << stats.continuous_area_samples << "\n";
        log << "area_clutter_mean_power=" << cfg.scene.area.mean_power << "\n";
        log << "area_clutter_texture_sigma=" << cfg.scene.area.texture_sigma << "\n";
        log << "strong_scatterer_count=" << cfg.scene.strong.count << "\n";
        log << "strong_rcs_db_min=" << cfg.scene.strong.rcs_db_min << "\n";
        log << "strong_rcs_db_max=" << cfg.scene.strong.rcs_db_max << "\n";
        log << "line_scatterer_count=" << cfg.scene.line.line_count << "\n";
        log << "line_points_per_line=" << cfg.scene.line.points_per_line << "\n";
        log << "line_rcs_db=" << cfg.scene.line.rcs_db << "\n";
        log << "noise_power=" << cfg.scene.noise.noise_power << "\n";
        log << "packets_written=" << stats.packets_written << "\n";
    }
    if (!writeStage2Report(joinPath(joinPath(run.output_dir, "reports"),
                                    "stage2_simulation_report.md"),
                           cfg, scene_opt, scatterers, stats, elapsed, err)) {
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

} // namespace

int main(int argc, char **argv)
{
    Stage2RunOptions opt;
    if (!parseStage2CommandLine(argc, argv, opt)) {
        std::cerr << "[stage2][ERR] " << opt.parse_error << "\n";
        return 1;
    }
    std::string err;
    Stage2RunConfig run;
    if (!opt.run_config.empty()) {
        if (!loadStage2RunConfig(opt.run_config, run, err)) {
            std::cerr << "[stage2][ERR] " << err << "\n";
            return 1;
        }
    } else {
        if (!makeLegacyStage2RunConfig(opt, run, err)) {
            std::cerr << "[stage2][ERR] " << err << "\n";
            return 1;
        }
        applyClutterOverrides(run.cfg, opt);
        run.global = makeTargetGlobal(run.cfg);
    }
    if (opt.validate && !validateStage2RunConfig(run, err) && !run.targets.empty()) {
        std::cerr << "[stage2][ERR] " << err << "\n";
        return 1;
    }
    return generateStage2Data(run);
}
