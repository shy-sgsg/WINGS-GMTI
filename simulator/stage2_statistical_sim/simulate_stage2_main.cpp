#include "stage2_config.h"
#include "stage2_lfm_forward.h"
#include "stage2_scene_generator.h"
#include "stage2_validator.h"

#include "../target_injection/lfm_echo_generator.h"
#include "../target_injection/truth_writer.h"

#include <chrono>
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
    if (!writeSceneTruth(joinPath(opt.output_dir, "truth"), scatterers, err)) {
        std::cerr << "[stage2][ERR] " << err << "\n";
        return 1;
    }

    TargetGlobalConfig global = makeTargetGlobal(cfg);
    TargetConfig target;
    bool target_loaded = false;
    if (opt.target_enabled && cfg.target.enabled && opt.scene_mode != "point_target_only" && opt.scene_mode != "noise_only") {
        if (!loadTargetConfig(cfg.target.target_config_path, cfg.radar, global, target, err)) {
            std::cerr << "[stage2][ERR] " << err << "\n";
            return 1;
        }
        global.amplitude_mode = cfg.target.amplitude_mode;
        global.target_snr_db = cfg.target.target_snr_db;
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
                if (cfg.scene.noise.enabled && opt.scene_mode != "point_target_only") {
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
