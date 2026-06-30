#include "target_config.h"
#include "truth_writer.h"

#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

using namespace gmti::target_injection;

namespace {

void printUsage()
{
    std::cout
        << "inject_cooperative_target --input-config outputs/stage1/config/temp_config_stage1_newsystem.xml\n"
        << "  --input-data-dir outputs/stage1/data --target-config targets.json\n"
        << "  --output-dir outputs/stage1_target --period-start 0 --period-count 1\n"
        << "  --visibility-mode hard_gate --amplitude-mode snr_db --target-snr-db 30\n";
}

uint64_t fileSize(const std::string &path)
{
    std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
    if (!in) return 0;
    return static_cast<uint64_t>(in.tellg());
}

bool fileExists(const std::string &path)
{
    std::ifstream in(path.c_str(), std::ios::binary);
    return static_cast<bool>(in);
}

void makeOutputDirs(const std::string &root)
{
    ensureDir(root);
    ensureDir(joinPath(root, "config"));
    ensureDir(joinPath(root, "data"));
    ensureDir(joinPath(root, "truth"));
    ensureDir(joinPath(root, "debug"));
    ensureDir(joinPath(root, "logs"));
    ensureDir(joinPath(root, "reports"));
}

std::string chooseInputFile(const RunOptions &run, const RadarConfig &radar)
{
    if (run.background_mode == "zero" || run.input_data_dir == "none") return "";
    const std::string candidate = joinPath(run.input_data_dir, "stage1_background_newprotocol.bin");
    if (fileExists(candidate)) return candidate;
    return radar.input_data_file;
}

void updateStatsFromPacket(const std::vector<uint8_t> &packet, const RadarConfig &radar, InjectionStats &stats)
{
    for (int n = 0; n < radar.pulse_len; ++n) {
        const size_t off = 256U + static_cast<size_t>(n) * 16U;
        const float v[4] = {
            loadF32LE(&packet[off]),
            loadF32LE(&packet[off + 4]),
            loadF32LE(&packet[off + 8]),
            loadF32LE(&packet[off + 12])
        };
        for (int k = 0; k < 4; ++k) {
            if (std::isnan(v[k])) stats.has_nan = true;
            if (std::isinf(v[k])) stats.has_inf = true;
            stats.max_amplitude = std::max(stats.max_amplitude, std::fabs(static_cast<double>(v[k])));
        }
    }
}

} // namespace

int main(int argc, char **argv)
{
    RunOptions run;
    if (!parseCommandLine(argc, argv, run)) {
        printUsage();
        return 0;
    }

    InjectionConfig cfg;
    cfg.run = run;
    std::string err;
    if (!loadRadarConfig(run.input_config, cfg.radar, err)) {
        std::cerr << "[target-inject][ERR] " << err << "\n";
        return 1;
    }
    if (!loadTargetConfig(run.target_config, cfg.radar, cfg.global, cfg.target, err)) {
        std::cerr << "[target-inject][ERR] " << err << "\n";
        return 1;
    }
    applyRunOverrides(run, cfg.global);
    makeOutputDirs(run.output_dir);

    const size_t packet_bytes = 256U + static_cast<size_t>(cfg.radar.pulse_len) * 16U;
    const int prts_per_period = cfg.radar.beam_count * cfg.radar.pulse_num;
    const uint64_t start_prt = static_cast<uint64_t>(run.period_start) * static_cast<uint64_t>(prts_per_period);
    const int beam_end = std::min(cfg.radar.beam_count,
                                  run.beam_start + (run.beam_count > 0 ? run.beam_count : cfg.radar.beam_count));
    const int pulse_end = std::min(cfg.radar.pulse_num,
                                   run.pulse_start + (run.pulse_count > 0 ? run.pulse_count : cfg.radar.pulse_num));
    const int beams_to_write = std::max(0, beam_end - run.beam_start);
    const int pulses_to_write = std::max(0, pulse_end - run.pulse_start);
    const uint64_t total_prts_to_write =
        static_cast<uint64_t>(run.period_count) * static_cast<uint64_t>(beams_to_write) *
        static_cast<uint64_t>(pulses_to_write);
    if (total_prts_to_write == 0U) {
        std::cerr << "[target-inject][ERR] empty generation range\n";
        return 1;
    }

    cfg.radar.input_data_file = chooseInputFile(run, cfg.radar);
    cfg.radar.output_data_file = joinPath(joinPath(run.output_dir, "data"), "stage1_background_with_target.bin");

    std::ifstream in;
    if (run.background_mode != "zero" && run.input_data_dir != "none") {
        const uint64_t input_bytes = fileSize(cfg.radar.input_data_file);
        if (input_bytes == 0 || input_bytes % packet_bytes != 0U) {
            std::cerr << "[target-inject][ERR] invalid input data file: "
                      << cfg.radar.input_data_file << "\n";
            return 1;
        }
        const uint64_t input_prts = input_bytes / packet_bytes;
        const uint64_t max_needed_prt =
            start_prt +
            static_cast<uint64_t>(std::max(0, run.period_count - 1)) * static_cast<uint64_t>(prts_per_period) +
            static_cast<uint64_t>(std::max(0, beam_end - 1)) * static_cast<uint64_t>(cfg.radar.pulse_num) +
            static_cast<uint64_t>(std::max(0, pulse_end - 1)) + 1U;
        if (max_needed_prt > input_prts) {
            std::cerr << "[target-inject][ERR] requested PRT range exceeds input file. input_prts="
                      << input_prts << " need_end=" << max_needed_prt << "\n";
            return 1;
        }
        in.open(cfg.radar.input_data_file.c_str(), std::ios::binary);
    }

    std::ofstream out(cfg.radar.output_data_file.c_str(), std::ios::binary);
    if (!out) {
        std::cerr << "[target-inject][ERR] failed to open output data: "
                  << cfg.radar.output_data_file << "\n";
        return 1;
    }

    TruthWriter truth;
    if (run.write_truth && !truth.open(joinPath(run.output_dir, "truth"), err)) {
        std::cerr << "[target-inject][ERR] " << err << "\n";
        return 1;
    }

    InjectionStats stats;
    std::vector<uint8_t> packet(packet_bytes, 0U);
    uint32_t prt_counter = 0;
    for (int p = 0; p < run.period_count; ++p) {
        const int period_id = run.period_start + p;
        for (int b = run.beam_start; b < beam_end; ++b) {
            for (int m = run.pulse_start; m < pulse_end; ++m) {
                if (run.background_mode == "zero" || run.input_data_dir == "none") {
                    const double t = pulseTimeSec(cfg.radar, period_id, b, m);
                    const double theta = cfg.radar.scan_min_deg + cfg.radar.scan_step_deg * static_cast<double>(b);
                    fillZeroPacketHeader(packet, cfg.radar, cfg.global, prt_counter, t, theta);
                } else {
                    const uint64_t input_prt =
                        static_cast<uint64_t>(period_id) * static_cast<uint64_t>(prts_per_period) +
                        static_cast<uint64_t>(b) * static_cast<uint64_t>(cfg.radar.pulse_num) +
                        static_cast<uint64_t>(m);
                    in.seekg(static_cast<std::streamoff>(input_prt * packet_bytes), std::ios::beg);
                    in.read(reinterpret_cast<char *>(&packet[0]), static_cast<std::streamsize>(packet_bytes));
                    if (in.gcount() != static_cast<std::streamsize>(packet_bytes)) {
                        std::cerr << "[target-inject][ERR] short read at period=" << period_id
                                  << " beam=" << b << " pulse=" << m << "\n";
                        return 1;
                    }
                    ++stats.packets_read;
                }

                PulseTruth pt = injectOnePulse(packet, cfg.radar, cfg.global, cfg.target, period_id, b, m);
                if (pt.injection_enabled) {
                    ++stats.pulses_injected;
                    stats.samples_injected += static_cast<uint64_t>(pt.injected_sample_count);
                }
                if (run.write_truth) truth.writePulse(pt);
                updateStatsFromPacket(packet, cfg.radar, stats);

                out.write(reinterpret_cast<const char *>(&packet[0]), static_cast<std::streamsize>(packet_bytes));
                if (!out) {
                    std::cerr << "[target-inject][ERR] failed to write output packet\n";
                    return 1;
                }
                ++stats.packets_written;
                ++prt_counter;
            }
        }
    }

    if (run.write_truth) {
        truth.writeSummary();
        truth.close();
    }

    if (run.write_config) {
        const std::string out_xml = joinPath(joinPath(run.output_dir, "config"), "temp_config_stage1_target.xml");
        if (!writeOutputConfig(run.input_config, out_xml, cfg.radar.output_data_file, err)) {
            std::cerr << "[target-inject][ERR] " << err << "\n";
            return 1;
        }
    }

    if (!writeTargetInjectionReport(joinPath(joinPath(run.output_dir, "reports"), "target_injection_report.md"),
                                    cfg, stats, "", err)) {
        std::cerr << "[target-inject][ERR] " << err << "\n";
        return 1;
    }

    std::cout << "[target-inject] output=" << cfg.radar.output_data_file
              << " packets=" << stats.packets_written
              << " pulses_injected=" << stats.pulses_injected
              << " samples_injected=" << stats.samples_injected
              << " max_amp=" << stats.max_amplitude << "\n";
    return (stats.has_nan || stats.has_inf) ? 2 : 0;
}
