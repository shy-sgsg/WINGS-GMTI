#include "stage1_validator.h"
#include "dbs/NewProtocolLayout.hpp"

#include <fstream>
#include <iomanip>
#include <algorithm>

namespace gmti {
namespace sim_stage1 {

bool writeParameterAudit(const Stage1OldSystemConfig &old_cfg,
                         const Stage1NewSystemConfig &new_cfg,
                         const FileAudit &audit,
                         const std::string &output_dir,
                         std::string &err)
{
    const std::string report_dir = pathJoin(output_dir, "reports");
    std::ofstream js(pathJoin(report_dir, "parameter_audit.json").c_str());
    std::ofstream md(pathJoin(report_dir, "parameter_audit.md").c_str());
    if (!js || !md) {
        err = "failed to open parameter audit outputs";
        return false;
    }
    js << std::fixed << std::setprecision(6);
    js << "{\n"
       << "  \"old\": {\n"
       << "    \"data_ch1\": \"" << old_cfg.data_ch1 << "\",\n"
       << "    \"data_ch2\": \"" << old_cfg.data_ch2 << "\",\n"
       << "    \"pulse_len\": " << old_cfg.pulse_len << ",\n"
       << "    \"rg_len\": " << old_cfg.rg_len << ",\n"
       << "    \"pulse_num\": " << old_cfg.pulse_num << ",\n"
       << "    \"read_pulse_num\": " << old_cfg.read_pulse_num << ",\n"
       << "    \"fc_ghz\": " << old_cfg.fc_ghz << ",\n"
       << "    \"br_mhz\": " << old_cfg.br_mhz << ",\n"
       << "    \"fs_mhz\": " << old_cfg.fs_mhz << ",\n"
       << "    \"tr_us\": " << old_cfg.tr_us << ",\n"
       << "    \"prf_hz\": " << old_cfg.prf_hz << ",\n"
       << "    \"scan_min_deg\": " << old_cfg.scan_min_deg << ",\n"
       << "    \"scan_max_deg\": " << old_cfg.scan_max_deg << ",\n"
       << "    \"az_count\": " << old_cfg.az_count << ",\n"
       << "    \"d_chan\": " << old_cfg.d_chan << "\n"
       << "  },\n"
       << "  \"new\": {\n"
       << "    \"scan_min_deg\": " << new_cfg.scan_min_deg << ",\n"
       << "    \"scan_max_deg\": " << new_cfg.scan_max_deg << ",\n"
       << "    \"scan_step_deg\": " << new_cfg.scan_step_deg << ",\n"
       << "    \"beam_count\": " << new_cfg.beam_count << ",\n"
       << "    \"pulse_num_new\": " << new_cfg.pulse_num_new << ",\n"
       << "    \"ddc_len_new\": " << new_cfg.ddc_len_new << ",\n"
       << "    \"fft_len_new\": " << new_cfg.fft_len_new << ",\n"
       << "    \"pc_crop_start\": " << new_cfg.pc_crop_start << ",\n"
       << "    \"pc_crop_len\": " << new_cfg.pc_crop_len << ",\n"
       << "    \"fc_new_ghz\": " << new_cfg.fc_new_ghz << ",\n"
       << "    \"br_new_mhz\": " << new_cfg.br_new_mhz << ",\n"
       << "    \"fs_new_mhz\": " << new_cfg.fs_new_mhz << ",\n"
       << "    \"tr_new_us\": " << new_cfg.tr_new_us << ",\n"
       << "    \"prf_new_hz\": " << new_cfg.prf_new_hz << ",\n"
       << "    \"sample_delay_us\": " << new_cfg.sample_delay_us << ",\n"
       << "    \"sample_window_us\": " << new_cfg.sample_window_us << ",\n"
       << "    \"period_time_sec\": " << new_cfg.periodTimeSec() << ",\n"
       << "    \"realtime_80_sec\": " << new_cfg.realtime80Sec() << ",\n"
       << "    \"strict_gate_sec\": 4.8\n"
       << "  },\n"
       << "  \"file_audit\": {\n"
       << "    \"complex_format\": \"" << audit.complex_format << "\",\n"
       << "    \"ch1_bytes\": " << audit.ch1_bytes << ",\n"
       << "    \"ch2_bytes\": " << audit.ch2_bytes << ",\n"
       << "    \"bytes_per_old_period_per_channel\": " << audit.bytes_per_old_period_per_channel << ",\n"
       << "    \"period_count\": " << audit.period_count << ",\n"
       << "    \"ok\": " << boolText(audit.ok) << ",\n"
       << "    \"message\": \"" << audit.message << "\"\n"
       << "  }\n"
       << "}\n";

    md << "# Stage 1 parameter audit\n\n";
    md << "## Old system\n\n";
    md << "- data_ch1: " << old_cfg.data_ch1 << "\n";
    md << "- data_ch2: " << old_cfg.data_ch2 << "\n";
    md << "- format: " << audit.complex_format << "\n";
    md << "- beams per period: " << old_cfg.az_count << "\n";
    md << "- pulses per beam: " << old_cfg.pulse_num << "\n";
    md << "- complex samples per pulse: " << old_cfg.pulse_len << "\n";
    md << "- PRF/fs/fc/Tr: " << old_cfg.prf_hz << " Hz / " << old_cfg.fs_mhz
       << " MHz / " << old_cfg.fc_ghz << " GHz / " << old_cfg.tr_us << " us\n\n";
    md << "## New system\n\n";
    md << "- scan: " << new_cfg.scan_min_deg << " to " << new_cfg.scan_max_deg
       << " deg, step " << new_cfg.scan_step_deg << ", beam_count " << new_cfg.beam_count << "\n";
    md << "- 61 x 130 / 1300 = " << new_cfg.periodTimeSec() << " s\n";
    md << "- 80 percent realtime gate: " << new_cfg.realtime80Sec() << " s\n";
    md << "- strict engineering gate: 4.8 s\n";
    md << "- sample_window_us x fs_new = " << new_cfg.sample_window_us << " x "
       << new_cfg.fs_new_mhz << " = " << (new_cfg.sample_window_us * new_cfg.fs_new_mhz) << " points\n\n";
    md << "## File audit\n\n";
    md << "- ch1 bytes: " << audit.ch1_bytes << "\n";
    md << "- ch2 bytes: " << audit.ch2_bytes << "\n";
    md << "- bytes per period per channel: " << audit.bytes_per_old_period_per_channel << "\n";
    md << "- period_count: " << audit.period_count << "\n";
    md << "- status: " << audit.message << "\n";
    return true;
}

bool writeDataStatsReport(const std::string &output_dir,
                          const GenerationStats &stats,
                          const Stage1RunOptions &opt,
                          const Stage1NewSystemConfig &new_cfg,
                          std::string &err)
{
    std::ofstream md(pathJoin(pathJoin(output_dir, "reports"), "data_statistics_report.md").c_str());
    if (!md) {
        err = "failed to open data_statistics_report.md";
        return false;
    }
    md << "# Stage 1 data statistics\n\n";
    md << "- pulse_resample_mode: " << opt.pulse_resample_mode << "\n";
    md << "- range_resize_mode: " << opt.range_resize_mode << "\n";
    md << "- packet order: period -> beam -> pulse -> range -> ch1 IQ + ch2 IQ\n";
    md << "- packets_written: " << stats.packets_written << "\n";
    md << "- output_bytes: " << stats.output_bytes << "\n";
    const std::size_t channel_count =
        static_cast<std::size_t>(std::max(2, new_cfg.new_protocol_channel_count));
    const std::string iq_type = new_cfg.iq_data_type.empty() ? "float32" : new_cfg.iq_data_type;
    const std::size_t packet_bytes =
        gmti::new_protocol::packetBytes(static_cast<std::size_t>(new_cfg.ddc_len_new),
                                        channel_count,
                                        iq_type);
    md << "- expected bytes per full period: " << static_cast<uint64_t>(new_cfg.beam_count) *
              static_cast<uint64_t>(new_cfg.pulse_num_new) *
              static_cast<uint64_t>(packet_bytes) << "\n";
    md << "- iq_data_type: " << iq_type << "\n";
    md << "- new_protocol_channel_count: " << channel_count << "\n";
    md << "- read channels: " << new_cfg.new_protocol_read_channel_1
       << ", " << new_cfg.new_protocol_read_channel_2 << "\n";
    md << "- ch1 mean_abs/rms/max_abs: " << stats.mean_abs_ch1 << " / " << stats.rms_ch1 << " / " << stats.max_abs_ch1 << "\n";
    md << "- ch2 mean_abs/rms/max_abs: " << stats.mean_abs_ch2 << " / " << stats.rms_ch2 << " / " << stats.max_abs_ch2 << "\n";
    md << "- has_nan/has_inf: " << boolText(stats.has_nan) << " / " << boolText(stats.has_inf) << "\n";
    return true;
}

bool writeAlgorithmIntegrationPlaceholder(const std::string &output_dir,
                                          const Stage1RunOptions &opt,
                                          const Stage1NewSystemConfig &new_cfg,
                                          const GenerationStats &stats,
                                          std::string &err)
{
    std::ofstream md(pathJoin(pathJoin(output_dir, "reports"), "algorithm_integration_report.md").c_str());
    if (!md) {
        err = "failed to open algorithm_integration_report.md";
        return false;
    }
    md << "# Stage 1 algorithm integration report\n\n";
    md << "- angle grid: " << new_cfg.scan_min_deg << " to " << new_cfg.scan_max_deg << "\n";
    md << "- pulse_resample_mode: " << opt.pulse_resample_mode << "\n";
    md << "- range_resize_mode: " << opt.range_resize_mode << "\n";
    md << "- input period_count generated: " << opt.period_count << "\n";
    md << "- output file size: " << stats.output_bytes << "\n";
    md << "- suggested smoke command: `./build/GMTI_core " << pathJoin(pathJoin(opt.output_dir, "config"), "temp_config_stage1_newsystem.xml") << "`\n";
    md << "- POS/velocity header: generated packets write UTC, interpolated lat/lon/alt, and vn/ve/vd derived from adjacent POS position differences. This avoids zero platform speed in squint estimation.\n";
    md << "- current status: data generation and static validation complete; runtime algorithm logs should be appended after running GMTI_core or GMTI_pipe_core.\n";
    md << "- stage conclusion: no cooperative target is injected in stage 1; out-of-range background beams are reused by forward cyclic angle wrapping for engineering link and statistics validation, not strict wide-angle geometry proof.\n";
    return true;
}

} // namespace sim_stage1
} // namespace gmti
