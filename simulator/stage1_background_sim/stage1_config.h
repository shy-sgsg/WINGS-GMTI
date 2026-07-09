#pragma once

#include <cstdint>
#include <string>

namespace gmti {
namespace sim_stage1 {

struct Stage1OldSystemConfig {
    std::string xml_path;
    std::string data_ch1;
    std::string data_ch2;
    std::string pos_path;
    std::string result_add;
    int info_len = 256;
    int pulse_len = 4096;
    int rg_len = 4096;
    int pulse_num = 256;
    int read_pulse_num = 256;
    int read_pulse_offset = -1;
    int pulse_dec = 1;
    int az_count = 26;
    int wavepos_st = 1;
    int wavepos_ed = 26;
    int wavepos_skip = 1;
    int skip_pulses = 0;
    int week_offset = 0;
    int sec_bias = 0;
    double fc_ghz = 17.0;
    double br_mhz = 50.0;
    double fs_mhz = 75.0;
    double tr_us = 38.0;
    double prf_hz = 2000.0;
    double scan_min_deg = -25.0;
    double scan_max_deg = 25.0;
    double d_chan = 0.17;
    double rmin_m = 13950.0;
};

struct Stage1NewSystemConfig {
    double scan_min_deg = -60.0;
    double scan_max_deg = 60.0;
    double scan_step_deg = 2.0;
    int beam_count = 61;
    int pulse_num_new = 130;
    std::string iq_data_type = "float32";
    int new_protocol_channel_count = 2;
    int new_protocol_read_channel_1 = 1;
    int new_protocol_read_channel_2 = 2;
    int ddc_len_new = 11820;
    int fft_len_new = 12288;
    int pc_crop_start = 3864;
    int pc_crop_len = 4096;
    double fc_new_ghz = 16.0;
    double br_new_mhz = 50.0;
    double fs_new_mhz = 60.0;
    double tr_new_us = 130.0;
    double prf_new_hz = 1300.0;
    double sample_delay_us = 488.0;
    double sample_window_us = 197.0;
    double beam_width_deg = 2.28;

    double periodTimeSec() const;
    double realtime80Sec() const;
};

struct Stage1RunOptions {
    std::string old_config_path = "temp_config.xml";
    std::string output_dir = "outputs/stage1";
    std::string pulse_resample_mode = "physical_time";
    std::string range_resize_mode = "fft_zero_pad";
    int period_start = 0;
    int period_count = 1;
    int beam_start = 0;
    int beam_count = 0;
    std::string channel_mode = "both";
    bool write_config = true;
    bool validate = true;
    bool generate_data = true;
};

bool parseOldConfigXml(const std::string &xml_path, Stage1OldSystemConfig &cfg, std::string &err);
bool parseCommandLine(int argc, char **argv, Stage1RunOptions &opt, Stage1NewSystemConfig &new_cfg);
bool ensureStage1Dirs(const std::string &output_dir, std::string &err);
std::string pathJoin(const std::string &a, const std::string &b);
std::string boolText(bool v);

} // namespace sim_stage1
} // namespace gmti
