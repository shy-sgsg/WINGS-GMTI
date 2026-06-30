#include "stage1_data_reader.h"

#include <cstdio>
#include <fstream>
#include <sstream>

namespace gmti {
namespace sim_stage1 {

namespace {

uint64_t fileSize(const std::string &path, bool &exists)
{
    std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
    exists = static_cast<bool>(in);
    return exists ? static_cast<uint64_t>(in.tellg()) : 0U;
}

static inline int16_t load_i16_le(const uint8_t *p)
{
    return static_cast<int16_t>(p[0] | (static_cast<uint16_t>(p[1]) << 8));
}

static inline uint32_t load_u32_le(const uint8_t *p)
{
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

} // namespace

bool auditOldFiles(const Stage1OldSystemConfig &cfg, FileAudit &audit)
{
    audit = FileAudit();
    audit.ch1_bytes = fileSize(cfg.data_ch1, audit.ch1_exists);
    audit.ch2_bytes = fileSize(cfg.data_ch2, audit.ch2_exists);
    audit.ch1_nonzero = audit.ch1_bytes > 0U;
    audit.ch2_nonzero = audit.ch2_bytes > 0U;
    const uint64_t bytes_iq = static_cast<uint64_t>(cfg.pulse_len) * 2U * sizeof(float);
    const uint64_t bytes_prt = static_cast<uint64_t>(cfg.info_len) + bytes_iq;
    audit.bytes_per_old_beam_per_channel = bytes_prt * static_cast<uint64_t>(cfg.pulse_num);
    audit.bytes_per_old_period_per_channel =
        audit.bytes_per_old_beam_per_channel * static_cast<uint64_t>(cfg.az_count);
    if (!audit.ch1_exists || !audit.ch2_exists || !audit.ch1_nonzero || !audit.ch2_nonzero) {
        audit.message = "old channel path missing or empty";
        return false;
    }
    if (audit.bytes_per_old_period_per_channel == 0U ||
        (audit.ch1_bytes % audit.bytes_per_old_period_per_channel) != 0U ||
        (audit.ch2_bytes % audit.bytes_per_old_period_per_channel) != 0U) {
        audit.message = "old file size is not an integer number of periods";
        return false;
    }
    audit.period_count_ch1 =
        static_cast<int>(audit.ch1_bytes / audit.bytes_per_old_period_per_channel);
    audit.period_count_ch2 =
        static_cast<int>(audit.ch2_bytes / audit.bytes_per_old_period_per_channel);
    if (audit.period_count_ch1 != audit.period_count_ch2) {
        audit.message = "two old channels have different period counts";
        return false;
    }
    audit.period_count = audit.period_count_ch1;
    audit.ok = true;
    audit.message = "ok";
    return true;
}

bool readOldBlock(const Stage1OldSystemConfig &cfg,
                  int period_id,
                  int old_beam_index,
                  int channel_id,
                  OldBlock &block,
                  std::string &err)
{
    if (old_beam_index < 0 || old_beam_index >= cfg.az_count) {
        err = "old beam index out of range";
        return false;
    }
    const std::string &path = (channel_id == 2) ? cfg.data_ch2 : cfg.data_ch1;
    FILE *fp = std::fopen(path.c_str(), "rb");
    if (!fp) {
        err = "failed to open old echo file: " + path;
        return false;
    }
    const uint64_t bytes_iq = static_cast<uint64_t>(cfg.pulse_len) * 2U * sizeof(float);
    const uint64_t bytes_prt = static_cast<uint64_t>(cfg.info_len) + bytes_iq;
    const uint64_t beam_index_global =
        static_cast<uint64_t>(period_id) * static_cast<uint64_t>(cfg.az_count) +
        static_cast<uint64_t>(old_beam_index);
    const uint64_t off = beam_index_global * bytes_prt * static_cast<uint64_t>(cfg.pulse_num);
    if (std::fseek(fp, static_cast<long>(off), SEEK_SET) != 0) {
        std::fclose(fp);
        err = "fseek failed for old echo block";
        return false;
    }
    block.samples.assign(static_cast<size_t>(cfg.pulse_num) * static_cast<size_t>(cfg.pulse_len),
                         std::complex<float>(0.0f, 0.0f));
    block.utc.assign(static_cast<size_t>(cfg.pulse_num), 0.0);
    block.angle_deg.assign(static_cast<size_t>(cfg.pulse_num), 0.0);
    std::vector<uint8_t> header(static_cast<size_t>(cfg.info_len));
    std::vector<float> iq(static_cast<size_t>(cfg.pulse_len) * 2U);
    for (int p = 0; p < cfg.pulse_num; ++p) {
        if (std::fread(header.data(), 1, header.size(), fp) != header.size()) {
            std::fclose(fp);
            err = "short read in old PRT header";
            return false;
        }
        if (header.size() > 219U) {
            block.angle_deg[static_cast<size_t>(p)] =
                static_cast<double>(load_i16_le(&header[218])) / 100.0;
        }
        if (header.size() >= 44U) {
            auto bcd = [](uint8_t x) { return static_cast<int>(x) - 6 * (static_cast<int>(x) / 16); };
            const int hh = bcd(header[36]);
            const int mm = bcd(header[37]);
            const int ss = bcd(header[38]);
            block.utc[static_cast<size_t>(p)] =
                hh * 3600.0 + mm * 60.0 + ss +
                static_cast<double>(load_u32_le(&header[39])) / 100e6 +
                cfg.week_offset * 24.0 * 3600.0 + cfg.sec_bias;
        }
        if (std::fread(iq.data(), sizeof(float), iq.size(), fp) != iq.size()) {
            std::fclose(fp);
            err = "short read in old IQ payload";
            return false;
        }
        std::complex<float> *row =
            &block.samples[static_cast<size_t>(p) * static_cast<size_t>(cfg.pulse_len)];
        for (int n = 0; n < cfg.pulse_len; ++n) {
            row[n] = std::complex<float>(iq[static_cast<size_t>(2 * n)],
                                         iq[static_cast<size_t>(2 * n + 1)]);
        }
    }
    std::fclose(fp);
    return true;
}

} // namespace sim_stage1
} // namespace gmti

