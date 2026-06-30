#pragma once

#include "stage1_config.h"

#include <complex>
#include <cstdint>
#include <string>
#include <vector>

namespace gmti {
namespace sim_stage1 {

struct OldBlock {
    std::vector<std::complex<float> > samples;
    std::vector<double> utc;
    std::vector<double> angle_deg;
};

struct FileAudit {
    bool ch1_exists = false;
    bool ch2_exists = false;
    bool ch1_nonzero = false;
    bool ch2_nonzero = false;
    uint64_t ch1_bytes = 0;
    uint64_t ch2_bytes = 0;
    uint64_t bytes_per_old_beam_per_channel = 0;
    uint64_t bytes_per_old_period_per_channel = 0;
    int period_count_ch1 = 0;
    int period_count_ch2 = 0;
    int period_count = 0;
    bool ok = false;
    std::string complex_format = "float32 IQ with 256-byte PRT header per channel";
    std::string message;
};

bool auditOldFiles(const Stage1OldSystemConfig &cfg, FileAudit &audit);
bool readOldBlock(const Stage1OldSystemConfig &cfg,
                  int period_id,
                  int old_beam_index,
                  int channel_id,
                  OldBlock &block,
                  std::string &err);

} // namespace sim_stage1
} // namespace gmti

