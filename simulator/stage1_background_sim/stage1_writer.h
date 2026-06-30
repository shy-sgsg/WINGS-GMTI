#pragma once

#include "stage1_beam_mapper.h"
#include "stage1_config.h"
#include "stage1_data_reader.h"
#include "stage1_pulse_resampler.h"
#include "stage1_range_resizer.h"

#include <complex>
#include <string>
#include <vector>

namespace gmti {
namespace sim_stage1 {

struct GenerationStats {
    uint64_t packets_written = 0;
    uint64_t output_bytes = 0;
    double mean_abs_ch1 = 0.0;
    double rms_ch1 = 0.0;
    double max_abs_ch1 = 0.0;
    double mean_abs_ch2 = 0.0;
    double rms_ch2 = 0.0;
    double max_abs_ch2 = 0.0;
    bool has_nan = false;
    bool has_inf = false;
};

bool generateStage1Data(const Stage1OldSystemConfig &old_cfg,
                        const Stage1NewSystemConfig &new_cfg,
                        const Stage1RunOptions &opt,
                        const std::vector<BeamMapEntry> &beam_map,
                        const std::vector<PulseMapEntry> &pulse_map,
                        RangeFftZeroPadResizer &resizer,
                        const std::string &out_file,
                        GenerationStats &stats,
                        std::string &err);

bool writeStage1ConfigXml(const Stage1OldSystemConfig &old_cfg,
                          const Stage1NewSystemConfig &new_cfg,
                          const Stage1RunOptions &opt,
                          const std::string &data_file,
                          std::string &err);

} // namespace sim_stage1
} // namespace gmti

