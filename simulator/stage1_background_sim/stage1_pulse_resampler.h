#pragma once

#include "stage1_config.h"

#include <string>
#include <vector>

namespace gmti {
namespace sim_stage1 {

struct PulseMapEntry {
    int new_pulse_index = 0;
    double old_pulse_pos = 0.0;
    int old_left_index = 0;
    int old_right_index = 0;
    double weight = 0.0;
};

bool buildPulseMap(const Stage1OldSystemConfig &old_cfg,
                   const Stage1NewSystemConfig &new_cfg,
                   const std::string &mode,
                   std::vector<PulseMapEntry> &map,
                   std::string &err);
bool writePulseMap(const std::string &output_dir,
                   const std::string &mode,
                   const std::vector<PulseMapEntry> &map,
                   std::string &err);

} // namespace sim_stage1
} // namespace gmti

