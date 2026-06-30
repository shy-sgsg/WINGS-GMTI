#pragma once

#include "stage1_config.h"

#include <string>
#include <vector>

namespace gmti {
namespace sim_stage1 {

struct BeamMapEntry {
    int new_beam_index = 0;
    double theta_new_deg = 0.0;
    double background_ref_angle_deg = 0.0;
    std::string source_mode;
    int source_left_beam_index = 0;
    double source_left_angle_deg = 0.0;
    int source_right_beam_index = 0;
    double source_right_angle_deg = 0.0;
    double interp_weight = 0.0;
    bool is_wrapped = false;
};

std::vector<double> buildOldAngles();
std::vector<double> buildNewAngles(const Stage1NewSystemConfig &cfg);
bool buildBeamMap(const Stage1NewSystemConfig &cfg, std::vector<BeamMapEntry> &map, std::string &err);
bool writeBeamTables(const std::string &output_dir, const std::vector<BeamMapEntry> &map, std::string &err);

} // namespace sim_stage1
} // namespace gmti
