#pragma once

#include "stage1_beam_mapper.h"
#include "stage1_config.h"
#include "stage1_data_reader.h"
#include "stage1_pulse_resampler.h"
#include "stage1_writer.h"

#include <string>
#include <vector>

namespace gmti {
namespace sim_stage1 {

bool writeParameterAudit(const Stage1OldSystemConfig &old_cfg,
                         const Stage1NewSystemConfig &new_cfg,
                         const FileAudit &audit,
                         const std::string &output_dir,
                         std::string &err);

bool writeDataStatsReport(const std::string &output_dir,
                          const GenerationStats &stats,
                          const Stage1RunOptions &opt,
                          const Stage1NewSystemConfig &new_cfg,
                          std::string &err);

bool writeAlgorithmIntegrationPlaceholder(const std::string &output_dir,
                                          const Stage1RunOptions &opt,
                                          const Stage1NewSystemConfig &new_cfg,
                                          const GenerationStats &stats,
                                          std::string &err);

} // namespace sim_stage1
} // namespace gmti

