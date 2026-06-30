#pragma once

#include "lfm_echo_generator.h"

#include <fstream>
#include <map>

namespace gmti {
namespace target_injection {

struct BeamSummaryAccumulator {
    int period_id = 0;
    int beam_id = 0;
    int target_id = 0;
    int rows = 0;
    int visible_pulse_count = 0;
    double sum_range_m = 0.0;
    double sum_range_sample = 0.0;
    double sum_target_azimuth_deg = 0.0;
    double sum_beam_gain = 0.0;
    double sum_radial_velocity_mps = 0.0;
    double sum_target_amplitude = 0.0;
    int injected_sample_count = 0;
};

class TruthWriter {
public:
    bool open(const std::string &truth_dir, std::string &err);
    void writePulse(const PulseTruth &t);
    void writeSummary();
    void close();

private:
    std::ofstream pulse_;
    std::ofstream summary_;
    std::map<std::pair<int, int>, BeamSummaryAccumulator> acc_;
};

bool writeTargetInjectionReport(const std::string &path,
                                const InjectionConfig &cfg,
                                const InjectionStats &stats,
                                const std::string &notes,
                                std::string &err);

} // namespace target_injection
} // namespace gmti

