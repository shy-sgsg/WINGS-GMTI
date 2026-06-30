#pragma once

#include "beam_visibility.h"
#include "channel_phase_model.h"

#include <vector>

namespace gmti {
namespace target_injection {

struct PulseTruth {
    int period_id = 0;
    int beam_id = 0;
    int pulse_id = 0;
    int target_id = 0;
    std::string target_name;
    GeometrySample geom;
    double beam_gain = 0.0;
    bool visible_by_beam = false;
    bool target_period_enabled = false;
    bool injection_enabled = false;
    double local_background_rms = 0.0;
    double target_amplitude = 0.0;
    double delta_phi_ch_rad = 0.0;
    int injected_sample_count = 0;
};

struct InjectionStats {
    uint64_t packets_read = 0;
    uint64_t packets_written = 0;
    uint64_t pulses_injected = 0;
    uint64_t samples_injected = 0;
    double max_amplitude = 0.0;
    bool has_nan = false;
    bool has_inf = false;
};

float loadF32LE(const uint8_t *p);
void storeF32LE(uint8_t *p, float v);
void fillZeroPacketHeader(std::vector<uint8_t> &packet,
                          const RadarConfig &radar,
                          const TargetGlobalConfig &global,
                          uint32_t prt_counter,
                          double utc,
                          double theta_deg);

PulseTruth injectOnePulse(std::vector<uint8_t> &packet,
                          const RadarConfig &radar,
                          const TargetGlobalConfig &global,
                          const TargetConfig &target,
                          int period_id,
                          int beam_id,
                          int pulse_id);

} // namespace target_injection
} // namespace gmti

