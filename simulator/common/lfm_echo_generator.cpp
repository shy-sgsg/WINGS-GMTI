#include "lfm_echo_generator.h"

#include <cmath>

namespace gmti {
namespace sim_common {

namespace {
const double kPi = 3.14159265358979323846;
}

std::complex<float> lfmSample(double amplitude,
                              double beam_gain,
                              double phase_rad,
                              double dt_sec,
                              double pulse_width_sec,
                              double chirp_rate_hz_per_sec,
                              int chirp_phase_sign)
{
    if (dt_sec < 0.0 || dt_sec >= pulse_width_sec) return std::complex<float>(0.0f, 0.0f);
    const double phase = phase_rad + static_cast<double>(chirp_phase_sign) * kPi * chirp_rate_hz_per_sec * dt_sec * dt_sec;
    return std::complex<float>(static_cast<float>(amplitude * beam_gain * std::cos(phase)),
                               static_cast<float>(amplitude * beam_gain * std::sin(phase)));
}

} // namespace sim_common
} // namespace gmti
