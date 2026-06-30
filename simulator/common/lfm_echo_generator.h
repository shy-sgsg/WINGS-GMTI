#pragma once

#include <complex>

namespace gmti {
namespace sim_common {

std::complex<float> lfmSample(double amplitude,
                              double beam_gain,
                              double phase_rad,
                              double dt_sec,
                              double pulse_width_sec,
                              double chirp_rate_hz_per_sec,
                              int chirp_phase_sign);

} // namespace sim_common
} // namespace gmti
