#include "lfm_echo_generator.h"

#include <algorithm>
#include <cstring>

namespace gmti {
namespace target_injection {

namespace {

const size_t kHeaderBytes = 256;
const size_t kBytesPerSample = 16;

void storeU16LE(std::vector<uint8_t> &b, size_t off, uint16_t v)
{
    b[off] = static_cast<uint8_t>(v & 0xffU);
    b[off + 1] = static_cast<uint8_t>((v >> 8) & 0xffU);
}

void storeI16LE(std::vector<uint8_t> &b, size_t off, int16_t v)
{
    storeU16LE(b, off, static_cast<uint16_t>(v));
}

void storeU32LE(std::vector<uint8_t> &b, size_t off, uint32_t v)
{
    for (int i = 0; i < 4; ++i) b[off + static_cast<size_t>(i)] = static_cast<uint8_t>((v >> (8 * i)) & 0xffU);
}

void storeU64LE(std::vector<uint8_t> &b, size_t off, uint64_t v)
{
    for (int i = 0; i < 8; ++i) b[off + static_cast<size_t>(i)] = static_cast<uint8_t>((v >> (8 * i)) & 0xffU);
}

void storeF64LE(std::vector<uint8_t> &b, size_t off, double v)
{
    uint64_t raw = 0;
    std::memcpy(&raw, &v, sizeof(double));
    storeU64LE(b, off, raw);
}

int16_t satI16(double x)
{
    if (x > 32767.0) return 32767;
    if (x < -32768.0) return -32768;
    return static_cast<int16_t>(std::floor(x + (x >= 0.0 ? 0.5 : -0.5)));
}

std::complex<float> loadComplexCh(const std::vector<uint8_t> &packet, int pulse_len, int n, int ch)
{
    (void)pulse_len;
    const size_t off = kHeaderBytes + static_cast<size_t>(n) * kBytesPerSample + (ch == 1 ? 0U : 8U);
    return std::complex<float>(loadF32LE(&packet[off]), loadF32LE(&packet[off + 4]));
}

void storeComplexCh(std::vector<uint8_t> &packet, int pulse_len, int n, int ch, const std::complex<float> &v)
{
    (void)pulse_len;
    const size_t off = kHeaderBytes + static_cast<size_t>(n) * kBytesPerSample + (ch == 1 ? 0U : 8U);
    storeF32LE(&packet[off], v.real());
    storeF32LE(&packet[off + 4], v.imag());
}

double estimateLocalRms(const std::vector<uint8_t> &packet, const RadarConfig &radar, int center, double floor_value)
{
    const int half = 64;
    const int lo = std::max(0, center - half);
    const int hi = std::min(radar.pulse_len - 1, center + half);
    double power = 0.0;
    int count = 0;
    for (int n = lo; n <= hi; ++n) {
        const std::complex<float> c1 = loadComplexCh(packet, radar.pulse_len, n, 1);
        const std::complex<float> c2 = loadComplexCh(packet, radar.pulse_len, n, 2);
        power += static_cast<double>(std::norm(c1) + std::norm(c2)) * 0.5;
        ++count;
    }
    if (count <= 0) return floor_value;
    return std::max(floor_value, std::sqrt(power / static_cast<double>(count)));
}

} // namespace

float loadF32LE(const uint8_t *p)
{
    uint32_t raw = 0;
    raw |= static_cast<uint32_t>(p[0]);
    raw |= static_cast<uint32_t>(p[1]) << 8;
    raw |= static_cast<uint32_t>(p[2]) << 16;
    raw |= static_cast<uint32_t>(p[3]) << 24;
    float v = 0.0f;
    std::memcpy(&v, &raw, sizeof(float));
    return v;
}

void storeF32LE(uint8_t *p, float v)
{
    uint32_t raw = 0;
    std::memcpy(&raw, &v, sizeof(float));
    p[0] = static_cast<uint8_t>(raw & 0xffU);
    p[1] = static_cast<uint8_t>((raw >> 8) & 0xffU);
    p[2] = static_cast<uint8_t>((raw >> 16) & 0xffU);
    p[3] = static_cast<uint8_t>((raw >> 24) & 0xffU);
}

void fillZeroPacketHeader(std::vector<uint8_t> &packet,
                          const RadarConfig &radar,
                          const TargetGlobalConfig &global,
                          uint32_t prt_counter,
                          double utc,
                          double theta_deg)
{
    std::fill(packet.begin(), packet.end(), 0U);
    const uint32_t prt_len = static_cast<uint32_t>(packet.size());
    storeU64LE(packet, 0, 0x5A5A5A5A5A5A5A5AULL);
    packet[8] = 5;
    storeU32LE(packet, 9, prt_len);
    storeF32LE(&packet[16], static_cast<float>(utc));
    storeU32LE(packet, 20, prt_counter);
    packet[88] = 0x02;
    packet[90] = 0x40;
    packet[92] = 0x0B;
    packet[93] = 0x01;
    storeI16LE(packet, 94, 0);
    storeU32LE(packet, 96, static_cast<uint32_t>(std::floor(utc * 1000.0 + 0.5)));
    storeF64LE(packet, 104, 0.0);
    storeF64LE(packet, 112, 0.0);
    storeF64LE(packet, 120, global.platform_height_m);
    storeF32LE(&packet[128], 0.0f);
    storeF32LE(&packet[132], static_cast<float>(global.platform_speed_mps));
    storeF32LE(&packet[136], 0.0f);
    storeF32LE(&packet[140], static_cast<float>(global.platform_speed_mps));
    packet[208] = static_cast<uint8_t>(prt_counter & 0xffU);
    storeI16LE(packet, 218, satI16(theta_deg * 100.0));
    storeU64LE(packet, 248, 0x5B5B5B5B5B5B5B5BULL);
    (void)radar;
}

PulseTruth injectOnePulse(std::vector<uint8_t> &packet,
                          const RadarConfig &radar,
                          const TargetGlobalConfig &global,
                          const TargetConfig &target,
                          int period_id,
                          int beam_id,
                          int pulse_id)
{
    PulseTruth truth;
    truth.period_id = period_id;
    truth.beam_id = beam_id;
    truth.pulse_id = pulse_id;
    truth.target_id = target.id;
    truth.target_name = target.name;
    truth.geom = evaluateGeometry(radar, global, target, period_id, beam_id, pulse_id);
    const VisibilityResult vr = evaluateVisibility(radar, global, truth.geom.angle_error_deg);
    truth.beam_gain = vr.beam_gain;
    truth.visible_by_beam = vr.visible;
    truth.target_period_enabled = target.enabled &&
                                  period_id >= target.start_period &&
                                  period_id <= target.end_period;
    truth.delta_phi_ch_rad = channelPhaseRad(radar, global, truth.geom);
    truth.local_background_rms = estimateLocalRms(packet, radar, truth.geom.range_sample_int, global.rms_floor);
    if (global.amplitude_mode == "direct_amplitude") {
        truth.target_amplitude = global.direct_amplitude;
    } else {
        truth.target_amplitude = truth.local_background_rms * std::pow(10.0, global.target_snr_db / 20.0);
    }
    truth.injection_enabled = truth.target_period_enabled &&
                              truth.geom.in_range_window &&
                              truth.visible_by_beam;
    if (!truth.injection_enabled) return truth;

    const double kr = radar.br_hz / radar.tr_sec;
    const double lambda = kC / radar.fc_hz;
    const int np = static_cast<int>(std::floor(radar.tr_sec * radar.fs_hz + 0.5));
    const int n_start = static_cast<int>(std::floor(truth.geom.range_sample_float));
    const int n_end = std::min(radar.pulse_len, n_start + np);
    const std::complex<double> ch2_phase =
        std::exp(std::complex<double>(0.0, truth.delta_phi_ch_rad));
    const double carrier_phase =
        static_cast<double>(global.carrier_phase_sign) * 4.0 * kPi * truth.geom.range_m / lambda;
    const std::complex<double> carrier = std::exp(std::complex<double>(0.0, carrier_phase));
    for (int n = std::max(0, n_start); n < n_end; ++n) {
        const double t_fast = static_cast<double>(n) / radar.fs_hz;
        const double dt = t_fast - truth.geom.tau_rel_sec;
        if (dt < 0.0 || dt >= radar.tr_sec) continue;
        const double chirp_phase =
            static_cast<double>(global.chirp_phase_sign) * kPi * kr * dt * dt;
        const std::complex<double> echo =
            truth.target_amplitude * truth.beam_gain *
            std::exp(std::complex<double>(0.0, chirp_phase)) * carrier;
        const std::complex<float> e1(static_cast<float>(echo.real()), static_cast<float>(echo.imag()));
        const std::complex<double> echo2 = echo * ch2_phase;
        const std::complex<float> e2(static_cast<float>(echo2.real()), static_cast<float>(echo2.imag()));
        const std::complex<float> c1 = loadComplexCh(packet, radar.pulse_len, n, 1);
        const std::complex<float> c2 = loadComplexCh(packet, radar.pulse_len, n, 2);
        const std::complex<float> o1 = c1 + e1;
        const std::complex<float> o2 = c2 + e2;
        storeComplexCh(packet, radar.pulse_len, n, 1, o1);
        storeComplexCh(packet, radar.pulse_len, n, 2, o2);
        ++truth.injected_sample_count;
    }
    return truth;
}

} // namespace target_injection
} // namespace gmti
