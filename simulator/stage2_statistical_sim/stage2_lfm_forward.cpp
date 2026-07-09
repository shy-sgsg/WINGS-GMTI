#include "stage2_lfm_forward.h"

#include "../common/SimulationGeometry.h"
#include "../target_injection/beam_visibility.h"
#include "../target_injection/channel_phase_model.h"
#include "../target_injection/radar_geometry.h"
#include "dbs/NewProtocolLayout.hpp"

#include <algorithm>
#include <cmath>
#include <complex>

namespace gmti {
namespace stage2 {

using gmti::target_injection::GeometrySample;
using gmti::target_injection::PlatformState;
using gmti::target_injection::TargetGlobalConfig;
using gmti::target_injection::TargetState;
using gmti::target_injection::Vec3;
using gmti::target_injection::deg2rad;
using gmti::target_injection::dot;
using gmti::target_injection::evaluatePlatformState;
using gmti::target_injection::evaluateVisibility;
using gmti::target_injection::kC;
using gmti::target_injection::kPi;
using gmti::target_injection::norm;
using gmti::target_injection::pulseTimeSec;
using gmti::target_injection::rad2deg;
using gmti::target_injection::wrapTo180;

namespace {

const size_t kHeaderBytes = gmti::new_protocol::kHeaderBytes;

std::complex<float> loadCh(const std::vector<uint8_t> &packet,
                           int n,
                           int ch,
                           std::size_t channel_count,
                           const std::string &iq_type)
{
    const std::size_t iq_bytes = gmti::new_protocol::bytesPerIq(iq_type);
    const size_t off = kHeaderBytes +
                       static_cast<size_t>(n) * gmti::new_protocol::sampleBytes(channel_count, iq_type) +
                       gmti::new_protocol::channelOffset(static_cast<size_t>(ch), iq_type);
    return std::complex<float>(
        gmti::new_protocol::loadIqAsFloat(&packet[off], iq_type),
        gmti::new_protocol::loadIqAsFloat(&packet[off + iq_bytes], iq_type));
}

void storeCh(std::vector<uint8_t> &packet,
             int n,
             int ch,
             std::size_t channel_count,
             const std::string &iq_type,
             const std::complex<float> &v)
{
    const std::size_t iq_bytes = gmti::new_protocol::bytesPerIq(iq_type);
    const size_t off = kHeaderBytes +
                       static_cast<size_t>(n) * gmti::new_protocol::sampleBytes(channel_count, iq_type) +
                       gmti::new_protocol::channelOffset(static_cast<size_t>(ch), iq_type);
    gmti::new_protocol::storeIqFromFloat(&packet[off], iq_type, v.real());
    gmti::new_protocol::storeIqFromFloat(&packet[off + iq_bytes], iq_type, v.imag());
}

GeometrySample scattererGeometry(const gmti::target_injection::RadarConfig &radar,
                                 const TargetGlobalConfig &global,
                                 const Scatterer &s,
                                 int period_id,
                                 int beam_id,
                                 int pulse_id)
{
    GeometrySample g;
    g.time_sec = pulseTimeSec(radar, period_id, beam_id, pulse_id);
    g.theta_cmd_deg = radar.scan_min_deg + radar.scan_step_deg * static_cast<double>(beam_id);
    g.theta_true_deg = g.theta_cmd_deg;
    g.platform = evaluatePlatformState(global, g.time_sec);
    g.target.position = s.position;
    g.target.velocity = Vec3(0.0, 0.0, 0.0);
    g.los = g.target.position - g.platform.position;
    g.range_m = norm(g.los);
    if (g.range_m > 0.0) g.los_unit = g.los * (1.0 / g.range_m);
    g.tau_abs_sec = 2.0 * g.range_m / kC;
    g.tau_rel_sec = g.tau_abs_sec - radar.sample_delay_sec;
    g.range_sample_float = g.tau_rel_sec * radar.fs_hz;
    g.range_sample_int = static_cast<int>(std::floor(g.range_sample_float + 0.5));
    g.in_range_window = (g.range_sample_float >= 0.0 &&
                         g.range_sample_float < static_cast<double>(radar.pulse_len));
    const gmti::sim_geometry::ENUPoint platform_enu =
        gmti::sim_geometry::localToEnu(
            gmti::sim_geometry::LocalPoint(g.platform.position.x, g.platform.position.y, g.platform.position.z),
            global.geometry);
    const gmti::sim_geometry::ENUPoint target_enu =
        gmti::sim_geometry::localToEnu(
            gmti::sim_geometry::LocalPoint(g.target.position.x, g.target.position.y, g.target.position.z),
            global.geometry);
    const gmti::sim_geometry::ENUVelocity platform_vel =
        gmti::sim_geometry::localVelocityToEnu(
            gmti::sim_geometry::LocalVelocity(g.platform.velocity.x, g.platform.velocity.y, g.platform.velocity.z),
            global.geometry);
    const double de = target_enu.e - platform_enu.e;
    const double dn = target_enu.n - platform_enu.n;
    const double horizontal = std::sqrt(de * de + dn * dn);
    const gmti::sim_geometry::LookVectorEN look =
        gmti::sim_geometry::makeAlgorithmLookVectorEN(platform_vel.ve, platform_vel.vn, g.theta_true_deg, global.geometry);
    if (horizontal > 1.0e-9) {
        const double ue = de / horizontal;
        const double un = dn / horizontal;
        const double dot_en = std::max(-1.0, std::min(1.0, look.east * ue + look.north * un));
        const double cross_en = look.east * un - look.north * ue;
        g.angle_error_deg = rad2deg(std::atan2(cross_en, dot_en));
    } else {
        g.angle_error_deg = 0.0;
    }
    g.target_azimuth_deg = wrapTo180(g.theta_true_deg + g.angle_error_deg);
    const Vec3 rel_v = g.target.velocity - g.platform.velocity;
    g.radial_velocity_mps = dot(rel_v, g.los_unit);
    return g;
}

double channelPhaseRad(const gmti::target_injection::RadarConfig &radar, const GeometrySample &g)
{
    const double lambda = kC / radar.fc_hz;
    return 2.0 * kPi * radar.d_chan_m * std::sin(deg2rad(g.theta_true_deg)) / lambda;
}

uint64_t mix64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

double hashUnit(uint64_t seed,
                int period_id,
                int beam_id,
                int pulse_id,
                int sample_id,
                uint64_t salt)
{
    uint64_t x = seed;
    x ^= mix64(static_cast<uint64_t>(period_id) + 0x100000001b3ULL * salt);
    x ^= mix64(static_cast<uint64_t>(beam_id) + 0x9e3779b9ULL * salt);
    x ^= mix64(static_cast<uint64_t>(pulse_id) + 0x7f4a7c15ULL * salt);
    x ^= mix64(static_cast<uint64_t>(sample_id) + 0xbf58476dULL * salt);
    const uint64_t y = mix64(x);
    return (static_cast<double>((y >> 11) & ((1ULL << 53) - 1)) + 0.5) *
           (1.0 / static_cast<double>(1ULL << 53));
}

double hashGaussian(uint64_t seed,
                    int period_id,
                    int beam_id,
                    int pulse_id,
                    int sample_id,
                    uint64_t salt0,
                    uint64_t salt1)
{
    const double u1 = std::max(1.0e-12, hashUnit(seed, period_id, beam_id, pulse_id, sample_id, salt0));
    const double u2 = hashUnit(seed, period_id, beam_id, pulse_id, sample_id, salt1);
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
}

} // namespace

void addScatterersToPacket(std::vector<uint8_t> &packet,
                           const gmti::target_injection::RadarConfig &radar,
                           const TargetGlobalConfig &global,
                           const ScattererList &scatterers,
                           int period_id,
                           int beam_id,
                           int pulse_id,
                           double beam_gain_threshold,
                           Stage2Stats &stats)
{
    const double kr = radar.br_hz / radar.tr_sec;
    const double lambda = kC / radar.fc_hz;
    const int np = static_cast<int>(std::floor(radar.tr_sec * radar.fs_hz + 0.5));
    const std::size_t channel_count = static_cast<std::size_t>(std::max(2, radar.new_protocol_channel_count));
    const std::string iq_type = radar.iq_data_type.empty() ? "float32" : radar.iq_data_type;
    const int ch1 = radar.new_protocol_read_channel_1;
    const int ch2 = radar.new_protocol_read_channel_2;
    for (size_t i = 0; i < scatterers.size(); ++i) {
        const Scatterer &s = scatterers[i];
        const GeometrySample g = scattererGeometry(radar, global, s, period_id, beam_id, pulse_id);
        if (!g.in_range_window) continue;
        const gmti::target_injection::VisibilityResult vis =
            evaluateVisibility(radar, global, g.angle_error_deg);
        if (!vis.visible || vis.beam_gain < beam_gain_threshold) continue;
        const int n_start = static_cast<int>(std::floor(g.range_sample_float));
        const int n_end = std::min(radar.pulse_len, n_start + np);
        if (n_end <= 0 || n_start >= radar.pulse_len) continue;
        const std::complex<double> ch2_phase = std::exp(std::complex<double>(0.0, channelPhaseRad(radar, g)));
        const double carrier_phase =
            static_cast<double>(global.carrier_phase_sign) * 4.0 * kPi * g.range_m / lambda + s.phase_rad;
        const std::complex<double> carrier = std::exp(std::complex<double>(0.0, carrier_phase));
        int touched = 0;
        for (int n = std::max(0, n_start); n < n_end; ++n) {
            const double dt = static_cast<double>(n) / radar.fs_hz - g.tau_rel_sec;
            if (dt < 0.0 || dt >= radar.tr_sec) continue;
            const double chirp_phase =
                static_cast<double>(global.chirp_phase_sign) * kPi * kr * dt * dt;
            const std::complex<double> echo =
                s.amplitude * vis.beam_gain *
                std::exp(std::complex<double>(0.0, chirp_phase)) * carrier;
            const std::complex<float> e1(static_cast<float>(echo.real()), static_cast<float>(echo.imag()));
            const std::complex<double> echo2 = echo * ch2_phase;
            const std::complex<float> e2(static_cast<float>(echo2.real()), static_cast<float>(echo2.imag()));
            storeCh(packet, n, ch1, channel_count, iq_type, loadCh(packet, n, ch1, channel_count, iq_type) + e1);
            storeCh(packet, n, ch2, channel_count, iq_type, loadCh(packet, n, ch2, channel_count, iq_type) + e2);
            ++touched;
        }
        if (touched > 0) {
            ++stats.scatterer_echoes;
            stats.scatterer_samples += static_cast<uint64_t>(touched);
        }
    }
}

void addContinuousAreaClutter(std::vector<uint8_t> &packet,
                              const gmti::target_injection::RadarConfig &radar,
                              const TargetGlobalConfig &global,
                              const SceneConfig &scene,
                              uint32_t random_seed,
                              int period_id,
                              int beam_id,
                              int pulse_id,
                              Stage2Stats &stats)
{
    if (!scene.area.enabled) return;
    const bool continuous_mode =
        scene.area.model == "continuous_texture" ||
        scene.area.model == "continuous_surface" ||
        scene.area.model == "continuous_grid" ||
        scene.area.model == "grid_texture";
    if (!continuous_mode) return;

    const std::size_t channel_count = static_cast<std::size_t>(std::max(2, radar.new_protocol_channel_count));
    const std::string iq_type = radar.iq_data_type.empty() ? "float32" : radar.iq_data_type;
    const int ch1 = radar.new_protocol_read_channel_1;
    const int ch2 = radar.new_protocol_read_channel_2;
    const double lambda = kC / radar.fc_hz;
    const double theta_cmd_deg =
        radar.scan_min_deg + radar.scan_step_deg * static_cast<double>(beam_id);
    const double ch2_phase_rad = 2.0 * kPi * radar.d_chan_m * std::sin(deg2rad(theta_cmd_deg)) / lambda;
    const std::complex<float> ch2_phase(
        static_cast<float>(std::cos(ch2_phase_rad)),
        static_cast<float>(std::sin(ch2_phase_rad)));
    const double beam_gain = 1.0;
    const double sigma = std::sqrt(std::max(1.0e-12, scene.area.mean_power) / 2.0);
    const double t = pulseTimeSec(radar, period_id, beam_id, pulse_id);
    const PlatformState platform = evaluatePlatformState(global, t);
    (void)platform;

    const uint64_t seed = static_cast<uint64_t>(random_seed);
    int injected = 0;
    for (int n = 0; n < radar.pulse_len; ++n) {
        const double range_m =
            0.5 * kC * (radar.sample_delay_sec + static_cast<double>(n) / radar.fs_hz);
        if (range_m < scene.range_min_m || range_m > scene.range_max_m) {
            continue;
        }
        const double u_phase = hashUnit(seed, period_id, beam_id, pulse_id, n, 11ULL);
        const double g0 = hashGaussian(seed, period_id, beam_id, pulse_id, n, 23ULL, 29ULL);
        const double g1 = hashGaussian(seed, period_id, beam_id, pulse_id, n, 31ULL, 37ULL);
        double amp = sigma * std::sqrt(-2.0 * std::log(std::max(1.0e-12, u_phase)));
        if (scene.area.texture_sigma > 0.0) {
            amp *= std::exp(scene.area.texture_sigma * g0);
        }
        amp *= scene.clutter_amplitude_scale * beam_gain;
        const double carrier_phase = static_cast<double>(global.carrier_phase_sign) * 4.0 * kPi * range_m / lambda;
        const double clutter_phase = 2.0 * kPi * hashUnit(seed, period_id, beam_id, pulse_id, n, 41ULL) + g1 * 0.15;
        const std::complex<double> echo1 =
            amp * std::exp(std::complex<double>(0.0, carrier_phase + clutter_phase));
        const std::complex<double> echo2 = echo1 * std::complex<double>(ch2_phase.real(), ch2_phase.imag());
        const std::complex<float> e1(static_cast<float>(echo1.real()), static_cast<float>(echo1.imag()));
        const std::complex<float> e2(static_cast<float>(echo2.real()), static_cast<float>(echo2.imag()));
        storeCh(packet, n, ch1, channel_count, iq_type, loadCh(packet, n, ch1, channel_count, iq_type) + e1);
        storeCh(packet, n, ch2, channel_count, iq_type, loadCh(packet, n, ch2, channel_count, iq_type) + e2);
        ++injected;
    }
    if (injected > 0) {
        ++stats.continuous_area_packets;
        stats.continuous_area_samples += static_cast<uint64_t>(injected);
    }
}

void addThermalNoise(std::vector<uint8_t> &packet,
                     const gmti::target_injection::RadarConfig &radar,
                     double noise_power,
                     std::mt19937 &rng,
                     Stage2Stats &stats)
{
    if (noise_power <= 0.0) return;
    const double sigma = std::sqrt(noise_power / 2.0);
    std::normal_distribution<float> n01(0.0f, static_cast<float>(sigma));
    const std::size_t channel_count = static_cast<std::size_t>(std::max(2, radar.new_protocol_channel_count));
    const std::string iq_type = radar.iq_data_type.empty() ? "float32" : radar.iq_data_type;
    const int ch1 = radar.new_protocol_read_channel_1;
    const int ch2 = radar.new_protocol_read_channel_2;
    for (int n = 0; n < radar.pulse_len; ++n) {
        const std::complex<float> z1(n01(rng), n01(rng));
        const std::complex<float> z2(n01(rng), n01(rng));
        storeCh(packet, n, ch1, channel_count, iq_type, loadCh(packet, n, ch1, channel_count, iq_type) + z1);
        storeCh(packet, n, ch2, channel_count, iq_type, loadCh(packet, n, ch2, channel_count, iq_type) + z2);
        stats.sum_noise_power += static_cast<double>(std::norm(z1) + std::norm(z2));
        stats.noise_samples += 2U;
    }
}

void scanPacketStats(const std::vector<uint8_t> &packet,
                     const gmti::target_injection::RadarConfig &radar,
                     Stage2Stats &stats)
{
    const std::size_t channel_count = static_cast<std::size_t>(std::max(2, radar.new_protocol_channel_count));
    const std::string iq_type = radar.iq_data_type.empty() ? "float32" : radar.iq_data_type;
    const int ch1 = radar.new_protocol_read_channel_1;
    const int ch2 = radar.new_protocol_read_channel_2;
    for (int n = 0; n < radar.pulse_len; ++n) {
        const std::complex<float> c1 = loadCh(packet, n, ch1, channel_count, iq_type);
        const std::complex<float> c2 = loadCh(packet, n, ch2, channel_count, iq_type);
        const float vals[4] = {c1.real(), c1.imag(), c2.real(), c2.imag()};
        for (int k = 0; k < 4; ++k) {
            if (std::isnan(vals[k])) stats.has_nan = true;
            if (std::isinf(vals[k])) stats.has_inf = true;
            stats.max_abs_component = std::max(stats.max_abs_component, std::fabs(static_cast<double>(vals[k])));
        }
    }
}

} // namespace stage2
} // namespace gmti
