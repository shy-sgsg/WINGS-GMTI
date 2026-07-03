#include "stage2_lfm_forward.h"

#include "../common/SimulationGeometry.h"
#include "../target_injection/beam_visibility.h"
#include "../target_injection/channel_phase_model.h"
#include "../target_injection/radar_geometry.h"

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
using gmti::target_injection::loadF32LE;
using gmti::target_injection::norm;
using gmti::target_injection::pulseTimeSec;
using gmti::target_injection::rad2deg;
using gmti::target_injection::storeF32LE;
using gmti::target_injection::wrapTo180;

namespace {

const size_t kHeaderBytes = 256;
const size_t kBytesPerSample = 16;

std::complex<float> loadCh(const std::vector<uint8_t> &packet, int n, int ch)
{
    const size_t off = kHeaderBytes + static_cast<size_t>(n) * kBytesPerSample + (ch == 1 ? 0U : 8U);
    return std::complex<float>(loadF32LE(&packet[off]), loadF32LE(&packet[off + 4]));
}

void storeCh(std::vector<uint8_t> &packet, int n, int ch, const std::complex<float> &v)
{
    const size_t off = kHeaderBytes + static_cast<size_t>(n) * kBytesPerSample + (ch == 1 ? 0U : 8U);
    storeF32LE(&packet[off], v.real());
    storeF32LE(&packet[off + 4], v.imag());
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
            storeCh(packet, n, 1, loadCh(packet, n, 1) + e1);
            storeCh(packet, n, 2, loadCh(packet, n, 2) + e2);
            ++touched;
        }
        if (touched > 0) {
            ++stats.scatterer_echoes;
            stats.scatterer_samples += static_cast<uint64_t>(touched);
        }
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
    for (int n = 0; n < radar.pulse_len; ++n) {
        const std::complex<float> z1(n01(rng), n01(rng));
        const std::complex<float> z2(n01(rng), n01(rng));
        storeCh(packet, n, 1, loadCh(packet, n, 1) + z1);
        storeCh(packet, n, 2, loadCh(packet, n, 2) + z2);
        stats.sum_noise_power += static_cast<double>(std::norm(z1) + std::norm(z2));
        stats.noise_samples += 2U;
    }
}

void scanPacketStats(const std::vector<uint8_t> &packet,
                     const gmti::target_injection::RadarConfig &radar,
                     Stage2Stats &stats)
{
    for (int n = 0; n < radar.pulse_len; ++n) {
        const std::complex<float> c1 = loadCh(packet, n, 1);
        const std::complex<float> c2 = loadCh(packet, n, 2);
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
