#ifndef MOTION_COMP_HPP
#define MOTION_COMP_HPP

#include "config_structs.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

inline double wrapPiLocal(double phase)
{
    if (!std::isfinite(phase)) {
        return phase;
    }
    phase = std::fmod(phase + M_PI, 2.0 * M_PI);
    if (phase < 0.0) {
        phase += 2.0 * M_PI;
    }
    return phase - M_PI;
}

inline double staticPhaseTheoryRad(const Config &cfg, double theta_true_deg, double lambda)
{
    if (!(lambda > 0.0) && cfg.fc > 0.0) {
        lambda = C / (cfg.fc * 1.0e9);
    }
    if (!(lambda > 0.0) || !(cfg.d_channel > 0.0) || !std::isfinite(theta_true_deg)) {
        return cfg.ati_phase_bias_rad;
    }
    return 2.0 * M_PI * cfg.d_channel * std::sin(theta_true_deg * M_PI / 180.0) / lambda
           + cfg.ati_phase_bias_rad;
}

struct MotionCompResult {
    bool ok = false;
    bool used_motion_comp = false;
    bool fallback_legacy = false;
    double phi_meas = std::numeric_limits<double>::quiet_NaN();
    double phi_static = std::numeric_limits<double>::quiet_NaN();
    double phi_static_total = std::numeric_limits<double>::quiet_NaN();
    double phi_res = std::numeric_limits<double>::quiet_NaN();
    double phi_static_at_zero = std::numeric_limits<double>::quiet_NaN();
    double phi_res_at_zero = std::numeric_limits<double>::quiet_NaN();
    double phi_static_geometry = std::numeric_limits<double>::quiet_NaN();
    double af_phase = std::numeric_limits<double>::quiet_NaN();
    double af_total = std::numeric_limits<double>::quiet_NaN();
    double af_geometry = std::numeric_limits<double>::quiet_NaN();
    double af_motion = 0.0;
    double v_radial = std::numeric_limits<double>::quiet_NaN();
    double v_from_phase_raw = std::numeric_limits<double>::quiet_NaN();
    double v_from_phi_res = std::numeric_limits<double>::quiet_NaN();
    double phi_motion = std::numeric_limits<double>::quiet_NaN();
    double delta_t = std::numeric_limits<double>::quiet_NaN();
    double denom = std::numeric_limits<double>::quiet_NaN();
    double denom_without_k = std::numeric_limits<double>::quiet_NaN();
    double k = std::numeric_limits<double>::quiet_NaN();
    double b = std::numeric_limits<double>::quiet_NaN();
    int p38_theory_sign = 1;
    int motion_doppler_axis_sign = 1;
    int ati_phase_to_velocity_sign = 1;
    const char *p38_mode = "clutter_cancel_38_paper_1_p38_cuda";
    const char *geometry_calib_mode = "strict_dual_channel_forward";
    const char *solver = "old";
};

inline MotionCompResult compensateTargetMotionDoppler(const Config &cfg,
                                                      const GMTIOutput::Plane &plane,
                                                      double phase,
                                                      double af_total,
                                                      double k,
                                                      double b,
                                                      double lambda,
                                                      double theta_true_deg)
{
    MotionCompResult out;
    out.k = k;
    out.b = b;
    out.motion_doppler_axis_sign = (cfg.motion_doppler_axis_sign < 0) ? -1 : 1;
    out.ati_phase_to_velocity_sign = (cfg.ati_phase_to_velocity_sign < 0) ? -1 : 1;
    out.phi_meas = phase;
    if (std::abs(k) >= 1.0e-12) {
        out.af_phase = (phase - b) / k;
    } else {
        out.af_phase = 0.0;
    }
    out.af_total = std::isfinite(af_total) ? af_total : out.af_phase;
    out.af_geometry = out.af_phase;
    out.phi_static_at_zero = cfg.ati_phase_bias_rad;
    out.phi_res_at_zero = wrapPiLocal(phase - out.phi_static_at_zero);
    out.phi_static = k * out.af_phase + b + cfg.ati_phase_bias_rad;
    out.phi_static_geometry = out.phi_static;

    if (!cfg.motion_comp_enable || !cfg.motion_comp_analytic_enable) {
        return out;
    }

    if (!(lambda > 0.0) && cfg.fc > 0.0) {
        lambda = C / (cfg.fc * 1.0e9);
    }
    const double delta_t = (plane.V > 0.0) ? (cfg.d_channel / plane.V)
                                           : std::numeric_limits<double>::quiet_NaN();
    const bool valid = std::isfinite(phase) &&
                       std::isfinite(out.af_total) &&
                       std::isfinite(out.af_phase) &&
                       std::isfinite(k) &&
                       std::isfinite(b) &&
                       std::abs(k) >= 1.0e-12 &&
                       lambda > 0.0 &&
                       plane.V > 1.0e-6 &&
                       cfg.d_channel > 1.0e-6 &&
                       delta_t > 1.0e-9 &&
                       std::isfinite(delta_t);
    if (!valid) {
        out.fallback_legacy = true;
        return out;
    }

    out.used_motion_comp = true;
    out.solver = "analytic";
    out.delta_t = delta_t;
    const double s_phi = (cfg.ati_phase_to_velocity_sign < 0) ? -1.0 : 1.0;
    const double s_dop = (cfg.motion_doppler_axis_sign < 0) ? -1.0 : 1.0;

    out.phi_static_total = staticPhaseTheoryRad(cfg, theta_true_deg, lambda);
    out.phi_res = wrapPiLocal(phase - out.phi_static_total);
    out.denom_without_k = s_phi * 4.0 * M_PI * delta_t / lambda;
    out.denom = out.denom_without_k;
    if (std::isfinite(out.denom_without_k) &&
        std::abs(out.denom_without_k) >= std::max(0.0, cfg.motion_comp_denom_min)) {
        out.v_from_phase_raw = phase / out.denom_without_k;
    }

    if (!std::isfinite(out.denom) ||
        std::abs(out.denom) < std::max(0.0, cfg.motion_comp_denom_min)) {
        out.fallback_legacy = true;
        return out;
    }

    out.v_radial = out.phi_res / out.denom;
    out.v_from_phi_res = out.v_radial;
    if (cfg.ati_vmax_mps > 0.0 && std::isfinite(out.v_radial)) {
        out.v_radial = std::max(-cfg.ati_vmax_mps,
                                std::min(cfg.ati_vmax_mps, out.v_radial));
    }
    out.af_motion = s_dop * 2.0 * out.v_radial / lambda;
    out.af_geometry = out.af_total - out.af_motion;
    out.phi_static_geometry = staticPhaseTheoryRad(cfg, theta_true_deg, lambda);
    out.phi_static = out.phi_static_geometry;
    out.phi_motion = wrapPiLocal(phase - out.phi_static_geometry);

    out.ok = std::isfinite(out.af_geometry) &&
             std::isfinite(out.v_radial) &&
             std::isfinite(out.af_motion) &&
             std::isfinite(out.phi_motion);
    out.fallback_legacy = !out.ok;
    return out;
}

#endif // MOTION_COMP_HPP
