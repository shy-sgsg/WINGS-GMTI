#ifndef MOTION_COMP_HPP
#define MOTION_COMP_HPP

#include "config_structs.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <limits>
#include <string>

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

inline double unwrapNearLocal(double phase, double ref)
{
    return ref + wrapPiLocal(phase - ref);
}

inline std::string normalizeMotionCompSolver(std::string solver)
{
    std::string out;
    out.reserve(solver.size());
    for (unsigned char ch : solver) {
        if (!std::isspace(ch)) {
            out.push_back(static_cast<char>(std::tolower(ch)));
        }
    }
    if (out.empty()) {
        return "analytic";
    }
    if (out == "iter") {
        return "iterative";
    }
    if (out == "root") {
        return "root1d";
    }
    if (out == "dbg") {
        return "debug";
    }
    return out;
}

inline double staticPhaseModel(double af_geo,
                               double k,
                               double b,
                               double ati_phase_bias_rad)
{
    return k * af_geo + b + ati_phase_bias_rad;
}

inline double strictGeometryStaticPhaseModel(double af_geo,
                                             double delta_t,
                                             double ati_phase_bias_rad)
{
    return 2.0 * M_PI * delta_t * af_geo + ati_phase_bias_rad;
}

inline double motionStaticPhaseModel(double af_geo,
                                     double k,
                                     double b,
                                     double ati_phase_bias_rad,
                                     double delta_t,
                                     double static_phase_ref = std::numeric_limits<double>::quiet_NaN())
{
    if (std::isfinite(static_phase_ref)) {
        return static_phase_ref + ati_phase_bias_rad;
    }
    if (std::isfinite(delta_t) && delta_t > 0.0) {
        return strictGeometryStaticPhaseModel(af_geo, delta_t, ati_phase_bias_rad);
    }
    return staticPhaseModel(af_geo, k, b, ati_phase_bias_rad);
}

inline double motionStaticPhaseKeff(double delta_t,
                                    double fallback_k,
                                    double static_phase_ref = std::numeric_limits<double>::quiet_NaN())
{
    if (std::isfinite(static_phase_ref)) {
        return 0.0;
    }
    if (std::isfinite(delta_t) && delta_t > 0.0) {
        return 2.0 * M_PI * delta_t;
    }
    return fallback_k;
}

inline double beamStaticPhaseRef(const Config &cfg,
                                 const GMTIOutput::Plane &plane,
                                 double lambda,
                                 double theta_true_deg)
{
    if (!(lambda > 0.0) || !(plane.V > 0.0) || !(cfg.d_channel > 0.0) ||
        !std::isfinite(theta_true_deg)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double delta_t = cfg.d_channel / plane.V;
    const double theta = theta_true_deg * M_PI / 180.0;
    const double af_stationary = -2.0 * plane.V * std::sin(theta) / lambda;
    return -2.0 * M_PI * delta_t * af_stationary;
}

inline double calcStaticPhaseKeff(double af_total,
                                  double k,
                                  double b,
                                  double ati_phase_bias_rad,
                                  double fd_res_hint)
{
    const double df = std::max(1.0, 0.25 * std::max(1.0, fd_res_hint));
    const double phi0 = staticPhaseModel(af_total, k, b, ati_phase_bias_rad);
    double phip = staticPhaseModel(af_total + df, k, b, ati_phase_bias_rad);
    double phim = staticPhaseModel(af_total - df, k, b, ati_phase_bias_rad);
    phip = unwrapNearLocal(phip, phi0);
    phim = unwrapNearLocal(phim, phi0);
    return (phip - phim) / (2.0 * df);
}

struct MotionCompResult {
    bool ok = false;
    bool used_motion_comp = false;
    bool fallback_legacy = false;
    double phi_meas = std::numeric_limits<double>::quiet_NaN();
    double phi_static_model = std::numeric_limits<double>::quiet_NaN();
    double phi_static = std::numeric_limits<double>::quiet_NaN();
    double phi_static_total = std::numeric_limits<double>::quiet_NaN();
    double phi_res = std::numeric_limits<double>::quiet_NaN();
    double phi_static_at_zero = std::numeric_limits<double>::quiet_NaN();
    double phi_res_at_zero = std::numeric_limits<double>::quiet_NaN();
    double phi_static_geometry = std::numeric_limits<double>::quiet_NaN();
    double af_phase = std::numeric_limits<double>::quiet_NaN();
    double af_total = std::numeric_limits<double>::quiet_NaN();
    double af_geometry = std::numeric_limits<double>::quiet_NaN();
    double af_motion = std::numeric_limits<double>::quiet_NaN();
    double phi_motion = std::numeric_limits<double>::quiet_NaN();
    double delta_t = std::numeric_limits<double>::quiet_NaN();
    double denom = std::numeric_limits<double>::quiet_NaN();
    double denom_without_k = std::numeric_limits<double>::quiet_NaN();
    double k = std::numeric_limits<double>::quiet_NaN();
    double b = std::numeric_limits<double>::quiet_NaN();
    double C_ati = std::numeric_limits<double>::quiet_NaN();
    double k_eff_static_phase_df = std::numeric_limits<double>::quiet_NaN();
    double v_radial = std::numeric_limits<double>::quiet_NaN();
    double v_from_phase_raw = std::numeric_limits<double>::quiet_NaN();
    double v_from_phi_res = std::numeric_limits<double>::quiet_NaN();
    double v_iterative_mps = std::numeric_limits<double>::quiet_NaN();
    double v_analytic_mps = std::numeric_limits<double>::quiet_NaN();
    double v_root1d_mps = std::numeric_limits<double>::quiet_NaN();
    double v_old_mps = std::numeric_limits<double>::quiet_NaN();
    double af_geometry_old_hz = std::numeric_limits<double>::quiet_NaN();
    double af_geometry_iterative_hz = std::numeric_limits<double>::quiet_NaN();
    double af_geometry_analytic_hz = std::numeric_limits<double>::quiet_NaN();
    double af_geometry_root1d_hz = std::numeric_limits<double>::quiet_NaN();
    double root1d_cost = std::numeric_limits<double>::quiet_NaN();
    int p38_theory_sign = 1;
    int motion_doppler_axis_sign = 1;
    int ati_phase_to_velocity_sign = 1;
    std::string phi_static_model_name = "p38_linear_fit";
    std::string p38_mode = "clutter_cancel_38_paper_1_p38_cuda";
    std::string geometry_calib_mode = "linear_p38_phase_vs_doppler";
    std::string solver_requested = "analytic";
    std::string solver = "old";
};

inline MotionCompResult solveOldMotionComp(const Config &cfg,
                                           double phase,
                                           double af_total,
                                           double k,
                                           double b,
                                           double lambda,
                                           double static_phase_ref = std::numeric_limits<double>::quiet_NaN())
{
    MotionCompResult out;
    (void)static_phase_ref;
    out.k = k;
    out.b = b;
    out.phi_meas = phase;
    out.af_phase = (std::abs(k) >= 1.0e-12) ? ((phase - b) / k) : 0.0;
    out.af_total = std::isfinite(af_total) ? af_total : out.af_phase;
    out.af_geometry = out.af_phase;
    out.phi_static_model = staticPhaseModel(out.af_total, k, b, cfg.ati_phase_bias_rad);
    out.phi_static = out.phi_static_model;
    out.phi_static_total = out.phi_static_model;
    out.phi_static_at_zero = staticPhaseModel(0.0, k, b, cfg.ati_phase_bias_rad);
    out.phi_res_at_zero = wrapPiLocal(phase - out.phi_static_at_zero);
    out.phi_res = wrapPiLocal(phase - out.phi_static_total);
    out.phi_static_geometry = out.phi_static_model;
    out.phi_motion = out.phi_res;
    out.phi_static_model_name = "p38_linear_fit";
    out.geometry_calib_mode = "linear_p38_phase_vs_doppler";
    out.solver_requested = normalizeMotionCompSolver(cfg.motion_comp_solver);
    out.solver = "old";
    out.v_radial = 0.0;
    out.af_motion = 0.0;
    out.v_old_mps = 0.0;
    out.af_geometry_old_hz = out.af_geometry;
    out.used_motion_comp = false;
    if (!(lambda > 0.0) && cfg.fc > 0.0) {
        lambda = C / (cfg.fc * 1.0e9);
    }
    const double delta_t = std::numeric_limits<double>::quiet_NaN();
    out.delta_t = delta_t;
    out.C_ati = std::numeric_limits<double>::quiet_NaN();
    out.k_eff_static_phase_df = std::numeric_limits<double>::quiet_NaN();
    return out;
}

inline MotionCompResult solveIterativeMotionComp(const Config &cfg,
                                                 const GMTIOutput::Plane &plane,
                                                 double phase,
                                                 double af_total,
                                                 double k,
                                                 double b,
                                                 double lambda,
                                                 double static_phase_ref = std::numeric_limits<double>::quiet_NaN())
{
    MotionCompResult out;
    out.k = k;
    out.b = b;
    out.phi_meas = phase;
    out.af_phase = (std::abs(k) >= 1.0e-12) ? ((phase - b) / k) : 0.0;
    out.af_total = std::isfinite(af_total) ? af_total : out.af_phase;
    out.phi_static_model_name = "strict_two_channel_geometry";
    out.geometry_calib_mode = "strict_two_channel_forward";
    out.solver_requested = "iterative";
    out.solver = "iterative";
    out.p38_theory_sign = 1;
    out.motion_doppler_axis_sign = (cfg.motion_doppler_axis_sign < 0) ? -1 : 1;
    out.ati_phase_to_velocity_sign = (cfg.ati_phase_to_velocity_sign < 0) ? -1 : 1;
    out.used_motion_comp = true;

    if (!(lambda > 0.0) && cfg.fc > 0.0) {
        lambda = C / (cfg.fc * 1.0e9);
    }
    if (!(lambda > 0.0) || !(cfg.d_channel > 0.0) || !(cfg.ati_vmax_mps > 0.0)) {
        out.fallback_legacy = true;
        return out;
    }

    const double delta_t = (plane.V > 0.0) ? (cfg.d_channel / plane.V)
                                           : std::numeric_limits<double>::quiet_NaN();
    out.delta_t = delta_t;
    const double s_phi = (cfg.ati_phase_to_velocity_sign < 0) ? -1.0 : 1.0;
    const double s_dop = (cfg.motion_doppler_axis_sign < 0) ? -1.0 : 1.0;
    const double C_ati = s_phi * 4.0 * M_PI * delta_t / lambda;
    out.C_ati = C_ati;
    out.denom_without_k = C_ati;
    out.k_eff_static_phase_df = motionStaticPhaseKeff(delta_t, k, static_phase_ref);
    if (!std::isfinite(C_ati) || std::abs(C_ati) < std::max(0.0, cfg.motion_comp_denom_min)) {
        out.fallback_legacy = true;
        return out;
    }

    double v = 0.0;
    double last_v = v;
    const int max_iter = std::max(1, cfg.motion_comp_iter);
    const double vmax = (cfg.ati_vmax_mps > 0.0) ? cfg.ati_vmax_mps : std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < max_iter; ++iter) {
        const double af_motion = s_dop * 2.0 * v / lambda;
        const double af_geo = out.af_total - af_motion;
        const double phi_static = motionStaticPhaseModel(af_geo, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
        const double phi_res = wrapPiLocal(phase - phi_static);
        const double v_new = phi_res / C_ati;
        if (!std::isfinite(v_new) || std::abs(v_new) > vmax) {
            out.fallback_legacy = true;
            return out;
        }
        last_v = v;
        v = v_new;
        if (std::abs(v - last_v) < std::max(0.0, cfg.motion_comp_iter_tol_mps)) {
            break;
        }
    }

    out.v_radial = v;
    out.v_iterative_mps = v;
    out.af_motion = s_dop * 2.0 * v / lambda;
    out.af_geometry = out.af_total - out.af_motion;
    out.af_geometry_iterative_hz = out.af_geometry;
    out.phi_static_total = motionStaticPhaseModel(out.af_total, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_static_model = out.phi_static_total;
    out.phi_static = motionStaticPhaseModel(out.af_geometry, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_static_geometry = out.phi_static;
    out.phi_static_at_zero = motionStaticPhaseModel(0.0, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_res_at_zero = wrapPiLocal(phase - out.phi_static_at_zero);
    out.phi_res = wrapPiLocal(phase - out.phi_static_total);
    out.phi_motion = wrapPiLocal(phase - out.phi_static_geometry);
    out.v_from_phase_raw = phase / C_ati;
    out.v_from_phi_res = out.phi_res / C_ati;
    out.denom = C_ati;
    out.ok = std::isfinite(out.v_radial) &&
             std::isfinite(out.af_geometry) &&
             std::isfinite(out.af_motion) &&
             std::isfinite(out.phi_motion);
    out.fallback_legacy = !out.ok;
    return out;
}

inline MotionCompResult solveAnalyticMotionComp(const Config &cfg,
                                                const GMTIOutput::Plane &plane,
                                                double phase,
                                                double af_total,
                                                double k,
                                                double b,
                                                double lambda,
                                                double static_phase_ref = std::numeric_limits<double>::quiet_NaN())
{
    MotionCompResult out;
    out.k = k;
    out.b = b;
    out.phi_meas = phase;
    out.af_phase = (std::abs(k) >= 1.0e-12) ? ((phase - b) / k) : 0.0;
    out.af_total = std::isfinite(af_total) ? af_total : out.af_phase;
    out.phi_static_model_name = "strict_two_channel_geometry";
    out.geometry_calib_mode = "strict_two_channel_forward";
    out.solver_requested = "analytic";
    out.solver = "analytic";
    out.p38_theory_sign = 1;
    out.motion_doppler_axis_sign = (cfg.motion_doppler_axis_sign < 0) ? -1 : 1;
    out.ati_phase_to_velocity_sign = (cfg.ati_phase_to_velocity_sign < 0) ? -1 : 1;
    out.used_motion_comp = true;

    if (!(lambda > 0.0) && cfg.fc > 0.0) {
        lambda = C / (cfg.fc * 1.0e9);
    }
    const double delta_t = (plane.V > 0.0) ? (cfg.d_channel / plane.V)
                                           : std::numeric_limits<double>::quiet_NaN();
    out.delta_t = delta_t;
    if (!std::isfinite(phase) ||
        !(lambda > 0.0) ||
        !(plane.V > 1.0e-6) ||
        !(cfg.d_channel > 1.0e-6) ||
        !std::isfinite(delta_t) ||
        !(std::abs(k) >= 1.0e-12)) {
        out.fallback_legacy = true;
        return out;
    }

    const double s_phi = (cfg.ati_phase_to_velocity_sign < 0) ? -1.0 : 1.0;
    const double s_dop = (cfg.motion_doppler_axis_sign < 0) ? -1.0 : 1.0;
    const double C_ati = s_phi * 4.0 * M_PI * delta_t / lambda;
    const double phi_static_total = motionStaticPhaseModel(out.af_total, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    const double phi_res = wrapPiLocal(phase - phi_static_total);
    const double k_eff = motionStaticPhaseKeff(delta_t, k, static_phase_ref);
    const double denom = C_ati - k_eff * s_dop * 2.0 / lambda;

    out.phi_static_model = phi_static_total;
    out.phi_static_total = phi_static_total;
    out.phi_static_at_zero = motionStaticPhaseModel(0.0, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_res_at_zero = wrapPiLocal(phase - out.phi_static_at_zero);
    out.phi_res = phi_res;
    out.C_ati = C_ati;
    out.k_eff_static_phase_df = k_eff;
    out.denom_without_k = C_ati;
    out.denom = denom;
    out.v_from_phase_raw = phase / C_ati;
    out.v_from_phi_res = phi_res / C_ati;

    if (!std::isfinite(denom) || std::abs(denom) < std::max(0.0, cfg.motion_comp_denom_min)) {
        out.fallback_legacy = true;
        return out;
    }

    auto analytic_residual = [&](double v) -> double {
        const double af_motion = s_dop * 2.0 * v / lambda;
        const double af_geo = out.af_total - af_motion;
        const double phi_static = motionStaticPhaseModel(af_geo, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
        const double phi_pred = phi_static + C_ati * v;
        return wrapPiLocal(phase - phi_pred);
    };
    auto analytic_cost = [&](double v) -> double {
        const double r = analytic_residual(v);
        return 1.0 - std::cos(r);
    };
    const double vmax = (cfg.ati_vmax_mps > 0.0) ? cfg.ati_vmax_mps : 1.6;
    const double step = std::max(0.001, cfg.motion_comp_root_grid_step_mps);
    double best_v = std::max(-vmax, std::min(vmax, phi_res / denom));
    double best_cost = analytic_cost(best_v);
    for (double v = -vmax; v <= vmax + 0.5 * step; v += step) {
        const double c = analytic_cost(v);
        if (c < best_cost) {
            best_cost = c;
            best_v = v;
        }
    }
    double lo = std::max(-vmax, best_v - step);
    double hi = std::min(vmax, best_v + step);
    const double gr = 0.6180339887498949;
    double c = hi - gr * (hi - lo);
    double d = lo + gr * (hi - lo);
    double fc = analytic_cost(c);
    double fd = analytic_cost(d);
    for (int i = 0; i < 48; ++i) {
        if (fc < fd) {
            hi = d;
            d = c;
            fd = fc;
            c = hi - gr * (hi - lo);
            fc = analytic_cost(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + gr * (hi - lo);
            fd = analytic_cost(d);
        }
    }
    out.v_radial = (fc < fd) ? c : d;
    out.af_motion = s_dop * 2.0 * out.v_radial / lambda;
    out.af_geometry = out.af_total - out.af_motion;
    out.phi_static = motionStaticPhaseModel(out.af_geometry, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_static_geometry = out.phi_static;
    out.phi_motion = wrapPiLocal(phase - motionStaticPhaseModel(out.af_geometry, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref));
    out.ok = std::isfinite(out.v_radial) &&
             std::isfinite(out.af_geometry) &&
             std::isfinite(out.af_motion) &&
             std::isfinite(out.phi_motion);
    out.fallback_legacy = !out.ok;
    out.v_analytic_mps = out.v_radial;
    out.af_geometry_analytic_hz = out.af_geometry;
    return out;
}

inline MotionCompResult solveRoot1DMotionComp(const Config &cfg,
                                              const GMTIOutput::Plane &plane,
                                              double phase,
                                              double af_total,
                                              double k,
                                              double b,
                                              double lambda,
                                              double static_phase_ref = std::numeric_limits<double>::quiet_NaN())
{
    MotionCompResult out;
    out.k = k;
    out.b = b;
    out.phi_meas = phase;
    out.af_phase = (std::abs(k) >= 1.0e-12) ? ((phase - b) / k) : 0.0;
    out.af_total = std::isfinite(af_total) ? af_total : out.af_phase;
    out.phi_static_model_name = "strict_two_channel_geometry";
    out.geometry_calib_mode = "strict_two_channel_forward";
    out.solver_requested = "root1d";
    out.solver = "root1d";
    out.p38_theory_sign = 1;
    out.motion_doppler_axis_sign = (cfg.motion_doppler_axis_sign < 0) ? -1 : 1;
    out.ati_phase_to_velocity_sign = (cfg.ati_phase_to_velocity_sign < 0) ? -1 : 1;
    out.used_motion_comp = true;

    if (!(lambda > 0.0) && cfg.fc > 0.0) {
        lambda = C / (cfg.fc * 1.0e9);
    }
    const double delta_t = (plane.V > 0.0) ? (cfg.d_channel / plane.V)
                                           : std::numeric_limits<double>::quiet_NaN();
    out.delta_t = delta_t;
    if (!std::isfinite(phase) ||
        !(lambda > 0.0) ||
        !(plane.V > 1.0e-6) ||
        !(cfg.d_channel > 1.0e-6) ||
        !std::isfinite(delta_t)) {
        out.fallback_legacy = true;
        return out;
    }

    const double s_phi = (cfg.ati_phase_to_velocity_sign < 0) ? -1.0 : 1.0;
    const double s_dop = (cfg.motion_doppler_axis_sign < 0) ? -1.0 : 1.0;
    const double C_ati = s_phi * 4.0 * M_PI * delta_t / lambda;
    out.C_ati = C_ati;
    out.denom_without_k = C_ati;
    out.k_eff_static_phase_df = motionStaticPhaseKeff(delta_t, k, static_phase_ref);

    const double vmax = (cfg.ati_vmax_mps > 0.0) ? cfg.ati_vmax_mps : 1.6;
    const double step = std::max(0.001, cfg.motion_comp_root_grid_step_mps);

    auto residual = [&](double v) -> double {
        const double af_motion = s_dop * 2.0 * v / lambda;
        const double af_geo = out.af_total - af_motion;
        const double phi_static = motionStaticPhaseModel(af_geo, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
        const double phi_pred = phi_static + C_ati * v;
        return wrapPiLocal(phase - phi_pred);
    };
    auto cost = [&](double v) -> double {
        const double r = residual(v);
        return 1.0 - std::cos(r);
    };

    double best_v = 0.0;
    double best_cost = std::numeric_limits<double>::infinity();
    for (double v = -vmax; v <= vmax + 0.5 * step; v += step) {
        const double c = cost(v);
        if (c < best_cost) {
            best_cost = c;
            best_v = v;
        }
    }

    double lo = std::max(-vmax, best_v - step);
    double hi = std::min(vmax, best_v + step);
    const double gr = 0.6180339887498949;
    double c = hi - gr * (hi - lo);
    double d = lo + gr * (hi - lo);
    double fc = cost(c);
    double fd = cost(d);
    for (int i = 0; i < 48; ++i) {
        if (fc < fd) {
            hi = d;
            d = c;
            fd = fc;
            c = hi - gr * (hi - lo);
            fc = cost(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + gr * (hi - lo);
            fd = cost(d);
        }
    }
    double v = (fc < fd) ? c : d;
    const double final_cost = std::min(fc, fd);
    out.root1d_cost = final_cost;
    out.v_radial = v;
    out.v_root1d_mps = v;

    if (!std::isfinite(v) || final_cost > std::max(0.0, cfg.motion_comp_root_cost_max)) {
        out.fallback_legacy = true;
        return out;
    }

    out.af_motion = s_dop * 2.0 * v / lambda;
    out.af_geometry = out.af_total - out.af_motion;
    out.af_geometry_root1d_hz = out.af_geometry;
    out.phi_static_model = motionStaticPhaseModel(out.af_total, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_static_total = motionStaticPhaseModel(out.af_total, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_static_at_zero = motionStaticPhaseModel(0.0, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_res_at_zero = wrapPiLocal(phase - out.phi_static_at_zero);
    out.phi_res = wrapPiLocal(phase - out.phi_static_total);
    out.phi_static = motionStaticPhaseModel(out.af_geometry, k, b, cfg.ati_phase_bias_rad, delta_t, static_phase_ref);
    out.phi_static_geometry = out.phi_static;
    out.phi_motion = wrapPiLocal(phase - out.phi_static_geometry);
    out.v_from_phase_raw = phase / C_ati;
    out.v_from_phi_res = out.phi_res / C_ati;
    out.ok = std::isfinite(out.v_radial) &&
             std::isfinite(out.af_geometry) &&
             std::isfinite(out.af_motion) &&
             std::isfinite(out.phi_motion);
    out.fallback_legacy = !out.ok;
    return out;
}

inline MotionCompResult solveMotionCompensation(const Config &cfg,
                                                const GMTIOutput::Plane &plane,
                                                double phase,
                                                double af_total,
                                                double k,
                                                double b,
                                                double lambda,
                                                double theta_true_deg = std::numeric_limits<double>::quiet_NaN())
{
    if (!(lambda > 0.0) && cfg.fc > 0.0) {
        lambda = C / (cfg.fc * 1.0e9);
    }
    const double static_phase_ref = beamStaticPhaseRef(cfg, plane, lambda, theta_true_deg);
    const std::string solver = normalizeMotionCompSolver(cfg.motion_comp_solver);
    MotionCompResult out;
    if (!cfg.motion_comp_enable || solver == "old") {
        out = solveOldMotionComp(cfg, phase, af_total, k, b, lambda, static_phase_ref);
    } else if (solver == "iterative") {
        out = solveIterativeMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);
    } else if (solver == "root1d") {
        out = solveRoot1DMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);
    } else if (solver == "debug") {
        const MotionCompResult old_res = solveOldMotionComp(cfg, phase, af_total, k, b, lambda, static_phase_ref);
        const MotionCompResult iter_res = solveIterativeMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);
        const MotionCompResult ana_res = solveAnalyticMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);
        const MotionCompResult root_res = solveRoot1DMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);

        out = ana_res;
        out.solver_requested = "debug";
        out.solver = "debug";
        out.used_motion_comp = true;
        out.v_old_mps = old_res.v_radial;
        out.af_geometry_old_hz = old_res.af_geometry;
        out.v_iterative_mps = iter_res.v_radial;
        out.v_analytic_mps = ana_res.v_radial;
        out.v_root1d_mps = root_res.v_radial;
        out.af_geometry_iterative_hz = iter_res.af_geometry;
        out.af_geometry_analytic_hz = ana_res.af_geometry;
        out.af_geometry_root1d_hz = root_res.af_geometry;
        out.root1d_cost = root_res.root1d_cost;
        return out;
    } else {
        out = solveAnalyticMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);
    }

    const bool want_compare = (solver == "debug");
    if (want_compare && cfg.motion_comp_enable) {
        const MotionCompResult old_res = solveOldMotionComp(cfg, phase, af_total, k, b, lambda, static_phase_ref);
        const MotionCompResult iter = solveIterativeMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);
        const MotionCompResult ana = solveAnalyticMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);
        const MotionCompResult root = solveRoot1DMotionComp(cfg, plane, phase, af_total, k, b, lambda, static_phase_ref);
        out.v_old_mps = old_res.v_radial;
        out.af_geometry_old_hz = old_res.af_geometry;
        out.v_iterative_mps = iter.v_radial;
        out.v_analytic_mps = ana.v_radial;
        out.v_root1d_mps = root.v_radial;
        out.af_geometry_iterative_hz = iter.af_geometry;
        out.af_geometry_analytic_hz = ana.af_geometry;
        out.af_geometry_root1d_hz = root.af_geometry;
        out.root1d_cost = root.root1d_cost;
    }
    return out;
}

#endif // MOTION_COMP_HPP
