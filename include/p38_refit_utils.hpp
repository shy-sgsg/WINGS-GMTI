#ifndef P38_REFIT_UTILS_HPP
#define P38_REFIT_UTILS_HPP

#include "config_structs.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <complex>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

struct P38StageMetrics {
    std::array<double, 2> p38{{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    }};
    double rmse = std::numeric_limits<double>::quiet_NaN();
    int sample_count = 0;
    double inlier_ratio = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
    std::string source;
};

inline std::vector<uint8_t> buildP38ExcludeMask(const Config &cfg,
                                                int Na,
                                                int Nr,
                                                const std::vector<int> &cluster_rows,
                                                const std::vector<int> &cluster_cols,
                                                const std::vector<int> &det_rows,
                                                const std::vector<int> &det_cols,
                                                const std::vector<float> &power_map)
{
    const size_t total = static_cast<size_t>(std::max(0, Na)) * static_cast<size_t>(std::max(0, Nr));
    std::vector<uint8_t> mask(total, 0u);
    if (Na <= 0 || Nr <= 0 || power_map.size() != total) {
        return mask;
    }

    const int row_guard = std::max(0, cfg.p38_refit_row_guard_bins);
    const int range_guard = std::max(0, cfg.p38_refit_range_guard_bins);

    auto mark_rect = [&](int r0, int r1, int c0, int c1) {
        r0 = std::max(0, r0);
        c0 = std::max(0, c0);
        r1 = std::min(Na - 1, r1);
        c1 = std::min(Nr - 1, c1);
        if (r0 > r1 || c0 > c1) return;
        for (int r = r0; r <= r1; ++r) {
            const size_t off = static_cast<size_t>(r) * static_cast<size_t>(Nr);
            for (int c = c0; c <= c1; ++c) {
                mask[off + static_cast<size_t>(c)] = 1u;
            }
        }
    };

    for (size_t i = 0; i < cluster_rows.size() && i < cluster_cols.size(); ++i) {
        const int r = cluster_rows[i];
        const int c = cluster_cols[i];
        mark_rect(r - row_guard, r + row_guard, c - range_guard, c + range_guard);
    }
    for (size_t i = 0; i < det_rows.size() && i < det_cols.size(); ++i) {
        const int r = det_rows[i];
        const int c = det_cols[i];
        mark_rect(r - row_guard, r + row_guard, c - range_guard, c + range_guard);
    }

    for (int r = 0; r < Na; ++r) {
        for (int c = 0; c < Nr; ++c) {
            const size_t idx = static_cast<size_t>(r) * static_cast<size_t>(Nr) + static_cast<size_t>(c);
            const float p = power_map[idx];
            if (!std::isfinite(p)) {
                mask[idx] = 1u;
                continue;
            }
            if (r < row_guard || r >= Na - row_guard ||
                c < range_guard || c >= Nr - range_guard) {
                mask[idx] = 1u;
            }
        }
    }

    std::vector<float> finite_power;
    finite_power.reserve(total);
    for (float p : power_map) {
        if (std::isfinite(p)) finite_power.push_back(p);
    }
    if (!finite_power.empty() && cfg.p38_refit_top_power_frac > 0.0) {
        const double keep_frac = std::max(0.0, std::min(1.0, 1.0 - cfg.p38_refit_top_power_frac));
        const size_t nth = static_cast<size_t>(std::floor(keep_frac * static_cast<double>(finite_power.size() - 1)));
        std::nth_element(finite_power.begin(), finite_power.begin() + static_cast<std::ptrdiff_t>(nth), finite_power.end());
        const float threshold = finite_power[nth];
        for (size_t i = 0; i < power_map.size(); ++i) {
            if (std::isfinite(power_map[i]) && power_map[i] >= threshold) {
                mask[i] = 1u;
            }
        }
    }

    return mask;
}

inline P38StageMetrics evaluateP38FitMetrics(const std::vector<double> &phase_samples,
                                             const std::vector<double> &row_fa_samples,
                                             const std::array<double, 2> &p38)
{
    P38StageMetrics out;
    out.p38 = p38;
    if (phase_samples.size() != row_fa_samples.size() || phase_samples.empty()) {
        return out;
    }

    const double k = p38[0];
    const double b = p38[1];
    double sum_sq = 0.0;
    int inliers = 0;
    int valid = 0;
    const double inlier_thresh = 0.35;
    for (size_t i = 0; i < phase_samples.size(); ++i) {
        const double phi = phase_samples[i];
        const double x = row_fa_samples[i];
        if (!std::isfinite(phi) || !std::isfinite(x)) continue;
        ++valid;
        double residual = phi - (k * x + b);
        residual = std::fmod(residual + M_PI, 2.0 * M_PI);
        if (residual < 0.0) residual += 2.0 * M_PI;
        residual -= M_PI;
        sum_sq += residual * residual;
        if (std::abs(residual) <= inlier_thresh) ++inliers;
    }
    out.sample_count = valid;
    if (valid > 0) {
        out.rmse = std::sqrt(sum_sq / static_cast<double>(valid));
        out.inlier_ratio = static_cast<double>(inliers) / static_cast<double>(valid);
        out.valid = true;
    }
    return out;
}

template <typename ComplexT, typename P38T, typename EstimateFn>
inline P38StageMetrics refitP38WithMask(const Config &cfg,
                                        const std::vector<double> &fa_axis,
                                        const std::vector<ComplexT> &F1_pre,
                                        const std::vector<ComplexT> &F2_pre,
                                        int az_st,
                                        int rg_st,
                                        int az_ed,
                                        int rg_ed,
                                        const std::vector<int> &cluster_rows,
                                        const std::vector<int> &cluster_cols,
                                        const std::vector<int> &det_rows,
                                        const std::vector<int> &det_cols,
                                        const std::vector<float> &phase_map,
                                        const std::vector<float> &power_map,
                                        EstimateFn &&estimate_fn,
                                        std::array<P38T, 2> &p38_out,
                                        std::vector<double> *phase_samples_out = nullptr,
                                        std::vector<double> *row_samples_out = nullptr,
                                        std::vector<uint8_t> *mask_out = nullptr)
{
    P38StageMetrics metrics;
    if (!cfg.p38_refit_enable) {
        return metrics;
    }
    const int Na = static_cast<int>(fa_axis.size());
    const int Nr = Na > 0 ? static_cast<int>(F1_pre.size() / static_cast<size_t>(Na)) : 0;
    if (Na <= 0 || Nr <= 0 || F1_pre.size() != F2_pre.size() ||
        static_cast<size_t>(Na) * static_cast<size_t>(Nr) != F1_pre.size() ||
        phase_map.size() != F1_pre.size() || power_map.size() != F1_pre.size()) {
        return metrics;
    }

    std::vector<uint8_t> mask = buildP38ExcludeMask(cfg, Na, Nr, cluster_rows, cluster_cols,
                                                    det_rows, det_cols, power_map);
    if (mask_out) {
        *mask_out = mask;
    }

    std::vector<ComplexT> F1_masked = F1_pre;
    std::vector<ComplexT> F2_masked = F2_pre;
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
            F1_masked[i] = ComplexT(0.0, 0.0);
            F2_masked[i] = ComplexT(0.0, 0.0);
        }
    }

    std::vector<ComplexT> prosig_38;
    std::vector<double> phase_samples;
    std::vector<double> row_samples;
    if (!estimate_fn(fa_axis,
                     F1_masked,
                     F2_masked,
                     az_st, rg_st, az_ed, rg_ed,
                     cfg,
                     prosig_38,
                     p38_out,
                     phase_samples,
                     row_samples)) {
        return metrics;
    }

    if (phase_samples_out) {
        *phase_samples_out = phase_samples;
    }
    if (row_samples_out) {
        *row_samples_out = row_samples;
    }
    std::array<double, 2> p38_double{{static_cast<double>(p38_out[0]), static_cast<double>(p38_out[1])}};
    metrics.p38 = p38_double;
    metrics = evaluateP38FitMetrics(phase_samples, row_samples, p38_double);
    return metrics;
}

#endif // P38_REFIT_UTILS_HPP
