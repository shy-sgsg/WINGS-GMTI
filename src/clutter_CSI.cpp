#include "GMTIProcessor.hpp"
#include "trig_lut.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdint>
#include <complex>
#include <cmath>
#include <algorithm> // lower_bound
#include <numeric>   // accumulate

bool GMTIProcessor::clutter_cancel_38_paper_1_p38(
    const std::vector<float>& y_faAxis,
    const std::vector<std::complex<float>>& F1f,
    const std::vector<std::complex<float>>& F2f,
    int az_st, int rg_st, int az_ed, int rg_ed,
    const Config& cfg,
    std::array<float,2>& p_38
) {
    const size_t Na = static_cast<size_t>(effectivePulseNum(cfg));
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    const size_t total = Na * Nr;

    if (F1f.size() != total || F2f.size() != total) return false;
    if (y_faAxis.size() != Na) return false;

    az_st = std::max(0, std::min(az_st, int(Na) - 1));
    az_ed = std::max(0, std::min(az_ed, int(Na) - 1));
    if (az_st > az_ed) return false;

    std::vector<float> phase_tra_38(Na, 0.0f), row_fa(Na, 0.0f);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int rr = az_st; rr <= az_ed; ++rr) {
        const size_t off = static_cast<size_t>(rr) * Nr;
        float num_re = 0.0f, num_im = 0.0f, den = 0.0f;
        for (size_t c = 0; c < Nr; ++c) {
            const auto &a = F1f[off + c];
            const auto &b = F2f[off + c];
            num_re += (a.real() * b.real() + a.imag() * b.imag());
            num_im += (a.imag() * b.real() - a.real() * b.imag());
            den += std::norm(a);
        }
        const std::complex<float> num(num_re, num_im);
        phase_tra_38[static_cast<size_t>(rr)] = (den > 0.0f) ? std::arg(num / den) : std::arg(num);
        row_fa[static_cast<size_t>(rr)] = y_faAxis[static_cast<size_t>(rr)];
    }

    auto unwrap_inplace = [](std::vector<float> &v) {
        for (size_t i = 1; i < v.size(); ++i) {
            float d = v[i] - v[i - 1];
            if (d > static_cast<float>(M_PI)) v[i] -= static_cast<float>(2.0 * M_PI);
            if (d < static_cast<float>(-M_PI)) v[i] += static_cast<float>(2.0 * M_PI);
        }
    };

    int mid_num = std::max(1, std::min(int(std::round(double(az_ed - az_st + 1))), int(Na)));
    const float mid_phase = phase_tra_38[size_t(mid_num - 1)];
    unwrap_inplace(phase_tra_38);
    const float delta = mid_phase - phase_tra_38[size_t(mid_num - 1)];
    for (float &v : phase_tra_38) v += delta;

    const int M = az_ed - az_st + 1;
    std::vector<float> phase_tra_38_cut(M);
    std::vector<float> row_fa_cut(M);
    for (int i = 0; i < M; ++i) {
        phase_tra_38_cut[i] = phase_tra_38[size_t(az_st + i)];
        row_fa_cut[i] = row_fa[size_t(az_st + i)];
    }

    auto linfit = [](const std::vector<float> &x, const std::vector<float> &y, float &k, float &b) -> bool {
        if (x.size() != y.size() || x.size() < 2) return false;
        const size_t n = x.size();
        float Sx = 0.0f, Sy = 0.0f, Sxx = 0.0f, Sxy = 0.0f;
        for (size_t i = 0; i < n; ++i) {
            Sx += x[i];
            Sy += y[i];
            Sxx += x[i] * x[i];
            Sxy += x[i] * y[i];
        }
        float det = static_cast<float>(n) * Sxx - Sx * Sx;
        if (std::fabs(det) < 1e-8f) return false;
        k = (static_cast<float>(n) * Sxy - Sx * Sy) / det;
        b = (Sy * Sxx - Sx * Sxy) / det;
        return true;
    };

    float k = 0.0f, b = 0.0f;
    if (!linfit(row_fa_cut, phase_tra_38_cut, k, b)) return false;
    p_38 = {k, b};
    return true;
}

bool GMTIProcessor::clutter_cancel_38_paper_1(
    const std::vector<float>& y_faAxis,
    const std::vector<std::complex<float>>& F1f,
    const std::vector<std::complex<float>>& F2f,
    int az_st, int rg_st, int az_ed, int rg_ed,
    const Config& cfg,
    std::vector<std::complex<float>>& prosig_38,
    std::array<float,2>& p_38,
    std::vector<float>& phase_tra_38_cut,
    std::vector<float>& row_fa_cut
) {
    auto wall_start = std::chrono::high_resolution_clock::now();

    const size_t Na = static_cast<size_t>(effectivePulseNum(cfg));
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    const size_t total = Na * Nr;

    if (F1f.size() != total || F2f.size() != total) return false;
    if (y_faAxis.size() != Na) return false;

    az_st = std::max(0, std::min(az_st, int(Na) - 1));
    az_ed = std::max(0, std::min(az_ed, int(Na) - 1));
    if (az_st > az_ed) return false;

    if (!clutter_cancel_38_paper_1_p38(y_faAxis, F1f, F2f, az_st, rg_st, az_ed, rg_ed, cfg, p_38)) {
        return false;
    }

    std::vector<float> phase_tra_38(Na, 0.0f), row_fa(Na, 0.0f);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int rr = az_st; rr <= az_ed; ++rr) {
        const size_t off = static_cast<size_t>(rr) * Nr;
        float num_re = 0.0f, num_im = 0.0f, den = 0.0f;
        for (size_t c = 0; c < Nr; ++c) {
            const auto &a = F1f[off + c];
            const auto &b = F2f[off + c];
            num_re += (a.real() * b.real() + a.imag() * b.imag());
            num_im += (a.imag() * b.real() - a.real() * b.imag());
            den += std::norm(a);
        }
        const std::complex<float> num(num_re, num_im);
        phase_tra_38[static_cast<size_t>(rr)] = (den > 0.0f) ? std::arg(num / den) : std::arg(num);
        row_fa[static_cast<size_t>(rr)] = y_faAxis[static_cast<size_t>(rr)];
    }

    auto unwrap_inplace = [](std::vector<float> &v) {
        for (size_t i = 1; i < v.size(); ++i) {
            float d = v[i] - v[i - 1];
            if (d > static_cast<float>(M_PI)) v[i] -= static_cast<float>(2.0 * M_PI);
            if (d < static_cast<float>(-M_PI)) v[i] += static_cast<float>(2.0 * M_PI);
        }
    };
    int mid_num = std::max(1, std::min(int(std::round(double(az_ed - az_st + 1))), int(Na)));
    const float mid_phase = phase_tra_38[size_t(mid_num - 1)];
    unwrap_inplace(phase_tra_38);
    const float delta = mid_phase - phase_tra_38[size_t(mid_num - 1)];
    for (float &v : phase_tra_38) v += delta;

    const int M = az_ed - az_st + 1;
    phase_tra_38_cut.resize(M);
    row_fa_cut.resize(M);
    for (int i = 0; i < M; ++i) {
        phase_tra_38_cut[i] = phase_tra_38[size_t(az_st + i)];
        row_fa_cut[i] = row_fa[size_t(az_st + i)];
    }

    prosig_38.assign(total, std::complex<float>(0.0f, 0.0f));
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t r = 0; r < Na; ++r) {
        const float phi = p_38[0] * y_faAxis[r] + p_38[1];
        const float cs = static_cast<float>(std::cos(phi));
        const float sn = static_cast<float>(std::sin(phi));
        const std::complex<float> az_fai(cs, sn);
        const size_t off = r * Nr;
        for (size_t c = 0; c < Nr; ++c) {
            const auto &a = F1f[off + c];
            const auto b2 = az_fai * F2f[off + c];
            const float aa = std::abs(a), bb = std::abs(b2);
            const float m = (aa < bb) ? aa : bb;
            const auto a_eq = (aa > 0.0f) ? (m / aa) * a : std::complex<float>(0.0f, 0.0f);
            const auto b_eq = (bb > 0.0f) ? (m / bb) * b2 : std::complex<float>(0.0f, 0.0f);
            prosig_38[off + c] = a_eq - b_eq;
        }
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> wall_ms = wall_end - wall_start;
    (void)wall_ms;
    return true;
}
