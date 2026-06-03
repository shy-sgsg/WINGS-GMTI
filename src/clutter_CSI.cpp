#include "GMTIProcessor.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdint>
#include <complex>
#include <cmath>
#include <algorithm> // lower_bound
#include <numeric>   // accumulate

bool GMTIProcessor::clutter_38_V2(
    const std::vector<std::complex<double>>& data1,
    const std::vector<std::complex<double>>& data2,
    const std::vector<std::complex<double>>& rg_phi,   // len = cfg.rg_len
    int skip_num,
    int az_st, int rg_st,
    int az_ed, int rg_ed,
    const std::vector<double>& y_faAxis,               // len = cfg.pulse_num
    double fa2,
    const Config& cfg,                                  // ☆ 新增：显式传入 cfg
    std::vector<std::complex<double>>& CSI_result_38paper,
    std::array<double,2>& p_38,
    std::vector<double>& phase_38_phase,
    std::vector<double>& phase_38_fa
) {
    const size_t Na = static_cast<size_t>(cfg.pulse_num);
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    const size_t total = Na * Nr;

    if (data1.size() != total || data2.size() != total) return false;
    if (rg_phi.size() != Nr) return false;
    if (y_faAxis.size() != Na) return false;

    // 1) 对齐 + FFT(沿方位) + fftshift + DBS（此函数本身已从 cfg 取 Na/Nr）
    std::vector<std::complex<float>> data1_f(total);
    std::vector<std::complex<float>> data2_f(total);
    for (size_t i = 0; i < total; ++i) {
        data1_f[i] = std::complex<float>(static_cast<float>(data1[i].real()), static_cast<float>(data1[i].imag()));
        data2_f[i] = std::complex<float>(static_cast<float>(data2[i].real()), static_cast<float>(data2[i].imag()));
    }

    std::vector<std::complex<float>> F1f_f, F2f_f;
    if (!alignFFTAndDBS_CUDA(data1_f, data2_f, skip_num, static_cast<float>(fa2), cfg, F1f_f, F2f_f)) return false;

    std::vector<std::complex<double>> F1f(total), F2f(total);
    for (size_t i = 0; i < total; ++i) {
        const auto v1 = F1f_f[i];
        const auto v2 = F2f_f[i];
        F1f[i] = std::complex<double>(v1.real(), v1.imag());
        F2f[i] = std::complex<double>(v2.real(), v2.imag());
    }

    // 2) 列广播：F2f = F2f .* rg_phi
    for (size_t r = 0; r < Na; ++r) {
        const size_t off = r * Nr;
        for (size_t c = 0; c < Nr; ++c) F2f[off + c] *= rg_phi[c];
    }

    // 3) 38所方案对消
    return clutter_cancel_38_paper_1(
        y_faAxis, F1f, F2f,
        az_st, rg_st, az_ed, rg_ed,
        cfg,                                         // ☆ 传入 cfg
        CSI_result_38paper, p_38, phase_38_phase, phase_38_fa
    );
}

bool GMTIProcessor::clutter_cancel_38_paper_1_p38(
    const std::vector<double>& y_faAxis,
    const std::vector<std::complex<double>>& F1f,
    const std::vector<std::complex<double>>& F2f,
    int az_st, int rg_st, int az_ed, int rg_ed,
    const Config& cfg,
    std::array<double,2>& p_38
) {
    const size_t Na = static_cast<size_t>(cfg.pulse_num);
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    const size_t total = Na * Nr;

    if (F1f.size() != total || F2f.size() != total) return false;
    if (y_faAxis.size() != Na) return false;

    az_st = std::max(0, std::min(az_st, int(Na) - 1));
    az_ed = std::max(0, std::min(az_ed, int(Na) - 1));
    if (az_st > az_ed) return false;

    std::vector<double> phase_tra_38(Na, 0.0), row_fa(Na, 0.0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int rr = az_st; rr <= az_ed; ++rr) {
        const size_t off = static_cast<size_t>(rr) * Nr;

        double num_re = 0.0, num_im = 0.0, den = 0.0;
        for (size_t c = 0; c < Nr; ++c) {
            const auto& a = F1f[off + c];
            const auto& b = F2f[off + c];
            num_re += (a.real() * b.real() + a.imag() * b.imag());
            num_im += (a.imag() * b.real() - a.real() * b.imag());
            den += std::norm(a);
        }

        const std::complex<double> num(num_re, num_im);
        phase_tra_38[static_cast<size_t>(rr)] = (den > 0.0) ? std::arg(num / den) : std::arg(num);
        row_fa[static_cast<size_t>(rr)] = y_faAxis[static_cast<size_t>(rr)];
    }

    auto unwrap_inplace = [](std::vector<double>& v) {
        for (size_t i = 1; i < v.size(); ++i) {
            double d = v[i] - v[i - 1];
            if (d > M_PI) v[i] -= 2 * M_PI;
            if (d < -M_PI) v[i] += 2 * M_PI;
        }
    };

    int mid_num = std::max(1, std::min(int(std::round(double(az_ed - az_st + 1))), int(Na)));
    const double mid_phase = phase_tra_38[size_t(mid_num - 1)];
    unwrap_inplace(phase_tra_38);
    const double delta = mid_phase - phase_tra_38[size_t(mid_num - 1)];
    for (double& v : phase_tra_38) v += delta;

    const int M = az_ed - az_st + 1;
    std::vector<double> phase_tra_38_cut(M);
    std::vector<double> row_fa_cut(M);
    for (int i = 0; i < M; ++i) {
        phase_tra_38_cut[i] = phase_tra_38[size_t(az_st + i)];
        row_fa_cut[i] = row_fa[size_t(az_st + i)];
    }

    auto linfit = [](const std::vector<double>& x, const std::vector<double>& y, double& k, double& b) -> bool {
        if (x.size() != y.size() || x.size() < 2) return false;
        const size_t n = x.size();
        double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
        for (size_t i = 0; i < n; ++i) {
            Sx += x[i];
            Sy += y[i];
            Sxx += x[i] * x[i];
            Sxy += x[i] * y[i];
        }
        double det = n * Sxx - Sx * Sx;
        if (std::fabs(det) < 1e-12) return false;
        k = (n * Sxy - Sx * Sy) / det;
        b = (Sy * Sxx - Sx * Sxy) / det;
        return true;
    };

    double k = 0.0, b = 0.0;
    if (!linfit(row_fa_cut, phase_tra_38_cut, k, b)) return false;
    p_38 = {k, b};
    return true;
}

bool GMTIProcessor::clutter_cancel_38_paper_1(
    const std::vector<double>& y_faAxis,                  // len = cfg.pulse_num
    const std::vector<std::complex<double>>& F1f,         // size = Na*Nr
    const std::vector<std::complex<double>>& F2f,         // size = Na*Nr
    int az_st, int rg_st, int az_ed, int rg_ed,
    const Config& cfg,
    std::vector<std::complex<double>>& prosig_38,
    std::array<double,2>& p_38,
    std::vector<double>& phase_tra_38_cut,
    std::vector<double>& row_fa_cut
) {

    auto wall_start = std::chrono::high_resolution_clock::now();

    const size_t Na = static_cast<size_t>(cfg.pulse_num);
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    const size_t total = Na * Nr;

    if (F1f.size() != total || F2f.size() != total) return false;
    if (y_faAxis.size() != Na) return false;

    az_st = std::max(0, std::min(az_st, int(Na)-1));
    az_ed = std::max(0, std::min(az_ed, int(Na)-1));
    if (az_st > az_ed) return false;

    if (!clutter_cancel_38_paper_1_p38(y_faAxis, F1f, F2f, az_st, rg_st, az_ed, rg_ed, cfg, p_38)) {
        return false;
    }

    // —— 按行估计相位轨迹（加权：|F1|^2）⇢ 行并行
    std::vector<double> phase_tra_38(Na, 0.0), row_fa(Na, 0.0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int rr = az_st; rr <= az_ed; ++rr) {
        const size_t off = static_cast<size_t>(rr) * Nr;

        // 局部累加（避免使用复杂 reduction，C++11/OMP3.1 友好）
        double num_re = 0.0, num_im = 0.0, den = 0.0;
        for (size_t c = 0; c < Nr; ++c) {
            const auto& a = F1f[off + c];
            const auto& b = F2f[off + c];
            // a * conj(b)
            num_re += ( a.real()* b.real() + a.imag()* b.imag());
            num_im += ( a.imag()* b.real() - a.real()* b.imag());
            den    += std::norm(a);
        }
        const std::complex<double> num(num_re, num_im);
        phase_tra_38[static_cast<size_t>(rr)] = (den>0.0) ? std::arg(num/den) : std::arg(num);
        row_fa    [static_cast<size_t>(rr)] = y_faAxis[static_cast<size_t>(rr)];
    }

    // —— 解缠 + 基准平移（顺序执行）
    auto unwrap_inplace = [](std::vector<double>& v){
        for (size_t i=1;i<v.size();++i){
            double d=v[i]-v[i-1];
            if (d> M_PI) v[i]-=2*M_PI;
            if (d<-M_PI) v[i]+=2*M_PI;
        }
    };
    int mid_num = std::max(1, std::min(int(std::round(double(az_ed-az_st+1))), int(Na)));
    const double mid_phase = phase_tra_38[size_t(mid_num-1)];
    unwrap_inplace(phase_tra_38);
    const double delta = mid_phase - phase_tra_38[size_t(mid_num-1)];
    for (double& v: phase_tra_38) v += delta;

    // —— 线段拟合 phase ≈ k*y + b（顺序执行；数据量小）
    const int M = az_ed - az_st + 1;
    phase_tra_38_cut.resize(M);
    row_fa_cut.resize(M);
    for (int i=0;i<M;++i){
        phase_tra_38_cut[i] = phase_tra_38[size_t(az_st+i)];
        row_fa_cut[i]       = row_fa[size_t(az_st+i)];
    }
    auto linfit = [](const std::vector<double>& x,const std::vector<double>& y,double& k,double& b)->bool{
        if (x.size()!=y.size() || x.size()<2) return false;
        const size_t n=x.size();
        double Sx=0,Sy=0,Sxx=0,Sxy=0;
        for(size_t i=0;i<n;++i){ Sx+=x[i]; Sy+=y[i]; Sxx+=x[i]*x[i]; Sxy+=x[i]*y[i]; }
        double det = n*Sxx - Sx*Sx; if (std::fabs(det)<1e-12) return false;
        k = (n*Sxy - Sx*Sy)/det; b = (Sy*Sxx - Sx*Sxy)/det; return true;
    };
    double k=0,b=0; if (!linfit(row_fa_cut, phase_tra_38_cut, k, b)) return false;
    p_38 = {k,b};
#ifdef DEBUG
    {
        std::ofstream fout("debug_phase_38paper_linefit.bin", std::ios::binary);
        fout.write(reinterpret_cast<const char *>(phase_tra_38_cut.data()), phase_tra_38_cut.size() * sizeof(double));
        fout.close();

        std::ofstream fout2("debug_row_fa_linefit.bin", std::ios::binary);
        fout2.write(reinterpret_cast<const char *>(row_fa_cut.data()), row_fa_cut.size() * sizeof(double));
        fout2.close();
    }
#endif
    // —— 相位补偿 ch2 并幅度匹配 + 相消 ⇢ 行并行
    prosig_38.assign(total, std::complex<double>(0.0,0.0));

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t r=0;r<Na;++r){
        // exp(j*phi) 改为 cos/sin，少一次复指数
        const double phi = k*y_faAxis[r] + b;
        const double cs = std::cos(phi), sn = std::sin(phi);
        const std::complex<double> az_fai(cs, sn);
        const size_t off = r*Nr;

        for (size_t c=0;c<Nr;++c){
            const auto& a = F1f[off+c];
            const auto  b2= az_fai*F2f[off+c];
            const double aa= std::abs(a), bb= std::abs(b2);
            const double m = (aa<bb) ? aa : bb;

            const auto a_eq = (aa>0)? (m/aa)*a : std::complex<double>(0,0);
            const auto b_eq = (bb>0)? (m/bb)*b2: std::complex<double>(0,0);
            prosig_38[off+c] = a_eq - b_eq;
        }
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> wall_ms = wall_end - wall_start;
    // printf("[TIME] Total Wall Time (Inc. Memcpy): %.3f ms\n", wall_ms.count());

    return true;
}

