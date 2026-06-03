#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include "GMTIProcessor.hpp"

// ---- 小工具：3x3 线性方程求解 A x = b（高斯消元） ------------------------
static bool solve3x3(double A[3][3], double b[3], double x[3]) {
    // 增广矩阵
    double M[3][4] = {
        {A[0][0], A[0][1], A[0][2], b[0]},
        {A[1][0], A[1][1], A[1][2], b[1]},
        {A[2][0], A[2][1], A[2][2], b[2]}
    };
    // 消元
    for (int i = 0; i < 3; ++i) {
        // 选主元
        int piv = i;
        for (int r = i+1; r < 3; ++r)
            if (std::fabs(M[r][i]) > std::fabs(M[piv][i])) piv = r;
        if (std::fabs(M[piv][i]) < 1e-12) return false;
        if (piv != i) for (int c = i; c < 4; ++c) std::swap(M[i][c], M[piv][c]);

        // 归一化
        double inv = 1.0 / M[i][i];
        for (int c = i; c < 4; ++c) M[i][c] *= inv;

        // 消去其他行
        for (int r = 0; r < 3; ++r) if (r != i) {
            double f = M[r][i];
            for (int c = i; c < 4; ++c) M[r][c] -= f * M[i][c];
        }
    }
    x[0] = M[0][3]; x[1] = M[1][3]; x[2] = M[2][3];
    return true;
}

// ---- 二次多项式最小二乘拟合：y ≈ a x^2 + b x + c，x,y 为等长向量 -----------
static bool fitQuadraticLeastSquares(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     double coeff[3]) {
    const size_t n = x.size();
    if (n < 3 || y.size() != n) return false;

    // 构建正规方程 A * p = b
    double S0= n;
    double S1=0, S2=0, S3=0, S4=0;
    double T0=0, T1=0, T2=0; // b 向量：Σ(y), Σ(x*y), Σ(x^2*y)

    for (size_t i = 0; i < n; ++i) {
        double xi = x[i], yi = y[i];
        double x2 = xi*xi, x3 = x2*xi, x4 = x2*x2;
        S1 += xi; S2 += x2; S3 += x3; S4 += x4;
        T0 += yi; T1 += xi*yi; T2 += x2*yi;
    }

    // A 对应于 [x^2 x 1] 的内积
    double A[3][3] = {
        { S4, S3, S2 },
        { S3, S2, S1 },
        { S2, S1, S0 }
    };
    double b[3] = { T2, T1, T0 };
    return solve3x3(A, b, coeff);
}

// ---- RANSAC + 二次拟合（等价 MATLAB 的 RANSAC2） -------------------------
static bool ransacQuadratic(const std::vector<int>&    xs_1based,   // 距离下标（1-based，仿 MATLAB）
                            const std::vector<double>& ys_phase,    // 相位
                            double thre,                            // 阈值(弧度)
                            double coeff_out[3])                    // 最终系数 a,b,c
{
    const int len_t = static_cast<int>(ys_phase.size());
    if (len_t < 3 || (int)xs_1based.size() != len_t) return false;

    const int final_count = 40;
    const int use_num     = std::min(10, len_t);

    std::mt19937 rng{ std::random_device{}() };
    std::uniform_int_distribution<int> uni(0, len_t - 1);

    double best_ratio = -1.0;
    std::vector<double> best_x, best_y;

    for (int it = 0; it < final_count; ++it) {
        // 随机子集
        std::vector<int>    idx(use_num);
        for (int k = 0; k < use_num; ++k) idx[k] = uni(rng);

        // 构建子集并拟合（X=[x^2,x,1]）
        std::vector<double> xs(use_num), ys(use_num);
        for (int k = 0; k < use_num; ++k) {
            xs[k] = static_cast<double>(xs_1based[idx[k]]);
            ys[k] = ys_phase[idx[k]];
        }

        double P[3];
        if (!fitQuadraticLeastSquares(xs, ys, P)) continue;

        // 统计内点
        std::vector<double> cur_x_in, cur_y_in;
        cur_x_in.reserve(len_t);
        cur_y_in.reserve(len_t);
        int inlier = 0;
        for (int i = 0; i < len_t; ++i) {
            double xi = static_cast<double>(xs_1based[i]);
            double yi = ys_phase[i];
            double yhat = P[0]*xi*xi + P[1]*xi + P[2];
            if (std::fabs(yi - yhat) < thre) {
                ++inlier;
                cur_x_in.push_back(xi);
                cur_y_in.push_back(yi);
            }
        }
        double ratio = static_cast<double>(inlier) / len_t;
        if (ratio > best_ratio) {
            best_ratio = ratio;
            best_x.swap(cur_x_in);
            best_y.swap(cur_y_in);
        }
    }

    // 若没有内点集，退化用全体点拟合
    if (best_x.size() < 3) {
        best_x.resize(len_t);
        best_y.resize(len_t);
        for (int i = 0; i < len_t; ++i) {
            best_x[i] = static_cast<double>(xs_1based[i]);
            best_y[i] = ys_phase[i];
        }
    }

    return fitQuadraticLeastSquares(best_x, best_y, coeff_out);
}

/// ---- 主函数：rg_correct ---------------------------------------------------
bool GMTIProcessor::rg_correct(const std::vector<std::complex<double>>& F1,
                               const std::vector<std::complex<double>>& F2,
                               const Config& cfg,
                               double th_rg,
                               std::vector<double>& phi_fit,
                               std::vector<double>& phi_diss_phase,
                               std::vector<int>&    phi_diss_range)
{
    const int Na = cfg.pulse_num;   // 使用类内 cfg_（你之前要求 readXmlParam 填充到私有成员）
    const int Nr = cfg.rg_len;
    const int rg_len = cfg.rg_len;

    if (Na <= 0 || Nr <= 0) return false;
    const size_t total = static_cast<size_t>(Na) * static_cast<size_t>(Nr);
    if (F1.size() != total || F2.size() != total) return false;

    // 范围校验 & 规范
    int rg_st = std::max(0, std::min(cfg.rg_st, Nr-1));
    int rg_ed = std::max(0, std::min(cfg.rg_ed, Nr-1));
    int az_st = std::max(0, std::min(cfg.az_st, Na-1));
    int az_ed = std::max(0, std::min(cfg.az_ed, Na-1));
    if (rg_st > rg_ed || az_st > az_ed) return false;

    if (th_rg <= 0.0) th_rg = 0.2 * M_PI;

    const int R = rg_ed - rg_st + 1;

    // printf("Na=%d, Nr=%d\n", Na, Nr);
    // printf("rg_st=%d, rg_ed=%d\n", rg_st, rg_ed);
    // printf("az_st=%d, az_ed=%d\n", az_st, az_ed);
    // printf("R=%d, th_rg=%.6f\n", R, th_rg);

    // 计算每个距离单元的相位：angle(sum_{az}(F1 .* conj(F2)))
    phi_diss_phase.assign(R, 0.0);
    phi_diss_range.resize(R);
    for (int c = rg_st, k = 0; c <= rg_ed; ++c, ++k) {
        std::complex<double> s(0.0, 0.0);
        for (int r = az_st; r <= az_ed; ++r) {
            const size_t off = static_cast<size_t>(r) * Nr + static_cast<size_t>(c);
            s += F1[off] * std::conj(F2[off]);
        }
        phi_diss_phase[k] = std::arg(s);
        // 为了与 MATLAB phi_diss(2,:) = rg_st:rg_ed 一致，用 1-based 存储
        phi_diss_range[k] = c + 1;
        // if(k<5) printf("phi_diss_phase[%d] = %f\n", k, phi_diss_phase[k]);
    }

    // RANSAC + 二次拟合，得到系数 a,b,c
    double P[3] = {0,0,0};
    if (!ransacQuadratic(phi_diss_range, phi_diss_phase, th_rg, P)) {
        // 拟合失败
        return false;
    }

    // 生成 phi_fit（与 MATLAB 保持：x=1..rg_len）
    phi_fit.resize(std::max(0, rg_len));
    for (int i = 0; i < rg_len; ++i) {
        double x = static_cast<double>(i + 1); // 1-based
        phi_fit[i] = P[0]*x*x + P[1]*x + P[2];
    }
    return true;
}

bool GMTIProcessor::rg_correct_CUDA(
    const Config& cfg,             // 配置参数
    double th_rg,                  // RANSAC 阈值
    std::vector<float>& phi_fit,   // 输出：拟合曲线
    std::vector<double>& phi_diss_phase, 
    std::vector<int>&    phi_diss_range)
{
    const int Na = cfg.pulse_num;
    const int Nr = cfg.rg_len;

    // 范围校验 & 规范
    int rg_st = std::max(0, std::min(cfg.rg_st, Nr-1));
    int rg_ed = std::max(0, std::min(cfg.rg_ed, Nr-1));
    int az_st = std::max(0, std::min(cfg.az_st, Na-1));
    int az_ed = std::max(0, std::min(cfg.az_ed, Na-1));
    if (rg_st > rg_ed || az_st > az_ed) return false;
    
    const int R = rg_ed - rg_st + 1;

    // printf("Na=%d, Nr=%d\n", Na, Nr);
    // printf("rg_st=%d, rg_ed=%d\n", rg_st, rg_ed);
    // printf("az_st=%d, az_ed=%d\n", az_st, az_ed);
    // printf("R=%d, th_rg=%.6f\n", R, th_rg);

    if (R <= 0) return false;

    // --- Step 1 & 2: GPU 累加与结果下载 ---
    // 调用我们刚才写的两个 GPU 辅助函数
    std::vector<std::complex<float>> h_sums;
    if (!cuda_stage_rg_sum_async(cfg)) return false;
    if (!cuda_download_rg_sums_sync(h_sums, Nr)) return false;

    // --- Step 3: CPU 计算相位 (std::arg) ---
    phi_diss_phase.assign(R, 0.0);
    phi_diss_range.resize(R);
    
    for (int c = cfg.rg_st, k = 0; c <= cfg.rg_ed; ++c, ++k) {
        // 直接对 GPU 传回的累加结果求相位
        phi_diss_phase[k] = static_cast<double>(std::arg(h_sums[c]));
        // 保持 1-based 坐标
        phi_diss_range[k] = c + 1;
        // if(k<5) printf("phi_diss_phase_gpu[%d] = %f\n", k, phi_diss_phase[k]);
    }

    // --- Step 4: RANSAC + 二次拟合 ---
    double P[3] = {0, 0, 0};
    if (!ransacQuadratic(phi_diss_range, phi_diss_phase, th_rg, P)) {
        // 如果拟合失败，通常是因为杂波区相位太乱
        return false;
    }

    // --- Step 5: 生成全距离向补偿曲线 ---
    phi_fit.resize(Nr);
    for (int i = 0; i < Nr; ++i) {
        double x = static_cast<double>(i + 1); // 1-based
        phi_fit[i] = static_cast<float>(P[0]*x*x + P[1]*x + P[2]);
    }

    return true;
}