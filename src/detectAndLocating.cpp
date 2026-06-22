#include <deque>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include "GMTIProcessor.hpp"
#include "trig_lut.hpp"
#include <iostream>

bool GMTIProcessor::target_select(const std::vector<std::complex<double>> &GMTI_dataf_1,
                                  const std::vector<std::complex<double>> &GMTI_dataf_2,
                                  const std::vector<int> &prow,
                                  const std::vector<int> &pcol,
                                  const Config &cfg,
                                  GMTIOutput::Detect &S_target)
{
    const size_t Na = static_cast<size_t>(effectivePulseNum(cfg));
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    const size_t total = Na * Nr;
    if (GMTI_dataf_1.size() != total || GMTI_dataf_2.size() != total)
        return false;
    if (prow.size() != pcol.size())
        return false;

    struct Item
    {
        int r, c;
        double A;
        std::complex<double> ch1, ch2;
    };
    std::vector<Item> items;
    items.reserve(prow.size());

    // 收集有效点
    for (size_t i = 0; i < prow.size(); ++i)
    {
        int r = prow[i];
        int c = pcol[i];
        if (r < 0 || c < 0 || static_cast<size_t>(r) >= Na || static_cast<size_t>(c) >= Nr)
        {
            continue; // 跳过越界点
        }
        const size_t off = static_cast<size_t>(r) * Nr + static_cast<size_t>(c);
        const auto v1 = GMTI_dataf_1[off];
        const auto v2 = GMTI_dataf_2[off];
        items.push_back(Item{r, c, std::abs(v1), v1, v2});
    }

    // 按幅度降序排序（等价 MATLAB 的 sortrows(...,3,'descend')）
    std::sort(items.begin(), items.end(),
              [](const Item &a, const Item &b)
              { return a.A > b.A; });

    // 填充输出
    S_target.prow.clear();
    S_target.pcol.clear();
    S_target.prow.reserve(items.size());
    S_target.pcol.reserve(items.size());

    for (const auto &it : items)
    {
        S_target.prow.push_back(it.r);
        S_target.pcol.push_back(it.c);
    }

    return true;
}

bool GMTIProcessor::dpca_cfar2_fast(const std::vector<std::complex<double>> &GMTI_new,
                                    double pf, int c_num, int b_num, const std::string &type,
                                    const Config &cfg,
                                    std::vector<double> &mydata,
                                    std::vector<int> &prow,
                                    std::vector<int> &pcol)
{
    const int H = effectivePulseNum(cfg); // 行（方位）
    const int W = cfg.rg_len;    // 列（距离）
    if (H <= 0 || W <= 0)
        return false;
    if ((int)GMTI_new.size() != H * W)
    {
        std::cerr << "dpca_cfar2: 输入尺寸与 cfg 不一致\n";
        return false;
    }
    if (!(pf > 0.0 && pf < 1.0))
    {
        std::cerr << "dpca_cfar2: pf 应在 (0,1) 内\n";
        return false;
    }

    const int g = std::max(0, c_num); // 保护单元半径
    const int b = std::max(1, b_num); // 背景带厚度
    const int R = g + b;              // 全窗半径
    const int win = 2 * R + 1;        // 全窗边长
    const int cut = 2 * g + 1;        // 保护(含CUT)边长
    const int total_bg = win * win - cut * cut;
    if (total_bg <= 0)
    {
        std::cerr << "dpca_cfar2: 背景单元数非法\n";
        return false;
    }

    // 统一 CFAR 类型
    std::string ty = type;
    for (auto &ch : ty)
        ch = (char)std::toupper((unsigned char)ch);
    if (ty != "CA" && ty != "GO" && ty != "OS")
        ty = "CA";
    if (ty == "OS")
    {
        std::cerr << "dpca_cfar2: 暂未实现 OS 的积分图优化，回退到 CA。\n";
        ty = "CA";
    }

    // —— 1) 构建积分图 I (H+1)*(W+1)，按行累计；功率用 double
    std::vector<double> I((size_t)(H + 1) * (W + 1), 0.0);
    auto atI = [&](int r, int c) -> double &
    { return I[(size_t)r * (W + 1) + (size_t)c]; };

    for (int r = 1; r <= H; ++r)
    {
        double rowsum = 0.0;
        const int row0 = r - 1;
        for (int c = 1; c <= W; ++c)
        {
            const int col0 = c - 1;
            const size_t idx = (size_t)row0 * W + (size_t)col0;
            const double re = GMTI_new[idx].real();
            const double im = GMTI_new[idx].imag();
            const double p = re * re + im * im; // double 功率
            rowsum += p;
            atI(r, c) = atI(r - 1, c) + rowsum;
        }
    }

    // 矩形求和：闭区间 [r0,r1]×[c0,c1]（0-based, 含端点）
    auto rect_sum = [&](int r0, int c0, int r1, int c1) -> double
    {
        return atI(r1 + 1, c1 + 1) - atI(r0, c1 + 1) - atI(r1 + 1, c0) + atI(r0, c0);
    };

    // —— 2) 预计算 alpha
    const double alpha = double(total_bg) * (std::pow(pf, -1.0 / double(total_bg)) - 1.0);

    // —— 3) 遍历检测（边缘跳过 R）
    mydata.assign((size_t)H * (size_t)W, 0.0);
    prow.clear();
    pcol.clear();
    prow.reserve((size_t)(H * W / 200));
    pcol.reserve((size_t)(H * W / 200));

    const int r_min = R, r_max = H - 1 - R;
    const int c_min = R, c_max = W - 1 - R;

    for (int r = r_min; r <= r_max; ++r)
    {
        for (int c = c_min; c <= c_max; ++c)
        {
            const size_t idx = (size_t)r * W + (size_t)c;
            const double re = GMTI_new[idx].real();
            const double im = GMTI_new[idx].imag();
            const double CUT = re * re + im * im; // 再算一次 CUT（double）

            // 全窗与保护窗边界（闭区间）
            const int r0 = r - R, r1 = r + R;
            const int c0 = c - R, c1 = c + R;

            const int rg0 = r - g, rg1 = r + g;
            const int cg0 = c - g, cg1 = c + g;

            double noise_level = 0.0;

            if (ty == "CA")
            {
                const double sum_full = rect_sum(r0, c0, r1, c1);
                const double sum_guard = rect_sum(rg0, cg0, rg1, cg1);
                const double bg_sum = sum_full - sum_guard;
                noise_level = bg_sum / double(total_bg);
            }
            else
            { // "GO"
                // 四向带（从全窗减去保护窗），取四个方向均值的最大值
                // 左带：行 [r0,r1], 列 [c0, cg0-1]
                double mleft = 0.0, mright = 0.0, mtop = 0.0, mbot = 0.0;

                const int lc0 = c0, lc1 = cg0 - 1;
                const int rc0 = cg1 + 1, rc1 = c1;
                const int tr0 = r0, tr1 = rg0 - 1;
                const int br0 = rg1 + 1, br1 = r1;

                if (lc0 <= lc1)
                {
                    const double s = rect_sum(r0, lc0, r1, lc1);
                    const int npx = (r1 - r0 + 1) * (lc1 - lc0 + 1);
                    if (npx > 0)
                        mleft = s / npx;
                }
                if (rc0 <= rc1)
                {
                    const double s = rect_sum(r0, rc0, r1, rc1);
                    const int npx = (r1 - r0 + 1) * (rc1 - rc0 + 1);
                    if (npx > 0)
                        mright = s / npx;
                }
                if (tr0 <= tr1)
                {
                    const double s = rect_sum(tr0, c0, tr1, c1);
                    const int npx = (tr1 - tr0 + 1) * (c1 - c0 + 1);
                    if (npx > 0)
                        mtop = s / npx;
                }
                if (br0 <= br1)
                {
                    const double s = rect_sum(br0, c0, br1, c1);
                    const int npx = (br1 - br0 + 1) * (c1 - c0 + 1);
                    if (npx > 0)
                        mbot = s / npx;
                }

                noise_level = std::max(std::max(mleft, mright), std::max(mtop, mbot));
            }

            const double thr = alpha * noise_level;
            if (CUT > thr)
            {
               if (r < 70 || r > 190)
                {
                    mydata[idx] = CUT; // ★ 初始策略：命中写功率
                    prow.push_back(r); // 0-based
                    pcol.push_back(c);
                }
            }
        }
    }
    return true;
}

bool GMTIProcessor::cluster_filter(const std::vector<double> &mydata,
                                   int min_points,
                                   const Config &cfg,
                                   std::vector<double> &refined_data,
                                   std::vector<int> &prow_new,
                                   std::vector<int> &pcol_new)
{
    const int H = effectivePulseNum(cfg);
    const int W = cfg.rg_len;
    if (H <= 0 || W <= 0)
        return false;
    const size_t total = static_cast<size_t>(H) * static_cast<size_t>(W);
    if (mydata.size() != total)
        return false;

    if (min_points < 1)
        min_points = 1;

    refined_data.assign(total, 0.0);
    prow_new.clear();
    pcol_new.clear();

    // 二值图：功率>0 视为候选
    std::vector<uint8_t> bin(total, 0);
    for (size_t i = 0; i < total; ++i)
        bin[i] = (mydata[i] > 0.0) ? 1 : 0;

    // 8 邻域偏移
    const int dr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const int dc[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

    std::vector<uint8_t> vis(total, 0);
    std::deque<size_t> q;
    std::vector<size_t> comp;
    comp.reserve(256);

    for (int r0 = 0; r0 < H; ++r0)
    {
        for (int c0 = 0; c0 < W; ++c0)
        {
            const size_t idx0 = static_cast<size_t>(r0) * W + static_cast<size_t>(c0);
            if (!bin[idx0] || vis[idx0])
                continue;

            // BFS 收集一个连通域
            comp.clear();
            vis[idx0] = 1;
            q.clear();
            q.push_back(idx0);

            while (!q.empty())
            {
                const size_t idx = q.front();
                q.pop_front();
                comp.push_back(idx);
                const int r = static_cast<int>(idx / W);
                const int c = static_cast<int>(idx % W);

                for (int k = 0; k < 8; ++k)
                {
                    const int rr = r + dr[k];
                    const int cc = c + dc[k];
                    if (rr < 0 || rr >= H || cc < 0 || cc >= W)
                        continue;
                    const size_t idn = static_cast<size_t>(rr) * W + static_cast<size_t>(cc);
                    if (!bin[idn] || vis[idn])
                        continue;
                    vis[idn] = 1;
                    q.push_back(idn);
                }
            }

            // 面积阈值
            if (static_cast<int>(comp.size()) >= min_points)
            {
                // 保留该区域原图
                for (size_t id : comp)
                    refined_data[id] = mydata[id];

                // 选择幅度最大的像素（mydata 为功率，最大功率即最大幅度）
                size_t best_id = comp[0];
                double best_val = mydata[best_id];
                for (size_t id : comp)
                {
                    if (mydata[id] > best_val)
                    {
                        best_val = mydata[id];
                        best_id = id;
                    }
                }
                const int br = static_cast<int>(best_id / W);
                const int bc = static_cast<int>(best_id % W);
                prow_new.push_back(br); // 0-based
                pcol_new.push_back(bc);
            }
        }
    }

    return true;
}

#include <vector>
#include <deque>
#include <cmath>    // cos, sin, atan2, sqrt

// 小工具：把角度包装到 (-pi, pi]
static inline double wrapAngleRad(double x)
{
    const double PI     = 3.14159265358979323846;
    const double TWO_PI = 6.28318530717958647692;
    while (x <= -PI) x += TWO_PI;
    while (x >  PI)  x -= TWO_PI;
    return x;
}

bool GMTIProcessor::cluster_filter_gap_phase(const std::vector<double> &mydata,
                                             const std::vector<double> &phase_map, // 干涉相位 (rad)
                                             int min_points,
                                             int max_gap,        // 允许的最大列向间断像元数
                                             double maxPhaseStd, // 相位标准差阈值 (rad)
                                             const Config &cfg,
                                             std::vector<double> &refined_data,
                                             std::vector<int> &prow_new,
                                             std::vector<int> &pcol_new,
                                             std::vector<double> &phaseStdList)
{
    const int H = effectivePulseNum(cfg);  // 行数（方位向）
    const int W = cfg.rg_len;     // 列数（距离向）
    if (H <= 0 || W <= 0)
        return false;

    const size_t total = static_cast<size_t>(H) * static_cast<size_t>(W);
    if (mydata.size() != total)
        return false;
    if (!phase_map.empty() && phase_map.size() != total)
        return false;

    if (min_points < 1)
        min_points = 1;
    if (max_gap < 1)
        max_gap = 1;  // max_gap=1 时退化为普通 8 邻域（列方向不扩展）

    refined_data.assign(total, 0.0);
    prow_new.clear();
    pcol_new.clear();
    phaseStdList.clear();

    // 二值图：功率>0 视为候选
    std::vector<uint8_t> bin(total, 0);
    for (size_t i = 0; i < total; ++i)
        bin[i] = (mydata[i] > 0.0) ? 1 : 0;

    std::vector<uint8_t> vis(total, 0);
    std::deque<size_t> q;
    std::vector<size_t> comp;
    comp.reserve(256);

    const bool usePhase = (!phase_map.empty() && maxPhaseStd > 0.0);

    for (int r0 = 0; r0 < H; ++r0)
    {
        for (int c0 = 0; c0 < W; ++c0)
        {
            const size_t idx0 = static_cast<size_t>(r0) * W + static_cast<size_t>(c0);
            if (!bin[idx0] || vis[idx0])
                continue;

            // ---------- BFS 收集一个“列向允许 max_gap 间断”的连通域 ----------
            comp.clear();
            vis[idx0] = 1;
            q.clear();
            q.push_back(idx0);

            while (!q.empty())
            {
                const size_t idx = q.front();
                q.pop_front();
                comp.push_back(idx);

                const int r = static_cast<int>(idx / W);
                const int c = static_cast<int>(idx % W);

                // 行方向 dr ∈ {-1, 0, 1}（不允许大间断）
                // 列方向 dc ∈ [-max_gap, max_gap]（允许最大 max_gap 间断）
                for (int dr = -1; dr <= 1; ++dr)
                {
                    const int rr = r + dr;
                    if (rr < 0 || rr >= H)
                        continue;

                    for (int dc = -max_gap; dc <= max_gap; ++dc)
                    {
                        if (dr == 0 && dc == 0)
                            continue;

                        const int cc = c + dc;
                        if (cc < 0 || cc >= W)
                            continue;

                        const size_t idn = static_cast<size_t>(rr) * W + static_cast<size_t>(cc);
                        if (!bin[idn] || vis[idn])
                            continue;

                        vis[idn] = 1;
                        q.push_back(idn);
                    }
                }
            }

            // ---------- (1) 面积阈值 ----------
            if (static_cast<int>(comp.size()) < min_points)
                continue;

            // ---------- (2) 干涉相位一致性判据 ----------
            double phaseStd = 0.0;
            if (usePhase)
            {
                // 圆统计：求平均相位方向
                double sumCos = 0.0;
                double sumSin = 0.0;

                for (size_t k = 0; k < comp.size(); ++k)
                {
                    const size_t id = comp[k];
                    const double phi = phase_map[id]; // 弧度
                    sumCos += gmti::trig_lut::cos(phi);
                    sumSin += gmti::trig_lut::sin(phi);
                }

                const double phi0 = gmti::trig_lut::atan2(sumSin, sumCos); // 平均相位方向

                // 把相位绕 phi0 展开，计算标准差
                double mean  = 0.0;
                double M2    = 0.0;
                size_t count = 0;

                for (size_t k = 0; k < comp.size(); ++k)
                {
                    const size_t id = comp[k];
                    double dphi = phase_map[id] - phi0;
                    dphi = wrapAngleRad(dphi);  // 映射到 (-pi, pi]

                    ++count;
                    const double delta  = dphi - mean;
                    mean += delta / static_cast<double>(count);
                    const double delta2 = dphi - mean;
                    M2 += delta * delta2;
                }

                if (count > 1)
                {
                    const double var = M2 / static_cast<double>(count - 1);
                    phaseStd = (var > 0.0) ? std::sqrt(var) : 0.0;
                }
                else
                {
                    phaseStd = 0.0;
                }

                // 相位波动太大 → 判为虚警
                if (phaseStd > maxPhaseStd)
                    continue;
            }

            // ---------- (3) 真正保留该目标 ----------
            // 把该连通域内的功率写回 refined_data
            for (size_t k = 0; k < comp.size(); ++k)
            {
                const size_t id = comp[k];
                refined_data[id] = mydata[id];
            }

            // 代表点：功率最大像元
            size_t best_id  = comp[0];
            double best_val = mydata[best_id];

            for (size_t k = 1; k < comp.size(); ++k)
            {
                const size_t id = comp[k];
                const double v  = mydata[id];
                if (v > best_val)
                {
                    best_val = v;
                    best_id  = id;
                }
            }

            const int br = static_cast<int>(best_id / W);
            const int bc = static_cast<int>(best_id % W);

            prow_new.push_back(br);       // 0-based
            pcol_new.push_back(bc);
            phaseStdList.push_back(phaseStd);
        }
    }

    return true;
}
