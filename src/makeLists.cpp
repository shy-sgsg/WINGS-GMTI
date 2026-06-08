#include "GMTIProcessor.hpp"
#include "geo/geoProj.hpp" // 包含 wrap180

#include <cmath>
#include <algorithm>
#include <vector>
#include <array>

bool GMTIProcessor::makePeriodList(const GMTIOutput::Plane& plane,
                                   const Config& cfg,
                                   std::vector<int>& periodList)
{
    if (!cfg.wavepos_use_roi) {
        return makeDefaultRangeList(cfg, periodList);
    } else {
        return makePeriodListROI51(cfg, plane, periodList);
    }
}

bool GMTIProcessor::makePeriodListROI51(const Config& cfg,
                                        const GMTIOutput::Plane& plane,
                                        std::vector<int>& periodList)
{
    periodList.clear();

    // 1) ROI 中心（四角均值；若跨±180°可先做经度unwrap，这里略）
    double lat_c = 0.0, lng_c = 0.0;

    lat_c = cfg.roi_ll_deg[0];
    lng_c = cfg.roi_ll_deg[1];

    double E_c=0.0, N_c=0.0;
    Gaussp3(lat_c, lng_c, cfg.L0, E_c, N_c);

    // 2) 观测方位角（飞机→ROI中心），坐标系：东为0°，逆时针为正
    double dE=0.0, dN=0.0;
    geoDiff_EN_m(plane.E, plane.N, E_c, N_c, dE, dN);
    const double bearing_deg = std::atan2(dN, dE) * 180.0 / 3.14159265358979323846;

    // 3) ★ 斜视角定义：以“飞行方向的右侧法线”为 0°，逆时针为正
    //    right-normal = heading - 90°
    //    squint = bearing - (heading - 90°) = bearing - heading + 90°
    const double squint_deg = wrap180(bearing_deg - plane.V_angle + 90.0);

    // 4) 扫描角域检查（默认 [-25,25]）
    const double a_min = std::min(cfg.scan_min_deg, cfg.scan_max_deg);
    const double a_max = std::max(cfg.scan_min_deg, cfg.scan_max_deg);
    const int p_first = std::min(cfg.wavepos_st, cfg.wavepos_ed);
    const int p_last = std::max(cfg.wavepos_st, cfg.wavepos_ed);
    const int p_count = p_last - p_first + 1;
    const int p_step = std::max(1, std::abs(cfg.wavepos_skip));
    if (a_max <= a_min || p_count <= 0) {
        return makeDefaultRangeList(cfg, periodList);
    }

    std::vector<int> candidates;
    candidates.reserve(static_cast<size_t>((p_count + p_step - 1) / p_step));
    for (int p = p_first; p <= p_last; p += p_step) {
        candidates.push_back(p);
    }
    if (candidates.empty()) {
        return false;
    }

    // 5) 按配置的扫描角范围和波位范围建立线性映射：
    //    p_first -> scan_min_deg，p_last -> scan_max_deg。
    const auto angleOfPeriod = [&](int p) -> double {
        if (p_count == 1) {
            return 0.5 * (a_min + a_max);
        }
        const double frac = static_cast<double>(p - p_first) / static_cast<double>(p_count - 1);
        return a_min + frac * (a_max - a_min);
    };

    // 6) 在按 wavepos_skip 抽样后的候选波位中，取与 ROI 斜视角最近的波位及其左右邻居。
    size_t best = 0;
    double best_err = std::abs(angleOfPeriod(candidates[0]) - squint_deg);
    for (size_t i = 1; i < candidates.size(); ++i) {
        const double err = std::abs(angleOfPeriod(candidates[i]) - squint_deg);
        if (err < best_err) {
            best = i;
            best_err = err;
        }
    }

    const size_t begin = (best == 0) ? 0 : best - 1;
    const size_t end = std::min(candidates.size() - 1, best + 1);
    for (size_t i = begin; i <= end; ++i) {
        periodList.push_back(candidates[i]);
    }

    return true;
}

bool GMTIProcessor::makeDefaultRangeList(const Config& cfg, std::vector<int>& periodList)
{
    periodList.clear();
    const int st   = cfg.wavepos_st;
    const int ed   = cfg.wavepos_ed;
    const int step = (cfg.wavepos_skip == 0 ? 1 : cfg.wavepos_skip);
    if (ed < st || step <= 0) return false;

    for (int v = st; v <= ed; v += step) periodList.push_back(v);

    return true;
}
