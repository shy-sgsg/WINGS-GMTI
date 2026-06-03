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
        return makePeriodListROI51(cfg.roi_ll_deg, plane, periodList, -25.0, 25.0);
    }
}

bool GMTIProcessor::makePeriodListROI51(const std::array<double,4>& roi_ll_deg,
                                        const GMTIOutput::Plane& plane,
                                        std::vector<int>& periodList,
                                        double scan_min_deg,
                                        double scan_max_deg)
{
    periodList.clear();

    // 1) ROI 中心（四角均值；若跨±180°可先做经度unwrap，这里略）
    double lat_c = 0.0, lng_c = 0.0;

    lat_c = roi_ll_deg[0];
    lng_c = roi_ll_deg[1];

    double E_c=0.0, N_c=0.0;
    Gaussp3(lat_c, lng_c, 117.0, E_c, N_c); // 参考经度117

    // 2) 观测方位角（飞机→ROI中心），坐标系：东为0°，逆时针为正
    double dE=0.0, dN=0.0;
    geoDiff_EN_m(plane.E, plane.N, E_c, N_c, dE, dN);
    const double bearing_deg = std::atan2(dN, dE) * 180.0 / 3.14159265358979323846;

    // 3) ★ 斜视角定义：以“飞行方向的右侧法线”为 0°，逆时针为正
    //    right-normal = heading - 90°
    //    squint = bearing - (heading - 90°) = bearing - heading + 90°
    const double squint_deg = wrap180(bearing_deg - plane.V_angle + 90.0);

    // 4) 扫描角域检查（默认 [-25,25]）
    const double a_min = std::min(scan_min_deg, scan_max_deg);
    const double a_max = std::max(scan_min_deg, scan_max_deg);
    if (a_max <= a_min) {
        // 配置异常：回默认中间 3 个
        periodList = {13, 14};
        return true;
    }
    if (squint_deg < a_min || squint_deg > a_max) {
        // 超域：回默认中间 3 个
        periodList = {13, 14};
        return true;
    }

    // 5) 将连续角 squint 映射到 1..51 的离散波位
    //    线性映射：-25° -> 1， +25° -> 51，0° -> 26
    //    公式：idx = round( ((squint - a_min) / (a_max - a_min)) * 50 ) + 1
    double frac = (squint_deg - a_min) / (a_max - a_min);   // [0,1]
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;

    int idx = static_cast<int>(std::floor(frac * 25.0 + 0.5)) + 1; // 1..51
    if (idx < 1)  idx = 1;
    if (idx > 26) idx = 26;

    // 6) 取毗邻 3 个波位（边界裁剪）
    const int i0 = std::max(1, idx - 1);
    const int i1 = idx;
    const int i2 = std::min(26, idx + 1);
    periodList = { i0, i1 };

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
    const int M = static_cast<int>(periodList.size());

    return true;
}