#include "GMTIProcessor.hpp"


// 角度->弧度
static inline double deg2rad(double deg) { return deg * M_PI / 180.0; }

// 在已升序的 faAxis 上找与 val 最接近的下标（0-based）
// 若 faAxis 为空返回 -1
static inline int nearest_index(const std::vector<double>& faAxis, double val) {
    if (faAxis.empty()) return -1;
    auto it = std::lower_bound(faAxis.begin(), faAxis.end(), val);
    if (it == faAxis.begin()) return 0;
    if (it == faAxis.end())   return static_cast<int>(faAxis.size() - 1);
    // 夹在 prev(it) 与 it 之间，取更近者
    auto iR = static_cast<int>(it - faAxis.begin());
    auto iL = iR - 1;
    return (std::abs(faAxis[iL] - val) <= std::abs(faAxis[iR] - val)) ? iL : iR;
}

/**
 * 动态支撑域计算（inline）
 * 等价 MATLAB：
 *   BW_az = 2 * plane.V * sind(2) / cfg.lambda;
 *   fd_st = fa2 - BW_az; az_center = argmin|faAxis-fa2|;
 *   fd_ed = fa2 + BW_az; az_st/az_ed 同理
 *
 * 参数：
 *   faAxis   : 频率轴（长度 ~ pulse_num）
 *   fa2      : 多普勒中心
 *   plane    : 需包含 V (m/s)
 *   cfg      : 需包含 lambda, rg_st, rg_ed
 * 输出：
 *   az_center, az_st, az_ed : 方位索引(0-based)
 *   fd_st, fd_ed, BW_az     : 便于后续使用的频率量
 *   rg_st_out, rg_ed_out    : 直接透传 cfg 的距离支撑域（方便调用端）
 */
bool GMTIProcessor::computeDynamicSupportDomain(
        const std::vector<double>& faAxis,
        double fa2,
        const GMTIOutput::Plane& plane,
        const Config& cfg,
        int& az_center,
        int& az_st,
        int& az_ed,
        double& fd_st,
        double& fd_ed,
        double& BW_az,
        int& rg_st_out,
        int& rg_ed_out)
{
    if (faAxis.empty() || cfg.lambda == 0.0) return false;

    // 距离支撑域直接来自 cfg
    rg_st_out = cfg.rg_st;
    rg_ed_out = cfg.rg_ed;

    // 带宽（注意 sind(2) -> sin(2°)）
    BW_az = 2.0 * plane.V * std::sin(deg2rad(2.0)) / cfg.lambda;

    // 动态上下边界频率
    fd_st = fa2 - BW_az;
    fd_ed = fa2 + BW_az;

    // 最近邻索引（0-based）
    az_center = nearest_index(faAxis, fa2);
    az_st     = nearest_index(faAxis, fd_st);
    az_ed     = nearest_index(faAxis, fd_ed);

    // 容错：若任一未找到
    if (az_center < 0 || az_st < 0 || az_ed < 0) return false;

    return true;
}
