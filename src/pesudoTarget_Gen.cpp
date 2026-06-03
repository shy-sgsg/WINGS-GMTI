#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>
#include <geo/geoProj.hpp>
#include <config_structs.hpp> // 包含 Config / GMTIOutput::Result 声明
#include "GMTIProcessor.hpp"
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <strings.h> // strcasecmp
#include <cctype>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>

struct RegionBounds
{
    double minLat, maxLat, minLon, maxLon;
    RegionBounds() : minLat(0), maxLat(0), minLon(0), maxLon(0) {}
};

// --- 工具 ---
static inline bool ends_with_txt_ci(const char *name)
{
    size_t n = std::strlen(name);
    if (n < 4)
        return false;
    const char *ext = name + (n - 4);
    return (strcasecmp(ext, ".txt") == 0);
}
static std::string trim_copy(const std::string &s)
{
    size_t i = 0, j = s.size();
    while (i < j && std::isspace((unsigned char)s[i]))
        ++i;
    while (j > i && std::isspace((unsigned char)s[j - 1]))
        --j;
    return s.substr(i, j - i);
}
static bool join_path(const std::string &dir, const char *name, std::string &out)
{
    if (dir.empty())
        return false;
    out = (dir[dir.size() - 1] == '/') ? (dir + name) : (dir + "/" + name);
    return true;
}

// --- 找目录下“最近修改”的 .txt ---
static bool find_latest_txt(const std::string &dir, std::string &latest_path)
{
    DIR *dp = opendir(dir.c_str());
    if (!dp)
        return false;

    bool found = false;
    time_t best_mtime = 0;
    struct dirent *ent;
    while ((ent = readdir(dp)) != NULL)
    {
        if (ent->d_name[0] == '.')
            continue;
        if (!ends_with_txt_ci(ent->d_name))
            continue;

        std::string path;
        if (!join_path(dir, ent->d_name, path))
            continue;

        struct stat st;
        if (stat(path.c_str(), &st) != 0)
            continue;
        if (!S_ISREG(st.st_mode))
            continue;

        if (!found || st.st_mtime > best_mtime)
        {
            found = true;
            best_mtime = st.st_mtime;
            latest_path = path;
        }
    }
    closedir(dp);
    return found;
}

// --- 解析 B0..B3, L0..L3 -> 包围盒 ---
static bool parse_region_bounds_from_file(const std::string &path, RegionBounds &out)
{
    std::ifstream fin(path.c_str());
    if (!fin)
        return false;

    double B[4] = {0, 0, 0, 0}, L[4] = {0, 0, 0, 0};
    bool Bhit[4] = {0, 0, 0, 0}, Lhit[4] = {0, 0, 0, 0};

    std::string line;
    while (std::getline(fin, line))
    {
        size_t pos = line.find('=');
        if (pos == std::string::npos)
            continue;
        std::string k = trim_copy(line.substr(0, pos));
        std::string v = trim_copy(line.substr(pos + 1));
        if (k.size() == 2 && (k[0] == 'B' || k[0] == 'L') && std::isdigit((unsigned char)k[1]))
        {
            int idx = k[1] - '0';
            if (idx >= 0 && idx <= 3)
            {
                char *endp = NULL;
                double val = std::strtod(v.c_str(), &endp);
                if (endp && *endp == '\0')
                {
                    if (k[0] == 'B')
                    {
                        B[idx] = val;
                        Bhit[idx] = true;
                    }
                    else
                    {
                        L[idx] = val;
                        Lhit[idx] = true;
                    }
                }
            }
        }
    }
    for (int i = 0; i < 4; ++i)
        if (!Bhit[i] || !Lhit[i])
            return false;

    double minB = B[0], maxB = B[0], minL = L[0], maxL = L[0];
    for (int i = 1; i < 4; ++i)
    {
        if (B[i] < minB)
            minB = B[i];
        if (B[i] > maxB)
            maxB = B[i];
        if (L[i] < minL)
            minL = L[i];
        if (L[i] > maxL)
            maxL = L[i];
    }
    out.minLat = minB;
    out.maxLat = maxB;
    out.minLon = minL;
    out.maxLon = maxL;
    return true;
}

static inline bool point_in_bounds(double lat, double lon, const RegionBounds &R)
{
    return (lat >= R.minLat && lat <= R.maxLat && lon >= R.minLon && lon <= R.maxLon);
}

// --- 直接判断：目标经纬度是否在“最新区域”内 ---
bool is_target_in_latest_region(const std::string &dir,
                                double targetLat, double targetLon)
{
    std::string latest;
    if (!find_latest_txt(dir, latest))
        return false; // 读取失败按“不在”处理

    RegionBounds rb;
    if (!parse_region_bounds_from_file(latest, rb))
        return false;
    DBG("边界"<<" minLat="<<rb.minLat<<" maxLat="<<rb.maxLat
        <<" minLon="<<rb.minLon<<" maxLon="<<rb.maxLon);
    return point_in_bounds(targetLat, targetLon, rb);
}

// 当天 0:00 起的本地秒数 + 当前 UTC 秒
static inline void get_day_seconds_local_and_utc(double &sec_of_day_local, double &utc_seconds)
{
    using clock = std::chrono::system_clock;
    const auto now_tp = clock::now();
    const std::time_t now_tt = clock::to_time_t(now_tp);
    utc_seconds = static_cast<double>(now_tt);

    std::tm local_tm{};
#if defined(_WIN32)
    localtime_s(&local_tm, &now_tt);
#else
    localtime_r(&now_tt, &local_tm);
#endif
    sec_of_day_local = local_tm.tm_hour * 3600.0 + local_tm.tm_min * 60.0 + local_tm.tm_sec;
}

// ========== 主函数：给定速度、给定 L0 的往返运动 ==========
void pesudoTarget_Gen(GMTIOutput &resRef,
                      double lat_st, double lon_st,
                      double lat_ed, double lon_ed,
                      double v_mps, double L0_deg,
                      const Config &cfg)
{
    // 1) 起止点投影到高斯平面（x=Northing, y=Easting）
    double Ns, Es, Ne, Ee; // N=Northing, E=Easting
    Gaussp3(lat_st, lon_st, L0_deg, Ns, Es);
    Gaussp3(lat_ed, lon_ed, L0_deg, Ne, Ee);

    // 2) 线段长度（米）与当前时间
    const double dE = Ee - Es;
    const double dN = Ne - Ns;
    const double D = std::hypot(dE, dN);

    double t_local = 0.0, t_utc = 0.0;
    get_day_seconds_local_and_utc(t_local, t_utc);

    if (D <= 1e-6)
    {
        // 退化：起止重合
        resRef.MT.clear();
        resRef.MT.push_back(lat_st);
        resRef.MT.push_back(lon_st);
        resRef.MT.push_back(cfg.MT_nowz);
        resRef.MT.push_back(0.0); // xP
        resRef.MT.push_back(0.0); // yP
        resRef.MT.push_back(t_utc);
        return;
    }

    // 3) 等速往返：路程 s_m = |v| * t_local，在 [0, 2D) 上取模
    const double v = std::fabs(v_mps);
    const double twoD = 2.0 * D;
    double s_m = std::fmod(v * t_local, twoD);
    if (s_m < 0.0)
        s_m += twoD;

    // 4) 归一化位置 frac ∈ [0,1]（0 起点，1 终点，超出后镜像回退）
    const double frac = (s_m <= D) ? (s_m / D) : ((twoD - s_m) / D);

    // 5) 当前平面点（E,N）
    const double Enow = Es + frac * dE;
    const double Nnow = Ns + frac * dN;

    // 6) 反算回经纬度；并计算相对起点位移 xP/yP（E/N 增量）
    double lat_now = 0.0, lon_now = 0.0;
    if (!Gaussp3RV(/*x=N*/ Nnow, /*y=E*/ Enow, L0_deg, lat_now, lon_now))
        return;

    const double xP = Nnow; // Easting 增量
    const double yP = Enow; // Northing 增量

    // 7) 输出（若要保留历史轨迹，请删掉 clear）
    if (is_target_in_latest_region("/root/桌面/GMTI/result/", lat_now, lon_now))
    {
        // 在区域内
        DBG("伪目标在有效区域内：lat="<<lat_now<<" lon="<< lon_now);
        resRef.MT.push_back(lat_now);
        resRef.MT.push_back(lon_now);
        resRef.MT.push_back(cfg.MT_nowz);
        resRef.MT.push_back(xP);
        resRef.MT.push_back(yP);
        resRef.MT.push_back(t_utc);
    }
}
