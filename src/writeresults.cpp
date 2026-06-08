#include "GMTIProcessor.hpp"
#include "geo/geoProj.hpp" // 包含 Gaussp3RV

#include <vector>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <string>
#include <iostream>
#include <cerrno>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h> // opendir, readdir, closedir
#include <cstdio>   // std::snprintf

static inline bool starts_with(const char *s, const char *pfx)
{
    while (*pfx)
    {
        if (*s++ != *pfx++)
            return false;
    }
    return true;
}
static inline bool ends_with(const char *s, const char *sfx)
{
    const size_t ls = std::strlen(s), lr = std::strlen(sfx);
    return (ls >= lr) && (std::strcmp(s + (ls - lr), sfx) == 0);
}

static inline double wrap180_deg(double angle_deg)
{
    angle_deg = std::fmod(angle_deg + 180.0, 360.0);
    if (angle_deg < 0.0)
        angle_deg += 360.0;
    return angle_deg - 180.0;
}

bool GMTIProcessor::nextGMTIFileName(const std::string &dir,
                                     std::string &out_path,
                                     int &out_index) const
{
    out_path.clear();
    out_index = 0;

    // 1) 确保目录存在
    std::string d = dir;
    if (!d.empty() && d.back() != '/' && d.back() != '\\')
        d.push_back('/');
    if (!mkdir_p(d))
    {
        std::cerr << "nextGMTIFileName: 创建目录失败: " << d << "\n";
        return false;
    }

    // 2) 打开目录并扫描现有文件，匹配 "GMTIxx.bin"
    int max_idx = 0;

    DIR *dp = ::opendir(d.c_str());
    if (dp)
    {
        while (dirent *ent = ::readdir(dp))
        {
            const char *name = ent->d_name;
            // 过滤掉 "." ".."
            if (name[0] == '.')
                continue;

            // 格式严格：长度==10，前缀GMTI，后缀.bin，中间两位数字
            // "GMTI" + 2 + ".bin" = 4 + 2 + 4 = 10
            const size_t len = std::strlen(name);
            if (len != 10)
                continue;
            if (!starts_with(name, "GMTI"))
                continue;
            if (!ends_with(name, ".bin"))
                continue;
            if (name[4] < '0' || name[4] > '9')
                continue;
            if (name[5] < '0' || name[5] > '9')
                continue;

            int idx = (name[4] - '0') * 10 + (name[5] - '0');
            if (idx >= 1 && idx <= 99)
            {
                if (idx > max_idx)
                    max_idx = idx;
            }
        }
        ::closedir(dp);
    }
    else
    {
        // 目录不可读也不致命（可能刚创建），从 01 开始
    }

    // 3) 下一个序号
    const int next_idx = max_idx + 1;
    if (next_idx > 99)
    {
        std::cerr << "nextGMTIFileName: 已达到 99 个文件，无法继续编号。\n";
        return false;
    }

    // 4) 组装路径
    char fname[16];
    std::snprintf(fname, sizeof(fname), "GMTI%02d.bin", next_idx);
    out_path = d + fname;
    out_index = next_idx;
    return true;
}

static inline void put_u8(std::vector<uint8_t> &buf, size_t pos, uint8_t v)
{
    buf[pos] = v;
}

static inline void put_u16_le(std::vector<uint8_t> &buf, size_t pos, uint16_t v)
{
    buf[pos + 0] = static_cast<uint8_t>(v & 0xFF);
    buf[pos + 1] = static_cast<uint8_t>((v >> 8) & 0xFF);
}

static inline void put_f32_le(std::vector<uint8_t> &buf, size_t pos, float v)
{
    static_assert(sizeof(float) == 4, "float must be 4 bytes");
    uint32_t u;
    std::memcpy(&u, &v, sizeof(u));
    // 写为 little-endian
    buf[pos + 0] = static_cast<uint8_t>(u & 0xFF);
    buf[pos + 1] = static_cast<uint8_t>((u >> 8) & 0xFF);
    buf[pos + 2] = static_cast<uint8_t>((u >> 16) & 0xFF);
    buf[pos + 3] = static_cast<uint8_t>((u >> 24) & 0xFF);
}

static inline void put_f64_le(std::vector<uint8_t> &buf, size_t pos, double v)
{
    static_assert(sizeof(double) == sizeof(uint64_t), "double must be 8 bytes");
    uint64_t u = 0;
    std::memcpy(&u, &v, sizeof(u));
    for (int i = 0; i < 8; ++i)
    {
        buf[pos + static_cast<size_t>(i)] = static_cast<uint8_t>((u >> (8 * i)) & 0xFF);
    }
}

static inline void put_i32_le(std::vector<uint8_t> &buf, size_t pos, int32_t v);
static inline int32_t quant_deg_to_i32(double deg);

static inline uint16_t quant_speed_to_u16(double speed)
{
    const double q = std::round(std::max(0.0, speed) / 0.01);
    if (q > 65535.0)
        return 65535;
    return static_cast<uint16_t>(q);
}

bool GMTIProcessor::writeResult(const std::vector<double> &res, const Config &cfg)
{
    // 布局：2 + N*36
    // 单记录：id(u16) + lon(i32) + lat(i32) + speed(u16) + direction(f64) + range(f64) + utc(f64)
    constexpr size_t REC_BYTES = 36;
    constexpr size_t HEADER_BYTES = 2;
    const size_t n_targets = res.size() / 8;
    const size_t FILE_BYTES = HEADER_BYTES + n_targets * REC_BYTES;

    std::vector<uint8_t> buf(FILE_BYTES, 0u);

    // [0..1] 写总数（u16, LE）
    put_u16_le(buf, 0, static_cast<uint16_t>(n_targets));

    // 逐目标写 36 字节记录（含每目标 utc）
    for (size_t i = 0; i < n_targets; ++i)
    {
        const size_t base = HEADER_BYTES + i * REC_BYTES;

        const double lat = res[i * 8 + 0];
        const double lng = res[i * 8 + 1];
        const double target_utc = res[i * 8 + 5];
        const double direction = wrap180_deg(res[i * 8 + 6]);
        const double range = res[i * 8 + 7];

        // 0..1: id (u16, LE)
        put_u16_le(buf, base + 0, static_cast<uint16_t>(std::max<size_t>(1, std::min<size_t>(65535, i + 1))));

        // 2..5: 经度（int32）
        put_i32_le(buf, base + 2, quant_deg_to_i32(lng));
        // 6..9: 纬度（int32）
        put_i32_le(buf, base + 6, quant_deg_to_i32(lat));
        // 10..11: 速度（u16），当前周期检测结果暂写 0
        put_u16_le(buf, base + 10, 0);
        // 12..19: 方向（double）
        put_f64_le(buf, base + 12, direction);
        // 20..27: 距离（double）
        put_f64_le(buf, base + 20, range);
        // 28..35: 目标 utc（double）
        put_f64_le(buf, base + 28, target_utc);
    }

    // 目标路径：优先使用当前回波文件编号，避免扫描目录自动加一导致关联窗口错位。
    std::string outpath;
    int idx = 0;
    if (cfg.result_file_id > 0 && cfg.result_file_id <= 99) {
        std::string d = cfg.result_add;
        if (!d.empty() && d.back() != '/' && d.back() != '\\') {
            d.push_back('/');
        }
        if (!mkdir_p(d)) {
            std::cerr << "writeResult: 创建目录失败: " << d << "\n";
            return false;
        }
        char fname[16];
        std::snprintf(fname, sizeof(fname), "GMTI%02d.bin", cfg.result_file_id);
        outpath = d + fname;
        idx = cfg.result_file_id;
    } else if (!nextGMTIFileName(cfg.result_add, outpath, idx)) {
        return false;
    }

    std::ofstream ofs(outpath.c_str(), std::ios::binary);
    if (!ofs)
    {
        std::cerr << "writeResult: 无法打开输出文件: " << outpath << "\n";
        return false;
    }
    ofs.write(reinterpret_cast<const char *>(buf.data()),
              static_cast<std::streamsize>(buf.size()));
    if (!ofs)
    {
        std::cerr << "writeResult: 写文件失败: " << outpath << "\n";
        return false;
    }
    ofs.close();
    std::cout << "[GMTI] writeResult: 写入当前检测结果 "
              << outpath << "，目标数: " << n_targets << std::endl;
    return true;
}

// ---------- 小工具：递归创建目录 ----------
bool GMTIProcessor::mkdir_p(const std::string &dir)
{
    if (dir.empty())
        return true;
    std::string cur;
    size_t i = 0;
    if (dir[0] == '/')
    {
        cur = "/";
        i = 1;
    }
    for (; i < dir.size(); ++i)
    {
        cur.push_back(dir[i]);
        if (dir[i] == '/' || i == dir.size() - 1)
        {
            if (cur == "/" || cur == "./" || cur == ".")
                continue;
            if (::mkdir(cur.c_str(), 0755) != 0)
            {
                if (errno == EEXIST)
                    continue;
                if (errno == ENOENT)
                    continue; // 父目录尚未就绪，继续循环
                std::cerr << "mkdir_p failed: " << cur << " errno=" << errno << "\n";
                return false;
            }
        }
    }
    return true;
}

// ---------- 占位：xy(m) -> (lat,lng) 度 ----------

// ---------- 写入帮助：LE 编码 ----------
static inline void put_i32_le(std::vector<uint8_t> &buf, size_t pos, int32_t v)
{
    uint32_t u = static_cast<uint32_t>(v);
    buf[pos + 0] = static_cast<uint8_t>(u & 0xFF);
    buf[pos + 1] = static_cast<uint8_t>((u >> 8) & 0xFF);
    buf[pos + 2] = static_cast<uint8_t>((u >> 16) & 0xFF);
    buf[pos + 3] = static_cast<uint8_t>((u >> 24) & 0xFF);
}

// 量化：度 -> int32（LSB=8.38191e-8 度），四舍五入并饱和
static inline int32_t quant_deg_to_i32(double deg)
{
    const double LSB = 8.38191e-8; // deg / LSB
    double q = deg / LSB;
    // 四舍五入
    if (q >= 0.0)
        q = std::floor(q + 0.5);
    else
        q = std::ceil(q - 0.5);
    // 饱和到 int32
    if (q > 2147483647.0)
        q = 2147483647.0;
    if (q < -2147483648.0)
        q = -2147483648.0;
    return static_cast<int32_t>(q);
}

// ---------- 核心：写入所有 tracks 的所有 pos ----------
bool GMTIProcessor::writeTracksBinary(const std::vector<Track> &tracks,
                                      double utcMid,
                                      const GMTIOutput::Plane &plane,
                                      const Config &cfg)
{
    (void)utcMid;
    (void)plane;

    // 逐目标写 28 字节记录，文件尾追加 8 字节 utc
    const size_t REC_BYTES = 28;
    const size_t HEADER_BYTES = 2;
    std::vector<uint8_t> buf;

    // 1) 将所有轨迹的滤波点 kf 展开为点序列（按 track.id 升序；每条内按时间顺序）
    struct Rec
    {
        uint16_t id;
        double x;
        double y;
        double speed;
        double direction;
        double range;
    };
    std::vector<Rec> points;
    points.reserve(64);

    // 拷贝并排序轨迹（确保 id 升序；id<=0 的放后面）
    std::vector<const Track *> order;
    order.reserve(tracks.size());
    for (size_t i = 0; i < tracks.size(); ++i)
        order.push_back(&tracks[i]);
    std::stable_sort(order.begin(), order.end(),
                     [](const Track *a, const Track *b)
                     {
                         const int ia = (a->id > 0) ? a->id : 0x7fffffff;
                         const int ib = (b->id > 0) ? b->id : 0x7fffffff;
                         if (ia != ib)
                             return ia < ib;
                         return a < b;
                     });

    for (size_t ti = 0; ti < order.size(); ++ti)
    {
        const Track &tr = *order[ti];
        const uint16_t tid = (tr.id > 0 && tr.id <= 65535) ? static_cast<uint16_t>(tr.id)
                                                           : static_cast<uint16_t>(std::min<size_t>(65535, ti + 1));
        // 只写入每条轨迹的第一个滤波点
        if (tr.kf.empty())
            continue;
        double speed = 0.0;
        if (tr.pos.size() >= 2 && tr.time.size() >= 2)
        {
            const size_t n = tr.pos.size();
            const double dx = tr.pos[n - 1][0] - tr.pos[n - 2][0];
            const double dy = tr.pos[n - 1][1] - tr.pos[n - 2][1];
            const double dt = tr.time[n - 1] - tr.time[n - 2];
            speed = std::sqrt(dx * dx + dy * dy) / (dt > 1e-6 ? dt : 1.0);
        }
        Rec r;
        r.id = tid;
        r.x = tr.kf.back()[0];
        r.y = tr.kf.back()[1];
        r.speed = speed;
        r.direction = tr.direction;
        r.range = tr.range;
        points.push_back(r);
    }

    std::cout << "[GMTI] writeTracksBinary: 共 " << tracks.size() << " 个目标，"
              << points.size() << " 个点。\n";

    // 2 字节目标数 + N*28 字节记录 + 8 字节 utc
    buf.reserve(HEADER_BYTES + points.size() * REC_BYTES + 8);
    buf.resize(HEADER_BYTES, 0u);
    put_u16_le(buf, 0, static_cast<uint16_t>(std::min<size_t>(65535, points.size())));

    // 逐点写 28 字节记录
    for (size_t i = 0; i < points.size(); ++i)
    {
        const size_t base = HEADER_BYTES + i * REC_BYTES;
        buf.resize(base + REC_BYTES, 0u);

        // 2 byte id
        put_u16_le(buf, base + 0, points[i].id);

        // 坐标转换：x,y(m) → lat,lng(deg)
        double lat = 0.0, lng = 0.0;
        Gaussp3RV(points[i].x, points[i].y, cfg.L0, lat, lng);

        // 2..5 经度，6..9 纬度
        put_i32_le(buf, base + 2, quant_deg_to_i32(lng));
        put_i32_le(buf, base + 6, quant_deg_to_i32(lat));
        // 10..11 速度：轨迹末端速度估计
        put_u16_le(buf, base + 10, quant_speed_to_u16(points[i].speed));
        // 12..19 方向，20..27 距离
        put_f64_le(buf, base + 12, wrap180_deg(points[i].direction));
        put_f64_le(buf, base + 20, points[i].range);
    }

    put_f64_le(buf, HEADER_BYTES + points.size() * REC_BYTES, utcMid);

    // 输出到 result_add 目录
    std::string outpath;
    int idx = 0;
    if (!nextGMTIFileName(cfg.result_add, outpath, idx))
    {
        return false; // 或者回退到固定文件名
    }

    std::ofstream ofs(outpath.c_str(), std::ios::binary);
    if (!ofs)
    {
        std::cerr << "writeTracksBinary: 无法打开输出文件: " << outpath << "\n";
        return false;
    }
    ofs.write(reinterpret_cast<const char *>(buf.data()),
              static_cast<std::streamsize>(buf.size()));
    if (!ofs)
    {
        std::cerr << "writeTracksBinary: 写文件失败: " << outpath << "\n";
        return false;
    }
    std::cout << "writeTracksBinary: 写文件成功: " << outpath << "\n";
    return true;
}
