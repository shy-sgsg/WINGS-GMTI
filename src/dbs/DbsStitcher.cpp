#include "dbs/DbsStitcher.hpp"
#include "dbs/DbsIO.hpp"
#include <fstream>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include "dbs/lodepng.h" // PNG写出
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <vector>
#include <cstdint>
#include <iostream>
#include <limits>
#include <cerrno>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h> // opendir, readdir, closedir
#include "dbs/rotation_xy.hpp"

static inline float c0() { return 299792458.0f; }

static size_t dbs_max_mosaic_pixels()
{
  const char *env = std::getenv("DBS_MAX_MOSAIC_PIXELS");
  if (env && *env)
  {
    char *end = nullptr;
    const unsigned long long v = std::strtoull(env, &end, 10);
    if (end != env && v > 0)
      return static_cast<size_t>(v);
  }
  return 60000000ULL;
}

DbsStitcher::DbsStitcher()
{
}
DbsStitcher::~DbsStitcher()
{
}

// ---------- 小工具：递归创建目录 ----------
bool DbsStitcher::mkdir_p(const std::string &dir)
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

bool DbsStitcher::nextGMTIFileNamePNG(const std::string &dir,
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
      if (!ends_with(name, ".png"))
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
  std::snprintf(fname, sizeof(fname), "GMTI%02d.png", next_idx);
  out_path = d + fname;
  out_index = next_idx;
  return true;
}

bool DbsStitcher::nextGMTIFileNameTxt(const std::string &dir,
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
      if (!ends_with(name, ".txt"))
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
  std::snprintf(fname, sizeof(fname), "GMTI%02d.txt", next_idx);
  out_path = d + fname;
  out_index = next_idx;
  return true;
}

bool DbsStitcher::loadPosAndEstimateVelocity(const Params &P, PosData &POS)
{
  // 读取 POS: 每帧 7 个 double
  std::ifstream ifs(P.pos_path, std::ios::binary);
  if (!ifs)
    return false;
  const size_t elems = 7;
  std::vector<double> buf(elems);
  while (ifs.read(reinterpret_cast<char *>(buf.data()), elems * sizeof(double)))
  {
    POS.gpstime.push_back(buf[0]); // 时间戳
    POS.lat_rad.push_back(buf[1]); // 纬度 (弧度)
    POS.lon_rad.push_back(buf[2]); // 经度 (弧度)
    POS.alt_m.push_back(buf[3]);   // 高度 (米)
  }
  size_t N = POS.gpstime.size();
  POS.vx.resize(N);
  POS.vy.resize(N);
  POS.vz.resize(N);
  POS.speed.resize(N);
  POS.x_m.resize(N);
  POS.y_m.resize(N);

  for (size_t i = 0; i < N; i++)
  {
    double lat_deg = POS.lat_rad[i] * 180.0 / M_PI;
    double lon_deg = POS.lon_rad[i] * 180.0 / M_PI;
    proj_ll_to_xy(lat_deg, lon_deg, POS.x_m[i], POS.y_m[i]);
  }
  // 中心差分速度
  auto diff = [&](const std::vector<double> &v, std::vector<double> &dv)
  {
    size_t N = v.size();
    dv.resize(N);
    for (size_t i = 0; i < N; i++)
    {
      size_t i0 = (i > 4) ? i - 4 : i;
      size_t i1 = (i + 4 < N) ? i + 4 : i;
      double dt = POS.gpstime[i1] - POS.gpstime[i0];
      dv[i] = (dt != 0.0) ? (v[i1] - v[i0]) / dt : 0.0;
    }
  };
  diff(POS.x_m, POS.vx);
  diff(POS.y_m, POS.vy);
  diff(POS.alt_m, POS.vz);
  for (size_t i = 0; i < N; i++)
    POS.speed[i] = std::sqrt(POS.vx[i] * POS.vx[i] + POS.vy[i] * POS.vy[i] + POS.vz[i] * POS.vz[i]);

  return true;
}

  int DbsStitcher::estimateEffectiveAzBins(const Params &P, const PosData &POS)
  {
    // 与 MATLAB 逻辑一致的估算（略简化）
    double Rbin = c0() / (2.0 * P.fs_hz);
    double air_h = (std::accumulate(POS.alt_m.begin(), POS.alt_m.end(), 0.0) / POS.alt_m.size()) - P.mean_ground_h;
    double rmax = P.Rmin_m + P.range_samp_used * Rbin;
    double max_cosphi = std::sqrt(std::max(0.0, 1.0 - (air_h / rmax) * (air_h / rmax)));

    int azN = std::max(1, int(std::floor(P.scan_max_az_deg)));
    double max_delta_sin = 0;
    for (int i = 1; i <= azN; i++)
    {
      double d = std::abs(std::sin((i + (P.beamwidth_deg + 1) / 2) * M_PI / 180.0) - std::sin((i - (P.beamwidth_deg + 1) / 2) * M_PI / 180.0));
      max_delta_sin = std::max(max_delta_sin, d);
    }
    double lam = c0() / P.fc_hz;
    double v_med = 0.0;
    if (POS.speed.size() > 100) // 确保有足够的元素
    {
      // 从索引100到末尾
      auto max_iter = std::max_element(POS.speed.begin() + 100, POS.speed.end());
      if (max_iter != POS.speed.end())
      {
        v_med = *max_iter;
      }
    }

    double max_dopp = 2 * v_med / lam * max_delta_sin * max_cosphi;
    int nEff = int(std::ceil(max_dopp / (1.0 * P.PRF / P.pulses_per_beam)));
    nEff = int(std::floor(nEff * 1.4));
    if (nEff % 2)
      nEff++;

    return nEff;
  }

bool DbsStitcher::estimateMosaicExtent(const Config &cfg, const RDData &RD, const MetaPack &meta,
                                       Bounds &b, Grid &grid)
{
  // === 用你“精简版静止目标定位”逻辑：只更新 min/max ===
  const int B = (int)meta.beams.size(); // bowei_num
  const int M = cfg.rg_len;
  const int nEff = RD.nEff;
  if (B <= 0 || M <= 0 || nEff <= 0 || cfg.fc <= 0.0 || cfg.dbs_out_res_m <= 0.0 ||
      RD.fd_axis.size() < static_cast<size_t>(B) || RD.rg_axis.size() < static_cast<size_t>(B))
  {
    std::cerr << "[fusion][dbs][extent] invalid inputs: B=" << B
              << " M=" << M << " nEff=" << nEff
              << " fc=" << cfg.fc << " out_res=" << cfg.dbs_out_res_m << std::endl;
    return false;
  }

  float lam = float(c0() / cfg.fc);
  float minX = std::numeric_limits<float>::infinity();
  float maxX = -std::numeric_limits<float>::infinity();
  float minY = std::numeric_limits<float>::infinity();
  float maxY = -std::numeric_limits<float>::infinity();

  int step = std::max(1, cfg.dbs_range_skip);

  std::vector<float> vN(B), vE(B), V_feiji(B), sin_jiaodu(B), cos_jiaodu(B);
  for (int b = 0; b < B; ++b)
  {
    vN[b] = meta.beams[b].vN;
    vE[b] = meta.beams[b].vE;
    V_feiji[b] = std::sqrt(vN[b] * vN[b] + vE[b] * vE[b]);
    const float denom = (std::fabs(vE[b]) > 1e-12f) ? std::fabs(vE[b]) : 1e-12f;
    const float jiaodu = std::atan(std::fabs(vN[b] / denom));
    sin_jiaodu[b] = std::sin(jiaodu);
    cos_jiaodu[b] = std::cos(jiaodu);
  }

  for (int n = 0; n < B; n++)
  {
    if ((n % std::max(1, cfg.dbs_beam_skip)) != 0)
      continue;

    const auto &fdv = RD.fd_axis[n]; // 1 x nEff
    const auto &rgv = RD.rg_axis[n]; // 1 x M
    if (fdv.size() < static_cast<size_t>(nEff) || rgv.size() < static_cast<size_t>(M))
    {
      std::cerr << "[fusion][dbs][extent] skip beam " << n
                << " due to axis size mismatch: fd=" << fdv.size()
                << " rg=" << rgv.size() << std::endl;
      continue;
    }

    float V_f = V_feiji[n];
    if (!std::isfinite(V_f) || V_f <= 1e-3f)
    {
      std::cerr << "[fusion][dbs][extent] skip beam " << n
                << " due to invalid speed V=" << V_f << std::endl;
      continue;
    }
    float Height = meta.beams[n].z - float(cfg.MT_nowz);
    float start_x = meta.beams[n].x, start_y = meta.beams[n].y;

    int flag = infer_flight_dir_flag(vN, vE); // 1/2/3/4
    int squint_side = cfg.squint_side;           // 0 for left, 1 for right
    flag = flag + squint_side * 4;


    for (int i = 0; i < M; i += step)
    {
      float R = rgv[i];
      for (int jj = 0; jj < nEff; jj++)
      {
        float fd = fdv[jj];
        float x0 = R * lam * fd / (2 * V_f);
        float t2 = R * R - x0 * x0 - Height * Height;
        if (t2 <= 0)
          continue;
        float y0 = std::sqrt(t2);

        float tx, ty;
        rotation_xy(x0, y0,
                    flag,
                    cos_jiaodu[n], sin_jiaodu[n],
                    tx, ty);

        float gx = tx + start_x;
        float gy = ty + start_y;
        if (!std::isfinite(gx) || !std::isfinite(gy))
          continue;
        minX = std::min(minX, gx);
        maxX = std::max(maxX, gx);
        minY = std::min(minY, gy);
        maxY = std::max(maxY, gy);
      }
    }
  }

  b.minX = minX;
  b.maxX = maxX;
  b.minY = minY;
  b.maxY = maxY;
  if (!std::isfinite(minX) || !std::isfinite(maxX) ||
      !std::isfinite(minY) || !std::isfinite(maxY) ||
      !(maxX > minX) || !(maxY > minY))
  {
    std::cerr << "[fusion][dbs][extent] invalid bounds: X=[" << minX << ", " << maxX
              << "] Y=[" << minY << ", " << maxY << "]" << std::endl;
    return false;
  }

  // 生成规则网格
  double dx = std::max(1.0, cfg.dbs_out_res_m);
  const double spanX = static_cast<double>(maxX) - static_cast<double>(minX);
  const double spanY = static_cast<double>(maxY) - static_cast<double>(minY);
  auto axis_count = [](double span, double res) -> double {
    return std::floor(span / res) + 1.0;
  };

  double nxD = axis_count(spanX, dx);
  double nyD = axis_count(spanY, dx);
  double pixD = nxD * nyD;
  const size_t maxPixels = dbs_max_mosaic_pixels();
  if (!std::isfinite(nxD) || !std::isfinite(nyD) || !std::isfinite(pixD) ||
      nxD <= 0.0 || nyD <= 0.0)
  {
    std::cerr << "[fusion][dbs][extent] invalid grid estimate: spanX=" << spanX
              << " spanY=" << spanY << " dx=" << dx << std::endl;
    return false;
  }

  if (pixD > static_cast<double>(maxPixels))
  {
    const double oldDx = dx;
    const double oldPixD = pixD;
    dx *= std::sqrt(pixD / static_cast<double>(maxPixels));
    nxD = axis_count(spanX, dx);
    nyD = axis_count(spanY, dx);
    pixD = nxD * nyD;
    std::cerr << "[fusion][dbs][extent][warn] mosaic grid too large at requested resolution: "
              << "spanX=" << spanX << "m spanY=" << spanY << "m out_res=" << oldDx
              << "m pixels~" << static_cast<unsigned long long>(oldPixD)
              << ". Auto coarsen output resolution to " << dx
              << "m. Set DBS_MAX_MOSAIC_PIXELS to override the cap (current "
              << maxPixels << ")." << std::endl;
  }

  if (pixD > static_cast<double>(maxPixels) * 1.02)
  {
    std::cerr << "[fusion][dbs][extent] grid still exceeds cap after coarsening: pixels~"
              << pixD << " cap=" << maxPixels << std::endl;
    return false;
  }

  const size_t nx = static_cast<size_t>(nxD);
  const size_t ny = static_cast<size_t>(nyD);
  grid.x.assign(nx, 0.0f);
  grid.y.assign(ny, 0.0f);
  for (size_t i = 0; i < nx; ++i)
    grid.x[i] = static_cast<float>(static_cast<double>(minX) + static_cast<double>(i) * dx);
  for (size_t i = 0; i < ny; ++i)
    grid.y[i] = static_cast<float>(static_cast<double>(minY) + static_cast<double>(i) * dx);
  grid.dx = static_cast<float>(dx);
  grid.dy = static_cast<float>(dx);

  const double hostMosaicGiB = static_cast<double>(nx) * static_cast<double>(ny) *
                               static_cast<double>(sizeof(float) + sizeof(uint16_t)) /
                               (1024.0 * 1024.0 * 1024.0);
  std::cout << "[fusion][dbs][extent] B=" << B << " M=" << M << " nEff=" << nEff
            << " boundsX=[" << minX << ", " << maxX << "] boundsY=[" << minY << ", " << maxY
            << "] grid=" << nx << "x" << ny << " dx=" << grid.dx
            << " host_mosaic_min~" << hostMosaicGiB << " GiB" << std::endl;
  return true;
}

// ---- 小工具：转置、左右翻转、均值、PNG写出(8bit灰度) ----
template <typename T>
static Image2D<T> transpose_copy(const Image2D<T> &src)
{
  Image2D<T> dst(src.cols, src.rows);
  for (int r = 0; r < src.rows; ++r)
    for (int c = 0; c < src.cols; ++c)
      dst.at(c, r) = src.at(r, c);
  return dst;
}
template <typename T>
static void fliplr_inplace(Image2D<T> &img)
{
  for (int r = 0; r < img.rows; ++r)
  {
    for (int c = 0; c < img.cols / 2; ++c)
    {
      std::swap(img.at(r, c), img.at(r, img.cols - 1 - c));
    }
  }
}
static double mean_of_image(const Image2D<float> &img)
{
  if (img.empty())
    return 0.0;
  double s = 0.0;
  for (size_t i = 0; i < img.buf.size(); ++i)
    s += img.buf[i];
  return s / (double)img.buf.size();
}
// 使用 lodepng 写 8bit 灰度 PNG（data 按行主）
static bool write_png_gray8(const std::string &path,
                            const unsigned char *data,
                            unsigned width, unsigned height)
{
  unsigned err = lodepng_encode_file(path.c_str(), data, width, height, LCT_GREY, 8);
  return err == 0;
}

bool DbsStitcher::writeProducts(const Config &cfg,
                                const Grid &grid,
                                const Mosaic &mosaic,
                                const Bounds &b)
{
  (void)grid;
  // 1) 阈值截断 + 归一化到 16bit
  Image2D<float> amp = mosaic.amp; // 拷贝一份做处理
  const double th_src = mean_of_image(amp) * 13.0;
  const double th = (th_src > 1e-12) ? th_src : 1.0; // 防止除零

  for (size_t i = 0; i < amp.buf.size(); ++i)
  {
    double v = amp.buf[i];
    if (v > th)
      v = th;
    amp.buf[i] = (float)(v * (65500.0 / th));
  }

  // 转为 uint16_t
  Image2D<uint16_t> u16(amp.rows, amp.cols);
  for (size_t i = 0; i < amp.buf.size(); ++i)
  {
    double v = amp.buf[i];
    if (v < 0.0)
      v = 0.0;
    if (v > 65535.0)
      v = 65535.0;
    u16.buf[i] = (uint16_t)(v + 0.5); // 四舍五入
  }

  // 2) RAW/PNG：先 转置，再 左右翻转（与原 MATLAB 一致）
  Image2D<uint16_t> out16 = transpose_copy(u16);
  fliplr_inplace(out16);

#ifdef DEBUG
  Image2D<float> temp_src = mosaic.r_save;
  Image2D<float> temp_dst = transpose_copy(temp_src);
  fliplr_inplace(temp_dst);

  std::ofstream src2dst_cfs0("DBS_range.dat", std::ios::binary);
  if (!src2dst_cfs0)
    {
        std::cerr << "writeResult: 无法打开输出文件: " << "DBS_range.dat" << "\n";
        return false;
    }
  src2dst_cfs0.write(reinterpret_cast<const char *>(temp_dst.buf.data()),
            static_cast<std::streamsize>(temp_dst.buf.size() * 4));

  temp_src = mosaic.fd_save;
  temp_dst = transpose_copy(temp_src);
  fliplr_inplace(temp_dst);

  std::ofstream src2dst_cfs1("DBS_doppler.dat", std::ios::binary);
  if (!src2dst_cfs1)
    {
        std::cerr << "writeResult: 无法打开输出文件: " << "DBS_doppler.dat" << "\n";
        return false;
    }
  src2dst_cfs1.write(reinterpret_cast<const char *>(temp_dst.buf.data()),
            static_cast<std::streamsize>(temp_dst.buf.size() * 4));

#endif

  // 3) 写 PNG（8bit）：直接用 RAW 同一张（转置+翻转后的）数据 /256
  std::vector<unsigned char> u8((size_t)out16.rows * out16.cols);
  for (size_t i = 0; i < u8.size(); ++i)
    u8[i] = (unsigned char)(out16.buf[i] >> 8); // /256

  std::string outpath;
  int idx = 0;
  if (!nextGMTIFileNamePNG(cfg.result_add, outpath, idx))
  {
    return false; // 或者回退到固定文件名
  }

  if (!write_png_gray8(outpath, &u8[0], (unsigned)out16.cols, (unsigned)out16.rows))
    return false;

  // 4) 四角经纬度 .txt
  double B0, L0, B1, L1, B2, L2, B3, L3;
  proj_xy_to_ll(b.minX, b.minY, B0, L0);
  proj_xy_to_ll(b.maxX, b.maxY, B1, L1);
  proj_xy_to_ll(b.maxX, b.minY, B2, L2);
  proj_xy_to_ll(b.minX, b.maxY, B3, L3);

  outpath.clear();
  idx = 0;
  if (!nextGMTIFileNameTxt(cfg.result_add, outpath, idx))
  {
    return false; // 或者回退到固定文件名
  }

  std::ofstream cfs(outpath.c_str(), std::ios::binary);
  if (!cfs)
    return false;
  cfs.setf(std::ios::fixed);
  cfs.precision(6);
  cfs << "B0 = " << B0 << "\n";
  cfs << "B1 = " << B1 << "\n";
  cfs << "L0 = " << L0 << "\n";
  cfs << "L1 = " << L1 << "\n";
  cfs << "B2 = " << B2 << "\n";
  cfs << "B3 = " << B3 << "\n";
  cfs << "L2 = " << L2 << "\n";
  cfs << "L3 = " << L3 << "\n";
  cfs.close();

  return true;
}
