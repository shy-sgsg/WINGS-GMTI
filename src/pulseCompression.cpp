#include <vector>
#include <complex>
#include <cmath>
#include <cstring>
#include <mutex>
#include "rangeCompress.hpp"
#include <fftw3.h>
#include <cstdint>

static inline int16_t load_int16_le(const uint8_t* p) {
  return static_cast<int16_t>(p[0] | (p[1] << 8));
}
static inline uint32_t load_u32_le(const uint8_t* p) {
  return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

inline std::vector<std::complex<double>>
buildHfSpec(int N, double fs, double Kr)
{
    std::vector<std::complex<double>> Hf(N);
    if (N <= 0 || fs <= 0.0 || Kr == 0.0) return Hf;

    const double df = fs / double(N);
    // 直接用“未移位的 DFT 频率”公式：
    //   n = 0..N-1
    //   f(n) = (n <= N/2-1) ? n*df : (n-N)*df
    for (int n = 0; n < N; ++n) {
        double fn = (n <= (N/2 - 1)) ? n*df : (n - N)*df;
        double phase = M_PI * (fn*fn) / Kr;  // -pi*f^2/Kr
        Hf[n] = std::complex<double>(std::cos(phase), std::sin(phase));
    }
    return Hf;
}

inline std::vector<std::complex<double>>
buildHfSpecFromRefFunc(const Config &P){
  int N = P.pulse_len;
  FILE* fp = std::fopen(P.reffunc_add.c_str(), "rb");
  if (!fp) { std::perror("fopen"); return {}; }
  std::vector<std::complex<double>> Hf(N);
  std::vector<double> H_amp(N), H_phase(N);
  // 读取参考函数数据
  for (int n = 0; n < N; ++n)
  {
    float amp = 0.0f;
    if(std::fread(&amp, sizeof(float), 1, fp) != 1) {
      std::fclose(fp);
      return {};
    }
    H_amp[n] = 1.0*amp;
  }
  for (int n = 0; n < N; ++n)
  {
    double re = 0.0, im = 0.0;
    float phase = 0.0f;
    if(std::fread(&phase, sizeof(float), 1, fp) != 1) {
      std::fclose(fp);
      return {};
    }
    re = H_amp[n] * std::cos(1.0*phase);
    im = H_amp[n] * std::sin(1.0*phase);
    Hf[n] = std::complex<double>(re, im);
  }
  std::fclose(fp);
  return Hf;
}

inline std::vector<std::complex<double>>
buildHfSpecFromParams(const Config& P)
{
    const int    N   = P.pulse_len;
    const double fs  = P.fs;
    const double Kr  = (P.Tr > 0.0) ? (P.Br / P.Tr) : 0.0;
    return buildHfSpec(N, fs, Kr);
}

bool rangeCompressFFT(const Config& P,
                      const std::vector<std::complex<double>>& in,
                      std::vector<std::complex<double>>& rc_out)
{
  const int W    = effectivePulseNum(P);
  const int Lraw = P.pulse_len;
  const int M    = P.rg_len;

  std::vector<std::complex<double>> HfSpec;
  if(P.hasRefFunc){
    HfSpec = buildHfSpecFromRefFunc(P);
  }
  else{
    HfSpec = buildHfSpecFromParams(P);
  }
  if ((int)HfSpec.size() != Lraw) return false;

  // 准备 FFTW 缓冲：W×Lraw（行主）
  fftw_complex* buf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (size_t)W * Lraw);
  if (!buf) return false;

  // FFTW planner 不是线程安全的：并发 period 时先串行创建/销毁 plan
  static std::mutex fftw_plan_mutex;
  fftw_plan fwd = nullptr;
  fftw_plan inv = nullptr;
  {
    std::lock_guard<std::mutex> lock(fftw_plan_mutex);

    // plan_many：沿 Lraw 做 FFT（每行一个 FFT），howmany=W
    int rank=1, n[1]={Lraw}, howmany=W;
    int istride=1, ostride=1, idist=Lraw, odist=Lraw;
    int* inembed=n; int* onembed=n;

    fwd = fftw_plan_many_dft(rank, n, howmany, buf, inembed, istride, idist,
                              buf, onembed, ostride, odist,
                              FFTW_FORWARD, FFTW_MEASURE);
    if (!fwd) { fftw_free(buf); return false; }
    inv = fftw_plan_many_dft(rank, n, howmany, buf, inembed, istride, idist,
                              buf, onembed, ostride, odist,
                              FFTW_BACKWARD, FFTW_MEASURE);
    if (!inv) { fftw_destroy_plan(fwd); fftw_free(buf); return false; }
  }

  // 拷贝输入到 buf
  for (int k = 0; k < W; ++k) {
    const std::complex<double>* src = &in[(size_t)k * Lraw];
    fftw_complex* dst = &buf[(size_t)k * Lraw];
    std::memcpy(dst, src, sizeof(std::complex<double>) * (size_t)Lraw);
  }

  // 前向 FFT
  fftw_execute(fwd);

  // 频域乘以 conj(HfSpec)
  for (int k = 0; k < W; ++k) {
    fftw_complex* row = &buf[(size_t)k * Lraw];
    for (int f = 0; f < Lraw; ++f) {
      std::complex<double>& X = *reinterpret_cast<std::complex<double>*>(&row[f]);
      X *= HfSpec[(size_t)f];
    }
  }

  // 逆向 FFT
  fftw_execute(inv);

  // 归一化 + 裁剪到 M，输出到 rc_out（W×M 行主）
  const int Lraw2M = Lraw / M;
  rc_out.resize((size_t)W * M);
  const double scale = 1.0 / double(Lraw);
  for (int k = 0; k < W; ++k) {
    const fftw_complex* row = &buf[(size_t)k * Lraw];
    for (int m = 0; m < M; ++m) {
      rc_out[(size_t)k * M + m] =
        (*reinterpret_cast<const std::complex<double>*>(&row[m * Lraw2M])) * scale;
    }
  }

  {
    std::lock_guard<std::mutex> lock(fftw_plan_mutex);
    fftw_destroy_plan(fwd);
    fftw_destroy_plan(inv);
  }
  fftw_free(buf);
  return true;
}

// 读取一个波位的原始 IQ + 帧头字段（float版本，避免不必要的double转换）
bool readBeamRawFloat(const Config& P, 
  const char* filepath, 
  int beamIdx, 
  std::vector<std::complex<float>>& data,
  std::vector<double>& fw_angle_deg,
  std::vector<double>& utc)
{
  const int hdr  = P.info_len;
  const int Lraw = P.pulse_len;
  const int W    = effectivePulseNum(P);

  const size_t bytesIQperPulse = (size_t)Lraw * 2 /*I/Q*/ * sizeof(float);
  const size_t stride          = (size_t)hdr + bytesIQperPulse;
  const size_t off_bytes       = (size_t)(beamIdx - 1) * W * stride
                               + (size_t)P.skip_az_num * stride;

  FILE* fp = std::fopen(filepath, "rb");
  if (!fp) { std::perror("fopen"); return false; }
  if (std::fseek(fp, (long)off_bytes, SEEK_SET) != 0) {
    std::perror("fseek"); std::fclose(fp); return false;
  }

  data.resize((size_t)W * Lraw);
  fw_angle_deg.assign(W, 0.0f);
  utc.assign(W, 0.0);

  std::vector<uint8_t> header(hdr);
  std::vector<float>   iq_interleaved((size_t)Lraw * 2);

  for (int k = 0; k < W; ++k) {
    // 1) 读帧头
    size_t nr = std::fread(header.data(), 1, hdr, fp);
    if (nr != (size_t)hdr) { std::fclose(fp); return false; }

    // 波束方位角：info(219:220)-[218:219] -> int16 -> /100（度）
    if (hdr > 218) {
      int16_t ang = load_int16_le(&header[218]);
      fw_angle_deg[k] = (double)ang / 100.0;
    } else {
      fw_angle_deg[k] = 0.0;
    }

    // 2) 读 IQ（float32 交错）
    nr = std::fread(iq_interleaved.data(), sizeof(float), (size_t)Lraw * 2, fp);
    if (nr != (size_t)Lraw * 2) { std::fclose(fp); return false; }

    // 交错转复数（保持float，不转double）
    std::complex<float>* row = &data[(size_t)k * Lraw];
    const float* src = iq_interleaved.data();
    for (int n = 0; n < Lraw; ++n) {
      row[n] = std::complex<float>(src[2*n + 0], src[2*n + 1]);
    }

    // 3) UTC 时间解码（与你的 MATLAB 一致）
    // info(37:39) 经 BCD 修正: x - 6*floor(x/16)
    auto bcd_fix = [](uint8_t u)->int { return (int)u - 6 * ((int)u / 16); };
    int hh = 0, mm = 0, ss = 0;
    if (hdr >= 40) {
      hh = bcd_fix(header[36]);  // 注意：MATLAB 的下标从1起，这里从0起，因此(37)->[36]
      mm = bcd_fix(header[37]);  // (38)->[37]
      ss = bcd_fix(header[38]);  // (39)->[38]
    }
    double subsecs = 0.0;
    if (hdr >= 44) {
      uint32_t ds = load_u32_le(&header[39]); // (40:43)->[39..42]
      subsecs = ds / 100e6; // /100e6
    }
    // + 周偏移 + 19 秒
    double t_utc = hh*3600.0 + mm*60.0 + ss + subsecs + P.week * 24.0 * 3600.0 + P.secBias;
    utc[k] = t_utc;
  }

  std::fclose(fp);
  return true;
}

// 读取一个波位的原始 IQ + 帧头字段（与 MATLAB 完全对应）
bool readBeamRaw(const Config& P, 
  const char* filepath, 
  int beamIdx, 
  std::vector<std::complex<double>>& data,
  std::vector<double>& fw_angle_deg,
  std::vector<double>& utc)
{
  const int hdr  = P.info_len;
  const int Lraw = P.pulse_len;
  const int W    = effectivePulseNum(P);

  const size_t bytesIQperPulse = (size_t)Lraw * 2 /*I/Q*/ * sizeof(float);
  const size_t stride          = (size_t)hdr + bytesIQperPulse;
  const size_t off_bytes       = (size_t)(beamIdx - 1) * W * stride
                               + (size_t)P.skip_az_num * stride;

  FILE* fp = std::fopen(filepath, "rb");
  if (!fp) { std::perror("fopen"); return false; }
  if (std::fseek(fp, (long)off_bytes, SEEK_SET) != 0) {
    std::perror("fseek"); std::fclose(fp); return false;
  }

  data.resize((size_t)W * Lraw);
  fw_angle_deg.assign(W, 0.0f);
  utc.assign(W, 0.0);

  std::vector<uint8_t> header(hdr);
  std::vector<float>   iq_interleaved((size_t)Lraw * 2);

  for (int k = 0; k < W; ++k) {
    // 1) 读帧头
    size_t nr = std::fread(header.data(), 1, hdr, fp);
    if (nr != (size_t)hdr) { std::fclose(fp); return false; }

    // 波束方位角：info(219:220)-[218:219] -> int16 -> /100（度）
    if (hdr > 218) {
      int16_t ang = load_int16_le(&header[218]);
      fw_angle_deg[k] = (double)ang / 100.0;
    } else {
      fw_angle_deg[k] = 0.0;
    }

    // 2) 读 IQ（float32 交错）
    nr = std::fread(iq_interleaved.data(), sizeof(float), (size_t)Lraw * 2, fp);
    if (nr != (size_t)Lraw * 2) { std::fclose(fp); return false; }

    // 交错转复数
    std::complex<double>* row = &data[(size_t)k * Lraw];
    const float* src = iq_interleaved.data();
    for (int n = 0; n < Lraw; ++n) {
      row[n] = std::complex<double>(src[2*n + 0], src[2*n + 1]);
    }

    // 3) UTC 时间解码（与你的 MATLAB 一致）
    // info(37:39) 经 BCD 修正: x - 6*floor(x/16)
    auto bcd_fix = [](uint8_t u)->int { return (int)u - 6 * ((int)u / 16); };
    int hh = 0, mm = 0, ss = 0;
    if (hdr >= 40) {
      hh = bcd_fix(header[36]);  // 注意：MATLAB 的下标从1起，这里从0起，因此(37)->[36]
      mm = bcd_fix(header[37]);  // (38)->[37]
      ss = bcd_fix(header[38]);  // (39)->[38]
    }
    double subsecs = 0.0;
    if (hdr >= 44) {
      uint32_t ds = load_u32_le(&header[39]); // (40:43)->[39..42]
      subsecs = ds / 100e6; // /100e6
    }
    // + 周偏移 + 19 秒
    double t_utc = hh*3600.0 + mm*60.0 + ss + subsecs + P.week * 24.0 * 3600.0 + P.secBias;
    utc[k] = t_utc;
  }

  std::fclose(fp);
  return true;
}
