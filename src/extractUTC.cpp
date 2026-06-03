#include "config_structs.hpp"
#include "rangeCompress.hpp"
#include "GMTIProcessor.hpp"

#include <vector>
#include <cstdint>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>

static inline int16_t load_int16_le(const uint8_t* p) {
  return static_cast<int16_t>(p[0] | (p[1] << 8));
}
static inline uint32_t load_u32_le(const uint8_t* p) {
  return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

// 读取一个波位的原始 IQ + 帧头字段（与 MATLAB 完全对应）
bool readBeamUTC(const Config& P, 
  const char* filepath, 
  int beamIdx, 
  std::vector<double>& utc)
{
  const int hdr  = P.info_len;
  const int Lraw = P.pulse_len;
  const int W    = P.pulse_num;

  const size_t bytesIQperPulse = (size_t)Lraw * 2 /*I/Q*/ * sizeof(float);
  const size_t stride          = (size_t)hdr + bytesIQperPulse;
  const size_t off_bytes       = (size_t)(beamIdx - 1) * W * stride
                               + (size_t)P.skip_az_num * stride;

  FILE* fp = std::fopen(filepath, "rb");
  if (!fp) { std::perror("fopen"); return false; }
  if (std::fseek(fp, (long)off_bytes, SEEK_SET) != 0) {
    std::perror("fseek"); std::fclose(fp); return false;
  }

  utc.assign(W, 0.0);

  std::vector<uint8_t> header(hdr);
  std::vector<float>   iq_interleaved((size_t)Lraw * 2);

  for (int k = 0; k < W; ++k) {
    // 1) 读帧头
    size_t nr = std::fread(header.data(), 1, hdr, fp);
    if (nr != (size_t)hdr) { std::fclose(fp); return false; }


    // 2) 读 IQ（float32 交错）
    nr = std::fread(iq_interleaved.data(), sizeof(float), (size_t)Lraw * 2, fp);
    if (nr != (size_t)Lraw * 2) { std::fclose(fp); return false; }

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

bool GMTIProcessor::extractPlanePV(const std::vector<std::vector<double>> &POS,
                                    const Config &cfg,
                                    GMTIOutput::Plane &plane)
{
  if (cfg.INFO_Type) {
    return extractPlanePVFromEcho(cfg, plane);
  }

  // Legacy protocol: use external POS file
  std::vector<double> pos_utc;
  const char* filepath = cfg.GMTI_Data_add.c_str();
  if(! readBeamUTC(cfg, filepath, 25, pos_utc) ) {
    return false;
  }

  if(! extractPlanePos(pos_utc, POS, cfg, plane) ) {
    return false;
  }

  return true;
}