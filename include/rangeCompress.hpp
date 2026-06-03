#pragma once
#include "config_structs.hpp"

bool rangeCompressFFT(const Config &P,
                      const std::vector<std::complex<double>> &in,
                      std::vector<std::complex<double>> &rc_out);

// 读取一个波位的原始 IQ + 帧头字段（与 MATLAB 完全对应）
bool readBeamRaw(const Config &P,
                 const char *filepath,
                 int beamIdx,
                 std::vector<std::complex<double>> &data,
                 std::vector<double> &fw_angle_deg,
                 std::vector<double> &utc);

// 读取一个波位的原始 IQ + 帧头字段（float版本，避免不必要的double转换）
bool readBeamRawFloat(const Config &P,
                      const char *filepath,
                      int beamIdx,
                      std::vector<std::complex<float>> &data,
                      std::vector<double> &fw_angle_deg,
                      std::vector<double> &utc);

bool readBeamUTC(const Config &P,
                 const char *filepath,
                 int beamIdx,
                 std::vector<double> &utc);