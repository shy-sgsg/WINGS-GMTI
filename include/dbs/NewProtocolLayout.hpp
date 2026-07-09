#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>

namespace gmti {
namespace new_protocol {

static const std::size_t kHeaderBytes = 256;
static const std::size_t kBytesPerFloatIq = sizeof(float);
static const std::size_t kBytesPerInt16Iq = sizeof(int16_t);

static const std::size_t kOffMagicHead = 0;
static const std::size_t kOffVersion = 8;
static const std::size_t kOffPrtLen = 9;
static const std::size_t kOffUtc = 16;
static const std::size_t kOffPrtCounter = 20;
static const std::size_t kOffLatDeg = 104;
static const std::size_t kOffLonDeg = 112;
static const std::size_t kOffHeightM = 120;
static const std::size_t kOffVnMps = 128;
static const std::size_t kOffVeMps = 132;
static const std::size_t kOffVdMps = 136;
static const std::size_t kOffSpeedMps = 140;
static const std::size_t kOffPrtLowByte = 208;
static const std::size_t kOffThetaDegX100 = 218;
static const std::size_t kOffMagicTail = 248;

inline bool isInt16IqType(const std::string &iq_data_type)
{
    return iq_data_type == "int16" || iq_data_type == "iq_int16" ||
           iq_data_type == "short" || iq_data_type == "s16";
}

inline std::size_t bytesPerIq(const std::string &iq_data_type)
{
    return isInt16IqType(iq_data_type) ? kBytesPerInt16Iq : kBytesPerFloatIq;
}

inline std::size_t bytesPerChannel(const std::string &iq_data_type)
{
    return 2U * bytesPerIq(iq_data_type);
}

inline std::size_t sampleBytes(std::size_t channel_count, const std::string &iq_data_type)
{
    return channel_count * bytesPerChannel(iq_data_type);
}

inline std::size_t packetBytes(std::size_t samples_per_prt,
                               std::size_t channel_count,
                               const std::string &iq_data_type)
{
    return kHeaderBytes + samples_per_prt * sampleBytes(channel_count, iq_data_type);
}

inline std::size_t channelOffset(std::size_t channel_index_1based, const std::string &iq_data_type)
{
    return (channel_index_1based - 1U) * bytesPerChannel(iq_data_type);
}

inline uint16_t loadU16LE(const uint8_t *p)
{
    return static_cast<uint16_t>(p[0] | (static_cast<uint16_t>(p[1]) << 8));
}

inline int16_t loadI16LE(const uint8_t *p)
{
    return static_cast<int16_t>(loadU16LE(p));
}

inline uint32_t loadU32LE(const uint8_t *p)
{
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

inline uint64_t loadU64LE(const uint8_t *p)
{
    return static_cast<uint64_t>(p[0]) |
           (static_cast<uint64_t>(p[1]) << 8) |
           (static_cast<uint64_t>(p[2]) << 16) |
           (static_cast<uint64_t>(p[3]) << 24) |
           (static_cast<uint64_t>(p[4]) << 32) |
           (static_cast<uint64_t>(p[5]) << 40) |
           (static_cast<uint64_t>(p[6]) << 48) |
           (static_cast<uint64_t>(p[7]) << 56);
}

inline float loadF32LE(const uint8_t *p)
{
    const uint32_t raw = loadU32LE(p);
    float v = 0.0f;
    std::memcpy(&v, &raw, sizeof(float));
    return v;
}

inline float loadIqAsFloat(const uint8_t *p, const std::string &iq_data_type)
{
    return isInt16IqType(iq_data_type)
        ? static_cast<float>(loadI16LE(p))
        : loadF32LE(p);
}

inline int16_t satI16FromFloat(float x)
{
    if (x > 32767.0f) return 32767;
    if (x < -32768.0f) return -32768;
    return static_cast<int16_t>(std::floor(x + (x >= 0.0f ? 0.5f : -0.5f)));
}

inline double loadF64LE(const uint8_t *p)
{
    const uint64_t raw = loadU64LE(p);
    double v = 0.0;
    std::memcpy(&v, &raw, sizeof(double));
    return v;
}

inline void storeU16LE(uint8_t *p, uint16_t v)
{
    p[0] = static_cast<uint8_t>(v & 0xffU);
    p[1] = static_cast<uint8_t>((v >> 8) & 0xffU);
}

inline void storeI16LE(uint8_t *p, int16_t v)
{
    storeU16LE(p, static_cast<uint16_t>(v));
}

inline void storeU32LE(uint8_t *p, uint32_t v)
{
    for (int i = 0; i < 4; ++i) {
        p[i] = static_cast<uint8_t>((v >> (8 * i)) & 0xffU);
    }
}

inline void storeU64LE(uint8_t *p, uint64_t v)
{
    for (int i = 0; i < 8; ++i) {
        p[i] = static_cast<uint8_t>((v >> (8 * i)) & 0xffU);
    }
}

inline void storeF32LE(uint8_t *p, float v)
{
    uint32_t raw = 0;
    std::memcpy(&raw, &v, sizeof(float));
    storeU32LE(p, raw);
}

inline void storeIqFromFloat(uint8_t *p, const std::string &iq_data_type, float v)
{
    if (isInt16IqType(iq_data_type)) {
        storeI16LE(p, satI16FromFloat(v));
    } else {
        storeF32LE(p, v);
    }
}

inline void storeF64LE(uint8_t *p, double v)
{
    uint64_t raw = 0;
    std::memcpy(&raw, &v, sizeof(double));
    storeU64LE(p, raw);
}

struct HeaderSample {
    double utc = 0.0;
    double lat_deg = 0.0;
    double lon_deg = 0.0;
    double height_m = 0.0;
    double vn_mps = 0.0;
    double ve_mps = 0.0;
    double vd_mps = 0.0;
    double theta_cmd_deg = 0.0;
};

inline HeaderSample readHeaderSample(const uint8_t *hdr)
{
    HeaderSample s;
    s.utc = static_cast<double>(loadF32LE(hdr + kOffUtc));
    s.lat_deg = loadF64LE(hdr + kOffLatDeg);
    s.lon_deg = loadF64LE(hdr + kOffLonDeg);
    s.height_m = loadF64LE(hdr + kOffHeightM);
    s.vn_mps = static_cast<double>(loadF32LE(hdr + kOffVnMps));
    s.ve_mps = static_cast<double>(loadF32LE(hdr + kOffVeMps));
    s.vd_mps = static_cast<double>(loadF32LE(hdr + kOffVdMps));
    s.theta_cmd_deg = static_cast<double>(loadI16LE(hdr + kOffThetaDegX100)) / 100.0;
    return s;
}

} // namespace new_protocol
} // namespace gmti
