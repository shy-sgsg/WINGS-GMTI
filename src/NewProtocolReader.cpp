#include "NewProtocolReader.hpp"

#include <cmath>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <iostream>

namespace {

static inline uint16_t load_u16_le(const uint8_t *p)
{
    return static_cast<uint16_t>(p[0] | (static_cast<uint16_t>(p[1]) << 8));
}

static inline int16_t load_i16_le(const uint8_t *p)
{
    return static_cast<int16_t>(load_u16_le(p));
}

static inline uint32_t load_u32_le(const uint8_t *p)
{
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

static inline uint64_t load_u64_le(const uint8_t *p)
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

static inline float load_f32_le(const uint8_t *p)
{
    uint32_t raw = load_u32_le(p);
    float v = 0.0f;
    std::memcpy(&v, &raw, sizeof(float));
    return v;
}

static inline double load_f64_le(const uint8_t *p)
{
    uint64_t raw = load_u64_le(p);
    double v = 0.0;
    std::memcpy(&v, &raw, sizeof(double));
    return v;
}

} // namespace

bool readPulseBlockNewProtocol(const Config &cfg,
                               int beamskip,
                               std::vector<std::complex<double>> &data1,
                               std::vector<std::complex<double>> &data2,
                               std::vector<double> &utc,
                               double &theta_sq,
                               std::vector<std::vector<double>> &posRaw)
{
    constexpr size_t kHeaderBytes = 256;
    constexpr size_t kSamplesPerPrt = 4096;
    constexpr size_t kBytesPerSample = 16; // ch1(I,Q) float32 + ch2(I,Q) float32
    constexpr size_t kPrtBytes = kHeaderBytes + kSamplesPerPrt * kBytesPerSample;

    if (cfg.pulse_len != static_cast<int>(kSamplesPerPrt)) {
        std::cerr << "[ERR] 新协议要求 pulse_len=4096，当前 pulse_len=" << cfg.pulse_len << std::endl;
        return false;
    }

    const size_t W = static_cast<size_t>(cfg.pulse_num);

    const std::string &echoPath = cfg.GMTI_Data_new.empty() ? cfg.GMTI_Data_add : cfg.GMTI_Data_new;
    std::ifstream fp(echoPath, std::ios::binary);
    if (!fp) {
        std::cerr << "[ERR] 无法打开新协议回波文件: " << echoPath << std::endl;
        return false;
    }

    fp.seekg(0, std::ios::end);
    const std::streamsize file_size = fp.tellg();
    if (file_size <= 0 || (static_cast<uint64_t>(file_size) % kPrtBytes) != 0U) {
        std::cerr << "[ERR] 新协议回波文件大小非法，不是完整 PRT 包整数倍: " << echoPath << std::endl;
        return false;
    }

    const size_t total_prt = static_cast<size_t>(static_cast<uint64_t>(file_size) / kPrtBytes);
    const size_t start_with_skip = (static_cast<size_t>(beamskip) - 1) * W + static_cast<size_t>(cfg.skip_az_num);
    const size_t start_without_skip = (static_cast<size_t>(beamskip) - 1) * W;

    size_t start_prt = start_with_skip;
    if (start_with_skip + W > total_prt) {
        if (start_without_skip + W <= total_prt) {
            std::cerr << "[WARN] 新协议读取带 skip_pulses 越界，自动按已裁剪文件读取（忽略 skip_pulses）: period="
                      << beamskip << " skip_pulses=" << cfg.skip_az_num << std::endl;
            start_prt = start_without_skip;
        } else {
            std::cerr << "[ERR] 新协议回波文件大小不足，无法读取 period=" << beamskip
                      << " (total_prt=" << total_prt
                      << ", need_end_prt=" << (start_with_skip + W)
                      << ", alt_need_end_prt=" << (start_without_skip + W) << ")" << std::endl;
            return false;
        }
    }

    const size_t start_byte = start_prt * kPrtBytes;
    if (static_cast<uint64_t>(file_size) < start_byte + W * kPrtBytes) {
        std::cerr << "[ERR] 新协议回波文件大小不足，无法读取 period=" << beamskip << std::endl;
        return false;
    }
    fp.seekg(static_cast<std::streamoff>(start_byte), std::ios::beg);

    data1.resize(W * kSamplesPerPrt);
    data2.resize(W * kSamplesPerPrt);
    utc.resize(W);
    posRaw.assign(W, std::vector<double>(7, 0.0));

    std::vector<uint8_t> packet(kPrtBytes);
    std::vector<double> fw_angle_deg(W, 0.0);

    for (size_t k = 0; k < W; ++k) {
        fp.read(reinterpret_cast<char *>(packet.data()), static_cast<std::streamsize>(kPrtBytes));
        if (fp.gcount() != static_cast<std::streamsize>(kPrtBytes)) {
            std::cerr << "[ERR] 读取新协议 PRT 失败, period=" << beamskip << " pulse=" << k << std::endl;
            return false;
        }

        const uint8_t *hdr = packet.data();
        utc[k] = static_cast<double>(load_f32_le(hdr + 16));
        posRaw[k][0] = utc[k];
        posRaw[k][1] = load_f64_le(hdr + 104) * M_PI / 180.0;
        posRaw[k][2] = load_f64_le(hdr + 112) * M_PI / 180.0;
        posRaw[k][3] = load_f64_le(hdr + 120);
        posRaw[k][4] = static_cast<double>(load_f32_le(hdr + 128));
        posRaw[k][5] = static_cast<double>(load_f32_le(hdr + 132));
        posRaw[k][6] = static_cast<double>(load_f32_le(hdr + 136));
        fw_angle_deg[k] = static_cast<double>(load_i16_le(hdr + 218)) / 100.0;

        const uint8_t *payload = hdr + kHeaderBytes;
        for (size_t n = 0; n < kSamplesPerPrt; ++n) {
            const size_t off = n * kBytesPerSample;
            const float ch1_i = load_f32_le(payload + off + 0);
            const float ch1_q = load_f32_le(payload + off + 4);
            const float ch2_i = load_f32_le(payload + off + 8);
            const float ch2_q = load_f32_le(payload + off + 12);
            data1[k * kSamplesPerPrt + n] = std::complex<double>(ch1_i, ch1_q);
            data2[k * kSamplesPerPrt + n] = std::complex<double>(ch2_i, ch2_q);
        }
    }

    theta_sq = fw_angle_deg[fw_angle_deg.size() / 2];
    return true;
}
