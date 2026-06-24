#include "dbs/NewProtocolReader.hpp"

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
                               std::vector<std::complex<float>> &data1,
                               std::vector<std::complex<float>> &data2,
                               std::vector<double> &utc,
                               double &theta_sq,
                               std::vector<std::vector<double>> &posRaw)
{
    constexpr size_t kHeaderBytes = 256;
    constexpr size_t kBytesPerSample = 16; // ch1(I,Q) float32 + ch2(I,Q) float32
    if (cfg.pulse_len <= 0) {
        std::cerr << "[ERR] 新协议 pulse_len 非法: " << cfg.pulse_len << std::endl;
        return false;
    }
    const size_t samples_per_prt = static_cast<size_t>(cfg.pulse_len);
    const size_t prt_bytes =
        kHeaderBytes + samples_per_prt * kBytesPerSample;

    if (cfg.pulse_num <= 0) {
        std::cerr << "[ERR] 新协议 pulse_num 非法: pulse_num=" << cfg.pulse_num << std::endl;
        return false;
    }
    const size_t period_pulses = static_cast<size_t>(cfg.pulse_num);
    size_t read_pulses = period_pulses;
    if (cfg.read_pulse_num > 0) {
        read_pulses = static_cast<size_t>(cfg.read_pulse_num);
    }
    if (read_pulses == 0U || read_pulses > period_pulses) {
        std::cerr << "[ERR] 新协议 read_pulse_num 非法: read_pulse_num="
                  << cfg.read_pulse_num << " pulse_num=" << cfg.pulse_num << std::endl;
        return false;
    }

    size_t pulse_offset = 0U;
    if (cfg.read_pulse_offset >= 0) {
        pulse_offset = static_cast<size_t>(cfg.read_pulse_offset);
    } else {
        pulse_offset = (period_pulses - read_pulses) / 2U;
    }
    if (pulse_offset + read_pulses > period_pulses) {
        std::cerr << "[ERR] 新协议读取窗口越界: pulse_num=" << cfg.pulse_num
                  << " read_pulse_num=" << read_pulses
                  << " read_pulse_offset=" << pulse_offset << std::endl;
        return false;
    }

    const std::string &echoPath = cfg.GMTI_Data_new.empty() ? cfg.GMTI_Data_add : cfg.GMTI_Data_new;
    std::ifstream fp(echoPath, std::ios::binary);
    if (!fp) {
        std::cerr << "[ERR] 无法打开新协议回波文件: " << echoPath << std::endl;
        return false;
    }

    fp.seekg(0, std::ios::end);
    const std::streamsize file_size = fp.tellg();
    if (file_size <= 0 || (static_cast<uint64_t>(file_size) % prt_bytes) != 0U) {
        std::cerr << "[ERR] 新协议回波文件大小非法，不是完整 PRT 包整数倍: " << echoPath << std::endl;
        return false;
    }

    const size_t total_prt =
        static_cast<size_t>(static_cast<uint64_t>(file_size) / prt_bytes);
    const size_t period_start_with_skip = (static_cast<size_t>(beamskip) - 1) * period_pulses + static_cast<size_t>(cfg.skip_az_num);
    const size_t period_start_without_skip = (static_cast<size_t>(beamskip) - 1) * period_pulses;

    size_t start_prt = period_start_with_skip + pulse_offset;
    if (start_prt + read_pulses > total_prt) {
        const size_t alt_start_prt = period_start_without_skip + pulse_offset;
        if (alt_start_prt + read_pulses <= total_prt) {
            std::cerr << "[WARN] 新协议读取带 skip_pulses 越界，自动按已裁剪文件读取（忽略 skip_pulses）: period="
                      << beamskip << " skip_pulses=" << cfg.skip_az_num << std::endl;
            start_prt = alt_start_prt;
        } else {
            std::cerr << "[ERR] 新协议回波文件大小不足，无法读取 period=" << beamskip
                      << " (total_prt=" << total_prt
                      << ", need_end_prt=" << (start_prt + read_pulses)
                      << ", alt_need_end_prt=" << (alt_start_prt + read_pulses) << ")" << std::endl;
            return false;
        }
    }

    const size_t start_byte = start_prt * prt_bytes;
    if (static_cast<uint64_t>(file_size) <
        start_byte + read_pulses * prt_bytes) {
        std::cerr << "[ERR] 新协议回波文件大小不足，无法读取 period=" << beamskip << std::endl;
        return false;
    }
    fp.seekg(static_cast<std::streamoff>(start_byte), std::ios::beg);

    data1.resize(read_pulses * samples_per_prt);
    data2.resize(read_pulses * samples_per_prt);
    utc.resize(read_pulses);
    posRaw.assign(read_pulses, std::vector<double>(7, 0.0));

    std::vector<uint8_t> packet(prt_bytes);
    std::vector<double> fw_angle_deg(read_pulses, 0.0);

    for (size_t k = 0; k < read_pulses; ++k) {
        fp.read(reinterpret_cast<char *>(packet.data()),
                static_cast<std::streamsize>(prt_bytes));
        if (fp.gcount() != static_cast<std::streamsize>(prt_bytes)) {
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
        for (size_t n = 0; n < samples_per_prt; ++n) {
            const size_t off = n * kBytesPerSample;
            const float ch1_i = load_f32_le(payload + off + 0);
            const float ch1_q = load_f32_le(payload + off + 4);
            const float ch2_i = load_f32_le(payload + off + 8);
            const float ch2_q = load_f32_le(payload + off + 12);
            data1[k * samples_per_prt + n] =
                std::complex<float>(ch1_i, ch1_q);
            data2[k * samples_per_prt + n] =
                std::complex<float>(ch2_i, ch2_q);
        }
    }

    theta_sq = fw_angle_deg[fw_angle_deg.size() / 2];
    return true;
}
