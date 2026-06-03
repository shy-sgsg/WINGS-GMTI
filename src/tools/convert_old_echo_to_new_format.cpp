#include "config_structs.hpp"
#include "rangeCompress.hpp"
#include "tinyxml.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct ConvertConfig {
    std::string data_ch1;
    std::string data_ch2;
    std::string data_new;
    std::string pos_file;
    int info_len = 0;
    int pulse_len = 0;
    int pulse_num = 0;
    int skip_az_num = 0;
    int week = 0;
    int sec_bias = 0;
};

struct PosRow {
    double t = 0.0;
    double lat = 0.0;
    double lon = 0.0;
    double alt = 0.0;
    double vn = 0.0;
    double ve = 0.0;
    double vd = 0.0;
};

struct PosSample {
    double lat_deg = 0.0;
    double lon_deg = 0.0;
    double alt_m = 0.0;
    double vn = 0.0;
    double ve = 0.0;
    double vd = 0.0;
};

static inline std::string getText(TiXmlElement *parent, const char *name) {
    if (!parent) {
        return "";
    }
    TiXmlElement *elem = parent->FirstChildElement(name);
    if (!elem || !elem->GetText()) {
        return "";
    }
    return elem->GetText();
}

bool loadConvertConfig(const std::string &xml_path, ConvertConfig &cfg) {
    TiXmlDocument doc;
    if (!doc.LoadFile(xml_path.c_str())) {
        std::cerr << "Load XML failed: " << xml_path << std::endl;
        return false;
    }

    TiXmlElement *root = doc.FirstChildElement("GMTI");
    if (!root) {
        std::cerr << "XML missing <GMTI> root" << std::endl;
        return false;
    }

    TiXmlElement *param = root->FirstChildElement("GMTI_parameter");
    if (!param) {
        std::cerr << "XML missing <GMTI_parameter>" << std::endl;
        return false;
    }

    cfg.data_ch1 = getText(param, "GMTI_data");
    cfg.data_ch2 = getText(param, "GMTI_data2");
    cfg.data_new = getText(param, "GMTI_data_new");
    cfg.pos_file = getText(param, "Plane_POS");

    const std::string info_len = getText(param, "info_len");
    const std::string pulse_len = getText(param, "pulse_len");
    const std::string pulse_num = getText(param, "pulse_num");
    const std::string skip = getText(param, "skip_pulses");
    const std::string week = getText(param, "week_offset");
    const std::string sec_bias = getText(param, "secBias");

    if (cfg.data_ch1.empty() || cfg.data_ch2.empty() || cfg.pos_file.empty() ||
        info_len.empty() || pulse_len.empty() || pulse_num.empty()) {
        std::cerr << "XML missing required fields for conversion" << std::endl;
        return false;
    }

    cfg.info_len = std::stoi(info_len);
    cfg.pulse_len = std::stoi(pulse_len);
    cfg.pulse_num = std::stoi(pulse_num);
    cfg.skip_az_num = skip.empty() ? 0 : std::stoi(skip);
    cfg.week = week.empty() ? 0 : std::stoi(week);
    cfg.sec_bias = sec_bias.empty() ? 0 : std::stoi(sec_bias);

    return true;
}

bool readPosFileLikePOSDataread(const std::string &path, std::vector<PosRow> &rows) {
    const int pos_len = 7;

    std::ifstream fin(path, std::ios::binary);
    if (!fin) {
        std::cerr << "Unable to open POS file: " << path << std::endl;
        return false;
    }

    fin.seekg(0, std::ios::end);
    const std::streamsize file_size = fin.tellg();
    fin.seekg(0, std::ios::beg);

    const std::streamsize pkg_len = pos_len * static_cast<std::streamsize>(sizeof(double));
    if (file_size <= 0 || (file_size % pkg_len) != 0) {
        std::cerr << "POS file size mismatch with 7-double frame format" << std::endl;
        return false;
    }

    const int frame_num = static_cast<int>(file_size / pkg_len);
    rows.resize(static_cast<size_t>(frame_num));

    for (int i = 0; i < frame_num; ++i) {
        double vals[pos_len] = {0.0};
        for (int j = 0; j < pos_len; ++j) {
            fin.read(reinterpret_cast<char *>(&vals[j]), sizeof(double));
            if (!fin) {
                std::cerr << "POS read failed at frame " << i << std::endl;
                return false;
            }
        }

        PosRow row;
        row.t = vals[0];
        row.lat = vals[1];
        row.lon = vals[2];
        row.alt = vals[3];
        row.vn = vals[4];
        row.ve = vals[5];
        row.vd = vals[6];
        rows[static_cast<size_t>(i)] = row;
    }

    return true;
}

static inline double radToDegIfLikely(double x) {
    const double abs_x = std::fabs(x);
    if (abs_x <= 3.2) {
        return x * 180.0 / M_PI;
    }
    return x;
}

PosSample interpPosAtUtc(const std::vector<PosRow> &rows, double t) {
    PosSample out;
    if (rows.empty()) {
        return out;
    }

    auto eval = [&](const PosRow &r) {
        out.lat_deg = radToDegIfLikely(r.lat);
        out.lon_deg = radToDegIfLikely(r.lon);
        out.alt_m = r.alt;
        out.vn = r.vn;
        out.ve = r.ve;
        out.vd = r.vd;
    };

    if (t <= rows.front().t) {
        eval(rows.front());
        return out;
    }
    if (t >= rows.back().t) {
        eval(rows.back());
        return out;
    }

    auto it = std::lower_bound(rows.begin(), rows.end(), t, [](const PosRow &r, double x) {
        return r.t < x;
    });

    const size_t ir = static_cast<size_t>(it - rows.begin());
    const size_t il = ir - 1;

    const PosRow &l = rows[il];
    const PosRow &r = rows[ir];
    const double dt = r.t - l.t;
    const double ratio = (dt == 0.0) ? 0.0 : (t - l.t) / dt;

    auto lerp = [&](double a, double b) {
        return a + ratio * (b - a);
    };

    out.lat_deg = radToDegIfLikely(lerp(l.lat, r.lat));
    out.lon_deg = radToDegIfLikely(lerp(l.lon, r.lon));
    out.alt_m = lerp(l.alt, r.alt);
    out.vn = lerp(l.vn, r.vn);
    out.ve = lerp(l.ve, r.ve);
    out.vd = lerp(l.vd, r.vd);

    return out;
}

static inline uint8_t sat_u8(int x) {
    if (x < 0) {
        return 0;
    }
    if (x > 255) {
        return 255;
    }
    return static_cast<uint8_t>(x);
}

static inline int16_t sat_i16_from_double(double x) {
    if (x > 32767.0) {
        return 32767;
    }
    if (x < -32768.0) {
        return -32768;
    }
    return static_cast<int16_t>(std::lround(x));
}

static inline int16_t sat_i16_from_iq(double x, double scale) {
    const double s = x * scale;
    return sat_i16_from_double(s);
}

static inline void write_u16_le(std::vector<uint8_t> &buf, size_t off, uint16_t v) {
    buf[off + 0] = static_cast<uint8_t>((v >> 0) & 0xFF);
    buf[off + 1] = static_cast<uint8_t>((v >> 8) & 0xFF);
}

static inline void write_i16_le(std::vector<uint8_t> &buf, size_t off, int16_t v) {
    write_u16_le(buf, off, static_cast<uint16_t>(v));
}

static inline void write_u32_le(std::vector<uint8_t> &buf, size_t off, uint32_t v) {
    buf[off + 0] = static_cast<uint8_t>((v >> 0) & 0xFF);
    buf[off + 1] = static_cast<uint8_t>((v >> 8) & 0xFF);
    buf[off + 2] = static_cast<uint8_t>((v >> 16) & 0xFF);
    buf[off + 3] = static_cast<uint8_t>((v >> 24) & 0xFF);
}

static inline void write_u64_le(std::vector<uint8_t> &buf, size_t off, uint64_t v) {
    for (int i = 0; i < 8; ++i) {
        buf[off + static_cast<size_t>(i)] = static_cast<uint8_t>((v >> (8 * i)) & 0xFF);
    }
}

static inline void write_f32_le(std::vector<uint8_t> &buf, size_t off, float v) {
    uint32_t raw = 0;
    std::memcpy(&raw, &v, sizeof(float));
    write_u32_le(buf, off, raw);
}

static inline void write_f64_le(std::vector<uint8_t> &buf, size_t off, double v) {
    uint64_t raw = 0;
    std::memcpy(&raw, &v, sizeof(double));
    write_u64_le(buf, off, raw);
}

uint32_t calcBeamCount(const std::string &path, int info_len, int pulse_len, int pulse_num) {
    std::ifstream fin(path, std::ios::binary | std::ios::ate);
    if (!fin) {
        return 0;
    }

    const std::streamsize file_size = fin.tellg();
    const std::streamsize iq_bytes_per_pulse = static_cast<std::streamsize>(pulse_len) * 2 * static_cast<std::streamsize>(sizeof(float));
    const std::streamsize stride = static_cast<std::streamsize>(info_len) + iq_bytes_per_pulse;
    const std::streamsize beam_bytes = stride * static_cast<std::streamsize>(pulse_num);

    if (beam_bytes <= 0 || file_size <= 0 || (file_size % beam_bytes) != 0) {
        return 0;
    }

    return static_cast<uint32_t>(file_size / beam_bytes);
}

void fillHeaderAndIns(
    std::vector<uint8_t> &pkt,
    uint8_t mode,
    uint32_t prt_counter,
    double utc_sec,
    int week,
    const PosSample &pos,
    double fw_angle_deg,
    uint32_t prt_len)
{
    pkt.assign(prt_len, 0);

    write_u64_le(pkt, 0, 0x5A5A5A5A5A5A5A5AULL);
    pkt[8] = mode;
    write_u32_le(pkt, 9, prt_len);
    write_f32_le(pkt, 16, static_cast<float>(utc_sec));
    write_u32_le(pkt, 20, prt_counter);

    pkt[88] = 0x02;
    pkt[89] = 0x00;
    pkt[90] = 0x40;
    pkt[91] = 0;
    pkt[92] = 0x0B;
    pkt[93] = 0x01;
    write_i16_le(pkt, 94, static_cast<int16_t>(week));

    const double sec_of_day = utc_sec - std::floor(utc_sec / 86400.0) * 86400.0;
    write_u32_le(pkt, 96, static_cast<uint32_t>(std::llround(sec_of_day * 1000.0)));

    write_f64_le(pkt, 104, pos.lat_deg);
    write_f64_le(pkt, 112, pos.lon_deg);
    write_f64_le(pkt, 120, pos.alt_m);

    write_f32_le(pkt, 128, static_cast<float>(pos.vn));
    write_f32_le(pkt, 132, static_cast<float>(pos.ve));
    write_f32_le(pkt, 136, static_cast<float>(pos.vd));
    const double speed = std::sqrt(pos.vn * pos.vn + pos.ve * pos.ve + pos.vd * pos.vd);
    write_f32_le(pkt, 140, static_cast<float>(speed));

    pkt[208] = static_cast<uint8_t>(prt_counter & 0xFFU);
    write_i16_le(pkt, 218, sat_i16_from_double(fw_angle_deg * 100.0));

    write_u64_le(pkt, 248, 0x5B5B5B5B5B5B5B5BULL);
}

void fillGmtiIqPayload(std::vector<uint8_t> &pkt,
                       const std::vector<std::complex<double>> &ch1,
                       const std::vector<std::complex<double>> &ch2,
                       int pulse_idx,
                       int pulse_len,
                       double iq_scale) {
    const int payload_start = 256;
    const int payload_bytes = static_cast<int>(pkt.size()) - payload_start;
    const int bytes_per_sample_2ch = 16; // ch1(float I,Q)=8 + ch2(float I,Q)=8
    const int max_samples = payload_bytes / bytes_per_sample_2ch;

    const int copy_samples = std::min(max_samples, pulse_len);
    const size_t row_off = static_cast<size_t>(pulse_idx) * static_cast<size_t>(pulse_len);

    for (int n = 0; n < copy_samples; ++n) {
        const std::complex<double> &a = ch1[row_off + static_cast<size_t>(n)];
        const std::complex<double> &b = ch2[row_off + static_cast<size_t>(n)];
        const size_t off = static_cast<size_t>(payload_start + n * bytes_per_sample_2ch);

        write_f32_le(pkt, off + 0, static_cast<float>(a.real()));
        write_f32_le(pkt, off + 4, static_cast<float>(a.imag()));
        write_f32_le(pkt, off + 8, static_cast<float>(b.real()));
        write_f32_le(pkt, off + 12, static_cast<float>(b.imag()));
    }
}

bool convertOldToNew(const ConvertConfig &cfg,
                     const std::vector<PosRow> &pos_rows,
                     const std::string &out_file,
                     int beam_start,
                     int beam_count,
                     double iq_scale,
                     uint8_t mode) {
    Config rcfg;
    rcfg.GMTI_Data_add = cfg.data_ch1;
    rcfg.GMTI_Data_add2 = cfg.data_ch2;
    rcfg.info_len = cfg.info_len;
    rcfg.pulse_len = cfg.pulse_len;
    rcfg.pulse_num = cfg.pulse_num;
    rcfg.skip_az_num = cfg.skip_az_num;
    rcfg.week = cfg.week;
    rcfg.secBias = cfg.sec_bias;

    // Diagnostic info: print config and file sizes to help debug "No usable beam after skip_pulses"
    std::cerr << "[DBG] convertOldToNew: data_ch1='" << cfg.data_ch1 << "' data_ch2='" << cfg.data_ch2 << "' pos='" << cfg.pos_file << "'\n";
    std::cerr << "[DBG] cfg.info_len=" << cfg.info_len << " pulse_len=" << cfg.pulse_len << " pulse_num=" << cfg.pulse_num << " skip_az_num=" << cfg.skip_az_num << " week=" << cfg.week << " sec_bias=" << cfg.sec_bias << "\n";

    // report file sizes if available
    auto report_size = [&](const std::string &p, const char *label) {
        std::ifstream f(p, std::ios::binary | std::ios::ate);
        if (!f) {
            std::cerr << "[DBG] cannot open " << label << "='" << p << "'\n";
            return (std::streamsize) -1;
        }
        std::streamsize s = f.tellg();
        std::cerr << "[DBG] " << label << " size=" << s << " bytes\n";
        return s;
    };

    const std::streamsize sz1 = report_size(cfg.data_ch1, "ch1");
    const std::streamsize sz2 = report_size(cfg.data_ch2, "ch2");

    const size_t iq_bytes_per_pulse = (size_t)cfg.pulse_len * 2 /*I/Q*/ * sizeof(float);
    const size_t stride = (size_t)cfg.info_len + iq_bytes_per_pulse;
    const size_t beam_bytes = stride * (size_t)cfg.pulse_num;
    std::cerr << "[DBG] iq_bytes_per_pulse=" << iq_bytes_per_pulse << " stride=" << stride << " beam_bytes(per beam)=" << beam_bytes << "\n";

    const uint32_t total_beams = calcBeamCount(cfg.data_ch1, cfg.info_len, cfg.pulse_len, cfg.pulse_num);
    std::cerr << "[DBG] calcBeamCount returned total_beams=" << total_beams << "\n";
    if (total_beams == 0U) {
        std::cerr << "Cannot infer beam count from old echo file" << std::endl;
        return false;
    }

    // Interpret skip_az_num as global number of pulses to skip from file start
    uint64_t total_pulses = uint64_t(total_beams) * uint64_t(cfg.pulse_num);
    if (uint64_t(cfg.skip_az_num) >= total_pulses) {
        std::cerr << "No usable data: skip_pulses >= total pulses in file" << std::endl;
        if (sz1 > 0) {
            std::cerr << "[DBG] ch1 file bytes=" << sz1 << " beam_bytes(per beam)=" << beam_bytes << " => total_beams_est=" << (sz1 / (std::streamsize)beam_bytes) << "\n";
        }
        if (sz2 > 0) {
            std::cerr << "[DBG] ch2 file bytes=" << sz2 << " beam_bytes(per beam)=" << beam_bytes << " => total_beams_est=" << (sz2 / (std::streamsize)beam_bytes) << "\n";
        }
        std::cerr << "[DBG] Suggestion: ensure <skip_pulses> (global pulses skipped) is less than total_beams * pulse_num." << std::endl;
        return false;
    }

    const int usable_beams = static_cast<int>((total_pulses - uint64_t(cfg.skip_az_num)) / uint64_t(cfg.pulse_num));
    std::cerr << "[DBG] total_pulses=" << total_pulses << " skip_pulses=" << cfg.skip_az_num << " usable_beams=" << usable_beams << "\n";

    const int start = std::max(1, beam_start);
    const int wanted = (beam_count <= 0) ? (usable_beams - (start - 1)) : beam_count;
    const int end = std::min(usable_beams, start + wanted - 1);
    if (start > end) {
        std::cerr << "Invalid beam range" << std::endl;
        return false;
    }

    std::ofstream fout(out_file, std::ios::binary);
    if (!fout) {
        std::cerr << "Unable to open output file: " << out_file << std::endl;
        return false;
    }

    // New protocol: 4096 sample points per PRT, each point stores ch1 I/Q + ch2 I/Q as float32.
    // Header is 256 bytes, payload is 4096 * 16 bytes = 65536 bytes, total PRT length = 65792 bytes.
    const uint32_t prt_len = 256U + 4096U * 16U;
    std::cerr << "[DBG] using new PRT length=" << prt_len << " bytes (256 header + 4096 samples * 16 bytes/sample)\n";

    std::vector<std::complex<double>> data1;
    std::vector<std::complex<double>> data2;
    std::vector<double> utc;
    std::vector<double> fw_angle_deg;
    std::vector<uint8_t> packet(4096, 0);

    uint32_t prt_counter = 0;

    for (int beam_idx = start; beam_idx <= end; ++beam_idx) {
        if (!readBeamRaw(rcfg, cfg.data_ch2.c_str(), beam_idx, data2, fw_angle_deg, utc)) {
            std::cerr << "readBeamRaw failed on channel-2, beam=" << beam_idx << std::endl;
            return false;
        }

        std::vector<double> utc_ch1;
        std::vector<double> fw_angle_dummy;
        if (!readBeamRaw(rcfg, cfg.data_ch1.c_str(), beam_idx, data1, fw_angle_dummy, utc_ch1)) {
            std::cerr << "readBeamRaw failed on channel-1, beam=" << beam_idx << std::endl;
            return false;
        }

        if (data1.size() != data2.size()) {
            std::cerr << "Channel size mismatch at beam=" << beam_idx << std::endl;
            return false;
        }
        if (utc.size() != utc_ch1.size()) {
            std::cerr << "UTC size mismatch at beam=" << beam_idx << std::endl;
            return false;
        }

        for (int k = 0; k < cfg.pulse_num; ++k) {
            const double t = utc[static_cast<size_t>(k)];
            const PosSample pos = interpPosAtUtc(pos_rows, t);
            const double fw = (k < static_cast<int>(fw_angle_deg.size())) ? fw_angle_deg[static_cast<size_t>(k)] : 0.0;

            fillHeaderAndIns(packet, mode, prt_counter, t, cfg.week, pos, fw, prt_len);
            fillGmtiIqPayload(packet, data1, data2, k, cfg.pulse_len, iq_scale);

            fout.write(reinterpret_cast<const char *>(packet.data()), static_cast<std::streamsize>(packet.size()));
            if (!fout) {
                std::cerr << "Write failed at beam=" << beam_idx << " pulse=" << k << std::endl;
                return false;
            }
            ++prt_counter;
        }

        std::cout << "Converted beam " << beam_idx << " / " << end << std::endl;
    }

    std::cout << "Conversion done, total PRT packets: " << prt_counter << std::endl;
    return true;
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <xml_path> [output_new_echo.bin] [beam_start=1] [beam_count=0(all)] [iq_scale=1.0] [mode=5]"
                  << "\nIf output_new_echo.bin is omitted, <GMTI_data_new> from XML will be used."
                  << std::endl;
        return 1;
    }

    const std::string xml_path = argv[1];
    const bool has_output_arg = (argc >= 3);
    const int beam_start = (argc > 3) ? std::max(1, std::atoi(argv[3])) : 1;
    const int beam_count = (argc > 4) ? std::atoi(argv[4]) : 0;
    const double iq_scale = (argc > 5) ? std::atof(argv[5]) : 1.0;
    const int mode_i = (argc > 6) ? std::atoi(argv[6]) : 5;
    const uint8_t mode = sat_u8(mode_i);

    ConvertConfig cfg;
    if (!loadConvertConfig(xml_path, cfg)) {
        return 2;
    }

    std::string out_file;
    if (has_output_arg) {
        out_file = argv[2];
    } else {
        out_file = cfg.data_new;
    }
    if (out_file.empty()) {
        std::cerr << "Output path is empty. Provide argv[2] or set <GMTI_data_new> in XML." << std::endl;
        return 2;
    }

    std::vector<PosRow> pos_rows;
    if (!readPosFileLikePOSDataread(cfg.pos_file, pos_rows)) {
        return 3;
    }

    if (!convertOldToNew(cfg, pos_rows, out_file, beam_start, beam_count, iq_scale, mode)) {
        return 4;
    }

    return 0;
}
