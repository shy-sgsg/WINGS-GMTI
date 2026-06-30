#include "stage1_writer.h"

#include "tinyxml.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace gmti {
namespace sim_stage1 {

namespace {

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

static inline int16_t satI16(double x)
{
    if (x > 32767.0) return 32767;
    if (x < -32768.0) return -32768;
    return static_cast<int16_t>(std::lround(x));
}

static inline void putU16(std::vector<uint8_t> &b, size_t off, uint16_t v)
{
    b[off] = static_cast<uint8_t>(v & 0xffU);
    b[off + 1] = static_cast<uint8_t>((v >> 8) & 0xffU);
}

static inline void putI16(std::vector<uint8_t> &b, size_t off, int16_t v)
{
    putU16(b, off, static_cast<uint16_t>(v));
}

static inline void putU32(std::vector<uint8_t> &b, size_t off, uint32_t v)
{
    for (int i = 0; i < 4; ++i) b[off + static_cast<size_t>(i)] = static_cast<uint8_t>((v >> (8 * i)) & 0xffU);
}

static inline void putU64(std::vector<uint8_t> &b, size_t off, uint64_t v)
{
    for (int i = 0; i < 8; ++i) b[off + static_cast<size_t>(i)] = static_cast<uint8_t>((v >> (8 * i)) & 0xffU);
}

static inline void putF32(std::vector<uint8_t> &b, size_t off, float v)
{
    uint32_t raw = 0;
    std::memcpy(&raw, &v, sizeof(float));
    putU32(b, off, raw);
}

static inline void putF64(std::vector<uint8_t> &b, size_t off, double v)
{
    uint64_t raw = 0;
    std::memcpy(&raw, &v, sizeof(double));
    putU64(b, off, raw);
}

std::complex<float> lerp(const std::complex<float> &a, const std::complex<float> &b, double w)
{
    const float wf = static_cast<float>(w);
    return a * (1.0f - wf) + b * wf;
}

double lerpDouble(double a, double b, double w)
{
    return a + (b - a) * w;
}

double radToDegIfLikely(double x)
{
    return std::fabs(x) <= 3.2 ? x * 180.0 / M_PI : x;
}

double degToRadIfLikely(double x)
{
    return std::fabs(x) > M_PI ? x * M_PI / 180.0 : x;
}

bool readPosRows(const std::string &path, std::vector<PosRow> &rows, std::string &err)
{
    std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
    if (!in) {
        err = "failed to open POS file: " + path;
        return false;
    }
    const std::streamsize bytes = in.tellg();
    const std::streamsize row_bytes = static_cast<std::streamsize>(7 * sizeof(double));
    if (bytes <= 0 || (bytes % row_bytes) != 0) {
        err = "POS file size is not a 7-double row multiple: " + path;
        return false;
    }
    in.seekg(0, std::ios::beg);
    rows.resize(static_cast<size_t>(bytes / row_bytes));
    for (size_t i = 0; i < rows.size(); ++i) {
        double v[7] = {0.0};
        for (int j = 0; j < 7; ++j) {
            in.read(reinterpret_cast<char *>(&v[j]), sizeof(double));
            if (!in) {
                err = "failed to read POS row";
                return false;
            }
        }
        rows[i].t = v[0];
        rows[i].lat = v[1];
        rows[i].lon = v[2];
        rows[i].alt = v[3];
        rows[i].vn = v[4];
        rows[i].ve = v[5];
        rows[i].vd = v[6];
    }
    return true;
}

PosSample interpPosAtUtc(const std::vector<PosRow> &rows, double t)
{
    const double earth_radius_m = 6378137.0;
    PosSample out;
    if (rows.empty()) return out;
    auto eval = [&](const PosRow &r, size_t vel_left, size_t vel_right) {
        out.lat_deg = radToDegIfLikely(r.lat);
        out.lon_deg = radToDegIfLikely(r.lon);
        out.alt_m = r.alt;
        const PosRow &a = rows[vel_left];
        const PosRow &b = rows[vel_right];
        const double dt = b.t - a.t;
        if (std::fabs(dt) > 1e-9) {
            const double lat_a = degToRadIfLikely(a.lat);
            const double lat_b = degToRadIfLikely(b.lat);
            const double lon_a = degToRadIfLikely(a.lon);
            const double lon_b = degToRadIfLikely(b.lon);
            const double lat_mid = 0.5 * (lat_a + lat_b);
            out.vn = (lat_b - lat_a) * earth_radius_m / dt;
            out.ve = (lon_b - lon_a) * earth_radius_m * std::cos(lat_mid) / dt;
            out.vd = -(b.alt - a.alt) / dt;
        } else {
            out.vn = r.vn;
            out.ve = r.ve;
            out.vd = r.vd;
        }
    };
    if (t <= rows.front().t) {
        eval(rows.front(), 0, rows.size() > 1 ? 1 : 0);
        return out;
    }
    if (t >= rows.back().t) {
        eval(rows.back(), rows.size() > 1 ? rows.size() - 2 : 0, rows.size() - 1);
        return out;
    }
    size_t right = 1;
    while (right < rows.size() && rows[right].t < t) {
        ++right;
    }
    const size_t left = right - 1;
    const double denom = rows[right].t - rows[left].t;
    const double w = (std::fabs(denom) < 1e-9) ? 0.0 : (t - rows[left].t) / denom;
    out.lat_deg = radToDegIfLikely(lerpDouble(rows[left].lat, rows[right].lat, w));
    out.lon_deg = radToDegIfLikely(lerpDouble(rows[left].lon, rows[right].lon, w));
    out.alt_m = lerpDouble(rows[left].alt, rows[right].alt, w);
    const double dt = rows[right].t - rows[left].t;
    if (std::fabs(dt) > 1e-9) {
        const double lat_l = degToRadIfLikely(rows[left].lat);
        const double lat_r = degToRadIfLikely(rows[right].lat);
        const double lon_l = degToRadIfLikely(rows[left].lon);
        const double lon_r = degToRadIfLikely(rows[right].lon);
        const double lat_mid = 0.5 * (lat_l + lat_r);
        out.vn = (lat_r - lat_l) * earth_radius_m / dt;
        out.ve = (lon_r - lon_l) * earth_radius_m * std::cos(lat_mid) / dt;
        out.vd = -(rows[right].alt - rows[left].alt) / dt;
    } else {
        out.vn = lerpDouble(rows[left].vn, rows[right].vn, w);
        out.ve = lerpDouble(rows[left].ve, rows[right].ve, w);
        out.vd = lerpDouble(rows[left].vd, rows[right].vd, w);
    }
    return out;
}

void fillHeader(std::vector<uint8_t> &pkt,
                uint32_t prt_counter,
                double utc,
                int week,
                const PosSample &pos,
                double theta_new_deg,
                uint32_t prt_len)
{
    std::fill(pkt.begin(), pkt.end(), 0U);
    putU64(pkt, 0, 0x5A5A5A5A5A5A5A5AULL);
    pkt[8] = 5;
    putU32(pkt, 9, prt_len);
    putF32(pkt, 16, static_cast<float>(utc));
    putU32(pkt, 20, prt_counter);
    pkt[88] = 0x02;
    pkt[90] = 0x40;
    pkt[92] = 0x0B;
    pkt[93] = 0x01;
    putI16(pkt, 94, static_cast<int16_t>(week));
    const double sec_of_day = utc - std::floor(utc / 86400.0) * 86400.0;
    putU32(pkt, 96, static_cast<uint32_t>(std::llround(sec_of_day * 1000.0)));
    putF64(pkt, 104, pos.lat_deg);
    putF64(pkt, 112, pos.lon_deg);
    putF64(pkt, 120, pos.alt_m);
    putF32(pkt, 128, static_cast<float>(pos.vn));
    putF32(pkt, 132, static_cast<float>(pos.ve));
    putF32(pkt, 136, static_cast<float>(pos.vd));
    const double speed = std::sqrt(pos.vn * pos.vn + pos.ve * pos.ve + pos.vd * pos.vd);
    putF32(pkt, 140, static_cast<float>(speed));
    pkt[208] = static_cast<uint8_t>(prt_counter & 0xffU);
    putI16(pkt, 218, satI16(theta_new_deg * 100.0));
    putU64(pkt, 248, 0x5B5B5B5B5B5B5B5BULL);
}

void updateAgg(const std::vector<std::complex<float> > &x, double &sum_abs, double &sum_power,
               double &max_abs, bool &has_nan, bool &has_inf)
{
    for (size_t i = 0; i < x.size(); ++i) {
        const float re = x[i].real();
        const float im = x[i].imag();
        if (std::isnan(re) || std::isnan(im)) has_nan = true;
        if (std::isinf(re) || std::isinf(im)) has_inf = true;
        const double a = std::abs(x[i]);
        sum_abs += a;
        sum_power += a * a;
        if (a > max_abs) max_abs = a;
    }
}

} // namespace

bool generateStage1Data(const Stage1OldSystemConfig &old_cfg,
                        const Stage1NewSystemConfig &new_cfg,
                        const Stage1RunOptions &opt,
                        const std::vector<BeamMapEntry> &beam_map,
                        const std::vector<PulseMapEntry> &pulse_map,
                        RangeFftZeroPadResizer &resizer,
                        const std::string &out_file,
                        GenerationStats &stats,
                        std::string &err)
{
    std::ofstream out(out_file.c_str(), std::ios::binary);
    if (!out) {
        err = "failed to open output data file: " + out_file;
        return false;
    }
    std::ofstream resize_stats(pathJoin(pathJoin(opt.output_dir, "debug"),
                                         "range_resize_stats.csv").c_str(),
                               std::ios::out | std::ios::app);
    if (!resize_stats) {
        err = "failed to open range_resize_stats.csv for real data stats";
        return false;
    }
    const uint32_t prt_len = static_cast<uint32_t>(256 + new_cfg.ddc_len_new * 16);
    std::vector<uint8_t> packet(prt_len);
    std::vector<std::complex<float> > row_old_ch1(static_cast<size_t>(old_cfg.pulse_len));
    std::vector<std::complex<float> > row_old_ch2(static_cast<size_t>(old_cfg.pulse_len));
    std::vector<std::complex<float> > row_new_ch1;
    std::vector<std::complex<float> > row_new_ch2;
    double sum_abs1 = 0.0, sum_pow1 = 0.0, max1 = 0.0;
    double sum_abs2 = 0.0, sum_pow2 = 0.0, max2 = 0.0;
    uint64_t sample_count = 0;
    uint32_t prt_counter = 0;
    const int gen_beam_start = std::max(0, opt.beam_start);
    const int gen_beam_count = (opt.beam_count > 0) ? opt.beam_count : new_cfg.beam_count - gen_beam_start;
    const int gen_beam_end = std::min(new_cfg.beam_count, gen_beam_start + gen_beam_count);
    std::vector<PosRow> pos_rows;
    if (!readPosRows(old_cfg.pos_path, pos_rows, err)) {
        return false;
    }

    for (int pidx = 0; pidx < opt.period_count; ++pidx) {
        const int period = opt.period_start + pidx;
        for (int nb = gen_beam_start; nb < gen_beam_end; ++nb) {
            const BeamMapEntry &bm = beam_map[static_cast<size_t>(nb)];
            OldBlock l1, r1, l2, r2;
            if (!readOldBlock(old_cfg, period, bm.source_left_beam_index, 1, l1, err)) return false;
            if (!readOldBlock(old_cfg, period, bm.source_right_beam_index, 1, r1, err)) return false;
            if (!readOldBlock(old_cfg, period, bm.source_left_beam_index, 2, l2, err)) return false;
            if (!readOldBlock(old_cfg, period, bm.source_right_beam_index, 2, r2, err)) return false;
            for (int np = 0; np < new_cfg.pulse_num_new; ++np) {
                const PulseMapEntry &pm = pulse_map[static_cast<size_t>(np)];
                const size_t ll = static_cast<size_t>(pm.old_left_index) * static_cast<size_t>(old_cfg.pulse_len);
                const size_t rr = static_cast<size_t>(pm.old_right_index) * static_cast<size_t>(old_cfg.pulse_len);
                for (int n = 0; n < old_cfg.pulse_len; ++n) {
                    const size_t o = static_cast<size_t>(n);
                    const std::complex<float> a1 = lerp(l1.samples[ll + o], r1.samples[ll + o], bm.interp_weight);
                    const std::complex<float> b1 = lerp(l1.samples[rr + o], r1.samples[rr + o], bm.interp_weight);
                    const std::complex<float> a2 = lerp(l2.samples[ll + o], r2.samples[ll + o], bm.interp_weight);
                    const std::complex<float> b2 = lerp(l2.samples[rr + o], r2.samples[rr + o], bm.interp_weight);
                    row_old_ch1[o] = lerp(a1, b1, pm.weight);
                    row_old_ch2[o] = lerp(a2, b2, pm.weight);
                }
                if (!resizer.resize(row_old_ch1.data(), row_new_ch1) ||
                    !resizer.resize(row_old_ch2.data(), row_new_ch2)) {
                    err = "range resize failed while generating data";
                    return false;
                }
                if (pidx == 0 && (nb - gen_beam_start) < 3 && np < 3) {
                    const SignalStats in1 = computeStats(row_old_ch1);
                    const SignalStats out1 = computeStats(row_new_ch1);
                    const SignalStats in2 = computeStats(row_old_ch2);
                    const SignalStats out2 = computeStats(row_new_ch2);
                    resize_stats << "real_p" << period << "_b" << nb << "_q" << np << "_ch1,"
                                 << in1.mean_abs << "," << in1.rms << "," << in1.max_abs << ","
                                 << out1.mean_abs << "," << out1.rms << "," << out1.max_abs << ","
                                 << boolText(out1.has_nan) << "," << boolText(out1.has_inf) << "\n";
                    resize_stats << "real_p" << period << "_b" << nb << "_q" << np << "_ch2,"
                                 << in2.mean_abs << "," << in2.rms << "," << in2.max_abs << ","
                                 << out2.mean_abs << "," << out2.rms << "," << out2.max_abs << ","
                                 << boolText(out2.has_nan) << "," << boolText(out2.has_inf) << "\n";
                }
                const double utc = lerp(std::complex<float>(static_cast<float>(l1.utc[static_cast<size_t>(pm.old_left_index)]), 0.0f),
                                        std::complex<float>(static_cast<float>(l1.utc[static_cast<size_t>(pm.old_right_index)]), 0.0f),
                                        pm.weight).real();
                const PosSample pos = interpPosAtUtc(pos_rows, utc);
                fillHeader(packet, prt_counter, utc, old_cfg.week_offset, pos, bm.theta_new_deg, prt_len);
                for (int n = 0; n < new_cfg.ddc_len_new; ++n) {
                    const size_t poff = 256U + static_cast<size_t>(n) * 16U;
                    const std::complex<float> z1 = (opt.channel_mode == "ch2") ? std::complex<float>(0, 0) : row_new_ch1[static_cast<size_t>(n)];
                    const std::complex<float> z2 = (opt.channel_mode == "ch1") ? std::complex<float>(0, 0) : row_new_ch2[static_cast<size_t>(n)];
                    putF32(packet, poff + 0, z1.real());
                    putF32(packet, poff + 4, z1.imag());
                    putF32(packet, poff + 8, z2.real());
                    putF32(packet, poff + 12, z2.imag());
                }
                out.write(reinterpret_cast<const char *>(packet.data()), static_cast<std::streamsize>(packet.size()));
                if (!out) {
                    err = "write failed for stage1 data";
                    return false;
                }
                updateAgg(row_new_ch1, sum_abs1, sum_pow1, max1, stats.has_nan, stats.has_inf);
                updateAgg(row_new_ch2, sum_abs2, sum_pow2, max2, stats.has_nan, stats.has_inf);
                sample_count += static_cast<uint64_t>(new_cfg.ddc_len_new);
                ++stats.packets_written;
                ++prt_counter;
            }
            std::cout << "[stage1] period=" << period << " beam=" << nb << " done\n";
        }
    }
    stats.output_bytes = stats.packets_written * static_cast<uint64_t>(prt_len);
    if (sample_count > 0U) {
        stats.mean_abs_ch1 = sum_abs1 / static_cast<double>(sample_count);
        stats.rms_ch1 = std::sqrt(sum_pow1 / static_cast<double>(sample_count));
        stats.max_abs_ch1 = max1;
        stats.mean_abs_ch2 = sum_abs2 / static_cast<double>(sample_count);
        stats.rms_ch2 = std::sqrt(sum_pow2 / static_cast<double>(sample_count));
        stats.max_abs_ch2 = max2;
    }
    return true;
}

bool writeStage1ConfigXml(const Stage1OldSystemConfig &old_cfg,
                          const Stage1NewSystemConfig &new_cfg,
                          const Stage1RunOptions &opt,
                          const std::string &data_file,
                          std::string &err)
{
    const std::string path = pathJoin(pathJoin(opt.output_dir, "config"), "temp_config_stage1_newsystem.xml");
    TiXmlDocument doc;
    if (!doc.LoadFile(old_cfg.xml_path.c_str())) {
        err = "failed to load old XML template: " + old_cfg.xml_path;
        return false;
    }
    TiXmlElement *root = doc.FirstChildElement("GMTI");
    TiXmlElement *param = root ? root->FirstChildElement("GMTI_parameter") : nullptr;
    if (!param) {
        err = "old XML template missing GMTI_parameter";
        return false;
    }

    auto setNode = [&](const char *name, const std::string &value) {
        TiXmlElement *e = param->FirstChildElement(name);
        if (!e) {
            e = new TiXmlElement(name);
            param->LinkEndChild(e);
        }
        e->Clear();
        e->LinkEndChild(new TiXmlText(value.c_str()));
    };
    auto setInt = [&](const char *name, int value) {
        std::ostringstream oss;
        oss << value;
        setNode(name, oss.str());
    };
    auto setDouble = [&](const char *name, double value) {
        std::ostringstream oss;
        oss << value;
        setNode(name, oss.str());
    };

    setNode("result_add", pathJoin(opt.output_dir, "algorithm_result"));
    setNode("INFO_Type", "1");
    setNode("GMTI_data_new", data_file);
    setNode("isSeparated", "separate");
    setInt("info_len", 256);
    setInt("pulse_len", new_cfg.ddc_len_new);
    setInt("rg_len", new_cfg.pc_crop_len);
    setInt("range_fft_len", new_cfg.fft_len_new);
    setInt("range_crop_start", new_cfg.pc_crop_start);
    setInt("range_compress_len", new_cfg.pc_crop_len);
    setInt("pulse_num", new_cfg.pulse_num_new);
    setInt("read_pulse_num", new_cfg.pulse_num_new);
    setInt("read_pulse_offset", 0);
    setDouble("fc", new_cfg.fc_new_ghz);
    setDouble("Br", new_cfg.br_new_mhz);
    setDouble("fs", new_cfg.fs_new_mhz);
    setDouble("Tr", new_cfg.tr_new_us);
    setDouble("PRF", new_cfg.prf_new_hz);
    setInt("skip_pulses", 0);
    setInt("wavepos_st", 1);
    setInt("wavepos_ed", new_cfg.beam_count);
    setInt("wavepos_skip", 1);
    setDouble("scan_min_deg", new_cfg.scan_min_deg);
    setDouble("scan_max_deg", new_cfg.scan_max_deg);
    setInt("az_count", new_cfg.beam_count);
    setDouble("boshu", new_cfg.beam_width_deg);
    setDouble("loc_beam_gate_deg", new_cfg.beam_width_deg * 2.0);
    setDouble("max_theta", std::max(std::fabs(new_cfg.scan_min_deg), std::fabs(new_cfg.scan_max_deg)));
    setInt("period_first", 1);
    setInt("period_num", new_cfg.beam_count);
    setInt("rg_ed", new_cfg.pc_crop_len - 1);

    if (!doc.SaveFile(path.c_str())) {
        err = "failed to save generated XML: " + path;
        return false;
    }
    return true;
}

} // namespace sim_stage1
} // namespace gmti
