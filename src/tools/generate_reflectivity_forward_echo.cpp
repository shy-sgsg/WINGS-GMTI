#include "dbs/lodepng.h"
#include "geo/geoProj.hpp"
#include "tinyxml.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kC = 299792458.0;

struct Config {
    std::string pos_file;
    std::string output_file;
    int week = 0;
    int sec_bias = 0;
    int target_beam_count = 61;
    int pulse_num = 130;
    int ddc_len = 11820;
    double scan_min_deg = -59.0;
    double scan_max_deg = 61.0;
    double prf_hz = 1300.0;
    double fc_hz = 17.0e9;
    double fs_hz = 225.0e6;
    double r_min_m = 13950.0;
    double l0_deg = 117.0;
    double d_channel_m = 0.17;
    double beamwidth_deg = 2.0;
    double force_alt_m = 6000.0;
    double force_speed_mps = 60.0;
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
    double e = 0.0;
    double n = 0.0;
};

struct TargetMotion {
    double utc0 = 0.0;
    double lat0_deg = 0.0;
    double lon0_deg = 0.0;
    double vn_mps = 0.0;
    double ve_mps = 0.0;
    double alt_m = 0.0;
};

struct Bounds {
    double min_x = 0.0;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
};

struct Scatterer {
    double x = 0.0;
    double y = 0.0;
    double amp = 0.0;
};

std::string getText(TiXmlElement *parent, const char *name) {
    TiXmlElement *elem = parent ? parent->FirstChildElement(name) : nullptr;
    return (elem && elem->GetText()) ? elem->GetText() : "";
}

double parseDouble(TiXmlElement *parent, const char *name, double fallback) {
    const std::string text = getText(parent, name);
    return text.empty() ? fallback : std::atof(text.c_str());
}

int parseInt(TiXmlElement *parent, const char *name, int fallback) {
    const std::string text = getText(parent, name);
    return text.empty() ? fallback : std::atoi(text.c_str());
}

bool loadConfig(const std::string &xml_path, Config &cfg) {
    TiXmlDocument doc;
    if (!doc.LoadFile(xml_path.c_str())) {
        std::cerr << "Load XML failed: " << xml_path << std::endl;
        return false;
    }
    TiXmlElement *root = doc.FirstChildElement("GMTI");
    TiXmlElement *param = root ? root->FirstChildElement("GMTI_parameter") : nullptr;
    if (!param) {
        std::cerr << "XML must contain <GMTI><GMTI_parameter>" << std::endl;
        return false;
    }

    cfg.pos_file = getText(param, "Plane_POS");
    cfg.output_file = getText(param, "GMTI_data_new");
    cfg.week = parseInt(param, "week_offset", cfg.week);
    cfg.sec_bias = parseInt(param, "secBias", cfg.sec_bias);
    cfg.target_beam_count =
        parseInt(param, "period_num",
                 parseInt(param, "az_count", cfg.target_beam_count));
    cfg.pulse_num = parseInt(param, "pulse_num", cfg.pulse_num);
    cfg.ddc_len = parseInt(param, "pulse_len", cfg.ddc_len);
    cfg.scan_min_deg = parseDouble(param, "scan_min_deg", cfg.scan_min_deg);
    cfg.scan_max_deg = parseDouble(param, "scan_max_deg", cfg.scan_max_deg);
    cfg.prf_hz = parseDouble(param, "PRF", cfg.prf_hz);
    cfg.fc_hz = parseDouble(param, "fc", cfg.fc_hz / 1.0e9) * 1.0e9;
    cfg.fs_hz = parseDouble(param, "fs", cfg.fs_hz / 1.0e6) * 1.0e6;
    cfg.r_min_m = parseDouble(param, "Rmin", cfg.r_min_m);
    cfg.l0_deg = parseDouble(param, "ref_lon", cfg.l0_deg);
    cfg.d_channel_m = parseDouble(param, "d_chan", cfg.d_channel_m);
    cfg.beamwidth_deg = parseDouble(param, "boshu", cfg.beamwidth_deg);

    if (cfg.pos_file.empty() || cfg.output_file.empty() ||
        cfg.target_beam_count <= 0 || cfg.pulse_num <= 0 ||
        cfg.ddc_len <= 0 || cfg.prf_hz <= 0.0 || cfg.fc_hz <= 0.0 ||
        cfg.fs_hz <= 0.0) {
        std::cerr << "Invalid forward echo config" << std::endl;
        return false;
    }
    return true;
}

bool readPosRows(const std::string &path, std::vector<PosRow> &rows) {
    std::ifstream fin(path, std::ios::binary);
    if (!fin) {
        std::cerr << "Unable to open POS file: " << path << std::endl;
        return false;
    }
    fin.seekg(0, std::ios::end);
    const std::streamsize bytes = fin.tellg();
    fin.seekg(0, std::ios::beg);
    const std::streamsize row_bytes = 7 * static_cast<std::streamsize>(sizeof(double));
    if (bytes <= 0 || (bytes % row_bytes) != 0) {
        std::cerr << "POS file size mismatch" << std::endl;
        return false;
    }
    rows.resize(static_cast<size_t>(bytes / row_bytes));
    for (size_t i = 0; i < rows.size(); ++i) {
        double v[7] = {0.0};
        fin.read(reinterpret_cast<char *>(v), sizeof(v));
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

double radToDegIfLikely(double v) {
    return std::fabs(v) <= 3.2 ? v * 180.0 / M_PI : v;
}

PosSample interpPos(const std::vector<PosRow> &rows, double utc, double l0_deg) {
    PosSample out;
    if (rows.empty()) {
        return out;
    }
    const auto copy_row = [&](const PosRow &r) {
        out.lat_deg = radToDegIfLikely(r.lat);
        out.lon_deg = radToDegIfLikely(r.lon);
        out.alt_m = r.alt;
        out.vn = r.vn;
        out.ve = r.ve;
        out.vd = r.vd;
        Gaussp3(out.lat_deg, out.lon_deg, l0_deg, out.e, out.n);
    };
    if (utc <= rows.front().t) {
        copy_row(rows.front());
        return out;
    }
    if (utc >= rows.back().t) {
        copy_row(rows.back());
        return out;
    }
    const auto it = std::lower_bound(
        rows.begin(), rows.end(), utc,
        [](const PosRow &r, double value) { return r.t < value; });
    const PosRow &b = *it;
    const PosRow &a = *(it - 1);
    const double ratio = (b.t == a.t) ? 0.0 : (utc - a.t) / (b.t - a.t);
    const auto lerp = [ratio](double x, double y) {
        return x + ratio * (y - x);
    };
    out.lat_deg = radToDegIfLikely(lerp(a.lat, b.lat));
    out.lon_deg = radToDegIfLikely(lerp(a.lon, b.lon));
    out.alt_m = lerp(a.alt, b.alt);
    out.vn = lerp(a.vn, b.vn);
    out.ve = lerp(a.ve, b.ve);
    out.vd = lerp(a.vd, b.vd);
    Gaussp3(out.lat_deg, out.lon_deg, l0_deg, out.e, out.n);
    return out;
}

TargetMotion buildMotion(const std::vector<PosRow> &rows,
                         double utc0,
                         const Config &cfg) {
    const double earth_radius_m = 6378137.0;
    const PosSample origin = interpPos(rows, utc0, cfg.l0_deg);
    const PosSample before = interpPos(rows, utc0 - 0.5, cfg.l0_deg);
    const PosSample after = interpPos(rows, utc0 + 0.5, cfg.l0_deg);
    const double mean_lat_rad =
        0.5 * (before.lat_deg + after.lat_deg) * M_PI / 180.0;
    const double north_delta =
        (after.lat_deg - before.lat_deg) * M_PI / 180.0 * earth_radius_m;
    const double east_delta =
        (after.lon_deg - before.lon_deg) * M_PI / 180.0 *
        earth_radius_m * std::cos(mean_lat_rad);
    const double norm = std::hypot(north_delta, east_delta);

    TargetMotion motion;
    motion.utc0 = utc0;
    motion.lat0_deg = origin.lat_deg;
    motion.lon0_deg = origin.lon_deg;
    motion.alt_m = cfg.force_alt_m;
    if (norm > 1e-6) {
        motion.vn_mps = cfg.force_speed_mps * north_delta / norm;
        motion.ve_mps = cfg.force_speed_mps * east_delta / norm;
    } else {
        motion.vn_mps = 0.0;
        motion.ve_mps = cfg.force_speed_mps;
    }
    return motion;
}

PosSample evalMotion(const TargetMotion &motion, double utc, double l0_deg) {
    const double earth_radius_m = 6378137.0;
    const double dt = utc - motion.utc0;
    PosSample pos;
    pos.lat_deg = motion.lat0_deg +
        (motion.vn_mps * dt / earth_radius_m) * 180.0 / M_PI;
    const double cos_lat =
        std::max(1e-9, std::fabs(std::cos(motion.lat0_deg * M_PI / 180.0)));
    pos.lon_deg = motion.lon0_deg +
        (motion.ve_mps * dt / (earth_radius_m * cos_lat)) * 180.0 / M_PI;
    pos.alt_m = motion.alt_m;
    pos.vn = motion.vn_mps;
    pos.ve = motion.ve_mps;
    pos.vd = 0.0;
    Gaussp3(pos.lat_deg, pos.lon_deg, l0_deg, pos.e, pos.n);
    return pos;
}

bool parseCornerFile(const std::string &path, double l0_deg, Bounds &bounds) {
    std::ifstream fin(path);
    if (!fin) {
        std::cerr << "Unable to open corner file: " << path << std::endl;
        return false;
    }
    std::map<std::string, double> vals;
    std::string key, eq;
    double value = 0.0;
    while (fin >> key >> eq >> value) {
        if (!key.empty()) {
            vals[key] = value;
        }
    }
    const char *required[] = {"B0", "B1", "B2", "B3", "L0", "L1", "L2", "L3"};
    for (const char *name : required) {
        if (vals.find(name) == vals.end()) {
            std::cerr << "Corner file missing " << name << std::endl;
            return false;
        }
    }
    double xs[4] = {0.0}, ys[4] = {0.0};
    for (int i = 0; i < 4; ++i) {
        const std::string bi = "B" + std::to_string(i);
        const std::string li = "L" + std::to_string(i);
        Gaussp3(vals[bi], vals[li], l0_deg, xs[i], ys[i]);
    }
    bounds.min_x = *std::min_element(xs, xs + 4);
    bounds.max_x = *std::max_element(xs, xs + 4);
    bounds.min_y = *std::min_element(ys, ys + 4);
    bounds.max_y = *std::max_element(ys, ys + 4);
    return bounds.max_x > bounds.min_x && bounds.max_y > bounds.min_y;
}

bool loadScatterers(const std::string &png_path,
                    const Bounds &bounds,
                    int stride,
                    size_t max_scatterers,
                    std::vector<Scatterer> &scatterers) {
    unsigned char *raw = nullptr;
    unsigned w = 0, h = 0;
    const unsigned err =
        lodepng_decode_file(&raw, &w, &h, png_path.c_str(), LCT_GREY, 8);
    if (err != 0 || !raw || w == 0 || h == 0) {
        std::cerr << "PNG decode failed: " << png_path
                  << " err=" << err << std::endl;
        std::free(raw);
        return false;
    }
    stride = std::max(1, stride);
    std::vector<Scatterer> all;
    for (unsigned r = 0; r < h; r += static_cast<unsigned>(stride)) {
        const double x =
            bounds.min_x +
            (h <= 1 ? 0.0 :
             (static_cast<double>(r) / static_cast<double>(h - 1)) *
             (bounds.max_x - bounds.min_x));
        for (unsigned c = 0; c < w; c += static_cast<unsigned>(stride)) {
            const uint8_t pix = raw[static_cast<size_t>(r) * w + c];
            if (pix == 0) {
                continue;
            }
            const double y =
                bounds.min_y +
                (w <= 1 ? 0.0 :
                 (static_cast<double>(w - 1 - c) / static_cast<double>(w - 1)) *
                 (bounds.max_y - bounds.min_y));
            Scatterer s;
            s.x = x;
            s.y = y;
            const double g = static_cast<double>(pix) / 255.0;
            s.amp = g * g;
            all.push_back(s);
        }
    }
    std::free(raw);

    if (max_scatterers > 0 && all.size() > max_scatterers) {
        std::nth_element(
            all.begin(), all.begin() + static_cast<std::ptrdiff_t>(max_scatterers),
            all.end(),
            [](const Scatterer &a, const Scatterer &b) {
                return a.amp > b.amp;
            });
        all.resize(max_scatterers);
    }
    scatterers.swap(all);
    std::cout << "scatterers: " << scatterers.size()
              << " from " << png_path << std::endl;
    return !scatterers.empty();
}

int16_t sat_i16(double x) {
    if (x > 32767.0) return 32767;
    if (x < -32768.0) return -32768;
    return static_cast<int16_t>(std::lround(x));
}

void write_u16(std::vector<uint8_t> &buf, size_t off, uint16_t v) {
    buf[off] = static_cast<uint8_t>(v & 0xff);
    buf[off + 1] = static_cast<uint8_t>((v >> 8) & 0xff);
}

void write_i16(std::vector<uint8_t> &buf, size_t off, int16_t v) {
    write_u16(buf, off, static_cast<uint16_t>(v));
}

void write_u32(std::vector<uint8_t> &buf, size_t off, uint32_t v) {
    for (int i = 0; i < 4; ++i)
        buf[off + static_cast<size_t>(i)] =
            static_cast<uint8_t>((v >> (8 * i)) & 0xff);
}

void write_u64(std::vector<uint8_t> &buf, size_t off, uint64_t v) {
    for (int i = 0; i < 8; ++i)
        buf[off + static_cast<size_t>(i)] =
            static_cast<uint8_t>((v >> (8 * i)) & 0xff);
}

void write_f32(std::vector<uint8_t> &buf, size_t off, float v) {
    uint32_t raw = 0;
    std::memcpy(&raw, &v, sizeof(raw));
    write_u32(buf, off, raw);
}

void write_f64(std::vector<uint8_t> &buf, size_t off, double v) {
    uint64_t raw = 0;
    std::memcpy(&raw, &v, sizeof(raw));
    write_u64(buf, off, raw);
}

void fillHeader(std::vector<uint8_t> &pkt,
                uint8_t mode,
                uint32_t prt_counter,
                double utc,
                int week,
                const PosSample &pos,
                double angle_deg,
                uint32_t prt_len) {
    pkt.assign(prt_len, 0);
    write_u64(pkt, 0, 0x5A5A5A5A5A5A5A5AULL);
    pkt[8] = mode;
    write_u32(pkt, 9, prt_len);
    write_f32(pkt, 16, static_cast<float>(utc));
    write_u32(pkt, 20, prt_counter);
    pkt[88] = 0x02;
    pkt[90] = 0x40;
    pkt[92] = 0x0B;
    pkt[93] = 0x01;
    write_i16(pkt, 94, static_cast<int16_t>(week));
    const double sec_day = utc - std::floor(utc / 86400.0) * 86400.0;
    write_u32(pkt, 96, static_cast<uint32_t>(std::llround(sec_day * 1000.0)));
    write_f64(pkt, 104, pos.lat_deg);
    write_f64(pkt, 112, pos.lon_deg);
    write_f64(pkt, 120, pos.alt_m);
    write_f32(pkt, 128, static_cast<float>(pos.vn));
    write_f32(pkt, 132, static_cast<float>(pos.ve));
    write_f32(pkt, 136, static_cast<float>(pos.vd));
    write_f32(pkt, 140, static_cast<float>(std::hypot(pos.vn, pos.ve)));
    pkt[208] = static_cast<uint8_t>(prt_counter & 0xffU);
    write_i16(pkt, 218, sat_i16(angle_deg * 100.0));
    write_u64(pkt, 248, 0x5B5B5B5B5B5B5B5BULL);
}

double wrap180(double angle) {
    angle = std::fmod(angle + 180.0, 360.0);
    if (angle < 0.0) angle += 360.0;
    return angle - 180.0;
}

void synthPulse(const Config &cfg,
                const std::vector<Scatterer> &scatterers,
                const PosSample &pos,
                double beam_angle_deg,
                std::vector<std::complex<float>> &ch1,
                std::vector<std::complex<float>> &ch2) {
    ch1.assign(static_cast<size_t>(cfg.ddc_len), {0.0f, 0.0f});
    ch2.assign(static_cast<size_t>(cfg.ddc_len), {0.0f, 0.0f});
    const double lambda = kC / cfg.fc_hz;
    const double range_bin = kC / (2.0 * cfg.fs_hz);
    const double heading_deg = std::atan2(pos.vn, pos.ve) * 180.0 / M_PI;
    const double heading_rad = heading_deg * M_PI / 180.0;
    const double ue = std::cos(heading_rad);
    const double un = std::sin(heading_rad);
    const double sigma_deg = std::max(0.25, cfg.beamwidth_deg * 0.5);
    const double alt = pos.alt_m;

    for (const Scatterer &s : scatterers) {
        const double de = s.x - pos.e;
        const double dn = s.y - pos.n;
        const double bearing_deg = std::atan2(dn, de) * 180.0 / M_PI;
        const double squint_deg = wrap180(bearing_deg - heading_deg + 90.0);
        const double err_deg = wrap180(squint_deg - beam_angle_deg);
        if (std::fabs(err_deg) > cfg.beamwidth_deg * 1.5) {
            continue;
        }
        const double beam_w =
            std::exp(-0.5 * (err_deg / sigma_deg) * (err_deg / sigma_deg));
        const double r1 = std::sqrt(de * de + dn * dn + alt * alt);
        const double de2 = de - cfg.d_channel_m * ue;
        const double dn2 = dn - cfg.d_channel_m * un;
        const double r2 = std::sqrt(de2 * de2 + dn2 * dn2 + alt * alt);
        const double idx = (r1 - cfg.r_min_m) / range_bin;
        const int i0 = static_cast<int>(std::floor(idx));
        const double frac = idx - static_cast<double>(i0);
        if (i0 < 0 || i0 + 1 >= cfg.ddc_len) {
            continue;
        }
        const double amp = s.amp * beam_w;
        const std::complex<double> z1 =
            std::polar(amp, -4.0 * M_PI * r1 / lambda);
        const std::complex<double> z2 =
            std::polar(amp, -4.0 * M_PI * r2 / lambda);
        const float w0 = static_cast<float>(1.0 - frac);
        const float w1 = static_cast<float>(frac);
        ch1[static_cast<size_t>(i0)] +=
            std::complex<float>(z1.real() * w0, z1.imag() * w0);
        ch1[static_cast<size_t>(i0 + 1)] +=
            std::complex<float>(z1.real() * w1, z1.imag() * w1);
        ch2[static_cast<size_t>(i0)] +=
            std::complex<float>(z2.real() * w0, z2.imag() * w0);
        ch2[static_cast<size_t>(i0 + 1)] +=
            std::complex<float>(z2.real() * w1, z2.imag() * w1);
    }
}

void fillPayload(std::vector<uint8_t> &pkt,
                 const std::vector<std::complex<float>> &ch1,
                 const std::vector<std::complex<float>> &ch2) {
    const size_t start = 256;
    for (size_t n = 0; n < ch1.size(); ++n) {
        const size_t off = start + n * 16U;
        write_f32(pkt, off + 0, ch1[n].real());
        write_f32(pkt, off + 4, ch1[n].imag());
        write_f32(pkt, off + 8, ch2[n].real());
        write_f32(pkt, off + 12, ch2[n].imag());
    }
}

bool generateForwardEcho(const Config &cfg,
                         const std::vector<Scatterer> &scatterers,
                         const std::vector<PosRow> &pos_rows,
                         const std::string &out_file,
                         int beam_limit,
                         int pulse_limit) {
    const uint32_t prt_len = 256U + static_cast<uint32_t>(cfg.ddc_len) * 16U;
    std::ofstream fout(out_file, std::ios::binary);
    if (!fout) {
        std::cerr << "Unable to open output file: " << out_file << std::endl;
        return false;
    }

    const double utc0 = pos_rows.empty() ? 0.0 : pos_rows.front().t;
    const TargetMotion motion = buildMotion(pos_rows, utc0, cfg);
    std::cout << "forward motion vn=" << motion.vn_mps
              << " ve=" << motion.ve_mps
              << " speed=" << std::hypot(motion.vn_mps, motion.ve_mps)
              << std::endl;

    std::vector<uint8_t> pkt(prt_len);
    std::vector<std::complex<float>> ch1, ch2;
    uint32_t counter = 0;
    const double step_deg =
        (cfg.target_beam_count <= 1) ? 0.0 :
        (cfg.scan_max_deg - cfg.scan_min_deg) /
            static_cast<double>(cfg.target_beam_count - 1);
    const int beam_count =
        (beam_limit > 0) ? std::min(cfg.target_beam_count, beam_limit)
                         : cfg.target_beam_count;
    const int pulse_count =
        (pulse_limit > 0) ? std::min(cfg.pulse_num, pulse_limit)
                          : cfg.pulse_num;

    for (int b = 0; b < beam_count; ++b) {
        const double angle = cfg.scan_min_deg + step_deg * b;
        for (int p = 0; p < pulse_count; ++p) {
            const double utc = utc0 + static_cast<double>(counter) / cfg.prf_hz;
            const PosSample pos = evalMotion(motion, utc, cfg.l0_deg);
            fillHeader(pkt, 5, counter, utc, cfg.week, pos, angle, prt_len);
            synthPulse(cfg, scatterers, pos, angle, ch1, ch2);
            fillPayload(pkt, ch1, ch2);
            fout.write(reinterpret_cast<const char *>(pkt.data()),
                       static_cast<std::streamsize>(pkt.size()));
            if (!fout) {
                std::cerr << "Write failed at beam=" << (b + 1)
                          << " pulse=" << p << std::endl;
                return false;
            }
            ++counter;
        }
        std::cout << "forward beam " << (b + 1) << " / "
                  << beam_count << " angle=" << angle
                  << " deg" << std::endl;
    }
    return true;
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <temp_config_61beam.xml> <reflectivity.png> <corners.txt>"
                     " [output_echo.dat] [scatter_stride=4] [max_scatterers=20000]"
                     " [beam_limit=0] [pulse_limit=0]"
                  << std::endl;
        return 1;
    }

    Config cfg;
    if (!loadConfig(argv[1], cfg)) {
        return 2;
    }
    const std::string png_path = argv[2];
    const std::string corner_path = argv[3];
    const std::string out_file = (argc > 4) ? argv[4] : cfg.output_file;
    const int scatter_stride = (argc > 5) ? std::atoi(argv[5]) : 4;
    const size_t max_scatterers =
        (argc > 6) ? static_cast<size_t>(std::strtoull(argv[6], nullptr, 10))
                   : 20000ULL;
    const int beam_limit = (argc > 7) ? std::atoi(argv[7]) : 0;
    const int pulse_limit = (argc > 8) ? std::atoi(argv[8]) : 0;

    Bounds bounds;
    if (!parseCornerFile(corner_path, cfg.l0_deg, bounds)) {
        return 3;
    }

    std::vector<Scatterer> scatterers;
    if (!loadScatterers(png_path, bounds, scatter_stride,
                        max_scatterers, scatterers)) {
        return 4;
    }

    std::vector<PosRow> pos_rows;
    if (!readPosRows(cfg.pos_file, pos_rows)) {
        return 5;
    }

    try {
        if (!generateForwardEcho(cfg, scatterers, pos_rows, out_file,
                                 beam_limit, pulse_limit)) {
            return 6;
        }
    } catch (const std::exception &e) {
        std::cerr << "Forward echo generation failed: " << e.what()
                  << std::endl;
        return 6;
    }
    return 0;
}
