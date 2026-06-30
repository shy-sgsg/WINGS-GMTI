#include "config_structs.hpp"
#include "rangeCompress.hpp"
#include "tinyxml.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
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

struct TargetProfile {
    int target_beam_count = 61;

    int src_beam_count_expected = 25;
    int src_pulse_num_expected = 256;
    int src_range_len_expected = 4096;
    double src_scan_start_deg = -23.0;
    double src_scan_step_deg = 2.0;

    int target_pulse_num = 130;
    int target_ddc_len = 11820;
    int target_fft_len = 12288;

    // These describe the downstream pulse-compression extraction contract.
    // This generator writes DDC samples and does not apply this crop itself.
    int range_crop_start = 3864;
    int range_crop_len = 4096;

    double prf_hz = 1300.0;
    double scan_start_deg = -59.0;
    double scan_step_deg = 2.0;

    double carrier_hz = 17.0e9;
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
};

struct TargetMotion {
    double utc0 = 0.0;
    double lat0_deg = 0.0;
    double lon0_deg = 0.0;
    double vn_mps = 0.0;
    double ve_mps = 0.0;
    double alt_m = 0.0;
};

std::string getText(TiXmlElement *parent, const char *name) {
    if (!parent) {
        return "";
    }
    TiXmlElement *elem = parent->FirstChildElement(name);
    return (elem && elem->GetText()) ? elem->GetText() : "";
}

bool loadConvertConfig(const std::string &xml_path, ConvertConfig &cfg) {
    TiXmlDocument doc;
    if (!doc.LoadFile(xml_path.c_str())) {
        std::cerr << "Load XML failed: " << xml_path << std::endl;
        return false;
    }

    TiXmlElement *root = doc.FirstChildElement("GMTI");
    TiXmlElement *param = root ? root->FirstChildElement("GMTI_parameter") : nullptr;
    if (!root || !param) {
        std::cerr << "XML must contain <GMTI><GMTI_parameter>" << std::endl;
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
        std::cerr << "XML missing required conversion fields" << std::endl;
        return false;
    }

    try {
        cfg.info_len = std::stoi(info_len);
        cfg.pulse_len = std::stoi(pulse_len);
        cfg.pulse_num = std::stoi(pulse_num);
        cfg.skip_az_num = skip.empty() ? 0 : std::stoi(skip);
        cfg.week = week.empty() ? 0 : std::stoi(week);
        cfg.sec_bias = sec_bias.empty() ? 0 : std::stoi(sec_bias);
    } catch (const std::exception &e) {
        std::cerr << "Invalid numeric XML field: " << e.what() << std::endl;
        return false;
    }
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
    const std::streamsize frame_bytes =
        pos_len * static_cast<std::streamsize>(sizeof(double));
    if (file_size <= 0 || (file_size % frame_bytes) != 0) {
        std::cerr << "POS file size mismatch with 7-double frame format" << std::endl;
        return false;
    }

    rows.resize(static_cast<size_t>(file_size / frame_bytes));
    for (size_t i = 0; i < rows.size(); ++i) {
        double v[pos_len] = {0.0};
        fin.read(reinterpret_cast<char *>(v), sizeof(v));
        if (!fin) {
            std::cerr << "POS read failed at frame " << i << std::endl;
            return false;
        }
        PosRow row;
        row.t = v[0];
        row.lat = v[1];
        row.lon = v[2];
        row.alt = v[3];
        row.vn = v[4];
        row.ve = v[5];
        row.vd = v[6];
        rows[i] = row;
    }
    return true;
}

double radToDegIfLikely(double x) {
    return std::fabs(x) <= 3.2 ? x * 180.0 / M_PI : x;
}

PosSample interpPosAtUtc(const std::vector<PosRow> &rows, double t) {
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
    };

    if (t <= rows.front().t) {
        copy_row(rows.front());
        return out;
    }
    if (t >= rows.back().t) {
        copy_row(rows.back());
        return out;
    }

    const auto it = std::lower_bound(
        rows.begin(), rows.end(), t,
        [](const PosRow &r, double value) { return r.t < value; });
    const PosRow &right = *it;
    const PosRow &left = *(it - 1);
    const double dt = right.t - left.t;
    const double ratio = dt == 0.0 ? 0.0 : (t - left.t) / dt;
    const auto lerp = [ratio](double a, double b) {
        return a + ratio * (b - a);
    };

    out.lat_deg = radToDegIfLikely(lerp(left.lat, right.lat));
    out.lon_deg = radToDegIfLikely(lerp(left.lon, right.lon));
    out.alt_m = lerp(left.alt, right.alt);
    out.vn = lerp(left.vn, right.vn);
    out.ve = lerp(left.ve, right.ve);
    out.vd = lerp(left.vd, right.vd);
    return out;
}

void forcePosProfile(PosSample &pos, double alt_m, double speed_mps) {
    pos.alt_m = alt_m;
    const double horizontal_norm =
        std::sqrt(pos.vn * pos.vn + pos.ve * pos.ve);
    if (horizontal_norm > 1e-6) {
        const double scale = speed_mps / horizontal_norm;
        pos.vn *= scale;
        pos.ve *= scale;
    } else {
        pos.vn = 0.0;
        pos.ve = speed_mps;
    }
    // Target profile represents level flight.
    pos.vd = 0.0;
}

TargetMotion buildTargetMotion(const std::vector<PosRow> &rows,
                               double utc0,
                               double alt_m,
                               double speed_mps) {
    const double earth_radius_m = 6378137.0;
    const double half_window_s = 0.5;
    const PosSample origin = interpPosAtUtc(rows, utc0);
    const PosSample before = interpPosAtUtc(rows, utc0 - half_window_s);
    const PosSample after = interpPosAtUtc(rows, utc0 + half_window_s);

    const double mean_lat_rad =
        0.5 * (before.lat_deg + after.lat_deg) * M_PI / 180.0;
    const double north_delta =
        (after.lat_deg - before.lat_deg) * M_PI / 180.0 * earth_radius_m;
    const double east_delta =
        (after.lon_deg - before.lon_deg) * M_PI / 180.0 *
        earth_radius_m * std::cos(mean_lat_rad);
    const double horizontal_delta =
        std::sqrt(north_delta * north_delta + east_delta * east_delta);

    TargetMotion motion;
    motion.utc0 = utc0;
    motion.lat0_deg = origin.lat_deg;
    motion.lon0_deg = origin.lon_deg;
    motion.alt_m = alt_m;
    if (horizontal_delta > 1e-6) {
        motion.vn_mps = speed_mps * north_delta / horizontal_delta;
        motion.ve_mps = speed_mps * east_delta / horizontal_delta;
    } else {
        motion.vn_mps = 0.0;
        motion.ve_mps = speed_mps;
    }
    return motion;
}

PosSample evaluateTargetMotion(const TargetMotion &motion, double utc) {
    const double earth_radius_m = 6378137.0;
    const double dt = utc - motion.utc0;
    PosSample pos;
    pos.lat_deg =
        motion.lat0_deg +
        (motion.vn_mps * dt / earth_radius_m) * 180.0 / M_PI;
    const double cos_lat =
        std::max(1e-9, std::fabs(std::cos(motion.lat0_deg * M_PI / 180.0)));
    pos.lon_deg =
        motion.lon0_deg +
        (motion.ve_mps * dt / (earth_radius_m * cos_lat)) * 180.0 / M_PI;
    pos.alt_m = motion.alt_m;
    pos.vn = motion.vn_mps;
    pos.ve = motion.ve_mps;
    pos.vd = 0.0;
    return pos;
}

uint8_t sat_u8(int x) {
    return static_cast<uint8_t>(std::max(0, std::min(255, x)));
}

int16_t sat_i16_from_double(double x) {
    if (x > 32767.0) {
        return 32767;
    }
    if (x < -32768.0) {
        return -32768;
    }
    return static_cast<int16_t>(std::lround(x));
}

void write_u16_le(std::vector<uint8_t> &buf, size_t off, uint16_t v) {
    buf[off] = static_cast<uint8_t>(v & 0xffU);
    buf[off + 1] = static_cast<uint8_t>((v >> 8) & 0xffU);
}

void write_i16_le(std::vector<uint8_t> &buf, size_t off, int16_t v) {
    write_u16_le(buf, off, static_cast<uint16_t>(v));
}

void write_u32_le(std::vector<uint8_t> &buf, size_t off, uint32_t v) {
    for (int i = 0; i < 4; ++i) {
        buf[off + static_cast<size_t>(i)] =
            static_cast<uint8_t>((v >> (8 * i)) & 0xffU);
    }
}

void write_u64_le(std::vector<uint8_t> &buf, size_t off, uint64_t v) {
    for (int i = 0; i < 8; ++i) {
        buf[off + static_cast<size_t>(i)] =
            static_cast<uint8_t>((v >> (8 * i)) & 0xffU);
    }
}

void write_f32_le(std::vector<uint8_t> &buf, size_t off, float v) {
    uint32_t raw = 0;
    std::memcpy(&raw, &v, sizeof(raw));
    write_u32_le(buf, off, raw);
}

void write_f64_le(std::vector<uint8_t> &buf, size_t off, double v) {
    uint64_t raw = 0;
    std::memcpy(&raw, &v, sizeof(raw));
    write_u64_le(buf, off, raw);
}

uint16_t read_u16_le(const std::vector<uint8_t> &buf, size_t off) {
    return static_cast<uint16_t>(
        static_cast<uint16_t>(buf[off]) |
        (static_cast<uint16_t>(buf[off + 1]) << 8));
}

uint32_t read_u32_le(const std::vector<uint8_t> &buf, size_t off) {
    uint32_t value = 0;
    for (int i = 0; i < 4; ++i) {
        value |= static_cast<uint32_t>(buf[off + static_cast<size_t>(i)])
                 << (8 * i);
    }
    return value;
}

float read_f32_le(const std::vector<uint8_t> &buf, size_t off) {
    const uint32_t raw = read_u32_le(buf, off);
    float value = 0.0f;
    std::memcpy(&value, &raw, sizeof(value));
    return value;
}

void fillHeaderAndIns(std::vector<uint8_t> &pkt,
                      uint8_t mode,
                      uint32_t prt_counter,
                      double utc_sec,
                      int week,
                      const PosSample &pos,
                      double fw_angle_deg,
                      uint32_t prt_len) {
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

    const double sec_of_day =
        utc_sec - std::floor(utc_sec / 86400.0) * 86400.0;
    write_u32_le(pkt, 96,
                 static_cast<uint32_t>(std::llround(sec_of_day * 1000.0)));
    write_f64_le(pkt, 104, pos.lat_deg);
    write_f64_le(pkt, 112, pos.lon_deg);
    write_f64_le(pkt, 120, pos.alt_m);
    write_f32_le(pkt, 128, static_cast<float>(pos.vn));
    write_f32_le(pkt, 132, static_cast<float>(pos.ve));
    write_f32_le(pkt, 136, static_cast<float>(pos.vd));
    const double speed =
        std::sqrt(pos.vn * pos.vn + pos.ve * pos.ve + pos.vd * pos.vd);
    write_f32_le(pkt, 140, static_cast<float>(speed));
    pkt[208] = static_cast<uint8_t>(prt_counter & 0xffU);
    write_i16_le(pkt, 218, sat_i16_from_double(fw_angle_deg * 100.0));
    write_u64_le(pkt, 248, 0x5B5B5B5B5B5B5B5BULL);
}

void fillGmtiIqPayloadOnePulse(
    std::vector<uint8_t> &pkt,
    const std::vector<std::complex<float>> &ch1,
    const std::vector<std::complex<float>> &ch2) {
    const int payload_start = 256;
    const int bytes_per_sample_2ch = 16;
    const int max_samples =
        (static_cast<int>(pkt.size()) - payload_start) / bytes_per_sample_2ch;
    if (static_cast<int>(ch1.size()) > max_samples ||
        static_cast<int>(ch2.size()) > max_samples) {
        throw std::runtime_error(
            "payload vector is larger than packet payload capacity");
    }

    const int n_samples =
        std::min(static_cast<int>(ch1.size()), static_cast<int>(ch2.size()));
    for (int n = 0; n < n_samples; ++n) {
        const size_t off =
            static_cast<size_t>(payload_start + n * bytes_per_sample_2ch);
        write_f32_le(pkt, off, ch1[static_cast<size_t>(n)].real());
        write_f32_le(pkt, off + 4, ch1[static_cast<size_t>(n)].imag());
        write_f32_le(pkt, off + 8, ch2[static_cast<size_t>(n)].real());
        write_f32_le(pkt, off + 12, ch2[static_cast<size_t>(n)].imag());
    }
}

uint32_t calcBeamCount(const std::string &path,
                       int info_len,
                       int pulse_len,
                       int pulse_num) {
    std::ifstream fin(path, std::ios::binary | std::ios::ate);
    if (!fin) {
        return 0;
    }
    const std::streamsize file_size = fin.tellg();
    const std::streamsize stride =
        static_cast<std::streamsize>(info_len) +
        static_cast<std::streamsize>(pulse_len) * 2 *
            static_cast<std::streamsize>(sizeof(float));
    const std::streamsize beam_bytes =
        stride * static_cast<std::streamsize>(pulse_num);
    if (beam_bytes <= 0 || file_size <= 0 || (file_size % beam_bytes) != 0) {
        return 0;
    }
    return static_cast<uint32_t>(file_size / beam_bytes);
}

struct BeamReplicaMapping {
    int source_beam = 0;
    double source_angle_deg = 0.0;
    double target_angle_deg = 0.0;
};

struct DebugSelection {
    std::vector<int> source_beams;
    std::vector<double> target_angles_deg;
};

struct PulsePowerStats {
    double sum_power_ch1 = 0.0;
    double sum_power_ch2 = 0.0;
    uint64_t sample_count = 0;
};

int mapTargetPulseToSourcePulse(int target_pulse,
                                int source_pulse_count,
                                int target_pulse_count);

int positiveModulo(int value, int modulus) {
    const int remainder = value % modulus;
    return remainder < 0 ? remainder + modulus : remainder;
}

std::string trimCopy(const std::string &text) {
    const std::string whitespace = " \t\r\n";
    const size_t first = text.find_first_not_of(whitespace);
    if (first == std::string::npos) {
        return "";
    }
    const size_t last = text.find_last_not_of(whitespace);
    return text.substr(first, last - first + 1);
}

bool parseIntCsv(const std::string &text,
                 std::vector<int> &values,
                 const char *name) {
    values.clear();
    const std::string trimmed = trimCopy(text);
    if (trimmed.empty() || trimmed == "all" || trimmed == "*") {
        return true;
    }
    std::stringstream ss(trimmed);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = trimCopy(token);
        if (token.empty()) {
            continue;
        }
        try {
            values.push_back(std::stoi(token));
        } catch (const std::exception &e) {
            std::cerr << "Invalid " << name << " list item '" << token
                      << "': " << e.what() << std::endl;
            return false;
        }
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return true;
}

bool parseDoubleCsv(const std::string &text,
                    std::vector<double> &values,
                    const char *name) {
    values.clear();
    const std::string trimmed = trimCopy(text);
    if (trimmed.empty() || trimmed == "all" || trimmed == "*") {
        return true;
    }
    std::stringstream ss(trimmed);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = trimCopy(token);
        if (token.empty()) {
            continue;
        }
        try {
            values.push_back(std::stod(token));
        } catch (const std::exception &e) {
            std::cerr << "Invalid " << name << " list item '" << token
                      << "': " << e.what() << std::endl;
            return false;
        }
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return true;
}

bool containsInt(const std::vector<int> &values, int value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

bool containsAngle(const std::vector<double> &angles,
                   double angle_deg,
                   double tolerance_deg) {
    for (double selected : angles) {
        if (std::fabs(selected - angle_deg) <= tolerance_deg) {
            return true;
        }
    }
    return false;
}

std::string joinInts(const std::vector<int> &values) {
    if (values.empty()) {
        return "all";
    }
    std::ostringstream os;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            os << ',';
        }
        os << values[i];
    }
    return os.str();
}

std::string joinDoubles(const std::vector<double> &values) {
    if (values.empty()) {
        return "all";
    }
    std::ostringstream os;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            os << ',';
        }
        os << values[i];
    }
    return os.str();
}

BeamReplicaMapping mapTargetBeamToReplicatedSource(
    int target_beam,
    int source_start_beam,
    int source_beam_count,
    const TargetProfile &target) {
    if (source_beam_count <= 0) {
        throw std::invalid_argument("source beam count must be positive");
    }

    const double target_angle_deg =
        target.scan_start_deg + target.scan_step_deg * target_beam;
    const int target_angle_slot = static_cast<int>(std::lround(
        (target_angle_deg - target.src_scan_start_deg) /
        target.src_scan_step_deg));
    const int source_offset =
        positiveModulo(target_angle_slot, source_beam_count);

    BeamReplicaMapping mapping;
    mapping.source_beam = source_start_beam + source_offset;
    mapping.source_angle_deg =
        target.src_scan_start_deg +
        target.src_scan_step_deg * source_offset;
    mapping.target_angle_deg = target_angle_deg;
    return mapping;
}

double theoreticalDopplerHz(double angle_deg,
                            double speed_mps,
                            double carrier_hz) {
    const double speed_of_light_mps = 299792458.0;
    const double wavelength_m = speed_of_light_mps / carrier_hz;
    return -2.0 * speed_mps / wavelength_m *
           std::sin(angle_deg * M_PI / 180.0);
}

double dopplerShiftForReplicaHz(const BeamReplicaMapping &mapping,
                                const TargetProfile &target) {
    return theoreticalDopplerHz(
               mapping.target_angle_deg,
               target.force_speed_mps,
               target.carrier_hz) -
           theoreticalDopplerHz(
               mapping.source_angle_deg,
               target.force_speed_mps,
               target.carrier_hz);
}

double wrapFrequencyHz(double frequency_hz, double prf_hz) {
    if (!(prf_hz > 0.0)) {
        return frequency_hz;
    }
    double wrapped = std::fmod(frequency_hz + 0.5 * prf_hz, prf_hz);
    if (wrapped < 0.0) {
        wrapped += prf_hz;
    }
    return wrapped - 0.5 * prf_hz;
}

double estimateMappedPulseDopplerHz(
    const std::vector<std::complex<float>> &data,
    int source_pulse_count,
    int source_range_len,
    int target_pulse_count,
    double target_prf_hz) {
    if (source_pulse_count <= 1 || source_range_len <= 0 ||
        target_pulse_count <= 1 || !(target_prf_hz > 0.0)) {
        return 0.0;
    }

    std::complex<double> corr(0.0, 0.0);
    const int range_step = std::max(1, source_range_len / 512);
    for (int target_pulse = 1;
         target_pulse < target_pulse_count;
         ++target_pulse) {
        const int src_pulse = mapTargetPulseToSourcePulse(
            target_pulse, source_pulse_count, target_pulse_count);
        const int prev_src_pulse = mapTargetPulseToSourcePulse(
            target_pulse - 1, source_pulse_count, target_pulse_count);
        const size_t src_offset =
            static_cast<size_t>(src_pulse) *
            static_cast<size_t>(source_range_len);
        const size_t prev_offset =
            static_cast<size_t>(prev_src_pulse) *
            static_cast<size_t>(source_range_len);
        for (int n = 0; n < source_range_len; n += range_step) {
            corr += static_cast<std::complex<double>>(
                        data[src_offset + static_cast<size_t>(n)]) *
                    std::conj(static_cast<std::complex<double>>(
                        data[prev_offset + static_cast<size_t>(n)]));
        }
    }
    if (std::abs(corr) <= 0.0) {
        return 0.0;
    }
    return target_prf_hz / (2.0 * M_PI) * std::arg(corr);
}

double dopplerShiftForReplicaFromMeasuredSource(
    const BeamReplicaMapping &mapping,
    const TargetProfile &target,
    double measured_source_wrapped_hz) {
    if (std::fabs(mapping.target_angle_deg - mapping.source_angle_deg) < 1e-9) {
        return 0.0;
    }
    const double desired_wrapped_hz = wrapFrequencyHz(
        theoreticalDopplerHz(mapping.target_angle_deg,
                             target.force_speed_mps,
                             target.carrier_hz),
        target.prf_hz);
    return wrapFrequencyHz(
        desired_wrapped_hz - measured_source_wrapped_hz,
        target.prf_hz);
}

void applySlowTimeFrequencyShift(
    std::vector<std::complex<float>> &samples,
    int pulse_index,
    double frequency_shift_hz,
    double prf_hz) {
    if (frequency_shift_hz == 0.0) {
        return;
    }
    const double phase =
        2.0 * M_PI * frequency_shift_hz *
        static_cast<double>(pulse_index) / prf_hz;
    const std::complex<float> rotator(
        static_cast<float>(std::cos(phase)),
        static_cast<float>(std::sin(phase)));
    for (std::complex<float> &sample : samples) {
        sample *= rotator;
    }
}

void accumulatePowerStats(const std::vector<std::complex<float>> &ch1,
                          const std::vector<std::complex<float>> &ch2,
                          PulsePowerStats &stats) {
    const size_t n = std::min(ch1.size(), ch2.size());
    for (size_t i = 0; i < n; ++i) {
        stats.sum_power_ch1 += std::norm(ch1[i]);
        stats.sum_power_ch2 += std::norm(ch2[i]);
    }
    stats.sample_count += static_cast<uint64_t>(n);
}

double meanPower(const PulsePowerStats &stats, int channel) {
    if (stats.sample_count == 0) {
        return 0.0;
    }
    const double sum =
        channel == 1 ? stats.sum_power_ch1 : stats.sum_power_ch2;
    return sum / static_cast<double>(stats.sample_count);
}

int mapTargetPulseToSourcePulse(int target_pulse,
                               int source_pulse_count,
                               int target_pulse_count) {
    if (source_pulse_count <= 1 || target_pulse_count <= 1) {
        return 0;
    }

    // Uniform time-domain extraction over the complete source aperture.
    // For 256 -> 130 this maps target pulse 0 to source pulse 0 and target
    // pulse 129 to source pulse 255. Intermediate indices are rounded to the
    // nearest source pulse, so the extraction interval alternates between
    // one and two source PRTs instead of discarding both aperture edges.
    const int source_pulse = static_cast<int>(std::lround(
        target_pulse * static_cast<double>(source_pulse_count - 1) /
        static_cast<double>(target_pulse_count - 1)));
    return std::max(0, std::min(source_pulse_count - 1, source_pulse));
}

class RangeFftResampler {
public:
    RangeFftResampler(int src_len, int fft_len, int out_len)
        : src_len_(src_len),
          fft_len_(fft_len),
          out_len_(out_len),
          spec_shift_(static_cast<size_t>(src_len)),
          spec_pad_shift_(static_cast<size_t>(fft_len)),
          spec_pad_(static_cast<size_t>(fft_len)) {
        if (src_len_ <= 0 || fft_len_ < src_len_ || out_len_ > fft_len_ ||
            out_len_ <= 0 || (src_len_ % 2) != 0 || (fft_len_ % 2) != 0) {
            throw std::invalid_argument("invalid FFT resampler lengths");
        }

        in_ = static_cast<fftwf_complex *>(
            fftwf_malloc(sizeof(fftwf_complex) * static_cast<size_t>(src_len_)));
        spec_src_ = static_cast<fftwf_complex *>(
            fftwf_malloc(sizeof(fftwf_complex) * static_cast<size_t>(src_len_)));
        inv_in_ = static_cast<fftwf_complex *>(
            fftwf_malloc(sizeof(fftwf_complex) * static_cast<size_t>(fft_len_)));
        time_full_ = static_cast<fftwf_complex *>(
            fftwf_malloc(sizeof(fftwf_complex) * static_cast<size_t>(fft_len_)));
        if (!in_ || !spec_src_ || !inv_in_ || !time_full_) {
            cleanup();
            throw std::bad_alloc();
        }

        forward_plan_ =
            fftwf_plan_dft_1d(src_len_, in_, spec_src_, FFTW_FORWARD, FFTW_ESTIMATE);
        inverse_plan_ =
            fftwf_plan_dft_1d(fft_len_, inv_in_, time_full_, FFTW_BACKWARD, FFTW_ESTIMATE);
        if (!forward_plan_ || !inverse_plan_) {
            cleanup();
            throw std::runtime_error("failed to create FFTW plans");
        }
    }

    ~RangeFftResampler() {
        cleanup();
    }

    RangeFftResampler(const RangeFftResampler &) = delete;
    RangeFftResampler &operator=(const RangeFftResampler &) = delete;

    bool process(const std::complex<float> *src,
                 std::vector<std::complex<float>> &out) {
        if (!src || !forward_plan_ || !inverse_plan_) {
            return false;
        }
        for (int n = 0; n < src_len_; ++n) {
            in_[n][0] = src[n].real();
            in_[n][1] = src[n].imag();
        }
        fftwf_execute(forward_plan_);

        for (int k = 0; k < src_len_; ++k) {
            const int src_index = (k + src_len_ / 2) % src_len_;
            spec_shift_[static_cast<size_t>(k)] =
                std::complex<float>(spec_src_[src_index][0],
                                    spec_src_[src_index][1]);
        }

        std::fill(spec_pad_shift_.begin(), spec_pad_shift_.end(),
                  std::complex<float>(0.0f, 0.0f));
        const int insert = (fft_len_ - src_len_) / 2;
        std::copy(spec_shift_.begin(), spec_shift_.end(),
                  spec_pad_shift_.begin() + insert);

        for (int k = 0; k < fft_len_; ++k) {
            spec_pad_[static_cast<size_t>(k)] =
                spec_pad_shift_[static_cast<size_t>((k + fft_len_ / 2) % fft_len_)];
            inv_in_[k][0] = spec_pad_[static_cast<size_t>(k)].real();
            inv_in_[k][1] = spec_pad_[static_cast<size_t>(k)].imag();
        }
        fftwf_execute(inverse_plan_);

        out.resize(static_cast<size_t>(out_len_));
        // FFTW leaves both transforms unnormalized. Dividing the padded IFFT
        // by src_len (not fft_len) preserves the original sample amplitude;
        // equivalently this is a conventional 1/fft_len IFFT multiplied by
        // fft_len/src_len for interpolation gain compensation.
        const float scale = 1.0f / static_cast<float>(src_len_);
        for (int n = 0; n < out_len_; ++n) {
            out[static_cast<size_t>(n)] =
                std::complex<float>(time_full_[n][0] * scale,
                                    time_full_[n][1] * scale);
        }
        return true;
    }

private:
    void cleanup() {
        if (forward_plan_) {
            fftwf_destroy_plan(forward_plan_);
            forward_plan_ = nullptr;
        }
        if (inverse_plan_) {
            fftwf_destroy_plan(inverse_plan_);
            inverse_plan_ = nullptr;
        }
        fftwf_free(in_);
        fftwf_free(spec_src_);
        fftwf_free(inv_in_);
        fftwf_free(time_full_);
        in_ = nullptr;
        spec_src_ = nullptr;
        inv_in_ = nullptr;
        time_full_ = nullptr;
    }

    int src_len_;
    int fft_len_;
    int out_len_;
    fftwf_complex *in_ = nullptr;
    fftwf_complex *spec_src_ = nullptr;
    fftwf_complex *inv_in_ = nullptr;
    fftwf_complex *time_full_ = nullptr;
    fftwf_plan forward_plan_ = nullptr;
    fftwf_plan inverse_plan_ = nullptr;
    std::vector<std::complex<float>> spec_shift_;
    std::vector<std::complex<float>> spec_pad_shift_;
    std::vector<std::complex<float>> spec_pad_;
};

bool runSelfTest() {
    const TargetProfile target;
    const uint32_t prt_len =
        256U + static_cast<uint32_t>(target.target_ddc_len) * 16U;
    const uint64_t packets =
        static_cast<uint64_t>(target.target_beam_count) *
        static_cast<uint64_t>(target.target_pulse_num);
    if (prt_len != 189376U || packets != 7930U ||
        packets * static_cast<uint64_t>(prt_len) != 1501751680ULL) {
        std::cerr << "[SELFTEST] packet geometry failed" << std::endl;
        return false;
    }
    const BeamReplicaMapping left =
        mapTargetBeamToReplicatedSource(0, 1, 25, target);
    const BeamReplicaMapping center_first =
        mapTargetBeamToReplicatedSource(18, 1, 25, target);
    const BeamReplicaMapping center_last =
        mapTargetBeamToReplicatedSource(42, 1, 25, target);
    const BeamReplicaMapping right =
        mapTargetBeamToReplicatedSource(60, 1, 25, target);
    if (left.source_beam != 8 || left.source_angle_deg != -9.0 ||
        left.target_angle_deg != -59.0 ||
        center_first.source_beam != 1 ||
        center_first.source_angle_deg != -23.0 ||
        center_first.target_angle_deg != -23.0 ||
        center_last.source_beam != 25 ||
        center_last.source_angle_deg != 25.0 ||
        center_last.target_angle_deg != 25.0 ||
        right.source_beam != 18 || right.source_angle_deg != 11.0 ||
        right.target_angle_deg != 61.0 ||
        dopplerShiftForReplicaHz(center_first, target) != 0.0 ||
        dopplerShiftForReplicaHz(center_last, target) != 0.0) {
        std::cerr << "[SELFTEST] replicated beam mapping failed" << std::endl;
        return false;
    }
    const double test_shift_hz = 137.5;
    std::complex<double> adjacent_correlation(0.0, 0.0);
    std::complex<float> previous(1.0f, 0.0f);
    for (int pulse = 0; pulse < 16; ++pulse) {
        std::vector<std::complex<float>> shifted(1, {1.0f, 0.0f});
        applySlowTimeFrequencyShift(
            shifted, pulse, test_shift_hz, target.prf_hz);
        if (pulse > 0) {
            adjacent_correlation +=
                static_cast<std::complex<double>>(shifted[0]) *
                std::conj(static_cast<std::complex<double>>(previous));
        }
        previous = shifted[0];
    }
    const double estimated_shift_hz =
        target.prf_hz / (2.0 * M_PI) * std::arg(adjacent_correlation);
    if (std::fabs(estimated_shift_hz - test_shift_hz) > 1e-4) {
        std::cerr << "[SELFTEST] slow-time Doppler shift sign failed"
                  << std::endl;
        return false;
    }
    if (wrapFrequencyHz(5832.77, target.prf_hz) < 632.0 ||
        wrapFrequencyHz(5832.77, target.prf_hz) > 633.5 ||
        dopplerShiftForReplicaFromMeasuredSource(
            left, target,
            wrapFrequencyHz(theoreticalDopplerHz(
                left.source_angle_deg,
                target.force_speed_mps,
                target.carrier_hz),
                target.prf_hz)) == 0.0) {
        std::cerr << "[SELFTEST] measured-source Doppler shift failed"
                  << std::endl;
        return false;
    }
    int previous_source_pulse = -1;
    for (int target_pulse = 0;
         target_pulse < target.target_pulse_num;
         ++target_pulse) {
        const int source_pulse = mapTargetPulseToSourcePulse(
            target_pulse, target.src_pulse_num_expected,
            target.target_pulse_num);
        if (source_pulse <= previous_source_pulse) {
            std::cerr << "[SELFTEST] pulse extraction mapping is not monotonic"
                      << std::endl;
            return false;
        }
        previous_source_pulse = source_pulse;
    }
    if (mapTargetPulseToSourcePulse(
            0, target.src_pulse_num_expected, target.target_pulse_num) != 0 ||
        mapTargetPulseToSourcePulse(
            target.target_pulse_num - 1,
            target.src_pulse_num_expected,
            target.target_pulse_num) !=
            target.src_pulse_num_expected - 1) {
        std::cerr << "[SELFTEST] pulse extraction endpoints failed"
                  << std::endl;
        return false;
    }

    RangeFftResampler resampler(
        target.src_range_len_expected,
        target.target_fft_len,
        target.target_ddc_len);
    const std::complex<float> expected(1.25f, -0.5f);
    std::vector<std::complex<float>> src(
        static_cast<size_t>(target.src_range_len_expected), expected);
    std::vector<std::complex<float>> out;
    if (!resampler.process(src.data(), out) ||
        out.size() != static_cast<size_t>(target.target_ddc_len)) {
        std::cerr << "[SELFTEST] FFT resampler execution failed" << std::endl;
        return false;
    }
    float max_error = 0.0f;
    for (size_t i = 0; i < out.size(); ++i) {
        max_error = std::max(max_error, std::abs(out[i] - expected));
    }
    if (max_error > 1e-4f) {
        std::cerr << "[SELFTEST] FFT interpolation normalization failed, max_error="
                  << max_error << std::endl;
        return false;
    }

    PosSample pos;
    pos.vn = 3.0;
    pos.ve = 4.0;
    forcePosProfile(pos, target.force_alt_m, target.force_speed_mps);
    std::vector<uint8_t> packet(prt_len, 0);
    fillHeaderAndIns(packet, 5, 1, 1.0 / target.prf_hz, 0, pos,
                     target.scan_start_deg, prt_len);
    const int16_t angle_raw =
        static_cast<int16_t>(read_u16_le(packet, 218));
    const double speed =
        std::sqrt(static_cast<double>(read_f32_le(packet, 128)) *
                      read_f32_le(packet, 128) +
                  static_cast<double>(read_f32_le(packet, 132)) *
                      read_f32_le(packet, 132) +
                  static_cast<double>(read_f32_le(packet, 136)) *
                      read_f32_le(packet, 136));
    if (packet.size() != prt_len || read_u32_le(packet, 9) != prt_len ||
        read_u32_le(packet, 20) != 1U || angle_raw != -5900 ||
        std::fabs(read_f32_le(packet, 16) - 1.0 / target.prf_hz) > 1e-9 ||
        std::fabs(pos.alt_m - 6000.0) > 1e-9 ||
        std::fabs(pos.vd) > 1e-9 ||
        std::fabs(speed - 60.0) > 1e-4) {
        std::cerr << "[SELFTEST] target header fields failed" << std::endl;
        return false;
    }

    std::cout << "[SELFTEST] PASS"
              << "\n  FFT interpolation max error: " << max_error
              << "\n  packet length: " << prt_len
              << "\n  packet count: " << packets
              << "\n  output bytes: "
              << packets * static_cast<uint64_t>(prt_len)
              << "\n  beam mapping: original -23..25 deg beams remain unchanged"
              << "\n  replicated edges: target -59<-source -9, "
                 "target 61<-source 11"
              << "\n  pulse extraction: target 0->source 0, target 129->source 255"
              << "\n  first angle: -59 deg"
              << "\n  UTC interval: " << 1.0 / target.prf_hz << " s"
              << "\n  forced speed: " << speed << " m/s"
              << std::endl;
    return true;
}

void printTargetProfile(const TargetProfile &target,
                        int source_pulse_count,
                        uint32_t prt_len) {
    const uint64_t packets =
        static_cast<uint64_t>(target.target_beam_count) *
        static_cast<uint64_t>(target.target_pulse_num);
    const uint64_t bytes = packets * static_cast<uint64_t>(prt_len);
    std::cout << "[TargetProfile]\n"
              << "NOTE: synthetic test data; the original 25 beams remain "
                 "unchanged and complete copies are Doppler-shifted to extend "
                 "the scan,\n"
              << "      pulses are uniformly extracted over the full source aperture,\n"
              << "      range samples are FFT zero-pad interpolated, and UTC/"
                 "angle/POS are rewritten.\n"
              << "source expected: " << target.src_beam_count_expected
              << " beams, " << target.src_pulse_num_expected
              << " pulses/beam, " << target.src_range_len_expected
              << " range samples\n"
              << "target beams: " << target.target_beam_count << '\n'
              << "target scan: " << target.scan_start_deg << " to "
              << target.scan_start_deg +
                     target.scan_step_deg * (target.target_beam_count - 1)
              << " deg, step " << target.scan_step_deg << " deg\n"
              << "source scan assumed: " << target.src_scan_start_deg
              << " to "
              << target.src_scan_start_deg +
                     target.src_scan_step_deg *
                         (target.src_beam_count_expected - 1)
              << " deg, step " << target.src_scan_step_deg << " deg\n"
              << "replica Doppler model: -2*v/lambda*sin(angle), fc="
              << target.carrier_hz / 1.0e9 << " GHz, v="
              << target.force_speed_mps << " m/s\n"
              << "target pulses/beam: " << target.target_pulse_num << '\n'
              << "azimuth pulse conversion: uniform time-domain extraction "
              << source_pulse_count << " -> " << target.target_pulse_num
              << " (full aperture, endpoints retained)\n"
              << "target PRF: " << target.prf_hz << " Hz\n"
              << "target total pulses: " << packets << '\n'
              << "target duration: "
              << static_cast<double>(packets) / target.prf_hz << " s\n"
              << "target DDC length: " << target.target_ddc_len << '\n'
              << "target FFT length: " << target.target_fft_len << '\n'
              << "downstream range crop: start=" << target.range_crop_start
              << ", length=" << target.range_crop_len << '\n'
              << "target packet length: " << prt_len << " bytes\n"
              << "expected output bytes: " << bytes << std::endl;
}

bool convertToTargetProfile(const ConvertConfig &cfg,
                            const std::vector<PosRow> &pos_rows,
                            const std::string &out_file,
                            int beam_start,
                            int beam_count,
                            double iq_scale,
                            uint8_t mode,
                            int debug_source_beam,
                            double debug_target_angle_deg,
                            const DebugSelection &debug_selection) {
    const TargetProfile target;
    if (iq_scale != 1.0) {
        std::cerr << "[WARN] iq_scale is ignored because payload IQ is float32"
                  << std::endl;
    }
    if (cfg.pulse_len != target.src_range_len_expected) {
        std::cerr << "[WARN] cfg.pulse_len=" << cfg.pulse_len
                  << ", expected " << target.src_range_len_expected
                  << " for this target profile" << std::endl;
    }
    if (cfg.pulse_num != target.src_pulse_num_expected) {
        std::cerr << "[WARN] cfg.pulse_num=" << cfg.pulse_num
                  << ", expected " << target.src_pulse_num_expected
                  << " for this target profile" << std::endl;
    }
    if (cfg.pulse_num < target.target_pulse_num) {
        std::cerr << "source pulse_num must be >= target pulse_num" << std::endl;
        return false;
    }

    const uint32_t beams1 =
        calcBeamCount(cfg.data_ch1, cfg.info_len, cfg.pulse_len, cfg.pulse_num);
    const uint32_t beams2 =
        calcBeamCount(cfg.data_ch2, cfg.info_len, cfg.pulse_len, cfg.pulse_num);
    if (beams1 == 0U || beams2 == 0U) {
        std::cerr << "Cannot infer beam count from one or both source files"
                  << std::endl;
        return false;
    }
    if (beams1 != beams2) {
        std::cerr << "[WARN] source channel beam counts differ: ch1=" << beams1
                  << ", ch2=" << beams2 << std::endl;
    }

    const uint64_t total_beams = std::min(beams1, beams2);
    if (total_beams ==
        static_cast<uint64_t>(target.src_beam_count_expected + 1)) {
        std::cout << "[INFO] source file contains " << total_beams
                  << " complete old-format beams; target profile will use "
                  << target.src_beam_count_expected << std::endl;
    } else if (total_beams !=
               static_cast<uint64_t>(target.src_beam_count_expected)) {
        std::cerr << "[WARN] source file contains " << total_beams
                  << " complete beams, expected "
                  << target.src_beam_count_expected
                  << "; the selected source window will still default to 25"
                  << std::endl;
    }
    const uint64_t total_pulses =
        total_beams * static_cast<uint64_t>(cfg.pulse_num);
    if (cfg.skip_az_num < 0 ||
        static_cast<uint64_t>(cfg.skip_az_num) >= total_pulses) {
        std::cerr << "Invalid skip_pulses for source file" << std::endl;
        return false;
    }
    const int usable_beams = static_cast<int>(
        (total_pulses - static_cast<uint64_t>(cfg.skip_az_num)) /
        static_cast<uint64_t>(cfg.pulse_num));
    const int skipped_full_beams = cfg.skip_az_num / cfg.pulse_num;
    const int skipped_partial_pulses = cfg.skip_az_num % cfg.pulse_num;
    std::cout << "source old-format beams: total=" << total_beams
              << ", usable=" << usable_beams
              << ", skip_pulses=" << cfg.skip_az_num << std::endl;
    // Auto mode discards the first beam when the old-format source contains
    // 26 complete beams, leaving source beams 2..26 as the expected 25-beam
    // window. An explicit command-line source_beam_start overrides this.
    const int source_start_beam =
        usable_beams == 26 && beam_start <= 1
            ? 2
            : (beam_start > 0 ? beam_start : 1);
    if (usable_beams == 26 && beam_start <= 1) {
        std::cout << "[INFO] source contains 26 beams; dropping source beam 1 "
                     "and using beams 2..26"
                  << std::endl;
    }
    const int available_from_start = usable_beams - source_start_beam + 1;
    if (available_from_start <= 0) {
        std::cerr << "No source beam available from beam_start="
                  << source_start_beam << std::endl;
        return false;
    }

    // A zero beam_count means the target profile's expected 25-beam source
    // window, clamped only when the source file has fewer usable beams.
    const int requested_source_beams =
        beam_count > 0 ? beam_count : target.src_beam_count_expected;
    const int source_beam_count =
        std::min(requested_source_beams, available_from_start);
    if (source_beam_count != target.src_beam_count_expected) {
        std::cerr << "[WARN] source beam count is " << source_beam_count
                  << ", expected " << target.src_beam_count_expected
                  << " for target profile" << std::endl;
    }
    if (beam_count > available_from_start) {
        std::cerr << "[WARN] requested beam_count=" << beam_count
                  << " clamped to " << source_beam_count << std::endl;
    }

    Config rcfg{};
    rcfg.GMTI_Data_add = cfg.data_ch1;
    rcfg.GMTI_Data_add2 = cfg.data_ch2;
    rcfg.info_len = cfg.info_len;
    rcfg.pulse_len = cfg.pulse_len;
    rcfg.pulse_num = cfg.pulse_num;
    rcfg.skip_az_num = cfg.skip_az_num;
    rcfg.week = cfg.week;
    rcfg.secBias = cfg.sec_bias;

    const uint32_t prt_len =
        256U + static_cast<uint32_t>(target.target_ddc_len) * 16U;
    printTargetProfile(target, cfg.pulse_num, prt_len);
    std::cout << "source beam window: start=" << source_start_beam
              << ", count=" << source_beam_count;
    if (skipped_partial_pulses == 0) {
        const int physical_start =
            skipped_full_beams + source_start_beam;
        const int physical_end =
            physical_start + source_beam_count - 1;
        std::cout << " (physical old-format beams "
                  << physical_start << ".." << physical_end << ")";
        if (physical_start == 2 && physical_end == 26) {
            std::cout << " [source beam 1 dropped]";
        }
    } else {
        std::cout << " (skip_pulses is not beam-aligned)";
    }
    std::cout << std::endl;

    for (int selected_source : debug_selection.source_beams) {
        if (selected_source < 1 || selected_source > source_beam_count) {
            std::cerr << "selected_source_beams item " << selected_source
                      << " is outside current source window 1.."
                      << source_beam_count << std::endl;
            return false;
        }
    }

    std::vector<int> selected_target_beams;
    selected_target_beams.reserve(static_cast<size_t>(target.target_beam_count));
    const double angle_tolerance =
        std::max(1e-6, 0.25 * std::fabs(target.scan_step_deg));
    for (int target_beam = 0;
         target_beam < target.target_beam_count;
         ++target_beam) {
        BeamReplicaMapping mapping =
            mapTargetBeamToReplicatedSource(
                target_beam, source_start_beam, source_beam_count, target);
        const int relative_source_beam =
            mapping.source_beam - source_start_beam + 1;
        const bool source_selected =
            debug_selection.source_beams.empty() ||
            containsInt(debug_selection.source_beams, relative_source_beam);
        const bool target_angle_selected =
            debug_selection.target_angles_deg.empty() ||
            containsAngle(debug_selection.target_angles_deg,
                          mapping.target_angle_deg,
                          angle_tolerance);
        if (source_selected && target_angle_selected) {
            selected_target_beams.push_back(target_beam);
        }
    }
    if (selected_target_beams.empty()) {
        std::cerr << "No target beams selected. selected_source_beams="
                  << joinInts(debug_selection.source_beams)
                  << ", selected_target_angles_deg="
                  << joinDoubles(debug_selection.target_angles_deg)
                  << std::endl;
        return false;
    }

    const uint64_t data_packets =
        static_cast<uint64_t>(selected_target_beams.size()) *
        static_cast<uint64_t>(target.target_pulse_num);
    const uint64_t expected_packets =
        static_cast<uint64_t>(target.target_beam_count) *
        static_cast<uint64_t>(target.target_pulse_num);
    std::cout << "[DEBUG_SELECT] selected_source_beams="
              << joinInts(debug_selection.source_beams)
              << ", selected_target_angles_deg="
              << joinDoubles(debug_selection.target_angles_deg)
              << ", data target beams="
              << selected_target_beams.size() << " / "
              << target.target_beam_count
              << ", zero-filled target beams="
              << target.target_beam_count -
                     static_cast<int>(selected_target_beams.size())
              << ", data packets=" << data_packets
              << ", full output packets=" << expected_packets
              << ", full output bytes="
              << expected_packets * static_cast<uint64_t>(prt_len)
              << std::endl;
    std::cout << "[DEBUG_SELECT] mapping preview:" << std::endl;
    for (int target_beam : selected_target_beams) {
        const BeamReplicaMapping mapping =
            mapTargetBeamToReplicatedSource(
                target_beam, source_start_beam, source_beam_count, target);
        const int relative_source_beam =
            mapping.source_beam - source_start_beam + 1;
        std::cout << "  out_target_beam=" << (target_beam + 1)
                  << ", target_angle=" << mapping.target_angle_deg
                  << " deg <- rel_source_beam=" << relative_source_beam
                  << ", file_source_beam=" << mapping.source_beam
                  << ", source_angle=" << mapping.source_angle_deg
                  << " deg" << std::endl;
    }

    std::ofstream fout(out_file, std::ios::binary);
    if (!fout) {
        std::cerr << "Unable to open output file: " << out_file << std::endl;
        return false;
    }

    RangeFftResampler resampler(
        cfg.pulse_len, target.target_fft_len, target.target_ddc_len);
    std::vector<uint8_t> packet(prt_len, 0);
    std::vector<std::complex<float>> data1;
    std::vector<std::complex<float>> data2;
    std::vector<std::complex<float>> ch1_out;
    std::vector<std::complex<float>> ch2_out;
    std::vector<double> utc;
    std::vector<double> utc_ch1;
    std::vector<double> fw_angle_unused;
    std::vector<double> fw_angle_unused_ch1;

    uint32_t prt_counter = 0;
    if (!readBeamRawFloat(rcfg, cfg.data_ch2.c_str(), source_start_beam,
                          data2, fw_angle_unused, utc)) {
        std::cerr << "readBeamRawFloat failed while probing UTC, beam="
                  << source_start_beam << std::endl;
        return false;
    }
    if (utc.empty()) {
        std::cerr << "Source UTC vector is empty while probing beam="
                  << source_start_beam << std::endl;
        return false;
    }
    const double utc0 = utc.front();
    const TargetMotion target_motion = buildTargetMotion(
        pos_rows, utc0, target.force_alt_m, target.force_speed_mps);
    std::cout << "target motion: vn=" << target_motion.vn_mps
              << " m/s, ve=" << target_motion.ve_mps
              << " m/s, vd=0 m/s, horizontal speed="
              << std::sqrt(target_motion.vn_mps * target_motion.vn_mps +
                           target_motion.ve_mps * target_motion.ve_mps)
              << " m/s" << std::endl;

    size_t selected_index = 0;
    for (int target_beam = 0;
         target_beam < target.target_beam_count;
         ++target_beam) {
        const bool write_real_data =
            containsInt(selected_target_beams, target_beam);
        BeamReplicaMapping beam_mapping =
            mapTargetBeamToReplicatedSource(
                target_beam, source_start_beam, source_beam_count, target);
        const double target_angle_deg = beam_mapping.target_angle_deg;
        if (!write_real_data) {
            int16_t first_header_angle_raw = 0;
            uint32_t first_header_prt_len = 0;
            for (int target_pulse = 0;
                 target_pulse < target.target_pulse_num;
                 ++target_pulse) {
                const double utc_out =
                    utc0 + static_cast<double>(prt_counter) / target.prf_hz;
                const PosSample pos =
                    evaluateTargetMotion(target_motion, utc_out);
                fillHeaderAndIns(packet, mode, prt_counter, utc_out, cfg.week,
                                 pos, target_angle_deg, prt_len);
                if (target_pulse == 0) {
                    first_header_prt_len = read_u32_le(packet, 9);
                    first_header_angle_raw =
                        static_cast<int16_t>(read_u16_le(packet, 218));
                }
                fout.write(reinterpret_cast<const char *>(packet.data()),
                           static_cast<std::streamsize>(packet.size()));
                if (!fout) {
                    std::cerr << "Write zero-fill failed at target_beam="
                              << target_beam
                              << " target_pulse=" << target_pulse
                              << std::endl;
                    return false;
                }
                ++prt_counter;
            }
            std::cout << "Zero-filled target beam " << target_beam + 1
                      << " / " << target.target_beam_count
                      << ", target_angle=" << target_angle_deg
                      << " deg, header_angle_raw="
                      << first_header_angle_raw
                      << ", header_prt_len=" << first_header_prt_len
                      << std::endl;
            continue;
        }

        const bool debug_override =
            debug_source_beam > 0 &&
            debug_source_beam <= source_beam_count &&
            std::isfinite(debug_target_angle_deg) &&
            std::fabs(beam_mapping.target_angle_deg -
                      debug_target_angle_deg) <=
                0.5 * std::fabs(target.scan_step_deg);
        if (debug_override) {
            const int source_offset = debug_source_beam - 1;
            beam_mapping.source_beam = source_start_beam + source_offset;
            beam_mapping.source_angle_deg =
                target.src_scan_start_deg +
                target.src_scan_step_deg * source_offset;
            std::cout << "[DEBUG] target beam " << target_beam + 1
                      << " at " << beam_mapping.target_angle_deg
                      << " deg is forced to source beam "
                      << debug_source_beam << " (source angle "
                      << beam_mapping.source_angle_deg << " deg)"
                      << std::endl;
        }
        const int src_beam = beam_mapping.source_beam;
        const int relative_source_beam =
            src_beam - source_start_beam + 1;

        if (!readBeamRawFloat(rcfg, cfg.data_ch2.c_str(), src_beam,
                              data2, fw_angle_unused, utc)) {
            std::cerr << "readBeamRawFloat failed on channel-2, beam="
                      << src_beam << std::endl;
            return false;
        }
        if (!readBeamRawFloat(rcfg, cfg.data_ch1.c_str(), src_beam,
                              data1, fw_angle_unused_ch1, utc_ch1)) {
            std::cerr << "readBeamRawFloat failed on channel-1, beam="
                      << src_beam << std::endl;
            return false;
        }

        const size_t expected_samples =
            static_cast<size_t>(cfg.pulse_num) *
            static_cast<size_t>(cfg.pulse_len);
        if (data1.size() != expected_samples ||
            data2.size() != expected_samples ||
            utc.size() != static_cast<size_t>(cfg.pulse_num) ||
            utc_ch1.size() != static_cast<size_t>(cfg.pulse_num)) {
            std::cerr << "Source beam data size mismatch at beam="
                      << src_beam << std::endl;
            return false;
        }
        const double measured_source_wrapped_hz =
            estimateMappedPulseDopplerHz(
                data1, cfg.pulse_num, cfg.pulse_len,
                target.target_pulse_num, target.prf_hz);
        const double doppler_shift_hz =
            dopplerShiftForReplicaFromMeasuredSource(
                beam_mapping, target, measured_source_wrapped_hz);
        const double desired_wrapped_hz = wrapFrequencyHz(
            theoreticalDopplerHz(target_angle_deg,
                                 target.force_speed_mps,
                                 target.carrier_hz),
            target.prf_hz);
        const double theoretical_source_wrapped_hz = wrapFrequencyHz(
            theoreticalDopplerHz(beam_mapping.source_angle_deg,
                                 target.force_speed_mps,
                                 target.carrier_hz),
            target.prf_hz);
        PulsePowerStats stats_before;
        PulsePowerStats stats_after;
        const int first_src_pulse = mapTargetPulseToSourcePulse(
            0, cfg.pulse_num, target.target_pulse_num);
        const int last_src_pulse = mapTargetPulseToSourcePulse(
            target.target_pulse_num - 1, cfg.pulse_num,
            target.target_pulse_num);
        int16_t first_header_angle_raw = 0;
        uint32_t first_header_prt_len = 0;

        for (int target_pulse = 0;
             target_pulse < target.target_pulse_num;
             ++target_pulse) {
            const int src_pulse = mapTargetPulseToSourcePulse(
                target_pulse, cfg.pulse_num, target.target_pulse_num);
            const size_t src_offset =
                static_cast<size_t>(src_pulse) *
                static_cast<size_t>(cfg.pulse_len);
            if (!resampler.process(data1.data() + src_offset, ch1_out)) {
                std::cerr << "range resample failed for ch1 target_beam="
                          << target_beam << " target_pulse=" << target_pulse
                          << std::endl;
                return false;
            }
            if (!resampler.process(data2.data() + src_offset, ch2_out)) {
                std::cerr << "range resample failed for ch2 target_beam="
                          << target_beam << " target_pulse=" << target_pulse
                          << std::endl;
                return false;
            }
            accumulatePowerStats(ch1_out, ch2_out, stats_before);
            applySlowTimeFrequencyShift(
                ch1_out, target_pulse, doppler_shift_hz, target.prf_hz);
            applySlowTimeFrequencyShift(
                ch2_out, target_pulse, doppler_shift_hz, target.prf_hz);
            accumulatePowerStats(ch1_out, ch2_out, stats_after);

            const double utc_out =
                utc0 + static_cast<double>(prt_counter) / target.prf_hz;
            const PosSample pos =
                evaluateTargetMotion(target_motion, utc_out);
            fillHeaderAndIns(packet, mode, prt_counter, utc_out, cfg.week,
                             pos, target_angle_deg, prt_len);
            if (target_pulse == 0) {
                first_header_prt_len = read_u32_le(packet, 9);
                first_header_angle_raw =
                    static_cast<int16_t>(read_u16_le(packet, 218));
            }
            fillGmtiIqPayloadOnePulse(packet, ch1_out, ch2_out);
            fout.write(reinterpret_cast<const char *>(packet.data()),
                       static_cast<std::streamsize>(packet.size()));
            if (!fout) {
                std::cerr << "Write failed at target_beam=" << target_beam
                          << " target_pulse=" << target_pulse << std::endl;
                return false;
            }
            ++prt_counter;
        }

        std::cout << "Converted selected beam " << (selected_index + 1)
                  << " / " << selected_target_beams.size()
                  << ", out_target_beam=" << target_beam + 1
                  << " / " << target.target_beam_count
                  << ", rel_source_beam=" << relative_source_beam
                  << ", file_source_beam=" << src_beam
                  << ", src_angle=" << beam_mapping.source_angle_deg
                  << " deg, target_angle=" << target_angle_deg
                  << " deg, theoretical_source_wrapped_fd="
                  << theoretical_source_wrapped_hz
                  << " Hz, measured_source_fd=" << measured_source_wrapped_hz
                  << " Hz, desired_wrapped_fd=" << desired_wrapped_hz
                  << " Hz, Doppler shift=" << doppler_shift_hz
                  << " Hz, pulse_map=" << first_src_pulse << ".."
                  << last_src_pulse
                  << ", mean_power_before(ch1,ch2)="
                  << meanPower(stats_before, 1) << ","
                  << meanPower(stats_before, 2)
                  << ", mean_power_after(ch1,ch2)="
                  << meanPower(stats_after, 1) << ","
                  << meanPower(stats_after, 2)
                  << ", header_angle_raw=" << first_header_angle_raw
                  << ", header_prt_len=" << first_header_prt_len
                  << std::endl;
        ++selected_index;
    }

    if (prt_counter != expected_packets) {
        std::cerr << "Internal packet-count mismatch: wrote " << prt_counter
                  << ", expected " << expected_packets << std::endl;
        return false;
    }
    fout.close();
    if (!fout) {
        std::cerr << "Failed while closing output file" << std::endl;
        return false;
    }

    std::cout << "Conversion done, total PRT packets: " << prt_counter
              << "\nOutput file: " << out_file
              << "\nOutput bytes: "
              << expected_packets * static_cast<uint64_t>(prt_len)
              << std::endl;
    return true;
}

} // namespace

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        try {
            return runSelfTest() ? 0 : 5;
        } catch (const std::exception &e) {
            std::cerr << "[SELFTEST] exception: " << e.what() << std::endl;
            return 5;
        }
    }
    if (argc < 2) {
        std::cerr
            << "Usage: " << argv[0]
            << " <xml_path> [output_target_echo.bin] [source_beam_start=auto]"
               " [source_beam_count=25] [iq_scale=1.0(ignored)] [mode=5]"
               " [debug_source_beam=0] [debug_target_angle_deg=nan]"
               " [selected_source_beams=all] [selected_target_angles_deg=all]\n"
            << "       " << argv[0] << " --self-test\n"
            << "If output is omitted, <GMTI_data_new> from XML is used.\n"
            << "selected_source_beams uses relative source indices in the "
               "current 25-beam window, e.g. 1,2,19.\n"
            << "selected_target_angles_deg filters output target angles, "
               "e.g. -23,27. Use all or * to disable a filter."
            << std::endl;
        return 1;
    }

    const std::string xml_path = argv[1];
    const bool has_output_arg = argc >= 3;
    const int beam_start =
        argc > 3 ? std::max(1, std::atoi(argv[3])) : 0;
    const int beam_count = argc > 4 ? std::atoi(argv[4]) : 0;
    const double iq_scale = argc > 5 ? std::atof(argv[5]) : 1.0;
    const uint8_t mode = sat_u8(argc > 6 ? std::atoi(argv[6]) : 5);
    const int debug_source_beam = argc > 7 ? std::atoi(argv[7]) : 0;
    const double debug_target_angle_deg =
        argc > 8 ? std::atof(argv[8]) :
                   std::numeric_limits<double>::quiet_NaN();
    DebugSelection debug_selection;
    if (argc > 9 &&
        !parseIntCsv(argv[9], debug_selection.source_beams,
                     "selected_source_beams")) {
        return 1;
    }
    if (argc > 10 &&
        !parseDoubleCsv(argv[10], debug_selection.target_angles_deg,
                        "selected_target_angles_deg")) {
        return 1;
    }
    if (beam_count < 0) {
        std::cerr << "source_beam_count must be >= 0" << std::endl;
        return 1;
    }

    ConvertConfig cfg;
    if (!loadConvertConfig(xml_path, cfg)) {
        return 2;
    }
    const std::string out_file =
        has_output_arg ? argv[2] : cfg.data_new;
    if (out_file.empty()) {
        std::cerr
            << "Output path is empty. Provide argv[2] or <GMTI_data_new>."
            << std::endl;
        return 2;
    }
    if (out_file == cfg.data_ch1 || out_file == cfg.data_ch2) {
        std::cerr << "Output path must not overwrite a source echo file"
                  << std::endl;
        return 2;
    }

    std::vector<PosRow> pos_rows;
    if (!readPosFileLikePOSDataread(cfg.pos_file, pos_rows)) {
        return 3;
    }

    try {
        if (!convertToTargetProfile(cfg, pos_rows, out_file, beam_start,
                                    beam_count, iq_scale, mode,
                                    debug_source_beam,
                                    debug_target_angle_deg,
                                    debug_selection)) {
            return 4;
        }
    } catch (const std::exception &e) {
        std::cerr << "Target-profile generation failed: " << e.what()
                  << std::endl;
        return 4;
    }
    return 0;
}
