#include "trackModule.hpp"
#include "GMTIProcessor.hpp"
#include "TrackManager.hpp"
#include "geo/geoProj.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <dirent.h>

// 手动实现 clamp，替代 std::clamp
template<typename T>
T clamp(T val, T min_val, T max_val) {
    if (val < min_val) return min_val;
    if (val > max_val) return max_val;
    return val;
}

// 前向声明
GMTIResult read_mt_result_bin(const std::string& binFile);
struct CurrentTargetPacket {
    uint16_t id;
    double lon;
    double lat;
    double speed;
    double direction;
    double range;
    double utc;
};
void writeCurrentCyclePackets(const std::vector<CurrentTargetPacket>& targets, const std::string& fileName);

namespace {

struct SlidingTrack {
    int id = 0;
    double e = 0.0;
    double n = 0.0;
    double utc = 0.0;
    double speed = 0.0;
    double direction = 0.0;
    double range = 0.0;
    int lastFrame = -1;
    std::deque<int> matchedFrames;
};

struct AssocEdge {
    double d;
    int trackIdx;
    int detIdx;
};

static inline void put_u16_le(std::vector<uint8_t>& buf, size_t pos, uint16_t v) {
    buf[pos + 0] = static_cast<uint8_t>(v & 0xFF);
    buf[pos + 1] = static_cast<uint8_t>((v >> 8) & 0xFF);
}

static inline void put_i32_le(std::vector<uint8_t>& buf, size_t pos, int32_t v) {
    const uint32_t u = static_cast<uint32_t>(v);
    buf[pos + 0] = static_cast<uint8_t>(u & 0xFF);
    buf[pos + 1] = static_cast<uint8_t>((u >> 8) & 0xFF);
    buf[pos + 2] = static_cast<uint8_t>((u >> 16) & 0xFF);
    buf[pos + 3] = static_cast<uint8_t>((u >> 24) & 0xFF);
}

static inline void put_f64_le(std::vector<uint8_t>& buf, size_t pos, double v) {
    static_assert(sizeof(double) == sizeof(uint64_t), "double must be 8 bytes");
    uint64_t u = 0;
    std::memcpy(&u, &v, sizeof(u));
    for (int i = 0; i < 8; ++i) {
        buf[pos + static_cast<size_t>(i)] = static_cast<uint8_t>((u >> (8 * i)) & 0xFF);
    }
}

static inline int32_t quant_deg_to_i32(double deg) {
    double q = deg / LSB;
    q = (q >= 0.0) ? std::floor(q + 0.5) : std::ceil(q - 0.5);
    if (q > 2147483647.0) q = 2147483647.0;
    if (q < -2147483648.0) q = -2147483648.0;
    return static_cast<int32_t>(q);
}

static inline uint16_t quant_speed_to_u16(double speed) {
    const double q = std::round(std::max(0.0, speed) / 0.01);
    if (q > 65535.0) return 65535;
    return static_cast<uint16_t>(q);
}

static inline double wrap180_deg(double angle_deg) {
    angle_deg = std::fmod(angle_deg + 180.0, 360.0);
    if (angle_deg < 0.0) angle_deg += 360.0;
    return angle_deg - 180.0;
}

static void printCurrentTargets(const std::vector<CurrentTargetPacket>& targets) {
    std::cout << "\n[TRACK] 目标信息列表" << std::endl;
    std::cout << "-----------------------------------------------------------------------------------------------" << std::endl;
    std::cout << std::left
              << std::setw(8)  << "TrackID"
              << std::setw(14) << "Lon(deg)"
              << std::setw(14) << "Lat(deg)"
              << std::setw(12) << "Speed(m/s)"
              << std::setw(12) << "Dir(deg)"
              << std::setw(12) << "Range(m)"
              << std::setw(14) << "UTC(s)"
              << std::endl;
    std::cout << "-----------------------------------------------------------------------------------------------" << std::endl;

    if (targets.empty()) {
        std::cout << "(empty)" << std::endl;
    } else {
        for (const auto& t : targets) {
            std::cout << std::left
                      << std::setw(8)  << t.id
                      << std::setw(14) << std::fixed << std::setprecision(6) << t.lon
                      << std::setw(14) << std::fixed << std::setprecision(6) << t.lat
                      << std::setw(12) << std::fixed << std::setprecision(2) << t.speed
                      << std::setw(12) << std::fixed << std::setprecision(2) << t.direction
                      << std::setw(12) << std::fixed << std::setprecision(2) << t.range
                      << std::setw(14) << std::fixed << std::setprecision(3) << t.utc
                      << std::endl;
        }
    }
    std::cout << "-----------------------------------------------------------------------------------------------\n" << std::endl;
}

static int resolveCurrentResultId(const Config& cfg)
{
    if (cfg.result_file_id > 0) {
        return cfg.result_file_id;
    }
    if (!cfg.track_idx_range.empty()) {
        return cfg.track_idx_range.back();
    }

    int latest = -1;
    DIR* dir = opendir(cfg.result_add.c_str());
    if (!dir) {
        return latest;
    }
    while (dirent* ent = readdir(dir)) {
        int id = -1;
        char tail = '\0';
        if (std::sscanf(ent->d_name, "GMTI%02d.bin%c", &id, &tail) == 1 && id > latest) {
            latest = id;
        }
    }
    closedir(dir);
    return latest;
}

static std::string makeGmtiBinPath(const std::string& result_dir, int result_id)
{
    char name[256];
    std::snprintf(name, sizeof(name), "GMTI%02d.bin", result_id);
    return result_dir + "/" + name;
}

static std::string makeGmtiTrackPath(const std::string& result_dir, int result_id)
{
    char name[256];
    std::snprintf(name, sizeof(name), "GMTI%02d_track.bin", result_id);
    return result_dir + "/" + name;
}

static void dumpTrackPointsCsv(const Config& cfg,
                               int result_id,
                               const std::vector<GMTIDetection>& outputs)
{
    const std::string path = cfg.result_add + "/track_points.csv";

    bool need_header = true;
    {
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        if (in.is_open()) {
            need_header = (in.tellg() <= 0);
        }
    }

    std::ofstream out(path, std::ios::out | std::ios::app);
    if (!out.is_open()) {
        std::cout << "[TRACK_ONLINE][WARN] cannot open track points csv: "
                  << path << std::endl;
        return;
    }

    if (need_header) {
        out << "result_id,track_id,e,n,lat,lon\n";
        if (!out.good()) {
            std::cout << "[TRACK_ONLINE][WARN] cannot write track points csv header: "
                      << path << std::endl;
            return;
        }
    }

    for (const auto& det : outputs) {
        out << result_id << ","
            << det.id << ","
            << std::fixed << std::setprecision(6) << det.e << ","
            << std::fixed << std::setprecision(6) << det.n << ","
            << std::fixed << std::setprecision(12) << det.lat << ","
            << std::fixed << std::setprecision(12) << det.lon << "\n";
    }

    if (!out.good()) {
        std::cout << "[TRACK_ONLINE][WARN] failed writing track points csv: "
                  << path << std::endl;
    }
}

static bool readCurrentDetections(const std::string& input_path,
                                  std::vector<GMTIDetection>& dets,
                                  double& utc_global)
{
    try {
        GMTIResult result = read_mt_result_bin(input_path);
        dets = result.targets;
        utc_global = result.utcGlobal;
        return true;
    } catch (const std::exception& ex) {
        std::cout << "[TRACK_ONLINE][WARN] read failed: " << input_path
                  << " error=" << ex.what() << std::endl;
    } catch (...) {
        std::cout << "[TRACK_ONLINE][WARN] read failed: " << input_path << std::endl;
    }
    dets.clear();
    utc_global = 0.0;
    return false;
}

static double estimateFrameUtc(const std::vector<GMTIDetection>& dets, double utc_global)
{
    if (std::isfinite(utc_global) && utc_global > 0.0) {
        return utc_global;
    }
    for (const auto& det : dets) {
        if (std::isfinite(det.utcMid) && det.utcMid > 0.0) {
            return det.utcMid;
        }
    }
    return 0.0;
}

static std::vector<CurrentTargetPacket> makeCurrentPackets(const std::vector<GMTIDetection>& detections)
{
    std::vector<CurrentTargetPacket> packets;
    packets.reserve(detections.size());
    for (const auto& det : detections) {
        packets.push_back(CurrentTargetPacket{
            det.id,
            det.lon,
            det.lat,
            det.speed,
            det.direction,
            det.range,
            det.utcMid
        });
    }
    return packets;
}

} // namespace

// Note: GMTIProcessor constructor is defined in the main header; do not redefine it here.

std::vector<GMTIDetection> trackModule(const Config& cfg) {
    return trackModule(cfg, nullptr);
}

std::vector<GMTIDetection> trackModule(const Config& cfg, TrackManager* manager) {
    const std::string& result_dir = cfg.result_add;
    const std::vector<int>& idx_range = cfg.track_idx_range;
    std::vector<GMTIDetection> latest_targets;

    if (idx_range.empty()) {
        std::cout << "idx_range 为空，无法进行航迹关联。" << std::endl;
        return latest_targets;
    }

    std::cout << "[TRACK] 关联输入目录: " << result_dir << std::endl;
    std::cout << "[TRACK] 关联周期 idx_range:";
    for (const int idx : idx_range) {
        std::cout << ' ' << idx;
    }
    std::cout << "  (truth_threshold=" << cfg.track_truth_threshold
              << ", gate_m=" << cfg.track_gate_m
              << ", v_max=" << cfg.track_v_max << ")" << std::endl;

    std::vector<std::vector<std::array<double,2>>> frame_dets;
    std::vector<std::vector<double>> frame_dirs;
    std::vector<std::vector<double>> frame_ranges;
    std::vector<std::vector<double>> frame_target_utc;
    frame_dets.reserve(idx_range.size());
    frame_dirs.reserve(idx_range.size());
    frame_ranges.reserve(idx_range.size());
    frame_target_utc.reserve(idx_range.size());
    int readable_frames = 0;
    int nonempty_frames = 0;
    int valid_utc_frames = 0;

    for (const int idx : idx_range) {
        char filename[256];
        sprintf(filename, "GMTI%02d.bin", idx);
        const std::string binPath = result_dir + "/" + filename;

        try {
            GMTIResult out = read_mt_result_bin(binPath);
            ++readable_frames;
            if (out.count > 0) {
                ++nonempty_frames;
            }
            if (std::isfinite(out.utcGlobal) && out.utcGlobal > 0.0) {
                ++valid_utc_frames;
            } else {
                std::cout << "[TRACK][WARN] idx=" << idx
                          << " UTC 非法或缺失，速度估计将使用帧序兜底: utc="
                          << out.utcGlobal << std::endl;
            }
            std::cout << "[TRACK] read idx=" << idx
                      << " file=" << binPath
                      << " count=" << out.count
                      << " utc=" << out.utcGlobal << std::endl;

            std::vector<std::array<double,2>> dets;
            std::vector<double> dirs;
            std::vector<double> ranges;
            std::vector<double> target_utcs;
            dets.reserve(static_cast<size_t>(out.count));
            dirs.reserve(static_cast<size_t>(out.count));
            ranges.reserve(static_cast<size_t>(out.count));
            target_utcs.reserve(static_cast<size_t>(out.count));
            for (int i = 0; i < out.count; ++i) {
                double E = 0.0, N = 0.0;
                Gaussp3(out.targets[i].lat, out.targets[i].lon, cfg.L0, E, N);

                dets.push_back({{E, N}});
                dirs.push_back(wrap180_deg(out.targets[i].direction));
                ranges.push_back(out.targets[i].range);
                target_utcs.push_back(out.targets[i].utcMid);
            }
            frame_dets.push_back(dets);
            frame_dirs.push_back(dirs);
            frame_ranges.push_back(ranges);
            frame_target_utc.push_back(target_utcs);
        } catch (const std::exception& ex) {
            std::cout << "[TRACK][WARN] 跳过序号 " << idx
                      << ": 文件处理出错: " << ex.what() << std::endl;
            frame_dets.push_back(std::vector<std::array<double,2>>());
            frame_dirs.push_back(std::vector<double>());
            frame_ranges.push_back(std::vector<double>());
            frame_target_utc.push_back(std::vector<double>());
        } catch (...) {
            std::cout << "[TRACK][WARN] 跳过序号 " << idx
                      << ": 文件处理出错。" << std::endl;
            frame_dets.push_back(std::vector<std::array<double,2>>());
            frame_dirs.push_back(std::vector<double>());
            frame_ranges.push_back(std::vector<double>());
            frame_target_utc.push_back(std::vector<double>());
        }
    }

    if (frame_dets.empty()) {
        std::cout << "[TRACK][WARN] 没有可用周期数据。" << std::endl;
        return latest_targets;
    }

    const int effective_window = std::max(1, static_cast<int>(frame_dets.size()));
    const int min_hits = std::max(1, cfg.track_truth_threshold);
    const double gate_m = std::max(0.0, cfg.track_gate_m);
    const double v_max = std::max(0.0, cfg.track_v_max);
    if (readable_frames < static_cast<int>(idx_range.size())) {
        std::cout << "[TRACK][WARN] 关联窗口存在缺失文件: readable_frames="
                  << readable_frames << "/" << idx_range.size() << std::endl;
    }
    if (valid_utc_frames < readable_frames) {
        std::cout << "[TRACK][WARN] 关联窗口存在 UTC 非法帧: valid_utc_frames="
                  << valid_utc_frames << "/" << readable_frames << std::endl;
    }
    if (nonempty_frames < min_hits) {
        std::cout << "[TRACK][WARN] 有效非空帧不足: nonempty_frames="
                  << nonempty_frames << ", truth_threshold=" << min_hits
                  << "，可能无法形成当前周期真目标。" << std::endl;
    }

    std::vector<SlidingTrack> tracks;
    tracks.reserve(128);
    int next_track_id = 1;

    std::vector<int> latest_det_to_track;
    const int last_k = static_cast<int>(frame_dets.size()) - 1;

    for (int k = 0; k < static_cast<int>(frame_dets.size()); ++k) {
        const auto& dets = frame_dets[k];
        std::vector<int> det_to_track(dets.size(), -1);
        std::vector<char> track_used(tracks.size(), 0);
        std::vector<char> det_used(dets.size(), 0);

        std::vector<AssocEdge> edges;
        edges.reserve(tracks.size() * std::max<size_t>(1, dets.size()));

        for (size_t ti = 0; ti < tracks.size(); ++ti) {
            const SlidingTrack& tr = tracks[ti];
            if (tr.lastFrame < 0 || k - tr.lastFrame > (effective_window - 1)) {
                continue;
            }
            for (size_t dj = 0; dj < dets.size(); ++dj) {
                const double dx = dets[dj][0] - tr.e;
                const double dy = dets[dj][1] - tr.n;
                const double d = std::sqrt(dx * dx + dy * dy);
                if (d < gate_m) {
                    edges.push_back({d, static_cast<int>(ti), static_cast<int>(dj)});
                }
            }
        }

        std::sort(edges.begin(), edges.end(), [](const AssocEdge& a, const AssocEdge& b) {
            return a.d < b.d;
        });

        for (const auto& e : edges) {
            if (track_used[e.trackIdx] || det_used[e.detIdx]) {
                continue;
            }

            SlidingTrack& tr = tracks[e.trackIdx];
            const double prev_e = tr.e;
            const double prev_n = tr.n;
            const double prev_utc = tr.utc;
            const double cur_utc = (e.detIdx < static_cast<int>(frame_target_utc[k].size()))
                ? frame_target_utc[k][e.detIdx]
                : prev_utc;
            const double dist = std::sqrt((dets[e.detIdx][0] - prev_e) * (dets[e.detIdx][0] - prev_e)
                                        + (dets[e.detIdx][1] - prev_n) * (dets[e.detIdx][1] - prev_n));

            double dt = cur_utc - prev_utc;
            if (dt <= 1e-3) {
                dt = static_cast<double>(std::max(1, k - tr.lastFrame));
            }

            const double speed = dist / dt;
            if (v_max > 0.0 && speed > v_max) {
                continue;
            }

            tr.speed = speed;
            tr.e = dets[e.detIdx][0];
            tr.n = dets[e.detIdx][1];
            tr.utc = cur_utc;
            tr.direction = wrap180_deg(frame_dirs[k][e.detIdx]);
            tr.range = frame_ranges[k][e.detIdx];
            tr.lastFrame = k;
            if (tr.matchedFrames.empty() || tr.matchedFrames.back() != k) {
                tr.matchedFrames.push_back(k);
            }

            det_to_track[e.detIdx] = e.trackIdx;
            track_used[e.trackIdx] = 1;
            det_used[e.detIdx] = 1;
        }

        for (size_t dj = 0; dj < dets.size(); ++dj) {
            if (det_used[dj]) {
                continue;
            }

            SlidingTrack tr;
            tr.id = next_track_id++;
            tr.e = dets[dj][0];
            tr.n = dets[dj][1];
            tr.utc = (dj < frame_target_utc[k].size()) ? frame_target_utc[k][dj] : 0.0;
            tr.speed = 0.0;
            tr.direction = wrap180_deg(frame_dirs[k][dj]);
            tr.range = frame_ranges[k][dj];
            tr.lastFrame = k;
            tr.matchedFrames.push_back(k);
            tracks.push_back(tr);

            det_to_track[dj] = static_cast<int>(tracks.size()) - 1;
        }

        const int keep_from = std::max(0, k - effective_window + 1);
        for (auto& tr : tracks) {
            while (!tr.matchedFrames.empty() && tr.matchedFrames.front() < keep_from) {
                tr.matchedFrames.pop_front();
            }
        }

        if (k == last_k) {
            latest_det_to_track.swap(det_to_track);
        }
    }

    std::vector<GMTIDetection> current_targets;
    const auto& latest_dets = frame_dets.back();
    current_targets.reserve(latest_dets.size());

    for (size_t j = 0; j < latest_dets.size(); ++j) {
        const int tr_idx = (j < latest_det_to_track.size()) ? latest_det_to_track[j] : -1;
        if (tr_idx < 0 || tr_idx >= static_cast<int>(tracks.size())) {
            continue;
        }
        const SlidingTrack& tr = tracks[tr_idx];
        if (static_cast<int>(tr.matchedFrames.size()) < min_hits) {
            continue;
        }

        double lat = 0.0;
        double lon = 0.0;
        Gaussp3RV(latest_dets[j][0], latest_dets[j][1], cfg.L0, lat, lon);

        const int safe_id = std::max(1, std::min(65535, tr.id));
        GMTIDetection det{};
        det.id = static_cast<uint16_t>(safe_id);
        det.lon = lon;
        det.lat = lat;
        det.speed = tr.speed;
        det.direction = tr.direction;
        det.range = tr.range;
        det.utcMid = tr.utc;
        det.e = latest_dets[j][0];
        det.n = latest_dets[j][1];
        current_targets.push_back(det);
    }

    if (manager != nullptr) {
        current_targets = manager->update(cfg, current_targets);
    }

    std::vector<CurrentTargetPacket> track_packets;
    track_packets.reserve(current_targets.size());
    for (const auto& det : current_targets) {
        track_packets.push_back(CurrentTargetPacket{
            det.id,
            det.lon,
            det.lat,
            det.speed,
            det.direction,
            det.range,
            det.utcMid
        });
    }

    std::sort(track_packets.begin(), track_packets.end(),
              [](const CurrentTargetPacket& a, const CurrentTargetPacket& b) {
                  return a.id < b.id;
              });

    std::cout << "多周期关联完成，当前周期真目标数: " << track_packets.size()
              << " (关联文件数=" << effective_window << ", 判真阈值=" << min_hits << ")" << std::endl;

    char out_name[256];
    sprintf(out_name, "GMTI%02d_track.bin", idx_range.back());
    const std::string save_path = result_dir + "/" + out_name;
    printCurrentTargets(track_packets);
    dumpTrackPointsCsv(cfg, idx_range.back(), current_targets);
    writeCurrentCyclePackets(track_packets, save_path);

    return current_targets;
}

std::vector<GMTIDetection> trackModuleOnline(const Config& cfg, TrackManager* manager)
{
    if (!manager) {
        return trackModule(cfg);
    }

    const int result_id = resolveCurrentResultId(cfg);
    if (result_id <= 0) {
        std::cout << "[TRACK_ONLINE][WARN] 无法解析当前 GMTI 结果编号，回退到旧 trackModule。" << std::endl;
        return trackModule(cfg, manager);
    }

    const std::string input_path = makeGmtiBinPath(cfg.result_add, result_id);

    std::vector<GMTIDetection> dets;
    double utc_global = 0.0;
    const bool ok = readCurrentDetections(input_path, dets, utc_global);
    if (!ok) {
        dets.clear();
    }

    for (auto& d : dets) {
        double E = 0.0;
        double N = 0.0;
        Gaussp3(d.lat, d.lon, cfg.L0, E, N);
        d.e = E;
        d.n = N;
    }

    const double frame_utc = estimateFrameUtc(dets, utc_global);
    std::vector<GMTIDetection> outputs =
        manager->updateRawDetections(cfg, dets, result_id, frame_utc);

    for (auto& out : outputs) {
        double lat = 0.0;
        double lon = 0.0;
        Gaussp3RV(out.e, out.n, cfg.L0, lat, lon);
        out.lat = lat;
        out.lon = lon;
    }

    std::sort(outputs.begin(), outputs.end(),
              [](const GMTIDetection& a, const GMTIDetection& b) {
                  return a.id < b.id;
              });

    std::vector<CurrentTargetPacket> packets = makeCurrentPackets(outputs);
    std::cout << "[TRACK_ONLINE] 当前周期确认输出目标数: " << packets.size()
              << " (result_id=" << result_id
              << ", raw_dets=" << dets.size() << ")" << std::endl;
    printCurrentTargets(packets);
    dumpTrackPointsCsv(cfg, result_id, outputs);
    writeCurrentCyclePackets(packets, makeGmtiTrackPath(cfg.result_add, result_id));
    return outputs;
}

GMTIResult read_mt_result_bin(const std::string& binFile) {
    GMTIResult out;
    std::ifstream fs(binFile, std::ios::binary);
    if (!fs.is_open()) {
        throw std::runtime_error("Cannot open file: " + binFile);
    }

    fs.seekg(0, std::ios::end);
    const std::streamoff file_size = fs.tellg();
    fs.seekg(0, std::ios::beg);

    // 1. 读取目标总数 (uint8)
    uint16_t N_raw;
    fs.read(reinterpret_cast<char*>(&N_raw), 2);
    int N = std::min((int)N_raw, MAX_TGT);
    out.count = N;

    const std::streamoff new_size_36 = 2 + static_cast<std::streamoff>(N) * 36;
    const std::streamoff new_size_28 = 2 + static_cast<std::streamoff>(N) * 28 + 8;
    const std::streamoff legacy_size_28_no_tail = 2 + static_cast<std::streamoff>(N) * 28;
    const std::streamoff old_size_27 = 2 + static_cast<std::streamoff>(N) * 27 + 4;
    const bool use_new_layout_36 = (file_size == new_size_36);
    const bool use_new_layout_28 = (file_size == new_size_28);
    const bool use_legacy_layout_28_no_tail = (file_size == legacy_size_28_no_tail);
    const bool use_old_layout = (!use_new_layout_36 && !use_new_layout_28 && !use_legacy_layout_28_no_tail && file_size == old_size_27);

    // 2. 循环读取每个目标
    for (int i = 0; i < N; ++i) {
        GMTIDetection det{};
        int32_t lon_int = 0, lat_int = 0;

        if (use_new_layout_36) {
            fs.read(reinterpret_cast<char*>(&det.id), 2);
            fs.read(reinterpret_cast<char*>(&lon_int), 4);
            fs.read(reinterpret_cast<char*>(&lat_int), 4);
            uint16_t speed_q = 0;
            fs.read(reinterpret_cast<char*>(&speed_q), 2);
            det.speed = static_cast<double>(speed_q) * 0.01;
            fs.read(reinterpret_cast<char*>(&det.direction), 8);
            fs.read(reinterpret_cast<char*>(&det.range), 8);
            fs.read(reinterpret_cast<char*>(&det.utcMid), 8);
        } else if (use_new_layout_28 || use_legacy_layout_28_no_tail) {
            fs.read(reinterpret_cast<char*>(&det.id), 2);
            fs.read(reinterpret_cast<char*>(&lon_int), 4);
            fs.read(reinterpret_cast<char*>(&lat_int), 4);
            uint16_t speed_q = 0;
            fs.read(reinterpret_cast<char*>(&speed_q), 2);
            det.speed = static_cast<double>(speed_q) * 0.01;
            fs.read(reinterpret_cast<char*>(&det.direction), 8);
            fs.read(reinterpret_cast<char*>(&det.range), 8);
            det.utcMid = 0.0;
        } else if (use_old_layout) {
            uint16_t old_id = 0;
            fs.read(reinterpret_cast<char*>(&old_id), 2);
            det.id = old_id;
            fs.read(reinterpret_cast<char*>(&lon_int), 4);
            fs.read(reinterpret_cast<char*>(&lat_int), 4);
            float speed_f = 0.0f;
            fs.read(reinterpret_cast<char*>(&speed_f), 4);
            det.speed = static_cast<double>(speed_f);
            det.direction = 0.0;
            det.range = 0.0;
            det.utcMid = 0.0;
            fs.seekg(4, std::ios::cur); // 跳过 4 字节保留位
        } else {
            throw std::runtime_error("Unsupported GMTI binary layout: " + binFile);
        }

        det.lon = lon_int * LSB;
        det.lat = lat_int * LSB;
        out.targets.push_back(det);
    }

    // 3. 读取全局 utcMid（仅兼容旧布局，不参与新布局速度计算）
    if (use_new_layout_36) {
        out.utcGlobal = out.targets.empty() ? 0.0 : out.targets.front().utcMid;
    } else if (use_new_layout_28) {
        fs.seekg(2 + static_cast<std::streamoff>(N) * 28, std::ios::beg);
        fs.read(reinterpret_cast<char*>(&out.utcGlobal), 8);
        for (auto& det : out.targets) {
            det.utcMid = out.utcGlobal;
        }
    } else if (use_legacy_layout_28_no_tail) {
        out.utcGlobal = 0.0;
        for (auto& det : out.targets) {
            det.utcMid = 0.0;
        }
    } else if (use_old_layout) {
        fs.seekg(2 + static_cast<std::streamoff>(N) * 27, std::ios::beg);
        float old_utc = 0.0f;
        fs.read(reinterpret_cast<char*>(&old_utc), 4);
        out.utcGlobal = static_cast<double>(old_utc);
        for (auto& det : out.targets) {
            det.utcMid = out.utcGlobal;
        }
    } else {
        out.utcGlobal = 0.0;
    }

    fs.close();
    return out;
}

void writeCurrentCyclePackets(const std::vector<CurrentTargetPacket>& targets, const std::string& fileName) {
    std::ofstream ofs(fileName, std::ios::binary);
    if (!ofs.is_open()) return;

    constexpr size_t REC_BYTES = 28;
    std::vector<uint8_t> pkt;

    for (const auto& t : targets) {
        const size_t base = pkt.size();
        pkt.resize(base + REC_BYTES, 0u);
        put_u16_le(pkt, base + 0, t.id);
        put_i32_le(pkt, base + 2, quant_deg_to_i32(t.lon));
        put_i32_le(pkt, base + 6, quant_deg_to_i32(t.lat));
        put_u16_le(pkt, base + 10, quant_speed_to_u16(t.speed));
        put_f64_le(pkt, base + 12, wrap180_deg(t.direction));
        put_f64_le(pkt, base + 20, t.range);
    }

    ofs.write(reinterpret_cast<const char*>(pkt.data()), static_cast<std::streamsize>(pkt.size()));

    ofs.close();
    std::cout << "当前周期 GMTI 目标信息包已保存: " << fileName
              << "，目标数: " << targets.size() << std::endl;
}


// ===== main 函数 =====
// NOTE: main function moved to main.cpp for unified GMTI_core target
// This function is now called from main.cpp after GMTI detection/location phase
// 
// int main(int argc, char* argv[]) {
//     if (argc < 3) {
//         std::cout << "Usage: GMTI_track <result_dir> <idx1> <idx2> ... [--win N]" << std::endl;
//         return 1;
//     }
//
//     std::string result_dir = argv[1];
//     std::vector<int> idx_range;
//     int win_len = 5;
//
//     for (int i = 2; i < argc; ++i) {
//         const std::string arg = argv[i];
//         if (arg == "--win") {
//             if (i + 1 >= argc) {
//                 std::cerr << "--win 缺少窗长参数。" << std::endl;
//                 return 1;
//             }
//             win_len = std::stoi(argv[++i]);
//             continue;
//         }
//         if (arg.rfind("--win=", 0) == 0) {
//             win_len = std::stoi(arg.substr(6));
//             continue;
//         }
//         idx_range.push_back(std::stoi(arg));
//     }
//
//     if (idx_range.empty()) {
//         std::cerr << "未提供周期索引。" << std::endl;
//         return 1;
//     }
//
//     win_len = clamp(win_len, 5, 10);
//     if (static_cast<int>(idx_range.size()) < win_len) {
//         std::cout << "输入周期数少于窗长，自动使用有效窗长: " << idx_range.size() << std::endl;
//     }
//
//     trackModule(result_dir, idx_range, win_len);
//     return 0;
// }
