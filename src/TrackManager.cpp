#include "TrackManager.hpp"
#include "geo/geoProj.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

namespace {

constexpr double kInvalidCost = 1.0e12;
constexpr double kMissCost = 25.0;
constexpr double kMahalanobisGate = 16.0;
constexpr double kMinMeasurementStdM = 20.0;
constexpr double kProcessAccelMps2 = 5.0;
constexpr double kPi = 3.14159265358979323846;
constexpr double kSmallDistanceM = 1.0;
constexpr double kSmallSpeedMps = 0.5;

struct DetectionPoint {
    GMTIDetection det;
    double e = 0.0;
    double n = 0.0;
};

struct AssocScore {
    bool valid = false;
    double cost = kInvalidCost;
    double dist = 0.0;
    double instant_speed = 0.0;
    double mahalanobis_d2 = kInvalidCost;
};

static bool validUtc(double utc)
{
    return std::isfinite(utc) && utc > 0.0;
}

static double currentUtc(const std::vector<DetectionPoint>& detections)
{
    for (const auto& det : detections) {
        if (validUtc(det.det.utcMid)) {
            return det.det.utcMid;
        }
    }
    return 0.0;
}

static double safeDefaultDt(const Config& cfg)
{
    return std::max(1e-3, cfg.track_default_dt);
}

static double frameDtLegacy(const ManagedTrack& tr, double utc, int result_id)
{
    if (validUtc(utc) && validUtc(tr.utc)) {
        const double dt = utc - tr.utc;
        if (dt > 1e-3) {
            return dt;
        }
    }
    if (result_id > 0 && tr.last_result_id > 0) {
        return static_cast<double>(std::max(1, result_id - tr.last_result_id));
    }
    return 1.0;
}

static double computeDtForPredict(const Config& cfg, const ManagedTrack& tr, double frame_utc)
{
    if (validUtc(frame_utc) && validUtc(tr.utc) && frame_utc > tr.utc) {
        return frame_utc - tr.utc;
    }
    return safeDefaultDt(cfg);
}

static double computeDtForDetection(const Config& cfg,
                                    const ManagedTrack& tr,
                                    const GMTIDetection& det,
                                    double frame_utc)
{
    if (validUtc(det.utcMid) && validUtc(tr.utc) && det.utcMid > tr.utc) {
        return det.utcMid - tr.utc;
    }
    if (validUtc(frame_utc) && validUtc(tr.utc) && frame_utc > tr.utc) {
        return frame_utc - tr.utc;
    }
    return safeDefaultDt(cfg);
}

static void initCovarianceLegacy(ManagedTrack& tr, double gate_m)
{
    tr.P.fill(0.0);
    const double pos_var = std::max(kMinMeasurementStdM, gate_m * 0.5);
    const double vel_std = std::max(10.0, gate_m * 0.05);
    tr.P[0] = pos_var * pos_var;
    tr.P[5] = pos_var * pos_var;
    tr.P[10] = vel_std * vel_std;
    tr.P[15] = vel_std * vel_std;
}

static void initCovarianceOnline(ManagedTrack& tr)
{
    tr.P.fill(0.0);
    tr.P[0] = 100.0;
    tr.P[5] = 100.0;
    tr.P[10] = 1000.0;
    tr.P[15] = 1000.0;
}

static void kalmanPredict(ManagedTrack& tr, const Config& cfg, double dt)
{
    dt = std::max(1e-3, dt);
    tr.e += tr.ve * dt;
    tr.n += tr.vn * dt;

    const double p00 = tr.P[0],  p01 = tr.P[1],  p02 = tr.P[2],  p03 = tr.P[3];
    const double p10 = tr.P[4],  p11 = tr.P[5],  p12 = tr.P[6],  p13 = tr.P[7];
    const double p20 = tr.P[8],  p21 = tr.P[9],  p22 = tr.P[10], p23 = tr.P[11];
    const double p30 = tr.P[12], p31 = tr.P[13], p32 = tr.P[14], p33 = tr.P[15];

    tr.P[0] = p00 + dt * (p20 + p02) + dt * dt * p22 + cfg.track_process_noise_pos;
    tr.P[1] = p01 + dt * (p21 + p03) + dt * dt * p23;
    tr.P[2] = p02 + dt * p22;
    tr.P[3] = p03 + dt * p23;

    tr.P[4] = p10 + dt * (p30 + p12) + dt * dt * p32;
    tr.P[5] = p11 + dt * (p31 + p13) + dt * dt * p33 + cfg.track_process_noise_pos;
    tr.P[6] = p12 + dt * p32;
    tr.P[7] = p13 + dt * p33;

    tr.P[8] = p20 + dt * p22;
    tr.P[9] = p21 + dt * p23;
    tr.P[10] = p22 + cfg.track_process_noise_vel;
    tr.P[11] = p23;

    tr.P[12] = p30 + dt * p32;
    tr.P[13] = p31 + dt * p33;
    tr.P[14] = p32;
    tr.P[15] = p33 + cfg.track_process_noise_vel;
}

static void predictTrackLegacy(ManagedTrack& tr, double dt)
{
    Config cfg;
    cfg.track_process_noise_pos = 0.25 * dt * dt * dt * dt * kProcessAccelMps2 * kProcessAccelMps2;
    cfg.track_process_noise_vel = dt * dt * kProcessAccelMps2 * kProcessAccelMps2;
    kalmanPredict(tr, cfg, dt);
}

static double mahalanobis2(const ManagedTrack& tr,
                           double det_e,
                           double det_n,
                           double measurement_var,
                           double& dist_out)
{
    const double rx = det_e - tr.e;
    const double ry = det_n - tr.n;
    dist_out = std::sqrt(rx * rx + ry * ry);

    const double s00 = tr.P[0] + measurement_var;
    const double s01 = tr.P[1];
    const double s10 = tr.P[4];
    const double s11 = tr.P[5] + measurement_var;
    const double det_s = s00 * s11 - s01 * s10;
    if (std::fabs(det_s) < 1e-9) {
        return kInvalidCost;
    }

    const double inv00 = s11 / det_s;
    const double inv01 = -s01 / det_s;
    const double inv10 = -s10 / det_s;
    const double inv11 = s00 / det_s;
    return rx * (inv00 * rx + inv01 * ry) + ry * (inv10 * rx + inv11 * ry);
}

static void kalmanUpdate(ManagedTrack& tr,
                         double det_e,
                         double det_n,
                         double measurement_var,
                         double* residual_d2)
{
    const double y0 = det_e - tr.e;
    const double y1 = det_n - tr.n;

    const double s00 = tr.P[0] + measurement_var;
    const double s01 = tr.P[1];
    const double s10 = tr.P[4];
    const double s11 = tr.P[5] + measurement_var;
    const double det_s = s00 * s11 - s01 * s10;
    if (std::fabs(det_s) < 1e-9) {
        tr.e = det_e;
        tr.n = det_n;
        if (residual_d2) {
            *residual_d2 = 0.0;
        }
        return;
    }

    const double inv00 = s11 / det_s;
    const double inv01 = -s01 / det_s;
    const double inv10 = -s10 / det_s;
    const double inv11 = s00 / det_s;
    if (residual_d2) {
        *residual_d2 = y0 * (inv00 * y0 + inv01 * y1) + y1 * (inv10 * y0 + inv11 * y1);
    }

    double K[8] = {0.0};
    for (int r = 0; r < 4; ++r) {
        const double p_r0 = tr.P[r * 4 + 0];
        const double p_r1 = tr.P[r * 4 + 1];
        K[r * 2 + 0] = p_r0 * inv00 + p_r1 * inv10;
        K[r * 2 + 1] = p_r0 * inv01 + p_r1 * inv11;
    }

    tr.e += K[0] * y0 + K[1] * y1;
    tr.n += K[2] * y0 + K[3] * y1;
    tr.ve += K[4] * y0 + K[5] * y1;
    tr.vn += K[6] * y0 + K[7] * y1;

    double oldP[16];
    for (int i = 0; i < 16; ++i) {
        oldP[i] = tr.P[i];
    }
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            tr.P[r * 4 + c] = oldP[r * 4 + c]
                - K[r * 2 + 0] * oldP[c]
                - K[r * 2 + 1] * oldP[4 + c];
        }
    }
}

static double median(std::deque<double> values)
{
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const size_t mid = values.size() / 2;
    if (values.size() % 2 == 0) {
        return 0.5 * (values[mid - 1] + values[mid]);
    }
    return values[mid];
}

static double absAngleDiff(double a, double b)
{
    double d = std::fmod(a - b + kPi, 2.0 * kPi);
    if (d < 0.0) {
        d += 2.0 * kPi;
    }
    return std::fabs(d - kPi);
}

static double circularMean(const std::deque<double>& angles)
{
    if (angles.empty()) {
        return 0.0;
    }
    double s = 0.0;
    double c = 0.0;
    for (double a : angles) {
        s += std::sin(a);
        c += std::cos(a);
    }
    return std::atan2(s, c);
}

static void pushBounded(std::deque<int>& values, int value, int window)
{
    values.push_back(value);
    while (static_cast<int>(values.size()) > std::max(1, window)) {
        values.pop_front();
    }
}

static void pushBounded(std::deque<double>& values, double value, int window)
{
    if (!std::isfinite(value)) {
        return;
    }
    values.push_back(value);
    while (static_cast<int>(values.size()) > std::max(1, window)) {
        values.pop_front();
    }
}

static void pushPointHistory(ManagedTrack& tr, double e, double n, double utc, int result_id, const Config& cfg)
{
    tr.point_history.push_back(TrackPoint{e, n, utc, result_id});
    const int keep = std::max(3, cfg.track_linearity_window);
    while (static_cast<int>(tr.point_history.size()) > keep) {
        tr.point_history.pop_front();
    }
}

static int recentHits(const ManagedTrack& tr, int window)
{
    int count = 0;
    const int n = std::min<int>(std::max(1, window), tr.hit_history.size());
    for (int i = 0; i < n; ++i) {
        count += tr.hit_history[tr.hit_history.size() - 1 - i];
    }
    return count;
}

static double distanceEN(double e0, double n0, double e1, double n1)
{
    const double de = e1 - e0;
    const double dn = n1 - n0;
    return std::sqrt(de * de + dn * dn);
}

static double computeLinearity(const ManagedTrack& tr, int window)
{
    const int n = std::min<int>(std::max(3, window), tr.point_history.size());
    if (n < 3) {
        return 1.0;
    }
    const int start = static_cast<int>(tr.point_history.size()) - n;
    const TrackPoint& first = tr.point_history[start];
    const TrackPoint& last = tr.point_history.back();
    const double straight = distanceEN(first.e, first.n, last.e, last.n);
    double path = 0.0;
    for (int i = start + 1; i < static_cast<int>(tr.point_history.size()); ++i) {
        path += distanceEN(tr.point_history[i - 1].e, tr.point_history[i - 1].n,
                           tr.point_history[i].e, tr.point_history[i].n);
    }
    if (path <= 1e-6) {
        return 1.0;
    }
    return straight / path;
}

static bool linearityOK(const ManagedTrack& tr, const Config& cfg)
{
    if (tr.point_history.size() < 3) {
        return true;
    }
    return computeLinearity(tr, cfg.track_linearity_window) >= cfg.track_min_linearity_confirm;
}

static bool shouldPromoteToConfirmed(const Config& cfg, const ManagedTrack& tr)
{
    if (recentHits(tr, cfg.track_confirm_window) < std::max(1, cfg.track_confirm_hits)) {
        return false;
    }
    return linearityOK(tr, cfg);
}

static bool shouldOutputThisFrame(const Config& cfg, const ManagedTrack& tr, int result_id)
{
    if (tr.state != TrackState::Confirmed) {
        return false;
    }
    if (!tr.matched_this_frame) {
        return false;
    }
    if (tr.last_result_id != result_id) {
        return false;
    }
    return recentHits(tr, cfg.track_confirm_window) >= std::max(1, cfg.track_confirm_hits);
}

static const char* stateName(TrackState state)
{
    switch (state) {
    case TrackState::Tentative:
        return "Tentative";
    case TrackState::Confirmed:
        return "Confirmed";
    case TrackState::Coasted:
        return "Coasted";
    case TrackState::Deleted:
        return "Deleted";
    }
    return "Unknown";
}

static void updateMotionHistoryFromEN(ManagedTrack& tr,
                                      const GMTIDetection& det,
                                      double dt,
                                      const Config& cfg)
{
    if (!tr.point_history.empty()) {
        const TrackPoint& last = tr.point_history.back();
        dt = std::max(1e-3, dt);
        const double de = det.e - last.e;
        const double dn = det.n - last.n;
        const double dist = std::sqrt(de * de + dn * dn);
        tr.ve = de / dt;
        tr.vn = dn / dt;
        tr.speed = std::sqrt(tr.ve * tr.ve + tr.vn * tr.vn);
        pushBounded(tr.speed_history, tr.speed, cfg.track_linearity_window);
        if (dist >= kSmallDistanceM && tr.speed >= kSmallSpeedMps) {
            pushBounded(tr.heading_history, std::atan2(dn, de), cfg.track_linearity_window);
        }
    }
}

static AssocScore computeAssocScore(const Config& cfg,
                                    const ManagedTrack& tr,
                                    const GMTIDetection& det,
                                    double frame_utc,
                                    double measurement_var)
{
    AssocScore score;
    score.mahalanobis_d2 = mahalanobis2(tr, det.e, det.n, measurement_var, score.dist);
    double gate_m = cfg.track_gate_m;
    double chi2_gate = cfg.track_chi2_gate;
    if (tr.state == TrackState::Tentative) {
        gate_m *= std::max(1.0, cfg.track_tentative_gate_scale);
        chi2_gate *= std::max(1.0, cfg.track_tentative_chi2_scale);
    }
    if (gate_m > 0.0 && score.dist > gate_m) {
        return score;
    }
    if (score.mahalanobis_d2 > chi2_gate) {
        return score;
    }

    const double dt = computeDtForDetection(cfg, tr, det, frame_utc);
    const TrackPoint* last_hit = tr.point_history.empty() ? nullptr : &tr.point_history.back();
    const double motion_dist = last_hit ? distanceEN(last_hit->e, last_hit->n, det.e, det.n) : score.dist;
    score.instant_speed = motion_dist / std::max(1e-3, dt);
    if (cfg.track_v_max > 0.0 && score.instant_speed > cfg.track_v_max) {
        return score;
    }

    const bool early_track = (tr.hit_count < 2 || tr.point_history.size() < 2);
    double speed_penalty = 0.0;
    if (!early_track && tr.speed_history.size() >= 2) {
        const double v_ref = median(tr.speed_history);
        speed_penalty = std::fabs(score.instant_speed - v_ref) / std::max(v_ref, 1.0);
    }

    double heading_penalty = 0.0;
    if (!early_track && last_hit && tr.heading_history.size() >= 2 &&
        motion_dist >= kSmallDistanceM && tr.speed >= kSmallSpeedMps) {
        const double heading_inst = std::atan2(det.n - last_hit->n, det.e - last_hit->e);
        const double heading_ref = circularMean(tr.heading_history);
        heading_penalty = absAngleDiff(heading_inst, heading_ref) / kPi;
    }

    score.valid = true;
    score.cost = score.mahalanobis_d2
        + cfg.track_speed_smooth_weight * speed_penalty
        + cfg.track_heading_weight * heading_penalty;
    return score;
}

static std::vector<int> solveAssignment(const std::vector<std::vector<double>>& cost)
{
    const int n = static_cast<int>(cost.size());
    if (n == 0) {
        return std::vector<int>();
    }
    const int m = static_cast<int>(cost.front().size());
    if (m == 0 || n > m) {
        return std::vector<int>(n, -1);
    }

    std::vector<double> u(n + 1, 0.0), v(m + 1, 0.0);
    std::vector<int> p(m + 1, 0), way(m + 1, 0);

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(m + 1, kInvalidCost);
        std::vector<char> used(m + 1, 0);
        do {
            used[j0] = 1;
            const int i0 = p[j0];
            double delta = kInvalidCost;
            int j1 = 0;
            for (int j = 1; j <= m; ++j) {
                if (used[j]) {
                    continue;
                }
                const double cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                if (cur < minv[j]) {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1 = j;
                }
            }
            if (delta >= kInvalidCost * 0.5) {
                break;
            }
            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);

        do {
            const int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    std::vector<int> assignment(n, -1);
    for (int j = 1; j <= m; ++j) {
        if (p[j] > 0 && p[j] <= n) {
            assignment[p[j] - 1] = j - 1;
        }
    }
    return assignment;
}

static bool ensureDirectory(const std::string& path)
{
    if (path.empty()) {
        return false;
    }

    std::string current;
    current.reserve(path.size());
    for (size_t i = 0; i < path.size(); ++i) {
        current.push_back(path[i]);
        if (path[i] != '/' && i + 1 != path.size()) {
            continue;
        }
        while (current.size() > 1 && current.back() == '/') {
            current.pop_back();
        }
        if (current.empty() || current == "/") {
            if (path[i] == '/' && current != "/") {
                current.push_back('/');
            }
            continue;
        }
        if (::mkdir(current.c_str(), 0777) != 0 && errno != EEXIST) {
            std::cerr << "[TRACK][WARN] cannot create debug dir: "
                      << current << " error=" << std::strerror(errno) << std::endl;
            return false;
        }
        if (i + 1 != path.size() && path[i] == '/' && current.back() != '/') {
            current.push_back('/');
        }
    }
    return true;
}

static std::string debugDir(const Config& cfg)
{
    if (!cfg.track_debug_dir.empty()) {
        return cfg.track_debug_dir;
    }
    if (cfg.result_add.empty()) {
        return "track_debug";
    }
    return cfg.result_add + "/track_debug";
}

static bool fileNeedsHeader(const std::string& path)
{
    std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
    return !in.is_open() || in.tellg() <= 0;
}

static std::ofstream openCsvAppend(const std::string& path, const char* header)
{
    const bool need_header = fileNeedsHeader(path);
    std::ofstream out(path.c_str(), std::ios::out | std::ios::app);
    if (!out.is_open()) {
        std::cerr << "[TRACK][WARN] cannot open debug csv: " << path << std::endl;
        return out;
    }
    if (need_header) {
        out << header << "\n";
    }
    out << std::fixed;
    return out;
}

} // namespace

std::vector<GMTIDetection> TrackManager::update(
    const Config& cfg,
    const std::vector<GMTIDetection>& current_targets)
{
    const int result_id = currentResultId(cfg);
    if (result_id > 0 && last_update_result_id_ > 0 && result_id < last_update_result_id_) {
        std::cout << "[TRACK][WARN] TrackManager result id moved backward: current="
                  << result_id << ", previous=" << last_update_result_id_
                  << ". Reset persistent tracks." << std::endl;
        reset();
    } else if (result_id > 0 && result_id == last_update_result_id_) {
        std::cout << "[TRACK][WARN] TrackManager received duplicate result id: "
                  << result_id << ". Persistent ids will be updated in place." << std::endl;
    }

    std::vector<DetectionPoint> detections;
    detections.reserve(current_targets.size());
    for (const auto& det : current_targets) {
        DetectionPoint pt;
        pt.det = det;
        Gaussp3(det.lat, det.lon, cfg.L0, pt.e, pt.n);
        detections.push_back(pt);
    }

    const double gate_m = std::max(0.0, cfg.track_gate_m);
    const double v_max = std::max(0.0, cfg.track_v_max);
    const int max_miss = maxMiss(cfg);
    const double measurement_std = std::max(kMinMeasurementStdM, gate_m * 0.25);
    const double measurement_var = measurement_std * measurement_std;
    const double update_utc = currentUtc(detections);

    std::vector<int> active_tracks;
    active_tracks.reserve(tracks_.size());
    for (size_t ti = 0; ti < tracks_.size(); ++ti) {
        ManagedTrack& tr = tracks_[ti];
        if (tr.state == TrackState::Deleted || tr.miss_count > max_miss) {
            continue;
        }
        if (result_id > 0 && tr.last_result_id > 0 &&
            result_id - tr.last_result_id > max_miss + 1) {
            continue;
        }
        if (tr.hit_count <= 0) {
            initCovarianceLegacy(tr, gate_m);
        }
        predictTrackLegacy(tr, frameDtLegacy(tr, update_utc, result_id));
        active_tracks.push_back(static_cast<int>(ti));
    }

    std::vector<char> track_used(tracks_.size(), 0);
    std::vector<char> det_used(detections.size(), 0);
    const int row_count = static_cast<int>(active_tracks.size());
    const int det_count = static_cast<int>(detections.size());
    if (row_count > 0) {
        const int col_count = det_count + row_count;
        std::vector<std::vector<double>> cost(row_count, std::vector<double>(col_count, kMissCost));

        for (int row = 0; row < row_count; ++row) {
            const ManagedTrack& tr = tracks_[active_tracks[row]];
            for (int di = 0; di < det_count; ++di) {
                double dist = 0.0;
                const double d2 = mahalanobis2(tr, detections[di].e, detections[di].n,
                                               measurement_var, dist);
                if (gate_m > 0.0 && dist > gate_m) {
                    cost[row][di] = kInvalidCost;
                    continue;
                }

                const double dt = frameDtLegacy(tr, detections[di].det.utcMid, result_id);
                const double speed = dist / std::max(1e-3, dt);
                if (v_max > 0.0 && speed > v_max) {
                    cost[row][di] = kInvalidCost;
                    continue;
                }
                if (d2 > kMahalanobisGate) {
                    cost[row][di] = kInvalidCost;
                    continue;
                }

                const double norm_dist = (gate_m > 0.0) ? dist / gate_m : 0.0;
                cost[row][di] = d2 + 0.5 * norm_dist;
            }
        }

        const std::vector<int> assignment = solveAssignment(cost);
        for (int row = 0; row < row_count; ++row) {
            const int ti = active_tracks[row];
            const int col = (row < static_cast<int>(assignment.size())) ? assignment[row] : -1;
            if (col < 0 || col >= det_count || cost[row][col] >= kMissCost) {
                continue;
            }

            ManagedTrack& tr = tracks_[ti];
            const DetectionPoint& pt = detections[col];
            const double dt = frameDtLegacy(tr, pt.det.utcMid, result_id);
            const double prev_e = tr.e;
            const double prev_n = tr.n;
            kalmanUpdate(tr, pt.e, pt.n, measurement_var, nullptr);
            const double dist = distanceEN(pt.e, pt.n, prev_e, prev_n);
            const double speed = dist / std::max(1e-3, dt);

            detections[col].det.id = tr.id;
            tr.e = pt.e;
            tr.n = pt.n;
            tr.utc = pt.det.utcMid;
            tr.speed = pt.det.speed > 0.0 ? pt.det.speed : speed;
            tr.direction = pt.det.direction;
            tr.range = pt.det.range;
            tr.last_result_id = result_id;
            tr.hit_count += 1;
            tr.miss_count = 0;
            tr.state = TrackState::Confirmed;
            tr.recent_result_ids.push_back(result_id);
            while (tr.recent_result_ids.size() > static_cast<size_t>(cfg.track_idx_window)) {
                tr.recent_result_ids.pop_front();
            }

            track_used[ti] = 1;
            det_used[col] = 1;
        }
    }

    for (size_t ti = 0; ti < tracks_.size(); ++ti) {
        if (tracks_[ti].state == TrackState::Deleted) {
            continue;
        }
        if (ti >= track_used.size() || !track_used[ti]) {
            tracks_[ti].miss_count += 1;
            tracks_[ti].state = tracks_[ti].miss_count > max_miss
                ? TrackState::Deleted
                : TrackState::Coasted;
        }
    }

    for (size_t di = 0; di < detections.size(); ++di) {
        if (det_used[di]) {
            continue;
        }
        ManagedTrack tr;
        tr.id = allocateId();
        tr.state = TrackState::Confirmed;
        tr.e = detections[di].e;
        tr.n = detections[di].n;
        tr.ve = 0.0;
        tr.vn = 0.0;
        initCovarianceLegacy(tr, gate_m);
        tr.utc = detections[di].det.utcMid;
        tr.speed = detections[di].det.speed;
        tr.direction = detections[di].det.direction;
        tr.range = detections[di].det.range;
        tr.last_result_id = result_id;
        tr.hit_count = 1;
        tr.miss_count = 0;
        tr.recent_result_ids.push_back(result_id);
        detections[di].det.id = tr.id;
        tracks_.push_back(tr);
    }

    pruneDeleted();
    if (result_id > 0) {
        last_update_result_id_ = result_id;
    }

    std::vector<GMTIDetection> out;
    out.reserve(detections.size());
    for (const auto& pt : detections) {
        out.push_back(pt.det);
    }
    std::sort(out.begin(), out.end(), [](const GMTIDetection& a, const GMTIDetection& b) {
        return a.id < b.id;
    });
    return out;
}

std::vector<GMTIDetection> TrackManager::updateRawDetections(
    const Config& cfg,
    const std::vector<GMTIDetection>& current_detections,
    int result_id,
    double frame_utc)
{
    if (result_id > 0 && last_update_result_id_ > 0 && result_id < last_update_result_id_) {
        std::cout << "[TRACK][WARN] TrackManager result id moved backward: current="
                  << result_id << ", previous=" << last_update_result_id_
                  << ". Reset persistent tracks." << std::endl;
        reset();
    }

    for (auto& tr : tracks_) {
        tr.matched_this_frame = false;
        tr.matched_det_index = -1;
        tr.is_output_this_frame = false;
    }

    std::vector<GMTIDetection> dets = current_detections;
    for (auto& tr : tracks_) {
        if (tr.state == TrackState::Deleted) {
            continue;
        }
        kalmanPredict(tr, cfg, computeDtForPredict(cfg, tr, frame_utc));
    }

    std::vector<int> active_tracks;
    active_tracks.reserve(tracks_.size());
    for (size_t ti = 0; ti < tracks_.size(); ++ti) {
        if (tracks_[ti].state != TrackState::Deleted) {
            active_tracks.push_back(static_cast<int>(ti));
        }
    }

    std::vector<char> det_used(dets.size(), 0);
    std::vector<int> det_to_track_id(dets.size(), -1);
    std::vector<AssocScore> accepted_scores(tracks_.size());
    int num_matched_tracks = 0;
    int num_new_tracks = 0;
    const int row_count = static_cast<int>(active_tracks.size());
    const int det_count = static_cast<int>(dets.size());
    const double dummy_cost = std::max(1.0, cfg.track_dummy_cost);
    const double measurement_std = std::max(1.0, cfg.track_measurement_noise_pos);
    const double measurement_var = measurement_std * measurement_std;

    if (row_count > 0 && det_count > 0) {
        const int col_count = det_count + row_count;
        std::vector<std::vector<double>> cost(row_count, std::vector<double>(col_count, dummy_cost));
        std::vector<std::vector<AssocScore>> scores(row_count, std::vector<AssocScore>(det_count));

        for (int row = 0; row < row_count; ++row) {
            const ManagedTrack& tr = tracks_[active_tracks[row]];
            for (int di = 0; di < det_count; ++di) {
                scores[row][di] = computeAssocScore(cfg, tr, dets[di], frame_utc, measurement_var);
                if (scores[row][di].valid) {
                    cost[row][di] = scores[row][di].cost;
                } else if (cfg.track_debug_level >= 3) {
                    std::cout << "[TRACK_ONLINE] cost invalid track=" << tr.id
                              << " det=" << di
                              << " d=" << scores[row][di].dist
                              << " d2=" << scores[row][di].mahalanobis_d2
                              << " v=" << scores[row][di].instant_speed << std::endl;
                }
            }
        }

        const std::vector<int> assignment = solveAssignment(cost);
        for (int row = 0; row < row_count; ++row) {
            const int ti = active_tracks[row];
            const int col = (row < static_cast<int>(assignment.size())) ? assignment[row] : -1;
            if (col < 0 || col >= det_count || cost[row][col] >= dummy_cost) {
                continue;
            }

            ManagedTrack& tr = tracks_[ti];
            GMTIDetection& det = dets[col];
            const AssocScore& score = scores[row][col];
            if (!score.valid) {
                continue;
            }

            const double dt = computeDtForDetection(cfg, tr, det, frame_utc);
            double residual_d2 = 0.0;
            kalmanUpdate(tr, det.e, det.n, measurement_var, &residual_d2);
            updateMotionHistoryFromEN(tr, det, dt, cfg);

            tr.e = det.e;
            tr.n = det.n;
            tr.matched_this_frame = true;
            tr.matched_det_index = col;
            tr.age++;
            tr.hit_count++;
            tr.miss_count = 0;
            tr.consecutive_hits++;
            tr.last_result_id = result_id;
            tr.utc = validUtc(det.utcMid) ? det.utcMid : frame_utc;
            tr.range = det.range;
            tr.direction = det.direction;
            pushBounded(tr.hit_history, 1, cfg.track_confirm_window);
            pushPointHistory(tr, det.e, det.n, tr.utc, result_id, cfg);
            pushBounded(tr.residual_history, residual_d2, cfg.track_linearity_window);
            if (tr.state == TrackState::Coasted) {
                tr.state = TrackState::Confirmed;
            }

            accepted_scores[ti] = score;
            det_used[col] = 1;
            det_to_track_id[col] = tr.id;
            num_matched_tracks++;

            if (cfg.track_debug_level >= 2) {
                std::cout << "[TRACK_ONLINE] assign track=" << tr.id
                          << " det=" << col
                          << " cost=" << score.cost
                          << " d=" << score.dist
                          << " v=" << score.instant_speed
                          << " state=" << (tr.state == TrackState::Confirmed ? "Confirmed" : "Tentative")
                          << std::endl;
            }
        }
    }

    const int max_missed = std::max(0, cfg.track_max_missed);
    const int tentative_max_missed = std::max(0, cfg.track_tentative_max_missed);
    for (auto& tr : tracks_) {
        if (tr.state == TrackState::Deleted || tr.matched_this_frame) {
            continue;
        }
        tr.age++;
        tr.miss_count++;
        tr.consecutive_hits = 0;
        pushBounded(tr.hit_history, 0, cfg.track_confirm_window);

        if (tr.state == TrackState::Tentative) {
            if (tr.miss_count > tentative_max_missed) {
                tr.state = TrackState::Deleted;
                if (cfg.track_debug_level >= 2) {
                    std::cout << "[TRACK_ONLINE] delete tentative id=" << tr.id
                              << " miss=" << tr.miss_count << std::endl;
                }
            }
            continue;
        }

        if (tr.state == TrackState::Confirmed || tr.state == TrackState::Coasted) {
            if (tr.miss_count > max_missed) {
                tr.state = TrackState::Deleted;
                if (cfg.track_debug_level >= 2) {
                    std::cout << "[TRACK_ONLINE] delete id=" << tr.id
                              << " miss=" << tr.miss_count << std::endl;
                }
            } else if (tr.miss_count > 0) {
                tr.state = TrackState::Coasted;
                if (cfg.track_debug_level >= 2) {
                    std::cout << "[TRACK_ONLINE] coast id=" << tr.id
                              << " miss=" << tr.miss_count << std::endl;
                }
            }
        } else if (tr.miss_count > max_missed) {
            tr.state = TrackState::Deleted;
            if (cfg.track_debug_level >= 2) {
                std::cout << "[TRACK_ONLINE] delete id=" << tr.id
                          << " miss=" << tr.miss_count << std::endl;
            }
        }
    }

    for (size_t di = 0; di < dets.size(); ++di) {
        if (det_used[di]) {
            continue;
        }
        const GMTIDetection& det = dets[di];
        ManagedTrack tr;
        tr.id = allocateNextId();
        tr.state = TrackState::Tentative;
        tr.e = det.e;
        tr.n = det.n;
        tr.ve = 0.0;
        tr.vn = 0.0;
        initCovarianceOnline(tr);
        tr.utc = validUtc(det.utcMid) ? det.utcMid : frame_utc;
        tr.speed = 0.0;
        tr.direction = det.direction;
        tr.range = det.range;
        tr.last_result_id = result_id;
        tr.age = 1;
        tr.hit_count = 1;
        tr.miss_count = 0;
        tr.consecutive_hits = 1;
        tr.matched_this_frame = true;
        tr.matched_det_index = static_cast<int>(di);
        tr.hit_history.push_back(1);
        tr.point_history.push_back(TrackPoint{det.e, det.n, tr.utc, result_id});
        tracks_.push_back(tr);
        num_new_tracks++;

        if (cfg.track_debug_level >= 2) {
            std::cout << "[TRACK_ONLINE] new tentative id=" << tr.id
                      << " det=" << di << std::endl;
        }
    }

    for (auto& tr : tracks_) {
        if (tr.state == TrackState::Tentative && shouldPromoteToConfirmed(cfg, tr)) {
            tr.state = TrackState::Confirmed;
            if (cfg.track_debug_level >= 2) {
                std::cout << "[TRACK_ONLINE] confirm id=" << tr.id
                          << " hits=" << recentHits(tr, cfg.track_confirm_window)
                          << "/" << cfg.track_confirm_window
                          << " linearity=" << computeLinearity(tr, cfg.track_linearity_window)
                          << std::endl;
            }
        }
    }

    std::vector<GMTIDetection> out;
    for (auto& tr : tracks_) {
        tr.is_output_this_frame = false;
    }
    for (const auto& tr : tracks_) {
        if (!shouldOutputThisFrame(cfg, tr, result_id)) {
            continue;
        }
        GMTIDetection det{};
        det.id = tr.id;
        det.e = tr.e;
        det.n = tr.n;
        det.speed = tr.speed;
        det.direction = tr.direction;
        det.range = tr.range;
        det.utcMid = tr.utc;
        out.push_back(det);
    }
    for (const auto& det : out) {
        ManagedTrack* tr = findTrackById(det.id);
        if (tr) {
            tr->is_output_this_frame = true;
        }
    }

    const int num_unmatched_detections = static_cast<int>(
        std::count(det_used.begin(), det_used.end(), 0));
    const int log_level = std::max(cfg.track_debug_level, cfg.track_debug_dump_level);
    if (log_level >= 1) {
        int num_tentative = 0;
        int num_confirmed = 0;
        int num_coasted = 0;
        int num_deleted = 0;
        for (const auto& tr : tracks_) {
            if (tr.state == TrackState::Tentative) ++num_tentative;
            if (tr.state == TrackState::Confirmed) ++num_confirmed;
            if (tr.state == TrackState::Coasted) ++num_coasted;
            if (tr.state == TrackState::Deleted) ++num_deleted;
        }
        const int active_count = static_cast<int>(tracks_.size()) - num_deleted;
        std::cout << "[TRACK] result_id=" << result_id
                  << " dets=" << dets.size()
                  << " active=" << active_count
                  << " tentative=" << num_tentative
                  << " confirmed=" << num_confirmed
                  << " coasted=" << num_coasted
                  << " deleted=" << num_deleted
                  << " outputs=" << out.size()
                  << " new=" << num_new_tracks
                  << " matched=" << num_matched_tracks
                  << " unmatched_dets=" << num_unmatched_detections
                  << std::endl;
        if (log_level >= 2) {
            for (const auto& tr : tracks_) {
                std::cout << "[TRACK] id=" << tr.id
                          << " state=" << stateName(tr.state)
                          << " hit=" << tr.hit_count
                          << " miss=" << tr.miss_count
                          << " recent=" << recentHits(tr, cfg.track_confirm_window)
                          << "/" << cfg.track_confirm_window
                          << " matched=" << (tr.matched_this_frame ? 1 : 0)
                          << " det=" << tr.matched_det_index
                          << " speed=" << tr.speed
                          << std::endl;
            }
        }
    }

    dumpDebugSnapshot(cfg, result_id, frame_utc, dets, det_to_track_id,
                      static_cast<int>(out.size()), num_new_tracks,
                      num_matched_tracks, num_unmatched_detections);

    pruneDeleted();
    if (result_id > 0) {
        last_update_result_id_ = result_id;
    }

    std::sort(out.begin(), out.end(), [](const GMTIDetection& a, const GMTIDetection& b) {
        return a.id < b.id;
    });
    return out;
}

void TrackManager::reset()
{
    tracks_.clear();
    next_id_ = 1;
    last_update_result_id_ = 0;
}

uint16_t TrackManager::allocateId()
{
    return allocateNextId();
}

uint16_t TrackManager::allocateNextId()
{
    uint32_t candidate = next_id_;
    for (uint32_t count = 0; count < 65535U; ++count) {
        if (candidate == 0U || candidate > 65535U) {
            candidate = 1U;
        }
        const uint16_t id = static_cast<uint16_t>(candidate);
        if (!isActiveId(id)) {
            next_id_ = static_cast<uint16_t>((candidate >= 65535U) ? 1U : candidate + 1U);
            return id;
        }
        ++candidate;
    }

    std::cout << "[TRACK][WARN] TrackManager id pool exhausted. Reusing 65535." << std::endl;
    return 65535;
}

bool TrackManager::isActiveId(uint16_t id) const
{
    for (const auto& tr : tracks_) {
        if (tr.state != TrackState::Deleted && tr.id == id) {
            return true;
        }
    }
    return false;
}

ManagedTrack* TrackManager::findTrackById(uint16_t id)
{
    for (auto& tr : tracks_) {
        if (tr.id == id) {
            return &tr;
        }
    }
    return nullptr;
}

int TrackManager::currentResultId(const Config& cfg) const
{
    if (cfg.result_file_id > 0) {
        return cfg.result_file_id;
    }
    if (!cfg.track_idx_range.empty()) {
        return cfg.track_idx_range.back();
    }
    return 0;
}

int TrackManager::maxMiss(const Config& cfg) const
{
    const int window = std::max(1, cfg.track_idx_window);
    const int truth = std::max(1, cfg.track_truth_threshold);
    return std::max(1, window - truth + 1);
}

void TrackManager::dumpDebugSnapshot(const Config& cfg,
                                     int result_id,
                                     double frame_utc,
                                     const std::vector<GMTIDetection>& detections,
                                     const std::vector<int>& det_to_track_id,
                                     int num_outputs,
                                     int num_new_tracks,
                                     int num_matched_tracks,
                                     int num_unmatched_detections) const
{
    if (!cfg.track_debug_dump) {
        return;
    }

    const std::string dir = debugDir(cfg);
    if (!ensureDirectory(dir)) {
        return;
    }

    int num_tentative = 0;
    int num_confirmed = 0;
    int num_coasted = 0;
    int num_deleted = 0;
    for (const auto& tr : tracks_) {
        switch (tr.state) {
        case TrackState::Tentative:
            ++num_tentative;
            break;
        case TrackState::Confirmed:
            ++num_confirmed;
            break;
        case TrackState::Coasted:
            ++num_coasted;
            break;
        case TrackState::Deleted:
            ++num_deleted;
            break;
        }
    }

    {
        std::ofstream frames = openCsvAppend(
            dir + "/track_frames.csv",
            "result_id,frame_utc,num_detections,num_tracks,num_tentative,num_confirmed,num_coasted,num_deleted,num_outputs,num_new_tracks,num_matched_tracks,num_unmatched_detections");
        if (frames.is_open()) {
            frames << std::setprecision(12)
                   << result_id << ","
                   << frame_utc << ","
                   << detections.size() << ","
                   << tracks_.size() << ","
                   << num_tentative << ","
                   << num_confirmed << ","
                   << num_coasted << ","
                   << num_deleted << ","
                   << num_outputs << ","
                   << num_new_tracks << ","
                   << num_matched_tracks << ","
                   << num_unmatched_detections << "\n";
        }
    }

    {
        std::ofstream dets = openCsvAppend(
            dir + "/track_detections.csv",
            "result_id,det_index,matched,matched_track_id,e,n,lat,lon,utc,range,direction");
        if (dets.is_open()) {
            for (size_t i = 0; i < detections.size(); ++i) {
                const int matched_track_id =
                    (i < det_to_track_id.size()) ? det_to_track_id[i] : -1;
                const bool matched = matched_track_id >= 0;
                dets << std::setprecision(12)
                     << result_id << ","
                     << i << ","
                     << (matched ? 1 : 0) << ","
                     << matched_track_id << ","
                     << std::setprecision(6) << detections[i].e << ","
                     << detections[i].n << ","
                     << std::setprecision(12) << detections[i].lat << ","
                     << detections[i].lon << ","
                     << detections[i].utcMid << ","
                     << detections[i].range << ","
                     << detections[i].direction << "\n";
            }
        }
    }

    {
        std::ofstream states = openCsvAppend(
            dir + "/track_states.csv",
            "result_id,track_id,state,matched_this_frame,matched_det_index,e,n,lat,lon,ve,vn,speed,heading,utc,range,direction,age,hit_count,miss_count,consecutive_hits,recent_hits,linearity,is_output");
        if (states.is_open()) {
            for (const auto& tr : tracks_) {
                double lat = 0.0;
                double lon = 0.0;
                Gaussp3RV(tr.e, tr.n, cfg.L0, lat, lon);
                const double heading = std::atan2(tr.vn, tr.ve);
                const int recent = recentHits(tr, cfg.track_confirm_window);
                const double linearity = computeLinearity(tr, cfg.track_linearity_window);
                states << std::setprecision(12)
                       << result_id << ","
                       << tr.id << ","
                       << stateName(tr.state) << ","
                       << (tr.matched_this_frame ? 1 : 0) << ","
                       << tr.matched_det_index << ","
                       << std::setprecision(6) << tr.e << ","
                       << tr.n << ","
                       << std::setprecision(12) << lat << ","
                       << lon << ","
                       << tr.ve << ","
                       << tr.vn << ","
                       << tr.speed << ","
                       << heading << ","
                       << tr.utc << ","
                       << tr.range << ","
                       << tr.direction << ","
                       << tr.age << ","
                       << tr.hit_count << ","
                       << tr.miss_count << ","
                       << tr.consecutive_hits << ","
                       << recent << ","
                       << linearity << ","
                       << (tr.is_output_this_frame ? 1 : 0) << "\n";
            }
        }
    }
}

void TrackManager::pruneDeleted()
{
    tracks_.erase(std::remove_if(tracks_.begin(), tracks_.end(),
                                 [](const ManagedTrack& tr) {
                                     return tr.state == TrackState::Deleted;
                                 }),
                  tracks_.end());
}
