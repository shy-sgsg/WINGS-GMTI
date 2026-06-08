#include "TrackManager.hpp"
#include "geo/geoProj.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

constexpr double kInvalidCost = 1.0e12;
constexpr double kMissCost = 25.0;
constexpr double kMahalanobisGate = 16.0;
constexpr double kMinMeasurementStdM = 20.0;
constexpr double kProcessAccelMps2 = 5.0;

struct DetectionPoint {
    GMTIDetection det;
    double e = 0.0;
    double n = 0.0;
};

static bool validUtc(double utc)
{
    return std::isfinite(utc) && utc > 0.0;
}

static double wrap180(double angle_deg)
{
    angle_deg = std::fmod(angle_deg + 180.0, 360.0);
    if (angle_deg < 0.0) {
        angle_deg += 360.0;
    }
    return angle_deg - 180.0;
}

static double directionCost(double a, double b)
{
    if (!std::isfinite(a) || !std::isfinite(b)) {
        return 0.0;
    }
    const double diff = std::fabs(wrap180(a - b));
    return std::min(4.0, diff / 45.0);
}

static double rangeCost(double a, double b)
{
    if (!std::isfinite(a) || !std::isfinite(b) || a <= 0.0 || b <= 0.0) {
        return 0.0;
    }
    const double denom = std::max(100.0, 0.5 * (std::fabs(a) + std::fabs(b)));
    return std::min(4.0, std::fabs(a - b) / denom);
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

static double frameDt(const ManagedTrack& tr, double utc, int result_id)
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

static void initCovariance(ManagedTrack& tr, double gate_m)
{
    tr.P.fill(0.0);
    const double pos_var = std::max(kMinMeasurementStdM, gate_m * 0.5);
    const double vel_std = std::max(10.0, gate_m * 0.05);
    tr.P[0] = pos_var * pos_var;
    tr.P[5] = pos_var * pos_var;
    tr.P[10] = vel_std * vel_std;
    tr.P[15] = vel_std * vel_std;
}

static void predictTrack(ManagedTrack& tr, double dt)
{
    dt = std::max(1e-3, dt);
    tr.e += tr.ve * dt;
    tr.n += tr.vn * dt;

    const double p00 = tr.P[0],  p01 = tr.P[1],  p02 = tr.P[2],  p03 = tr.P[3];
    const double p10 = tr.P[4],  p11 = tr.P[5],  p12 = tr.P[6],  p13 = tr.P[7];
    const double p20 = tr.P[8],  p21 = tr.P[9],  p22 = tr.P[10], p23 = tr.P[11];
    const double p30 = tr.P[12], p31 = tr.P[13], p32 = tr.P[14], p33 = tr.P[15];

    tr.P[0] = p00 + dt * (p20 + p02) + dt * dt * p22;
    tr.P[1] = p01 + dt * (p21 + p03) + dt * dt * p23;
    tr.P[2] = p02 + dt * p22;
    tr.P[3] = p03 + dt * p23;

    tr.P[4] = p10 + dt * (p30 + p12) + dt * dt * p32;
    tr.P[5] = p11 + dt * (p31 + p13) + dt * dt * p33;
    tr.P[6] = p12 + dt * p32;
    tr.P[7] = p13 + dt * p33;

    tr.P[8] = p20 + dt * p22;
    tr.P[9] = p21 + dt * p23;
    tr.P[10] = p22;
    tr.P[11] = p23;

    tr.P[12] = p30 + dt * p32;
    tr.P[13] = p31 + dt * p33;
    tr.P[14] = p32;
    tr.P[15] = p33;

    const double q = kProcessAccelMps2 * kProcessAccelMps2;
    const double dt2 = dt * dt;
    const double dt3 = dt2 * dt;
    const double dt4 = dt2 * dt2;
    tr.P[0] += 0.25 * dt4 * q;
    tr.P[5] += 0.25 * dt4 * q;
    tr.P[2] += 0.5 * dt3 * q;
    tr.P[8] += 0.5 * dt3 * q;
    tr.P[7] += 0.5 * dt3 * q;
    tr.P[13] += 0.5 * dt3 * q;
    tr.P[10] += dt2 * q;
    tr.P[15] += dt2 * q;
}

static double mahalanobis2(const ManagedTrack& tr,
                           const DetectionPoint& det,
                           double measurement_var,
                           double& dist_out)
{
    const double rx = det.e - tr.e;
    const double ry = det.n - tr.n;
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

static void updateTrackKalman(ManagedTrack& tr,
                              const DetectionPoint& det,
                              double measurement_var,
                              double dt)
{
    const double z0 = det.e;
    const double z1 = det.n;
    const double y0 = z0 - tr.e;
    const double y1 = z1 - tr.n;

    const double s00 = tr.P[0] + measurement_var;
    const double s01 = tr.P[1];
    const double s10 = tr.P[4];
    const double s11 = tr.P[5] + measurement_var;
    const double det_s = s00 * s11 - s01 * s10;
    if (std::fabs(det_s) < 1e-9) {
        tr.e = det.e;
        tr.n = det.n;
        return;
    }

    const double inv00 = s11 / det_s;
    const double inv01 = -s01 / det_s;
    const double inv10 = -s10 / det_s;
    const double inv11 = s00 / det_s;

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

    dt = std::max(1e-3, dt);
    const double measured_ve = (det.e - (tr.e - tr.ve * dt)) / dt;
    const double measured_vn = (det.n - (tr.n - tr.vn * dt)) / dt;
    if (std::isfinite(measured_ve) && std::isfinite(measured_vn)) {
        tr.ve = 0.7 * tr.ve + 0.3 * measured_ve;
        tr.vn = 0.7 * tr.vn + 0.3 * measured_vn;
    }
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
            initCovariance(tr, gate_m);
        }
        predictTrack(tr, frameDt(tr, update_utc, result_id));
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
                const double d2 = mahalanobis2(tr, detections[di], measurement_var, dist);
                if (gate_m > 0.0 && dist > gate_m) {
                    cost[row][di] = kInvalidCost;
                    continue;
                }

                const double dt = frameDt(tr, detections[di].det.utcMid, result_id);
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
                cost[row][di] = d2
                    + 0.5 * norm_dist
                    + 0.25 * directionCost(tr.direction, detections[di].det.direction)
                    + 0.25 * rangeCost(tr.range, detections[di].det.range);
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
            const double dt = frameDt(tr, pt.det.utcMid, result_id);
            const double prev_e = tr.e;
            const double prev_n = tr.n;
            updateTrackKalman(tr, pt, measurement_var, dt);
            const double dist = std::sqrt((pt.e - prev_e) * (pt.e - prev_e)
                                        + (pt.n - prev_n) * (pt.n - prev_n));
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
        initCovariance(tr, gate_m);
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

void TrackManager::reset()
{
    tracks_.clear();
    next_id_ = 1;
    last_update_result_id_ = 0;
}

uint16_t TrackManager::allocateId()
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

void TrackManager::pruneDeleted()
{
    tracks_.erase(std::remove_if(tracks_.begin(), tracks_.end(),
                                 [](const ManagedTrack& tr) {
                                     return tr.state == TrackState::Deleted;
                                 }),
                  tracks_.end());
}
