#pragma once

#include "trackModule.hpp"
#include <array>
#include <cstdint>
#include <deque>
#include <string>
#include <vector>

enum class TrackState {
    Tentative,
    Confirmed,
    Coasted,
    Deleted
};

struct TrackPoint {
    double e = 0.0;
    double n = 0.0;
    double utc = 0.0;
    int result_id = -1;

    TrackPoint() {}
    TrackPoint(double e_in, double n_in, double utc_in, int result_id_in)
        : e(e_in), n(n_in), utc(utc_in), result_id(result_id_in) {}
};

struct ManagedTrack {
    uint16_t id = 0;
    TrackState state = TrackState::Tentative;
    double e = 0.0;
    double n = 0.0;
    double ve = 0.0;
    double vn = 0.0;
    std::array<double, 16> P{};
    double utc = 0.0;
    double speed = 0.0;
    double direction = 0.0;
    double range = 0.0;
    int last_result_id = -1;
    int age = 0;
    int hit_count = 0;
    int miss_count = 0;
    int consecutive_hits = 0;
    bool matched_this_frame = false;
    int matched_det_index = -1;
    bool is_output_this_frame = false;
    std::deque<int> hit_history;
    std::deque<TrackPoint> point_history;
    std::deque<double> speed_history;
    std::deque<double> heading_history;
    std::deque<double> residual_history;
    std::deque<int> recent_result_ids;
};

class TrackManager {
public:
    std::vector<GMTIDetection> update(const Config& cfg,
                                      const std::vector<GMTIDetection>& current_targets);
    std::vector<GMTIDetection> updateRawDetections(
        const Config& cfg,
        const std::vector<GMTIDetection>& current_detections,
        int result_id,
        double frame_utc);
    void reset();

private:
    uint16_t allocateId();
    uint16_t allocateNextId();
    bool isActiveId(uint16_t id) const;
    ManagedTrack* findTrackById(uint16_t id);
    int currentResultId(const Config& cfg) const;
    int maxMiss(const Config& cfg) const;
    void dumpDebugSnapshot(const Config& cfg,
                           int result_id,
                           double frame_utc,
                           const std::vector<GMTIDetection>& detections,
                           const std::vector<int>& det_to_track_id,
                           int num_outputs,
                           int num_new_tracks,
                           int num_matched_tracks,
                           int num_unmatched_detections) const;
    void pruneDeleted();

    uint16_t next_id_ = 1;
    int last_update_result_id_ = 0;
    std::vector<ManagedTrack> tracks_;
};
