#pragma once

#include "trackModule.hpp"
#include <array>
#include <cstdint>
#include <deque>
#include <vector>

enum class TrackState {
    Confirmed,
    Coasted,
    Deleted
};

struct ManagedTrack {
    uint16_t id = 0;
    TrackState state = TrackState::Confirmed;
    double e = 0.0;
    double n = 0.0;
    double ve = 0.0;
    double vn = 0.0;
    std::array<double, 16> P{};
    double utc = 0.0;
    double speed = 0.0;
    double direction = 0.0;
    double range = 0.0;
    int last_result_id = 0;
    int hit_count = 0;
    int miss_count = 0;
    std::deque<int> recent_result_ids;
};

class TrackManager {
public:
    std::vector<GMTIDetection> update(const Config& cfg,
                                      const std::vector<GMTIDetection>& current_targets);
    void reset();

private:
    uint16_t allocateId();
    bool isActiveId(uint16_t id) const;
    int currentResultId(const Config& cfg) const;
    int maxMiss(const Config& cfg) const;
    void pruneDeleted();

    uint16_t next_id_ = 1;
    int last_update_result_id_ = 0;
    std::vector<ManagedTrack> tracks_;
};
