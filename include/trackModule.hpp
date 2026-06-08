#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "config_structs.hpp"

// 常量定义，需与MATLAB一致
const double LSB = 8.38191e-8;
const int MAX_TGT = 1000;

struct GMTIDetection {
    uint16_t id;
    double lon;
    double lat;
    double speed;
    double direction;
    double range;
    double utcMid;
};

struct GMTIResult {
    int count;
    std::vector<GMTIDetection> targets;
    double utcGlobal;
};

class TrackManager;

std::vector<GMTIDetection> trackModule(const Config& cfg);
std::vector<GMTIDetection> trackModule(const Config& cfg, TrackManager* manager);
