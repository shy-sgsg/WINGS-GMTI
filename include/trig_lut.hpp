#pragma once

#include <cstddef>

namespace gmti {
namespace trig_lut {

constexpr int kTrigLutSize = 16384;

enum class Mode {
    Lut = 0,
    Math = 1,
    Compare = 2,
};

bool initialize(bool uploadDeviceTable = true);
bool isInitialized();
bool setMode(Mode mode, bool uploadDeviceMode = true);
bool setModeFromString(const char* modeName, bool uploadDeviceMode = true);
Mode mode();
const char* modeName();
void configureFromEnv();
void benchmark(std::size_t iterations = 1u << 24);
void printCompareSummary();

double sin(double x);
double cos(double x);
double tan(double x);
double asin(double x);
double acos(double x);
double atan(double x);
double atan2(double y, double x);

float sinf(float x);
float cosf(float x);
float tanf(float x);
float asinf(float x);
float acosf(float x);
float atanf(float x);
float atan2f(float y, float x);

} // namespace trig_lut
} // namespace gmti
