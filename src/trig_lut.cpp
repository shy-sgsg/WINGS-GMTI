#include "trig_lut.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTwoPi = 2.0 * kPi;
constexpr int kSize = gmti::trig_lut::kTrigLutSize;
constexpr int kMask = kSize - 1;

} // namespace

extern "C" bool gmtiUploadTrigLutToDevice(const float* sinTable,
                                          const float* cosTable,
                                          const float* asinTable,
                                          const float* atanTable,
                                          int tableSize);
extern "C" bool gmtiSetDeviceTrigMode(int mode);

namespace {

struct CompareStats {
    const char* name;
    std::size_t count;
    double sum_abs_err;
    double max_abs_err;
    double math_ns;
    double lut_ns;
};

std::array<float, kSize> g_sin_table;
std::array<float, kSize> g_cos_table;
std::array<float, kSize + 1> g_asin_table;
std::array<float, kSize + 1> g_atan_table;
bool g_initialized = false;
gmti::trig_lut::Mode g_mode = gmti::trig_lut::Mode::Lut;
std::mutex g_stats_mutex;
CompareStats g_stats[] = {
    {"sin", 0, 0.0, 0.0, 0.0, 0.0},
    {"cos", 0, 0.0, 0.0, 0.0, 0.0},
    {"tan", 0, 0.0, 0.0, 0.0, 0.0},
    {"asin", 0, 0.0, 0.0, 0.0, 0.0},
    {"acos", 0, 0.0, 0.0, 0.0, 0.0},
    {"atan", 0, 0.0, 0.0, 0.0, 0.0},
    {"atan2", 0, 0.0, 0.0, 0.0, 0.0},
};
bool g_summary_registered = false;

enum StatIndex {
    kStatSin = 0,
    kStatCos,
    kStatTan,
    kStatAsin,
    kStatAcos,
    kStatAtan,
    kStatAtan2,
};

inline int angleIndex(double x)
{
    const double scaled = x * (static_cast<double>(kSize) / kTwoPi);
    const long long idx = static_cast<long long>(scaled + (scaled >= 0.0 ? 0.5 : -0.5));
    return static_cast<int>(idx) & kMask;
}

inline int unitIndex(double x)
{
    x = std::max(-1.0, std::min(1.0, x));
    return static_cast<int>((x + 1.0) * (0.5 * static_cast<double>(kSize)) + 0.5);
}

inline int zeroOneIndex(double x)
{
    x = std::max(0.0, std::min(1.0, x));
    return static_cast<int>(x * static_cast<double>(kSize) + 0.5);
}

inline double lutSin(double x)
{
    return static_cast<double>(g_sin_table[angleIndex(x)]);
}

inline double lutCos(double x)
{
    return static_cast<double>(g_cos_table[angleIndex(x)]);
}

inline double lutTan(double x)
{
    const int idx = angleIndex(x);
    return static_cast<double>(g_sin_table[idx]) / static_cast<double>(g_cos_table[idx]);
}

inline double lutAsin(double x)
{
    return static_cast<double>(g_asin_table[unitIndex(x)]);
}

inline double lutAcos(double x)
{
    return 0.5 * kPi - lutAsin(x);
}

inline double lutAtan(double x)
{
    if (x > 1.0) {
        return 0.5 * kPi - static_cast<double>(g_atan_table[zeroOneIndex(1.0 / x)]);
    }
    if (x < -1.0) {
        return -0.5 * kPi + static_cast<double>(g_atan_table[zeroOneIndex(-1.0 / x)]);
    }
    const double ax = std::fabs(x);
    const double a = static_cast<double>(g_atan_table[zeroOneIndex(ax)]);
    return x < 0.0 ? -a : a;
}

inline double lutAtan2(double y, double x)
{
    if (x > 0.0) return lutAtan(y / x);
    if (x < 0.0) return (y >= 0.0 ? kPi : -kPi) + lutAtan(y / x);
    if (y > 0.0) return 0.5 * kPi;
    if (y < 0.0) return -0.5 * kPi;
    return 0.0;
}

inline void recordCompare(StatIndex idx, double mathValue, double lutValue,
                          double mathNs, double lutNs)
{
    const double err = std::fabs(mathValue - lutValue);
    std::lock_guard<std::mutex> lock(g_stats_mutex);
    CompareStats& s = g_stats[idx];
    ++s.count;
    s.sum_abs_err += err;
    s.math_ns += mathNs;
    s.lut_ns += lutNs;
    if (err > s.max_abs_err) {
        s.max_abs_err = err;
    }
}

template <typename MathFn, typename LutFn>
inline double selectValue(StatIndex idx, MathFn mathFn, LutFn lutFn)
{
    if (g_mode == gmti::trig_lut::Mode::Math) return mathFn();
    if (g_mode == gmti::trig_lut::Mode::Lut) return lutFn();

    auto mathStart = std::chrono::high_resolution_clock::now();
    const double mathValue = mathFn();
    auto mathEnd = std::chrono::high_resolution_clock::now();

    auto lutStart = std::chrono::high_resolution_clock::now();
    const double lutValue = lutFn();
    auto lutEnd = std::chrono::high_resolution_clock::now();

    if (g_mode == gmti::trig_lut::Mode::Compare) {
        const double mathNs = std::chrono::duration<double, std::nano>(mathEnd - mathStart).count();
        const double lutNs = std::chrono::duration<double, std::nano>(lutEnd - lutStart).count();
        recordCompare(idx, mathValue, lutValue, mathNs, lutNs);
        return mathValue;
    }
    return mathValue;
}

} // namespace

namespace gmti {
namespace trig_lut {

bool initialize(bool uploadDeviceTable)
{
    configureFromEnv();
    if (!g_initialized) {
        for (int i = 0; i < kSize; ++i) {
            const double angle = kTwoPi * static_cast<double>(i) / static_cast<double>(kSize);
            g_sin_table[i] = static_cast<float>(std::sin(angle));
            g_cos_table[i] = static_cast<float>(std::cos(angle));
        }
        for (int i = 0; i <= kSize; ++i) {
            const double unit = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(kSize);
            const double zeroOne = static_cast<double>(i) / static_cast<double>(kSize);
            g_asin_table[i] = static_cast<float>(std::asin(unit));
            g_atan_table[i] = static_cast<float>(std::atan(zeroOne));
        }
        g_initialized = true;
    }

    if (!g_summary_registered) {
        std::atexit(printCompareSummary);
        g_summary_registered = true;
    }

    if (uploadDeviceTable && g_mode != Mode::Math) {
        if (!gmtiUploadTrigLutToDevice(g_sin_table.data(), g_cos_table.data(),
                                       g_asin_table.data(), g_atan_table.data(), kSize)) {
            return false;
        }
    }
    return !uploadDeviceTable || gmtiSetDeviceTrigMode(static_cast<int>(g_mode));
}

bool isInitialized()
{
    return g_initialized;
}

bool setMode(Mode mode, bool uploadDeviceMode)
{
    g_mode = mode;
    std::cout << "[TRIG-LUT] mode=" << modeName() << std::endl;
    if (uploadDeviceMode) {
        return gmtiSetDeviceTrigMode(static_cast<int>(g_mode));
    }
    return true;
}

bool setModeFromString(const char* modeNameText, bool uploadDeviceMode)
{
    if (!modeNameText || modeNameText[0] == '\0') {
        return true;
    }
    if (std::strcmp(modeNameText, "lut") == 0 || std::strcmp(modeNameText, "lookup") == 0) {
        return setMode(Mode::Lut, uploadDeviceMode);
    }
    if (std::strcmp(modeNameText, "math") == 0 || std::strcmp(modeNameText, "compute") == 0) {
        return setMode(Mode::Math, uploadDeviceMode);
    }
    if (std::strcmp(modeNameText, "compare") == 0) {
        return setMode(Mode::Compare, uploadDeviceMode);
    }
    std::cerr << "[TRIG-LUT][ERR] unknown mode: " << modeNameText
              << " (expected lut|math|compare)" << std::endl;
    return false;
}

Mode mode()
{
    return g_mode;
}

const char* modeName()
{
    switch (g_mode) {
    case Mode::Lut: return "lut";
    case Mode::Math: return "math";
    case Mode::Compare: return "compare";
    }
    return "unknown";
}

void configureFromEnv()
{
    static bool env_loaded = false;
    if (env_loaded) return;
    env_loaded = true;
    const char* env = std::getenv("GMTI_TRIG_MODE");
    if (env && env[0] != '\0') {
        setModeFromString(env, false);
    }
}

double sin(double x)
{
    if (!g_initialized) initialize(false);
    return selectValue(kStatSin, [&]() { return std::sin(x); }, [&]() { return lutSin(x); });
}

double cos(double x)
{
    if (!g_initialized) initialize(false);
    return selectValue(kStatCos, [&]() { return std::cos(x); }, [&]() { return lutCos(x); });
}

double tan(double x)
{
    if (!g_initialized) initialize(false);
    return selectValue(kStatTan, [&]() { return std::tan(x); }, [&]() { return lutTan(x); });
}

double asin(double x)
{
    if (!g_initialized) initialize(false);
    return selectValue(kStatAsin, [&]() { return std::asin(x); }, [&]() { return lutAsin(x); });
}

double acos(double x)
{
    if (!g_initialized) initialize(false);
    return selectValue(kStatAcos, [&]() { return std::acos(x); }, [&]() { return lutAcos(x); });
}

double atan(double x)
{
    if (!g_initialized) initialize(false);
    return selectValue(kStatAtan, [&]() { return std::atan(x); }, [&]() { return lutAtan(x); });
}

double atan2(double y, double x)
{
    if (!g_initialized) initialize(false);
    return selectValue(kStatAtan2, [&]() { return std::atan2(y, x); }, [&]() { return lutAtan2(y, x); });
}

float sinf(float x) { return static_cast<float>(sin(static_cast<double>(x))); }
float cosf(float x) { return static_cast<float>(cos(static_cast<double>(x))); }
float tanf(float x) { return static_cast<float>(tan(static_cast<double>(x))); }
float asinf(float x) { return static_cast<float>(asin(static_cast<double>(x))); }
float acosf(float x) { return static_cast<float>(acos(static_cast<double>(x))); }
float atanf(float x) { return static_cast<float>(atan(static_cast<double>(x))); }
float atan2f(float y, float x) { return static_cast<float>(atan2(static_cast<double>(y), static_cast<double>(x))); }

void benchmark(std::size_t iterations)
{
    initialize(false);
    volatile double sink = 0.0;
    const double step = 0.000317;
    double maxErr = 0.0;
    double sumErr = 0.0;

    auto t0 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < iterations; ++i) {
        const double x = static_cast<double>(i) * step;
        sink += std::sin(x) + std::cos(x);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < iterations; ++i) {
        const double x = static_cast<double>(i) * step;
        const double lutValue = lutSin(x) + lutCos(x);
        sink += lutValue;
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    for (std::size_t i = 0; i < iterations; ++i) {
        const double x = static_cast<double>(i) * step;
        const double lutValue = lutSin(x) + lutCos(x);
        const double mathValue = std::sin(x) + std::cos(x);
        const double err = std::fabs(mathValue - lutValue);
        sumErr += err;
        if (err > maxErr) maxErr = err;
    }

    const double mathMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double lutMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << std::setprecision(8)
              << "[TRIG-LUT][CPU] mode=" << modeName()
              << " iterations=" << iterations
              << " std_sin_cos_ms=" << mathMs
              << " lut_sin_cos_ms=" << lutMs
              << " speedup=" << (lutMs > 0.0 ? mathMs / lutMs : 0.0)
              << " mean_abs_err=" << (iterations > 0 ? sumErr / static_cast<double>(iterations) : 0.0)
              << " max_abs_err=" << maxErr
              << " sink=" << sink << std::endl;
}

void printCompareSummary()
{
    if (g_mode != Mode::Compare) return;
    std::lock_guard<std::mutex> lock(g_stats_mutex);
    std::cout << "[TRIG-LUT][COMPARE] precision/time summary";
    bool any = false;
    for (const auto& s : g_stats) {
        if (s.count == 0) continue;
        any = true;
        const double mathMs = s.math_ns / 1.0e6;
        const double lutMs = s.lut_ns / 1.0e6;
        std::cout << "\n  " << s.name
                  << ": count=" << s.count
                  << " mean_abs_err=" << (s.sum_abs_err / static_cast<double>(s.count))
                  << " max_abs_err=" << s.max_abs_err
                  << " math_ms=" << mathMs
                  << " lut_ms=" << lutMs
                  << " speedup=" << (lutMs > 0.0 ? mathMs / lutMs : 0.0);
    }
    if (!any) {
        std::cout << ": no trig calls recorded";
    }
    std::cout << std::endl;
}

} // namespace trig_lut
} // namespace gmti
