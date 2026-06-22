#pragma once

#include <cuda_runtime.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gmti {
namespace trig_lut_device {

constexpr int kTrigLutSize = 16384;
constexpr int kTrigLutMask = kTrigLutSize - 1;
constexpr int kModeLut = 0;
constexpr int kModeMath = 1;
constexpr int kModeCompare = 2;

struct TrigLutConfig {
    const float* sin_table;
    const float* cos_table;
    const float* asin_table;
    const float* atan_table;
    int mode;
};

__device__ __forceinline__ int angleIndex(float x)
{
    constexpr float kScale = static_cast<float>(kTrigLutSize) / (2.0f * static_cast<float>(M_PI));
    int idx = static_cast<int>(x * kScale + (x >= 0.0f ? 0.5f : -0.5f));
    return idx & kTrigLutMask;
}

__device__ __forceinline__ int unitIndex(float x)
{
    x = fminf(1.0f, fmaxf(-1.0f, x));
    return static_cast<int>((x + 1.0f) * (0.5f * static_cast<float>(kTrigLutSize)) + 0.5f);
}

__device__ __forceinline__ int zeroOneIndex(float x)
{
    x = fminf(1.0f, fmaxf(0.0f, x));
    return static_cast<int>(x * static_cast<float>(kTrigLutSize) + 0.5f);
}

__device__ __forceinline__ float lut_sinf(float x, const TrigLutConfig& cfg)
{
    return cfg.sin_table[angleIndex(x)];
}

__device__ __forceinline__ float lut_cosf(float x, const TrigLutConfig& cfg)
{
    return cfg.cos_table[angleIndex(x)];
}

__device__ __forceinline__ float sinf(float x, const TrigLutConfig& cfg)
{
    return cfg.mode == kModeLut ? lut_sinf(x, cfg) : ::sinf(x);
}

__device__ __forceinline__ float cosf(float x, const TrigLutConfig& cfg)
{
    return cfg.mode == kModeLut ? lut_cosf(x, cfg) : ::cosf(x);
}

__device__ __forceinline__ void sincosf(float x, float* s, float* c,
                                        const TrigLutConfig& cfg)
{
    if (cfg.mode != kModeLut) {
        ::sincosf(x, s, c);
        return;
    }
    const int idx = angleIndex(x);
    *s = cfg.sin_table[idx];
    *c = cfg.cos_table[idx];
}

__device__ __forceinline__ float tanf(float x, const TrigLutConfig& cfg)
{
    if (cfg.mode != kModeLut) return ::tanf(x);
    const int idx = angleIndex(x);
    return cfg.sin_table[idx] / cfg.cos_table[idx];
}

__device__ __forceinline__ float lut_asinf(float x, const TrigLutConfig& cfg)
{
    return cfg.asin_table[unitIndex(x)];
}

__device__ __forceinline__ float asinf(float x, const TrigLutConfig& cfg)
{
    return cfg.mode == kModeLut ? lut_asinf(x, cfg) : ::asinf(x);
}

__device__ __forceinline__ float acosf(float x, const TrigLutConfig& cfg)
{
    return cfg.mode == kModeLut
        ? 1.57079632679489661923f - lut_asinf(x, cfg)
        : ::acosf(x);
}

__device__ __forceinline__ float lut_atanf(float x, const TrigLutConfig& cfg)
{
    if (x > 1.0f) {
        return 1.57079632679489661923f - cfg.atan_table[zeroOneIndex(1.0f / x)];
    }
    if (x < -1.0f) {
        return -1.57079632679489661923f + cfg.atan_table[zeroOneIndex(-1.0f / x)];
    }
    const float ax = fabsf(x);
    const float a = cfg.atan_table[zeroOneIndex(ax)];
    return x < 0.0f ? -a : a;
}

__device__ __forceinline__ float atanf(float x, const TrigLutConfig& cfg)
{
    return cfg.mode == kModeLut ? lut_atanf(x, cfg) : ::atanf(x);
}

__device__ __forceinline__ float lut_atan2f(float y, float x,
                                           const TrigLutConfig& cfg)
{
    if (x > 0.0f) return lut_atanf(y / x, cfg);
    if (x < 0.0f) {
        return (y >= 0.0f ? static_cast<float>(M_PI) : -static_cast<float>(M_PI))
            + lut_atanf(y / x, cfg);
    }
    if (y > 0.0f) return 1.57079632679489661923f;
    if (y < 0.0f) return -1.57079632679489661923f;
    return 0.0f;
}

__device__ __forceinline__ float atan2f(float y, float x,
                                       const TrigLutConfig& cfg)
{
    return cfg.mode == kModeLut ? lut_atan2f(y, x, cfg) : ::atan2f(y, x);
}

} // namespace trig_lut_device
} // namespace gmti
