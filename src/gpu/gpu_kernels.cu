#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/functional.h>
#include <thrust/count.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>
#include <cstdint>
#include <iostream>
#include <chrono>
#include "GMTIProcessor.hpp"
#include "trig_lut.hpp"
#include "trig_lut_device.cuh"

extern "C" gmti::trig_lut_device::TrigLutConfig gmtiGetDeviceTrigLutConfig();

// Use cuFFT / CUDA types for device kernels
using cudacd = cuFloatComplex;

#define CUDA_CHECK(call) \
do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return false; \
    } \
} while (0)

#define CUFFT_CHECK(call) \
do { \
    cufftResult err = call; \
    if (err != CUFFT_SUCCESS) { \
        std::cerr << "cuFFT error: " << err << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return false; \
    } \
} while (0)

__global__ void init_labels_kernel(const float* mydata, int* labels, int total);
__global__ void propagate_labels_kernel(const int* labels_in, int* labels_out,
                                        int H, int W, int max_gap, int* changed);
__global__ void accumulate_stats_kernel(const float* mydata, const float* phase_map,
                                        const int* labels, int total,
                                        int* counts, float* sum_cos, float* sum_sin,
                                        float* max_power, int* max_idx,
                                        bool use_phase,
                                        gmti::trig_lut_device::TrigLutConfig trig_cfg);
__global__ void compute_mean_phase_kernel(const int* counts, const float* sum_cos,
                                          const float* sum_sin, float* mean_phi, int total,
                                          gmti::trig_lut_device::TrigLutConfig trig_cfg);
__global__ void accumulate_phase_var_kernel(const float* phase_map, const int* labels,
                                            const float* mean_phi, int total,
                                            float* sum_sq);
__global__ void compute_phase_std_kernel(const int* counts, const float* sum_sq,
                                         float* phase_std, int total);
__global__ void compute_amp_kernel(const cuFloatComplex* F1,
                                   const int* prow,
                                   const int* pcol,
                                   float* amp,
                                   int total_points,
                                   int H,
                                   int W);

__device__ float wrap_angle_rad_device(float x)
{
    const float PI = 3.14159265358979323846f;
    const float TWO_PI = 6.28318530717958647692f;
    while (x <= -PI) x += TWO_PI;
    while (x > PI) x -= TWO_PI;
    return x;
}

__device__ float atomic_max_float(float* addr, float value)
{
    int* addr_as_int = reinterpret_cast<int*>(addr);
    int old = *addr_as_int;
    int assumed;
    do {
        assumed = old;
        float old_val = __int_as_float(assumed);
        if (old_val >= value) break;
        old = atomicCAS(addr_as_int, assumed, __float_as_int(value));
    } while (assumed != old);
    return __int_as_float(old);
}

__device__ float atomic_add_float(float* addr, float value)
{
    return atomicAdd(addr, value);
}

struct NonNegativeAmp
{
    __host__ __device__ bool operator()(float v) const { return v >= 0.0f; }
};


/**
 * @brief 距离向相位补偿核函数
 * 每个线程负责一个距离门 (Column)，纵向处理所有脉冲
 */
__global__ void apply_phase_correction_kernel(
    cuFloatComplex* F1, 
    cuFloatComplex* F2,
    const float* phi_fit, 
    int Na, int Nr,
    gmti::trig_lut_device::TrigLutConfig trig_cfg)
{
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= Nr) return;

    float angle = phi_fit[c];
    float s, cr;
    gmti::trig_lut_device::sincosf(angle, &s, &cr, trig_cfg);
    cuFloatComplex phi_factor = make_cuFloatComplex(cr, s);

    for (int r = 0; r < Na; ++r) {
        size_t idx = (size_t)r * Nr + c;
        F2[idx] = cuCmulf(F2[idx], phi_factor);
    }
}

// 计算两通道共轭相乘的纵向累加
__global__ void compute_phase_sum_kernel(
    const cuFloatComplex* F1, 
    const cuFloatComplex* F2,
    cuFloatComplex* out_sums,
    int az_st, int az_ed,
    int Na, int Nr,
    int rg_st, int rg_ed) 
{
    // 每个线程处理一个距离向索引 c
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    
    // 越界检查：只处理用户指定的距离向范围 [rg_st, rg_ed]
    if (c < rg_st || c > rg_ed) return;

    cuFloatComplex sum = make_cuFloatComplex(0.0f, 0.0f);

    // 纵向累加指定的多普勒范围 [az_st, az_ed]
    for (int r = az_st; r <= az_ed; ++r) {
        size_t idx = (size_t)r * Nr + c;
        
        cuFloatComplex val1 = F1[idx];
        cuFloatComplex val2 = F2[idx];
        
        // 计算 F1 * conj(F2)
        // cuCmulf(a, cuConjf(b))
        cuFloatComplex res = cuCmulf(val1, cuConjf(val2));
        
        sum = cuCaddf(sum, res);
    }

    // 将结果写入输出向量（长度为 Nr 的显存空间）
    out_sums[c] = sum;
}

__global__ void compute_row_phase_kernel(
    const cuFloatComplex* F1,
    const cuFloatComplex* F2,
    float* phase_out,
    int Na,
    int Nr,
    int az_st,
    int az_ed,
    gmti::trig_lut_device::TrigLutConfig trig_cfg)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r < az_st || r > az_ed || r >= Na) return;

    const size_t off = static_cast<size_t>(r) * static_cast<size_t>(Nr);
    float num_re = 0.0f;
    float num_im = 0.0f;
    float den = 0.0f;

    for (int c = 0; c < Nr; ++c) {
        const cuFloatComplex a = F1[off + c];
        const cuFloatComplex b = F2[off + c];
        num_re += (cuCrealf(a) * cuCrealf(b) + cuCimagf(a) * cuCimagf(b));
        num_im += (cuCimagf(a) * cuCrealf(b) - cuCrealf(a) * cuCimagf(b));
        den += cuCrealf(a) * cuCrealf(a) + cuCimagf(a) * cuCimagf(a);
    }

    if (den > 0.0f) {
        num_re /= den;
        num_im /= den;
    }
    phase_out[r] = gmti::trig_lut_device::atan2f(num_im, num_re, trig_cfg);
}

__global__ void phase_map_kernel(
    const cuFloatComplex* F1,
    const cuFloatComplex* F2,
    float* out,
    size_t total,
    gmti::trig_lut_device::TrigLutConfig trig_cfg)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const cuFloatComplex a = F1[idx];
    const cuFloatComplex b = F2[idx];
    float num_re = (cuCrealf(a) * cuCrealf(b) + cuCimagf(a) * cuCimagf(b));
    float num_im = (cuCimagf(a) * cuCrealf(b) - cuCrealf(a) * cuCimagf(b));
    out[idx] = gmti::trig_lut_device::atan2f(num_im, num_re, trig_cfg);
}

__global__ void mix_detect_data_kernel(
    const cuFloatComplex* f2_in,
    const cuFloatComplex* csi_in,
    cuFloatComplex* out,
    int H,
    int W,
    int band_st,
    int band_ed)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = static_cast<size_t>(H) * static_cast<size_t>(W);
    if (idx >= total) return;

    int r = static_cast<int>(idx / W);
    if (r >= band_st && r <= band_ed) {
        out[idx] = csi_in[idx];
    } else {
        out[idx] = f2_in[idx];
    }
}

__global__ void power_map_kernel(
    const cuFloatComplex* in,
    float* power,
    size_t total)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const cuFloatComplex v = in[idx];
    const float re = cuCrealf(v);
    const float im = cuCimagf(v);
    power[idx] = re * re + im * im;
}

__global__ void row_prefix_kernel(
    const float* power,
    float* row_prefix,
    int H,
    int W)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= H) return;

    float sum = 0.0f;
    size_t base = static_cast<size_t>(r) * static_cast<size_t>(W);
    for (int c = 0; c < W; ++c) {
        sum += power[base + static_cast<size_t>(c)];
        row_prefix[base + static_cast<size_t>(c)] = sum;
    }
}

__global__ void col_prefix_to_integral_kernel(
    const float* row_prefix,
    float* integral,
    int H,
    int W)
{
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= W) return;

    float sum = 0.0f;
    for (int r = 0; r < H; ++r) {
        sum += row_prefix[static_cast<size_t>(r) * static_cast<size_t>(W) + static_cast<size_t>(c)];
        size_t idx = static_cast<size_t>(r + 1) * static_cast<size_t>(W + 1) + static_cast<size_t>(c + 1);
        integral[idx] = sum;
    }
}

__device__ float rect_sum_device(const float* integral, int W, int r0, int c0, int r1, int c1)
{
    size_t a = static_cast<size_t>(r1 + 1) * static_cast<size_t>(W + 1) + static_cast<size_t>(c1 + 1);
    size_t b0 = static_cast<size_t>(r0) * static_cast<size_t>(W + 1) + static_cast<size_t>(c1 + 1);
    size_t c0i = static_cast<size_t>(r1 + 1) * static_cast<size_t>(W + 1) + static_cast<size_t>(c0);
    size_t d = static_cast<size_t>(r0) * static_cast<size_t>(W + 1) + static_cast<size_t>(c0);
    return integral[a] - integral[b0] - integral[c0i] + integral[d];
}

__global__ void cfar_detect_kernel(
    const cuFloatComplex* detect,
    const float* integral,
    float* out,
    int H,
    int W,
    int g,
    int b,
    int total_bg,
    float alpha,
    bool use_go)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = static_cast<size_t>(H) * static_cast<size_t>(W);
    if (idx >= total) return;

    int r = static_cast<int>(idx / W);
    int c = static_cast<int>(idx % W);

    const int R = g + b;
    const int r_min = R;
    const int r_max = H - 1 - R;
    const int c_min = R;
    const int c_max = W - 1 - R;

    if (r < r_min || r > r_max || c < c_min || c > c_max) {
        out[idx] = 0.0f;
        return;
    }

    const cuFloatComplex v = detect[idx];
    const float re = cuCrealf(v);
    const float im = cuCimagf(v);
    const float CUT = re * re + im * im;

    const int r0 = r - R, r1 = r + R;
    const int c0 = c - R, c1 = c + R;

    const int rg0 = r - g, rg1 = r + g;
    const int cg0 = c - g, cg1 = c + g;

    float noise_level = 0.0f;
    if (!use_go) {
        const float sum_full = rect_sum_device(integral, W, r0, c0, r1, c1);
        const float sum_guard = rect_sum_device(integral, W, rg0, cg0, rg1, cg1);
        const float bg_sum = sum_full - sum_guard;
        noise_level = bg_sum / static_cast<float>(total_bg);
    } else {
        float mleft = 0.0f, mright = 0.0f, mtop = 0.0f, mbot = 0.0f;

        const int lc0 = c0, lc1 = cg0 - 1;
        const int rc0 = cg1 + 1, rc1 = c1;
        const int tr0 = r0, tr1 = rg0 - 1;
        const int br0 = rg1 + 1, br1 = r1;

        if (lc0 <= lc1) {
            const float s = rect_sum_device(integral, W, r0, lc0, r1, lc1);
            const int npx = (r1 - r0 + 1) * (lc1 - lc0 + 1);
            if (npx > 0) mleft = s / npx;
        }
        if (rc0 <= rc1) {
            const float s = rect_sum_device(integral, W, r0, rc0, r1, rc1);
            const int npx = (r1 - r0 + 1) * (rc1 - rc0 + 1);
            if (npx > 0) mright = s / npx;
        }
        if (tr0 <= tr1) {
            const float s = rect_sum_device(integral, W, tr0, c0, tr1, c1);
            const int npx = (tr1 - tr0 + 1) * (c1 - c0 + 1);
            if (npx > 0) mtop = s / npx;
        }
        if (br0 <= br1) {
            const float s = rect_sum_device(integral, W, br0, c0, br1, c1);
            const int npx = (br1 - br0 + 1) * (c1 - c0 + 1);
            if (npx > 0) mbot = s / npx;
        }

        float m1 = (mleft > mright) ? mleft : mright;
        float m2 = (mtop > mbot) ? mtop : mbot;
        noise_level = (m1 > m2) ? m1 : m2;
    }

    const float thr = alpha * noise_level;
    if (CUT > thr) {
        if (r < 70 || r > 190) {
            out[idx] = CUT;
        } else {
            out[idx] = 0.0f;
        }
    } else {
        out[idx] = 0.0f;
    }
}

__global__ void init_labels_kernel(const float* mydata, int* labels, int total)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;
    labels[idx] = (mydata[idx] > 0.0f) ? idx : -1;
}

__global__ void propagate_labels_kernel(const int* labels_in, int* labels_out,
                                        int H, int W, int max_gap, int* changed)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = H * W;
    if (idx >= total) return;

    int label = labels_in[idx];
    if (label < 0) {
        labels_out[idx] = -1;
        return;
    }

    int r = idx / W;
    int c = idx - r * W;
    int min_label = label;

    for (int dr = -1; dr <= 1; ++dr) {
        int rr = r + dr;
        if (rr < 0 || rr >= H) continue;
        int row_base = rr * W;
        for (int dc = -max_gap; dc <= max_gap; ++dc) {
            if (dr == 0 && dc == 0) continue;
            int cc = c + dc;
            if (cc < 0 || cc >= W) continue;
            int nidx = row_base + cc;
            int nlabel = labels_in[nidx];
            if (nlabel >= 0 && nlabel < min_label) min_label = nlabel;
        }
    }

    labels_out[idx] = min_label;
    if (min_label != label) atomicExch(changed, 1);
}

__global__ void accumulate_stats_kernel(const float* mydata, const float* phase_map,
                                        const int* labels, int total,
                                        int* counts, float* sum_cos, float* sum_sin,
                                        float* max_power, int* max_idx,
                                        bool use_phase,
                                        gmti::trig_lut_device::TrigLutConfig trig_cfg)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int label = labels[idx];
    if (label < 0) return;

    atomicAdd(&counts[label], 1);

    if (use_phase) {
        float phi = phase_map[idx];
        atomic_add_float(&sum_cos[label], gmti::trig_lut_device::cosf(phi, trig_cfg));
        atomic_add_float(&sum_sin[label], gmti::trig_lut_device::sinf(phi, trig_cfg));
    }

    float p = mydata[idx];
    float old = atomic_max_float(&max_power[label], p);
    if (p > old) {
        atomicExch(&max_idx[label], idx);
    }
}

__global__ void compute_mean_phase_kernel(const int* counts, const float* sum_cos,
                                          const float* sum_sin, float* mean_phi, int total,
                                          gmti::trig_lut_device::TrigLutConfig trig_cfg)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;
    if (counts[idx] > 0) {
        mean_phi[idx] = gmti::trig_lut_device::atan2f(
            sum_sin[idx], sum_cos[idx], trig_cfg);
    } else {
        mean_phi[idx] = 0.0f;
    }
}

__global__ void accumulate_phase_var_kernel(const float* phase_map, const int* labels,
                                            const float* mean_phi, int total,
                                            float* sum_sq)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int label = labels[idx];
    if (label < 0) return;

    float dphi = phase_map[idx] - mean_phi[label];
    dphi = wrap_angle_rad_device(dphi);
    atomic_add_float(&sum_sq[label], dphi * dphi);
}

__global__ void compute_phase_std_kernel(const int* counts, const float* sum_sq,
                                         float* phase_std, int total)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;
    if (counts[idx] > 1) {
        phase_std[idx] = sqrtf(sum_sq[idx] / static_cast<float>(counts[idx] - 1));
    } else {
        phase_std[idx] = 0.0f;
    }
}

__global__ void compute_amp_kernel(const cuFloatComplex* F1,
                                   const int* prow,
                                   const int* pcol,
                                   float* amp,
                                   int total_points,
                                   int H,
                                   int W)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_points) return;

    int r = prow[idx];
    int c = pcol[idx];
    if (r < 0 || c < 0 || r >= H || c >= W) {
        amp[idx] = -1.0f;
        return;
    }

    size_t off = static_cast<size_t>(r) * static_cast<size_t>(W) + static_cast<size_t>(c);
    const cuFloatComplex v = F1[off];
    float re = cuCrealf(v);
    float im = cuCimagf(v);
    amp[idx] = sqrtf(re * re + im * im);
}

__global__ void clutter_cancel_kernel(
    const cuFloatComplex* F1,
    const cuFloatComplex* F2,
    const float* y_faAxis,
    cuFloatComplex* out,
    int Na,
    int Nr,
    float k,
    float b,
    gmti::trig_lut_device::TrigLutConfig trig_cfg)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = static_cast<size_t>(Na) * static_cast<size_t>(Nr);
    if (idx >= total) return;

    int r = static_cast<int>(idx / Nr);
    const cuFloatComplex a = F1[idx];
    const cuFloatComplex f2 = F2[idx];

    float phi = k * y_faAxis[r] + b;
    float s, c;
    gmti::trig_lut_device::sincosf(phi, &s, &c, trig_cfg);
    const cuFloatComplex az_fai = make_cuFloatComplex(c, s);
    const cuFloatComplex b2 = cuCmulf(az_fai, f2);

    float aa = cuCabsf(a);
    float bb = cuCabsf(b2);
    float m = (aa < bb) ? aa : bb;

    cuFloatComplex a_eq = make_cuFloatComplex(0.0f, 0.0f);
    cuFloatComplex b_eq = make_cuFloatComplex(0.0f, 0.0f);
    if (aa > 0.0f) {
        float scale = m / aa;
        a_eq = make_cuFloatComplex(cuCrealf(a) * scale, cuCimagf(a) * scale);
    }
    if (bb > 0.0f) {
        float scale = m / bb;
        b_eq = make_cuFloatComplex(cuCrealf(b2) * scale, cuCimagf(b2) * scale);
    }

    out[idx] = cuCsubf(a_eq, b_eq);
}

// CUDA kernel for channel alignment
// *** wichtig: 与 CPU 版本对齐，应用 (-1)^row 预乘因子（模拟 FFTshift）***
__global__ void align_two_channels_fast_cuda(const cudacd* d1,
                                             const cudacd* d2,
                                             cudacd* a1,
                                             cudacd* a2,
                                             int skip,
                                             size_t Na,
                                             size_t Nr)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = Na * Nr;
    if (idx >= total) return;

    size_t row = idx / Nr;
    size_t col = idx % Nr;
    
    // 边界检查
    if (row >= Na || col >= Nr) return;
    
    // 计算 (-1)^row 预乘因子，这是为了对齐 CPU 的 FFTW + fftshift 逻辑
    float s = (row & 1) ? -1.0f : 1.0f;

    // 默认填零
    cudacd res1 = make_cuFloatComplex(0.0f, 0.0f);
    cudacd res2 = make_cuFloatComplex(0.0f, 0.0f);

    if (skip > 0) {
        // CPU 逻辑：d1 从行 (skip-1) 开始复制到 a1 开头，前面的行填零
        // size_t src_row0 = (size_t)(skip - 1);
        size_t src_row0 = (size_t)(skip);
        if (row < (Na - src_row0)) {
            // 有效行：从 d1[src_row0 + row] 复制
            size_t src_idx = (src_row0 + row) * Nr + col;
            res1 = d1[src_idx];
        } else {
            // 无效行：填零
            res1 = make_cuFloatComplex(0.0f, 0.0f);
        }
        // d2 保持不变
        res2 = d2[idx];
    } 
    else if (skip < 0) {
        // CPU 逻辑：d2 从行 (-skip-1) 开始复制到 a2 开头，前面的行填零
        // size_t src_row0 = (size_t)(-skip - 1);
        size_t src_row0 = (size_t)(-skip);
        if (row < (Na - src_row0)) {
            // 有效行：从 d2[src_row0 + row] 复制
            size_t src_idx = (src_row0 + row) * Nr + col;
            res2 = d2[src_idx];
        } else {
            // 无效行：填零
            res2 = make_cuFloatComplex(0.0f, 0.0f);
        }
        // d1 保持不变
        res1 = d1[idx];
    }
    else {
        // skip == 0，直接复制
        res1 = d1[idx];
        res2 = d2[idx];
    }

    // 统一应用预乘因子 s 并写入
    a1[idx] = make_cuFloatComplex(cuCrealf(res1) * s, cuCimagf(res1) * s);
    a2[idx] = make_cuFloatComplex(cuCrealf(res2) * s, cuCimagf(res2) * s);
}

// CUDA kernel for DBS center (circular shift rows)
__global__ void dbs_center_by_fa_fast_cuda(const cudacd* inA,
                                           const cudacd* inB,
                                           cudacd* outA,
                                           cudacd* outB,
                                           size_t Na,
                                           size_t Nr,
                                           size_t k)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = Na * Nr;
    if (idx >= total) return;

    size_t row = idx / Nr;
    size_t col = idx % Nr;
    
    // 边界检查
    if (row >= Na || col >= Nr) return;
    
    // 鲁棒的取模运算：确保结果永远为正
    size_t shift = (Na - k) % Na;
    size_t new_row = (row + shift) % Na;
    size_t new_idx = new_row * Nr + col;

    outA[new_idx] = inA[idx];
    outB[new_idx] = inB[idx];
}

__global__ void az_decimate_two_channels_kernel(const cudacd* in1,
                                                const cudacd* in2,
                                                cudacd* out1,
                                                cudacd* out2,
                                                int W_orig,
                                                int M,
                                                int dec,
                                                int W_new)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = static_cast<size_t>(W_new) * static_cast<size_t>(M);
    if (idx >= total) return;

    int k_new = static_cast<int>(idx / M);
    int m = static_cast<int>(idx % M);
    size_t in_base = static_cast<size_t>(k_new) * static_cast<size_t>(dec) * static_cast<size_t>(M);
    cudacd acc1 = make_cuFloatComplex(0.0f, 0.0f);
    cudacd acc2 = make_cuFloatComplex(0.0f, 0.0f);
    for (int d = 0; d < dec; ++d) {
        size_t in_idx = in_base + static_cast<size_t>(d) * static_cast<size_t>(M) + static_cast<size_t>(m);
        if (in_idx < static_cast<size_t>(W_orig) * static_cast<size_t>(M)) {
            acc1 = cuCaddf(acc1, in1[in_idx]);
            acc2 = cuCaddf(acc2, in2[in_idx]);
        }
    }
    out1[idx] = acc1;
    out2[idx] = acc2;
}

__global__ void scale_complex_kernel(cudacd* data, size_t total, float coef)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;
    cudacd v = data[idx];
    data[idx] = make_cuFloatComplex(v.x * coef, v.y * coef);
}

__global__ void adjacent_corr_reduce_kernel(const cudacd* data,
                                            cudacd* block_sums,
                                            int k,
                                            int Na,
                                            int Nr)
{
    extern __shared__ cudacd shared[];
    int tid = threadIdx.x;
    size_t total = static_cast<size_t>(Na - k) * static_cast<size_t>(Nr);
    cudacd acc = make_cuFloatComplex(0.0f, 0.0f);

    for (size_t idx = blockIdx.x * blockDim.x + tid;
         idx < total;
         idx += static_cast<size_t>(blockDim.x) * static_cast<size_t>(gridDim.x)) {
        int row = static_cast<int>(idx / Nr) + k;
        int col = static_cast<int>(idx % Nr);
        size_t cur = static_cast<size_t>(row) * static_cast<size_t>(Nr) + static_cast<size_t>(col);
        size_t prev = static_cast<size_t>(row - k) * static_cast<size_t>(Nr) + static_cast<size_t>(col);
        acc = cuCaddf(acc, cuCmulf(data[cur], cuConjf(data[prev])));
    }

    shared[tid] = acc;
    __syncthreads();

    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            shared[tid] = cuCaddf(shared[tid], shared[tid + stride]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        block_sums[blockIdx.x] = shared[0];
    }
}

// Production-grade CUDA implementation with pre-allocated workspace
bool GMTIProcessor::alignFFTAndDBS_CUDA(const std::vector<std::complex<float>> &data1,
                                        const std::vector<std::complex<float>> &data2,
                                        int skip,
                                        float fa2,
                                        const Config &cfg,
                                        std::vector<std::complex<float>> &out1,
                                        std::vector<std::complex<float>> &out2)
{
    // --- 准备 CUDA Event 用于测量内核执行耗时 ---
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    auto wall_start = std::chrono::high_resolution_clock::now();

    const size_t Na = static_cast<size_t>(effectivePulseNum(cfg));
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    const size_t total = Na * Nr;
    if (Na == 0 || Nr == 0) return false;
    if (data1.size() != total || data2.size() != total) return false;

    // ===== 1) 内存管理：按需分配或复用 =====
    // 需要 6 个缓冲：d_d1, d_d2, d_a1, d_a2, d_t1, d_t2
    size_t needed_bytes = 6 * total * sizeof(cudacd);
    
    if (d_workspace == nullptr || d_workspace_bytes < needed_bytes) {
        // 释放旧内存
        if (d_workspace != nullptr) {
            CUDA_CHECK(cudaFree(d_workspace));
        }
        // 新分配
        CUDA_CHECK(cudaMalloc(&d_workspace, needed_bytes));
        d_workspace_bytes = needed_bytes;
        DBG("CUDA workspace allocated: " << (needed_bytes / (1024*1024)) << " MB");
    }

    // 指针分散：将 workspace 分成 6 个部分
    cudacd* d_d1 = reinterpret_cast<cudacd*>(d_workspace) + 0 * total;
    cudacd* d_d2 = reinterpret_cast<cudacd*>(d_workspace) + 1 * total;
    cudacd* d_a1 = reinterpret_cast<cudacd*>(d_workspace) + 2 * total;
    cudacd* d_a2 = reinterpret_cast<cudacd*>(d_workspace) + 3 * total;
    cudacd* d_t1 = reinterpret_cast<cudacd*>(d_workspace) + 4 * total;
    cudacd* d_t2 = reinterpret_cast<cudacd*>(d_workspace) + 5 * total;

    // Host to Device
    CUDA_CHECK(cudaMemcpy(d_d1, reinterpret_cast<const cudacd*>(data1.data()), total * sizeof(cudacd), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_d2, reinterpret_cast<const cudacd*>(data2.data()), total * sizeof(cudacd), cudaMemcpyHostToDevice));

    // 初始化中间缓冲区为零，避免残留数据导致的NaN
    CUDA_CHECK(cudaMemset(d_t1, 0, total * sizeof(cudacd)));
    CUDA_CHECK(cudaMemset(d_t2, 0, total * sizeof(cudacd)));

    cudaEventRecord(start);
    // ===== 2) 通道对齐 =====
    int threads = 256;
    int blocks = static_cast<int>((total + threads - 1) / threads);
    align_two_channels_fast_cuda<<<blocks, threads>>>(d_d1, d_d2, d_a1, d_a2, skip, Na, Nr);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // ===== 3) 列向 FFT （使用缓存的计划）=====
    // 检查是否需要重建计划
    if (cached_Na_ != (int)Na || cached_Nr_ != (int)Nr) {
        // 销毁旧计划
        if (cufft_plan_ != -1) {
            CUFFT_CHECK(cufftDestroy((cufftHandle)cufft_plan_));
        }

        // 创建新计划
        cufftHandle plan;
        int n[] = {(int)Na}; 
        int howmany = (int)Nr;

        // 输入布局：
        // 相邻元素间距 istride = Nr (跨过一行取同列下一个点)
        // 相邻批次间距 idist = 1   (第一列算完，下一批从第一行第二个点开始)
        int istride = (int)Nr;
        int idist = 1;

        // 输出布局：
        // 如果你希望输出保持 HxW 布局（不转置），则 ostride/odist 必须与输入一致
        int ostride = (int)Nr;
        int odist = 1;

        CUFFT_CHECK(cufftPlanMany(&plan, 1, n,
                      n, istride, idist,
                      n, ostride, odist,
                      CUFFT_C2C, howmany));

        // 缓存计划（cufftHandle 本质是 int）
        cufft_plan_ = (int)plan;
        cached_Na_ = static_cast<int>(Na);
        cached_Nr_ = static_cast<int>(Nr);
        DBG("cuFFT plan created and cached for Na=" << Na << ", Nr=" << Nr);
    }

    // 使用缓存的计划执行 FFT
    cufftHandle plan = (cufftHandle)cufft_plan_;
    CUFFT_CHECK(cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_a1), 
                             reinterpret_cast<cufftComplex*>(d_a1), CUFFT_FORWARD));
    CUFFT_CHECK(cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_a2), 
                             reinterpret_cast<cufftComplex*>(d_a2), CUFFT_FORWARD));

    // ===== 4) DBS 中心化 =====
    int width = static_cast<int>(Na);
    if (!(cfg.PRF > 0.0)) return false;
    int center_num = static_cast<int>(std::floor((fa2 + 0.5f * static_cast<float>(cfg.PRF)) /
                                                 static_cast<float>(cfg.PRF) * width)) + 1;
    int cstart = center_num - width / 2;
    int k_1b = cstart;
    while (k_1b < 1) k_1b += width;
    while (k_1b > width) k_1b -= width;
    size_t k = static_cast<size_t>(k_1b - 1);

    // 使用明确的最终指针，避免局部指针交换导致的内存访问错误
    cudacd* res_ptr1 = d_a1;
    cudacd* res_ptr2 = d_a2;

    if (k != 0) {
        dbs_center_by_fa_fast_cuda<<<blocks, threads>>>(d_a1, d_a2, d_t1, d_t2, Na, Nr, k);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        // 明确指向 DBS 后的结果
        res_ptr1 = d_t1;
        res_ptr2 = d_t2;
    }
    cudaEventRecord(stop);

    // ===== 5) Device to Host =====
    out1.resize(total);
    out2.resize(total);
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<cudacd*>(out1.data()), res_ptr1, total * sizeof(cudacd), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaGetLastError()); // 检查可能的内存越界错误
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<cudacd*>(out2.data()), res_ptr2, total * sizeof(cudacd), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaGetLastError()); // 检查可能的内存越界错误

    auto wall_end = std::chrono::high_resolution_clock::now();

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    std::chrono::duration<float, std::milli> wall_ms = wall_end - wall_start;

    // printf("[TIME] GPU Kernel Pure Exec: %.3f ms\n", milliseconds);
    // printf("[TIME] Total Wall Time (Inc. Memcpy): %.3f ms\n", wall_ms.count());

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return true;
}

bool GMTIProcessor::cuda_stage_align_async(int skip, size_t Na, size_t Nr) {
    if (gpu_ptrs_.d1 == nullptr) {
        std::cerr << "CRITICAL: gpu_ptrs_.d1 is NULL before kernel launch!" << std::endl;
        return false;
    }
    
    int total = static_cast<int>(Na * Nr);
    int threads = 256;
    int blocks = (total + threads - 1) / threads;

    // 执行对齐核函数：从 d1/d2 到 a1/a2
    align_two_channels_fast_cuda<<<blocks, threads, 0, stream_compute_>>>(
        (cudacd*)gpu_ptrs_.d1, (cudacd*)gpu_ptrs_.d2, 
        (cudacd*)gpu_ptrs_.a1, (cudacd*)gpu_ptrs_.a2, 
        skip, (int)Na, (int)Nr
    );
    return true;
}

bool GMTIProcessor::cuda_stage_fft_async(size_t Na, size_t Nr) {
    // 执行 FFT：a1/a2 原地 (In-place) 计算
    cufftExecC2C(cufft_plan_, 
                 (cufftComplex*)gpu_ptrs_.a1, 
                 (cufftComplex*)gpu_ptrs_.a1, 
                 CUFFT_FORWARD);
                 
    cufftExecC2C(cufft_plan_, 
                 (cufftComplex*)gpu_ptrs_.a2, 
                 (cufftComplex*)gpu_ptrs_.a2, 
                 CUFFT_FORWARD);
    return true;
}

bool GMTIProcessor::cuda_stage_dbs_async(float fa2, float prf, size_t Na, size_t Nr) {
    // 1. 计算偏移 k
    if (!(prf > 0.0f) || Na == 0 || Nr == 0) {
        return false;
    }
    int width = (int)Na;
    int center_num = (int)std::floor((fa2 + 0.5f * prf) / prf * width) + 1;
    int k_1b = center_num - width / 2;
    while (k_1b < 1) k_1b += width;
    while (k_1b > width) k_1b -= width;
    size_t k = (size_t)(k_1b - 1);

    // 2. 启动 DBS 核函数：从 a1/a2 到 t1/t2
    int total = (int)(Na * Nr);
    int threads = 256;
    int blocks = (total + threads - 1) / threads;

    dbs_center_by_fa_fast_cuda<<<blocks, threads, 0, stream_compute_>>>(
        (cudacd*)gpu_ptrs_.a1, (cudacd*)gpu_ptrs_.a2, 
        (cudacd*)gpu_ptrs_.t1, (cudacd*)gpu_ptrs_.t2, 
        (int)Na, (int)Nr, k
    );
    return true;
}

bool GMTIProcessor::cuda_stage_az_decimate_async(int W_orig, int M, int dec) {
    if (W_orig <= 0 || M <= 0 || dec <= 0 || gpu_ptrs_.d1 == nullptr || gpu_ptrs_.d2 == nullptr ||
        gpu_ptrs_.a1 == nullptr || gpu_ptrs_.a2 == nullptr) {
        return false;
    }
    if (W_orig % dec != 0) {
        return false;
    }
    const int W_new = W_orig / dec;
    if (dec == 1) {
        return true;
    }

    const size_t total_out = static_cast<size_t>(W_new) * static_cast<size_t>(M);
    const int threads = 256;
    const int blocks = static_cast<int>((total_out + threads - 1) / threads);
    az_decimate_two_channels_kernel<<<blocks, threads, 0, stream_compute_>>>(
        (const cudacd*)gpu_ptrs_.d1, (const cudacd*)gpu_ptrs_.d2,
        (cudacd*)gpu_ptrs_.a1, (cudacd*)gpu_ptrs_.a2,
        W_orig, M, dec, W_new);
    if (cudaGetLastError() != cudaSuccess) {
        ERR("az_decimate_two_channels_kernel launch failed");
        return false;
    }

    const size_t bytes = total_out * sizeof(cudacd);
    CUDA_CHECK(cudaMemcpyAsync(gpu_ptrs_.d1, gpu_ptrs_.a1, bytes, cudaMemcpyDeviceToDevice, stream_compute_));
    CUDA_CHECK(cudaMemcpyAsync(gpu_ptrs_.d2, gpu_ptrs_.a2, bytes, cudaMemcpyDeviceToDevice, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));
    return true;
}

bool GMTIProcessor::cuda_scale_channel2_async(float coef, size_t total) {
    if (gpu_ptrs_.d2 == nullptr || total == 0) {
        return false;
    }
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    scale_complex_kernel<<<blocks, threads, 0, stream_compute_>>>((cudacd*)gpu_ptrs_.d2, total, coef);
    if (cudaGetLastError() != cudaSuccess) {
        ERR("scale_complex_kernel launch failed");
        return false;
    }
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));
    return true;
}

bool GMTIProcessor::cuda_compute_fd_ctr_from_d1(int k, const Config& cfg, double& fd_ctr) {
    fd_ctr = 0.0;
    const int Na = effectivePulseNum(cfg);
    const int Nr = cfg.rg_len;
    if (k <= 0 || Na <= k || Nr <= 0 || !(cfg.PRF > 0.0) || gpu_ptrs_.d1 == nullptr) {
        return false;
    }

    const int threads = 256;
    const int blocks = 256;
    cudacd* d_block_sums = nullptr;
    if (cudaMalloc(&d_block_sums, static_cast<size_t>(blocks) * sizeof(cudacd)) != cudaSuccess) {
        return false;
    }
    adjacent_corr_reduce_kernel<<<blocks, threads, threads * sizeof(cudacd), stream_compute_>>>(
        (const cudacd*)gpu_ptrs_.d1, d_block_sums, k, Na, Nr);
    if (cudaGetLastError() != cudaSuccess) {
        cudaFree(d_block_sums);
        ERR("adjacent_corr_reduce_kernel launch failed");
        return false;
    }

    std::vector<cudacd> h_sums(blocks);
    CUDA_CHECK(cudaMemcpyAsync(h_sums.data(), d_block_sums,
                               static_cast<size_t>(blocks) * sizeof(cudacd),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));
    cudaFree(d_block_sums);

    double re = 0.0;
    double im = 0.0;
    for (const auto& v : h_sums) {
        re += static_cast<double>(v.x);
        im += static_cast<double>(v.y);
    }
    fd_ctr = (cfg.PRF / (2.0 * M_PI)) * gmti::trig_lut::atan2(im, re);
    return true;
}

bool GMTIProcessor::cuda_stage_rg_sum_async(const Config& cfg) {
    int Nr = cfg.rg_len;
    int threads = 256;
    int blocks = (Nr + threads - 1) / threads;

    // 使用刚才在 initcuFFTPlans 中算好的指针
    // 注意：输入通常是 a1, a2 (对齐+FFT+DBS 后的结果)
    compute_phase_sum_kernel<<<blocks, threads, 0, stream_compute_>>>(
        (const cuFloatComplex*)gpu_ptrs_.t1,
        (const cuFloatComplex*)gpu_ptrs_.t2,
        (cuFloatComplex*)gpu_ptrs_.d_rg_sums,
        cfg.az_st, cfg.az_ed,
        effectivePulseNum(cfg), cfg.rg_len,
        cfg.rg_st, cfg.rg_ed
    );

    // 检查核函数启动是否有误
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        ERR("Kernel launch failed (rg_sum): " << cudaGetErrorString(err));
        return false;
    }
    
    return true;
}

/**
 * @brief 距离向相位校正 (GPU 异步实现)
 * @param phi_fit CPU 拟合出的相位曲线 (长度为 Nr)
 * @param Na 脉冲数
 * @param Nr 距离门数
 */
bool GMTIProcessor::cuda_apply_rg_correction_async(const std::vector<float>& phi_fit, int Na, int Nr) {
    if (phi_fit.empty() || d_phi_fit_ == nullptr) return false;

    // 1. 拷贝到专属区域，绝无冲突
    cudaMemcpyAsync(d_phi_fit_, phi_fit.data(), Nr * sizeof(float), 
                    cudaMemcpyHostToDevice, stream_compute_);

    // 2. 启动 Kernel
    int threads = 256;
    int blocks = (Nr + threads - 1) / threads;
    const gmti::trig_lut_device::TrigLutConfig trig_cfg =
        gmtiGetDeviceTrigLutConfig();
    apply_phase_correction_kernel<<<blocks, threads, 0, stream_compute_>>>(
        (cuFloatComplex*)gpu_ptrs_.a1, 
        (cuFloatComplex*)gpu_ptrs_.a2,
        d_phi_fit_, // 使用专属指针
        Na, Nr,
        trig_cfg
    );
    return true;
}

bool GMTIProcessor::cuda_download_rg_sums_sync(std::vector<std::complex<float>>& h_sums, int Nr) {
    h_sums.resize(Nr);
    
    // 异步拷贝
    cudaMemcpyAsync(h_sums.data(), 
                    gpu_ptrs_.d_rg_sums, 
                    Nr * sizeof(std::complex<float>), 
                    cudaMemcpyDeviceToHost, 
                    stream_compute_);
    
    // 强制同步流，确保 CPU 拿到完整数据后再进行后续的 RANSAC 拟合
    cudaStreamSynchronize(stream_compute_);
    
    return true;
}

bool GMTIProcessor::cuda_download_phase_map(std::vector<float>& phase_map, size_t total) {
    if (total == 0) return false;

    float* d_phase = nullptr;
    CUDA_CHECK(cudaMalloc(&d_phase, total * sizeof(float)));

    int threads = 256;
    int blocks = static_cast<int>((total + threads - 1) / threads);
    const gmti::trig_lut_device::TrigLutConfig trig_cfg =
        gmtiGetDeviceTrigLutConfig();
    phase_map_kernel<<<blocks, threads, 0, stream_compute_>>>(
        (const cuFloatComplex*)gpu_ptrs_.t1,
        (const cuFloatComplex*)gpu_ptrs_.t2,
        d_phase,
        total,
        trig_cfg);
    CUDA_CHECK(cudaGetLastError());

    phase_map.resize(total);
    CUDA_CHECK(cudaMemcpyAsync(phase_map.data(), d_phase, total * sizeof(float),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));

    CUDA_CHECK(cudaFree(d_phase));
    return true;
}

bool GMTIProcessor::cuda_download_csi_sync(std::vector<std::complex<float>> &out, size_t total)
{
    if (total == 0) return false;
    if (gpu_ptrs_.csi == nullptr) return false;

    out.resize(total);
    CUDA_CHECK(cudaMemcpyAsync(out.data(), gpu_ptrs_.csi, total * sizeof(cudacd),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));
    return true;
}

bool GMTIProcessor::dpca_cfar2_fast_cuda(const std::vector<std::complex<float>> &CSI_out,
                                         int band_st, int band_ed,
                                         float pf, int c_num, int b_num, const std::string &type,
                                         const Config &cfg,
                                         std::vector<float> &mydata)
{
    const int H = effectivePulseNum(cfg);
    const int W = cfg.rg_len;
    if (H <= 0 || W <= 0) return false;

    const size_t total = static_cast<size_t>(H) * static_cast<size_t>(W);
    if (!CSI_out.empty() && CSI_out.size() != total) {
        std::cerr << "dpca_cfar2_cuda: CSI_out size mismatch\n";
        return false;
    }
    if (!(pf > 0.0f && pf < 1.0f)) {
        std::cerr << "dpca_cfar2_cuda: pf should be in (0,1)\n";
        return false;
    }

    const int g = std::max(0, c_num);
    const int b = std::max(1, b_num);
    const int R = g + b;
    const int win = 2 * R + 1;
    const int cut = 2 * g + 1;
    const int total_bg = win * win - cut * cut;
    if (total_bg <= 0) {
        std::cerr << "dpca_cfar2_cuda: invalid background count\n";
        return false;
    }

    std::string ty = type;
    for (auto &ch : ty) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    if (ty != "CA" && ty != "GO" && ty != "OS") ty = "CA";
    if (ty == "OS") {
        std::cerr << "dpca_cfar2_cuda: OS not implemented, fallback to CA.\n";
        ty = "CA";
    }
    const bool use_go = (ty == "GO");

    const float alpha = static_cast<float>(total_bg) * (std::pow(pf, -1.0f / static_cast<float>(total_bg)) - 1.0f);

    cuFloatComplex* d_csi = nullptr;
    bool owns_csi = false;
    cuFloatComplex* d_detect = nullptr;
    float* d_power = nullptr;
    float* d_row_prefix = nullptr;
    float* d_integral = nullptr;
    float* d_mydata = nullptr;

    if (gpu_ptrs_.csi != nullptr) {
        d_csi = reinterpret_cast<cuFloatComplex*>(gpu_ptrs_.csi);
    } else {
        CUDA_CHECK(cudaMalloc(&d_csi, total * sizeof(cudacd)));
        owns_csi = true;
    }
    CUDA_CHECK(cudaMalloc(&d_detect, total * sizeof(cudacd)));
    CUDA_CHECK(cudaMalloc(&d_power, total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_row_prefix, total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_integral, (static_cast<size_t>(H) + 1) * (static_cast<size_t>(W) + 1) * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_mydata, total * sizeof(float)));

    if (owns_csi) {
        if (CSI_out.empty()) {
            std::cerr << "dpca_cfar2_cuda: CSI_out is empty and no device CSI buffer is available.\n";
            CUDA_CHECK(cudaFree(d_csi));
            return false;
        }
        CUDA_CHECK(cudaMemcpyAsync(d_csi, CSI_out.data(), total * sizeof(cudacd),
                                   cudaMemcpyHostToDevice, stream_compute_));
    }

    band_st = std::max(0, std::min(band_st, H - 1));
    band_ed = std::max(0, std::min(band_ed, H - 1));

    int threads = 256;
    int blocks = static_cast<int>((total + threads - 1) / threads);
    mix_detect_data_kernel<<<blocks, threads, 0, stream_compute_>>>(
        (const cuFloatComplex*)gpu_ptrs_.t2,
        d_csi,
        d_detect,
        H,
        W,
        band_st,
        band_ed);
    CUDA_CHECK(cudaGetLastError());

    power_map_kernel<<<blocks, threads, 0, stream_compute_>>>(
        d_detect,
        d_power,
        total);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemsetAsync(d_integral, 0, (static_cast<size_t>(H) + 1) * (static_cast<size_t>(W) + 1) * sizeof(float), stream_compute_));

    int row_blocks = (H + threads - 1) / threads;
    row_prefix_kernel<<<row_blocks, threads, 0, stream_compute_>>>(
        d_power,
        d_row_prefix,
        H,
        W);
    CUDA_CHECK(cudaGetLastError());

    int col_blocks = (W + threads - 1) / threads;
    col_prefix_to_integral_kernel<<<col_blocks, threads, 0, stream_compute_>>>(
        d_row_prefix,
        d_integral,
        H,
        W);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemsetAsync(d_mydata, 0, total * sizeof(float), stream_compute_));

    cfar_detect_kernel<<<blocks, threads, 0, stream_compute_>>>(
        d_detect,
        d_integral,
        d_mydata,
        H,
        W,
        g,
        b,
        total_bg,
        alpha,
        use_go);
    CUDA_CHECK(cudaGetLastError());

    mydata.resize(total);
    CUDA_CHECK(cudaMemcpyAsync(mydata.data(), d_mydata, total * sizeof(float),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));

    if (owns_csi) {
        CUDA_CHECK(cudaFree(d_csi));
    }
    CUDA_CHECK(cudaFree(d_detect));
    CUDA_CHECK(cudaFree(d_power));
    CUDA_CHECK(cudaFree(d_row_prefix));
    CUDA_CHECK(cudaFree(d_integral));
    CUDA_CHECK(cudaFree(d_mydata));

    return true;
}

bool GMTIProcessor::cluster_filter_gap_phase_cuda(const std::vector<float> &mydata,
                                                  const std::vector<float> &phase_map,
                                                  int min_points,
                                                  int max_gap,
                                                  float maxPhaseStd,
                                                  const Config &cfg,
                                                  std::vector<float> &refined_data,
                                                  std::vector<int> &prow_new,
                                                  std::vector<int> &pcol_new,
                                                  std::vector<float> &phaseStdList)
{
    const int H = effectivePulseNum(cfg);
    const int W = cfg.rg_len;
    if (H <= 0 || W <= 0) return false;

    const int total = H * W;
    if ((int)mydata.size() != total) return false;
    if (!phase_map.empty() && (int)phase_map.size() != total) return false;

    if (min_points < 1) min_points = 1;
    if (max_gap < 1) max_gap = 1;

    const bool use_phase = (!phase_map.empty() && maxPhaseStd > 0.0f);

    float* d_mydata = nullptr;
    float* d_phase = nullptr;
    int* d_labels_a = nullptr;
    int* d_labels_b = nullptr;
    int* d_changed = nullptr;
    int* d_counts = nullptr;
    float* d_sum_cos = nullptr;
    float* d_sum_sin = nullptr;
    float* d_max_power = nullptr;
    int* d_max_idx = nullptr;
    float* d_mean_phi = nullptr;
    float* d_sum_sq = nullptr;
    float* d_phase_std = nullptr;

    CUDA_CHECK(cudaMalloc(&d_mydata, total * sizeof(float)));
    CUDA_CHECK(cudaMemcpyAsync(d_mydata, mydata.data(), total * sizeof(float),
                               cudaMemcpyHostToDevice, stream_compute_));

    if (use_phase) {
        CUDA_CHECK(cudaMalloc(&d_phase, total * sizeof(float)));
        CUDA_CHECK(cudaMemcpyAsync(d_phase, phase_map.data(), total * sizeof(float),
                                   cudaMemcpyHostToDevice, stream_compute_));
    }

    CUDA_CHECK(cudaMalloc(&d_labels_a, total * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_labels_b, total * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_changed, sizeof(int)));

    int threads = 256;
    int blocks = (total + threads - 1) / threads;

    init_labels_kernel<<<blocks, threads, 0, stream_compute_>>>(d_mydata, d_labels_a, total);
    CUDA_CHECK(cudaGetLastError());

    int max_iters = H + W;
    for (int it = 0; it < max_iters; ++it) {
        CUDA_CHECK(cudaMemsetAsync(d_changed, 0, sizeof(int), stream_compute_));
        propagate_labels_kernel<<<blocks, threads, 0, stream_compute_>>>(
            d_labels_a, d_labels_b, H, W, max_gap, d_changed);
        CUDA_CHECK(cudaGetLastError());

        int h_changed = 0;
        CUDA_CHECK(cudaMemcpyAsync(&h_changed, d_changed, sizeof(int),
                                   cudaMemcpyDeviceToHost, stream_compute_));
        CUDA_CHECK(cudaStreamSynchronize(stream_compute_));

        std::swap(d_labels_a, d_labels_b);
        if (!h_changed) break;
    }

    CUDA_CHECK(cudaMalloc(&d_counts, total * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_sum_cos, total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_sum_sin, total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_max_power, total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_max_idx, total * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_mean_phi, total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_sum_sq, total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_phase_std, total * sizeof(float)));

    CUDA_CHECK(cudaMemsetAsync(d_counts, 0, total * sizeof(int), stream_compute_));
    CUDA_CHECK(cudaMemsetAsync(d_sum_cos, 0, total * sizeof(float), stream_compute_));
    CUDA_CHECK(cudaMemsetAsync(d_sum_sin, 0, total * sizeof(float), stream_compute_));
    CUDA_CHECK(cudaMemsetAsync(d_sum_sq, 0, total * sizeof(float), stream_compute_));
    CUDA_CHECK(cudaMemsetAsync(d_max_power, 0, total * sizeof(float), stream_compute_));
    CUDA_CHECK(cudaMemsetAsync(d_max_idx, 0xff, total * sizeof(int), stream_compute_));

    const gmti::trig_lut_device::TrigLutConfig trig_cfg =
        gmtiGetDeviceTrigLutConfig();
    accumulate_stats_kernel<<<blocks, threads, 0, stream_compute_>>>(
        d_mydata,
        use_phase ? d_phase : d_mydata,
        d_labels_a,
        total,
        d_counts,
        d_sum_cos,
        d_sum_sin,
        d_max_power,
        d_max_idx,
        use_phase,
        trig_cfg);
    CUDA_CHECK(cudaGetLastError());

    if (use_phase) {
        compute_mean_phase_kernel<<<blocks, threads, 0, stream_compute_>>>(
            d_counts, d_sum_cos, d_sum_sin, d_mean_phi, total, trig_cfg);
        CUDA_CHECK(cudaGetLastError());

        accumulate_phase_var_kernel<<<blocks, threads, 0, stream_compute_>>>(
            d_phase, d_labels_a, d_mean_phi, total, d_sum_sq);
        CUDA_CHECK(cudaGetLastError());

        compute_phase_std_kernel<<<blocks, threads, 0, stream_compute_>>>(
            d_counts, d_sum_sq, d_phase_std, total);
        CUDA_CHECK(cudaGetLastError());
    } else {
        CUDA_CHECK(cudaMemsetAsync(d_phase_std, 0, total * sizeof(float), stream_compute_));
    }

    std::vector<int> labels(total);
    std::vector<int> counts(total);
    std::vector<int> max_idx(total);
    std::vector<float> phase_std(total);

    CUDA_CHECK(cudaMemcpyAsync(labels.data(), d_labels_a, total * sizeof(int),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaMemcpyAsync(counts.data(), d_counts, total * sizeof(int),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaMemcpyAsync(max_idx.data(), d_max_idx, total * sizeof(int),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaMemcpyAsync(phase_std.data(), d_phase_std, total * sizeof(float),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));

    std::vector<uint8_t> valid_label(total, 0);
    for (int i = 0; i < total; ++i) {
        if (counts[i] < min_points) continue;
        if (use_phase && phase_std[i] > maxPhaseStd) continue;
        valid_label[i] = 1;
    }

    refined_data.assign(static_cast<size_t>(total), 0.0f);
    prow_new.clear();
    pcol_new.clear();
    phaseStdList.clear();

    for (int i = 0; i < total; ++i) {
        int label = labels[i];
        if (label >= 0 && valid_label[label]) {
            refined_data[static_cast<size_t>(i)] = mydata[static_cast<size_t>(i)];
        }
    }

    for (int i = 0; i < total; ++i) {
        if (!valid_label[i]) continue;
        int best = max_idx[i];
        if (best < 0 || best >= total) continue;
        int r = best / W;
        int c = best - r * W;
        prow_new.push_back(r);
        pcol_new.push_back(c);
        if (use_phase) phaseStdList.push_back(phase_std[i]);
    }

    cudaFree(d_mydata);
    if (d_phase) cudaFree(d_phase);
    cudaFree(d_labels_a);
    cudaFree(d_labels_b);
    cudaFree(d_changed);
    cudaFree(d_counts);
    cudaFree(d_sum_cos);
    cudaFree(d_sum_sin);
    cudaFree(d_max_power);
    cudaFree(d_max_idx);
    cudaFree(d_mean_phi);
    cudaFree(d_sum_sq);
    cudaFree(d_phase_std);

    return true;
}

bool GMTIProcessor::target_select_cuda(const std::vector<int> &prow,
                                       const std::vector<int> &pcol,
                                       const Config &cfg,
                                       GMTIOutput::Detect &S_target)
{
    const int H = effectivePulseNum(cfg);
    const int W = cfg.rg_len;
    const size_t total = static_cast<size_t>(H) * static_cast<size_t>(W);
    if (H <= 0 || W <= 0) return false;
    if (prow.size() != pcol.size()) return false;
    if (gpu_ptrs_.t1 == nullptr) return false;

    const int n = static_cast<int>(prow.size());
    S_target.prow.clear();
    S_target.pcol.clear();
    if (n == 0) return true;

    thrust::device_vector<int> d_prow(prow.begin(), prow.end());
    thrust::device_vector<int> d_pcol(pcol.begin(), pcol.end());
    thrust::device_vector<float> d_amp(n);

    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    compute_amp_kernel<<<blocks, threads, 0, stream_compute_>>>(
        (const cuFloatComplex*)gpu_ptrs_.t1,
        thrust::raw_pointer_cast(d_prow.data()),
        thrust::raw_pointer_cast(d_pcol.data()),
        thrust::raw_pointer_cast(d_amp.data()),
        n,
        H,
        W);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));

    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d_prow.begin(), d_pcol.begin()));
    auto zip_end = thrust::make_zip_iterator(thrust::make_tuple(d_prow.end(), d_pcol.end()));
    thrust::sort_by_key(d_amp.begin(), d_amp.end(), zip_begin, thrust::greater<float>());

    auto valid_count = thrust::count_if(d_amp.begin(), d_amp.end(), NonNegativeAmp());

    if (valid_count == 0) return true;

    S_target.prow.resize(valid_count);
    S_target.pcol.resize(valid_count);
    thrust::copy(d_prow.begin(), d_prow.begin() + valid_count, S_target.prow.begin());
    thrust::copy(d_pcol.begin(), d_pcol.begin() + valid_count, S_target.pcol.begin());

    return true;
}

bool GMTIProcessor::clutter_cancel_38_paper_1_p38_cuda(
    const std::vector<double>& y_faAxis,
    int az_st, int rg_st, int az_ed, int rg_ed,
    const Config& cfg,
    std::array<double,2>& p_38,
    std::vector<double>* phase_tra_38)
{
    (void)rg_st;
    (void)rg_ed;
    const size_t Na = static_cast<size_t>(effectivePulseNum(cfg));
    const size_t Nr = static_cast<size_t>(cfg.rg_len);

    if (y_faAxis.size() != Na) return false;

    az_st = std::max(0, std::min(az_st, int(Na) - 1));
    az_ed = std::max(0, std::min(az_ed, int(Na) - 1));
    if (az_st > az_ed) return false;

    // 确认是否需要
    if (!initcuFFTPlans(cfg)) return false;

    std::vector<double> phase_tra_local;
    std::vector<double>& phase_tra = phase_tra_38 ? *phase_tra_38 : phase_tra_local;

    float* d_phase = nullptr;
    CUDA_CHECK(cudaMalloc(&d_phase, Na * sizeof(float)));
    CUDA_CHECK(cudaMemsetAsync(d_phase, 0, Na * sizeof(float), stream_compute_));

    int threads = 256;
    int blocks = (static_cast<int>(Na) + threads - 1) / threads;
    const gmti::trig_lut_device::TrigLutConfig trig_cfg =
        gmtiGetDeviceTrigLutConfig();
    compute_row_phase_kernel<<<blocks, threads, 0, stream_compute_>>>(
        (const cuFloatComplex*)gpu_ptrs_.t1,
        (const cuFloatComplex*)gpu_ptrs_.t2,
        d_phase,
        (int)Na,
        (int)Nr,
        az_st,
        az_ed,
        trig_cfg);
    CUDA_CHECK(cudaGetLastError());

    std::vector<float> phase_tra_f(Na, 0.0f);
    CUDA_CHECK(cudaMemcpyAsync(phase_tra_f.data(), d_phase, Na * sizeof(float),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));
    CUDA_CHECK(cudaFree(d_phase));

    phase_tra.resize(Na);
    for (size_t i = 0; i < Na; ++i) phase_tra[i] = static_cast<double>(phase_tra_f[i]);

    auto unwrap_inplace = [](std::vector<double>& v) {
        for (size_t i = 1; i < v.size(); ++i) {
            double d = v[i] - v[i - 1];
            if (d > M_PI) v[i] -= 2 * M_PI;
            if (d < -M_PI) v[i] += 2 * M_PI;
        }
    };

    int mid_num = std::max(1, std::min(int(std::round(double(az_ed - az_st + 1))), int(Na)));
    const double mid_phase = phase_tra[size_t(mid_num - 1)];
    unwrap_inplace(phase_tra);
    const double delta = mid_phase - phase_tra[size_t(mid_num - 1)];
    for (double& v : phase_tra) v += delta;

    const int M = az_ed - az_st + 1;
    std::vector<double> phase_tra_38_cut(M);
    std::vector<double> row_fa_cut(M);
    for (int i = 0; i < M; ++i) {
        phase_tra_38_cut[i] = phase_tra[size_t(az_st + i)];
        row_fa_cut[i] = y_faAxis[size_t(az_st + i)];
    }

    auto linfit = [](const std::vector<double>& x, const std::vector<double>& y, double& k, double& b) -> bool {
        if (x.size() != y.size() || x.size() < 2) return false;
        const size_t n = x.size();
        double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
        for (size_t i = 0; i < n; ++i) {
            Sx += x[i];
            Sy += y[i];
            Sxx += x[i] * x[i];
            Sxy += x[i] * y[i];
        }
        double det = n * Sxx - Sx * Sx;
        if (std::fabs(det) < 1e-12) return false;
        k = (n * Sxy - Sx * Sy) / det;
        b = (Sy * Sxx - Sx * Sxy) / det;
        return true;
    };

    double k = 0.0, b = 0.0;
    if (!linfit(row_fa_cut, phase_tra_38_cut, k, b)) return false;
    p_38 = {k, b};
    return true;
}

bool GMTIProcessor::clutter_cancel_38_paper_1_cuda(
    const std::vector<double>& y_faAxis,
    int az_st, int rg_st, int az_ed, int rg_ed,
    const Config& cfg,
    std::vector<std::complex<double>>& prosig_38,
    std::array<double,2>& p_38,
    std::vector<double>& phase_tra_38_cut,
    std::vector<double>& row_fa_cut)
{
    (void)rg_st;
    (void)rg_ed;
    const size_t Na = static_cast<size_t>(effectivePulseNum(cfg));
    const size_t Nr = static_cast<size_t>(cfg.rg_len);
    const size_t total = Na * Nr;

    if (y_faAxis.size() != Na) return false;

    az_st = std::max(0, std::min(az_st, int(Na) - 1));
    az_ed = std::max(0, std::min(az_ed, int(Na) - 1));
    if (az_st > az_ed) return false;

    std::vector<double> phase_tra_38;
    if (!clutter_cancel_38_paper_1_p38_cuda(y_faAxis, az_st, rg_st, az_ed, rg_ed, cfg, p_38, &phase_tra_38)) {
        return false;
    }

    const int M = az_ed - az_st + 1;
    phase_tra_38_cut.resize(M);
    row_fa_cut.resize(M);
    for (int i = 0; i < M; ++i) {
        phase_tra_38_cut[i] = phase_tra_38[size_t(az_st + i)];
        row_fa_cut[i] = y_faAxis[size_t(az_st + i)];
    }

    std::vector<float> y_faAxis_f(Na);
    for (size_t i = 0; i < Na; ++i) y_faAxis_f[i] = static_cast<float>(y_faAxis[i]);

    float* d_fa = nullptr;
    CUDA_CHECK(cudaMalloc(&d_fa, Na * sizeof(float)));
    CUDA_CHECK(cudaMemcpyAsync(d_fa, y_faAxis_f.data(), Na * sizeof(float),
                               cudaMemcpyHostToDevice, stream_compute_));

    if (gpu_ptrs_.csi == nullptr) {
        std::cerr << "CSI device buffer is not initialized.\n";
        CUDA_CHECK(cudaFree(d_fa));
        return false;
    }

    int threads = 256;
    int blocks = (static_cast<int>(total) + threads - 1) / threads;
    const gmti::trig_lut_device::TrigLutConfig trig_cfg =
        gmtiGetDeviceTrigLutConfig();
    clutter_cancel_kernel<<<blocks, threads, 0, stream_compute_>>>(
        (const cuFloatComplex*)gpu_ptrs_.t1,
        (const cuFloatComplex*)gpu_ptrs_.t2,
        d_fa,
        (cuFloatComplex*)gpu_ptrs_.csi,
        (int)Na,
        (int)Nr,
        static_cast<float>(p_38[0]),
        static_cast<float>(p_38[1]),
        trig_cfg);
    CUDA_CHECK(cudaGetLastError());

    std::vector<std::complex<float>> prosig_38_f(total);
    CUDA_CHECK(cudaMemcpyAsync(prosig_38_f.data(), gpu_ptrs_.csi, total * sizeof(cudacd),
                               cudaMemcpyDeviceToHost, stream_compute_));
    CUDA_CHECK(cudaStreamSynchronize(stream_compute_));
    prosig_38.resize(total);
    for (size_t i = 0; i < total; ++i) {
        prosig_38[i] = std::complex<double>(prosig_38_f[i].real(), prosig_38_f[i].imag());
    }

    CUDA_CHECK(cudaFree(d_fa));
    return true;
}

// 清理 CUDA 资源（在析构或程序结束时调用）
void GMTIProcessor::cleanupCUDAResources()
{
    if (d_workspace != nullptr) {
        cudaFree(d_workspace);
        d_workspace = nullptr;
        d_workspace_bytes = 0;
        DBG("CUDA workspace freed");
    }

    if (cufft_plan_ != -1) {
        cufftDestroy((cufftHandle)cufft_plan_);
        cufft_plan_ = -1;
        cached_Na_ = cached_Nr_ = 0;
        DBG("cuFFT plan destroyed");
    }
}
