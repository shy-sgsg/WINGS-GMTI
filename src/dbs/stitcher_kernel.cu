#include "dbs/stitcher_kernel.hpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <stdio.h>

__device__ inline float sample2D_linear(
    const float* amp,
    int width,
    int height,
    float x,
    float y)
{
    if (!amp || width <= 1 || height <= 1) {
        return 0.0f;
    }
    if (x < 0.0f || y < 0.0f || x > (float)(width - 1) || y > (float)(height - 1)) {
        return 0.0f;
    }

    x = fmaxf(0.0f, fminf(x, (float)(width - 1)));
    y = fmaxf(0.0f, fminf(y, (float)(height - 1)));

    int x0 = (int)floorf(x);
    int y0 = (int)floorf(y);
    int x1 = (x0 + 1 < width) ? (x0 + 1) : (width - 1);
    int y1 = (y0 + 1 < height) ? (y0 + 1) : (height - 1);
    float wx = x - (float)x0;
    float wy = y - (float)y0;

    float v00 = amp[y0 * width + x0];
    float v01 = amp[y0 * width + x1];
    float v10 = amp[y1 * width + x0];
    float v11 = amp[y1 * width + x1];

    return (1.0f - wy) * (1.0f - wx) * v00
         + (1.0f - wy) * wx * v01
         + wy * (1.0f - wx) * v10
         + wy * wx * v11;
}

__device__ inline float clamp_unit_device(float x)
{
    return x < -1.0f ? -1.0f : (x > 1.0f ? 1.0f : x);
}

__device__ inline float wrap180_device(float angle_deg)
{
    angle_deg = fmodf(angle_deg + 180.0f, 360.0f);
    if (angle_deg < 0.0f) angle_deg += 360.0f;
    return angle_deg - 180.0f;
}

__device__ inline bool in_beam_angle_gate_device(
    float fd, float lambda, float velocity,
    float beam_angle_deg, float gate_deg)
{
    if (!(velocity > 0.0f) || !(lambda > 0.0f) || !(gate_deg > 0.0f)) return true;
    const float ratio = clamp_unit_device(-fd * lambda / (2.0f * velocity));
    const float pixel_angle_deg = asinf(ratio) * 180.0f / 3.14159265358979323846f;
    return fabsf(wrap180_device(pixel_angle_deg - beam_angle_deg)) <= gate_deg;
}

__global__ void buildMosaicFullKernel(
    float* d_amp_mosaic,
    uint16_t* d_which_beam,
    const float* d_all_amps,
    const float* d_grid_x,
    const float* d_grid_y,
    const BeamDevParams* d_beams,
    int nx, int ny, int B, int M, int nEff,
    float R_min, float R_bin, float lambda, float Height, int flag_base,
#if DBS_USE_CUDA_TEXTURE_OBJECTS
    bool useTexInterp,
    cudaTextureObject_t texObj,
    int tex_beam_idx)
#else
    bool useTexInterp)
#endif
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < 0 || i >= nx || j < 0 || j >= ny) {
        // printf("[KERNEL] Thread out of bounds: i=%d j=%d nx=%d ny=%d\n", i, j, nx, ny);
        return;
    }

    float x_label = d_grid_x[i];
    float y_label = d_grid_y[j];

    float best_amp = 0.0f;
    uint16_t best_beam_idx = 0;

    for (int b = 0; b < B; ++b) {
#if DBS_USE_CUDA_TEXTURE_OBJECTS
        if (useTexInterp && b != tex_beam_idx) continue;
#endif
        BeamDevParams bp = d_beams[b];
        float dx = x_label - bp.x;
        float dy = y_label - bp.y;
        float c = bp.cos_j;
        float s = bp.sin_j;
        int flag = bp.flag;
        float x_tmp = 0.0f, y_tmp = 0.0f;
        switch (flag)
        {
        case 1: x_tmp = -dx * c - dy * s; y_tmp = -dx * s + dy * c; break;
        case 2: x_tmp = dx * c + dy * s; y_tmp = dx * s - dy * c; break;
        case 3: x_tmp = -dx * c + dy * s; y_tmp = dx * s + dy * c; break;
        case 4: x_tmp = dx * c - dy * s; y_tmp = -dx * s - dy * c; break;
        case 6: x_tmp = dx * c + dy * s; y_tmp = -dx * s + dy * c; break;
        case 5: x_tmp = -dx * c - dy * s; y_tmp = dx * s - dy * c; break;
        case 8: x_tmp = dx * c - dy * s; y_tmp = dx * s + dy * c; break;
        case 7: x_tmp = -dx * c + dy * s; y_tmp = -dx * s - dy * c; break;
        default: x_tmp = dx; y_tmp = dy; break;
        }
        float R = sqrtf(x_tmp * x_tmp + y_tmp * y_tmp + Height * Height);
        float fd = (R > 1e-6f) ? (2.0f * bp.V_feiji * x_tmp / (R * lambda)) : 0.0f;
        if (!in_beam_angle_gate_device(fd, lambda, bp.V_feiji,
                                       bp.angle_deg, bp.angle_gate_deg)) {
            continue;
        }
        // 严格物理范围判断
        float R_max = R_min + (M - 1) * R_bin;
        float fd_max = bp.min_fd + (nEff - 1) * bp.delta_fd;
        if (R < R_min || R > R_max || fd < bp.min_fd || fd > fd_max) {
            continue;
        }
        float R_idx = (R - R_min) / R_bin;
        float fd_idx = (fd - bp.min_fd) / (bp.delta_fd > 1e-9f ? bp.delta_fd : 1e-9f);
        // clamp 到合法范围
        R_idx = fmaxf(0.0f, fminf(R_idx, (float)(M - 1)));
        fd_idx = fmaxf(0.0f, fminf(fd_idx, (float)(nEff - 1)));
        if (R_idx >= 0.0f && R_idx < (float)(M - 1) && fd_idx >= 0.0f && fd_idx < (float)(nEff - 1)) {
            float amp_val = 0.0f;
            if (useTexInterp) {
#if DBS_USE_CUDA_TEXTURE_OBJECTS
                amp_val = tex2D<float>(texObj, R_idx + 0.5f, fd_idx + 0.5f);
#else
                const float* b_data = d_all_amps + (size_t)bp.offset;
                amp_val = sample2D_linear(b_data, M, nEff, R_idx, fd_idx);
#endif
            } else {
                const float* b_data = d_all_amps + (size_t)bp.offset;
                amp_val = sample2D_linear(b_data, M, nEff, R_idx, fd_idx);
            }
            if (fabsf(amp_val) > fabsf(best_amp)) {
                best_amp = amp_val;
                best_beam_idx = (uint16_t)(b + 1);
            }
        }
    }

    int out_idx = j * nx + i;
    if (out_idx < 0 || out_idx >= nx * ny) {
        // printf("[KERNEL] out_idx out of bounds: out_idx=%d nx=%d ny=%d\n", out_idx, nx, ny);
        return;
    }
    d_amp_mosaic[out_idx] = best_amp;
    d_which_beam[out_idx] = best_beam_idx;
}

// 包装函数实现
void launchBuildMosaicKernel(
    float* d_amp_mosaic, uint16_t* d_which_beam,
    const float* d_all_amps,
    const float* d_grid_x, const float* d_grid_y,
    const void* d_beams,
    int nx, int ny, int B, int M, int nEff,
    float R_min, float R_bin, float lambda, float Height, int flag_base,
    bool useTexInterp)
{
    dim3 block(16, 16);
    dim3 gsize((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
        // printf("[HOST] launchBuildMosaicKernel: nx=%d ny=%d B=%d M=%d nEff=%d block=(%d,%d) grid=(%d,%d) useTexInterp=%d\n", nx, ny, B, M, nEff, block.x, block.y, gsize.x, gsize.y, useTexInterp);
        // printf("[HOST] d_amp_mosaic=%p d_which_beam=%p d_all_amps=%p d_grid_x=%p d_grid_y=%p d_beams=%p\n", d_amp_mosaic, d_which_beam, d_all_amps, d_grid_x, d_grid_y, d_beams);
        cudaError_t err;
    if (!useTexInterp) {
        buildMosaicFullKernel<<<gsize, block>>>(
            d_amp_mosaic, d_which_beam, d_all_amps,
            d_grid_x, d_grid_y, (const BeamDevParams*)d_beams,
            nx, ny, B, M, nEff,
            R_min, R_bin, lambda, Height, flag_base,
            false
#if DBS_USE_CUDA_TEXTURE_OBJECTS
            , 0,
            0 // 非纹理路径，tex_beam_idx=0
#endif
        );
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            printf("[HOST] cudaDeviceSynchronize error: %s\n", cudaGetErrorString(err));
        }
        return;
    }
    // 多波位纹理对象数组插值
    extern void launchBuildMosaicKernelTexArray(
        float* d_amp_mosaic, uint16_t* d_which_beam,
        const float* d_all_amps,
        const float* d_grid_x, const float* d_grid_y,
        const void* d_beams,
        int nx, int ny, int B, int M, int nEff,
        float R_min, float R_bin, float lambda, float Height, int flag_base
    );
    launchBuildMosaicKernelTexArray(
        d_amp_mosaic, d_which_beam,
        d_all_amps,
        d_grid_x, d_grid_y,
        d_beams,
        nx, ny, B, M, nEff,
        R_min, R_bin, lambda, Height, flag_base
    );
}
