#pragma once
#include "dbs/stitcher_kernel.hpp"
#include <cuda_runtime.h>

#ifndef DBS_USE_CUDA_TEXTURE_OBJECTS
#define DBS_USE_CUDA_TEXTURE_OBJECTS 0
#endif

// 多波位 kernel 启动。CoreX 默认使用 const float** 指针数组；
// 需要保留旧纹理对象路径时，配置 DBS_USE_CUDA_TEXTURE_OBJECTS=1。
__global__ void buildMosaicFullKernelTexArray(
    float* d_amp_mosaic,
    uint16_t* d_which_beam,
    const float* d_all_amps,
    const float* d_grid_x,
    const float* d_grid_y,
    const BeamDevParams* d_beams,
    int nx, int ny, int B, int M, int nEff,
    float R_min, float R_bin, float lambda, float Height, int flag_base,
#if DBS_USE_CUDA_TEXTURE_OBJECTS
    cudaTextureObject_t* texObjArr // B 个纹理对象
#else
    const float** ampArr // B 个 device amp 指针
#endif
);

// host 启动器
void launchBuildMosaicKernelTexArray(
    float* d_amp_mosaic, uint16_t* d_which_beam,
    const float* d_all_amps,
    const float* d_grid_x, const float* d_grid_y,
    const void* d_beams,
    int nx, int ny, int B, int M, int nEff,
    float R_min, float R_bin, float lambda, float Height, int flag_base
);
