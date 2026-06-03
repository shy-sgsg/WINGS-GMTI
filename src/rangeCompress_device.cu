// Device-only range compression: read raw data from gpu_ptrs_.d1 (size W x Lraw)
// and write compressed result into gpu_ptrs_.a1 (size W x M). Uses cuFFT and device kernels.

#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include "GMTIProcessor.hpp"

__global__ void mul_Hf_dev(cuFloatComplex* data, const cuFloatComplex* Hf, int Lraw, int W) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)W * (size_t)Lraw;
    if (idx >= total) return;
    int f = idx % Lraw;
    cuFloatComplex a = data[idx];
    cuFloatComplex h = Hf[f];
    data[idx] = cuCmulf(a, h);
}

__global__ void decimate_and_scale(const cuFloatComplex* in, cuFloatComplex* out, int Lraw, int M, int W, int Lraw2M, float scale) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total_out = (size_t)W * (size_t)M;
    if (idx >= total_out) return;
    int k = idx / M;
    int m = idx % M;
    size_t in_idx = (size_t)k * (size_t)Lraw + (size_t)m * (size_t)Lraw2M;
    cuFloatComplex v = in[in_idx];
    out[idx] = make_cuFloatComplex(v.x * scale, v.y * scale);
}

bool GMTIProcessor::rangeCompressCUFFT_device(int Lraw, int M, const Config &cfg) {
    const int W = cfg.pulse_num;
    if (W <= 0 || Lraw <= 0 || M <= 0) return false;
    const size_t total_in = (size_t)W * (size_t)Lraw;
    const size_t total_out = (size_t)W * (size_t)M;

    // device input pointer (raw) and output pointer (compressed)
    cuFloatComplex* d_in_src = reinterpret_cast<cuFloatComplex*>(gpu_ptrs_.d1);
    cuFloatComplex* d_out = reinterpret_cast<cuFloatComplex*>(gpu_ptrs_.d1); // write compressed result back into d1
    if (!d_in_src || !d_out) return false;

    // allocate temporary buffer and copy source raw into it (device->device)
    cuFloatComplex* d_buf = nullptr;
    size_t total_in_bytes = total_in * sizeof(cuFloatComplex);
    if (cudaMalloc((void**)&d_buf, total_in_bytes) != cudaSuccess) return false;
    if (cudaMemcpyAsync(d_buf, d_in_src, total_in_bytes, cudaMemcpyDeviceToDevice, stream_compute_) != cudaSuccess) { cudaFree(d_buf); return false; }

    // build Hf on host and copy to device
    const double Kr = (cfg.Tr > 0.0) ? (cfg.Br / cfg.Tr) : 0.0;
    const double fs = cfg.fs;
    if (Kr == 0.0 || fs <= 0.0) return false;
    const double df = fs / double(Lraw);
    std::vector<cuFloatComplex> h_Hf(Lraw);
    for (int n=0;n<Lraw;++n) {
        double fn = (n <= (Lraw/2 - 1)) ? n*df : (n - Lraw)*df;
        double phase = M_PI * (fn*fn) / Kr;
        float re = (float)std::cos(phase);
        float im = (float)std::sin(phase);
        h_Hf[n] = make_cuFloatComplex(re, im);
    }
    cuFloatComplex* d_Hf = nullptr;
    if (cudaMalloc((void**)&d_Hf, (size_t)Lraw * sizeof(cuFloatComplex)) != cudaSuccess) return false;
    if (cudaMemcpyAsync(d_Hf, h_Hf.data(), (size_t)Lraw * sizeof(cuFloatComplex), cudaMemcpyHostToDevice, stream_compute_) != cudaSuccess) { cudaFree(d_Hf); return false; }

    // plan many
    cufftHandle plan;
    int nvec[1] = { Lraw };
    if (cufftPlanMany(&plan, 1, nvec, &Lraw, 1, Lraw, &Lraw, 1, Lraw, CUFFT_C2C, W) != CUFFT_SUCCESS) { cudaFree(d_Hf); cudaFree(d_buf); return false; }
    CUFFT_CHECK(cufftSetStream(plan, stream_compute_));

    if (cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_buf), reinterpret_cast<cufftComplex*>(d_buf), CUFFT_FORWARD) != CUFFT_SUCCESS) { cufftDestroy(plan); cudaFree(d_Hf); cudaFree(d_buf); return false; }

    int threads = 256; int blocks = (total_in + threads - 1) / threads;
    mul_Hf_dev<<<blocks, threads, 0, stream_compute_>>>(d_buf, d_Hf, Lraw, W);
    if (cudaGetLastError() != cudaSuccess) { 
        ERR("[GPU] mul_Hf_dev kernel launch failed");
        cufftDestroy(plan); cudaFree(d_Hf); cudaFree(d_buf); 
        return false; 
    }

    cufftResult ifft_result = cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_buf), reinterpret_cast<cufftComplex*>(d_buf), CUFFT_INVERSE);
    if (ifft_result != CUFFT_SUCCESS) { 
        ERR("[GPU] cuFFT inverse FFT failed (error=" << ifft_result << ")");
        cufftDestroy(plan); cudaFree(d_Hf); cudaFree(d_buf); 
        return false; 
    }

    // decimate and scale into d_out (which is gpu_ptrs_.d1)
    int Lraw2M = Lraw / M;
    if (Lraw2M <= 0) {
        ERR("[GPU] Invalid decimation factor: Lraw=" << Lraw << " M=" << M);
        cufftDestroy(plan); cudaFree(d_Hf); cudaFree(d_buf); 
        return false;
    }
    float scale = 1.0f / (float)Lraw;
    size_t total_out_elems = (size_t)W * (size_t)M;
    int blocks2 = (total_out_elems + threads - 1) / threads;
    decimate_and_scale<<<blocks2, threads, 0, stream_compute_>>>(d_buf, d_out, Lraw, M, W, Lraw2M, scale);
    if (cudaGetLastError() != cudaSuccess) { 
        ERR("[GPU] decimate_and_scale kernel launch failed");
        cufftDestroy(plan); cudaFree(d_Hf); cudaFree(d_buf); 
        return false; 
    }

    cudaError_t sync_err = cudaStreamSynchronize(stream_compute_);
    if (sync_err != cudaSuccess) {
        ERR("[GPU] Stream synchronize failed: " << cudaGetErrorString(sync_err));
        cufftDestroy(plan); cudaFree(d_Hf); cudaFree(d_buf);
        return false;
    }

    // cleanup
    cufftDestroy(plan);
    cudaFree(d_Hf);
    cudaFree(d_buf);
    return true;
}
