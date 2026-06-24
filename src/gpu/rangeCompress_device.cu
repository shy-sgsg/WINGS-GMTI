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
#include "trig_lut.hpp"

__global__ void mul_Hf_dev(cuFloatComplex* data, const cuFloatComplex* Hf, int Lraw, int W) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)W * (size_t)Lraw;
    if (idx >= total) return;
    int f = idx % Lraw;
    cuFloatComplex a = data[idx];
    cuFloatComplex h = Hf[f];
    data[idx] = cuCmulf(a, h);
}

__global__ void pad_range_rows(const cuFloatComplex* in,
                               cuFloatComplex* out,
                               int Lraw,
                               int Nfft,
                               int W) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)W * (size_t)Nfft;
    if (idx >= total) return;
    int k = idx / Nfft;
    int n = idx % Nfft;
    out[idx] = n < Lraw
        ? in[(size_t)k * (size_t)Lraw + (size_t)n]
        : make_cuFloatComplex(0.0f, 0.0f);
}

__global__ void crop_and_scale(const cuFloatComplex* in,
                               cuFloatComplex* out,
                               int Nfft,
                               int M,
                               int W,
                               int crop,
                               int sample_step,
                               float scale) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total_out = (size_t)W * (size_t)M;
    if (idx >= total_out) return;
    int k = idx / M;
    int m = idx % M;
    size_t in_idx = (size_t)k * (size_t)Nfft +
                    (size_t)crop + (size_t)m * (size_t)sample_step;
    cuFloatComplex v = in[in_idx];
    out[idx] = make_cuFloatComplex(v.x * scale, v.y * scale);
}

bool GMTIProcessor::rangeCompressCUFFT_device(int Lraw, int M, const Config &cfg) {
    const int W = effectivePulseNum(cfg);
    const int Nfft = effectiveRangeFftLen(cfg);
    const int crop = cfg.range_crop_start;
    M = effectiveRangeCompressLen(cfg);
    if (W <= 0 || Lraw <= 0 || Nfft < Lraw || M <= 0 ||
        crop < 0 || crop + M > Nfft) return false;
    const int sample_step = usesRangeCropWindow(cfg) ? 1 : (Nfft / M);
    const size_t total_fft = (size_t)W * (size_t)Nfft;

    // device input pointer (raw) and output pointer (compressed)
    cuFloatComplex* d_in_src = reinterpret_cast<cuFloatComplex*>(gpu_ptrs_.d1);
    cuFloatComplex* d_out = reinterpret_cast<cuFloatComplex*>(gpu_ptrs_.d1); // write compressed result back into d1
    if (!d_in_src || !d_out) return false;

    // allocate temporary buffer and copy source raw into it (device->device)
    cuFloatComplex* d_buf = nullptr;
    size_t total_fft_bytes = total_fft * sizeof(cuFloatComplex);
    if (cudaMalloc((void**)&d_buf, total_fft_bytes) != cudaSuccess) return false;

    // build Hf on host and copy to device
    const double Kr = (cfg.Tr > 0.0) ? (cfg.Br / cfg.Tr) : 0.0;
    const double fs = cfg.fs;
    if (Kr == 0.0 || fs <= 0.0) return false;
    const double df = fs / double(Nfft);
    std::vector<cuFloatComplex> h_Hf(Nfft);
    for (int n=0;n<Nfft;++n) {
        double fn = (n <= (Nfft/2 - 1)) ? n*df : (n - Nfft)*df;
        double phase = M_PI * (fn*fn) / Kr;
        float re = (float)gmti::trig_lut::cos(phase);
        float im = (float)gmti::trig_lut::sin(phase);
        h_Hf[n] = make_cuFloatComplex(re, im);
    }
    cuFloatComplex* d_Hf = nullptr;
    if (cudaMalloc((void**)&d_Hf, (size_t)Nfft * sizeof(cuFloatComplex)) != cudaSuccess) return false;
    if (cudaMemcpyAsync(d_Hf, h_Hf.data(), (size_t)Nfft * sizeof(cuFloatComplex), cudaMemcpyHostToDevice, stream_compute_) != cudaSuccess) { cudaFree(d_Hf); return false; }

    // plan many
    cufftHandle plan;
    int nvec[1] = { Nfft };
    int inembed[1] = { Nfft };
    int onembed[1] = { Nfft };
    if (cufftPlanMany(&plan, 1, nvec, inembed, 1, Nfft,
                      onembed, 1, Nfft, CUFFT_C2C, W) != CUFFT_SUCCESS) {
        cudaFree(d_Hf); cudaFree(d_buf); return false;
    }
    CUFFT_CHECK(cufftSetStream(plan, stream_compute_));

    int threads = 256;
    int blocks = (total_fft + threads - 1) / threads;
    pad_range_rows<<<blocks, threads, 0, stream_compute_>>>(
        d_in_src, d_buf, Lraw, Nfft, W);
    if (cudaGetLastError() != cudaSuccess) {
        cufftDestroy(plan); cudaFree(d_Hf); cudaFree(d_buf); return false;
    }
    if (cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_buf), reinterpret_cast<cufftComplex*>(d_buf), CUFFT_FORWARD) != CUFFT_SUCCESS) { cufftDestroy(plan); cudaFree(d_Hf); cudaFree(d_buf); return false; }

    mul_Hf_dev<<<blocks, threads, 0, stream_compute_>>>(d_buf, d_Hf, Nfft, W);
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

    // Crop the configured contiguous IFFT window into d_out.
    float scale = 1.0f / (float)Nfft;
    size_t total_out_elems = (size_t)W * (size_t)M;
    int blocks2 = (total_out_elems + threads - 1) / threads;
    crop_and_scale<<<blocks2, threads, 0, stream_compute_>>>(
        d_buf, d_out, Nfft, M, W, crop, sample_step, scale);
    if (cudaGetLastError() != cudaSuccess) { 
        ERR("[GPU] crop_and_scale kernel launch failed");
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

bool GMTIProcessor::rangeCompressCUFFT_device_pair(int Lraw, int M, const Config &cfg) {
    const int W = effectivePulseNum(cfg);
    const int Nfft = effectiveRangeFftLen(cfg);
    const int crop = cfg.range_crop_start;
    M = effectiveRangeCompressLen(cfg);
    if (W <= 0 || Lraw <= 0 || Nfft < Lraw || M <= 0 ||
        crop < 0 || crop + M > Nfft) return false;
    const int sample_step = usesRangeCropWindow(cfg) ? 1 : (Nfft / M);
    const size_t total_fft = (size_t)W * (size_t)Nfft;
    const size_t total_out = (size_t)W * (size_t)M;

    cuFloatComplex* d_ch1 = reinterpret_cast<cuFloatComplex*>(gpu_ptrs_.d1);
    cuFloatComplex* d_ch2 = reinterpret_cast<cuFloatComplex*>(gpu_ptrs_.d2);
    if (!d_ch1 || !d_ch2) return false;

    const double Kr = (cfg.Tr > 0.0) ? (cfg.Br / cfg.Tr) : 0.0;
    const double fs = cfg.fs;
    if (Kr == 0.0 || fs <= 0.0) return false;

    cuFloatComplex* d_buf = nullptr;
    cuFloatComplex* d_Hf = nullptr;
    if (cudaMalloc((void**)&d_buf, total_fft * sizeof(cuFloatComplex)) != cudaSuccess) return false;
    if (cudaMalloc((void**)&d_Hf, (size_t)Nfft * sizeof(cuFloatComplex)) != cudaSuccess) {
        cudaFree(d_buf);
        return false;
    }

    const double df = fs / double(Nfft);
    std::vector<cuFloatComplex> h_Hf(Nfft);
    for (int n = 0; n < Nfft; ++n) {
        double fn = (n <= (Nfft / 2 - 1)) ? n * df : (n - Nfft) * df;
        double phase = M_PI * (fn * fn) / Kr;
        h_Hf[n] = make_cuFloatComplex((float)gmti::trig_lut::cos(phase), (float)gmti::trig_lut::sin(phase));
    }
    if (cudaMemcpyAsync(d_Hf, h_Hf.data(), (size_t)Nfft * sizeof(cuFloatComplex),
                        cudaMemcpyHostToDevice, stream_compute_) != cudaSuccess) {
        cudaFree(d_Hf);
        cudaFree(d_buf);
        return false;
    }

    cufftHandle plan;
    int nvec[1] = { Nfft };
    int inembed[1] = { Nfft };
    int onembed[1] = { Nfft };
    if (cufftPlanMany(&plan, 1, nvec, inembed, 1, Nfft,
                      onembed, 1, Nfft,
                      CUFFT_C2C, W) != CUFFT_SUCCESS) {
        cudaFree(d_Hf);
        cudaFree(d_buf);
        return false;
    }
    CUFFT_CHECK(cufftSetStream(plan, stream_compute_));

    const int threads = 256;
    const int blocks_in = (total_fft + threads - 1) / threads;
    const int blocks_out = (total_out + threads - 1) / threads;
    const float scale = 1.0f / (float)Nfft;

    auto process_one_channel = [&](cuFloatComplex* d_in_out) -> bool {
        pad_range_rows<<<blocks_in, threads, 0, stream_compute_>>>(
            d_in_out, d_buf, Lraw, Nfft, W);
        if (cudaGetLastError() != cudaSuccess) {
            return false;
        }
        if (cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_buf),
                         reinterpret_cast<cufftComplex*>(d_buf),
                         CUFFT_FORWARD) != CUFFT_SUCCESS) {
            return false;
        }
        mul_Hf_dev<<<blocks_in, threads, 0, stream_compute_>>>(d_buf, d_Hf, Nfft, W);
        if (cudaGetLastError() != cudaSuccess) {
            ERR("[GPU] mul_Hf_dev kernel launch failed");
            return false;
        }
        if (cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_buf),
                         reinterpret_cast<cufftComplex*>(d_buf),
                         CUFFT_INVERSE) != CUFFT_SUCCESS) {
            return false;
        }
        crop_and_scale<<<blocks_out, threads, 0, stream_compute_>>>(
            d_buf, d_in_out, Nfft, M, W, crop, sample_step, scale);
        if (cudaGetLastError() != cudaSuccess) {
            ERR("[GPU] crop_and_scale kernel launch failed");
            return false;
        }
        return true;
    };

    const bool ok = process_one_channel(d_ch1) && process_one_channel(d_ch2);
    cudaError_t sync_err = cudaStreamSynchronize(stream_compute_);
    if (sync_err != cudaSuccess) {
        ERR("[GPU] Stream synchronize failed: " << cudaGetErrorString(sync_err));
    }

    cufftDestroy(plan);
    cudaFree(d_Hf);
    cudaFree(d_buf);
    return ok && sync_err == cudaSuccess;
}
