// GPU range compression using cuFFT (float precision)
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include "GMTIProcessor.hpp"
#include "trig_lut.hpp"

static inline cuFloatComplex to_cu(const std::complex<float>& x) {
    return make_cuFloatComplex(x.real(), x.imag());
}

__global__ void mul_Hf_kernel(cuFloatComplex* data, const cuFloatComplex* Hf, int Lraw, int W)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)W * (size_t)Lraw;
    if (idx >= total) return;
    int f = idx % Lraw;
    cuFloatComplex a = data[idx];
    cuFloatComplex h = Hf[f];
    // multiply by Hf (no conjugate): data[idx] = data[idx] * h
    data[idx] = cuCmulf(a, h);
}

bool GMTIProcessor::rangeCompressCUFFT(const std::vector<std::complex<float>> &in,
                                       std::vector<std::complex<float>> &rc_out,
                                       const Config &cfg)
{
    const int W = effectivePulseNum(cfg);
    const int Lraw = cfg.pulse_len;
    const int M = cfg.rg_len;
    if (W <= 0 || Lraw <= 0 || M <= 0) return false;

    const size_t total = (size_t)W * (size_t)Lraw;

    // build HfSpec (float) from params: HfSpec[n] = exp(j * pi * f(n)^2 / Kr)
    const double Kr = (cfg.Tr > 0.0) ? (cfg.Br / cfg.Tr) : 0.0;
    const double fs = cfg.fs;
    std::vector<cuFloatComplex> h_buf(total);
    std::vector<cuFloatComplex> h_Hf(Lraw);
    if (Kr == 0.0 || fs <= 0.0) return false;
    const double df = fs / double(Lraw);
    for (int n=0;n<Lraw;++n) {
        double fn = (n <= (Lraw/2 - 1)) ? n*df : (n - Lraw)*df;
        double phase = M_PI * (fn*fn) / Kr;
        float re = (float)gmti::trig_lut::cos(phase);
        float im = (float)gmti::trig_lut::sin(phase);
        h_Hf[n] = make_cuFloatComplex(re, im);
    }
    for (int k=0;k<W;++k) {
        for (int n=0;n<Lraw;++n) {
            const auto &v = in[(size_t)k * Lraw + n];
            h_buf[(size_t)k * Lraw + n] = make_cuFloatComplex(v.real(), v.imag());
        }
    }

    cuFloatComplex *d_buf = nullptr;
    cuFloatComplex *d_Hf = nullptr;
    cudaError_t cerr;
    cerr = cudaMalloc((void**)&d_buf, total * sizeof(cuFloatComplex));
    if (cerr != cudaSuccess) { std::cerr<<"cudaMalloc d_buf failed"<<std::endl; return false; }
    cerr = cudaMalloc((void**)&d_Hf, (size_t)Lraw * sizeof(cuFloatComplex));
    if (cerr != cudaSuccess) { cudaFree(d_buf); std::cerr<<"cudaMalloc d_Hf failed"<<std::endl; return false; }

    cudaMemcpy(d_buf, h_buf.data(), total * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Hf, h_Hf.data(), (size_t)Lraw * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);

    // cufft plan: 1-D FFT length Lraw, howmany = W
    cufftHandle plan;
    int n[1] = { Lraw };
    int rank = 1;
    int istride = 1, ostride = 1, idist = Lraw, odist = Lraw;
    int inembed[1] = { Lraw }, onembed[1] = { Lraw };
    if (cufftPlanMany(&plan, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, W) != CUFFT_SUCCESS) {
        cudaFree(d_buf); cudaFree(d_Hf); std::cerr<<"cufftPlanMany failed"<<std::endl; return false;
    }
    CUFFT_CHECK(cufftSetStream(plan, stream_compute_));

    // forward
    if (cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_buf), reinterpret_cast<cufftComplex*>(d_buf), CUFFT_FORWARD) != CUFFT_SUCCESS) {
        cufftDestroy(plan); cudaFree(d_buf); cudaFree(d_Hf); std::cerr<<"cufftExecC2C forward failed"<<std::endl; return false;
    }

    // multiply by Hf (match CPU: X *= Hf[f])
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    mul_Hf_kernel<<<blocks, threads, 0, stream_compute_>>>(d_buf, d_Hf, Lraw, W);
    cudaError_t err = cudaGetLastError(); if (err != cudaSuccess) { cufftDestroy(plan); cudaFree(d_buf); cudaFree(d_Hf); std::cerr<<"mul_Hf_kernel failed"<<std::endl; return false; }

    // inverse
    if (cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(d_buf), reinterpret_cast<cufftComplex*>(d_buf), CUFFT_INVERSE) != CUFFT_SUCCESS) {
        cufftDestroy(plan); cudaFree(d_buf); cudaFree(d_Hf); std::cerr<<"cufftExecC2C inverse failed"<<std::endl; return false;
    }

    // download
    std::vector<cuFloatComplex> h_out(total);
    cudaMemcpyAsync(h_out.data(), d_buf, total * sizeof(cuFloatComplex), cudaMemcpyDeviceToHost, stream_compute_);
    cudaStreamSynchronize(stream_compute_);

    // scale and decimate to M
    rc_out.resize((size_t)W * M);
    int Lraw2M = Lraw / M;
    float scale = 1.0f / (float)Lraw;
    for (int k=0;k<W;++k) {
        for (int m=0;m<M;++m) {
            cuFloatComplex v = h_out[(size_t)k * Lraw + (size_t)m * Lraw2M];
            float re = v.x * scale;
            float im = v.y * scale;
            rc_out[(size_t)k * M + m] = std::complex<float>(re, im);
        }
    }

    // cleanup
    cufftDestroy(plan);
    cudaFree(d_buf);
    cudaFree(d_Hf);
    return true;
}
