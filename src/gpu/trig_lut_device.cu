#include "trig_lut_device.cuh"

#include <cuda_runtime.h>
#include <iostream>

namespace {

gmti::trig_lut_device::TrigLutConfig g_device_config = {
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    gmti::trig_lut_device::kModeLut
};

} // namespace

extern "C" bool gmtiUploadTrigLutToDevice(const float* sinTable,
                                          const float* cosTable,
                                          const float* asinTable,
                                          const float* atanTable,
                                          int tableSize)
{
    if (tableSize != gmti::trig_lut_device::kTrigLutSize) {
        std::cerr << "[TRIG-LUT][ERR] unexpected table size: " << tableSize << std::endl;
        return false;
    }

    cudaError_t err = cudaSuccess;
    if (!g_device_config.sin_table) {
        float* ptr = nullptr;
        err = cudaMalloc(&ptr, tableSize * sizeof(float));
        if (err == cudaSuccess) g_device_config.sin_table = ptr;
    }
    if (err == cudaSuccess && !g_device_config.cos_table) {
        float* ptr = nullptr;
        err = cudaMalloc(&ptr, tableSize * sizeof(float));
        if (err == cudaSuccess) g_device_config.cos_table = ptr;
    }
    if (err == cudaSuccess && !g_device_config.asin_table) {
        float* ptr = nullptr;
        err = cudaMalloc(&ptr, (tableSize + 1) * sizeof(float));
        if (err == cudaSuccess) g_device_config.asin_table = ptr;
    }
    if (err == cudaSuccess && !g_device_config.atan_table) {
        float* ptr = nullptr;
        err = cudaMalloc(&ptr, (tableSize + 1) * sizeof(float));
        if (err == cudaSuccess) g_device_config.atan_table = ptr;
    }
    if (err == cudaSuccess) {
        err = cudaMemcpy(const_cast<float*>(g_device_config.sin_table), sinTable,
                         tableSize * sizeof(float), cudaMemcpyHostToDevice);
    }
    if (err == cudaSuccess) {
        err = cudaMemcpy(const_cast<float*>(g_device_config.cos_table), cosTable,
                         tableSize * sizeof(float), cudaMemcpyHostToDevice);
    }
    if (err == cudaSuccess) {
        err = cudaMemcpy(const_cast<float*>(g_device_config.asin_table), asinTable,
                         (tableSize + 1) * sizeof(float), cudaMemcpyHostToDevice);
    }
    if (err == cudaSuccess) {
        err = cudaMemcpy(const_cast<float*>(g_device_config.atan_table), atanTable,
                         (tableSize + 1) * sizeof(float), cudaMemcpyHostToDevice);
    }
    if (err != cudaSuccess) {
        std::cerr << "[TRIG-LUT][ERR] device table upload failed: "
                  << cudaGetErrorString(err) << std::endl;
        return false;
    }

    std::cout << "[TRIG-LUT] uploaded " << tableSize
              << "-entry sin/cos/asin/atan tables to device memory" << std::endl;
    return true;
}

extern "C" bool gmtiSetDeviceTrigMode(int mode)
{
    g_device_config.mode = mode;
    return true;
}

extern "C" gmti::trig_lut_device::TrigLutConfig gmtiGetDeviceTrigLutConfig()
{
    return g_device_config;
}
