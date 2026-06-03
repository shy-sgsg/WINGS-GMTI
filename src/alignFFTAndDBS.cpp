#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include <fftw3.h>
#include <iostream>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "GMTIProcessor.hpp"
typedef std::complex<double> cd;

// —— 工具：计算字节数
template <typename T>
static inline std::size_t bytes(std::size_t n) { return n * sizeof(T); }

static inline bool align_two_channels_fast(const std::vector<cd> &d1,
                                           const std::vector<cd> &d2,
                                           int skip,
                                           std::size_t Na, std::size_t Nr,
                                           std::vector<cd> &a1,
                                           std::vector<cd> &a2)
{
    const std::size_t total = Na * Nr;
    if (d1.size() != total || d2.size() != total)
        return false;

    a1.resize(total);
    a2.resize(total);

    // 默认完整复制
    std::memcpy(a1.data(), d1.data(), bytes<cd>(total));
    std::memcpy(a2.data(), d2.data(), bytes<cd>(total));

    if (skip == 0)
        return true;

    if (skip > 0)
    {
        const std::size_t s = static_cast<std::size_t>(skip);
        if (s > Na)
            return true;
        const std::size_t rows_copy = Na - s + 1;
        const std::size_t src_row0 = s - 1;

        std::memset(a1.data(), 0, bytes<cd>(total));
        const cd *src = d1.data() + src_row0 * Nr;
        cd *dst = a1.data();
        std::memcpy(dst, src, bytes<cd>(rows_copy * Nr));
    }
    else
    {
        const std::size_t s = static_cast<std::size_t>(-skip);
        if (s > Na)
            return true;
        const std::size_t rows_copy = Na - s + 1;
        const std::size_t src_row0 = s - 1;

        std::memset(a2.data(), 0, bytes<cd>(total));
        const cd *src = d2.data() + src_row0 * Nr;
        cd *dst = a2.data();
        std::memcpy(dst, src, bytes<cd>(rows_copy * Nr));
    }
    return true;
}

static inline void dbs_center_by_fa_fast(std::vector<cd> &A,
                                         std::vector<cd> &B,
                                         std::size_t Na, std::size_t Nr,
                                         double f_d)
{
    if (Na == 0 || Nr == 0)
        return;

    const int width = static_cast<int>(Na);
    const int center_num = static_cast<int>(std::floor((f_d + 1000.0) / 2000.0 * width)) + 1;

    int cstart = center_num - width / 2; // 1-based
    int k_1b = cstart;
    while (k_1b < 1)
        k_1b += width;
    while (k_1b > width)
        k_1b -= width;

    const std::size_t k = static_cast<std::size_t>(k_1b - 1); // 0-based
    if (k == 0)
        return;

    const std::size_t top_rows = Na - k; // [k .. Na-1]
    const std::size_t bot_rows = k;      // [0 .. k-1]

    std::vector<cd> tmpA(A.size()), tmpB(B.size());

    #ifdef _OPENMP
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            std::memcpy(tmpA.data(),                 A.data() + k * Nr, bytes<cd>(top_rows * Nr));
            std::memcpy(tmpA.data() + top_rows * Nr, A.data(),          bytes<cd>(bot_rows * Nr));
        }
        #pragma omp section
        {
            std::memcpy(tmpB.data(),                 B.data() + k * Nr, bytes<cd>(top_rows * Nr));
            std::memcpy(tmpB.data() + top_rows * Nr, B.data(),          bytes<cd>(bot_rows * Nr));
        }
    }
    #else
        std::memcpy(tmpA.data(),                 A.data() + k * Nr, bytes<cd>(top_rows * Nr));
        std::memcpy(tmpA.data() + top_rows * Nr, A.data(),          bytes<cd>(bot_rows * Nr));
        std::memcpy(tmpB.data(),                 B.data() + k * Nr, bytes<cd>(top_rows * Nr));
        std::memcpy(tmpB.data() + top_rows * Nr, B.data(),          bytes<cd>(bot_rows * Nr));
    #endif

    A.swap(tmpA);
    B.swap(tmpB);
}

// 主函数：对齐+FFT(FFTW)+fftshift+DBS
bool GMTIProcessor::alignFFTAndDBS(const std::vector<std::complex<double>> &data1,
                                   const std::vector<std::complex<double>> &data2,
                                   int skip,
                                   double fa2,
                                   const Config &cfg,
                                   std::vector<std::complex<double>> &out1,
                                   std::vector<std::complex<double>> &out2)
{
    const std::size_t Na = static_cast<std::size_t>(cfg.pulse_num);
    const std::size_t Nr = static_cast<std::size_t>(cfg.rg_len);
    const std::size_t total = Na * Nr;
    if (Na == 0 || Nr == 0) return false;
    if (data1.size() != total || data2.size() != total) return false;

    // 确保 ColFFTTranspose 已 init(Na, Nr)
    // 若你不想在外面手动调用 initFFT，可在此处懒初始化：
    // static bool inited = false; if (!inited) { colfft_.init(Na, Nr); inited = true; }
    // 这里按更显式的风格，要求外部先 initFFT(Na, Nr)。

    // 1) 通道对齐（fast）
    align_two_channels_fast(data1, data2, skip, Na, Nr, a1_, a2_);

    // 2) 列向 FFT（通过“转置 + 行向 batched FFT + 反转置”，已包含 (-1)^r 预乘等效 fftshift）
    out1.resize(total);
    out2.resize(total);
    if (!colfft_.execute(a1_.data(), out1.data())) return false;
    if (!colfft_.execute(a2_.data(), out2.data())) return false;

    // 3) DBS 中心化（两段 memcpy）
    dbs_center_by_fa_fast(out1, out2, Na, Nr, fa2);

    return true;
}

bool GMTIProcessor::initFFTPlans(const Config &cfg)
{
    const int H = cfg.pulse_num / cfg.pulse_dec;
    const int W = cfg.rg_len;
    int threads = std::max(1, (int)std::thread::hardware_concurrency());
    DBG("initFFTPlans: Using " << threads << " threads for FFTW");
    if (!colfft_.init(H, W, true, threads))
    {
        ERR("initFFTPlans: FFTW plan init failed (H=" << H << ", W=" << W << ")");
        return false;
    }
    return true;
}

bool GMTIProcessor::initcuFFTPlans(const Config &cfg) {
    // 使用 size_t 避免 32位 溢出
    size_t Na = static_cast<size_t>(cfg.pulse_num);
    size_t Nr = static_cast<size_t>(cfg.rg_len);
    size_t total = Na * Nr;
    
    // 显存中每个复数占 8 字节 (float2)
    size_t single_buf_bytes = total * sizeof(cudacd);
    size_t needed_bytes = 7 * single_buf_bytes + Nr * sizeof(cudacd); // 7 个缓冲区(CSI) + 距离向相位校正中间结果

    // --- 1. 预分配显存 (仅在必要时重新分配) ---
    if (d_workspace == nullptr || d_workspace_bytes < needed_bytes) {
        if (d_workspace) {
            cudaFree(d_workspace);
        }
        // 使用标准的 CUDA 错误检查
        cudaError_t err = cudaMalloc(&d_workspace, needed_bytes);
        if (err != cudaSuccess) {
            ERR("GPU memory allocation failed: " << cudaGetErrorString(err));
            return false;
        }
        d_workspace_bytes = needed_bytes;
        
        // 精确计算偏移指针
        uint8_t* base = reinterpret_cast<uint8_t*>(d_workspace);
        gpu_ptrs_.d1 = base + 0 * single_buf_bytes;
        gpu_ptrs_.d2 = base + 1 * single_buf_bytes;
        gpu_ptrs_.a1 = base + 2 * single_buf_bytes;
        gpu_ptrs_.a2 = base + 3 * single_buf_bytes;
        gpu_ptrs_.t1 = base + 4 * single_buf_bytes;
        gpu_ptrs_.t2 = base + 5 * single_buf_bytes;
        gpu_ptrs_.csi = base + 6 * single_buf_bytes;
        gpu_ptrs_.d_rg_sums = base + 7 * single_buf_bytes; // 紧跟在主缓冲区后面

        cudaMalloc(&d_phi_fit_, Nr * sizeof(float));
        
        DBG("CUDA Workspace re-allocated. Total: " << needed_bytes / (1024*1024) << " MB");
    }

    // --- 2. cuFFT 计划更新检查 ---
    if (cufft_plan_ != 0 && (int)Na == cached_Na_ && (int)Nr == cached_Nr_) {
        return true; 
    }

    if (cufft_plan_ != 0) {
        cufftDestroy(cufft_plan_);
        cufft_plan_ = 0;
    }

    // --- 3. 配置列向 FFT 布局 ---
    int rank = 1;           
    int n[] = { (int)Na };       
    int istride = (int)Nr;       
    int idist = 1;          
    int ostride = (int)Nr;       
    int odist = 1;
    
    // 注意：inembed/onembed 设置为 NULL 通常表示连续内存，
    // 但在 stride != 1 的情况下，传入 n 是更安全的做法。
    int inembed[] = { (int)Na }; 
    int onembed[] = { (int)Na };

    CUFFT_CHECK(cufftPlanMany(&cufft_plan_, 
                              rank, 
                              n, 
                              inembed, istride, idist, 
                              onembed, ostride, odist, 
                              CUFFT_C2C, 
                              (int)Nr)); 

    // --- 4. 绑定流 (确保异步) ---
    CUFFT_CHECK(cufftSetStream(cufft_plan_, stream_compute_));

    cached_Na_ = (int)Na;
    cached_Nr_ = (int)Nr;
    
    return true;
}
