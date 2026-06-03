#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// #define DEBUG // 调试开关cd 
# define ENABLE_TIMING // 启用计时

#include <iostream>
#include <cuComplex.h>
#include <cufft.h>
#include "config_structs.hpp"
#include "dbs/DbsFusionTypes.hpp"
#include <vector>
#include <string>
#include <cstdint>
#include <limits>
#include "tinyxml.h"
#include <array>
#include <cstddef>
#include <cuda_runtime.h>
#include "transposeFFT.hpp"

#ifdef DEBUG
  #define DBG(msg)  std::cout << "[DEBUG] " << __FILE__ << ":" << __LINE__ << " " << msg << std::endl
#else
  #define DBG(msg)  ((void)0)
#endif

#define ERR(msg)  std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << " " << msg << std::endl

#ifdef _OPENMP
#include <omp.h>         // OpenMP 支持
#endif

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

// ============= 计时宏 =============
// 编译时开启计时：g++ -DENABLE_TIMING ...
#ifdef ENABLE_TIMING
    #include <chrono>
    #define TIMING_START(name) \
        auto _timer_start_##name = std::chrono::high_resolution_clock::now()
  
    #define TIMING_END(name) \
        do { \
                auto _timer_end_##name = std::chrono::high_resolution_clock::now(); \
                auto _duration = std::chrono::duration_cast<std::chrono::milliseconds>(_timer_end_##name - _timer_start_##name); \
                std::cout << "[TIMING] " << #name << ": " << _duration.count() << " ms" << std::endl; \
        } while (0)
  
    #define TIMING_POINT(name) \
        do { \
                static auto _last_time = std::chrono::high_resolution_clock::now(); \
                auto _now = std::chrono::high_resolution_clock::now(); \
                auto _elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(_now - _last_time); \
                std::cout << "[TIMING] Checkpoint " << #name << ": " << _elapsed.count() << " ms from last point" << std::endl; \
                _last_time = _now; \
        } while (0)
  
    #define TIMING_SCOPE(name) \
        struct _TimingScope { \
                const char* _name; \
                std::chrono::high_resolution_clock::time_point _start; \
                _TimingScope(const char* n) : _name(n), _start(std::chrono::high_resolution_clock::now()) {} \
                ~_TimingScope() { \
                        auto _end = std::chrono::high_resolution_clock::now(); \
                        auto _duration = std::chrono::duration_cast<std::chrono::milliseconds>(_end - _start); \
                        std::cout << "[TIMING-SCOPE] " << _name << ": " << _duration.count() << " ms" << std::endl; \
                } \
        } _timing_obj_##name(#name)
  
#else
    #define TIMING_START(name) ((void)0)
    #define TIMING_END(name) ((void)0)
    #define TIMING_POINT(name) ((void)0)
    #define TIMING_SCOPE(name) ((void)0)
#endif
class GMTIProcessor
{
public:
    // 构造函数，初始化配置
    GMTIProcessor() : d_workspace(nullptr), d_workspace_bytes(0), cufft_plan_(0), 
                      cached_Na_(0), cached_Nr_(0), stream_compute_(0)
    {
        // 初始化时创建流
        cudaStreamCreate(&stream_compute_);
        // 初始化指针结构体为 null
        gpu_ptrs_ = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
    }

    ~GMTIProcessor() {
        cleanupCUDAResources();
        if(stream_compute_) cudaStreamDestroy(stream_compute_);
    }

    // 在主流程开始前调用一次
    bool initFFTPlans(const Config& cfg);
    bool initcuFFTPlans(const Config& cfg);

    // Debug: compare CPU vs GPU range compression for a single period
    bool debug_compare_range(const Config &cfg, int periodIdx);

    // 读取配置文件并更新私有成员 cfg_
    bool readXmlParam(const std::string &xmlFile, Config &cfg);

    // 读取 POS 数据
    bool POS_dataread(const std::string &posFile, std::vector<std::vector<double>> &POS_data, int &POS_num);

    // 处理一个周期的数据
    bool processOnePeriod(int periodIdx, const Config &cfg, const std::vector<std::vector<double>> &posRaw, GMTIOutput &out);
    bool processOnePeriod3(int periodIdx, const Config &cfg, const std::vector<std::vector<double>> &posRaw, GMTIOutput &out);
    bool processOnePeriodFusionCache(int periodIdx,
                                     const Config &cfg,
                                     const std::vector<std::vector<double>> &posRaw,
                                     size_t slot,
                                     FusionGroupContext &ctx);

    // 写入结果到文件
    bool writeResult(const std::vector<double> &res, const Config &cfg);

    // ★ 便捷封装：直接做跟踪
    bool trackFromMT(
        const std::vector<std::vector<double>>& MT_frames,
        const Config& cfg,
        double v_max, double sigma_thresh, int min_len, int max_gap,
        std::vector<Track>& tracks_out, bool bypassTracking = false
    );

    // 写入“航迹二进制”文件
    bool writeTracksBinary(const std::vector<Track>& tracks,
                           double utcMid,           // 秒
                           const GMTIOutput::Plane& plane,
                           const Config& cfg);
    
    // 根据 cfg 选择策略，生成 periodList
    bool makePeriodList(const GMTIOutput::Plane& plane,
                        const Config& cfg,
                        std::vector<int>& periodList);

    bool extractPlanePV(const std::vector<std::vector<double>> &POS,
                                    const Config &cfg,
                                    GMTIOutput::Plane &plane);
    bool extractPlanePVFromEcho(const Config &cfg,
                              GMTIOutput::Plane &plane);

    // 生成下一个文件名：dir 末尾可有/无'/'，返回完整路径与编号
    bool nextGMTIFileName(const std::string& dir,
                          std::string& out_path,
                          int& out_index) const;

    bool alignFFTAndDBS(const std::vector<std::complex<double>> &data1,
                        const std::vector<std::complex<double>> &data2,
                        int skip,
                        double fa2,
                        const Config &cfg,
                        std::vector<std::complex<double>> &out1,
                        std::vector<std::complex<double>> &out2);
    bool alignFFTAndDBS_CUDA(const std::vector<std::complex<float>> &data1,
                             const std::vector<std::complex<float>> &data2,
                             int skip,
                             float fa2,
                             const Config &cfg,
                             std::vector<std::complex<float>> &out1,
                             std::vector<std::complex<float>> &out2);

    // 并行处理多个 period（在单 GPU 上使用多个 GMTIProcessor 实例并发运行）
    // 输出 results 与输入 periodList 顺序一致
    bool processPeriodsParallel(const std::vector<int> &periodList,
                                const Config &cfg,
                                const std::vector<std::vector<double>> &posRaw,
                                std::vector<GMTIOutput> &results);
    bool processPeriodsParallelFusion(const std::vector<int> &periodList,
                                      const Config &cfg,
                                      const std::vector<std::vector<double>> &posRaw,
                                      FusionGroupContext &ctx);
    bool processPeriodsParallelFusion(const std::vector<int> &periodList,
                                      const Config &cfg,
                                      const std::vector<std::vector<double>> &posRaw,
                                      FusionGroupContext &ctx,
                                      std::vector<GMTIOutput> &results);

    // Compute dataset-level squint by reading the center period(s) and
    // estimating fd_ctr from center-window data. Returns true on success
    // and writes degrees into out_squint. Implemented in processOnePeriod.cpp
    bool computeDatasetSquintFromCenter(const std::vector<int> &periodList,
                                        const Config &cfg,
                                        const std::vector<std::vector<double>> &posRaw,
                                        double &out_squint);

    // Global squint handling: allow computing/setting a dataset-wide squint angle
    static double estimateSquintAngleDeg(const std::vector<std::complex<double>> &data,
                                         const GMTIOutput::Plane &plane,
                                         const Config &cfg);
    static double estimateSquintAngleDeg(const GMTIOutput::Plane &plane, const Config &cfg, double fd_ctr = std::numeric_limits<double>::quiet_NaN());
    void setGlobalSquint(double deg) { global_squint_deg_ = deg; have_global_squint_ = true; }
    bool hasGlobalSquint() const { return have_global_squint_; }
    double getGlobalSquint() const { return global_squint_deg_; }

private:
    // 私有成员函数
    bool openEchoFiles(const Config &cfg, std::ifstream &fid1, std::ifstream &fid2);
    bool readPulseBlock(const Config &cfg,
                        int beamskip,
                        std::vector<std::complex<double>> &data1,
                        std::vector<std::complex<double>> &data2,
                        std::vector<double> &utc,
                        double& theta_sq);
    bool pulseCompression(std::vector<std::complex<double>> &data1, std::vector<std::complex<double>> &data2, const Config &cfg);
    // GPU implementation of range compression using cuFFT (float kernels)
    bool rangeCompressCUFFT(const std::vector<std::complex<float>> &in,
                            std::vector<std::complex<float>> &rc_out,
                            const Config &cfg);
    // Device-only variant: operate on buffers already uploaded to GPU (gpu_ptrs_.d1 -> gpu_ptrs_.a1)
    bool rangeCompressCUFFT_device(int Lraw, int M, const Config &cfg);
    // Debug helper: download compressed data from device buffer d1 (float complex)
    bool debug_download_d1(std::vector<std::complex<float>> &out, size_t total);
    bool extractUTC(const std::vector<uint8_t> &headers, const Config &cfg, std::vector<double> &utc);
    bool extractPlanePos(const std::vector<double> &t_utc,
                         const std::vector<std::vector<double>> &POS,
                         const Config &cfg,
                         GMTIOutput::Plane &plane);
    bool computeDynamicSupportDomain(
        const std::vector<double> &faAxis,
        double fa2,
        const GMTIOutput::Plane &plane,
        const Config &cfg,
        int &az_center,
        int &az_st,
        int &az_ed,
        double &fd_st,
        double &fd_ed,
        double &BW_az,
        int &rg_st_out,
        int &rg_ed_out);

    // Float版本：用于GPU路径，避免不必要的double转换
    bool readPulseBlockFloat(const Config &cfg,
                             int beamskip,
                             std::vector<std::complex<float>> &data1,
                             std::vector<std::complex<float>> &data2,
                             std::vector<double> &utc,
                             double& theta_sq);

    // 估计距离向相位依赖（二次项），并输出拟合值与原始相位-距离散点
    // F1, F2: 频域（已 FFT+DBS）数据，尺寸 [Na x Nr]，行优先一维存储
    // rg_st..rg_ed, az_st..az_ed：支撑域（0-based，闭区间）
    // rg_len：输出 phi_fit 的长度（通常等于 cfg.rg_len）
    // thre_rg：RANSAC内点阈值（相位弧度），<=0 则取 0.2*pi
    // 输出：
    //   phi_fit        : 长度 rg_len，按 x=1..rg_len 的二次多项式值（与 MATLAB 对齐）
    //   phi_diss_phase : 长度 (rg_ed-rg_st+1)，对应各距离单元的相位 angle(sum(F1.*conj(F2), az))
    //   phi_diss_range : 同长度，存储对应的“距离下标”，与 MATLAB 近似（这里用 1-based: rg_st+1 .. rg_ed+1）
    bool rg_correct(const std::vector<std::complex<double>> &F1,
                    const std::vector<std::complex<double>> &F2,
                    const Config &cfg,
                    double th_rg,
                    std::vector<double> &phi_fit,
                    std::vector<double> &phi_diss_phase,
                    std::vector<int> &phi_diss_range);
    bool rg_correct_CUDA(const Config &cfg,
                         double th_rg,
                         std::vector<float> &phi_fit,
                         std::vector<double> &phi_diss_phase,
                         std::vector<int> &phi_diss_range);

    bool clutter_cancel_38_paper_1(
        const std::vector<double> &y_faAxis,          // len = cfg.pulse_num
        const std::vector<std::complex<double>> &F1f, // size = Na*Nr
        const std::vector<std::complex<double>> &F2f, // size = Na*Nr
        int az_st, int rg_st, int az_ed, int rg_ed,
        const Config &cfg, // ☆ 新增：显式传入 cfg
        std::vector<std::complex<double>> &prosig_38,
        std::array<double, 2> &p_38,
        std::vector<double> &phase_tra_38_cut,
        std::vector<double> &row_fa_cut);

    bool clutter_cancel_38_paper_1_p38(
        const std::vector<double> &y_faAxis,
        const std::vector<std::complex<double>> &F1f,
        const std::vector<std::complex<double>> &F2f,
        int az_st, int rg_st, int az_ed, int rg_ed,
        const Config &cfg,
        std::array<double, 2> &p_38);

    bool clutter_cancel_38_paper_1_p38_cuda(
        const std::vector<double> &y_faAxis,
        int az_st, int rg_st, int az_ed, int rg_ed,
        const Config &cfg,
        std::array<double, 2> &p_38,
        std::vector<double> *phase_tra_38 = nullptr);

    bool clutter_cancel_38_paper_1_cuda(
        const std::vector<double> &y_faAxis,
        int az_st, int rg_st, int az_ed, int rg_ed,
        const Config &cfg,
        std::vector<std::complex<double>> &prosig_38,
        std::array<double, 2> &p_38,
        std::vector<double> &phase_tra_38_cut,
        std::vector<double> &row_fa_cut);

    bool clutter_cancel_38_paper_2(
        const std::vector<double>& y_faAxis,
        const std::vector<std::complex<double>>& F1f,
        const std::vector<std::complex<double>>& F2f,
        int az_st, int rg_st, int az_ed, int rg_ed,
        const Config& cfg,
        std::vector<std::complex<double>>& prosig_38,
        std::array<double,2>& p_38,
        std::vector<double>& phase_tra_38_cut,
        std::vector<double>& row_fa_cut
    );

    bool clutter_38_V2(
        const std::vector<std::complex<double>> &data1,
        const std::vector<std::complex<double>> &data2,
        const std::vector<std::complex<double>> &rg_phi, // len = cfg.rg_len
        int skip_num,
        int az_st, int rg_st,
        int az_ed, int rg_ed,
        const std::vector<double> &y_faAxis, // len = cfg.pulse_num
        double fa2,
        const Config &cfg, // ☆ 新增：显式传入 cfg
        std::vector<std::complex<double>> &CSI_result_38paper,
        std::array<double, 2> &p_38,
        std::vector<double> &phase_38_phase,
        std::vector<double> &phase_38_fa);

    bool dpca_cfar2_fast(const std::vector<std::complex<double>> &GMTI_new,
                    double pf, int c_num, int b_num, const std::string &type,
                    const Config &cfg,
                    std::vector<double> &mydata,
                    std::vector<int> &prow,
                    std::vector<int> &pcol);
    bool dpca_cfar2_fast_cuda(const std::vector<std::complex<float>> &CSI_out,
                              int band_st, int band_ed,
                              float pf, int c_num, int b_num, const std::string &type,
                              const Config &cfg,
                              std::vector<float> &mydata);
    bool target_select(const std::vector<std::complex<double>> &GMTI_dataf_1,
                       const std::vector<std::complex<double>> &GMTI_dataf_2,
                       const std::vector<int> &prow,
                       const std::vector<int> &pcol,
                       const Config &cfg,
                       GMTIOutput::Detect &S_target);
    bool target_select_cuda(const std::vector<int> &prow,
                            const std::vector<int> &pcol,
                            const Config &cfg,
                            GMTIOutput::Detect &S_target);

    bool cluster_filter(const std::vector<double> &mydata,
                        int min_points,
                        const Config &cfg,
                        std::vector<double> &refined_data,
                        std::vector<int> &prow_new,
                        std::vector<int> &pcol_new);

    bool cluster_filter_gap_phase(const std::vector<double> &mydata,
                                             const std::vector<double> &phase_map, // 干涉相位 (rad)
                                             int min_points,
                                             int max_gap,        // 允许的最大列向间断像元数
                                             double maxPhaseStd, // 相位标准差阈值 (rad)
                                             const Config &cfg,
                                             std::vector<double> &refined_data,
                                             std::vector<int> &prow_new,
                                             std::vector<int> &pcol_new,
                                             std::vector<double> &phaseStdList);
    bool cluster_filter_gap_phase_cuda(const std::vector<float> &mydata,
                                       const std::vector<float> &phase_map,
                                       int min_points,
                                       int max_gap,
                                       float maxPhaseStd,
                                       const Config &cfg,
                                       std::vector<float> &refined_data,
                                       std::vector<int> &prow_new,
                                       std::vector<int> &pcol_new,
                                       std::vector<float> &phaseStdList);

    // 已实现：卡尔曼跟踪（前面我给过实现）
    bool kalmanTrack(
        const std::vector<std::vector<std::array<double,2>>>& MT_pos,
        const std::vector<std::vector<double>>& MT_dir,
        const std::vector<std::vector<double>>& MT_range,
        double T, double v_max, double sigma_thresh, int min_len, int max_gap,
        const Config& cfg, std::vector<Track>& tracks_out, bool bypassTracking = false);

    // ★ 新增：从若干帧的 res.MT（每帧 stride=6）构造 MT_pos，并用 utcMid 自动估计 T
    bool buildMTposAndT(
        const std::vector<std::vector<double>>& MT_frames, // 每个元素是该帧的 res.MT 扁平数组
        const Config& cfg,
        std::vector<std::vector<std::array<double,2>>>& MT_pos, // 输出：每帧一组{x,y}
        std::vector<std::vector<double>>& MT_dir,              // 输出：每帧方向
        std::vector<std::vector<double>>& MT_range,            // 输出：每帧距离
        double& T_out                                          // 输出：秒
    );

    // ==== 小型内联线代工具（C++11，避免外部依赖）====
    static inline void mat4_eye(std::array<double,16>& A) {
        for (int i=0;i<16;++i) A[i]=0.0;
        A[0]=A[5]=A[10]=A[15]=1.0;
    }
    static inline void mat4_mul(const std::array<double,16>& A,
                                const std::array<double,16>& B,
                                std::array<double,16>& C) {
        for (int r=0;r<4;++r){
            for (int c=0;c<4;++c){
                double s=0.0;
                for (int k=0;k<4;++k) s += A[r*4+k]*B[k*4+c];
                C[r*4+c]=s;
            }
        }
    }
    static inline void mat4_add_inplace(std::array<double,16>& A,
                                        const std::array<double,16>& B) {
        for (int i=0;i<16;++i) A[i]+=B[i];
    }
    static inline void mat4_mul_vec4(const std::array<double,16>& A,
                                     const std::array<double,4>& x,
                                     std::array<double,4>& y) {
        for (int r=0;r<4;++r){
            double s=0.0; for (int k=0;k<4;++k) s += A[r*4+k]*x[k];
            y[r]=s;
        }
    }
    static inline void vec4_add_inplace(std::array<double,4>& x,
                                        const std::array<double,4>& dx) {
        for (int i=0;i<4;++i) x[i]+=dx[i];
    }

    // 递归建目录（CentOS7 无 <filesystem>）
    static bool mkdir_p(const std::string& dir);

    // —— 工具：角度归一化到 [-180,180)
    static inline double wrap180(double deg) {
        while (deg >= 180.0) deg -= 360.0;
        while (deg <  -180.0) deg += 360.0;
        return deg;
    }

    // 平面近似：以飞机点为原点的 ENU（米）向量（经纬度单位：度）
    // ★ 直接对 EN 坐标计算差分
    static inline void geoDiff_EN_m(double E0, double N0,
                                    double E,  double N,
                                    double& dE_m, double& dN_m) {
        dE_m = E - E0;
        dN_m = N - N0;
    }

    // 默认策略：按 st..ed..skip 取“中间 3 个”（或尽量多）
    static bool makeDefaultRangeList(const Config& cfg, std::vector<int>& periodList);

    // ROI: [lat1,lng1, lat2,lng2, lat3,lng3, lat4,lng4] (单位：度)
    // plane: 需包含 lat_deg, lng_deg（飞机位置，度）与 V_angle_deg（航向角，度；东为0°，逆时针为正）
    // 输出 periodList：3 个波位号（1..51）
    bool makePeriodListROI51(const std::array<double,4>& roi_ll_deg,
                             const GMTIOutput::Plane& plane,
                             std::vector<int>& periodList,
                             double scan_min_deg = -25.0,
                             double scan_max_deg =  25.0);

    Config cfg_; // 私有成员 cfg_ 用来存储配置
    ColFFTTranspose colfft_; // 用于列 FFT 的对象
    std::vector<std::complex<double>> a1_, a2_;  // 复用缓冲，减少分配

    // ===== CUDA 核心资源 (持久化) =====
    // Optional dataset-wide squint estimate (degrees)
    double global_squint_deg_ = 0.0;
    bool have_global_squint_ = false;
    cudaStream_t stream_compute_;       // 专门的计算流，实现异步派发

    // Device memory workspace (预分配，复用)
    void* d_workspace = nullptr;        // CUDA unified buffer: d_d1, d_d2, d_a1, d_a2, d_t1, d_t2
    size_t d_workspace_bytes = 0;       // 已分配的字节数
    
    // 显存池内的偏移指针 (仅作为引用，不单独 malloc/free)
    // 这样做是为了在 align, FFT, DBS 之间传递数据时，不需要在函数里重新计算偏移
    struct {
        void* d1; // 原始通道1
        void* d2; // 原始通道2
        void* a1; // 对齐/FFT后通道1
        void* a2; // 对齐/FFT后通道2
        void* t1; // DBS临时缓冲1
        void* t2; // DBS临时缓冲2
        void* csi; // 对消结果 (CSI) 缓冲
        void* d_rg_sums; // 距离向相位校正的中间结果（每距离单元的 az 维和）
    } gpu_ptrs_;

    float* d_phi_fit_ = nullptr; // 距离向相位校正的拟合值（长度 cfg.rg_len）

    std::vector<float> cached_phi_fit_;

    // cuFFT 计划缓存
    // 将原有的 int 改为标准的 cufftHandle 语义更清晰
    cufftHandle cufft_plan_ = 0; 
    int cached_Na_ = 0, cached_Nr_ = 0;
    
    // --- 第一阶段：资源检查与数据上传 ---
    // 负责检查 d_workspace 尺寸是否够，并异步把数据丢进 GPU
    bool cuda_upload_async(const std::vector<std::complex<float>> &data1,
                           const std::vector<std::complex<float>> &data2,
                           size_t Na, size_t Nr);

    // --- 第二阶段：通道对齐 Kernel ---
    bool cuda_stage_align_async(int skip, size_t Na, size_t Nr);

    // --- 第三阶段：多维 FFT ---
    bool cuda_stage_fft_async(size_t Na, size_t Nr);

    // --- 第四阶段：DBS 中心化 ---
    bool cuda_stage_dbs_async(float fa2, size_t Na, size_t Nr);

    bool exportDbsCacheAfterRecenter(const Config& cfg,
                                     const FusionBeamMeta& beamMeta,
                                     size_t slot,
                                     size_t Na,
                                     size_t Nr,
                                     RDData& rd,
                                     MetaPack& meta);

    // --- 第五阶段：数据下载与同步 ---
    bool cuda_download_sync(std::vector<std::complex<float>> &out1,
                            std::vector<std::complex<float>> &out2,
                            size_t total);

    bool cuda_stage_rg_sum_async(const Config& cfg);
    bool cuda_download_rg_sums_sync(std::vector<std::complex<float>>& h_sums, int Nr);
    bool cuda_apply_rg_correction_async(const std::vector<float>& phi_fit, int Na, int Nr);
    bool cuda_download_phase_map(std::vector<float> &phase_map, size_t total);
    bool cuda_download_csi_sync(std::vector<std::complex<float>> &out, size_t total);

    // ===== 清理函数（在析构函数或明确调用时释放 GPU 内存）=====
    void cleanupCUDAResources();
};

#endif // FUNCTIONS_H
