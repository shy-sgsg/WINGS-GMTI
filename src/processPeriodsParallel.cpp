// processPeriodsParallel.cpp
// 并行派发多个 period 到多个独立的 GMTIProcessor 实例（单 GPU，多 stream/多 workspace）

#include "GMTIProcessor.hpp"
#include "dbs/DbsFusion.hpp"
#include <thread>
#include <atomic>
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

bool GMTIProcessor::processPeriodsParallel(const std::vector<int> &periodList,
                                           const Config &cfg,
                                           const std::vector<std::vector<double>> &posRaw,
                                           std::vector<GMTIOutput> &results)
{
    if (periodList.empty()) return true;

    // 查询可用显存
    size_t freeBytes = 0, totalBytes = 0;
    cudaError_t cerr = cudaMemGetInfo(&freeBytes, &totalBytes);
    if (cerr != cudaSuccess) {
        std::cerr << "cudaMemGetInfo failed: " << cudaGetErrorString(cerr) << std::endl;
        return false;
    }
    std::cout << "[parallel] GPU free/total MB: " << (freeBytes/1024/1024) << " / " << (totalBytes/1024/1024) << std::endl;

    // 每个实例所需的显存基于当前对象的 d_workspace_bytes（在 initcuFFTPlans() 后设置）
    size_t per_instance_bytes = this->d_workspace_bytes;
    if (per_instance_bytes == 0) {
        // 保守回退值（约200MB），在未初始化的情况下使用
        per_instance_bytes = 100ull * 1024ull * 1024ull;
    }

    // std::cout << "[parallel] per_instance_bytes MB: " << (per_instance_bytes/1024/1024) << std::endl;

    // 额外开销（CFAR / phase map / 临时 device malloc），按单 period 的实际峰值做保守预算
    const size_t overhead = 40ull * 1024ull * 1024ull;
    size_t per_needed = per_instance_bytes + overhead;

    size_t max_by_mem = std::max<size_t>(1, static_cast<int>(freeBytes / per_needed));
    size_t requested = periodList.size();
    // 限制并发实例数，避免过多小片段导致上下文切换，设置合理上限
    const size_t HARD_CAP = 8;
    size_t instances = std::min<size_t>({requested, max_by_mem, HARD_CAP});

    std::cout << "[parallel] requested=" << requested << " max_by_mem=" << max_by_mem
              << " hard_cap=" << HARD_CAP << " per_needed_mb=" << (per_needed/1024/1024)
              << " -> instances=" << instances << std::endl;

    if (instances == 0) instances = 1;

        // --- Compute dataset-wide squint once using DBS-like center-wavepos method.
        //     If disabled, use the XML value directly. If enabled but estimation fails,
        //     fall back to XML cfg.squint_angle; final fallback 0.
        double final_squint = 0.0;
        if (!cfg.estimate_error_angle) {
            final_squint = std::isfinite(cfg.squint_angle) ? cfg.squint_angle : 0.0;
            std::cout << "[parallel] estimate_error_angle disabled, using XML squint="
                      << final_squint << " deg" << std::endl;
        } else {
            double computed = 0.0;
            bool ok = this->computeDatasetSquintFromCenter(periodList, cfg, posRaw, computed);
            if (!ok || !std::isfinite(computed) || std::abs(computed) > 30.0) {
                // try XML-provided value
                if (std::isfinite(cfg.squint_angle)) {
                    final_squint = cfg.squint_angle;
                    std::cout << "[parallel] computeDatasetSquintFromCenter failed, using XML squint="
                              << final_squint << " deg" << std::endl;
                } else {
                    final_squint = 0.0;
                    std::cout << "[parallel] computeDatasetSquintFromCenter failed, no XML fallback, using 0 deg" << std::endl;
                }
            } else {
                final_squint = computed;
                std::cout << "[parallel] estimated squint angle = " << computed << " deg" << std::endl;
            }
        }
        std::cout << "[parallel] final squint angle set to = " << final_squint << " deg" << std::endl;

    // 创建 worker 实例
    // 如果用户在 XML 中关闭了波位并行，则改为顺序逐波位计算
    if (!cfg.wavepos_parallel) {
        std::cout << "[parallel] wavepos_parallel disabled: running sequentially" << std::endl;
        // 需要为当前实例初始化 FFT/cuFFT plans
        if (!this->initFFTPlans(cfg)) {
            std::cerr << "initFFTPlans failed for sequential processing" << std::endl;
            return false;
        }
        if (!this->initcuFFTPlans(cfg)) {
            std::cerr << "initcuFFTPlans failed for sequential processing" << std::endl;
            return false;
        }
        this->setGlobalSquint(final_squint);

        results.clear();
        results.resize(periodList.size());
        bool all_ok = true;
        for (size_t i = 0; i < periodList.size(); ++i) {
            int per = periodList[i];
            GMTIOutput out;
            bool s = this->processOnePeriod(per, cfg, posRaw, out);
            results[i] = std::move(out);
            all_ok = all_ok && s;
            std::cout << "[sequential] period " << per << " -> " << (s?"OK":"FAIL") << std::endl;
        }
        return all_ok;
    }

    std::vector<std::unique_ptr<GMTIProcessor>> workers;
    workers.reserve(instances);
    for (size_t i = 0; i < instances; ++i) {
        workers.emplace_back(new GMTIProcessor());
        // 每个 worker 需要初始化 FFT/cuFFT plans 与显存
        std::cout << "[parallel] Initializing worker " << i << std::endl;
        if (!workers.back()->initFFTPlans(cfg)) {
            std::cerr << "Worker initFFTPlans failed for instance " << i << std::endl;
            return false;
        }
        if (!workers.back()->initcuFFTPlans(cfg)) {
            std::cerr << "Worker initcuFFTPlans failed for instance " << i << std::endl;
            return false;
        }
        // Propagate dataset-wide squint to each worker
        workers.back()->setGlobalSquint(final_squint);
        // std::cout << "[parallel] Worker " << i << " initialized (workspace bytes=" << workers.back()->d_workspace_bytes/1024/1024 << " MB)" << std::endl;
    }

    // 准备输出容器
    const size_t N = periodList.size();
    results.clear();
    results.resize(N);
    std::vector<char> ok(N, 0);

    std::atomic_size_t next_idx(0);

    // 线程池：每个 worker 对应一个线程，抢占式分配 period
    std::vector<std::thread> threads;
    threads.reserve(instances);

    for (size_t w = 0; w < instances; ++w) {
        threads.emplace_back([w, &workers, &periodList, &cfg, &posRaw, &results, &ok, &next_idx]() {
            GMTIProcessor *worker = workers[w].get();
            while (true) {
                size_t i = next_idx.fetch_add(1);
                if (i >= periodList.size()) break;
                int per = periodList[i];
                GMTIOutput out;
                bool s = worker->processOnePeriod(per, cfg, posRaw, out);
                results[i] = std::move(out);
                ok[i] = s ? 1 : 0;
                std::cout << "[parallel] period " << per << " done by worker " << w << " -> " << (s?"OK":"FAIL") << std::endl;
            }
        });
    }

    for (auto &t : threads) if (t.joinable()) t.join();

    // 汇总结果一致性检查
    bool all_ok = std::all_of(ok.begin(), ok.end(), [](char v){return v != 0;});
    if (!all_ok) std::cerr << "Warning: some periods failed in parallel processing." << std::endl;

    return all_ok;
}

bool GMTIProcessor::processPeriodsParallelFusion(const std::vector<int> &periodList,
                                                 const Config &cfg,
                                                 const std::vector<std::vector<double>> &posRaw,
                                                 FusionGroupContext &ctx)
{
    if (periodList.empty()) return true;

    ctx.reset(periodList);

    size_t freeBytes = 0, totalBytes = 0;
    cudaError_t cerr = cudaMemGetInfo(&freeBytes, &totalBytes);
    if (cerr != cudaSuccess) {
        std::cerr << "cudaMemGetInfo failed: " << cudaGetErrorString(cerr) << std::endl;
        return false;
    }

    const size_t requested = periodList.size();
    const size_t HARD_CAP = 8;
    const size_t instances = std::max<size_t>(1, std::min(requested, HARD_CAP));
    const size_t pulseDec = static_cast<size_t>(std::max(1, cfg.pulse_dec));
    const int procPulseNum = effectivePulseNum(cfg);
    const size_t rdRows = (static_cast<size_t>(std::max(0, procPulseNum)) + pulseDec - 1) / pulseDec;
    const size_t rdCols = static_cast<size_t>(std::max(0, cfg.rg_len));
    const double rdCacheGiB = static_cast<double>(periodList.size()) * static_cast<double>(rdRows) *
                              static_cast<double>(rdCols) * static_cast<double>(sizeof(float)) /
                              (1024.0 * 1024.0 * 1024.0);
    std::cout << "[fusion] group start: periods=" << periodList.size()
              << " workers=" << instances
              << " cuda_free=" << (static_cast<double>(freeBytes) / (1024.0 * 1024.0 * 1024.0))
              << "GiB/" << (static_cast<double>(totalBytes) / (1024.0 * 1024.0 * 1024.0))
              << "GiB pulse_num=" << cfg.pulse_num
              << " read_pulse_num=" << cfg.read_pulse_num
              << " process_pulse_num=" << procPulseNum
              << " pulse_dec=" << cfg.pulse_dec
              << " range_input_len=" << cfg.pulse_len
              << " range_fft_len=" << effectiveRangeFftLen(cfg)
              << " range_crop_start=" << cfg.range_crop_start
              << " range_output_len=" << effectiveRangeCompressLen(cfg)
              << " rg_len=" << cfg.rg_len
              << " estimated_dbs_amp_cache~" << rdCacheGiB << "GiB" << std::endl;

    std::vector<std::unique_ptr<GMTIProcessor>> workers;
    workers.reserve(instances);
    for (size_t i = 0; i < instances; ++i) {
        workers.emplace_back(new GMTIProcessor());
        if (!workers.back()->initFFTPlans(cfg)) {
            return false;
        }
        if (!workers.back()->initcuFFTPlans(cfg)) {
            return false;
        }
    }

    std::atomic_size_t next_idx(0);
    std::vector<char> ok(periodList.size(), 0);
    std::vector<std::thread> threads;
    threads.reserve(instances);

    for (size_t w = 0; w < instances; ++w) {
        threads.emplace_back([w, &workers, &periodList, &cfg, &posRaw, &ctx, &ok, &next_idx]() {
            GMTIProcessor *worker = workers[w].get();
            while (true) {
                const size_t slot = next_idx.fetch_add(1);
                if (slot >= periodList.size()) break;
                const int per = periodList[slot];
                const bool s = worker->processOnePeriodFusionCache(per, cfg, posRaw, slot, ctx);
                ok[slot] = s ? 1 : 0;
                // std::cout << "[fusion] period " << per << " slot " << slot
                //           << " done by worker " << w << " -> " << (s ? "OK" : "FAIL") << std::endl;
            }
        });
    }

    for (auto &t : threads) if (t.joinable()) t.join();

    const bool all_ok = std::all_of(ok.begin(), ok.end(), [](char v){ return v != 0; });
    if (!all_ok) {
        return false;
    }

    std::vector<int> activePeriods;
    std::vector<FusionBeamMeta> activeBeamMeta;
    activePeriods.reserve(periodList.size());
    activeBeamMeta.reserve(ctx.beam_meta.size());
    for (size_t slot = 0; slot < periodList.size(); ++slot) {
        if (!fusionSlotHasSignal(ctx, slot)) {
            continue;
        }
        activePeriods.push_back(periodList[slot]);
        activeBeamMeta.push_back(ctx.beam_meta[slot]);
    }
    if (activePeriods.empty()) {
        std::cerr << "[fusion] no active DBS beam after zero-fill filtering" << std::endl;
        return false;
    }
    std::cout << "[fusion] active DBS beams for bias:";
    for (int p : activePeriods) {
        std::cout << ' ' << p;
    }
    std::cout << std::endl;

    double biasDeg = 0.0;
    if (!cfg.estimate_error_angle) {
        biasDeg = std::isfinite(cfg.squint_angle) ? cfg.squint_angle : 0.0;
        std::cout << "[fusion] estimate_error_angle disabled, using XML beam_pointing_bias="
                  << biasDeg << " deg" << std::endl;
    } else {
        if (!estimateBeamPointingBiasByCenterBeams(activePeriods, activeBeamMeta, biasDeg)) {
            std::cerr << "[fusion] estimateBeamPointingBiasByCenterBeams failed" << std::endl;
            return false;
        }
    }
    if (!applyBeamPointingBiasToFusionContext(biasDeg, ctx)) {
        std::cerr << "[fusion] applyBeamPointingBiasToFusionContext failed" << std::endl;
        return false;
    }

    std::cout << "[fusion] beam_pointing_bias=" << biasDeg << " deg" << std::endl;
    return true;
}

bool GMTIProcessor::processPeriodsParallelFusion(const std::vector<int> &periodList,
                                                 const Config &cfg,
                                                 const std::vector<std::vector<double>> &posRaw,
                                                 FusionGroupContext &ctx,
                                                 std::vector<GMTIOutput> &results)
{
    if (!processPeriodsParallelFusion(periodList, cfg, posRaw, ctx)) {
        return false;
    }
    if (!relocateFusionDetections(ctx, cfg, results)) {
        std::cerr << "[fusion] relocateFusionDetections failed" << std::endl;
        return false;
    }
    return true;
}
