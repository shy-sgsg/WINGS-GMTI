#include "stage2_validator.h"

#include "../target_injection/target_common.h"

#include <fstream>

namespace gmti {
namespace stage2 {

bool writeStage2Report(const std::string &path,
                       const Stage2Config &cfg,
                       const Stage2RunOptions &opt,
                       const ScattererList &scatterers,
                       const Stage2Stats &stats,
                       double elapsed_sec,
                       std::string &err)
{
    std::ofstream out(path.c_str());
    if (!out) {
        err = "failed to write stage2 report";
        return false;
    }
    size_t area = 0, strong = 0, line = 0, single = 0;
    for (size_t i = 0; i < scatterers.size(); ++i) {
        if (scatterers[i].type == "area") ++area;
        else if (scatterers[i].type == "strong") ++strong;
        else if (scatterers[i].type == "line") ++line;
        else ++single;
    }
    out << "# 阶段二统计杂波纯仿真报告\n\n";
    out << "## 新系统参数\n\n";
    out << "- fc: " << cfg.radar.fc_hz / 1.0e9 << " GHz\n";
    out << "- Br: " << cfg.radar.br_hz / 1.0e6 << " MHz\n";
    out << "- fs: " << cfg.radar.fs_hz / 1.0e6 << " MHz\n";
    out << "- Tr: " << cfg.radar.tr_sec * 1.0e6 << " us\n";
    out << "- PRF: " << cfg.radar.prf_hz << " Hz\n";
    out << "- DDC_len: " << cfg.radar.pulse_len << "\n";
    out << "- FFT_len: " << cfg.radar.range_fft_len << "\n";
    out << "- pc_crop_start: " << cfg.radar.range_crop_start << "\n";
    out << "- pc_crop_len: " << cfg.radar.range_crop_len << "\n";
    out << "- beam_count: " << cfg.radar.beam_count << "\n";
    out << "- pulse_num: " << cfg.radar.pulse_num << "\n";
    out << "- scan_min_deg: " << cfg.radar.scan_min_deg << "\n";
    out << "- scan_step_deg: " << cfg.radar.scan_step_deg << "\n";
    out << "- beam_width_deg: " << cfg.radar.beam_width_deg << "\n\n";
    out << "- platform_height_m: " << cfg.platform_height_m << "\n";
    out << "- platform_speed_mps: " << cfg.platform_speed_mps << "\n\n";
    out << "## 场景配置\n\n";
    out << "- scene_mode: " << opt.scene_mode << "\n";
    out << "- random_seed: " << cfg.sim.random_seed << "\n";
    out << "- area_scatterers: " << area << "\n";
    out << "- strong_scatterers: " << strong << "\n";
    out << "- line_scatterers: " << line << "\n";
    out << "- single_or_other_scatterers: " << single << "\n";
    out << "- thermal_noise_enabled: " << (cfg.scene.noise.enabled ? "true" : "false") << "\n";
    out << "- thermal_noise_power: " << cfg.scene.noise.noise_power << "\n";
    out << "- target_enabled: " << (opt.target_enabled && cfg.target.enabled ? "true" : "false") << "\n\n";
    out << "## 生成统计\n\n";
    out << "- packets_written: " << stats.packets_written << "\n";
    out << "- scatterer_echoes: " << stats.scatterer_echoes << "\n";
    out << "- scatterer_samples: " << stats.scatterer_samples << "\n";
    out << "- target_pulses_injected: " << stats.target_pulses_injected << "\n";
    out << "- target_samples_injected: " << stats.target_samples_injected << "\n";
    out << "- max_abs_component: " << stats.max_abs_component << "\n";
    out << "- mean_noise_power_per_complex_sample: "
        << (stats.noise_samples ? stats.sum_noise_power / static_cast<double>(stats.noise_samples) : 0.0) << "\n";
    out << "- has_nan: " << (stats.has_nan ? "true" : "false") << "\n";
    out << "- has_inf: " << (stats.has_inf ? "true" : "false") << "\n";
    out << "- generation_time_ms: " << elapsed_sec * 1000.0 << "\n\n";
    out << "## 当前局限\n\n";
    out << "- 第一版平台模型为 ideal_straight，真实 POS/姿态尚未接入。\n";
    out << "- LFM rect 延迟采用 `0 <= dt < Tr`，仍需用算法脉压结果标定峰值固定偏移。\n";
    out << "- 检测评价和航迹评价文件已预留，当前尚未自动解析 GMTI 检测输出。\n";
    return true;
}

} // namespace stage2
} // namespace gmti
