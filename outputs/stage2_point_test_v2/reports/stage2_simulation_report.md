# 阶段二统计杂波纯仿真报告

## 新系统参数

- fc: 16 GHz
- Br: 50 MHz
- fs: 60 MHz
- Tr: 130 us
- PRF: 1300 Hz
- DDC_len: 11820
- FFT_len: 12288
- pc_crop_start: 3864
- pc_crop_len: 4096
- beam_count: 61
- pulse_num: 130
- scan_min_deg: -60
- scan_step_deg: 2
- beam_width_deg: 2.28

- platform_height_m: 6000
- platform_speed_mps: 60

## 场景配置

- scene_mode: point_target_only
- random_seed: 202606
- area_scatterers: 0
- strong_scatterers: 0
- line_scatterers: 0
- single_or_other_scatterers: 1
- thermal_noise_enabled: true
- thermal_noise_power: 0.01
- target_enabled: false

## 生成统计

- packets_written: 7930
- scatterer_echoes: 390
- scatterer_samples: 2788119
- target_pulses_injected: 0
- target_samples_injected: 0
- max_abs_component: 1
- mean_noise_power_per_complex_sample: 0
- has_nan: false
- has_inf: false
- generation_time_ms: 1453.81

## 当前局限

- 第一版平台模型为 ideal_straight，真实 POS/姿态尚未接入。
- LFM rect 延迟采用 `0 <= dt < Tr`，仍需用算法脉压结果标定峰值固定偏移。
- 检测评价和航迹评价文件已预留，当前尚未自动解析 GMTI 检测输出。
