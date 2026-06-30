# 合作目标 DDC 域注入报告

## 输入输出

- input_config: `outputs/stage1/config/temp_config_stage1_newsystem.xml`
- input_data_dir: `outputs/stage1/data`
- output_dir: `outputs/stage1_target`
- output_data_file: `outputs/stage1_target/data/stage1_background_with_target.bin`

## 新系统参数

- fc: 16 GHz
- Br: 50 MHz
- fs: 60 MHz
- Tr: 130 us
- PRF: 1300 Hz
- DDC_len: 11820
- beam_count: 61
- pulse_num: 130
- scan_min_deg: -61
- scan_step_deg: 2
- d_chan_m: 0.17

## 目标与注入模式

- target_id: 1
- target_name: strong_cooperative_target
- init_mode: range_azimuth
- p0_m: (84788, 0, 0)
- v_mps: (10, 0, 0)
- visibility_mode: hard_gate
- amplitude_mode: snr_db
- target_snr_db: 30
- direct_amplitude: 1
- chirp_phase_sign: 1
- carrier_phase_sign: -1
- channel_phase_mode: baseline_approx

## 生成统计

- packets_read: 1
- packets_written: 1
- pulses_injected: 1
- samples_injected: 7136
- max_amplitude: 153302
- has_nan: false
- has_inf: false

## 当前说明

- 合作目标注入位置是 DDC 原始快时间域 `range_sample_float`，不是脉压后 4096 点 bin。
- 第一版 LFM 延迟约定采用 `0 <= t_fast - tau_rel < Tr`，后续需用单点脉压测试标定固定偏移和符号。
- 第一版平台模型默认为 `ideal_platform`，真实 POS/导航接入留给下一步。
