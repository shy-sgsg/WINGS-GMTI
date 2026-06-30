# GMTI 阶段二一键仿真测试报告

## 1. 测试结论

- 输出目录：`outputs/stage2_oneclick_point`
- 回波数据：`outputs/stage2_oneclick_point/data/stage2_statistical_newprotocol.bin`，大小 1.40 GiB
- 仿真 PRT 数：7930
- 散射点总数：1，其中面杂波 0，强散射点 0，线状散射点 0
- 合作目标 truth：pulse 0 行，beam summary 0 行
- NaN/Inf：has_nan=false，has_inf=false
- 算法运行：未执行或未发现 `logs/algorithm_run.log`，成像/对消/检测/航迹指标仅报告数据侧准备情况。

## 2. 仿真参数设置

| 参数 | 值 |
|---|---:|
| fc_ghz | 16.0 |
| bandwidth_mhz | 50.0 |
| fs_mhz | 60.0 |
| pulse_width_us | 130.0 |
| prf_hz | 1300.0 |
| sample_delay_us | 488.0 |
| sample_window_us | 197.0 |
| ddc_len | 11820 |
| fft_len | 12288 |
| pc_crop_start | 3864 |
| pc_crop_len | 4096 |
| scan_min_deg | -60.0 |
| scan_step_deg | 2.0 |
| beam_count | 61 |
| beam_width_deg | 2.28 |
| pulse_num | 130 |
| platform_height_m | 6000.0 |
| platform_speed_mps | 60.0 |
| d_chan_m | 0.17 |

- random_seed：`202606`
- scene_mode：`point_target_only`
- target_enabled：`false`

## 3. 回波数据生成

- 数据格式：新协议 PRT，256 字节头 + `11820 * 16` 字节双通道 float32 IQ。
- `scatterer_echoes`：390
- `scatterer_samples`：2788119
- `max_abs_component`：1
- `mean_noise_power_per_complex_sample`：0
- 抽样 DDC RMS：0，抽样平均复功率：0，抽样最大分量幅度：0

## 4. 合作目标注入

- target_config_path：`targets.json`
- target_snr_db：25.0
- amplitude_mode：`snr_db`
- target_pulses_injected：0
- target_samples_injected：0
- 可见 beam 数：0，已注入 beam 数：0，示例 beam：无

## 5. 成像效果

- 未发现算法成像产品。当前报告只能确认 DDC 回波和 truth 已生成；请用生成的 XML 跑 `GMTI_core` 后再刷新报告。

## 6. 对消效果

- 未运行算法，无法评价对消效果。数据侧已生成统计杂波、强散射点、线状散射体和噪声。

## 7. 检测效果

- 未发现 `GMTIxx.bin` 检测结果，检测率、虚警数和漏检数暂不可评估。

## 8. 最小可检测速度评估

- 当前没有自动 SNR/速度梯度实验结果，因此这里给出体制尺度估计，不能替代真实检测门限。
- 波长：0.018737 m
- 单波位驻留时间：0.1 s
- 多普勒 bin 宽度：10 Hz
- 对应速度分辨尺度：0.0936851 m/s
- 第一盲速尺度：12.1791 m/s
- 建议后续自动生成 0, 1, 2, 3, 5, 10 m/s 目标速度组，并以 90% 检测率定义最小可检测速度。

## 9. 定位精度

- 未发现检测结果，无法计算距离误差、方位误差和位置误差。
- truth 文件已经提供 `range_m`、`range_sample_float`、`target_azimuth_deg`、`theta_true_deg`，可作为后续匹配基准。

## 10. 航迹关联

- 未发现 `GMTIxx_track.bin`，confirmed_track_count、id_switch_count、track_break_count 暂不可评估。

## 11. 输出文件索引

- `outputs/stage2_oneclick_point/config/stage2_config.json`：存在
- `outputs/stage2_oneclick_point/config/temp_config_stage2_newsystem.xml`：存在
- `outputs/stage2_oneclick_point/data/stage2_statistical_newprotocol.bin`：存在
- `outputs/stage2_oneclick_point/truth/scene_truth.csv`：存在
- `outputs/stage2_oneclick_point/truth/target_truth_pulse.csv`：未生成
- `outputs/stage2_oneclick_point/truth/target_truth_beam_summary.csv`：未生成
- `outputs/stage2_oneclick_point/reports/stage2_simulation_report.md`：存在
- `outputs/stage2_oneclick_point/logs/simulate_stage2.log`：存在
- `outputs/stage2_oneclick_point/logs/algorithm_run.log`：未生成

## 12. 下一步建议

1. 用单点场景跑 `GMTI_core`，标定 LFM 符号、rect 定义和脉压峰固定偏移。
2. 增加 `GMTIxx.bin` 和 `GMTIxx_track.bin` 解析器，自动计算检测率、虚警率、定位误差和航迹稳定性。
3. 增加速度/SNR 梯度批量生成，给出真实最小可检测速度曲线。
4. 对 30000 点场景增加 beam/range 索引和并行化，降低生成时间。
