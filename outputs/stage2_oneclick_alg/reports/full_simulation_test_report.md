# GMTI 阶段二一键仿真测试报告

## 1. 测试结论

- 输出目录：`outputs/stage2_oneclick_alg`
- 回波数据：`outputs/stage2_oneclick_alg/data/stage2_statistical_newprotocol.bin`，大小 1.40 GiB
- 仿真 PRT 数：7930
- 散射点总数：1070，其中面杂波 1000，强散射点 20，线状散射点 50
- 合作目标 truth：pulse 7930 行，beam summary 61 行
- NaN/Inf：has_nan=false，has_inf=false
- 算法运行：log 已生成，normal_exit=True，错误 1，警告 0

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
- scene_mode：`full`
- target_enabled：`true`

## 3. 回波数据生成

- 数据格式：新协议 PRT，256 字节头 + `11820 * 16` 字节双通道 float32 IQ。
- `scatterer_echoes`：157418
- `scatterer_samples`：929388612
- `max_abs_component`：64.8316
- `mean_noise_power_per_complex_sample`：0.00999972
- 抽样 DDC RMS：3.07148，抽样平均复功率：9.43399，抽样最大分量幅度：11.9355

## 4. 合作目标注入

- target_config_path：`targets.json`
- target_snr_db：25.0
- amplitude_mode：`snr_db`
- target_pulses_injected：130
- target_samples_injected：927789
- 可见 beam 数：1，已注入 beam 数：1，示例 beam：30
- truth 平均径向速度：-49.8748 m/s
- truth 平均 DDC range_sample：4682.66

## 5. 成像效果

- 算法图像/文本产品数量：PNG 5，TXT 1
- 示例文件：GMTI02.bin, GMTI02_track.bin, GMTI02.png, GMTI02_angle_labels.png, GMTI02_beam19_source1_before.png, GMTI02_beam44_source1_after.png, GMTI02_source1_before_after.png, GMTI02.txt

## 6. 对消效果

- 算法日志存在，但未解析到明确对消指标。建议后续在算法中输出对消前后功率或 SCR。

## 7. 检测效果

- 检测结果文件数量：1
- 当前脚本尚未解析 `GMTIxx.bin` 二进制字段，检测率/虚警率需下一步接入结果解析。

## 8. 最小可检测速度评估

- 当前没有自动 SNR/速度梯度实验结果，因此这里给出体制尺度估计，不能替代真实检测门限。
- 波长：0.018737 m
- 单波位驻留时间：0.1 s
- 多普勒 bin 宽度：10 Hz
- 对应速度分辨尺度：0.0936851 m/s
- 第一盲速尺度：12.1791 m/s
- 建议后续自动生成 0, 1, 2, 3, 5, 10 m/s 目标速度组，并以 90% 检测率定义最小可检测速度。

## 9. 定位精度

- 检测结果存在，但当前未解析检测点坐标，尚未与 target truth 匹配。
- truth 文件已经提供 `range_m`、`range_sample_float`、`target_azimuth_deg`、`theta_true_deg`，可作为后续匹配基准。

## 10. 航迹关联

- 航迹文件数量：1
- 当前脚本尚未解析 28 字节 track 包字段，TrackID 稳定性需要下一步接入解析。

## 11. 输出文件索引

- `outputs/stage2_oneclick_alg/config/stage2_config.json`：存在
- `outputs/stage2_oneclick_alg/config/temp_config_stage2_newsystem.xml`：存在
- `outputs/stage2_oneclick_alg/data/stage2_statistical_newprotocol.bin`：存在
- `outputs/stage2_oneclick_alg/truth/scene_truth.csv`：存在
- `outputs/stage2_oneclick_alg/truth/target_truth_pulse.csv`：存在
- `outputs/stage2_oneclick_alg/truth/target_truth_beam_summary.csv`：存在
- `outputs/stage2_oneclick_alg/reports/stage2_simulation_report.md`：存在
- `outputs/stage2_oneclick_alg/logs/simulate_stage2.log`：存在
- `outputs/stage2_oneclick_alg/logs/algorithm_run.log`：存在

## 12. 下一步建议

1. 用单点场景跑 `GMTI_core`，标定 LFM 符号、rect 定义和脉压峰固定偏移。
2. 增加 `GMTIxx.bin` 和 `GMTIxx_track.bin` 解析器，自动计算检测率、虚警率、定位误差和航迹稳定性。
3. 增加速度/SNR 梯度批量生成，给出真实最小可检测速度曲线。
4. 对 30000 点场景增加 beam/range 索引和并行化，降低生成时间。
