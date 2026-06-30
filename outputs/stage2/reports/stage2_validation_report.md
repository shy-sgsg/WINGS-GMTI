# 阶段二统计杂波纯仿真验证报告

## 编译

已执行：

```bash
cmake -S . -B build
cmake --build build --target simulate_stage2_statistical -j2
```

结果：编译通过。

## 测试 1：单点静止散射点

命令：

```bash
./build/simulate_stage2_statistical --stage2-config stage2_config.json --output-dir outputs/stage2_point_test --scene-mode point_target_only --target-enabled false --single-scatterer-range 85000 --single-scatterer-azimuth 0 --single-scatterer-amplitude 1.0 --period-count 1 --validate true
```

结果：

- packets: 7930
- scatterers: 1
- scatterer_echoes: 390
- target_pulses: 0
- elapsed_ms: 约 1611

结论：单点静止散射点 DDC 正演可生成完整单周期数据。脉压峰位置尚需接入算法做标定。

## 测试 2：纯噪声

命令：

```bash
./build/simulate_stage2_statistical --stage2-config stage2_config.json --output-dir outputs/stage2_noise_only --scene-mode noise_only --target-enabled false --period-count 1 --validate true
```

结果：

- packets: 7930
- scatterers: 0
- elapsed_ms: 约 6348
- 报告中 mean_noise_power_per_complex_sample 接近配置值 0.01。

结论：热噪声功率控制正常，无 NaN/Inf。

## 测试 3：小规模 full 场景 + 合作目标

命令：

```bash
./build/simulate_stage2_statistical --stage2-config stage2_config.json --target-config targets.json --output-dir outputs/stage2 --scene-mode full --target-enabled true --period-count 1 --validate true
```

结果：

- packets_written: 7930
- scatterers: 1070
- scatterer_echoes: 157418
- scatterer_samples: 929388612
- target_pulses_injected: 130
- target_samples_injected: 927789
- max_abs_component: 64.8316
- mean_noise_power_per_complex_sample: 0.00999972
- has_nan: false
- has_inf: false
- generation_time_ms: 34479

输出：

- `outputs/stage2/data/stage2_statistical_newprotocol.bin`
- `outputs/stage2/truth/scene_truth.csv`
- `outputs/stage2/truth/area_clutter_scatterers.csv`
- `outputs/stage2/truth/strong_scatterers.csv`
- `outputs/stage2/truth/line_scatterers.csv`
- `outputs/stage2/truth/target_truth_pulse.csv`
- `outputs/stage2/truth/target_truth_beam_summary.csv`
- `outputs/stage2/reports/stage2_simulation_report.md`

## 尚未完成

- 尚未调用 CUDA GMTI 算法做单点脉压峰标定。
- 检测评价和航迹评价 CSV 目前是占位输出，尚未自动解析 GMTI 检测结果。
- 30000 点全量杂波正演尚未优化，当前建议先用小规模参数调试。
- `use_existing_nav` 和真实姿态/地理坐标转换尚未接入。

