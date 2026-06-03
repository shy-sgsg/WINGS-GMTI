# DBS/GMTI 融合修改记录

## 修改目标

将 `DBS_CUDA` 中的 DBS 成像流程融入 GMTI 检测定位工程，以 GMTI 工程为主体，统一使用新协议回波文件。融合模式下，单波位阶段只做数据读取、脉压/抽取、GMTI 检测和 DBS 中间数据缓存；所有波位完成后，再估计波束指向偏差、统一解多普勒模糊，并重新定位目标。

## 关键原则

- 新协议为主：双通道回波、UTC、POS 均从同一个新协议回波文件读取。
- GMTI 主流程为主体：DBS 作为融合成像输出模块挂接。
- 单波位 worker 内不做最终解模糊定位。
- 每个波位只保存未解模糊多普勒中心 `fd_ctr_wrapped`。
- 所有波位完成后，复用 DBS `updateFdCtrEstimates` 的中心波位思想估计全局 `beam_pointing_bias`。
- 得到每波位真实角度 `theta_true = theta_sq + beam_pointing_bias` 后，再调用 `unwrap_prf_to_model` 得到 `fa2`，并用 `fa_shift = fa2 - fa_ctr` 修正目标多普勒后定位。

## 新增/修改文件

### 新增文件

- `include/dbs/NewProtocolReader.hpp`
- `src/dbs/NewProtocolReader.cpp`
- `include/DbsFusionTypes.hpp`
- `include/DbsFusion.hpp`
- `src/dbs/DbsFusion.cpp`
- `DBS_GMTI融合修改记录.md`

### 修改文件

- `CMakeLists.txt`
- `include/config_structs.hpp`
- `include/GMTIProcessor.hpp`
- `src/loadXML.cpp`
- `src/main.cpp`
- `src/pipe/MainCtrl.cpp`
- `src/processOnePeriod.cpp`
- `src/processPeriodsParallel.cpp`
- `src/dbs/DbsStitcher.cpp`
- `src/dbs/buildMosaic.cpp`
- `temp_config.xml`

## 配置变化

XML 增加开关：

```xml
<enable_dbs_fusion>1</enable_dbs_fusion>
```

- `true/1`：启用 DBS/GMTI 融合流程。
- `false/0` 或不配置：保持原 GMTI 流程。

## 数据结构变化

### FusionGroupContext

`FusionGroupContext` 是一个扫描周期级上下文，进入并行前按 `periodList` 预分配，每个 worker 按 slot 写入自己的波位缓存：

- `periodList`：完整扫描波位集合。
- `rd`：DBS 拼接需要的 `RDData`。
- `meta`：DBS 拼接需要的 `MetaPack`。
- `beam_meta`：每波位元数据，包括 `theta_sq`、未解模糊 `fd_ctr_wrapped`、最终 `theta_true/fa2` 等。
- `detections`：每波位 GMTI 检测原始点。
- `done`：每个 slot 的完成标记。

### DetectionRaw

每个检测点缓存：

- `prow/pcol`
- `range_m`
- `af_wrapped`：未解模糊、未加 `fa_shift` 的目标多普勒估计。
- `phase`：相对参考相位解缠后的目标相位。
- `amplitude`
- `utc_mid`

## 算法流程

### 总流程

```text
读取 XML
  |
  |-- enable_dbs_fusion = false
  |       |
  |       `-- 原 GMTI processPeriodsParallel 流程
  |
  `-- enable_dbs_fusion = true
          |
          |-- 创建 FusionGroupContext，按 periodList 预分配 slot
          |
          |-- 多 worker 并行处理每个波位
          |
          |-- 所有波位完成后估计 beam_pointing_bias
          |
          |-- 逐波位计算 theta_true 和 fa2
          |
          |-- 用 fa_shift 统一重定位 GMTI 目标
          |
          |-- 使用缓存 RDData/MetaPack 生成 DBS 拼接图
          |
          `-- 原 GMTI 结果写盘与航迹关联
```

### 单波位 worker 流程

```text
读取新协议回波 block
  |
  |-- 脉压或 PC 抽取
  |
  |-- pulse_dec 抽取
  |
  |-- 提取飞机位姿
  |
  |-- computeDoppler 得到未解模糊 fa_ctr
  |
  |-- 动态支撑域计算
  |
  |-- GPU 上传双通道数据
  |
  |-- skip=0 方位 FFT/DBS recenter
  |       |
  |       `-- 缓存 DBS RDData/MetaPack
  |
  |-- 按原 GMTI 双通道流程做对齐、FFT、DBS、相位校正、CSI、CFAR、聚类、目标筛选
  |
  `-- 缓存目标原始检测点 DetectionRaw
```

注意：此阶段不调用 `unwrap_prf_to_model`，不做最终定位。

### 组级偏差估计与解模糊

```text
选择 periodList 中心波位
  |
  |-- 使用 asin(-fd_ctr_wrapped * lambda / (2V)) 得到中心波位估计角
  |
  |-- 与 theta_sq 比较，得到 beam_pointing_bias
  |
  |-- 对每个波位：
  |       theta_true = theta_sq + beam_pointing_bias
  |       fa2 = unwrap_prf_to_model(fd_ctr_wrapped, PRF, theta_true, V, fc)
  |       fa_shift = fa2 - fd_ctr_wrapped
  |
  `-- 更新 DBS meta.fd_ctr 和 RD.fd_axis
```

### GMTI 最终重定位

对每个目标点：

```cpp
double af_ransac;
if (std::abs(k) < 1e-12) {
    af_ransac = cached.af_wrapped;
} else {
    af_ransac = (cached.phase - b) / k;
}
af_ransac += fa_shift;
```

随后沿用原 GMTI 几何定位公式：

```cpp
const double sinA = af_ransac * lambda / (2.0 * plane.V);
const double py = Rg * sinA;
const double px2 = Rg * Rg - py * py - dz * dz;
rotation_xy(py, px, flag, cosT, sinT, xP, yP);
Gaussp3RV(xP, yP, cfg.L0, lat, lng);
```

输出 `GMTIOutput::MT` 仍保持：

```text
[lat, lng, MT_nowz, xP, yP, utc, relative_dir, range]
```

### DBS 成像输出

融合模式下，所有波位处理完成后，直接使用上下文中的：

- `ctx.rd`
- `ctx.meta`

调用 DBS 拼接：

```text
estimateMosaicExtent
  |
  `-- buildMosaicGPU
          |
          `-- writeProducts
```

DBS 输出目录使用 `cfg.result_add`。

## 验证记录

已通过编译验证：

```bash
cmake --build build --target GMTI_core -j 4
cmake --build build --target GMTI_pipe_core -j 4
```

当前验证仅覆盖编译和链接。仍需要使用真实新协议回波数据做端到端运行验证，重点检查：

- `beam_pointing_bias` 数值是否合理。
- 每波位 `fa_shift` 是否符合预期。
- GMTI 重定位结果与旧流程/人工基准是否一致。
- DBS 拼接图是否正常输出、方向和范围是否正确。

## 2026-06-03 内存问题排查与保护

### 初步判断

本次内存爆掉的高风险点不是递归或无限循环，而是 DBS 拼图阶段输出网格过大：

```text
estimateMosaicExtent 生成 grid.x/grid.y
  |
  `-- buildMosaicGPU 按 nx * ny 分配整幅 mosaic
```

原流程按 `P.out_res_m = 1m` 生成规则网格，没有对 `nx * ny` 设置上限。如果地理覆盖范围较大，容易产生几千万到几亿像素。旧 GPU 拼图路径还会额外分配 `r_save/fd_save/x_rel_save/y_rel_save` 四张普通运行中不使用的大图，进一步放大主机内存占用。

### 本次修改

- `processPeriodsParallelFusion` 增加组级日志：
  - period 数量
  - worker 数量
  - 当前 CUDA free/total
  - `pulse_num/pulse_dec/rg_len`
  - DBS 幅度缓存估算
- `runDbsFusionImaging` 增加 DBS 成像入口日志：
  - beam 数量
  - RD 缓存尺寸
  - RD 幅度缓存估算
  - PRF/Rmin/fs/out_res/useGpu
- `estimateMosaicExtent` 增加输入合法性、速度合法性、bounds 合法性检查。
- `estimateMosaicExtent` 在创建网格前估算 `nx * ny`，默认限制为 `60000000` 像素。
- 如果请求分辨率下网格超过限制，自动放粗 DBS 输出分辨率，并打印 warning。
- `buildMosaicGPU` 在分配前再次检查像素数，并打印 host/device/amp pack 内存估算。
- `buildMosaicGPU` 普通编译下不再分配未使用的 `r_save/fd_save/x_rel_save/y_rel_save` 四张大图。
- `buildMosaicGPU` 对 `cudaMalloc` 增加失败检查，失败时释放已经申请的显存并返回 false。

### 2026-06-03 参数来源修正

融合后 GMTI 是主流程，DBS 成像不再单独读取 DBS XML 的 `Params P`。DBS 拼图链路改为直接使用 GMTI 已读取的 `Config cfg`：

- `raw_fenbianlv` -> `cfg.dbs_out_res_m`
- `n_tiaoguo` -> `cfg.dbs_beam_skip`
- `len_tiaoguo` -> `cfg.dbs_range_skip`
- 可选 `dbs_interp_mode` -> `cfg.dbs_interp_mode`

涉及函数：

```text
runDbsFusionImaging(cfg)
  |
  |-- estimateMosaicExtent(cfg, ...)
  |-- buildMosaicGPU(cfg, ...)
  `-- writeProducts(cfg, ...)
```

因此 DBS 融合成像中的输出分辨率不再使用 `Params::out_res_m` 默认值，而是使用 GMTI XML 中的 `raw_fenbianlv`。

### 可调参数

默认最大 DBS 拼图像素数：

```bash
DBS_MAX_MOSAIC_PIXELS=60000000
```

如 32G 内存设备需要更高 DBS 输出分辨率，可在运行前临时提高，例如：

```bash
export DBS_MAX_MOSAIC_PIXELS=100000000
```

提高该值会增加主机内存、显存、输出转置和 PNG 写图阶段的峰值内存。

### 后续实测重点日志

运行时重点观察：

```text
[fusion] group start: ...
[fusion][dbs] imaging start: ...
[fusion][dbs][extent] ...
[fusion][dbs][extent][warn] ...
[fusion][dbs][buildGPU] memory estimate: ...
```

如果仍然爆内存，优先根据最后一条成功打印的日志判断卡在哪个阶段。
