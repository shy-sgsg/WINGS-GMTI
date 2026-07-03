# GMTI TrackManager 可视化工具说明

## 1. 功能简介

`tools/visualize_track_manager.py` 用于读取在线 `TrackManager` 生成的调试 CSV，并可视化航迹管理过程。它可以帮助判断短航迹 ID 过多、断轨、漏检保留、输出点过滤等问题来源。

主要功能包括：

1. 航迹全局路径；
2. 单周期快照；
3. 航迹生命周期；
4. 速度变化；
5. 短航迹统计；
6. 每周期航迹生长图；
7. GIF 动画；
8. ROI 局部区域查看。

## 2. 输入文件说明

输入目录通常是：

```text
result_final/track_debug
```

该目录由在线 `TrackManager` 调试快照生成，包含三个 CSV。

### track_frames.csv

每个处理周期一行，用于查看整体统计。

字段：

```text
result_id
frame_utc
num_detections
num_tracks
num_tentative
num_confirmed
num_coasted
num_deleted
num_outputs
num_new_tracks
num_matched_tracks
num_unmatched_detections
```

含义：

- `result_id`：当前周期编号，例如 16、17、18。
- `frame_utc`：当前帧时间。
- `num_detections`：当前周期原始检测点数量。
- `num_tracks`：当前周期 TrackManager 内部航迹总数，包括即将删除的 Deleted。
- `num_tentative`：Tentative 航迹数。
- `num_confirmed`：Confirmed 航迹数。
- `num_coasted`：Coasted 航迹数。
- `num_deleted`：Deleted 航迹数，删除前会保留一行快照。
- `num_outputs`：最终输出目标数，即写入 `GMTIxx_track.bin` 的数量。
- `num_new_tracks`：当前周期新建 Tentative 数量。
- `num_matched_tracks`：当前周期匹配到检测的航迹数。
- `num_unmatched_detections`：当前周期未被任何航迹关联的检测点数量。

### track_detections.csv

每个周期、每个原始检测点一行，用于观察检测点是否被航迹吸收。

字段：

```text
result_id
det_index
matched
matched_track_id
e
n
lat
lon
utc
range
direction
```

含义：

- `result_id`：当前周期编号。
- `det_index`：当前周期检测点索引。
- `matched`：1 表示该检测点被某条 TrackManager 航迹使用。
- `matched_track_id`：匹配到的航迹 ID；未匹配时为 -1。
- `e,n`：高斯投影平面坐标，单位 m。
- `lat,lon`：经纬度。
- `utc`：检测点时间。
- `range`：目标到飞机平台距离，只保存用于诊断，不参与关联代价。
- `direction`：相对飞机航向角方向，只保存用于诊断，不参与关联代价。

### track_states.csv

这是核心文件，每个周期、每条 TrackManager 内部航迹一行，用于查看完整生命周期。

字段：

```text
result_id
track_id
state
matched_this_frame
matched_det_index
e
n
lat
lon
ve
vn
speed
heading
utc
range
direction
age
hit_count
miss_count
consecutive_hits
recent_hits
linearity
is_output
```

含义：

- `track_id`：TrackManager 持久航迹 ID。
- `state`：TrackManager 内部状态，取值为 `Tentative`、`Confirmed`、`Coasted`、`Deleted`。
- `matched_this_frame=1`：当前周期匹配到了真实检测。
- `matched_det_index`：匹配到的检测点索引，未匹配时为 -1。
- `e,n,lat,lon`：航迹当前坐标。
- `ve,vn`：EN 平面速度分量。
- `speed`：速度，仅用于诊断和显示。
- `heading`：EN 平面运动航向。
- `utc`：航迹当前时间。
- `range,direction`：最近一次检测属性，只保存和输出，不参与关联代价。
- `age`：航迹生命周期周期数。
- `hit_count`：累计命中次数。
- `miss_count`：连续漏检次数。
- `consecutive_hits`：连续命中次数。
- `recent_hits`：确认窗口内命中次数。
- `linearity`：近期轨迹直线度。
- `is_output=1`：该点是最终写入 `GMTIxx_track.bin` 的目标。

重点说明：

- `track_states.csv` 是诊断 TrackManager 的核心文件。
- `state` 表示 TrackManager 内部状态，不等同于最终输出。
- `is_output=1` 表示最终写入 `GMTIxx_track.bin` 的目标。
- `matched_this_frame=1` 表示当前周期匹配到了真实检测。
- `Coasted` 是预测保留状态，不代表当前帧真实检测，也不会输出到 `GMTIxx_track.bin`。

## 3. 输出文件说明

脚本默认输出到：

```text
<debug-dir>/figs
```

可用 `--out-dir` 指定其他目录。

输出文件：

```text
track_global_paths.<format>
track_snapshot_last.<format>
track_outputs_global_paths.<format>
track_outputs_snapshot_last.<format>
track_state_timeline.<format>
track_speed_timeline.<format>
track_short_track_summary.<format>
frames/frame_000016.<format>
output_frames/output_frame_000016.<format>
track_growth.gif
track_outputs_growth.gif
```

用途：

- `track_global_paths`：查看全局航迹路径、状态点、最终输出点和原始检测点。
- `track_snapshot_last`：查看最后一个周期或指定周期处理完后的 TrackManager 全局状态。
- `track_outputs_global_paths`：只查看最终确认并输出的点，以及由这些点连成的输出航迹。
- `track_outputs_snapshot_last`：只查看最后一个周期或指定周期的最终输出点和输出航迹历史。
- `track_state_timeline`：查看每条航迹随周期变化的状态，适合观察频繁新建 ID、短 Tentative、过快 Deleted。
- `track_speed_timeline`：查看每条航迹命中检测时的速度变化，默认不画 Coasted 预测速度。
- `track_short_track_summary`：统计短航迹数量，用于对比优化前后短 ID 是否减少。
- `frames/frame_000016.<format>`：每个周期一张快照，用于逐帧检查航迹生长过程。
- `output_frames/output_frame_000016.<format>`：每个周期一张最终输出专用快照，只显示 `is_output=1` 的点和航迹。
- `track_growth.gif`：由 PNG 帧合成的航迹生长动画。
- `track_outputs_growth.gif`：由 PNG 帧合成的最终输出航迹生长动画。

## 4. 常用命令示例

### 全局可视化

```bash
python tools/visualize_track_manager.py \
  --debug-dir result_final/track_debug \
  --coord en \
  --format pdf
```

### 过滤短航迹

```bash
python tools/visualize_track_manager.py \
  --debug-dir result_final/track_debug \
  --coord en \
  --min-len 3 \
  --format pdf
```

### 指定周期范围

```bash
python tools/visualize_track_manager.py \
  --debug-dir result_final/track_debug \
  --coord en \
  --start-id 16 \
  --end-id 40 \
  --format pdf
```

### 指定局部 ROI

```bash
python tools/visualize_track_manager.py \
  --debug-dir result_final/track_debug \
  --coord en \
  --xmin 12000 --xmax 18000 \
  --ymin 35000 --ymax 42000 \
  --min-len 3 \
  --format pdf
```

### 生成每周期图片和 GIF

```bash
python tools/visualize_track_manager.py \
  --debug-dir result_final/track_debug \
  --coord en \
  --min-len 3 \
  --make-frames \
  --make-gif \
  --format png
```

### 中心点半径局部查看

```bash
python tools/visualize_track_manager.py \
  --debug-dir result_final/track_debug \
  --coord en \
  --center-x 15000 \
  --center-y 38000 \
  --radius 3000 \
  --min-len 3 \
  --make-frames \
  --make-gif \
  --format png
```

## 5. 图中标注说明

图中符号含义：

空间航迹图，包括 `track_global_paths`、`track_snapshot_last`、逐周期 `frames` 和 GIF：

1. 不同颜色表示不同 `track_id`。
2. marker 形状表示 `state`。
3. 圆点表示 `Tentative`。
4. 方块表示 `Confirmed`。
5. 三角表示 `Coasted`。
6. 叉号表示 `Deleted`。
7. 黑边或更大的点表示当前周期最终输出目标，即 `is_output=1`。
8. 灰色点表示原始检测。
9. 空心点表示未匹配检测。
10. 末端文字 `ID=xx` 表示航迹编号。
11. `v=xx m/s` 表示速度，仅供参考。
12. `Coasted` 点是预测点，不是当前真实检测。

状态时间线 `track_state_timeline`：

1. 颜色表示 `state`。
2. marker 形状也按 `state` 区分。
3. 黑边表示该周期该航迹点最终输出到 `GMTIxx_track.bin`。

速度时间线 `track_speed_timeline`：

1. 颜色表示 `track_id`。
2. 默认只画 `matched_this_frame=1` 的点，避免把 Coasted 预测速度误认为真实检测速度。

## 6. 样式控制

空间图中颜色只表示 `track_id`，形状表示 `state`，黑边表示 final output。线宽、marker 大小、字体大小都是绘图显示单位，不是物理距离单位。

ROI 局部放大后，如果图面太挤，可以使用样式 preset：

```bash
--style-preset roi
```

或更紧凑的：

```bash
--style-preset dense
```

也可以手动调节：

```bash
--track-linewidth 0.6
--marker-scale 0.5
--detection-marker-size 4
--label-fontsize 6
```

示例：

```bash
python tools/visualize_track_manager.py \
  --debug-dir result_final/track_debug \
  --coord en \
  --xmin 652500 --xmax 655000 \
  --ymin 4346000 --ymax 4347000 \
  --min-len 4 \
  --format pdf \
  --style-preset dense \
  --track-linewidth 0.6 \
  --marker-scale 0.5 \
  --detection-marker-size 4 \
  --label-fontsize 6
```

如果指定了 ROI 参数并且没有显式指定 `--style-preset`，脚本会自动使用 `roi` 样式。可见航迹数超过 `--max-label-tracks` 时，会自动关闭 ID 标注并打印提示，避免文字遮挡。

逐周期图片和 GIF 默认使用更克制的 frame 样式。可用 `--hide-labels-in-frames` 关闭逐帧 ID 标注。

## 7. 常见问题

### 为什么图上很多短航迹？

可能原因：

1. 虚警多；
2. 关联门限太严；
3. 目标漏检；
4. TrackManager 还在 Tentative 阶段；
5. 只看最终输出点会隐藏内部航迹过程。

建议先使用：

```bash
--min-len 3
```

并重点查看：

```text
track_short_track_summary
track_state_timeline
```

### 为什么有些点是 Coasted？

`Coasted` 表示航迹当前周期未匹配到检测，但 TrackManager 暂时保留预测状态。它不会输出到 `GMTIxx_track.bin`，也不代表当前帧真实检测。

### 为什么 GIF 必须用 PNG？

GIF 是栅格动画，PDF/SVG 是矢量静态图。生成 GIF 时脚本会额外生成 PNG 帧，再把 PNG 帧按 `result_id` 升序合成为动画。

## 8. 参数列表

```text
--debug-dir
--coord
--out-dir
--format
--result-id
--start-id
--end-id
--min-len
--show-detections
--hide-detections
--show-speed
--label-tracks
--no-label-tracks
--make-frames
--make-gif
--gif-duration-ms
--keep-frame-png
--xmin
--xmax
--ymin
--ymax
--roi-mode
--center-x
--center-y
--radius
--padding-ratio
--track-linewidth
--track-alpha
--marker-scale
--detection-marker-size
--detection-alpha
--output-marker-scale
--output-edge-width
--label-fontsize
--legend-fontsize
--title-fontsize
--axis-label-fontsize
--tick-fontsize
--hide-labels-in-frames
--max-label-tracks
--style-preset
--fig-width
--fig-height
--frame-fig-width
--frame-fig-height
```

说明：

- `--debug-dir`：包含三份调试 CSV 的目录。
- `--coord`：`en` 或 `geo`，默认 `en`。
- `--out-dir`：输出图目录，默认 `<debug-dir>/figs`。
- `--format`：`png`、`pdf`、`svg`。
- `--result-id`：只画某个周期。
- `--start-id / --end-id`：周期范围。
- `--min-len`：按 `matched_this_frame=1` 次数过滤短航迹。
- `--show-detections / --hide-detections`：显示或隐藏原始检测点。
- `--show-speed`：在快照 ID 标注中附带速度。
- `--label-tracks / --no-label-tracks`：显示或隐藏航迹 ID。
- `--make-frames`：为每个周期生成快照图。
- `--make-gif`：生成 `track_growth.gif`。
- `--gif-duration-ms`：GIF 每帧时长，默认 500 ms。
- `--keep-frame-png`：是否保留 GIF 使用的 PNG 帧，默认 true，可传 `false`。
- `--xmin / --xmax / --ymin / --ymax`：矩形 ROI。
- `--roi-mode any`：航迹任一点进入 ROI，则保留整条历史。
- `--roi-mode points`：只显示 ROI 内点和线段。
- `--center-x / --center-y / --radius`：EN 坐标下按中心点和半径生成 ROI。
- `--padding-ratio`：坐标范围四周 padding，默认 0.05。
- `--track-linewidth`：航迹历史线宽，默认 1.2。
- `--track-alpha`：航迹线透明度，默认 0.85。
- `--marker-scale`：统一缩放航迹状态 marker 大小，默认 1.0。
- `--detection-marker-size`：检测点 marker 大小，默认 10。
- `--detection-alpha`：检测点透明度，默认 0.35。
- `--output-marker-scale`：final output 点相对普通 Confirmed 点的放大倍数，默认 1.5。
- `--output-edge-width`：final output 黑边线宽，默认 1.2。
- `--label-fontsize`：航迹 ID 和速度标签字体大小，默认 8。
- `--legend-fontsize`：图例字体大小，默认 8。
- `--title-fontsize`：标题字体大小，默认 11。
- `--axis-label-fontsize`：坐标轴字体大小，默认 10。
- `--tick-fontsize`：刻度字体大小，默认 8。
- `--hide-labels-in-frames`：生成逐周期 frames/GIF 时隐藏航迹 ID 标签。
- `--max-label-tracks`：最多给多少条航迹标注 ID，默认 80。
- `--style-preset`：样式 preset，可选 `default`、`global`、`roi`、`dense`。
- `--fig-width / --fig-height`：全局图尺寸，默认 12 x 10。
- `--frame-fig-width / --frame-fig-height`：逐周期帧图尺寸，默认 10 x 8。
