#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
BUILD_DIR="$ROOT/build"
OUT_ROOT="$ROOT/outputs/stage2_motion_comp_matrix_$(date +%Y%m%d_%H%M%S)"

# 默认先跑最关键的 5 个，确认链路没问题后再改成：
# CASES_TO_RUN="C0 C1 C2 C3 C4 C5 C6 C7 C8 C9"
CASES_TO_RUN="C0 C1 C2 C3 C4"

mkdir -p "$OUT_ROOT"

echo "[INFO] ROOT=$ROOT"
echo "[INFO] OUT_ROOT=$OUT_ROOT"

echo "[INFO] Build..."
cmake --build "$BUILD_DIR" -j2

SIM_EXE="$BUILD_DIR/simulate_stage2_statistical"
GMTI_EXE="$BUILD_DIR/GMTI_core"

if [ ! -x "$SIM_EXE" ]; then
  echo "[ERR] not found executable: $SIM_EXE"
  exit 1
fi

if [ ! -x "$GMTI_EXE" ]; then
  echo "[ERR] not found executable: $GMTI_EXE"
  exit 1
fi

# 自动寻找 stage2_config.json
STAGE2_CONFIG=""
for p in \
  "$ROOT/outputs/stage2/config/stage2_config.json" \
  "$ROOT/simulator/stage2/config/stage2_config.json" \
  "$ROOT/simulator/stage2/stage2_config.json" \
  "$ROOT/config/stage2_config.json" \
  "$ROOT/stage2_config.json"
do
  if [ -f "$p" ]; then
    STAGE2_CONFIG="$p"
    break
  fi
done

if [ -z "$STAGE2_CONFIG" ]; then
  STAGE2_CONFIG="$(find "$ROOT" -name 'stage2_config.json' | head -n 1 || true)"
fi

if [ -z "$STAGE2_CONFIG" ] || [ ! -f "$STAGE2_CONFIG" ]; then
  echo "[ERR] stage2_config.json not found"
  echo "Please set STAGE2_CONFIG manually in this script."
  exit 1
fi

echo "[INFO] STAGE2_CONFIG=$STAGE2_CONFIG"

# 自动寻找目标配置；如果 stage2_config 内部已经写了 target_config_path，也可以不传。
TARGET_CONFIG=""
for p in \
  "$ROOT/outputs/stage2/config/cooperative_target.json" \
  "$ROOT/simulator/target_injection/config/cooperative_target.json" \
  "$ROOT/simulator/target_injection/cooperative_target.json" \
  "$ROOT/config/cooperative_target.json" \
  "$ROOT/cooperative_target.json"
do
  if [ -f "$p" ]; then
    TARGET_CONFIG="$p"
    break
  fi
done

if [ -z "$TARGET_CONFIG" ]; then
  TARGET_CONFIG="$(find "$ROOT" \( -name '*target*.json' -o -name '*cooperative*.json' \) | head -n 1 || true)"
fi

if [ -n "$TARGET_CONFIG" ] && [ -f "$TARGET_CONFIG" ]; then
  echo "[INFO] TARGET_CONFIG=$TARGET_CONFIG"
else
  echo "[WARN] target config not found. If stage2_config already contains target_config_path, this may still work."
  TARGET_CONFIG=""
fi

patch_xml_for_gmti() {
  local xml_in="$1"
  local result_dir="$2"
  local motion_enable="$3"
  local fusion_enable="$4"

  python3 - "$xml_in" "$result_dir" "$motion_enable" "$fusion_enable" <<'PY'
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

xml_path = Path(sys.argv[1])
result_dir = sys.argv[2]
motion_enable = sys.argv[3]
fusion_enable = sys.argv[4]

tree = ET.parse(xml_path)
root = tree.getroot()

def set_tag(name, value):
    node = root.find(".//" + name)
    if node is None:
        node = ET.SubElement(root, name)
    node.text = str(value)

# 输出目录
set_tag("result_add", result_dir)

# 运动补偿相关
set_tag("motion_comp_enable", motion_enable)
set_tag("motion_comp_use_row_doppler", "1")
set_tag("motion_comp_debug", "1")

# 解析解模式，如果你的代码没有这个字段，不影响；如果有则启用
set_tag("motion_comp_solver", "analytic")

# 符号：当前正向注入的 angle(F1*conj(F2)) 运动相位随正径向速度为正，
# 但 Doppler 轴上同一运动量表现为负频移。
set_tag("ati_phase_to_velocity_sign", "1")
set_tag("motion_doppler_axis_sign", "-1")

# 兼容旧字段
set_tag("ati_velocity_sign", "1")

# 速度限幅，先设成略小于或接近无模糊范围，避免绕相位的目标误导
set_tag("ati_vmax_mps", "2.0")
set_tag("ati_phase_bias_rad", "0.0")

# 融合路径开关
set_tag("enable_dbs_fusion", fusion_enable)

tree.write(xml_path, encoding="utf-8", xml_declaration=True)
PY
}

run_one_case() {
  local case_id="$1"
  local scene_mode="$2"
  local target_enabled="$3"
  local vr_mps="$4"
  local area_count="$5"
  local area_power="$6"
  local strong_count="$7"
  local strong_min="$8"
  local strong_max="$9"
  local line_count="${10}"
  local line_points="${11}"
  local line_rcs="${12}"
  local noise_power="${13}"
  local motion_enable="${14}"
  local fusion_enable="${15}"

  local case_dir="$OUT_ROOT/$case_id"
  local sim_dir="$case_dir/sim"
  local result_dir="$case_dir/gmti_result"
  local log_dir="$case_dir/logs"

  mkdir -p "$sim_dir" "$result_dir" "$log_dir"

  echo ""
  echo "============================================================"
  echo "[CASE] $case_id"
  echo "scene=$scene_mode target=$target_enabled vr=$vr_mps area=$area_count strong=$strong_count line=$line_count noise=$noise_power motion=$motion_enable fusion=$fusion_enable"
  echo "============================================================"

  local target_args=()
  if [ "$target_enabled" = "1" ]; then
    target_args+=(--target-enabled 1)
    if [ -n "$TARGET_CONFIG" ]; then
      target_args+=(--target-config "$TARGET_CONFIG")
    fi
    target_args+=(--moving-target-beam-id 31)
    target_args+=(--moving-target-expected-bin 880)
    target_args+=(--moving-target-velocity-mode radial)
    target_args+=(--moving-target-radial-speed-mps "$vr_mps")
    target_args+=(--moving-target-tangential-speed-mps 0)
    target_args+=(--moving-target-rcs-db 40)
  else
    target_args+=(--target-enabled 0)
  fi

  "$SIM_EXE" \
    --stage2-config "$STAGE2_CONFIG" \
    --output-dir "$sim_dir" \
    --scene-mode "$scene_mode" \
    --period-start 0 \
    --period-count 1 \
    --clutter-amplitude-scale 1 \
    --area-clutter-scatterer-count "$area_count" \
    --area-clutter-mean-power "$area_power" \
    --area-clutter-texture-sigma 0.5 \
    --strong-scatterer-count "$strong_count" \
    --strong-rcs-db-min "$strong_min" \
    --strong-rcs-db-max "$strong_max" \
    --line-scatterer-count "$line_count" \
    --line-points-per-line "$line_points" \
    --line-rcs-db "$line_rcs" \
    --noise-power "$noise_power" \
    "${target_args[@]}" \
    > "$log_dir/simulate_stage2.stdout.log" \
    2> "$log_dir/simulate_stage2.stderr.log"

  local xml="$sim_dir/config/temp_config_stage2_newsystem.xml"
  if [ ! -f "$xml" ]; then
    echo "[ERR] XML not generated: $xml"
    exit 1
  fi

  patch_xml_for_gmti "$xml" "$result_dir" "$motion_enable" "$fusion_enable"

  "$GMTI_EXE" "$xml" \
    > "$log_dir/GMTI_core.stdout.log" \
    2> "$log_dir/GMTI_core.stderr.log" || {
      echo "[ERR] GMTI_core failed for $case_id"
      tail -n 80 "$log_dir/GMTI_core.stderr.log" || true
      exit 1
    }

  # 尝试寻找 detection_results.csv
  local det_csv=""
  for p in \
    "$result_dir/detection_results.csv" \
    "$sim_dir/detection_results.csv" \
    "$case_dir/detection_results.csv"
  do
    if [ -f "$p" ]; then
      det_csv="$p"
      break
    fi
  done

  if [ -z "$det_csv" ]; then
    det_csv="$(find "$case_dir" -name 'detection_results.csv' | head -n 1 || true)"
  fi

  if [ -n "$det_csv" ] && [ -f "$det_csv" ]; then
    cp "$det_csv" "$case_dir/detection_results.csv"
    echo "[INFO] detection_results.csv: $det_csv"
  else
    echo "[WARN] detection_results.csv not found for $case_id"
  fi

  echo "[DONE] $case_id"
}

# ------------------------------------------------------------
# Case 定义
# 参数：
# case_id scene_mode target_enabled vr area_count area_power strong_count strong_min strong_max line_count line_points line_rcs noise motion_enable fusion_enable
# ------------------------------------------------------------

for case_name in $CASES_TO_RUN; do
  case "$case_name" in
    C0)
      # 冒烟：极简杂波，无目标
      run_one_case C0_smoke_clutter_only \
        clutter_only 0 0 \
        2 1.0 \
        1 10 20 \
        0 0 0 \
        0 \
        0 0
      ;;
    C1)
      # 纯杂波：检查 p38 稳定性
      run_one_case C1_clutter_p38_check \
        clutter_only 0 0 \
        100 0.5 \
        0 10 20 \
        0 0 0 \
        0 \
        0 0
      ;;
    C2)
      # 静止目标：开启补偿后不应乱补
      run_one_case C2_static_target_motion_on \
        clutter_only 1 0 \
        100 0.5 \
        0 10 20 \
        0 0 0 \
        0 \
        1 0
      ;;
    C3)
      # +1 m/s：核心正速度测试
      run_one_case C3_pos1_target_motion_on \
        clutter_only 1 1.0 \
        100 0.5 \
        0 10 20 \
        0 0 0 \
        0 \
        1 0
      ;;
    C4)
      # -1 m/s：核心负速度测试
      run_one_case C4_neg1_target_motion_on \
        clutter_only 1 -1.0 \
        100 0.5 \
        0 10 20 \
        0 0 0 \
        0 \
        1 0
      ;;
    C5)
      # +0.5 m/s：低速测试
      run_one_case C5_pos05_target_motion_on \
        clutter_only 1 0.5 \
        100 0.5 \
        0 10 20 \
        0 0 0 \
        0 \
        1 0
      ;;
    C6)
      # +1.3 m/s：接近 ATI 无模糊边界
      run_one_case C6_pos13_target_motion_on \
        clutter_only 1 1.3 \
        100 0.5 \
        0 10 20 \
        0 0 0 \
        0 \
        1 0
      ;;
    C7)
      # 面杂波增强
      run_one_case C7_pos1_area300_motion_on \
        clutter_only 1 1.0 \
        300 1.0 \
        0 10 20 \
        0 0 0 \
        0 \
        1 0
      ;;
    C8)
      # 面杂波 + 强散射点
      run_one_case C8_pos1_area300_strong3_motion_on \
        clutter_only 1 1.0 \
        300 1.0 \
        3 10 20 \
        0 0 0 \
        0 \
        1 0
      ;;
    C9)
      # 面杂波 + 强散射点 + 线状散射
      run_one_case C9_pos1_complex_motion_on \
        clutter_only 1 1.0 \
        300 1.0 \
        3 10 20 \
        1 20 10 \
        0 \
        1 0
      ;;
    *)
      echo "[ERR] unknown case: $case_name"
      exit 1
      ;;
  esac
done

# ------------------------------------------------------------
# 汇总结果
# ------------------------------------------------------------
python3 - "$OUT_ROOT" <<'PY'
import csv
import math
import sys
from pathlib import Path

root = Path(sys.argv[1])
summary_path = root / "summary_motion_comp.csv"

fields = [
    "case",
    "det_count",
    "best_det_id",
    "beam_id",
    "row",
    "range_bin",
    "phase_rad",
    "phi_res_rad",
    "denom",
    "v_radial_mps",
    "af_total_hz",
    "af_motion_hz",
    "af_geometry_hz",
    "motion_comp_enable",
    "motion_comp_valid",
    "motion_comp_used",
    "loc_used_mode",
    "old_e",
    "old_n",
    "new_e",
    "new_n",
]

def fnum(x):
    try:
        if x is None or x == "":
            return math.nan
        return float(x)
    except Exception:
        return math.nan

def pick_best(rows):
    if not rows:
        return {}
    # 优先选 range_bin 最接近 880 的检测点；如果没有 range_bin，就取第一行
    def score(r):
        rb = fnum(r.get("range_bin"))
        if math.isfinite(rb):
            return abs(rb - 880)
        return 1e18
    return min(rows, key=score)

with summary_path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()

    for case_dir in sorted(root.iterdir()):
        if not case_dir.is_dir():
            continue
        det_csv = case_dir / "detection_results.csv"
        rows = []
        if det_csv.exists():
            with det_csv.open("r", encoding="utf-8", newline="") as df:
                rows = list(csv.DictReader(df))
        best = pick_best(rows)

        out = {"case": case_dir.name, "det_count": len(rows)}
        for k in fields:
            if k in ("case", "det_count"):
                continue
            out[k] = best.get(k, "") if best else ""
        writer.writerow(out)

print(f"[SUMMARY] {summary_path}")
print(summary_path.read_text(encoding="utf-8"))
PY

echo ""
echo "[ALL DONE]"
echo "Output root:"
echo "$OUT_ROOT"
echo ""
echo "Summary:"
echo "$OUT_ROOT/summary_motion_comp.csv"
