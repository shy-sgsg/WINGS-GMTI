#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

STAGE2_CONFIG="stage2_config.json"
TARGET_CONFIG="targets.json"
OUTPUT_DIR="outputs/stage2_oneclick"
SCENE_MODE="full"
TARGET_ENABLED="true"
PERIOD_COUNT="1"
RUN_ALGORITHM="false"
BUILD_TYPE=""
SINGLE_POINT_BEAM_ID=""
SINGLE_POINT_EXPECTED_BIN=""
MOVING_TARGET_BEAM_ID=""
MOVING_TARGET_EXPECTED_BIN=""
MOVING_TARGET_SPEED_MPS=""
MOVING_TARGET_VELOCITY_MODE=""
MOVING_TARGET_VE_MPS=""
MOVING_TARGET_VN_MPS=""
MOVING_TARGET_RADIAL_SPEED_MPS=""
MOVING_TARGET_TANGENTIAL_SPEED_MPS=""
MOVING_TARGET_RCS_DB=""
CLUTTER_AMPLITUDE_SCALE=""
AREA_CLUTTER_SCATTERER_COUNT=""
AREA_CLUTTER_MEAN_POWER=""
AREA_CLUTTER_TEXTURE_SIGMA=""
STRONG_SCATTERER_COUNT=""
STRONG_RCS_DB_MIN=""
STRONG_RCS_DB_MAX=""
LINE_SCATTERER_COUNT=""
LINE_POINTS_PER_LINE=""
LINE_RCS_DB=""
NOISE_POWER=""

usage() {
  cat <<'EOF'
Usage:
  scripts/run_stage2_one_click_test.sh [options]

Options:
  --stage2-config PATH       Default: stage2_config.json
  --target-config PATH       Default: targets.json
  --output-dir DIR           Default: outputs/stage2_oneclick
  --scene-mode MODE          full | clutter_only | noise_only | point_target_only | area_clutter_only | target_only
                             aliases: cooperative_target_only, empty_target
  --target-enabled true|false
  --period-count N           Default: 1
  --single-point-beam-id B     1-based, same as algorithm detection beam_id
  --single-point-expected-bin N
  --moving-target-beam-id B    1-based, same as algorithm detection beam_id
  --moving-target-expected-bin N
  --moving-target-speed-mps V
  --moving-target-velocity-mode enu|radial|tangential
  --moving-target-ve-mps V
  --moving-target-vn-mps V
  --moving-target-radial-speed-mps V
  --moving-target-tangential-speed-mps V
  --moving-target-rcs-db RCS
  --clutter-amplitude-scale S
  --area-clutter-scatterer-count N
  --area-clutter-mean-power P
  --area-clutter-texture-sigma S
  --strong-scatterer-count N
  --strong-rcs-db-min DB
  --strong-rcs-db-max DB
  --line-scatterer-count N
  --line-points-per-line N
  --line-rcs-db DB
  --noise-power P
  --run-algorithm true|false Default: false
  --build-type TYPE          Optional CMake build type

Examples:
  scripts/run_stage2_one_click_test.sh --scene-mode point_target_only --target-enabled false --output-dir outputs/stage2_point_ci
  scripts/run_stage2_one_click_test.sh --scene-mode full --target-enabled true --period-count 1 --run-algorithm true
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stage2-config) STAGE2_CONFIG="$2"; shift 2 ;;
    --target-config) TARGET_CONFIG="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --scene-mode) SCENE_MODE="$2"; shift 2 ;;
    --target-enabled) TARGET_ENABLED="$2"; shift 2 ;;
    --period-count) PERIOD_COUNT="$2"; shift 2 ;;
    --single-point-beam-id) SINGLE_POINT_BEAM_ID="$2"; shift 2 ;;
    --single-point-expected-bin) SINGLE_POINT_EXPECTED_BIN="$2"; shift 2 ;;
    --moving-target-beam-id) MOVING_TARGET_BEAM_ID="$2"; shift 2 ;;
    --moving-target-expected-bin) MOVING_TARGET_EXPECTED_BIN="$2"; shift 2 ;;
    --moving-target-speed-mps) MOVING_TARGET_SPEED_MPS="$2"; shift 2 ;;
    --moving-target-velocity-mode) MOVING_TARGET_VELOCITY_MODE="$2"; shift 2 ;;
    --moving-target-ve-mps) MOVING_TARGET_VE_MPS="$2"; shift 2 ;;
    --moving-target-vn-mps) MOVING_TARGET_VN_MPS="$2"; shift 2 ;;
    --moving-target-radial-speed-mps) MOVING_TARGET_RADIAL_SPEED_MPS="$2"; shift 2 ;;
    --moving-target-tangential-speed-mps) MOVING_TARGET_TANGENTIAL_SPEED_MPS="$2"; shift 2 ;;
    --moving-target-rcs-db) MOVING_TARGET_RCS_DB="$2"; shift 2 ;;
    --clutter-amplitude-scale) CLUTTER_AMPLITUDE_SCALE="$2"; shift 2 ;;
    --area-clutter-scatterer-count) AREA_CLUTTER_SCATTERER_COUNT="$2"; shift 2 ;;
    --area-clutter-mean-power) AREA_CLUTTER_MEAN_POWER="$2"; shift 2 ;;
    --area-clutter-texture-sigma) AREA_CLUTTER_TEXTURE_SIGMA="$2"; shift 2 ;;
    --strong-scatterer-count) STRONG_SCATTERER_COUNT="$2"; shift 2 ;;
    --strong-rcs-db-min) STRONG_RCS_DB_MIN="$2"; shift 2 ;;
    --strong-rcs-db-max) STRONG_RCS_DB_MAX="$2"; shift 2 ;;
    --line-scatterer-count) LINE_SCATTERER_COUNT="$2"; shift 2 ;;
    --line-points-per-line) LINE_POINTS_PER_LINE="$2"; shift 2 ;;
    --line-rcs-db) LINE_RCS_DB="$2"; shift 2 ;;
    --noise-power) NOISE_POWER="$2"; shift 2 ;;
    --run-algorithm) RUN_ALGORITHM="$2"; shift 2 ;;
    --build-type) BUILD_TYPE="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[oneclick][ERR] unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "$OUTPUT_DIR/logs" "$OUTPUT_DIR/reports"

echo "[oneclick] configure/build"
if [[ -n "$BUILD_TYPE" ]]; then
  cmake -S . -B build -DCMAKE_BUILD_TYPE="$BUILD_TYPE" | tee "$OUTPUT_DIR/logs/cmake_configure.log"
else
  cmake -S . -B build | tee "$OUTPUT_DIR/logs/cmake_configure.log"
fi
cmake --build build --target simulate_stage2_statistical -j2 | tee "$OUTPUT_DIR/logs/build_stage2.log"

if [[ "$RUN_ALGORITHM" == "true" ]]; then
  cmake --build build --target GMTI_core -j2 | tee "$OUTPUT_DIR/logs/build_gmti_core.log"
fi

echo "[oneclick] generate stage2 data: scene=${SCENE_MODE}, target=${TARGET_ENABLED}, periods=${PERIOD_COUNT}"
SIM_ARGS=(
  ./build/simulate_stage2_statistical
  --stage2-config "$STAGE2_CONFIG"
  --target-config "$TARGET_CONFIG"
  --output-dir "$OUTPUT_DIR"
  --scene-mode "$SCENE_MODE"
  --target-enabled "$TARGET_ENABLED"
  --period-count "$PERIOD_COUNT"
  --validate true
)
if [[ -n "$SINGLE_POINT_BEAM_ID" ]]; then
  SIM_ARGS+=(--single-point-beam-id "$SINGLE_POINT_BEAM_ID")
fi
if [[ -n "$SINGLE_POINT_EXPECTED_BIN" ]]; then
  SIM_ARGS+=(--single-point-expected-bin "$SINGLE_POINT_EXPECTED_BIN")
fi
if [[ -n "$MOVING_TARGET_BEAM_ID" ]]; then
  SIM_ARGS+=(--moving-target-beam-id "$MOVING_TARGET_BEAM_ID")
fi
if [[ -n "$MOVING_TARGET_EXPECTED_BIN" ]]; then
  SIM_ARGS+=(--moving-target-expected-bin "$MOVING_TARGET_EXPECTED_BIN")
fi
if [[ -n "$MOVING_TARGET_SPEED_MPS" ]]; then
  SIM_ARGS+=(--moving-target-speed-mps "$MOVING_TARGET_SPEED_MPS")
fi
if [[ -n "$MOVING_TARGET_VELOCITY_MODE" ]]; then
  SIM_ARGS+=(--moving-target-velocity-mode "$MOVING_TARGET_VELOCITY_MODE")
fi
if [[ -n "$MOVING_TARGET_VE_MPS" ]]; then
  SIM_ARGS+=(--moving-target-ve-mps "$MOVING_TARGET_VE_MPS")
fi
if [[ -n "$MOVING_TARGET_VN_MPS" ]]; then
  SIM_ARGS+=(--moving-target-vn-mps "$MOVING_TARGET_VN_MPS")
fi
if [[ -n "$MOVING_TARGET_RADIAL_SPEED_MPS" ]]; then
  SIM_ARGS+=(--moving-target-radial-speed-mps "$MOVING_TARGET_RADIAL_SPEED_MPS")
fi
if [[ -n "$MOVING_TARGET_TANGENTIAL_SPEED_MPS" ]]; then
  SIM_ARGS+=(--moving-target-tangential-speed-mps "$MOVING_TARGET_TANGENTIAL_SPEED_MPS")
fi
if [[ -n "$MOVING_TARGET_RCS_DB" ]]; then
  SIM_ARGS+=(--moving-target-rcs-db "$MOVING_TARGET_RCS_DB")
fi
if [[ -n "$CLUTTER_AMPLITUDE_SCALE" ]]; then
  SIM_ARGS+=(--clutter-amplitude-scale "$CLUTTER_AMPLITUDE_SCALE")
fi
if [[ -n "$AREA_CLUTTER_SCATTERER_COUNT" ]]; then
  SIM_ARGS+=(--area-clutter-scatterer-count "$AREA_CLUTTER_SCATTERER_COUNT")
fi
if [[ -n "$AREA_CLUTTER_MEAN_POWER" ]]; then
  SIM_ARGS+=(--area-clutter-mean-power "$AREA_CLUTTER_MEAN_POWER")
fi
if [[ -n "$AREA_CLUTTER_TEXTURE_SIGMA" ]]; then
  SIM_ARGS+=(--area-clutter-texture-sigma "$AREA_CLUTTER_TEXTURE_SIGMA")
fi
if [[ -n "$STRONG_SCATTERER_COUNT" ]]; then
  SIM_ARGS+=(--strong-scatterer-count "$STRONG_SCATTERER_COUNT")
fi
if [[ -n "$STRONG_RCS_DB_MIN" ]]; then
  SIM_ARGS+=(--strong-rcs-db-min "$STRONG_RCS_DB_MIN")
fi
if [[ -n "$STRONG_RCS_DB_MAX" ]]; then
  SIM_ARGS+=(--strong-rcs-db-max "$STRONG_RCS_DB_MAX")
fi
if [[ -n "$LINE_SCATTERER_COUNT" ]]; then
  SIM_ARGS+=(--line-scatterer-count "$LINE_SCATTERER_COUNT")
fi
if [[ -n "$LINE_POINTS_PER_LINE" ]]; then
  SIM_ARGS+=(--line-points-per-line "$LINE_POINTS_PER_LINE")
fi
if [[ -n "$LINE_RCS_DB" ]]; then
  SIM_ARGS+=(--line-rcs-db "$LINE_RCS_DB")
fi
if [[ -n "$NOISE_POWER" ]]; then
  SIM_ARGS+=(--noise-power "$NOISE_POWER")
fi
"${SIM_ARGS[@]}" 2>&1 | tee "$OUTPUT_DIR/logs/simulate_stage2_console.log"

python3 - "$OUTPUT_DIR" "$TARGET_ENABLED" <<'PY'
import csv
import math
import sys
from pathlib import Path

out_dir = Path(sys.argv[1])
target_enabled = sys.argv[2].lower() in {"1", "true", "yes", "on"}
path = out_dir / "truth" / "target_truth_beam_summary.csv"

def to_float(v, default=math.nan):
    try:
        return float(v)
    except Exception:
        return default

def fmt(v):
    if v is None:
        return ""
    return str(v)

if not path.exists():
    print(f"[oneclick][target] summary not found: {path}")
    if target_enabled:
        print("[oneclick][ERR] target-enabled=true but target truth summary is missing.", file=sys.stderr)
        sys.exit(1)
    sys.exit(0)

with path.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))

injected = [r for r in rows if to_float(r.get("injected_sample_count"), 0.0) > 0.0]
print(f"[oneclick][target] injected rows = {len(injected)}")

for r in injected[:5]:
    fields = [
        "beam_id",
        "beam_id_0based",
        "beam_id_1based",
        "expected_range_bin",
        "visible_pulse_count",
        "injected_sample_count",
        "moving_target_speed_mps",
        "rcs_db",
        "target_ve_mps",
        "target_vn_mps",
        "target_vr_self_mps",
        "target_vt_self_mps",
        "af_motion_truth_hz",
    ]
    print("[oneclick][target] " + " ".join(f"{k}={fmt(r.get(k, ''))}" for k in fields))

if target_enabled and not injected:
    print("[oneclick][ERR] target-enabled=true but no cooperative target was injected.", file=sys.stderr)
    sys.exit(1)
PY

if [[ "$RUN_ALGORITHM" == "true" ]]; then
  echo "[oneclick] run GMTI_core"
  mkdir -p "$OUTPUT_DIR/algorithm_result"
  ./build/GMTI_core "$OUTPUT_DIR/config/temp_config_stage2_newsystem.xml" \
    > "$OUTPUT_DIR/logs/algorithm_run.log" 2>&1 || {
      echo "[oneclick][WARN] GMTI_core exited non-zero; report will include the log." >&2
    }
else
  echo "[oneclick] skip GMTI_core; pass --run-algorithm true to include imaging/cancellation/detection/tracking products"
fi

echo "[oneclick] generate report"
REPORT_PATH="$(python3 simulator/reporting/generate_simulation_report.py \
  --output-dir "$OUTPUT_DIR" \
  --stage2-config "$STAGE2_CONFIG" \
  --target-config "$TARGET_CONFIG" \
  --scene-mode "$SCENE_MODE" \
  --target-enabled "$TARGET_ENABLED")"

echo "[oneclick] report: $REPORT_PATH"
