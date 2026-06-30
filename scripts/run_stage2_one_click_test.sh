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

usage() {
  cat <<'EOF'
Usage:
  scripts/run_stage2_one_click_test.sh [options]

Options:
  --stage2-config PATH       Default: stage2_config.json
  --target-config PATH       Default: targets.json
  --output-dir DIR           Default: outputs/stage2_oneclick
  --scene-mode MODE          full | clutter_only | noise_only | point_target_only | area_clutter_only
  --target-enabled true|false
  --period-count N           Default: 1
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
./build/simulate_stage2_statistical \
  --stage2-config "$STAGE2_CONFIG" \
  --target-config "$TARGET_CONFIG" \
  --output-dir "$OUTPUT_DIR" \
  --scene-mode "$SCENE_MODE" \
  --target-enabled "$TARGET_ENABLED" \
  --period-count "$PERIOD_COUNT" \
  --validate true \
  2>&1 | tee "$OUTPUT_DIR/logs/simulate_stage2_console.log"

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

