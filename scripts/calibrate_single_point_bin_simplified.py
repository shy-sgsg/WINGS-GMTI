#!/usr/bin/env python3
"""Simplified GMTI calibration / localization checker.

目标：少传参数、自动找文件、输出更详细的诊断结果。

最常用运行方式：

    python3 calibrate_single_point_bin_simplified.py /path/to/case_or_result_dir

脚本会自动递归查找：
    runtime_config_dump.json
    detection_results.csv
    scene_truth.csv
    target_truth_pulse.csv
    target_truth_beam_summary.csv

输出目录默认：
    <输入目录>/calibration_check

说明：门限类参数已经固定在本文件顶部的 FIXED 配置中。一般测试不需要命令行传一堆参数。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

C_MPS = 299_792_458.0

# ===============================
# 固定参数区：正常测试只改这里
# ===============================
FIXED = {
    # bin 映射：expected_pc_bin = round(range_sample_float) - pc_crop_start + pc_bin_offset
    "pc_bin_offset": 0,

    # 波位匹配：truth 与 detection 允许相邻 1 个波位；若有角度字段，则优先用角度差判断。
    "max_beam_diff": 1,
    "max_theta_diff_deg": 1.0,

    # 单点脉压峰标定门限
    "single_max_range_error_m": 1000.0,
    "single_max_range_bin_error": 50.0,

    # 合作动目标定位匹配门限
    "moving_max_range_error_m": 300.0,
    "moving_max_range_bin_error": 20.0,
    "moving_max_position_error_m": 500.0,
    "moving_max_velocity_error_mps": 20.0,

    # detection_results.csv 里的 beam_id 基准：auto / zero / one
    # auto：优先根据角度判断，无法判断时按原值理解为 1-based。
    "det_beam_id_base": "auto",

    # 多 truth 时是否只选一个最强 truth 做单点标定
    "single_truth_only": True,

    # 明细输出中保留前多少个候选
    "top_candidates_in_report": 20,
}

RUNTIME_CONFIG_NAMES = [
    "runtime_config_dump.json",
    "runtime_config.json",
    "run_config.json",
]
DETECTION_NAMES = [
    "detection_results.csv",
    "detections.csv",
]
SCENE_TRUTH_NAMES = [
    "scene_truth.csv",
    "sar_scatterers_truth.csv",
]
TRUTH_PULSE_NAMES = [
    "target_truth_pulse.csv",
    "truth_targets_by_pulse.csv",
    "target_truth.csv",
    "target_truth_pulse_detail.csv",
]
TRUTH_BEAM_NAMES = [
    "target_truth_beam_summary.csv",
    "target_truth_by_beam.csv",
    "truth_targets_by_beam.csv",
    "truth_targets_by_beam_summary.csv",
]

# 调试字段：如果 detection_results.csv 有这些列，会原样带入输出，便于检查运动补偿。
DETECTION_DEBUG_FIELDS = [
    "phase_rad",
    "af_phase_hz",
    "af_total_hz",
    "af_geometry_hz",
    "af_motion_hz",
    "v_radial_mps",
    "phi_motion_rad",
    "delta_t_s",
    "motion_comp_valid",
    "motion_comp_enable",
    "old_e",
    "old_n",
    "new_e",
    "new_n",
]


def to_float(value: Any, default: float = math.nan) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except Exception:
        return default


def to_int(value: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except Exception:
        return default


def truth_bool(row: Dict[str, Any], name: str) -> bool:
    value = str(row.get(name, "")).strip().lower()
    return value in {"1", "true", "yes", "y", "on"}


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def read_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def nested_get(data: Dict[str, Any], paths: Sequence[Sequence[str]], default: Any = None) -> Any:
    for keys in paths:
        cur: Any = data
        ok = True
        for key in keys:
            if not isinstance(cur, dict) or key not in cur:
                ok = False
                break
            cur = cur[key]
        if ok:
            return cur
    return default


def first_value(row: Dict[str, Any], names: Sequence[str], default: str = "") -> str:
    for name in names:
        value = row.get(name, "")
        if value not in (None, ""):
            return str(value)
    return default


def first_float(row: Dict[str, Any], names: Sequence[str], default: float = math.nan) -> float:
    return to_float(first_value(row, names, ""), default)


def first_int(row: Dict[str, Any], names: Sequence[str], default: Optional[int] = None) -> Optional[int]:
    return to_int(first_value(row, names, ""), default)


def truth_target_id(row: Dict[str, Any]) -> str:
    return first_value(row, ["target_id", "target_name", "scatterer_id", "id"], "")


def truth_period_id(row: Dict[str, Any]) -> int:
    return first_int(row, ["period_id"], 0) or 0


def truth_beam_id(row: Dict[str, Any]) -> Optional[int]:
    b = first_int(row, ["beam_id", "beam_id_1based"], None)
    if b is None:
        b0 = first_int(row, ["beam_id_0based"], None)
        if b0 is not None:
            b = b0 + 1
    return b


def det_beam_id_1based(det: Dict[str, Any], truth: Dict[str, Any]) -> Optional[int]:
    db = first_int(det, ["beam_id"], None)
    if db is None:
        return None
    base = FIXED["det_beam_id_base"]
    if base == "zero":
        return db + 1
    if base == "one":
        return db

    # auto：如果角度接近，则优先尊重 detection 原始 beam_id；否则也按 1-based 处理。
    tt = first_float(truth, ["theta_cmd_deg", "theta_true_deg", "target_azimuth_deg", "initial_azimuth_deg"])
    dt = first_float(det, ["theta_cmd_deg", "theta_true_deg"])
    tb = truth_beam_id(truth)
    if math.isfinite(tt) and math.isfinite(dt) and abs(dt - tt) <= FIXED["max_theta_diff_deg"]:
        if tb is not None and db == tb:
            return db
        if tb is not None and db + 1 == tb:
            return db + 1
    return db


def runtime_pc_crop_start(runtime_config: Dict[str, Any]) -> Optional[int]:
    value = nested_get(
        runtime_config,
        [
            ["pulse_compression", "pc_crop_start"],
            ["pulse_compression", "range_crop_start"],
            ["data_layout", "range_crop_start"],
        ],
    )
    return to_int(value, None)


def runtime_fs_hz(runtime_config: Dict[str, Any]) -> float:
    value = nested_get(
        runtime_config,
        [
            ["waveform", "fs_hz"],
            ["waveform", "fs"],
            ["radar", "fs_hz"],
            ["radar", "fs"],
        ],
    )
    return to_float(value)


def runtime_sample_delay_s(runtime_config: Dict[str, Any]) -> float:
    # 兼容 sample_delay_us / sample_delay_s / sample_delay 三种写法。
    value_us = nested_get(
        runtime_config,
        [["waveform", "sample_delay_us"], ["data_layout", "sample_delay_us"]],
    )
    if value_us is not None:
        return to_float(value_us) * 1.0e-6

    value_s = nested_get(
        runtime_config,
        [["waveform", "sample_delay_s"], ["data_layout", "sample_delay_s"]],
    )
    if value_s is not None:
        return to_float(value_s)

    value = nested_get(
        runtime_config,
        [["waveform", "sample_delay"], ["data_layout", "sample_delay"]],
    )
    v = to_float(value)
    if not math.isfinite(v):
        return math.nan
    # 一般采样延迟如果 > 1e-3，很可能单位是 us；否则按 s。
    return v * 1.0e-6 if abs(v) > 1e-3 else v


def ensure_range_sample(row: Dict[str, Any], runtime_config: Dict[str, Any], notes: List[str]) -> Dict[str, Any]:
    out = dict(row)
    if first_value(out, ["range_sample_float"], "") != "":
        if first_value(out, ["range_sample_int"], "") == "":
            rs = first_float(out, ["range_sample_float"])
            if math.isfinite(rs):
                out["range_sample_int"] = str(int(round(rs)))
        return out

    fs_hz = runtime_fs_hz(runtime_config)
    sample_delay_s = runtime_sample_delay_s(runtime_config)
    range_m = first_float(out, ["range_m", "initial_range_m"])
    if math.isfinite(fs_hz) and math.isfinite(sample_delay_s) and math.isfinite(range_m):
        tau_abs_s = 2.0 * range_m / C_MPS
        rs = (tau_abs_s - sample_delay_s) * fs_hz
        out["range_sample_float"] = f"{rs:.12g}"
        out["range_sample_int"] = str(int(round(rs)))
        notes.append("truth 缺少 range_sample_float，已由 range_m、fs 和 sample_delay 自动计算。")
    return out


def expected_pc_bin_float(truth: Dict[str, Any], runtime_config: Dict[str, Any]) -> float:
    pc_crop_start = runtime_pc_crop_start(runtime_config)
    rs = first_float(truth, ["range_sample_float"])
    if pc_crop_start is None or not math.isfinite(rs):
        return math.nan
    return rs - pc_crop_start + float(FIXED["pc_bin_offset"])


def expected_pc_bin_round(truth: Dict[str, Any], runtime_config: Dict[str, Any]) -> Optional[int]:
    pc_crop_start = runtime_pc_crop_start(runtime_config)
    rs = first_float(truth, ["range_sample_float"])
    if pc_crop_start is None or not math.isfinite(rs):
        return None
    return int(round(rs)) - pc_crop_start + int(FIXED["pc_bin_offset"])


def geo_error_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    if not all(math.isfinite(v) for v in [lat1, lon1, lat2, lon2]):
        return math.nan
    r = 6_371_000.0
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2.0) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    return 2.0 * r * math.atan2(math.sqrt(a), math.sqrt(max(0.0, 1.0 - a)))


def csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value


def write_dict_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: Optional[List[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: csv_value(row.get(k, "")) for k in fieldnames})


def write_json(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
        f.write("\n")


@dataclass
class InputPaths:
    root_dir: Path
    runtime_config: Optional[Path]
    detections: Optional[Path]
    scene_truth: Optional[Path]
    truth_pulse: Optional[Path]
    truth_beam: Optional[Path]
    output_dir: Path


def find_first(root: Path, names: Sequence[str], output_dir: Optional[Path] = None) -> Optional[Path]:
    candidates: List[Path] = []
    names_lower = {n.lower() for n in names}
    for dirpath, dirnames, filenames in os.walk(root):
        d = Path(dirpath)
        if output_dir is not None:
            try:
                d.relative_to(output_dir)
                continue
            except ValueError:
                pass
        for filename in filenames:
            if filename.lower() in names_lower:
                candidates.append(d / filename)
    if not candidates:
        return None
    # 优先根目录浅层，其次按文件名顺序。
    candidates.sort(key=lambda p: (len(p.relative_to(root).parts), names.index(p.name) if p.name in names else 999, str(p)))
    return candidates[0]


def resolve_inputs(args: argparse.Namespace) -> InputPaths:
    root = Path(args.result_dir).resolve()
    out_dir = Path(args.output_dir).resolve() if args.output_dir else root / "calibration_check"
    runtime = Path(args.runtime_config).resolve() if args.runtime_config else find_first(root, RUNTIME_CONFIG_NAMES, out_dir)
    detections = Path(args.detections).resolve() if args.detections else find_first(root, DETECTION_NAMES, out_dir)
    scene_truth = Path(args.scene_truth).resolve() if args.scene_truth else find_first(root, SCENE_TRUTH_NAMES, out_dir)
    truth_pulse = Path(args.truth_pulse).resolve() if args.truth_pulse else find_first(root, TRUTH_PULSE_NAMES, out_dir)
    truth_beam = Path(args.truth_beam).resolve() if args.truth_beam else find_first(root, TRUTH_BEAM_NAMES, out_dir)
    return InputPaths(root, runtime, detections, scene_truth, truth_pulse, truth_beam, out_dir)


def validate_common_inputs(paths: InputPaths) -> None:
    missing = []
    if paths.runtime_config is None or not paths.runtime_config.exists():
        missing.append("runtime_config_dump.json")
    if paths.detections is None or not paths.detections.exists():
        missing.append("detection_results.csv")
    if missing:
        raise FileNotFoundError("缺少必要输入文件：" + ", ".join(missing))


# ===============================
# truth 选择
# ===============================

def moving_truth_injected(row: Dict[str, Any]) -> bool:
    if first_value(row, ["beam_injected"], "") != "":
        return truth_bool(row, "beam_injected")
    if first_value(row, ["injection_enabled"], "") != "":
        return truth_bool(row, "injection_enabled")
    if first_value(row, ["injected_sample_count"], "") != "":
        return first_float(row, ["injected_sample_count"], 0.0) > 0.0
    return False


def moving_truth_visible(row: Dict[str, Any]) -> bool:
    if first_value(row, ["visible_by_beam"], "") != "":
        return truth_bool(row, "visible_by_beam")
    if first_value(row, ["visible_pulse_count"], "") != "":
        return first_float(row, ["visible_pulse_count"], 0.0) > 0.0
    return False


def moving_truth_valid(row: Dict[str, Any]) -> bool:
    return moving_truth_injected(row) and moving_truth_visible(row)


def aggregate_pulse_truth(pulse_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    groups: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
    for row in pulse_rows:
        if not moving_truth_valid(row):
            continue
        tid = truth_target_id(row)
        period = first_value(row, ["period_id"], "0")
        beam = str(truth_beam_id(row) or "")
        groups.setdefault((tid, period, beam), []).append(row)

    numeric_fields = [
        "range_m", "range_sample_float", "range_sample_int",
        "e", "n", "target_e", "target_n",
        "lat", "lon", "target_lat", "target_lon",
        "theta_cmd_deg", "theta_true_deg", "target_azimuth_deg",
        "radial_velocity_mps", "vr_mps", "velocity_mps",
    ]
    out: List[Dict[str, Any]] = []
    for (tid, period, beam), rows in groups.items():
        merged = dict(rows[len(rows) // 2])
        merged["target_id"] = tid
        merged["period_id"] = period
        merged["beam_id"] = beam
        for field in numeric_fields:
            vals = [to_float(r.get(field)) for r in rows]
            vals = [v for v in vals if math.isfinite(v)]
            if vals:
                merged[field] = str(sum(vals) / len(vals))
        merged["visible_pulse_count"] = str(len(rows))
        out.append(merged)
    return out


def choose_scene_truths(scene_rows: List[Dict[str, Any]], runtime_config: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[str]]:
    notes = ["truth_source=scene_truth"]
    scored: List[Dict[str, Any]] = []
    for row in scene_rows:
        r = dict(row)
        if first_value(r, ["target_id"], "") == "" and first_value(r, ["scatterer_id"], "") != "":
            r["target_id"] = r.get("scatterer_id", "")
        if first_value(r, ["range_m"], "") == "":
            r["range_m"] = first_value(r, ["initial_range_m"], "")
        if first_value(r, ["theta_true_deg"], "") == "":
            r["theta_true_deg"] = first_value(r, ["initial_azimuth_deg", "target_azimuth_deg"], "")
        if truth_beam_id(r) is None and first_value(r, ["beam_id_0based"], "") != "":
            r["beam_id"] = str((first_int(r, ["beam_id_0based"], -1) or -1) + 1)
        r = ensure_range_sample(r, runtime_config, notes)

        type_score = 1.0 if str(r.get("type", "")).strip() in {"single_point", "point", "scatterer", "strong"} else 0.0
        visible_score = 1.0 if truth_bool(r, "visible_by_beam") else 0.0
        amp = first_float(r, ["amplitude", "target_amplitude"], math.nan)
        rcs = first_float(r, ["rcs_db"], math.nan)
        score = type_score * 1000.0 + visible_score * 100.0 + (amp if math.isfinite(amp) else 0.0) + (rcs if math.isfinite(rcs) else -999.0) * 0.001
        r["truth_score"] = score
        scored.append(r)

    scored.sort(key=lambda r: r.get("truth_score", -1e18), reverse=True)
    selected = scored[:1] if FIXED["single_truth_only"] and scored else scored
    if not scored:
        notes.append("scene_truth.csv 中没有可用散射点。")
    return selected, scored, notes


def choose_cooperative_truths(pulse_rows: List[Dict[str, Any]], beam_rows: List[Dict[str, Any]], runtime_config: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[str]]:
    notes = ["truth_source=cooperative_target"]

    # 优先使用逐波位 beam summary；没有则从 pulse truth 聚合。
    beam_valid = [dict(r) for r in beam_rows if moving_truth_valid(r)]
    if beam_valid:
        candidates = beam_valid
        notes.append("使用 target_truth_beam_summary 中 beam_injected/visible 记录。")
    else:
        candidates = aggregate_pulse_truth(pulse_rows)
        notes.append("beam summary 没有可用逐波位记录，已从 target_truth_pulse 聚合。")

    scored: List[Dict[str, Any]] = []
    for row in candidates:
        r = ensure_range_sample(dict(row), runtime_config, notes)
        injected = 1.0 if moving_truth_injected(r) else 0.0
        visible = 1.0 if moving_truth_visible(r) else 0.0
        injected_count = first_float(r, ["injected_sample_count"], 0.0)
        visible_count = first_float(r, ["visible_pulse_count"], 0.0)
        beam_gain = first_float(r, ["mean_beam_gain", "beam_gain"], math.nan)
        amp = first_float(r, ["target_amplitude", "amplitude"], math.nan)
        score = injected * 10000.0 + visible * 5000.0 + injected_count + visible_count + (beam_gain if math.isfinite(beam_gain) else 0.0) + (amp if math.isfinite(amp) else 0.0) * 0.001
        r["truth_score"] = score
        scored.append(r)

    scored.sort(key=lambda r: r.get("truth_score", -1e18), reverse=True)
    selected = scored[:1] if FIXED["single_truth_only"] and scored else scored
    if not scored:
        notes.append("合作目标 truth 中没有 injection/visible 有效记录。")
    return selected, scored, notes


# ===============================
# 候选计算和匹配
# ===============================

def beam_match(truth: Dict[str, Any], det: Dict[str, Any]) -> Tuple[bool, Optional[float], str]:
    tt = first_float(truth, ["theta_cmd_deg", "theta_true_deg", "target_azimuth_deg", "initial_azimuth_deg"])
    dt = first_float(det, ["theta_cmd_deg", "theta_true_deg"])
    if math.isfinite(tt) and math.isfinite(dt):
        diff = abs(dt - tt)
        return diff <= FIXED["max_theta_diff_deg"], diff, "theta_deg"

    tb = truth_beam_id(truth)
    db = det_beam_id_1based(det, truth)
    if tb is None or db is None:
        return False, None, "missing_beam_id"
    diff = abs(db - tb)
    return diff <= FIXED["max_beam_diff"], float(diff), "beam_id"


def make_candidate(truth: Dict[str, Any], det: Dict[str, Any], runtime_config: Dict[str, Any], *, mode: str) -> Dict[str, Any]:
    tp = truth_period_id(truth)
    dp = first_int(det, ["period_id"], None)
    period_ok = (tp <= 0 or dp is None or dp == tp)

    beam_ok, beam_metric, beam_mode = beam_match(truth, det)

    truth_range = first_float(truth, ["range_m", "initial_range_m"])
    det_range = first_float(det, ["range_m"])
    range_error = det_range - truth_range if math.isfinite(det_range) and math.isfinite(truth_range) else math.nan

    exp_bin_float = expected_pc_bin_float(truth, runtime_config)
    exp_bin_round = expected_pc_bin_round(truth, runtime_config)
    det_bin = first_int(det, ["range_bin", "col"], None)
    bin_error_float = det_bin - exp_bin_float if det_bin is not None and math.isfinite(exp_bin_float) else math.nan
    bin_error_round = det_bin - exp_bin_round if det_bin is not None and exp_bin_round is not None else None

    te = first_float(truth, ["target_e", "e", "x_m", "target_x"])
    tn = first_float(truth, ["target_n", "n", "y_m", "target_y"])
    de = first_float(det, ["e", "new_e"])
    dn = first_float(det, ["n", "new_n"])
    pos_error = math.hypot(de - te, dn - tn) if all(math.isfinite(v) for v in [te, tn, de, dn]) else math.nan

    tlat = first_float(truth, ["lat", "target_lat"])
    tlon = first_float(truth, ["lon", "target_lon"])
    dlat = first_float(det, ["lat"])
    dlon = first_float(det, ["lon"])
    geo_err = geo_error_m(tlat, tlon, dlat, dlon)

    tt = first_float(truth, ["theta_true_deg", "theta_cmd_deg", "target_azimuth_deg", "initial_azimuth_deg"])
    dt = first_float(det, ["theta_true_deg", "theta_cmd_deg"])
    angle_error = dt - tt if math.isfinite(dt) and math.isfinite(tt) else math.nan

    tv = first_float(truth, ["radial_velocity_mps", "vr_mps", "velocity_mps"])
    dv = first_float(det, ["radial_velocity_mps", "v_radial_mps"])
    vel_error = dv - tv if math.isfinite(dv) and math.isfinite(tv) else math.nan

    reasons: List[str] = []
    if not period_ok:
        reasons.append("period_mismatch")
    if not beam_ok:
        reasons.append(f"beam_mismatch_by_{beam_mode}")
    if not math.isfinite(det_range):
        reasons.append("missing_det_range_m")
    if not math.isfinite(truth_range):
        reasons.append("missing_truth_range_m")

    if mode == "single":
        max_range_error = FIXED["single_max_range_error_m"]
        max_bin_error = FIXED["single_max_range_bin_error"]
    else:
        max_range_error = FIXED["moving_max_range_error_m"]
        max_bin_error = FIXED["moving_max_range_bin_error"]

    if math.isfinite(range_error) and abs(range_error) > max_range_error:
        reasons.append("range_gate_failed")
    if math.isfinite(bin_error_float) and abs(bin_error_float) > max_bin_error:
        reasons.append("range_bin_gate_failed")
    if mode == "moving" and math.isfinite(pos_error) and pos_error > FIXED["moving_max_position_error_m"]:
        reasons.append("position_gate_failed")
    if mode == "moving" and math.isfinite(vel_error) and abs(vel_error) > FIXED["moving_max_velocity_error_mps"]:
        reasons.append("velocity_gate_failed")

    if math.isfinite(bin_error_float):
        cost = abs(bin_error_float) * 100.0 + (abs(range_error) if math.isfinite(range_error) else 0.0)
    elif math.isfinite(range_error):
        cost = abs(range_error) + (beam_metric if beam_metric is not None else 1000.0) * 100.0
    elif math.isfinite(pos_error):
        cost = pos_error
    else:
        cost = 1.0e18
        reasons.append("missing_all_cost_fields")

    out = {
        "target_id": truth_target_id(truth),
        "truth_period_id": tp,
        "truth_beam_id": truth_beam_id(truth),
        "det_id": first_value(det, ["det_id"], ""),
        "det_period_id": dp,
        "det_beam_id_raw": first_int(det, ["beam_id"], None),
        "det_beam_id_1based": det_beam_id_1based(det, truth),
        "beam_match_mode": beam_mode,
        "beam_metric": beam_metric,
        "truth_range_m": truth_range,
        "det_range_m": det_range,
        "range_error_m": range_error,
        "truth_range_sample_float": first_float(truth, ["range_sample_float"]),
        "truth_range_sample_int": first_int(truth, ["range_sample_int"], None),
        "pc_crop_start": runtime_pc_crop_start(runtime_config),
        "pc_bin_offset_used": FIXED["pc_bin_offset"],
        "truth_pc_bin_expected_float": exp_bin_float,
        "truth_pc_bin_expected_round": exp_bin_round,
        "det_range_bin": det_bin,
        "range_bin_error_float": bin_error_float,
        "range_bin_error_round": bin_error_round,
        "truth_e": te,
        "truth_n": tn,
        "det_e": de,
        "det_n": dn,
        "pos_error_m": pos_error,
        "truth_lat": tlat,
        "truth_lon": tlon,
        "det_lat": dlat,
        "det_lon": dlon,
        "geo_error_m": geo_err,
        "truth_theta_deg": tt,
        "det_theta_deg": dt,
        "angle_error_deg": angle_error,
        "truth_radial_velocity_mps": tv,
        "det_radial_velocity_mps": dv,
        "velocity_error_mps": vel_error,
        "total_cost": cost,
        "gate_passed": not reasons,
        "gate_fail_reason": ";".join(reasons),
    }

    # 附带检测中的运动补偿调试字段。
    for field in DETECTION_DEBUG_FIELDS:
        if field in det:
            out[field] = det.get(field, "")
    return out


def summarize_matches(matches: List[Dict[str, Any]], detections: List[Dict[str, Any]], truths: List[Dict[str, Any]], notes: List[str]) -> Dict[str, Any]:
    matched = [m for m in matches if m.get("match_status") == "matched"]

    def mean_abs(field: str) -> Optional[float]:
        vals = [abs(to_float(m.get(field))) for m in matched]
        vals = [v for v in vals if math.isfinite(v)]
        return sum(vals) / len(vals) if vals else None

    return {
        "status": "passed" if len(matched) == len(truths) and matched else ("partial" if matched else "failed"),
        "truth_count": len(truths),
        "detection_count": len(detections),
        "matched_count": len(matched),
        "unmatched_truth_count": len([m for m in matches if m.get("match_status") == "unmatched_truth"]),
        "mean_abs_range_error_m": mean_abs("range_error_m"),
        "mean_abs_range_bin_error_float": mean_abs("range_bin_error_float"),
        "mean_abs_pos_error_m": mean_abs("pos_error_m"),
        "mean_abs_velocity_error_mps": mean_abs("velocity_error_mps"),
        "fixed_parameters": FIXED,
        "notes": notes,
    }


# ===============================
# 单点 bin 标定
# ===============================

def run_single_point(paths: InputPaths) -> Dict[str, Any]:
    validate_common_inputs(paths)
    runtime = read_json(paths.runtime_config)  # type: ignore[arg-type]
    detections = read_csv(paths.detections)  # type: ignore[arg-type]

    if paths.scene_truth and paths.scene_truth.exists():
        truths, truth_candidates, notes = choose_scene_truths(read_csv(paths.scene_truth), runtime)
    elif paths.truth_pulse and paths.truth_beam and paths.truth_pulse.exists() and paths.truth_beam.exists():
        truths, truth_candidates, notes = choose_cooperative_truths(read_csv(paths.truth_pulse), read_csv(paths.truth_beam), runtime)
    else:
        raise FileNotFoundError("单点标定需要 scene_truth.csv，或 target_truth_pulse.csv + target_truth_beam_summary.csv。")

    out_dir = paths.output_dir / "single_point_bin"
    out_dir.mkdir(parents=True, exist_ok=True)

    truth_candidate_rows = []
    for r in truth_candidates:
        truth_candidate_rows.append({
            "target_id": truth_target_id(r),
            "period_id": truth_period_id(r),
            "beam_id": truth_beam_id(r),
            "range_m": first_float(r, ["range_m", "initial_range_m"]),
            "range_sample_float": first_float(r, ["range_sample_float"]),
            "range_sample_int": first_int(r, ["range_sample_int"], None),
            "theta_deg": first_float(r, ["theta_true_deg", "theta_cmd_deg", "target_azimuth_deg", "initial_azimuth_deg"]),
            "truth_score": r.get("truth_score", ""),
        })
    write_dict_csv(out_dir / "truth_candidates.csv", truth_candidate_rows)

    if not truths:
        summary = {
            "calibration_status": "failed",
            "reason": "no_truth_selected",
            "fixed_parameters": FIXED,
            "input_files": input_files_dict(paths),
            "notes": notes,
        }
        write_json(out_dir / "single_point_bin_calibration.json", summary)
        return summary

    truth = truths[0]
    candidates = [make_candidate(truth, det, runtime, mode="single") for det in detections]
    candidates.sort(key=lambda r: r.get("total_cost", 1.0e18))
    best = next((c for c in candidates if c.get("gate_passed")), candidates[0] if candidates else None)

    status = "success" if best and best.get("gate_passed") else "failed"
    result_row = dict(best or {})
    result_row["calibration_status"] = status
    result_row["case_id"] = nested_get(runtime, [["run_info", "case_id"]], "")
    result_row["run_id"] = nested_get(runtime, [["run_info", "run_id"]], "")
    result_row["result_id"] = nested_get(runtime, [["run_info", "result_id"]], "")
    result_row["truth_source"] = "scene_truth" if paths.scene_truth else "cooperative_target"

    write_dict_csv(out_dir / "single_point_candidate_detections.csv", candidates)
    write_dict_csv(out_dir / "single_point_bin_calibration.csv", [result_row])

    summary = {
        "calibration_status": status,
        "case_id": result_row.get("case_id", ""),
        "run_id": result_row.get("run_id", ""),
        "result_id": result_row.get("result_id", ""),
        "truth_source": result_row.get("truth_source", ""),
        "selected_truth": truth_candidate_rows[0] if truth_candidate_rows else {},
        "selected_detection": result_row,
        "candidate_count": len(candidates),
        "gate_passed_candidate_count": len([c for c in candidates if c.get("gate_passed")]),
        "fixed_parameters": FIXED,
        "input_files": input_files_dict(paths),
        "notes": notes,
    }
    write_json(out_dir / "single_point_bin_calibration.json", summary)
    write_single_report(out_dir / "single_point_calibration_report.md", summary, candidates)
    return summary


# ===============================
# 合作动目标定位检查
# ===============================

def run_moving(paths: InputPaths) -> Dict[str, Any]:
    validate_common_inputs(paths)
    if not (paths.truth_pulse and paths.truth_beam and paths.truth_pulse.exists() and paths.truth_beam.exists()):
        raise FileNotFoundError("合作动目标定位检查需要 target_truth_pulse.csv + target_truth_beam_summary.csv。")

    runtime = read_json(paths.runtime_config)  # type: ignore[arg-type]
    detections = read_csv(paths.detections)  # type: ignore[arg-type]
    pulse_rows = read_csv(paths.truth_pulse)
    beam_rows = read_csv(paths.truth_beam)
    truths, truth_candidates, notes = choose_cooperative_truths(pulse_rows, beam_rows, runtime)
    if not FIXED["single_truth_only"]:
        truths = truth_candidates

    out_dir = paths.output_dir / "moving_target_localization"
    out_dir.mkdir(parents=True, exist_ok=True)

    candidate_rows: List[Dict[str, Any]] = []
    best_by_truth: Dict[int, Optional[Dict[str, Any]]] = {}
    for ti, truth in enumerate(truths):
        local_candidates = [make_candidate(truth, det, runtime, mode="moving") for det in detections]
        for c in local_candidates:
            c["truth_index"] = ti
        candidate_rows.extend(local_candidates)
        valid = [c for c in local_candidates if c.get("gate_passed")]
        best_by_truth[ti] = min(valid, key=lambda c: c.get("total_cost", 1.0e18)) if valid else None

    # 简单一对一：如果多个 truth 选中同一 detection，只保留 cost 最低的。
    det_owner: Dict[Tuple[Any, Any, Any, Any], Tuple[int, Dict[str, Any]]] = {}
    for ti, cand in best_by_truth.items():
        if cand is None:
            continue
        det_key = (cand.get("det_id"), cand.get("det_period_id"), cand.get("det_beam_id_raw"), cand.get("det_range_bin"))
        old = det_owner.get(det_key)
        if old is None or cand.get("total_cost", 1.0e18) < old[1].get("total_cost", 1.0e18):
            det_owner[det_key] = (ti, cand)

    result_rows: List[Dict[str, Any]] = []
    for ti, truth in enumerate(truths):
        cand = best_by_truth.get(ti)
        if cand is None:
            row = make_candidate(truth, {}, runtime, mode="moving")
            row["match_status"] = "unmatched_truth"
            row["reason"] = "no_gate_passed_candidate"
            result_rows.append(row)
            continue
        det_key = (cand.get("det_id"), cand.get("det_period_id"), cand.get("det_beam_id_raw"), cand.get("det_range_bin"))
        row = dict(cand)
        if det_owner.get(det_key, (None, None))[0] == ti:
            row["match_status"] = "matched"
            row["reason"] = ""
        else:
            row["match_status"] = "duplicate_det"
            row["reason"] = "same_detection_selected_by_lower_cost_truth"
        result_rows.append(row)

    case_id = nested_get(runtime, [["run_info", "case_id"]], "")
    run_id = nested_get(runtime, [["run_info", "run_id"]], "")
    result_id = nested_get(runtime, [["run_info", "result_id"]], "")
    for row in result_rows:
        row["case_id"] = case_id
        row["run_id"] = run_id
        row["result_id"] = result_id

    candidate_rows.sort(key=lambda r: r.get("total_cost", 1.0e18))
    write_dict_csv(out_dir / "moving_target_matches.csv", result_rows)
    write_dict_csv(out_dir / "moving_target_candidate_matches.csv", candidate_rows)

    summary = summarize_matches(result_rows, detections, truths, notes)
    summary.update({
        "case_id": case_id,
        "run_id": run_id,
        "result_id": result_id,
        "candidate_count": len(candidate_rows),
        "input_files": input_files_dict(paths),
    })
    write_json(out_dir / "moving_target_localization_summary.json", summary)
    write_moving_report(out_dir / "moving_target_localization_report.md", summary, result_rows, candidate_rows)
    return summary


# ===============================
# 报告与 manifest
# ===============================

def input_files_dict(paths: InputPaths) -> Dict[str, Optional[str]]:
    return {
        "root_dir": str(paths.root_dir),
        "runtime_config": str(paths.runtime_config) if paths.runtime_config else None,
        "detections": str(paths.detections) if paths.detections else None,
        "scene_truth": str(paths.scene_truth) if paths.scene_truth else None,
        "truth_pulse": str(paths.truth_pulse) if paths.truth_pulse else None,
        "truth_beam": str(paths.truth_beam) if paths.truth_beam else None,
        "output_dir": str(paths.output_dir),
    }


def fmt(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if not math.isfinite(v):
            return ""
        return f"{v:.6g}"
    return str(v)


def write_single_report(path: Path, summary: Dict[str, Any], candidates: List[Dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("# 单点目标 range-bin 标定报告\n\n")
        f.write("## 1. 结论\n\n")
        f.write(f"- 标定状态：`{summary.get('calibration_status')}`\n")
        f.write(f"- case_id：`{summary.get('case_id', '')}`\n")
        f.write(f"- truth_source：`{summary.get('truth_source', '')}`\n")
        det = summary.get("selected_detection", {})
        f.write(f"- truth range_sample_float：`{fmt(det.get('truth_range_sample_float'))}`\n")
        f.write(f"- pc_crop_start：`{fmt(det.get('pc_crop_start'))}`\n")
        f.write(f"- expected_pc_bin_float：`{fmt(det.get('truth_pc_bin_expected_float'))}`\n")
        f.write(f"- expected_pc_bin_round：`{fmt(det.get('truth_pc_bin_expected_round'))}`\n")
        f.write(f"- detection range_bin：`{fmt(det.get('det_range_bin'))}`\n")
        f.write(f"- range_bin_error_float：`{fmt(det.get('range_bin_error_float'))}`\n")
        f.write(f"- range_error_m：`{fmt(det.get('range_error_m'))}`\n")
        f.write(f"- gate_fail_reason：`{det.get('gate_fail_reason', '')}`\n\n")

        f.write("## 2. 固定参数\n\n")
        for k, v in FIXED.items():
            f.write(f"- {k}: `{v}`\n")
        f.write("\n")

        f.write("## 3. 输入文件\n\n")
        for k, v in summary.get("input_files", {}).items():
            f.write(f"- {k}: `{v}`\n")
        f.write("\n")

        f.write("## 4. 最接近的候选检测点\n\n")
        f.write("| rank | det_id | det_period | det_beam | det_bin | det_range_m | bin_err | range_err_m | cost | gate | reason |\n")
        f.write("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n")
        for i, row in enumerate(candidates[: int(FIXED["top_candidates_in_report"])], 1):
            f.write(
                f"| {i} | {fmt(row.get('det_id'))} | {fmt(row.get('det_period_id'))} | "
                f"{fmt(row.get('det_beam_id_raw'))} | {fmt(row.get('det_range_bin'))} | "
                f"{fmt(row.get('det_range_m'))} | {fmt(row.get('range_bin_error_float'))} | "
                f"{fmt(row.get('range_error_m'))} | {fmt(row.get('total_cost'))} | "
                f"{row.get('gate_passed')} | {row.get('gate_fail_reason', '')} |\n"
            )
        f.write("\n")

        f.write("## 5. 说明\n\n")
        for note in summary.get("notes", []):
            f.write(f"- {note}\n")
        if summary.get("calibration_status") != "success":
            f.write("\n标定失败时，不建议继续做严肃的检测率、定位误差或 SCNR 目标 ROI 统计，应先查看 `single_point_candidate_detections.csv` 的失败原因。\n")


def write_moving_report(path: Path, summary: Dict[str, Any], matches: List[Dict[str, Any]], candidates: List[Dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("# 合作动目标定位检查报告\n\n")
        f.write("## 1. 汇总结论\n\n")
        for key in [
            "status", "truth_count", "detection_count", "candidate_count", "matched_count",
            "unmatched_truth_count", "mean_abs_range_error_m", "mean_abs_range_bin_error_float",
            "mean_abs_pos_error_m", "mean_abs_velocity_error_mps",
        ]:
            f.write(f"- {key}: `{fmt(summary.get(key))}`\n")
        f.write("\n")

        f.write("## 2. 匹配结果\n\n")
        f.write("| target | truth_beam | det_id | det_beam | bin_err | range_err_m | pos_err_m | vel_err_mps | status | reason |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---|---|\n")
        for row in matches:
            f.write(
                f"| {row.get('target_id', '')} | {fmt(row.get('truth_beam_id'))} | {fmt(row.get('det_id'))} | "
                f"{fmt(row.get('det_beam_id_raw'))} | {fmt(row.get('range_bin_error_float'))} | "
                f"{fmt(row.get('range_error_m'))} | {fmt(row.get('pos_error_m'))} | "
                f"{fmt(row.get('velocity_error_mps'))} | {row.get('match_status', '')} | {row.get('reason', row.get('gate_fail_reason', ''))} |\n"
            )
        f.write("\n")

        f.write("## 3. 固定参数\n\n")
        for k, v in FIXED.items():
            f.write(f"- {k}: `{v}`\n")
        f.write("\n")

        f.write("## 4. 输入文件\n\n")
        for k, v in summary.get("input_files", {}).items():
            f.write(f"- {k}: `{v}`\n")
        f.write("\n")

        f.write("## 5. 最接近的候选检测点\n\n")
        f.write("| rank | target | det_id | det_beam | bin_err | range_err_m | pos_err_m | vel_err_mps | cost | gate | reason |\n")
        f.write("|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---|\n")
        for i, row in enumerate(candidates[: int(FIXED["top_candidates_in_report"])], 1):
            f.write(
                f"| {i} | {row.get('target_id', '')} | {fmt(row.get('det_id'))} | {fmt(row.get('det_beam_id_raw'))} | "
                f"{fmt(row.get('range_bin_error_float'))} | {fmt(row.get('range_error_m'))} | "
                f"{fmt(row.get('pos_error_m'))} | {fmt(row.get('velocity_error_mps'))} | "
                f"{fmt(row.get('total_cost'))} | {row.get('gate_passed')} | {row.get('gate_fail_reason', '')} |\n"
            )


def write_manifest(paths: InputPaths, summaries: Dict[str, Dict[str, Any]]) -> None:
    manifest = {
        "tool": "calibrate_single_point_bin_simplified.py",
        "input_files": input_files_dict(paths),
        "fixed_parameters": FIXED,
        "summaries": summaries,
        "generated_dirs": {
            "single_point_bin": str(paths.output_dir / "single_point_bin"),
            "moving_target_localization": str(paths.output_dir / "moving_target_localization"),
        },
    }
    write_json(paths.output_dir / "calibration_manifest.json", manifest)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="少参数版 GMTI 单点 range-bin 标定 / 合作动目标定位检查脚本。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("result_dir", help="结果目录或 case 根目录。脚本会递归查找输入文件。")
    parser.add_argument("--mode", choices=["auto", "single", "moving", "both"], default="auto", help="auto 会根据找到的 truth 文件自动决定运行内容。")
    parser.add_argument("--output-dir", help="输出目录。默认 <result_dir>/calibration_check。")

    # 可选覆盖项：一般不用传，只有自动查找失败或同目录多份文件时才需要。
    parser.add_argument("--runtime-config", help="手动指定 runtime_config_dump.json。")
    parser.add_argument("--detections", help="手动指定 detection_results.csv。")
    parser.add_argument("--scene-truth", help="手动指定 scene_truth.csv。")
    parser.add_argument("--truth-pulse", help="手动指定 target_truth_pulse.csv。")
    parser.add_argument("--truth-beam", help="手动指定 target_truth_beam_summary.csv。")
    args = parser.parse_args()

    try:
        paths = resolve_inputs(args)
        validate_common_inputs(paths)
        paths.output_dir.mkdir(parents=True, exist_ok=True)

        summaries: Dict[str, Dict[str, Any]] = {}

        mode = args.mode
        if mode == "auto":
            has_scene = paths.scene_truth is not None and paths.scene_truth.exists()
            has_moving = (
                paths.truth_pulse is not None and paths.truth_pulse.exists()
                and paths.truth_beam is not None and paths.truth_beam.exists()
            )
            if has_scene and has_moving:
                mode = "both"
            elif has_scene:
                mode = "single"
            elif has_moving:
                mode = "moving"
            else:
                raise FileNotFoundError("未找到 scene_truth.csv，也未找到 target_truth_pulse.csv + target_truth_beam_summary.csv。")

        if mode in {"single", "both"}:
            print("[calibration] running single_point_bin ...")
            summaries["single_point_bin"] = run_single_point(paths)
            print("[calibration] single_point_bin status=", summaries["single_point_bin"].get("calibration_status"))

        if mode in {"moving", "both"}:
            print("[calibration] running moving_target_localization ...")
            summaries["moving_target_localization"] = run_moving(paths)
            print("[calibration] moving_target_localization status=", summaries["moving_target_localization"].get("status"))

        write_manifest(paths, summaries)
        print(f"[calibration] outputs: {paths.output_dir}")
        return 0
    except Exception as exc:
        print(f"[calibration][ERR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
