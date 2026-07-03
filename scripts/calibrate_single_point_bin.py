#!/usr/bin/env python3
"""Calibrate single-point truth range_sample_float against detection range_bin."""

import argparse
import csv
import json
import math
import os
import sys


def to_float(value, default=math.nan):
    try:
        if value is None or value == "":
            return default
        return float(value)
    except Exception:
        return default


def to_int(value, default=None):
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except Exception:
        return default


def truth_bool(row, name):
    value = str(row.get(name, "")).strip().lower()
    return value in {"1", "true", "yes", "y"}


def read_csv(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def read_runtime_config(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def nested_get(data, keys, default=None):
    cur = data
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def runtime_pc_crop_start(runtime_config):
    pc_crop_start = nested_get(runtime_config, ["pulse_compression", "pc_crop_start"])
    if pc_crop_start is None:
        pc_crop_start = nested_get(runtime_config, ["pulse_compression", "range_crop_start"])
    return to_int(pc_crop_start, None)


def runtime_sample_params(runtime_config):
    fs_hz = nested_get(runtime_config, ["waveform", "fs_hz"])
    sample_delay_us = nested_get(runtime_config, ["waveform", "sample_delay_us"])
    if sample_delay_us is None:
        sample_delay_us = nested_get(runtime_config, ["data_layout", "sample_delay_us"])
    return to_float(fs_hz), to_float(sample_delay_us)


def choose_truth(pulse_rows, beam_rows):
    candidates = [
        r for r in pulse_rows
        if truth_bool(r, "injection_enabled") and truth_bool(r, "visible_by_beam")
    ]
    notes = []
    if not candidates:
        candidates = [r for r in pulse_rows if truth_bool(r, "injection_enabled")]
        notes.append("未找到 visible_by_beam=true 的注入记录，退化为 injection_enabled=true。")
    if not candidates:
        return None, "truth 中没有 injection_enabled=true 的记录，无法做单点标定。", notes

    beam_gain_by_key = {}
    for row in beam_rows:
        key = (row.get("period_id", ""), row.get("beam_id", ""), row.get("target_id", ""))
        beam_gain_by_key[key] = to_float(row.get("mean_beam_gain"), math.nan)

    def score(row):
        row_beam_id = row.get("beam_id", row.get("beam_id_1based", ""))
        if row_beam_id == "":
            b0 = to_int(row.get("beam_id_0based"), None)
            row_beam_id = "" if b0 is None else str(b0 + 1)
        key = (
            row.get("period_id", ""),
            row_beam_id,
            row.get("target_id", ""),
        )
        beam_gain = beam_gain_by_key.get(key, to_float(row.get("beam_gain"), math.nan))
        amp = to_float(row.get("target_amplitude"), math.nan)
        return (
            beam_gain if math.isfinite(beam_gain) else -1.0,
            amp if math.isfinite(amp) else -1.0,
        )

    return max(candidates, key=score), "", notes


def choose_scene_truth(scene_rows, runtime_config):
    notes = ["truth_source=scene_truth"]
    candidates = [
        r for r in scene_rows
        if str(r.get("type", "")).strip() in {"single_point", "point", "scatterer", "strong"}
    ]
    if not candidates:
        candidates = list(scene_rows)
        notes.append("scene_truth 中未找到 single_point/point/scatterer/strong 类型，退化为全部散射点。")
    if not candidates:
        return None, "scene_truth.csv 中没有散射点记录，无法做单点标定。", notes

    def score(row):
        visible = 1.0 if truth_bool(row, "visible_by_beam") else 0.0
        amp = to_float(row.get("amplitude"), math.nan)
        rcs = to_float(row.get("rcs_db"), math.nan)
        return (
            visible,
            amp if math.isfinite(amp) else -1.0,
            rcs if math.isfinite(rcs) else -999.0,
        )

    truth = max(candidates, key=score)
    truth = dict(truth)
    if "target_id" not in truth and "scatterer_id" in truth:
        truth["target_id"] = truth.get("scatterer_id")
    if "range_m" not in truth or truth.get("range_m", "") == "":
        truth["range_m"] = truth.get("initial_range_m", "")
    if "theta_true_deg" not in truth or truth.get("theta_true_deg", "") == "":
        truth["theta_true_deg"] = truth.get("initial_azimuth_deg", truth.get("target_azimuth_deg", ""))
    if "beam_id" not in truth or truth.get("beam_id", "") == "":
        if truth.get("beam_id_1based", "") != "":
            truth["beam_id"] = truth.get("beam_id_1based")
        elif truth.get("beam_id_0based", "") != "":
            b0 = to_int(truth.get("beam_id_0based"), None)
            if b0 is not None:
                truth["beam_id"] = str(b0 + 1)

    if truth.get("range_sample_float", "") == "":
        fs_hz, sample_delay_us = runtime_sample_params(runtime_config)
        range_m = to_float(truth.get("range_m"))
        if math.isfinite(fs_hz) and math.isfinite(sample_delay_us) and math.isfinite(range_m):
            tau_abs_sec = 2.0 * range_m / 299792458.0
            tau_rel_sec = tau_abs_sec - sample_delay_us * 1.0e-6
            rs = tau_rel_sec * fs_hz
            truth["range_sample_float"] = str(rs)
            truth["range_sample_int"] = str(int(round(rs)))
            notes.append("scene_truth 缺少 range_sample_float，已使用 runtime_config 的 fs_hz/sample_delay_us 和 range_m 兜底计算。")
        else:
            return None, "scene_truth.csv 缺少 range_sample_float，且 runtime_config 不足以兜底计算。", notes
    elif truth.get("range_sample_int", "") == "":
        rs = to_float(truth.get("range_sample_float"))
        if math.isfinite(rs):
            truth["range_sample_int"] = str(int(round(rs)))

    if truth.get("period_id", "") == "":
        truth["period_id"] = "0"
    return truth, "", notes


def choose_detection(detections, truth, max_range_error_m):
    truth_period = to_int(truth.get("period_id"), -1)
    truth_beam = truth_beam_id(truth)
    truth_range = to_float(truth.get("range_m"))

    candidates = []
    for row in detections:
        det_range = to_float(row.get("range_m"))
        if not math.isfinite(det_range):
            continue
        period = to_int(row.get("period_id"), -999999)
        beam = to_int(row.get("beam_id"), -999999)
        # Many current truth files use period_id=0 for a whole scan while detections use
        # beam/period indices 1..61. Treat truth period 0 as scan-global.
        period_ok = (truth_period <= 0 or period == truth_period or period == -999999)
        beam_ok = (
            truth_beam < 0 or beam == truth_beam or
            (beam != -999999 and abs(beam - truth_beam) <= 1)
        )
        if period_ok and beam_ok:
            candidates.append(row)

    if not candidates:
        return None, "未找到同周期、同波位或相邻波位的检测点。"
    if not candidates:
        return None, "未找到可用于标定的检测点。"

    best = min(candidates, key=lambda r: abs(to_float(r.get("range_m")) - truth_range))
    best_range = to_float(best.get("range_m"))
    if (
        math.isfinite(best_range) and
        math.isfinite(truth_range) and
        max_range_error_m is not None and
        abs(best_range - truth_range) > max_range_error_m
    ):
        return None, (
            "候选检测点与 truth 距离差过大，不能用于单点 bin 标定："
            f"range_error_m={best_range - truth_range:.6g}, "
            f"max_range_error_m={max_range_error_m:.6g}。"
        )
    return best, ""


def make_result(truth, detection, runtime_config, notes):
    pc_crop_start = runtime_pc_crop_start(runtime_config)

    truth_sample = to_float(truth.get("range_sample_float"))
    truth_sample_int = to_int(truth.get("range_sample_int"), None)
    if truth_sample_int is None and math.isfinite(truth_sample):
        truth_sample_int = int(round(truth_sample))
    truth_range = to_float(truth.get("range_m"))
    det_range_bin = to_int(detection.get("range_bin"), None)
    det_range = to_float(detection.get("range_m"))

    mapping = {
        "range_bin_type": "unknown",
        "pc_crop_start": pc_crop_start,
        "expected_pc_bin_without_offset": None,
        "pc_bin_offset": None,
        "bin_error": None,
        "range_error_m": det_range - truth_range if math.isfinite(det_range) and math.isfinite(truth_range) else None,
        "direct_bin_compare_valid": False,
    }

    if det_range_bin is not None and truth_sample_int is not None:
        if pc_crop_start is not None:
            expected = truth_sample_int - pc_crop_start
            mapping["range_bin_type"] = "pc_crop_bin"
            mapping["expected_pc_bin_without_offset"] = expected
            mapping["pc_bin_offset"] = det_range_bin - expected
            mapping["bin_error"] = det_range_bin - expected
            mapping["direct_bin_compare_valid"] = True
        else:
            mapping["range_bin_type"] = "unknown_or_raw_ddc_bin"
            mapping["expected_pc_bin_without_offset"] = truth_sample_int
            mapping["pc_bin_offset"] = det_range_bin - truth_sample_int
            mapping["bin_error"] = det_range_bin - truth_sample

    return {
        "calibration_status": "success",
        "case_id": nested_get(runtime_config, ["run_info", "case_id"], ""),
        "truth_source": "scene_truth" if any("truth_source=scene_truth" in n for n in notes) else "cooperative_target",
        "truth": {
            "target_id": to_int(truth.get("target_id"), None),
            "period_id": to_int(truth.get("period_id"), None),
            "beam_id": truth_beam_id(truth),
            "range_sample_float": truth_sample,
            "range_sample_int": truth_sample_int,
            "range_m": truth_range,
            "theta_true_deg": to_float(truth.get("theta_true_deg")),
            "radial_velocity_mps": to_float(truth.get("radial_velocity_mps")),
        },
        "detection": {
            "det_id": to_int(detection.get("det_id"), None),
            "period_id": to_int(detection.get("period_id"), None),
            "beam_id": to_int(detection.get("beam_id"), None),
            "range_bin": det_range_bin,
            "range_m": det_range,
            "row": to_int(detection.get("row"), None),
            "col": to_int(detection.get("col"), None),
            "theta_true_deg": to_float(detection.get("theta_true_deg")),
            "amplitude": to_float(detection.get("amplitude")),
        },
        "mapping": mapping,
        "notes": " ".join(notes),
    }


def failure_result(message, runtime_config=None, notes=None):
    notes = notes or []
    return {
        "calibration_status": "failed",
        "case_id": nested_get(runtime_config or {}, ["run_info", "case_id"], ""),
        "truth_source": "scene_truth" if any("truth_source=scene_truth" in n for n in notes) else "cooperative_target",
        "truth": {},
        "detection": {},
        "mapping": {
            "range_bin_type": "unknown",
            "pc_crop_start": runtime_pc_crop_start(runtime_config or {}),
            "expected_pc_bin_without_offset": None,
            "pc_bin_offset": None,
            "bin_error": None,
            "range_error_m": None,
            "direct_bin_compare_valid": False,
        },
        "notes": " ".join(notes + [message]),
    }


def attach_failure_truth_mapping(result):
    truth = result.get("truth", {})
    mapping = result.setdefault("mapping", {})
    truth_sample_int = truth.get("range_sample_int")
    pc_crop_start = mapping.get("pc_crop_start")
    if truth_sample_int is None:
        truth_sample = truth.get("range_sample_float")
        if isinstance(truth_sample, (int, float)) and math.isfinite(truth_sample):
            truth_sample_int = int(round(truth_sample))
            truth["range_sample_int"] = truth_sample_int
    if truth_sample_int is not None and pc_crop_start is not None:
        mapping["range_bin_type"] = "pc_crop_bin"
        mapping["expected_pc_bin_without_offset"] = truth_sample_int - pc_crop_start


def write_outputs(result, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    json_path = os.path.join(output_dir, "single_point_bin_calibration.json")
    csv_path = os.path.join(output_dir, "single_point_bin_calibration.csv")
    report_path = os.path.join(output_dir, "single_point_calibration_report.md")

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
        f.write("\n")

    row = {
        "calibration_status": result.get("calibration_status", ""),
        "case_id": result.get("case_id", ""),
        "truth_source": result.get("truth_source", ""),
        "truth_target_id": result.get("truth", {}).get("target_id"),
        "truth_period_id": result.get("truth", {}).get("period_id"),
        "truth_beam_id": result.get("truth", {}).get("beam_id"),
        "truth_range_sample_float": result.get("truth", {}).get("range_sample_float"),
        "truth_range_sample_int": result.get("truth", {}).get("range_sample_int"),
        "truth_range_m": result.get("truth", {}).get("range_m"),
        "det_id": result.get("detection", {}).get("det_id"),
        "det_period_id": result.get("detection", {}).get("period_id"),
        "det_beam_id": result.get("detection", {}).get("beam_id"),
        "det_range_bin": result.get("detection", {}).get("range_bin"),
        "det_range_m": result.get("detection", {}).get("range_m"),
        "pc_crop_start": result.get("mapping", {}).get("pc_crop_start"),
        "expected_pc_bin_without_offset": result.get("mapping", {}).get("expected_pc_bin_without_offset"),
        "pc_bin_offset": result.get("mapping", {}).get("pc_bin_offset"),
        "bin_error": result.get("mapping", {}).get("bin_error"),
        "range_error_m": result.get("mapping", {}).get("range_error_m"),
        "direct_bin_compare_valid": result.get("mapping", {}).get("direct_bin_compare_valid"),
        "notes": result.get("notes", ""),
    }
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# 单点目标脉压峰标定报告\n\n")
        f.write(f"- 标定状态：`{result.get('calibration_status')}`\n")
        f.write(f"- case_id：`{result.get('case_id', '')}`\n")
        f.write(f"- truth_source：`{result.get('truth_source', '')}`\n")
        f.write(f"- truth range_sample_float：`{result.get('truth', {}).get('range_sample_float')}`\n")
        f.write(f"- detection range_bin：`{result.get('detection', {}).get('range_bin')}`\n")
        f.write(f"- direct_bin_compare_valid：`{result.get('mapping', {}).get('direct_bin_compare_valid')}`\n")
        f.write(f"- pc_crop_start：`{result.get('mapping', {}).get('pc_crop_start')}`\n")
        f.write(f"- expected_pc_bin_without_offset：`{result.get('mapping', {}).get('expected_pc_bin_without_offset')}`\n")
        f.write(f"- pc_bin_offset：`{result.get('mapping', {}).get('pc_bin_offset')}`\n")
        f.write(f"- bin_error：`{result.get('mapping', {}).get('bin_error')}`\n")
        f.write(f"- range_error_m：`{result.get('mapping', {}).get('range_error_m')}`\n")
        f.write(f"- 说明：{result.get('notes', '')}\n\n")
        if result.get("calibration_status") == "success":
            f.write("后续 P4 truth 匹配应优先使用本标定结果修正 expected bin；如果后续单点复核确认 offset 稳定，可将该 offset 作为同配置下的 range-bin 映射校正。\n")
        else:
            f.write("单点标定失败时，不建议继续做严肃检测率、定位误差或 SCNR 目标 ROI 统计。\n")

    return json_path, csv_path, report_path


def update_manifest(output_dir, paths):
    result_dir = os.path.dirname(os.path.abspath(output_dir))
    manifest_path = os.path.join(result_dir, "run_manifest.json")
    if not os.path.exists(manifest_path):
        return
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)
        generated = manifest.setdefault("generated_files", {})
        generated["single_point_bin_calibration_json"] = paths[0]
        generated["single_point_bin_calibration_csv"] = paths[1]
        generated["single_point_calibration_report"] = paths[2]
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, ensure_ascii=False, indent=2)
            f.write("\n")
    except Exception as exc:
        print(f"[calibrate_single_point_bin][WARN] update manifest failed: {exc}", file=sys.stderr)


def parse_bool_arg(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def first_value(row, names, default=""):
    for name in names:
        value = row.get(name, "")
        if value not in (None, ""):
            return value
    return default


def first_float(row, names, default=math.nan):
    return to_float(first_value(row, names, ""), default)


def first_int(row, names, default=None):
    return to_int(first_value(row, names, ""), default)


def moving_truth_injected(row):
    if "beam_injected" in row and row.get("beam_injected", "") != "":
        return truth_bool(row, "beam_injected")
    if "injection_enabled" in row and row.get("injection_enabled", "") != "":
        return truth_bool(row, "injection_enabled")
    if "injected_sample_count" in row and row.get("injected_sample_count", "") != "":
        return to_float(row.get("injected_sample_count"), 0.0) > 0.0
    return False


def moving_truth_visible(row):
    if "visible_by_beam" in row and row.get("visible_by_beam", "") != "":
        return truth_bool(row, "visible_by_beam")
    if "visible_pulse_count" in row and row.get("visible_pulse_count", "") != "":
        return to_float(row.get("visible_pulse_count"), 0.0) > 0.0
    return False


def moving_truth_valid(row):
    return moving_truth_injected(row) and moving_truth_visible(row)


def truth_visible(row):
    return moving_truth_valid(row)


def moving_truth_key(row):
    target_id = first_value(row, ["target_id", "target_name", "id"], "")
    period_id = first_value(row, ["period_id"], "0")
    beam_id = first_value(row, ["beam_id", "beam_id_1based"], "")
    if beam_id == "":
        b0 = to_int(row.get("beam_id_0based"), None)
        beam_id = "" if b0 is None else str(b0 + 1)
    return target_id, period_id, beam_id


def aggregate_pulse_truth(pulse_rows):
    groups = {}
    for row in pulse_rows:
        if not truth_visible(row):
            continue
        groups.setdefault(moving_truth_key(row), []).append(row)
    out = []
    numeric_fields = [
        "range_m", "range_sample_float", "range_sample_int", "e", "n",
        "target_e", "target_n", "lat", "lon", "target_lat", "target_lon",
        "theta_cmd_deg", "theta_true_deg", "target_azimuth_deg",
        "radial_velocity_mps", "vr_mps", "velocity_mps",
    ]
    for key, rows in groups.items():
        merged = dict(rows[len(rows) // 2])
        for field in numeric_fields:
            vals = [to_float(r.get(field)) for r in rows]
            vals = [v for v in vals if math.isfinite(v)]
            if vals:
                merged[field] = str(sum(vals) / len(vals))
        merged["target_id"], merged["period_id"], merged["beam_id"] = key
        out.append(merged)
    return out


def choose_moving_truths(pulse_rows, beam_rows, args):
    notes = []
    beam_candidates = [dict(r) for r in beam_rows if moving_truth_valid(r)]
    usable_beam = [
        r for r in beam_candidates
        if first_value(r, ["beam_id", "beam_id_1based", "beam_id_0based"], "") != ""
    ]
    if usable_beam:
        notes.append("truth_source=target_truth_beam_summary")
        truths = usable_beam
    else:
        notes.append("beam summary 缺少可用逐波位记录，已从 target_truth_pulse 聚合。")
        notes.append("truth_source=target_truth_pulse_aggregated")
        truths = aggregate_pulse_truth(pulse_rows)

    if args.target_id:
        truths = [r for r in truths if str(truth_target_id(r)) == str(args.target_id)]
        notes.append(f"filter_target_id={args.target_id}")
    if args.target_beam_id is not None:
        truths = [r for r in truths if truth_beam_id(r) == args.target_beam_id]
        notes.append(f"filter_target_beam_id={args.target_beam_id}")
    if args.single_truth_only and len(truths) > 1:
        truths = sorted(
            truths,
            key=lambda r: (
                -to_float(r.get("injected_sample_count"), 0.0),
                -to_float(r.get("visible_pulse_count"), 0.0),
                truth_beam_id(r) if truth_beam_id(r) is not None else 999999,
            ),
        )[:1]
        notes.append("single_truth_only=true，仅保留最强有效 truth。")
    return truths, notes


def truth_beam_id(row):
    b = first_int(row, ["beam_id", "beam_id_1based"], None)
    if b is None:
        b0 = first_int(row, ["beam_id_0based"], None)
        if b0 is not None:
            b = b0 + 1
    return b


def truth_period_id(row):
    return first_int(row, ["period_id"], 0)


def truth_target_id(row):
    return first_value(row, ["target_id", "target_name", "id"], "")


def det_beam_id_1based(det, truth, args):
    db = first_int(det, ["beam_id"], None)
    if db is None:
        return None
    if args.det_beam_id_base == "zero":
        return db + 1
    if args.det_beam_id_base == "one":
        return db
    tt = first_float(truth, ["theta_cmd_deg", "theta_true_deg", "target_azimuth_deg"])
    dt = first_float(det, ["theta_cmd_deg", "theta_true_deg"])
    tb = truth_beam_id(truth)
    if math.isfinite(tt) and math.isfinite(dt) and abs(dt - tt) <= 1.0:
        if tb is not None and db == tb:
            return db
        if tb is not None and db + 1 == tb:
            return db + 1
    return db


def same_beam_or_near(truth, det, args):
    tt = first_float(truth, ["theta_cmd_deg", "theta_true_deg", "target_azimuth_deg"])
    dt = first_float(det, ["theta_cmd_deg", "theta_true_deg"])
    if math.isfinite(tt) and math.isfinite(dt):
        angle_diff = abs(dt - tt)
        return angle_diff <= 1.0, angle_diff, "theta"
    tb = truth_beam_id(truth)
    db1 = det_beam_id_1based(det, truth, args)
    if tb is None or db1 is None:
        return False, None, "missing_beam_id"
    return abs(db1 - tb) <= args.max_beam_diff, abs(db1 - tb), "beam_id"


def expected_truth_bin(row, pc_crop_start, pc_bin_offset):
    sample = first_float(row, ["range_sample_float"], math.nan)
    if not math.isfinite(sample) or pc_crop_start is None:
        return None
    return int(round(sample)) - pc_crop_start + pc_bin_offset


def geo_error_m(lat1, lon1, lat2, lon2):
    if not all(math.isfinite(v) for v in [lat1, lon1, lat2, lon2]):
        return math.nan
    r = 6371000.0
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2.0) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    return 2.0 * r * math.atan2(math.sqrt(a), math.sqrt(max(0.0, 1.0 - a)))


def moving_candidate_row(truth, det, runtime_config, args):
    pc_crop_start = runtime_pc_crop_start(runtime_config)
    tid = truth_target_id(truth)
    tp = truth_period_id(truth)
    tb = truth_beam_id(truth)
    dp = first_int(det, ["period_id"], None)
    raw_db = first_int(det, ["beam_id"], None)
    db = det_beam_id_1based(det, truth, args)
    beam_ok, beam_metric, beam_match_mode = same_beam_or_near(truth, det, args)
    beam_diff = beam_metric if beam_match_mode == "beam_id" else (
        abs(db - tb) if db is not None and tb is not None else None
    )

    truth_range = first_float(truth, ["range_m"])
    det_range = first_float(det, ["range_m"])
    range_error = det_range - truth_range if math.isfinite(det_range) and math.isfinite(truth_range) else math.nan

    truth_bin = expected_truth_bin(truth, pc_crop_start, args.pc_bin_offset)
    det_bin = first_int(det, ["range_bin"], None)
    bin_error = det_bin - truth_bin if det_bin is not None and truth_bin is not None else None

    te = first_float(truth, ["target_e", "e", "x_m", "target_x"])
    tn = first_float(truth, ["target_n", "n", "y_m", "target_y"])
    de = first_float(det, ["e"])
    dn = first_float(det, ["n"])
    pos_error = math.hypot(de - te, dn - tn) if all(math.isfinite(v) for v in [te, tn, de, dn]) else math.nan

    tlat = first_float(truth, ["lat", "target_lat"])
    tlon = first_float(truth, ["lon", "target_lon"])
    dlat = first_float(det, ["lat"])
    dlon = first_float(det, ["lon"])
    geo_err = geo_error_m(tlat, tlon, dlat, dlon)

    tt = first_float(truth, ["theta_true_deg", "theta_cmd_deg", "target_azimuth_deg"])
    dt = first_float(det, ["theta_true_deg", "theta_cmd_deg"])
    angle_error = dt - tt if math.isfinite(dt) and math.isfinite(tt) else math.nan

    tv = first_float(truth, ["radial_velocity_mps", "vr_mps", "velocity_mps"])
    dv = first_float(det, ["radial_velocity_mps"])
    vel_error = dv - tv if math.isfinite(dv) and math.isfinite(tv) else math.nan

    reasons = []
    soft_reasons = []
    if tp is not None and dp is not None and tp > 0 and dp != tp:
        reasons.append("period_mismatch")
    if not beam_ok:
        reasons.append("beam_diff_exceeded")
    if bin_error is not None and abs(bin_error) > args.max_range_bin_error:
        reasons.append("range_bin_gate_failed")
    if args.match_by in {"range", "position"} and math.isfinite(range_error) and abs(range_error) > args.max_range_error_m:
        reasons.append("range_gate_failed")
    if math.isfinite(pos_error) and pos_error > args.max_position_error_m:
        if args.ignore_position_gate:
            soft_reasons.append("position_error_large")
        else:
            reasons.append("position_gate_failed")
    if math.isfinite(vel_error) and abs(vel_error) > args.max_velocity_error_mps:
        reasons.append("velocity_gate_failed")

    if args.match_by == "range_bin" and bin_error is not None:
        total_cost = abs(bin_error) * 100.0 + (abs(range_error) if math.isfinite(range_error) else 0.0)
    elif args.match_by == "range" and math.isfinite(range_error):
        total_cost = abs(range_error) + float(beam_diff or 0) * 100.0
    elif args.match_by == "position" and math.isfinite(pos_error):
        total_cost = pos_error
    elif math.isfinite(range_error):
        total_cost = abs(range_error) + float(beam_diff or 0) * 100.0
    elif bin_error is not None:
        total_cost = abs(bin_error) * 100.0
    else:
        total_cost = 1.0e18
        reasons.append("missing_position_range_and_bin")

    if not beam_ok and beam_diff is None:
        reasons.append("missing_beam_id")
    if not math.isfinite(truth_range):
        reasons.append("missing_truth_range_m")
    if not math.isfinite(det_range):
        reasons.append("missing_det_range_m")

    gate_passed = not reasons
    return {
        "target_id": tid,
        "period_id": tp,
        "beam_id": tb,
        "det_id": first_value(det, ["det_id"], ""),
        "det_period_id": dp,
        "det_beam_id": raw_db,
        "det_beam_id_1based": db,
        "beam_diff": beam_diff,
        "truth_range_m": truth_range,
        "det_range_m": det_range,
        "range_error_m": range_error,
        "truth_range_bin_expected": truth_bin,
        "det_range_bin": det_bin,
        "range_bin_error": bin_error,
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
        "total_cost": total_cost,
        "gate_passed": gate_passed,
        "gate_fail_reason": ";".join(reasons),
        "reason": ";".join(soft_reasons),
    }


def run_moving_target_localization(args):
    runtime_config = read_runtime_config(args.runtime_config)
    detections = read_csv(args.detections)
    pulse_rows = read_csv(args.truth_pulse)
    beam_rows = read_csv(args.truth_beam)
    truths, notes = choose_moving_truths(pulse_rows, beam_rows, args)

    candidate_rows = []
    best_by_truth = {}
    for ti, truth in enumerate(truths):
        best = None
        for det in detections:
            cand = moving_candidate_row(truth, det, runtime_config, args)
            cand["_truth_index"] = ti
            candidate_rows.append(cand)
            if cand["gate_passed"] and (best is None or cand["total_cost"] < best["total_cost"]):
                best = cand
        best_by_truth[ti] = best

    assigned = {}
    for ti, cand in best_by_truth.items():
        if cand is None:
            continue
        det_key = (cand.get("det_id"), cand.get("det_period_id"), cand.get("det_beam_id"), cand.get("det_range_bin"))
        prev = assigned.get(det_key)
        if prev is None or cand["total_cost"] < prev["total_cost"]:
            assigned[det_key] = cand

    result_rows = []
    used_truth = set()
    used_dets = set()
    for ti, truth in enumerate(truths):
        cand = best_by_truth.get(ti)
        if cand is None:
            base = moving_candidate_row(truth, {}, runtime_config, args)
            base["match_status"] = "unmatched_truth"
            base["reason"] = "no_gate_passed_candidate"
            result_rows.append(base)
            continue
        det_key = (cand.get("det_id"), cand.get("det_period_id"), cand.get("det_beam_id"), cand.get("det_range_bin"))
        if assigned.get(det_key) is not cand:
            cand = dict(cand)
            cand["match_status"] = "duplicate_det"
            cand["reason"] = "same_detection_selected_by_lower_cost_truth"
        else:
            cand = dict(cand)
            cand["match_status"] = "matched"
            cand["reason"] = cand.get("reason", "")
            used_truth.add(ti)
            used_dets.add(det_key)
        result_rows.append(cand)

    det_keys_all = set()
    for d in detections:
        det_keys_all.add((first_value(d, ["det_id"], ""), first_int(d, ["period_id"], None), first_int(d, ["beam_id"], None), first_int(d, ["range_bin"], None)))

    matched = [r for r in result_rows if r.get("match_status") == "matched"]
    def mean_field(rows, field):
        vals = [r.get(field) for r in rows]
        vals = [v for v in vals if isinstance(v, (int, float)) and math.isfinite(v)]
        return sum(vals) / len(vals) if vals else None

    status = "passed" if matched and len(matched) == len(truths) else ("partial" if matched else "failed")
    summary = {
        "status": status,
        "truth_count": len(truths),
        "detection_count": len(detections),
        "candidate_count": len(candidate_rows),
        "matched_count": len(matched),
        "unmatched_truth_count": len([r for r in result_rows if r.get("match_status") == "unmatched_truth"]),
        "false_detection_count": max(0, len(det_keys_all - used_dets)),
        "mean_range_error_m": mean_field(matched, "range_error_m"),
        "mean_range_bin_error": mean_field(matched, "range_bin_error"),
        "mean_pos_error_m": mean_field(matched, "pos_error_m"),
        "mean_velocity_error_mps": mean_field(matched, "velocity_error_mps"),
        "pc_bin_offset_used": args.pc_bin_offset,
        "notes": notes,
    }

    write_moving_outputs(args.output_dir, result_rows, candidate_rows, summary, args)
    return summary, result_rows


def csv_value(value):
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value


def write_dict_csv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: csv_value(row.get(k, "")) for k in fieldnames})


def write_moving_outputs(output_dir, result_rows, candidate_rows, summary, args):
    os.makedirs(output_dir, exist_ok=True)
    result_path = os.path.join(output_dir, "合作动目标定位标定结果.csv")
    detail_path = os.path.join(output_dir, "合作动目标候选匹配明细.csv")
    json_path = os.path.join(output_dir, "合作动目标定位标定汇总.json")
    report_path = os.path.join(output_dir, "合作动目标定位标定报告.md")

    result_fields = [
        "case_id", "run_id", "result_id",
        "target_id", "period_id", "beam_id",
        "det_id", "det_period_id", "det_beam_id",
        "truth_range_m", "det_range_m", "range_error_m",
        "truth_range_bin_expected", "det_range_bin", "range_bin_error",
        "truth_e", "truth_n", "det_e", "det_n", "pos_error_m",
        "truth_lat", "truth_lon", "det_lat", "det_lon", "geo_error_m",
        "truth_theta_deg", "det_theta_deg", "angle_error_deg",
        "truth_radial_velocity_mps", "det_radial_velocity_mps", "velocity_error_mps",
        "match_status", "reason",
    ]
    case_id = nested_get(read_runtime_config(args.runtime_config), ["run_info", "case_id"], "")
    run_id = nested_get(read_runtime_config(args.runtime_config), ["run_info", "run_id"], "")
    result_id = nested_get(read_runtime_config(args.runtime_config), ["run_info", "result_id"], "")
    for row in result_rows:
        row.setdefault("case_id", case_id)
        row.setdefault("run_id", run_id)
        row.setdefault("result_id", result_id)
        row.setdefault("reason", row.get("gate_fail_reason", ""))
    write_dict_csv(result_path, result_rows, result_fields)

    detail_fields = [
        "target_id", "period_id", "beam_id",
        "det_id", "det_period_id", "det_beam_id",
        "beam_diff",
        "truth_range_m", "det_range_m", "range_error_m",
        "truth_range_bin_expected", "det_range_bin", "range_bin_error",
        "pos_error_m",
        "velocity_error_mps",
        "total_cost",
        "gate_passed",
        "gate_fail_reason",
    ]
    write_dict_csv(detail_path, candidate_rows, detail_fields)

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
        f.write("\n")

    if args.write_docs:
        with open(report_path, "w", encoding="utf-8") as f:
            f.write("# 合作动目标定位标定报告\n\n")
            f.write("## 目的\n\n检查合作动目标 truth 与 `detection_results.csv` 的最终定位输出是否一致。\n\n")
            f.write("## 输入\n\n")
            f.write(f"- truth_pulse: `{args.truth_pulse}`\n")
            f.write(f"- truth_beam: `{args.truth_beam}`\n")
            f.write(f"- detections: `{args.detections}`\n")
            f.write(f"- runtime_config: `{args.runtime_config}`\n\n")
            f.write("## 门限\n\n")
            f.write(f"- pc_bin_offset: `{args.pc_bin_offset}`\n")
            f.write(f"- max_beam_diff: `{args.max_beam_diff}`\n")
            f.write(f"- match_by: `{args.match_by}`\n")
            f.write(f"- max_range_bin_error: `{args.max_range_bin_error}`\n")
            f.write(f"- max_range_error_m: `{args.max_range_error_m}`\n")
            f.write(f"- max_position_error_m: `{args.max_position_error_m}`\n")
            f.write(f"- ignore_position_gate: `{args.ignore_position_gate}`\n")
            f.write(f"- det_beam_id_base: `{args.det_beam_id_base}`\n")
            f.write(f"- max_velocity_error_mps: `{args.max_velocity_error_mps}`\n\n")
            f.write("## 汇总\n\n")
            for key in ["status", "truth_count", "detection_count", "candidate_count", "matched_count", "unmatched_truth_count", "false_detection_count"]:
                f.write(f"- {key}: `{summary.get(key)}`\n")
            f.write(f"- mean_range_error_m: `{summary.get('mean_range_error_m')}`\n")
            f.write(f"- mean_range_bin_error: `{summary.get('mean_range_bin_error')}`\n")
            f.write(f"- mean_pos_error_m: `{summary.get('mean_pos_error_m')}`\n")
            f.write(f"- mean_velocity_error_mps: `{summary.get('mean_velocity_error_mps')}`\n\n")
            if summary.get("matched_count", 0) > 0:
                f.write("结论：存在合作动目标匹配，可用本结果检查 range/bin/position/velocity 输出误差。是否进入 P3/P4 取决于误差是否满足当前实验门限。\n")
            else:
                f.write("结论：未形成有效匹配。请优先查看候选匹配明细中的 gate_fail_reason，再决定是否调整数据生成或检查检测输出。\n")


def main():
    parser = argparse.ArgumentParser(description="Calibrate single point truth bin to detection bin.")
    parser.add_argument("--mode", choices=["single_point_bin", "moving_target_localization"], default="single_point_bin")
    parser.add_argument("--truth-pulse")
    parser.add_argument("--truth-beam")
    parser.add_argument("--scene-truth")
    parser.add_argument("--detections", required=True)
    parser.add_argument("--runtime-config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--pc-bin-offset", type=int, default=0)
    parser.add_argument("--max-beam-diff", type=int, default=1)
    parser.add_argument("--target-id")
    parser.add_argument("--target-beam-id", type=int, default=None)
    parser.add_argument("--single-truth-only", default="true")
    parser.add_argument("--match-by", choices=["range_bin", "range", "position"], default="range_bin")
    parser.add_argument("--max-range-bin-error", type=float, default=20.0)
    parser.add_argument("--ignore-position-gate", default="true")
    parser.add_argument("--det-beam-id-base", choices=["auto", "zero", "one"], default="auto")
    parser.add_argument(
        "--max-range-error-m",
        type=float,
        default=None,
        help="Maximum allowed absolute range_m error. Defaults: 1000 for single_point_bin, 300 for moving_target_localization.",
    )
    parser.add_argument("--max-position-error-m", type=float, default=500.0)
    parser.add_argument("--max-velocity-error-mps", type=float, default=20.0)
    parser.add_argument("--write-docs", default="true")
    args = parser.parse_args()
    args.write_docs = parse_bool_arg(args.write_docs)
    args.single_truth_only = parse_bool_arg(args.single_truth_only)
    args.ignore_position_gate = parse_bool_arg(args.ignore_position_gate)
    if args.max_range_error_m is None:
        args.max_range_error_m = 300.0 if args.mode == "moving_target_localization" else 1000.0

    try:
        if args.mode == "moving_target_localization":
            if not args.truth_pulse or not args.truth_beam:
                raise ValueError("moving_target_localization 模式需要同时传 --truth-pulse 和 --truth-beam。")
            summary, _ = run_moving_target_localization(args)
            print(
                "[P2.5] status={status} truth={truth_count} det={detection_count} matched={matched_count}".format(**summary)
            )
            return 0

        runtime_config = read_runtime_config(args.runtime_config)
        det_rows = read_csv(args.detections)

        if args.scene_truth:
            if args.truth_pulse or args.truth_beam:
                raise ValueError("--scene-truth 模式不能同时传 --truth-pulse/--truth-beam。")
            truth, error, notes = choose_scene_truth(read_csv(args.scene_truth), runtime_config)
        else:
            if not args.truth_pulse or not args.truth_beam:
                raise ValueError("合作目标模式需要同时传 --truth-pulse 和 --truth-beam；静止散射点模式请传 --scene-truth。")
            truth_rows = read_csv(args.truth_pulse)
            beam_rows = read_csv(args.truth_beam)
            truth, error, notes = choose_truth(truth_rows, beam_rows)
        if truth is None:
            result = failure_result(error, runtime_config, notes)
        else:
            detection, detection_error = choose_detection(det_rows, truth, args.max_range_error_m)
            if detection is None:
                result = failure_result(detection_error, runtime_config, notes)
                result["truth"] = {
                    "target_id": to_int(truth.get("target_id"), None),
                    "period_id": to_int(truth.get("period_id"), None),
                    "beam_id": truth_beam_id(truth),
                    "range_sample_float": to_float(truth.get("range_sample_float")),
                    "range_sample_int": to_int(truth.get("range_sample_int"), None),
                    "range_m": to_float(truth.get("range_m")),
                }
                attach_failure_truth_mapping(result)
            else:
                result = make_result(truth, detection, runtime_config, notes)

        paths = write_outputs(result, args.output_dir)
        update_manifest(args.output_dir, paths)
    except Exception as exc:
        print(f"[calibrate_single_point_bin][ERR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
