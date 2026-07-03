#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path


def to_float(v, default=math.nan):
    try:
        if v is None or v == "":
            return default
        return float(v)
    except Exception:
        return default


def to_int(v, default=None):
    try:
        if v is None or v == "":
            return default
        return int(float(v))
    except Exception:
        return default


def read_csv(path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def choose_truth(rows, target_id=None, beam_id=None):
    candidates = [
        r for r in rows
        if to_float(r.get("injected_sample_count"), 0.0) > 0.0
    ]
    if target_id is not None:
        candidates = [r for r in candidates if str(r.get("target_id", "")) == str(target_id)]
    if beam_id is not None:
        candidates = [r for r in candidates if to_int(r.get("beam_id")) == beam_id]
    if not candidates:
        raise SystemExit("no injected truth row matched")
    candidates.sort(key=lambda r: -to_float(r.get("injected_sample_count"), 0.0))
    return candidates[0]


def choose_detection(rows, truth, beam_id=None, det_id=None):
    expected_bin = to_int(truth.get("expected_range_bin"), None)
    truth_beam = to_int(truth.get("beam_id"), None)
    if beam_id is not None:
        truth_beam = beam_id
    candidates = list(rows)
    if det_id is not None:
        candidates = [r for r in candidates if to_int(r.get("det_id")) == det_id]
    if truth_beam is not None:
        candidates = [r for r in candidates if to_int(r.get("beam_id")) == truth_beam]
    if expected_bin is not None:
        candidates.sort(key=lambda r: abs((to_int(r.get("range_bin"), 10**9) or 10**9) - expected_bin))
    if not candidates:
        raise SystemExit("no detection row matched")
    return candidates[0]


def choose_trace(rows, det):
    beam = to_int(det.get("beam_id"))
    rng = to_int(det.get("range_bin"))
    row = to_int(det.get("row"))
    det_id = to_int(det.get("det_id"))
    candidates = [
        r for r in rows
        if to_int(r.get("beam_id")) == beam
        and to_int(r.get("range_bin")) == rng
        and to_int(r.get("row")) == row
    ]
    if det_id is not None:
        exact = [r for r in candidates if to_int(r.get("det_id")) == det_id]
        if exact:
            candidates = exact
    if not candidates:
        return None
    return candidates[0]


def fmt(v):
    if isinstance(v, float):
        if math.isnan(v):
            return "nan"
        return f"{v:.12g}"
    return str(v)


def main():
    ap = argparse.ArgumentParser(
        description="Compare stage2 physical truth with GMTI fusion localization trace."
    )
    ap.add_argument("--case", required=True, help="stage2 output case directory")
    ap.add_argument("--target-id", default=None)
    ap.add_argument("--beam-id", type=int, default=None, help="1-based beam id")
    ap.add_argument("--det-id", type=int, default=None)
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    case = Path(args.case)
    truth_rows = read_csv(case / "truth" / "target_truth_beam_summary.csv")
    det_rows = read_csv(case / "algorithm_result" / "detection_results.csv")
    trace_rows = read_csv(case / "algorithm_result" / "fusion_localization_trace.csv")
    runtime_path = case / "algorithm_result" / "runtime_config_dump.json"
    runtime = {}
    if runtime_path.exists():
        runtime = json.loads(runtime_path.read_text(encoding="utf-8"))

    truth = choose_truth(truth_rows, args.target_id, args.beam_id)
    det = choose_detection(det_rows, truth, args.beam_id, args.det_id)
    trace = choose_trace(trace_rows, det) if trace_rows else None

    truth_e = to_float(truth.get("target_e"))
    truth_n = to_float(truth.get("target_n"))
    det_e = to_float(det.get("e"))
    det_n = to_float(det.get("n"))
    de = det_e - truth_e
    dn = det_n - truth_n
    pos = math.hypot(de, dn)

    lines = []
    lines.append("# Stage2 Fusion Localization Compare")
    lines.append("")
    lines.append("## Matched Rows")
    lines.append(f"- case: `{case}`")
    lines.append(f"- truth target_id: `{truth.get('target_id', '')}`")
    lines.append(f"- truth beam_id: `{truth.get('beam_id', '')}`")
    lines.append(f"- detection det_id: `{det.get('det_id', '')}`")
    lines.append(f"- detection beam_id: `{det.get('beam_id', '')}`")
    lines.append("")
    lines.append("## Physical Truth vs Final Detection")
    rows = [
        ("expected_bin", truth.get("expected_range_bin", "")),
        ("detection_range_bin", det.get("range_bin", "")),
        ("range_m", det.get("range_m", "")),
        ("theta_cmd_deg", det.get("theta_cmd_deg", "")),
        ("theta_true_deg", det.get("theta_true_deg", "")),
        ("physical_truth_e", truth_e),
        ("physical_truth_n", truth_n),
        ("final_output_e", det_e),
        ("final_output_n", det_n),
        ("de_m", de),
        ("dn_m", dn),
        ("pos_error_m", pos),
    ]
    for k, v in rows:
        lines.append(f"- {k}: `{fmt(v)}`")

    lines.append("")
    lines.append("## Fusion Localization Trace")
    if trace is None:
        lines.append("- trace: `missing`")
        lines.append("- action: rerun GMTI_core with runtime diagnostics enabled to create `algorithm_result/fusion_localization_trace.csv`.")
    else:
        trace_fields = [
            "period_id", "beam_id", "det_id", "range_bin", "range_m", "row", "col",
            "theta_cmd_deg", "theta_true_deg", "theta_used_deg",
            "ref_platform_e", "ref_platform_n", "ref_platform_lat", "ref_platform_lon",
            "ref_platform_v", "ref_platform_v_angle_deg", "ref_platform_ve", "ref_platform_vn",
            "af_wrapped", "af_ransac", "fa_shift", "af_used", "sinA", "px", "py",
            "ground_range_m", "look_e_used", "look_n_used",
            "direct_localized_e", "direct_localized_n",
            "fused_localized_e", "fused_localized_n",
            "final_output_e", "final_output_n",
            "beam_center_dir_deg", "target_direction_deg", "beam_dir_err_deg",
        ]
        for k in trace_fields:
            lines.append(f"- {k}: `{trace.get(k, '')}`")
        alg_e = to_float(trace.get("final_output_e"))
        alg_n = to_float(trace.get("final_output_n"))
        lines.append("")
        lines.append("## Algorithm Expected vs Physical Truth")
        lines.append(f"- algorithm_expected_e: `{fmt(alg_e)}`")
        lines.append(f"- algorithm_expected_n: `{fmt(alg_n)}`")
        lines.append(f"- algorithm_minus_physical_de_m: `{fmt(alg_e - truth_e)}`")
        lines.append(f"- algorithm_minus_physical_dn_m: `{fmt(alg_n - truth_n)}`")
        lines.append(f"- algorithm_minus_physical_pos_m: `{fmt(math.hypot(alg_e - truth_e, alg_n - truth_n))}`")

    if runtime:
        lines.append("")
        lines.append("## Runtime")
        for k in ["enable_dbs_fusion", "estimate_error_angle", "squint_angle", "loc_beam_gate_deg", "squint_side"]:
            if k in runtime:
                lines.append(f"- {k}: `{runtime[k]}`")

    text = "\n".join(lines) + "\n"
    if args.output:
        out = Path(args.output)
    else:
        out = case / "reports" / "fusion_localization_compare.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(text)
    print(f"[compare] wrote {out}")


if __name__ == "__main__":
    main()
