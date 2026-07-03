#!/usr/bin/env python3
"""Check stage2 simulation geometry consistency for a generated case."""

import argparse
import csv
import json
import math
import struct
import xml.etree.ElementTree as ET
from pathlib import Path


def to_float(v, default=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def to_int(v, default=-1):
    try:
        return int(round(float(v)))
    except (TypeError, ValueError):
        return default


def read_rows(path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def xml_text(root, name, default=""):
    elem = root.find(".//" + name)
    return elem.text.strip() if elem is not None and elem.text else default


def read_header(path, packet_bytes, prt_index):
    with path.open("rb") as f:
        f.seek(packet_bytes * prt_index)
        h = f.read(256)
    if len(h) != 256:
        raise RuntimeError(f"failed to read header prt={prt_index}")
    return {
        "utc": struct.unpack_from("<f", h, 16)[0],
        "lat": struct.unpack_from("<d", h, 104)[0],
        "lon": struct.unpack_from("<d", h, 112)[0],
        "height": struct.unpack_from("<d", h, 120)[0],
        "vn": struct.unpack_from("<f", h, 128)[0],
        "ve": struct.unpack_from("<f", h, 132)[0],
        "vd": struct.unpack_from("<f", h, 136)[0],
        "theta": struct.unpack_from("<h", h, 218)[0] / 100.0,
    }


def add_check(checks, name, passed, value=None, tolerance=None):
    checks.append({
        "name": name,
        "passed": bool(passed),
        "value": value,
        "tolerance": tolerance,
    })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, type=Path)
    ap.add_argument("--expected-bin", type=int, default=None)
    args = ap.parse_args()

    case = args.case
    reports = case / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    xml_path = case / "config" / "temp_config_stage2_newsystem.xml"
    truth_path = case / "truth" / "target_truth_beam_summary.csv"
    pulse_path = case / "truth" / "target_truth_pulse.csv"
    if not xml_path.exists():
        raise SystemExit(f"missing XML: {xml_path}")
    if not truth_path.exists():
        raise SystemExit(f"missing truth: {truth_path}")

    root = ET.parse(xml_path).getroot()
    pulse_len = to_int(xml_text(root, "pulse_len"))
    pulse_num = to_int(xml_text(root, "pulse_num"))
    beam_count = to_int(xml_text(root, "az_count"))
    range_crop_start = to_int(xml_text(root, "range_crop_start"))
    data_file = xml_text(root, "GMTI_data_new")
    data_path = Path(data_file)
    if not data_path.is_absolute():
        data_path = Path.cwd() / data_path

    rows = read_rows(truth_path)
    injected = [r for r in rows if to_float(r.get("injected_sample_count"), 0.0) > 0.0]
    checks = []
    add_check(checks, "one injected beam-summary row", len(injected) == 1, len(injected))
    if not injected:
        status = "FAIL"
        result = {"status": status, "checks": checks}
    else:
        t = injected[0]
        expected_range_bin = to_int(t.get("expected_range_bin"), None)
        if expected_range_bin is None or expected_range_bin < 0:
            expected_range_bin = int(round(to_float(t.get("range_sample_float"))) - range_crop_start)
        expected_bin = args.expected_bin if args.expected_bin is not None else expected_range_bin
        add_check(checks, "expected_range_bin", expected_range_bin == expected_bin,
                  {"actual": expected_range_bin, "expected": expected_bin})

        ref_e = to_float(t.get("ref_platform_e"))
        ref_n = to_float(t.get("ref_platform_n"))
        target_e = to_float(t.get("target_e"))
        target_n = to_float(t.get("target_n"))
        look_e = to_float(t.get("look_e"))
        look_n = to_float(t.get("look_n"))
        ground_range = to_float(t.get("ground_range_m"))
        pred_e = ref_e + ground_range * look_e
        pred_n = ref_n + ground_range * look_n
        en_err = math.hypot(target_e - pred_e, target_n - pred_n)
        add_check(checks, "target_e/n follows ref + ground_range * look", en_err < 1.0e-3,
                  en_err, "1e-3 m")
        look_norm_err = abs(math.hypot(look_e, look_n) - 1.0)
        add_check(checks, "look vector unit length", look_norm_err < 1.0e-9,
                  look_norm_err, "1e-9")

        inv_ground = math.hypot(target_e - ref_e, target_n - ref_n)
        add_check(checks, "target-ref distance equals ground_range",
                  abs(inv_ground - ground_range) < 1.0e-3,
                  abs(inv_ground - ground_range), "1e-3 m")

        ref_pulse_idx = to_int(t.get("ref_pulse_idx"))
        beam_id = to_int(t.get("beam_id_0based", t.get("beam_id")))
        prt_index = beam_id * pulse_num + ref_pulse_idx
        packet_bytes = 256 + pulse_len * 16
        header = read_header(data_path, packet_bytes, prt_index)
        lat_err = abs(header["lat"] - to_float(t.get("ref_platform_lat")))
        lon_err = abs(header["lon"] - to_float(t.get("ref_platform_lon")))
        add_check(checks, "POS header lat matches ref_platform_lat", lat_err < 1.0e-7,
                  lat_err, "1e-7 deg")
        add_check(checks, "POS header lon matches ref_platform_lon", lon_err < 1.0e-7,
                  lon_err, "1e-7 deg")
        add_check(checks, "POS header lat/lon valid",
                  abs(header["lat"]) > 1.0 and abs(header["lon"]) > 1.0,
                  {"lat": header["lat"], "lon": header["lon"]})
        add_check(checks, "POS header height near platform z",
                  abs(header["height"] - to_float(t.get("ref_platform_z"))) < 1.0e-6,
                  abs(header["height"] - to_float(t.get("ref_platform_z"))), "1e-6 m")
        if "ref_platform_ve" in t and "ref_platform_vn" in t:
            add_check(checks, "POS header ve matches truth",
                      abs(header["ve"] - to_float(t.get("ref_platform_ve"))) < 1.0e-5,
                      abs(header["ve"] - to_float(t.get("ref_platform_ve"))), "1e-5 m/s")
            add_check(checks, "POS header vn matches truth",
                      abs(header["vn"] - to_float(t.get("ref_platform_vn"))) < 1.0e-5,
                      abs(header["vn"] - to_float(t.get("ref_platform_vn"))), "1e-5 m/s")

        if pulse_path.exists():
            pulse_rows = read_rows(pulse_path)
            inj_pulses = [r for r in pulse_rows if to_int(r.get("injection_enabled"), 0) == 1]
            add_check(checks, "injected pulse rows present", len(inj_pulses) > 0, len(inj_pulses))

        status = "PASS" if all(c["passed"] for c in checks) else "FAIL"
        result = {
            "status": status,
            "case": str(case),
            "geometry_config_name": t.get("geometry_config_name", ""),
            "expected_range_bin": expected_range_bin,
            "truth": {
                "ref_platform_e": ref_e,
                "ref_platform_n": ref_n,
                "target_e": target_e,
                "target_n": target_n,
                "look_e": look_e,
                "look_n": look_n,
                "ground_range_m": ground_range,
            },
            "header": header,
            "checks": checks,
        }

    json_path = reports / "geometry_consistency_check.json"
    md_path = reports / "geometry_consistency_report.md"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = [
        "# Stage2 Geometry Consistency Report",
        "",
        f"- case: `{case}`",
        f"- status: `{result['status']}`",
        f"- json: `{json_path}`",
        "",
        "## Checks",
        "",
    ]
    for c in checks:
        lines.append(f"- [{'PASS' if c['passed'] else 'FAIL'}] {c['name']}: `{c.get('value')}`")
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(result["status"])
    print(md_path)
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
