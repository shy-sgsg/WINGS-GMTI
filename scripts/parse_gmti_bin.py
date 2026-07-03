#!/usr/bin/env python3
"""Parse GMTIxx.bin detection files.

Structure source: src/writeresults.cpp::GMTIProcessor::writeResult.
Current layout is little-endian:
  uint16 target_count
  repeated target_count times:
    uint16 id
    int32 lon_quantized, LSB = 8.38191e-8 deg
    int32 lat_quantized, LSB = 8.38191e-8 deg
    uint16 speed, scale = 0.01 m/s
    float64 direction_deg
    float64 range_m
    float64 utc
"""

import argparse
import csv
import json
import os
import struct
import sys

LSB_DEG = 8.38191e-8
HEADER_BYTES = 2
REC_BYTES = 36

FIELDS = [
    "case_id", "run_id", "result_id", "period_id", "beam_id", "det_id",
    "range_bin", "range_m", "theta_cmd_deg", "theta_true_deg", "row", "col",
    "e", "n", "lat", "lon", "utc", "amplitude", "power",
    "radial_velocity_mps", "source_file",
]


def read_runtime_context(path):
    context = {"case_id": "", "run_id": "", "result_id": ""}
    if not path:
        return context
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    run_info = data.get("run_info", {})
    context["case_id"] = str(run_info.get("case_id", ""))
    context["run_id"] = str(run_info.get("run_id", ""))
    context["result_id"] = str(run_info.get("result_id", ""))
    return context


def infer_result_id(path):
    stem = os.path.splitext(os.path.basename(path))[0]
    return stem if stem.startswith("GMTI") else ""


def parse_file(path):
    with open(path, "rb") as f:
        data = f.read()
    if len(data) < HEADER_BYTES:
        if len(data) == 0:
            return []
        raise ValueError(f"file too small for GMTI header: {len(data)} bytes")

    (count,) = struct.unpack_from("<H", data, 0)
    expected = HEADER_BYTES + count * REC_BYTES
    if len(data) != expected:
        raise ValueError(
            f"unsupported GMTI binary size: got {len(data)} bytes, "
            f"expected {expected} for count={count} and record={REC_BYTES}"
        )

    rows = []
    off = HEADER_BYTES
    for i in range(count):
        det_id, lon_q, lat_q, speed_q, direction, range_m, utc = struct.unpack_from(
            "<HiiHddd", data, off
        )
        off += REC_BYTES
        rows.append({
            "det_id": i,
            "range_m": range_m,
            "lat": lat_q * LSB_DEG,
            "lon": lon_q * LSB_DEG,
            "utc": utc,
            "radial_velocity_mps": speed_q * 0.01,
            "source_file": path,
        })
    return rows


def write_csv(path, rows, context, source_file):
    result_id = context.get("result_id") or infer_result_id(source_file)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            out = {k: "" for k in FIELDS}
            out.update({
                "case_id": context.get("case_id", ""),
                "run_id": context.get("run_id", ""),
                "result_id": result_id,
                "period_id": "",
                "beam_id": "",
                "range_bin": "",
                "theta_cmd_deg": "",
                "theta_true_deg": "",
                "row": "",
                "col": "",
                "e": "",
                "n": "",
                "amplitude": "",
                "power": "",
            })
            out.update(row)
            writer.writerow(out)


def main():
    parser = argparse.ArgumentParser(description="Parse GMTIxx.bin into CSV.")
    parser.add_argument("--input", required=True, help="Input GMTIxx.bin")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--config", default="", help="Optional runtime_config_dump.json")
    args = parser.parse_args()

    try:
        rows = parse_file(args.input)
        context = read_runtime_context(args.config)
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        write_csv(args.output, rows, context, args.input)
    except Exception as exc:
        print(f"[parse_gmti_bin][ERR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
