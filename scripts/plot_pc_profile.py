#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def find_profile(debug_dir: Path, beam_id: int | None):
    files = sorted(debug_dir.glob("pc_range_profile_beam*.csv"))
    if not files:
        raise FileNotFoundError(f"No pc_range_profile_beam*.csv found in {debug_dir}")

    if beam_id is None:
        return files[0]

    key1 = f"beam{beam_id}.csv"
    key2 = f"beam{beam_id:02d}.csv"

    for f in files:
        if key1 in f.name or key2 in f.name:
            return f

    raise FileNotFoundError(f"No profile found for beam {beam_id} in {debug_dir}")


def find_amp2d(debug_dir: Path, beam_id: int | None):
    files = sorted(debug_dir.glob("pc_amplitude_2d_beam*.csv"))
    if not files:
        return None

    if beam_id is None:
        return files[0]

    key1 = f"beam{beam_id}.csv"
    key2 = f"beam{beam_id:02d}.csv"

    for f in files:
        if key1 in f.name or key2 in f.name:
            return f

    return None


def plot_1d(profile_path: Path, expected_bin: int, window: int, out_dir: Path):
    df = pd.read_csv(profile_path)

    required = {"bin", "range_m", "power", "power_db"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{profile_path} missing columns: {missing}")

    df["bin"] = pd.to_numeric(df["bin"], errors="coerce")
    df["range_m"] = pd.to_numeric(df["range_m"], errors="coerce")
    df["power"] = pd.to_numeric(df["power"], errors="coerce")
    df["power_db"] = pd.to_numeric(df["power_db"], errors="coerce")
    df = df.dropna(subset=["bin", "power", "power_db"]).copy()

    global_peak = df.loc[df["power"].idxmax()]

    left = expected_bin - window
    right = expected_bin + window
    win = df[(df["bin"] >= left) & (df["bin"] <= right)].copy()

    if len(win) == 0:
        raise ValueError(f"No data in expected window [{left}, {right}]")

    local_peak = win.loc[win["power"].idxmax()]

    print("==== Pulse compression profile ====")
    print(f"profile: {profile_path}")
    print(f"expected_bin: {expected_bin}")
    print(
        "local_peak_bin:",
        int(local_peak["bin"]),
        "offset:",
        int(local_peak["bin"]) - expected_bin,
        "range_m:",
        float(local_peak["range_m"]),
        "power_db:",
        float(local_peak["power_db"]),
    )
    print(
        "global_peak_bin:",
        int(global_peak["bin"]),
        "range_m:",
        float(global_peak["range_m"]),
        "power_db:",
        float(global_peak["power_db"]),
    )

    print("\n==== Top peaks in expected window ====")
    print(
        win.sort_values("power", ascending=False)
        .head(20)[["bin", "range_m", "power_db"]]
        .to_string(index=False)
    )

    beam_name = profile_path.stem.replace("pc_range_profile_", "")

    # Full profile
    plt.figure(figsize=(12, 5))
    plt.plot(df["bin"], df["power_db"], linewidth=1.0)
    plt.axvline(expected_bin, linestyle="--", linewidth=1.0, label=f"expected={expected_bin}")
    plt.axvline(local_peak["bin"], linestyle=":", linewidth=1.2, label=f"local peak={int(local_peak['bin'])}")
    plt.axvline(global_peak["bin"], linestyle="-.", linewidth=1.2, label=f"global peak={int(global_peak['bin'])}")
    plt.xlabel("Range bin")
    plt.ylabel("Power (dB)")
    plt.title(f"Pulse-compressed range profile - {beam_name}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out1 = out_dir / f"{profile_path.stem}.png"
    plt.savefig(out1, dpi=160)
    plt.close()

    # Expected window zoom
    plt.figure(figsize=(12, 5))
    plt.plot(win["bin"], win["power_db"], marker="o", linewidth=1.0, markersize=3)
    plt.axvline(expected_bin, linestyle="--", linewidth=1.0, label=f"expected={expected_bin}")
    plt.axvline(local_peak["bin"], linestyle=":", linewidth=1.2, label=f"local peak={int(local_peak['bin'])}")
    plt.xlabel("Range bin")
    plt.ylabel("Power (dB)")
    plt.title(f"Expected-window zoom - {beam_name}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out2 = out_dir / f"{profile_path.stem}_expected_window.png"
    plt.savefig(out2, dpi=160)
    plt.close()

    print(f"\nSaved: {out1}")
    print(f"Saved: {out2}")

    return int(local_peak["bin"]), int(global_peak["bin"])


def plot_2d(amp2d_path: Path, expected_bin: int, window: int, out_dir: Path):
    df = pd.read_csv(amp2d_path)

    required = {"pulse_idx", "bin", "power_db"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{amp2d_path} missing columns: {missing}")

    df["pulse_idx"] = pd.to_numeric(df["pulse_idx"], errors="coerce")
    df["bin"] = pd.to_numeric(df["bin"], errors="coerce")
    df["power_db"] = pd.to_numeric(df["power_db"], errors="coerce")
    df = df.dropna(subset=["pulse_idx", "bin", "power_db"]).copy()

    df["pulse_idx"] = df["pulse_idx"].astype(int)
    df["bin"] = df["bin"].astype(int)

    left = expected_bin - window
    right = expected_bin + window
    win = df[(df["bin"] >= left) & (df["bin"] <= right)].copy()

    if len(win) == 0:
        print(f"[WARN] No 2D data in expected window [{left}, {right}]")
        return

    mat = win.pivot_table(
        index="pulse_idx",
        columns="bin",
        values="power_db",
        aggfunc="mean"
    )

    bins = mat.columns.to_numpy()
    pulses = mat.index.to_numpy()
    image = mat.to_numpy()

    beam_name = amp2d_path.stem.replace("pc_amplitude_2d_", "")

    plt.figure(figsize=(12, 6))
    plt.imshow(
        image,
        aspect="auto",
        origin="lower",
        extent=[bins.min(), bins.max(), pulses.min(), pulses.max()],
    )
    plt.axvline(expected_bin, linestyle="--", linewidth=1.0)
    plt.xlabel("Range bin")
    plt.ylabel("Pulse index")
    plt.title(f"Pulse-range power map near expected bin - {beam_name}")
    plt.colorbar(label="Power (dB)")
    plt.tight_layout()

    out = out_dir / f"{amp2d_path.stem}_expected_window.png"
    plt.savefig(out, dpi=160)
    plt.close()

    print(f"Saved: {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--debug-dir",
        default="outputs/p15_single_point/algorithm_result/debug",
        help="Directory containing pc_range_profile_beamXX.csv"
    )
    parser.add_argument("--beam-id", type=int, default=None)
    parser.add_argument("--expected-bin", type=int, default=880)
    parser.add_argument("--window", type=int, default=120)
    parser.add_argument("--plot-2d", action="store_true")
    args = parser.parse_args()

    debug_dir = Path(args.debug_dir)
    out_dir = debug_dir

    profile_path = find_profile(debug_dir, args.beam_id)
    plot_1d(profile_path, args.expected_bin, args.window, out_dir)

    if args.plot_2d:
        amp2d_path = find_amp2d(debug_dir, args.beam_id)
        if amp2d_path is None:
            print("[WARN] No pc_amplitude_2d_beamXX.csv found, skip 2D plot.")
        else:
            plot_2d(amp2d_path, args.expected_bin, args.window, out_dir)


if __name__ == "__main__":
    main()
