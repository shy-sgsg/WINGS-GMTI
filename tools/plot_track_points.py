#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "gmti_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "gmti_cache"))

import matplotlib.pyplot as plt
import pandas as pd


REQUIRED_COLUMNS = ["result_id", "track_id", "e", "n", "lat", "lon"]


def parse_args():
    parser = argparse.ArgumentParser(description="Plot GMTI track_points.csv.")
    parser.add_argument("--csv", required=True, help="Path to track_points.csv")
    parser.add_argument("--out", required=True, help="Output image path")
    parser.add_argument("--coord", choices=["en", "geo"], default="en",
                        help="Coordinate system to plot")
    parser.add_argument(
        "--min-len",
        type=int,
        default=3,
        help="Only plot tracks with at least this many points. Default: 3.",
    )
    return parser.parse_args()


def load_track_points(csv_path):
    if not os.path.exists(csv_path):
        raise ValueError(f"CSV does not exist: {csv_path}")
    if os.path.getsize(csv_path) == 0:
        raise ValueError(f"CSV is empty: {csv_path}")

    try:
        df = pd.read_csv(csv_path)
    except pd.errors.EmptyDataError:
        raise ValueError(f"CSV is empty or has no header: {csv_path}")

    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(
            "CSV is missing required columns: " + ", ".join(missing)
        )
    if df.empty:
        raise ValueError(f"CSV has no track point rows: {csv_path}")

    return df


def filter_short_tracks(df, min_len):
    track_lengths = df.groupby("track_id").size()
    keep_ids = track_lengths[track_lengths >= min_len].index

    print(f"Loaded {len(df)} points from {track_lengths.size} tracks")
    print(f"Keeping {len(keep_ids)} tracks with length >= {min_len}")
    print(f"Filtered out {track_lengths.size - len(keep_ids)} short tracks")

    return df[df["track_id"].isin(keep_ids)].copy()


def plot_track_points(df, out_path, coord, min_len):
    if coord == "en":
        xcol, ycol = "e", "n"
        xlabel, ylabel = "E (m)", "N (m)"
    else:
        xcol, ycol = "lon", "lat"
        xlabel, ylabel = "Longitude", "Latitude"

    fig, ax = plt.subplots(figsize=(12, 10))
    track_ids = sorted(df["track_id"].dropna().unique())
    if not track_ids:
        raise ValueError("CSV has no valid track_id values")

    cmap = plt.get_cmap("tab20")

    for i, tid in enumerate(track_ids):
        g = df[df["track_id"] == tid].sort_values("result_id")
        color = cmap(i % 20)
        ax.plot(
            g[xcol],
            g[ycol],
            marker="o",
            linewidth=0.5,
            markersize=1,
            color=color,
            label=f"ID {tid}",
        )

        # last = g.iloc[-1]
        # ax.text(last[xcol], last[ycol], f"ID {tid}", fontsize=8)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"GMTI track points, track length >= {min_len}")
    if coord == "en":
        ax.axis("equal")
    ax.grid(True)
    ax.legend(loc="best", fontsize=8)

    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    try:
        df = load_track_points(args.csv)
        df = filter_short_tracks(df, args.min_len)
        if df.empty:
            print(f"No tracks with length >= {args.min_len}. Nothing to plot.")
            return 0
        plot_track_points(df, args.out, args.coord, args.min_len)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
