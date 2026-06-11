#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
import os
import shutil
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "gmti_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "gmti_cache"))

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


STATE_MARKERS = {
    "Tentative": "o",
    "Confirmed": "s",
    "Coasted": "^",
    "Deleted": "x",
}

STATE_SIZES = {
    "Tentative": 25,
    "Confirmed": 45,
    "Coasted": 55,
    "Deleted": 70,
}

STATE_COLORS = {
    "Tentative": "#7f7f7f",
    "Confirmed": "#2ca02c",
    "Coasted": "#ff7f0e",
    "Deleted": "#d62728",
}


@dataclass
class PlotStyle:
    track_linewidth: float
    track_alpha: float
    marker_scale: float
    detection_marker_size: float
    detection_alpha: float
    output_marker_scale: float
    output_edge_width: float
    label_fontsize: float
    legend_fontsize: float
    title_fontsize: float
    axis_label_fontsize: float
    tick_fontsize: float


def parse_args():
    parser = argparse.ArgumentParser(description="Visualize TrackManager debug snapshots.")
    parser.add_argument("--debug-dir", required=True)
    parser.add_argument("--coord", choices=["en", "geo"], default="en")
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--format", choices=["png", "pdf", "svg"], default="svg")
    parser.add_argument("--result-id", type=int, default=None)
    parser.add_argument("--start-id", type=int, default=None)
    parser.add_argument("--end-id", type=int, default=None)
    parser.add_argument("--min-len", type=int, default=3)
    parser.add_argument("--show-detections", action="store_true", default=True)
    parser.add_argument("--hide-detections", action="store_true")
    parser.add_argument("--show-speed", action="store_true", default=True)
    parser.add_argument("--label-tracks", action="store_true", default=True)
    parser.add_argument("--no-label-tracks", action="store_true")
    parser.add_argument("--make-frames", action="store_true")
    parser.add_argument("--make-gif", action="store_true")
    parser.add_argument("--gif-duration-ms", type=int, default=500)
    parser.add_argument("--keep-frame-png", nargs="?", const="true", default="true",
                        choices=["true", "false", "1", "0"],
                        help="Keep PNG frames used for GIF. Default: true.")
    parser.add_argument("--xmin", type=float, default=None)
    parser.add_argument("--xmax", type=float, default=None)
    parser.add_argument("--ymin", type=float, default=None)
    parser.add_argument("--ymax", type=float, default=None)
    parser.add_argument("--roi-mode", choices=["any", "points"], default="any")
    parser.add_argument("--center-x", type=float, default=None)
    parser.add_argument("--center-y", type=float, default=None)
    parser.add_argument("--radius", type=float, default=None)
    parser.add_argument("--padding-ratio", type=float, default=0.05)
    parser.add_argument("--track-linewidth", type=float, default=None)
    parser.add_argument("--track-alpha", type=float, default=None)
    parser.add_argument("--marker-scale", type=float, default=None)
    parser.add_argument("--detection-marker-size", type=float, default=None)
    parser.add_argument("--detection-alpha", type=float, default=None)
    parser.add_argument("--output-marker-scale", type=float, default=None)
    parser.add_argument("--output-edge-width", type=float, default=None)
    parser.add_argument("--label-fontsize", type=float, default=None)
    parser.add_argument("--legend-fontsize", type=float, default=None)
    parser.add_argument("--title-fontsize", type=float, default=None)
    parser.add_argument("--axis-label-fontsize", type=float, default=None)
    parser.add_argument("--tick-fontsize", type=float, default=None)
    parser.add_argument("--hide-labels-in-frames", action="store_true")
    parser.add_argument("--max-label-tracks", type=int, default=80)
    parser.add_argument("--style-preset", choices=["default", "global", "roi", "dense"],
                        default=None)
    parser.add_argument("--fig-width", type=float, default=12)
    parser.add_argument("--fig-height", type=float, default=10)
    parser.add_argument("--frame-fig-width", type=float, default=10)
    parser.add_argument("--frame-fig-height", type=float, default=8)
    return parser.parse_args()


def bool_arg(value):
    return str(value).lower() in ("true", "1", "yes", "y")


def resolve_style(args, is_roi=False, is_frame=False, track_count=0):
    preset = args.style_preset
    if preset is None:
        if track_count > 100:
            preset = "dense"
            print("Visible track count exceeds 100; using dense style preset.")
        elif is_roi:
            preset = "roi"
        else:
            preset = "default"

    base = {
        "track_linewidth": 1.2,
        "track_alpha": 0.85,
        "marker_scale": 1.0,
        "detection_marker_size": 10.0,
        "detection_alpha": 0.35,
        "output_marker_scale": 1.5,
        "output_edge_width": 1.2,
        "label_fontsize": 8.0,
        "legend_fontsize": 8.0,
        "title_fontsize": 11.0,
        "axis_label_fontsize": 10.0,
        "tick_fontsize": 8.0,
    }

    if preset in ("default", "global"):
        pass
    elif preset == "roi":
        base.update({
            "track_linewidth": 0.8,
            "marker_scale": 0.75,
            "detection_marker_size": 6.0,
            "label_fontsize": 7.0,
            "legend_fontsize": 7.0,
        })
    elif preset == "dense":
        base.update({
            "track_linewidth": 0.6,
            "track_alpha": 0.75,
            "marker_scale": 0.55,
            "detection_marker_size": 4.0,
            "detection_alpha": 0.25,
            "label_fontsize": 6.0,
            "legend_fontsize": 6.0,
        })

    if is_frame and args.style_preset is None:
        base["track_linewidth"] = min(base["track_linewidth"], 0.8)
        base["marker_scale"] = min(base["marker_scale"], 0.8)
        base["label_fontsize"] = min(base["label_fontsize"], 7.0)
        base["legend_fontsize"] = min(base["legend_fontsize"], 7.0)

    overrides = {
        "track_linewidth": args.track_linewidth,
        "track_alpha": args.track_alpha,
        "marker_scale": args.marker_scale,
        "detection_marker_size": args.detection_marker_size,
        "detection_alpha": args.detection_alpha,
        "output_marker_scale": args.output_marker_scale,
        "output_edge_width": args.output_edge_width,
        "label_fontsize": args.label_fontsize,
        "legend_fontsize": args.legend_fontsize,
        "title_fontsize": args.title_fontsize,
        "axis_label_fontsize": args.axis_label_fontsize,
        "tick_fontsize": args.tick_fontsize,
    }
    for key, value in overrides.items():
        if value is not None:
            base[key] = value

    return PlotStyle(**base)


def read_csv_checked(path, required):
    if not os.path.exists(path):
        raise ValueError(f"missing CSV: {path}")
    if os.path.getsize(path) == 0:
        raise ValueError(f"empty CSV: {path}")
    df = pd.read_csv(path)
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {', '.join(missing)}")
    return df


def load_debug_csvs(debug_dir):
    frames = read_csv_checked(
        os.path.join(debug_dir, "track_frames.csv"),
        [
            "result_id", "num_detections", "num_tracks", "num_tentative",
            "num_confirmed", "num_coasted", "num_deleted", "num_outputs",
        ],
    )
    detections = read_csv_checked(
        os.path.join(debug_dir, "track_detections.csv"),
        [
            "result_id", "det_index", "matched", "matched_track_id",
            "e", "n", "lat", "lon",
        ],
    )
    states = read_csv_checked(
        os.path.join(debug_dir, "track_states.csv"),
        [
            "result_id", "track_id", "state", "matched_this_frame",
            "matched_det_index", "e", "n", "lat", "lon", "speed",
            "hit_count", "miss_count", "recent_hits", "is_output",
        ],
    )
    return frames, detections, states


def apply_period_filter(states, detections, frames, start_id, end_id, result_id):
    if result_id is not None:
        start_id = result_id
        end_id = result_id
    for_df = []
    for df in (states, detections, frames):
        out = df
        if start_id is not None:
            out = out[out["result_id"] >= start_id]
        if end_id is not None:
            out = out[out["result_id"] <= end_id]
        for_df.append(out.copy())
    return for_df[0], for_df[1], for_df[2]


def get_coord_columns(coord):
    if coord == "en":
        return "e", "n", "E (m)", "N (m)"
    return "lon", "lat", "Longitude", "Latitude"


def compute_roi_from_args(args):
    xmin, xmax, ymin, ymax = args.xmin, args.xmax, args.ymin, args.ymax
    if args.radius is not None:
        if args.coord != "en":
            print("--radius is currently intended for EN coordinates.", file=sys.stderr)
        elif args.center_x is None or args.center_y is None:
            raise ValueError("--radius requires --center-x and --center-y")
        else:
            xmin = args.center_x - args.radius
            xmax = args.center_x + args.radius
            ymin = args.center_y - args.radius
            ymax = args.center_y + args.radius
    if all(v is None for v in (xmin, xmax, ymin, ymax)):
        return None
    return {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax}


def in_roi(df, roi, xcol, ycol):
    mask = pd.Series(True, index=df.index)
    if roi is None:
        return mask
    if roi["xmin"] is not None:
        mask &= df[xcol] >= roi["xmin"]
    if roi["xmax"] is not None:
        mask &= df[xcol] <= roi["xmax"]
    if roi["ymin"] is not None:
        mask &= df[ycol] >= roi["ymin"]
    if roi["ymax"] is not None:
        mask &= df[ycol] <= roi["ymax"]
    return mask


def apply_roi_filter(states, detections, roi, coord, roi_mode):
    if roi is None:
        return states.copy(), detections.copy()
    xcol, ycol, _, _ = get_coord_columns(coord)
    state_mask = in_roi(states, roi, xcol, ycol)
    det_mask = in_roi(detections, roi, xcol, ycol)
    detections = detections[det_mask].copy()
    if roi_mode == "points":
        states = states[state_mask].copy()
    else:
        keep_ids = set(states.loc[state_mask, "track_id"].unique())
        states = states[states["track_id"].isin(keep_ids)].copy()
    return states, detections


def filter_by_min_len(states, min_len):
    matched = states[states["matched_this_frame"] == 1]
    lengths = matched.groupby("track_id").size()
    all_ids = sorted(states["track_id"].dropna().unique())
    keep_ids = set(lengths[lengths >= min_len].index)
    print(f"Loaded {len(states)} points from {len(all_ids)} tracks")
    print(f"Keeping {len(keep_ids)} tracks with length >= {min_len}")
    print(f"Filtered out {len(all_ids) - len(keep_ids)} short tracks")
    return states[states["track_id"].isin(keep_ids)].copy(), keep_ids, lengths


def get_track_colors(track_ids):
    cmap = plt.get_cmap("tab20")
    return {tid: cmap(i % 20) for i, tid in enumerate(sorted(track_ids))}


def roi_title(roi, coord):
    if roi is None:
        return ""
    xname, yname = ("E", "N") if coord == "en" else ("lon", "lat")
    def fmt(v):
        return "*" if v is None else f"{v:g}"
    return (
        f"ROI: {xname}=[{fmt(roi['xmin'])}, {fmt(roi['xmax'])}], "
        f"{yname}=[{fmt(roi['ymin'])}, {fmt(roi['ymax'])}]"
    )


def padded_limits(values, lo=None, hi=None, padding_ratio=0.05):
    series = pd.Series(values).dropna()
    if lo is None:
        lo = float(series.min()) if not series.empty else 0.0
    if hi is None:
        hi = float(series.max()) if not series.empty else lo + 1.0
    if hi == lo:
        pad = max(abs(lo) * padding_ratio, 1.0)
    else:
        pad = (hi - lo) * padding_ratio
    return lo - pad, hi + pad


def compute_axis_limits(states, detections, roi, coord, padding_ratio):
    xcol, ycol, _, _ = get_coord_columns(coord)
    xvals = pd.concat([states[xcol], detections[xcol]], ignore_index=True)
    yvals = pd.concat([states[ycol], detections[ycol]], ignore_index=True)
    xlim = padded_limits(xvals, None if roi is None else roi["xmin"],
                         None if roi is None else roi["xmax"], padding_ratio)
    ylim = padded_limits(yvals, None if roi is None else roi["ymin"],
                         None if roi is None else roi["ymax"], padding_ratio)
    return xlim, ylim


def apply_axes(ax, coord, style, xlim=None, ylim=None):
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if coord == "en":
        ax.set_aspect("equal", adjustable="box")
    ax.tick_params(labelsize=style.tick_fontsize)
    ax.grid(True, alpha=0.3)


def save_figure(fig, path, fmt):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if fmt == "png":
        fig.savefig(path, dpi=300, bbox_inches="tight")
    else:
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path}")


def save_named_figure(fig, out_dir, stem, fmt):
    save_figure(fig, os.path.join(out_dir, f"{stem}.{fmt}"), fmt)


def add_marker_legend(ax, style, include_tracks=True, include_track_id_legend=False):
    handles = [
        Line2D([0], [0], marker="o", color="#4c78a8", markerfacecolor="#4c78a8",
               markeredgecolor="#4c78a8", label="Tentative", linestyle="none"),
        Line2D([0], [0], marker="s", color="#4c78a8", markerfacecolor="#4c78a8",
               markeredgecolor="#4c78a8", label="Confirmed", linestyle="none"),
        Line2D([0], [0], marker="^", color="#4c78a8", markerfacecolor="#4c78a8",
               markeredgecolor="#4c78a8", label="Coasted", linestyle="none"),
        Line2D([0], [0], marker="x", color="#4c78a8", label="Deleted", linestyle="none"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor="none",
               label="black edge = final output", linestyle="none"),
        Line2D([0], [0], marker=".", color="#bdbdbd", label="gray points = detections", linestyle="none"),
        Line2D([0], [0], marker="o", color="#969696", markerfacecolor="none",
               label="unmatched detection", linestyle="none"),
    ]
    if include_tracks:
        handles.insert(0, Line2D([0], [0], color="#4c78a8", label="track color = track_id"))
        handles.insert(1, Line2D([0], [0], color="none", label="marker shape = track state"))
    if include_track_id_legend:
        handles.append(Line2D([0], [0], color="none", label="individual IDs use fixed colors"))
    ax.legend(handles=handles, loc="best", fontsize=style.legend_fontsize)


def plot_detections(ax, dets, xcol, ycol, style):
    if dets.empty:
        return
    matched = dets[dets["matched"] == 1]
    unmatched = dets[dets["matched"] != 1]
    if not matched.empty:
        ax.scatter(matched[xcol], matched[ycol], s=style.detection_marker_size,
                   c="#bdbdbd", alpha=style.detection_alpha, marker=".")
    if not unmatched.empty:
        ax.scatter(unmatched[xcol], unmatched[ycol],
                   s=style.detection_marker_size * 1.8,
                   facecolors="none", edgecolors="#737373",
                   alpha=min(1.0, style.detection_alpha + 0.35),
                   marker="o", linewidths=max(0.6, style.output_edge_width * 0.7))


def should_label_tracks(args, track_count, is_frame=False):
    if is_frame and args.hide_labels_in_frames:
        return False
    if not args.label_tracks:
        return False
    if track_count > args.max_label_tracks:
        return False
    return True


def marker_size(state, is_output, style):
    size = STATE_SIZES.get(state, 40) * style.marker_scale
    if is_output:
        size *= style.output_marker_scale
    return size


def scatter_state_points(ax, data, xcol, ycol, state, color, style):
    if data.empty:
        return
    marker = STATE_MARKERS.get(state, "o")
    alpha = 0.7 if state == "Coasted" else 0.85
    output = data[data["is_output"] == 1]
    normal = data[data["is_output"] != 1]
    if not normal.empty:
        if marker == "x":
            ax.scatter(normal[xcol], normal[ycol],
                       s=marker_size(state, False, style),
                       marker=marker, color=color, alpha=alpha,
                       linewidths=max(0.8, style.output_edge_width * 0.8))
        else:
            ax.scatter(normal[xcol], normal[ycol],
                       s=marker_size(state, False, style),
                       marker=marker, color=color, edgecolors=color,
                       alpha=alpha, linewidths=0.5)
    if not output.empty:
        ax.scatter(output[xcol], output[ycol],
                   s=marker_size(state, True, style),
                   marker=marker, color=color, edgecolors="black",
                   alpha=0.95, linewidths=style.output_edge_width)


def plot_track_paths(ax, states, colors, args, style, current_id=None,
                     label_ends=True, is_frame=False):
    xcol, ycol, _, _ = get_coord_columns(args.coord)
    track_count = states["track_id"].nunique()
    do_labels = should_label_tracks(args, track_count, is_frame=is_frame)
    for tid in sorted(states["track_id"].unique()):
        g = states[states["track_id"] == tid].sort_values("result_id")
        color = colors.get(tid, "#4c78a8")
        ax.plot(g[xcol], g[ycol], color=color,
                linewidth=style.track_linewidth, alpha=style.track_alpha)
        for state, marker in STATE_MARKERS.items():
            s = g[g["state"] == state]
            scatter_state_points(ax, s, xcol, ycol, state, color, style)
        if current_id is not None:
            c = g[g["result_id"] == current_id]
            if not c.empty:
                row = c.iloc[-1]
                state = row["state"]
                marker = STATE_MARKERS.get(state, "o")
                is_output = row["is_output"] == 1
                if marker == "x":
                    ax.scatter([row[xcol]], [row[ycol]],
                               s=marker_size(state, is_output, style) * 1.25,
                               marker=marker, color="black" if is_output else color,
                               linewidths=style.output_edge_width if is_output else 0.7,
                               alpha=0.95 if is_output else 0.85)
                else:
                    ax.scatter([row[xcol]], [row[ycol]],
                               s=marker_size(state, is_output, style) * 1.25,
                               marker=marker, color=color,
                               edgecolors="black" if is_output else color,
                               linewidths=style.output_edge_width if is_output else 0.7,
                               alpha=0.95 if is_output else (0.7 if state == "Coasted" else 0.85))
                label = f"ID={tid}"
                if args.show_speed:
                    label += f" v={row['speed']:.1f} m/s"
                if do_labels:
                    ax.text(row[xcol], row[ycol], label, fontsize=style.label_fontsize)
        elif label_ends and do_labels:
            last = g.iloc[-1]
            ax.text(last[xcol], last[ycol], f"ID={tid}", fontsize=style.label_fontsize)


def title_suffix(args, roi):
    pieces = [f"min_len >= {args.min_len}"]
    roi_text = roi_title(roi, args.coord)
    if roi_text:
        pieces.append(roi_text)
    return "\n".join(pieces)


def plot_global_paths(states, detections, frames, args, out_dir, colors, roi, axis_limits, style):
    xcol, ycol, xlabel, ylabel = get_coord_columns(args.coord)
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
    if args.show_detections:
        plot_detections(ax, detections, xcol, ycol, style)
    plot_track_paths(ax, states, colors, args, style)
    ax.set_xlabel(xlabel, fontsize=style.axis_label_fontsize)
    ax.set_ylabel(ylabel, fontsize=style.axis_label_fontsize)
    rid_min = int(states["result_id"].min())
    rid_max = int(states["result_id"].max())
    ax.set_title(
        f"Track global paths result_id {rid_min}-{rid_max}\n{title_suffix(args, roi)}",
        fontsize=style.title_fontsize,
    )
    apply_axes(ax, args.coord, style, *axis_limits)
    add_marker_legend(ax, style)
    save_named_figure(fig, out_dir, "track_global_paths", args.format)


def frame_stats_text(frames, rid):
    row = frames[frames["result_id"] == rid]
    if row.empty:
        return "dets=?, outputs=?, Tentative=?, Confirmed=?, Coasted=?, Deleted=?"
    r = row.iloc[-1]
    return (
        f"dets={int(r['num_detections'])}, outputs={int(r['num_outputs'])}, "
        f"Tentative={int(r['num_tentative'])}, Confirmed={int(r['num_confirmed'])}, "
        f"Coasted={int(r['num_coasted'])}, Deleted={int(r['num_deleted'])}"
    )


def plot_snapshot(states, detections, frames, args, out_dir, colors, roi, axis_limits,
                  style, result_id=None, stem="track_snapshot_last", fmt=None,
                  is_frame=False):
    xcol, ycol, xlabel, ylabel = get_coord_columns(args.coord)
    rid = int(result_id if result_id is not None else frames["result_id"].max())
    hist = states[states["result_id"] <= rid]
    current_dets = detections[detections["result_id"] == rid]

    if is_frame:
        figsize = (args.frame_fig_width, args.frame_fig_height)
    else:
        figsize = (args.fig_width, args.fig_height)
    fig, ax = plt.subplots(figsize=figsize)
    if args.show_detections:
        plot_detections(ax, current_dets, xcol, ycol, style)
    plot_track_paths(ax, hist, colors, args, style, current_id=rid,
                     label_ends=False, is_frame=is_frame)

    ax.set_xlabel(xlabel, fontsize=style.axis_label_fontsize)
    ax.set_ylabel(ylabel, fontsize=style.axis_label_fontsize)
    if is_frame:
        title = f"TrackManager snapshot, result {rid}\n{frame_stats_text(frames, rid)}"
    else:
        title = f"TrackManager snapshot, result {rid}\n{frame_stats_text(frames, rid)}\n{title_suffix(args, roi)}"
    ax.set_title(title, fontsize=style.title_fontsize)
    apply_axes(ax, args.coord, style, *axis_limits)
    add_marker_legend(ax, style)
    save_named_figure(fig, out_dir, stem, fmt or args.format)


def plot_state_timeline(states, args, out_dir, style):
    fig, ax = plt.subplots(figsize=(args.fig_width, max(6, args.fig_height * 0.8)))
    for state, color in STATE_COLORS.items():
        s = states[states["state"] == state]
        if s.empty:
            continue
        ax.scatter(s["result_id"], s["track_id"],
                   s=STATE_SIZES.get(state, 35) * style.marker_scale,
                   marker=STATE_MARKERS.get(state, "o"), color=color, label=state)
    out = states[states["is_output"] == 1]
    if not out.empty:
        ax.scatter(out["result_id"], out["track_id"],
                   s=STATE_SIZES["Confirmed"] * style.marker_scale * style.output_marker_scale,
                   facecolors="none", edgecolors="black",
                   linewidths=style.output_edge_width, label="output")
    ax.set_xlabel("result_id", fontsize=style.axis_label_fontsize)
    ax.set_ylabel("track_id", fontsize=style.axis_label_fontsize)
    ax.set_title(f"Track state timeline, min_len >= {args.min_len}",
                 fontsize=style.title_fontsize)
    ax.tick_params(labelsize=style.tick_fontsize)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=style.legend_fontsize)
    save_named_figure(fig, out_dir, "track_state_timeline", args.format)


def plot_speed_timeline(states, args, out_dir, colors, style):
    matched = states[states["matched_this_frame"] == 1]
    fig, ax = plt.subplots(figsize=(args.fig_width, max(6, args.fig_height * 0.8)))
    for tid in sorted(matched["track_id"].unique()):
        g = matched[matched["track_id"] == tid].sort_values("result_id")
        ax.plot(g["result_id"], g["speed"], marker="o",
                linewidth=style.track_linewidth, alpha=style.track_alpha,
                markersize=max(2.0, 3.0 * style.marker_scale),
                color=colors.get(tid), label=f"ID {tid}")
    ax.set_xlabel("result_id", fontsize=style.axis_label_fontsize)
    ax.set_ylabel("speed", fontsize=style.axis_label_fontsize)
    ax.set_title(f"Matched-track speed timeline, min_len >= {args.min_len}",
                 fontsize=style.title_fontsize)
    ax.tick_params(labelsize=style.tick_fontsize)
    ax.grid(True, alpha=0.3)
    if matched["track_id"].nunique() <= 20:
        ax.legend(loc="best", fontsize=style.legend_fontsize)
    save_named_figure(fig, out_dir, "track_speed_timeline", args.format)


def plot_short_track_summary(all_states, args, out_dir, style):
    matched_len = all_states[all_states["matched_this_frame"] == 1].groupby("track_id").size()
    output_len = all_states[all_states["is_output"] == 1].groupby("track_id").size()
    life_len = all_states.groupby("track_id").size()
    all_ids = sorted(all_states["track_id"].unique())
    matched_len = matched_len.reindex(all_ids, fill_value=0)
    output_len = output_len.reindex(all_ids, fill_value=0)
    life_len = life_len.reindex(all_ids, fill_value=0)
    total = len(all_ids)
    bins = {
        "1": int((matched_len == 1).sum()),
        "2": int((matched_len == 2).sum()),
        "3-5": int(((matched_len >= 3) & (matched_len <= 5)).sum()),
        "6-10": int(((matched_len >= 6) & (matched_len <= 10)).sum()),
        ">10": int((matched_len > 10).sum()),
    }
    summary = [
        f"Total tracks: {total}",
        f"Tracks len=1: {bins['1']}",
        f"Tracks len=2: {bins['2']}",
        f"Tracks len>=3: {int((matched_len >= 3).sum())}",
        f"Median length: {matched_len.median() if total else 0:g}",
        f"Max length: {matched_len.max() if total else 0:g}",
        f"Median lifecycle: {life_len.median() if total else 0:g}",
        f"Max output hits: {output_len.max() if total else 0:g}",
    ]
    print("\n".join(summary))
    fig, ax = plt.subplots(figsize=(max(8, args.fig_width * 0.8), max(6, args.fig_height * 0.7)))
    ax.bar(list(bins.keys()), list(bins.values()), color="#4c78a8")
    ax.set_xlabel("matched_this_frame count", fontsize=style.axis_label_fontsize)
    ax.set_ylabel("track count", fontsize=style.axis_label_fontsize)
    ax.set_title("Short track summary", fontsize=style.title_fontsize)
    ax.tick_params(labelsize=style.tick_fontsize)
    ax.grid(axis="y", alpha=0.3)
    ax.text(0.98, 0.95, "\n".join(summary), transform=ax.transAxes,
            va="top", ha="right", fontsize=style.label_fontsize,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "#cccccc"})
    save_named_figure(fig, out_dir, "track_short_track_summary", args.format)


def make_frame_sequence(states, detections, frames, args, out_dir, colors, roi, axis_limits,
                        style, fmt=None, frame_dir=None):
    frame_fmt = fmt or args.format
    frame_dir = frame_dir or os.path.join(out_dir, "frames")
    os.makedirs(frame_dir, exist_ok=True)
    result_ids = sorted(frames["result_id"].unique())
    paths = []
    for rid in result_ids:
        stem = f"frame_{int(rid):06d}"
        plot_snapshot(states, detections, frames, args, frame_dir, colors, roi,
                      axis_limits, style, result_id=int(rid), stem=stem,
                      fmt=frame_fmt, is_frame=True)
        paths.append(os.path.join(frame_dir, f"{stem}.{frame_fmt}"))
    return paths


def make_gif_from_frames(frame_paths, gif_path, duration_ms):
    try:
        from PIL import Image
    except ImportError:
        try:
            import imageio.v2 as imageio
        except ImportError:
            print("GIF generation requires imageio. Please install it with: pip install imageio",
                  file=sys.stderr)
            return False
        try:
            images = [imageio.imread(path) for path in frame_paths]
            imageio.mimsave(gif_path, images, duration=duration_ms / 1000.0)
            print(f"Wrote {gif_path}")
            return True
        except Exception as exc:
            print(f"GIF generation failed: {exc}", file=sys.stderr)
            return False

    try:
        images = [Image.open(path).convert("P", palette=Image.ADAPTIVE) for path in frame_paths]
        if not images:
            print("No PNG frames available for GIF generation.", file=sys.stderr)
            return False
        os.makedirs(os.path.dirname(gif_path), exist_ok=True)
        images[0].save(gif_path, save_all=True, append_images=images[1:],
                       duration=duration_ms, loop=0)
        for image in images:
            image.close()
        print(f"Wrote {gif_path}")
        return True
    except Exception as exc:
        print(f"GIF generation failed: {exc}", file=sys.stderr)
        return False


def main():
    args = parse_args()
    args.show_detections = args.show_detections and not args.hide_detections
    args.label_tracks = args.label_tracks and not args.no_label_tracks
    keep_frame_png = bool_arg(args.keep_frame_png)
    out_dir = args.out_dir or os.path.join(args.debug_dir, "figs")

    try:
        frames, detections, states = load_debug_csvs(args.debug_dir)
        states, detections, frames = apply_period_filter(
            states, detections, frames, args.start_id, args.end_id, args.result_id)
        if states.empty:
            print("No track states in the selected result_id range.")
            return 0

        roi = compute_roi_from_args(args)
        states, detections = apply_roi_filter(states, detections, roi, args.coord, args.roi_mode)
        if states.empty:
            print("No tracks remain after ROI filtering.")
            return 0

        filtered_states, keep_ids, _ = filter_by_min_len(states, args.min_len)
        if filtered_states.empty:
            print(f"No tracks with length >= {args.min_len}. Nothing to plot.")
            return 0
        filtered_detections = detections.copy()
        track_count = filtered_states["track_id"].nunique()
        if args.label_tracks and track_count > args.max_label_tracks:
            print(
                f"Visible track count {track_count} exceeds --max-label-tracks "
                f"{args.max_label_tracks}; disabling track ID labels."
            )
            args.label_tracks = False
        is_roi_view = roi is not None
        style = resolve_style(args, is_roi=is_roi_view, is_frame=False,
                              track_count=track_count)
        frame_style = resolve_style(args, is_roi=is_roi_view, is_frame=True,
                                    track_count=track_count)
        colors = get_track_colors(filtered_states["track_id"].unique())
        axis_limits = compute_axis_limits(
            filtered_states, filtered_detections, roi, args.coord, args.padding_ratio)

        plot_global_paths(filtered_states, filtered_detections, frames, args, out_dir,
                          colors, roi, axis_limits, style)
        plot_snapshot(filtered_states, filtered_detections, frames, args, out_dir,
                      colors, roi, axis_limits, style)
        plot_state_timeline(filtered_states, args, out_dir, style)
        plot_speed_timeline(filtered_states, args, out_dir, colors, style)
        plot_short_track_summary(states, args, out_dir, style)

        frame_paths = None
        if args.make_frames:
            frame_paths = make_frame_sequence(filtered_states, filtered_detections, frames, args, out_dir,
                                             colors, roi, axis_limits, frame_style)

        if args.make_gif:
            gif_png_dir = os.path.join(out_dir, "frames")
            temp_dir = None
            if args.format == "png" and frame_paths is not None:
                png_paths = frame_paths
            else:
                temp_dir = os.path.join(out_dir, "_gif_png_frames")
                gif_png_dir = temp_dir
                png_paths = make_frame_sequence(filtered_states, filtered_detections, frames, args,
                                                out_dir, colors, roi, axis_limits,
                                                frame_style, fmt="png", frame_dir=gif_png_dir)
            made = make_gif_from_frames(
                png_paths, os.path.join(out_dir, "track_growth.gif"), args.gif_duration_ms)
            if made and not keep_frame_png and temp_dir:
                shutil.rmtree(temp_dir, ignore_errors=True)
            elif made and not keep_frame_png and args.format == "png" and frame_paths is None:
                for path in png_paths:
                    try:
                        os.remove(path)
                    except OSError:
                        pass
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
