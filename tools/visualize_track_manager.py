#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
import os
import re
import shutil
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "gmti_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "gmti_cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
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
    parser.add_argument("--debug-dir", default="/home/shy/AIR/小长/GMTI_Data/Mission068/result_final/track_debug")
    parser.add_argument("--coord", choices=["en", "geo"], default="geo")
    parser.add_argument("--out-dir", default="/home/shy/AIR/小长/GMTI_Data/Mission068/result_final/track_debug_out")
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
    parser.add_argument("--make-gif", action="store_true", default=True)
    parser.add_argument("--make-output-gif", action="store_true")
    parser.add_argument("--gif-duration-ms", type=int, default=500)
    parser.add_argument("--keep-frame-png", nargs="?", const="true", default="false",
                        choices=["true", "false", "1", "0"],
                        help="Keep PNG frames used for GIF. Default: false.")
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
    parser.add_argument("--fig-width", type=float, default=12)
    parser.add_argument("--fig-height", type=float, default=10)
    parser.add_argument("--frame-fig-width", type=float, default=10)
    parser.add_argument("--frame-fig-height", type=float, default=8)
    parser.add_argument("--show-dbs-bg", action="store_true", default=True,
                        help="Overlay DBS image as background on spatial track figures.")
    parser.add_argument("--no-dbs-bg", action="store_true",
                        help="Disable DBS background overlay.")
    parser.add_argument("--dbs-image", default="/home/shy/AIR/小长/GMTI_Data/Mission068/result_final/GMTI13.png",
                        help="Path to DBS PNG image.")
    parser.add_argument("--dbs-corner-file", default="/home/shy/AIR/小长/GMTI_Data/Mission068/result_final/GMTI13.txt",
                        help="Path to DBS corner lon/lat txt file.")
    parser.add_argument(
        "--dbs-bg-mode",
        choices=["warp", "image"],
        default="image",
        help=(
            "DBS background mode. "
            "warp: warp DBS image to regular lon/lat grid. "
            "image: keep original DBS PNG unchanged and project lon/lat points to image pixels."
        ),
    )
    parser.add_argument("--dbs-alpha", type=float, default=1.0,
                        help="DBS background alpha.")
    parser.add_argument(
        "--dbs-color-mode",
        choices=["original", "gray"],
        default="original",
        help=(
            "DBS background color mode. "
            "original: keep original PNG colors. "
            "gray: convert DBS background to grayscale."
        ),
    )
    parser.add_argument("--dbs-warp-width", type=int, default=4096,
                        help="Output width for warp mode.")
    parser.add_argument("--dbs-warp-height", type=int, default=6144,
                        help="Output height for warp mode.")
    parser.add_argument("--dbs-show-pixel-axis", action="store_true", default=False,
                        help="In image mode, show DBS image pixel axes.")
    parser.add_argument(
        "--dbs-native-resolution",
        action="store_true",
        default=True,
        help=(
            "In DBS image mode, save spatial PNG outputs at the original DBS PNG "
            "pixel resolution. This also applies to PNG frames used for GIF."
        ),
    )
    parser.add_argument(
        "--dbs-native-dpi",
        type=int,
        default=100,
        help="DPI used to map original DBS pixel size to Matplotlib figure inches.",
    )
    parser.add_argument("--dbs-debug", action="store_true", default=False,
                        help="Write DBS transform debug figures.")
    parser.add_argument("--dbs-clip-to-image", action="store_true", default=True,
                        help="In image mode, clip transformed points outside DBS image bounds.")
    parser.add_argument("--no-dbs-clip-to-image", action="store_true")
    parser.add_argument(
        "--dbs-image-corner-mode",
        choices=["direct", "cpp_transpose_fliplr"],
        default="cpp_transpose_fliplr",
        help=(
            "Corner mapping for DBS image/warp modes. "
            "direct maps txt logical corners directly to PNG corners. "
            "cpp_transpose_fliplr considers the C++ transpose_copy + fliplr_inplace before PNG output."
        ),
    )
    parser.add_argument("--axis-from", choices=["tracks", "dbs", "union"], default="tracks")
    return parser.parse_args()


def bool_arg(value):
    return str(value).lower() in ("true", "1", "yes", "y")


def resolve_style(args):
    base = {
        "track_linewidth": 0.8,
        "track_alpha": 0.8,
        "marker_scale": 0.5,
        "detection_marker_size": 3.0,
        "detection_alpha": 0.8,
        "output_marker_scale": 1.2,
        "output_edge_width": 0.6,
        "label_fontsize": 7.0,
        "legend_fontsize": 7.0,
        "title_fontsize": 12.0,
        "axis_label_fontsize": 10.0,
        "tick_fontsize": 7.0,
    }

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


def parse_dbs_corner_txt(path):
    vals = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            match = re.match(r"(\w+)\s*=\s*([+-]?\d+(?:\.\d+)?)", line)
            if match:
                vals[match.group(1)] = float(match.group(2))

    required = ["B0", "B1", "B2", "B3", "L0", "L1", "L2", "L3"]
    missing = [key for key in required if key not in vals]
    if missing:
        raise ValueError(f"{path} missing DBS corner keys: {missing}")

    corners = {
        0: {"lon": vals["L0"], "lat": vals["B0"]},
        1: {"lon": vals["L1"], "lat": vals["B1"]},
        2: {"lon": vals["L2"], "lat": vals["B2"]},
        3: {"lon": vals["L3"], "lat": vals["B3"]},
    }

    print("Loaded DBS corners:")
    print("  corner 0: top-left     lon={}, lat={}".format(corners[0]["lon"], corners[0]["lat"]))
    print("  corner 1: bottom-right lon={}, lat={}".format(corners[1]["lon"], corners[1]["lat"]))
    print("  corner 2: top-right    lon={}, lat={}".format(corners[2]["lon"], corners[2]["lat"]))
    print("  corner 3: bottom-left  lon={}, lat={}".format(corners[3]["lon"], corners[3]["lat"]))

    return corners


def read_dbs_image(path):
    try:
        import matplotlib.image as mpimg
        return mpimg.imread(path)
    except Exception as exc:
        raise ValueError(f"failed to read DBS image {path}: {exc}") from exc


def dbs_extent_from_corners(corners):
    lons = [corners[idx]["lon"] for idx in [0, 1, 2, 3]]
    lats = [corners[idx]["lat"] for idx in [0, 1, 2, 3]]
    return [
        float(min(lons)),
        float(max(lons)),
        float(min(lats)),
        float(max(lats)),
    ]


def dbs_geo_corners_image_order(corners, image_corner_mode="direct"):
    if image_corner_mode == "cpp_transpose_fliplr":
        return np.float32([
            [corners[3]["lon"], corners[3]["lat"]],
            [corners[0]["lon"], corners[0]["lat"]],
            [corners[2]["lon"], corners[2]["lat"]],
            [corners[1]["lon"], corners[1]["lat"]],
        ])
    return np.float32([
        [corners[0]["lon"], corners[0]["lat"]],
        [corners[2]["lon"], corners[2]["lat"]],
        [corners[1]["lon"], corners[1]["lat"]],
        [corners[3]["lon"], corners[3]["lat"]],
    ])


def build_geo_to_image_transform(corners, image_shape, image_corner_mode="direct"):
    try:
        import cv2
    except ImportError as exc:
        raise RuntimeError(
            "DBS image mode requires OpenCV. Please install it with: "
            "pip install opencv-python"
        ) from exc

    h, w = image_shape[:2]
    src_geo = dbs_geo_corners_image_order(corners, image_corner_mode)
    dst_pix = np.float32([
        [0, 0],
        [w - 1, 0],
        [w - 1, h - 1],
        [0, h - 1],
    ])
    matrix = cv2.getPerspectiveTransform(src_geo, dst_pix)

    print(f"DBS image mode geo->pixel transform, corner_mode={image_corner_mode}:")
    print("src_geo [lon, lat] =")
    print(src_geo)
    print("dst_pix [x, y] =")
    print(dst_pix)
    print("H_geo_to_pix =")
    print(matrix)
    return matrix, src_geo, dst_pix


def transform_lonlat_to_pixel(lon, lat, matrix):
    lon = np.asarray(lon, dtype=np.float64)
    lat = np.asarray(lat, dtype=np.float64)
    ones = np.ones_like(lon)
    pts = np.vstack([lon, lat, ones])
    out = matrix @ pts
    x = out[0, :] / out[2, :]
    y = out[1, :] / out[2, :]
    return x, y


def apply_dbs_color_mode(img, args):
    if args.dbs_color_mode == "original":
        return img

    if args.dbs_color_mode == "gray":
        arr = np.asarray(img)
        if arr.ndim == 2:
            return arr
        if arr.ndim == 3:
            rgb = arr[..., :3].astype(np.float32)
            if rgb.max() > 1.0:
                rgb = rgb / 255.0
            gray = (
                0.299 * rgb[..., 0] +
                0.587 * rgb[..., 1] +
                0.114 * rgb[..., 2]
            )
            return gray
        raise ValueError(f"Unsupported DBS image shape for gray mode: {arr.shape}")

    raise ValueError(f"Unsupported dbs_color_mode: {args.dbs_color_mode}")


def build_dbs_background_image(args):
    if args.coord != "geo":
        raise ValueError(
            "DBS image mode requires --coord geo, because lon/lat are needed "
            "to project tracks and detections into DBS image pixels."
        )
    img_raw = plt.imread(args.dbs_image)
    h, w = img_raw.shape[:2]
    img = apply_dbs_color_mode(img_raw, args)
    corners = parse_dbs_corner_txt(args.dbs_corner_file)
    matrix, src_geo, dst_pix = build_geo_to_image_transform(
        corners, img_raw.shape, args.dbs_image_corner_mode)
    return {
        "mode": "image",
        "image": img,
        "original_image": img_raw,
        "width": w,
        "height": h,
        "extent": [0, w - 1, h - 1, 0],
        "origin": "upper",
        "alpha": args.dbs_alpha,
        "color_mode": args.dbs_color_mode,
        "H_geo_to_pix": matrix,
        "src_geo": src_geo,
        "dst_pix": dst_pix,
        "image_corner_mode": args.dbs_image_corner_mode,
        "corners": corners,
    }


def build_dbs_background_warp(args):
    try:
        import cv2
    except ImportError as exc:
        raise RuntimeError(
            "DBS warp mode requires OpenCV. Please install it with: "
            "pip install opencv-python"
        ) from exc

    img_raw = plt.imread(args.dbs_image)
    h, w = img_raw.shape[:2]
    img = apply_dbs_color_mode(img_raw, args)
    corners = parse_dbs_corner_txt(args.dbs_corner_file)

    src_pix = np.float32([
        [0, 0],
        [w - 1, 0],
        [w - 1, h - 1],
        [0, h - 1],
    ])
    geo = dbs_geo_corners_image_order(corners, args.dbs_image_corner_mode)
    lon_min = float(geo[:, 0].min())
    lon_max = float(geo[:, 0].max())
    lat_min = float(geo[:, 1].min())
    lat_max = float(geo[:, 1].max())
    if lon_max == lon_min or lat_max == lat_min:
        raise ValueError("DBS corner extent has zero width or height")

    out_w = int(args.dbs_warp_width)
    out_h = int(args.dbs_warp_height)
    dst_pix = np.zeros((4, 2), dtype=np.float32)
    dst_pix[:, 0] = (geo[:, 0] - lon_min) / (lon_max - lon_min) * (out_w - 1)
    dst_pix[:, 1] = (lat_max - geo[:, 1]) / (lat_max - lat_min) * (out_h - 1)
    matrix = cv2.getPerspectiveTransform(src_pix, dst_pix)
    warped = cv2.warpPerspective(np.asarray(img), matrix, (out_w, out_h))

    print(f"DBS warp mode pixel->regular-lonlat transform, corner_mode={args.dbs_image_corner_mode}:")
    print("src_pix [x, y] =")
    print(src_pix)
    print("geo in PNG corner order [lon, lat] =")
    print(geo)
    print("dst_pix on regular lon/lat grid [x, y] =")
    print(dst_pix)
    print("H_pix_to_warp =")
    print(matrix)

    return {
        "mode": "warp",
        "image": warped,
        "original_image": img_raw,
        "width": out_w,
        "height": out_h,
        "extent": [lon_min, lon_max, lat_min, lat_max],
        "origin": "upper",
        "alpha": args.dbs_alpha,
        "color_mode": args.dbs_color_mode,
        "src_pix": src_pix,
        "dst_pix": dst_pix,
        "H_pix_to_warp": matrix,
        "geo": geo,
        "image_corner_mode": args.dbs_image_corner_mode,
        "corners": corners,
    }


def build_dbs_background(args):
    if not args.show_dbs_bg:
        return None
    if args.coord != "geo":
        if args.dbs_bg_mode == "image":
            raise ValueError(
                "DBS image mode requires --coord geo because tracks/detections must "
                "be mapped from lon/lat to DBS image pixels."
            )
        raise ValueError(
            "DBS background currently requires --coord geo because the DBS corner file "
            "is given in lon/lat. Please run with --coord geo, or add lon/lat to EN "
            "projection conversion later."
        )
    if not args.dbs_image:
        raise ValueError("--show-dbs-bg requires --dbs-image")
    if not args.dbs_corner_file:
        raise ValueError("--show-dbs-bg requires --dbs-corner-file")
    if not os.path.exists(args.dbs_image):
        raise ValueError(f"DBS image not found: {args.dbs_image}")
    if not os.path.exists(args.dbs_corner_file):
        raise ValueError(f"DBS corner file not found: {args.dbs_corner_file}")
    if args.dbs_bg_mode == "image":
        return build_dbs_background_image(args)
    if args.dbs_bg_mode == "warp":
        return build_dbs_background_warp(args)
    raise ValueError(f"Unsupported DBS mode: {args.dbs_bg_mode}")


def add_dbs_background(ax, dbs_bg):
    if dbs_bg is None:
        return
    origin = "upper" if dbs_bg["mode"] == "image" else dbs_bg["origin"]
    kwargs = {
        "extent": dbs_bg["extent"],
        "origin": origin,
        "alpha": dbs_bg["alpha"],
        "zorder": 0,
    }
    if dbs_bg.get("color_mode") == "gray":
        kwargs["cmap"] = "gray"
        kwargs["vmin"] = 0.0
        kwargs["vmax"] = 1.0
    ax.imshow(dbs_bg["image"], **kwargs)


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


def get_plot_coord_columns(args):
    if args.show_dbs_bg and args.dbs_bg_mode == "image":
        return "dbs_x", "dbs_y", "DBS image x (pixel)", "DBS image y (pixel)"
    return get_coord_columns(args.coord)


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


def compute_axis_limits_with_dbs(states, detections, roi, coord, padding_ratio, dbs_bg, axis_from):
    if dbs_bg is not None and dbs_bg.get("mode") == "image":
        return (0, dbs_bg["width"] - 1), (dbs_bg["height"] - 1, 0)
    if dbs_bg is not None and dbs_bg.get("mode") == "warp":
        extent = dbs_bg["extent"]
        return (extent[0], extent[1]), (extent[2], extent[3])
    if dbs_bg is None or axis_from == "tracks":
        return compute_axis_limits(states, detections, roi, coord, padding_ratio)
    extent = dbs_bg["extent"]
    xcol, ycol, _, _ = get_coord_columns(coord)
    xvals = pd.concat([
        states[xcol],
        detections[xcol],
        pd.Series([extent[0], extent[1]]),
    ], ignore_index=True)
    yvals = pd.concat([
        states[ycol],
        detections[ycol],
        pd.Series([extent[2], extent[3]]),
    ], ignore_index=True)
    return padded_limits(xvals, None, None, padding_ratio), padded_limits(yvals, None, None, padding_ratio)


def apply_axes(ax, args, style, xlim=None, ylim=None):
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    image_mode = args.show_dbs_bg and args.dbs_bg_mode == "image"
    if args.coord == "en" or image_mode:
        ax.set_aspect("equal", adjustable="box")
    ax.tick_params(labelsize=style.tick_fontsize)
    ax.grid(True, alpha=0.3)
    if image_mode and not args.dbs_show_pixel_axis:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)


def add_dbs_pixel_columns(states, detections, dbs_bg):
    matrix = dbs_bg["H_geo_to_pix"]
    states = states.copy()
    detections = detections.copy()
    states["dbs_x"], states["dbs_y"] = transform_lonlat_to_pixel(
        states["lon"], states["lat"], matrix)
    detections["dbs_x"], detections["dbs_y"] = transform_lonlat_to_pixel(
        detections["lon"], detections["lat"], matrix)
    return states, detections


def clip_to_dbs_image(df, dbs_bg, xcol="dbs_x", ycol="dbs_y"):
    if df.empty:
        return df
    w = dbs_bg["width"]
    h = dbs_bg["height"]
    mask = (
        (df[xcol] >= 0) &
        (df[xcol] <= w - 1) &
        (df[ycol] >= 0) &
        (df[ycol] <= h - 1)
    )
    return df[mask].copy()


def use_dbs_native_resolution(args, dbs_bg, fmt):
    return (
        args.dbs_native_resolution
        and dbs_bg is not None
        and dbs_bg.get("mode") == "image"
        and fmt == "png"
    )


def create_spatial_figure(args, dbs_bg=None, is_frame=False, fmt=None):
    out_fmt = fmt or args.format
    if use_dbs_native_resolution(args, dbs_bg, out_fmt):
        dpi = max(1, int(args.dbs_native_dpi))
        width = int(dbs_bg["width"])
        height = int(dbs_bg["height"])
        fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.set_position([0, 0, 1, 1])
        fig._dbs_native_png = True
        fig._dbs_native_dpi = dpi
        return fig, ax

    if is_frame:
        figsize = (args.frame_fig_width, args.frame_fig_height)
    else:
        figsize = (args.fig_width, args.fig_height)
    fig, ax = plt.subplots(figsize=figsize)
    fig._dbs_native_png = False
    return fig, ax


def save_figure(fig, path, fmt):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if fmt == "png":
        if getattr(fig, "_dbs_native_png", False):
            fig.savefig(
                path,
                dpi=getattr(fig, "_dbs_native_dpi", 100),
                bbox_inches=None,
                pad_inches=0,
            )
        else:
            fig.savefig(path, dpi=300, bbox_inches="tight")
    else:
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path}")


def save_named_figure(fig, out_dir, stem, fmt):
    save_figure(fig, os.path.join(out_dir, f"{stem}.{fmt}"), fmt)


def corner_pixels_for_dbs_image(image_shape, image_corner_mode="direct"):
    h, w = image_shape[:2]
    if image_corner_mode == "cpp_transpose_fliplr":
        return {
            3: (0, 0),
            0: (w - 1, 0),
            2: (w - 1, h - 1),
            1: (0, h - 1),
        }
    return {
        0: (0, 0),
        2: (w - 1, 0),
        1: (w - 1, h - 1),
        3: (0, h - 1),
    }


def save_dbs_debug_figures(args, img, corners, dbs_bg, out_dir, detections, states, style):
    if not args.dbs_debug or dbs_bg is None:
        return
    debug_dir = os.path.join(out_dir, "dbs_debug")
    os.makedirs(debug_dir, exist_ok=True)

    image_corner_mode = dbs_bg.get("image_corner_mode", "direct")
    corner_pixels = corner_pixels_for_dbs_image(img.shape, image_corner_mode)
    display_img = apply_dbs_color_mode(img, args)
    imshow_kwargs = {"origin": "upper"}
    if dbs_bg.get("color_mode") == "gray":
        imshow_kwargs.update({"cmap": "gray", "vmin": 0.0, "vmax": 1.0})
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(display_img, **imshow_kwargs)
    if image_corner_mode == "cpp_transpose_fliplr":
        labels = {
            3: "corner 3\nPNG left-top",
            0: "corner 0\nPNG right-top",
            2: "corner 2\nPNG right-bottom",
            1: "corner 1\nPNG left-bottom",
        }
        draw_order = [3, 0, 2, 1]
    else:
        labels = {
            0: "corner 0\nPNG left-top",
            2: "corner 2\nPNG right-top",
            1: "corner 1\nPNG right-bottom",
            3: "corner 3\nPNG left-bottom",
        }
        draw_order = [0, 2, 1, 3]
    for idx in draw_order:
        x, y = corner_pixels[idx]
        ax.scatter([x], [y], s=80, c="yellow", edgecolors="black", zorder=3)
        row = corners[idx]
        ax.text(x, y,
                f"{labels[idx]}\nlon={row['lon']:.6f}\nlat={row['lat']:.6f}",
                color="yellow", fontsize=8, zorder=4)
    ax.set_title(f"DBS original corners, image_corner_mode={image_corner_mode}")
    save_figure(fig, os.path.join(debug_dir, "dbs_original_with_corner_labels.png"), "png")

    fig, ax = plt.subplots(figsize=(8, 10))
    add_dbs_background(ax, dbs_bg)
    ax.set_xlim(dbs_bg["extent"][0], dbs_bg["extent"][1])
    ax.set_ylim(dbs_bg["extent"][2], dbs_bg["extent"][3])
    if dbs_bg.get("mode") == "image":
        ax.set_xlabel("DBS image x (pixel)")
        ax.set_ylabel("DBS image y (pixel)")
        ax.set_aspect("equal", adjustable="box")
    else:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
    ax.set_title(f"DBS {dbs_bg['mode']} preview")
    save_figure(fig, os.path.join(debug_dir, "dbs_warped_preview.png"), "png")

    fig, ax = plt.subplots(figsize=(8, 10))
    add_dbs_background(ax, dbs_bg)
    if dbs_bg.get("mode") == "image" and not detections.empty and "dbs_x" in detections.columns:
        plot_detections(ax, detections, "dbs_x", "dbs_y", style)
    elif not detections.empty:
        plot_detections(ax, detections, "lon", "lat", style)
    ax.set_xlim(dbs_bg["extent"][0], dbs_bg["extent"][1])
    ax.set_ylim(dbs_bg["extent"][2], dbs_bg["extent"][3])
    if dbs_bg.get("mode") == "image":
        ax.set_xlabel("DBS image x (pixel)")
        ax.set_ylabel("DBS image y (pixel)")
        ax.set_aspect("equal", adjustable="box")
    else:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
    ax.set_title("DBS overlay preview with detections")
    save_figure(fig, os.path.join(debug_dir, "dbs_overlay_preview.png"), "png")

    if dbs_bg.get("mode") == "image":
        fig, ax = plt.subplots(figsize=(8, 8))
        add_dbs_background(ax, dbs_bg)
        if not detections.empty and "dbs_x" in detections.columns:
            plot_detections(ax, detections, "dbs_x", "dbs_y", style)
        if not states.empty and "dbs_x" in states.columns:
            sample = states[states["matched_this_frame"] == 1]
            if not sample.empty:
                ax.scatter(sample["dbs_x"], sample["dbs_y"], s=12,
                           c="#1f77b4", alpha=0.8, marker="s", zorder=4)
        ax.set_xlim(0, dbs_bg["width"] - 1)
        ax.set_ylim(dbs_bg["height"] - 1, 0)
        ax.set_aspect("equal", adjustable="box")
        if not args.dbs_show_pixel_axis:
            ax.set_xticks([])
            ax.set_yticks([])
        ax.set_title("DBS image-mode overlay preview")
        save_figure(fig, os.path.join(debug_dir, "dbs_image_mode_overlay_preview.png"), "png")


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


def add_output_legend(ax, style):
    handles = [
        Line2D([0], [0], color="#4c78a8", label="track color = track_id"),
        Line2D([0], [0], marker="s", color="#4c78a8", markerfacecolor="#4c78a8",
               markeredgecolor="black", label="final output point", linestyle="none"),
        Line2D([0], [0], color="#4c78a8", linewidth=style.track_linewidth,
               label="final output track"),
    ]
    ax.legend(handles=handles, loc="best", fontsize=style.legend_fontsize)


def plot_detections(ax, dets, xcol, ycol, style):
    if dets.empty:
        return
    matched = dets[dets["matched"] == 1]
    unmatched = dets[dets["matched"] != 1]
    if not matched.empty:
        ax.scatter(matched[xcol], matched[ycol], s=style.detection_marker_size,
                   c="#bdbdbd", alpha=style.detection_alpha, marker=".", zorder=2)
    if not unmatched.empty:
        ax.scatter(unmatched[xcol], unmatched[ycol],
                   s=style.detection_marker_size * 1.8,
                   facecolors="none", edgecolors="#737373",
                   alpha=min(1.0, style.detection_alpha + 0.35),
                   marker="o", linewidths=max(0.6, style.output_edge_width * 0.7),
                   zorder=2)


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
                       linewidths=max(0.8, style.output_edge_width * 0.8),
                       zorder=4)
        else:
            ax.scatter(normal[xcol], normal[ycol],
                       s=marker_size(state, False, style),
                       marker=marker, color=color, edgecolors=color,
                       alpha=alpha, linewidths=0.5, zorder=4)
    if not output.empty:
        ax.scatter(output[xcol], output[ycol],
                   s=marker_size(state, True, style),
                   marker=marker, color=color, edgecolors="black",
                   alpha=0.95, linewidths=style.output_edge_width, zorder=4)


def plot_track_paths(ax, states, colors, args, style, current_id=None,
                     label_ends=True, is_frame=False):
    xcol, ycol, _, _ = get_plot_coord_columns(args)
    track_count = states["track_id"].nunique()
    do_labels = should_label_tracks(args, track_count, is_frame=is_frame)
    for tid in sorted(states["track_id"].unique()):
        g = states[states["track_id"] == tid].sort_values("result_id")
        color = colors.get(tid, "#4c78a8")
        ax.plot(g[xcol], g[ycol], color=color,
                linewidth=style.track_linewidth, alpha=style.track_alpha,
                zorder=3)
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
                               alpha=0.95 if is_output else 0.85, zorder=4)
                else:
                    ax.scatter([row[xcol]], [row[ycol]],
                               s=marker_size(state, is_output, style) * 1.25,
                               marker=marker, color=color,
                               edgecolors="black" if is_output else color,
                               linewidths=style.output_edge_width if is_output else 0.7,
                               alpha=0.95 if is_output else (0.7 if state == "Coasted" else 0.85),
                               zorder=4)
                label = f"ID={tid}"
                if args.show_speed:
                    label += f" v={row['speed']:.1f} m/s"
                if do_labels:
                    ax.text(row[xcol], row[ycol], label, fontsize=style.label_fontsize,
                            zorder=5)
        elif label_ends and do_labels:
            last = g.iloc[-1]
            ax.text(last[xcol], last[ycol], f"ID={tid}",
                    fontsize=style.label_fontsize, zorder=5)


def title_suffix(args, roi):
    pieces = [f"min_len >= {args.min_len}"]
    roi_text = roi_title(roi, args.coord)
    if roi_text:
        pieces.append(roi_text)
    return "\n".join(pieces)


def plot_output_track_paths(ax, output_states, colors, args, style,
                            current_id=None, label_ends=True, is_frame=False):
    xcol, ycol, _, _ = get_plot_coord_columns(args)
    if output_states.empty:
        return
    track_count = output_states["track_id"].nunique()
    do_labels = should_label_tracks(args, track_count, is_frame=is_frame)
    for tid in sorted(output_states["track_id"].unique()):
        g = output_states[output_states["track_id"] == tid].sort_values("result_id")
        color = colors.get(tid, "#4c78a8")
        ax.plot(g[xcol], g[ycol], color=color,
                linewidth=style.track_linewidth * 1.3,
                alpha=min(1.0, style.track_alpha + 0.15), zorder=3)
        ax.scatter(g[xcol], g[ycol],
                   s=marker_size("Confirmed", True, style),
                   marker=STATE_MARKERS["Confirmed"], color=color,
                   edgecolors="black", linewidths=style.output_edge_width,
                   alpha=0.95, zorder=4)
        if current_id is not None:
            c = g[g["result_id"] == current_id]
            if not c.empty:
                row = c.iloc[-1]
                ax.scatter([row[xcol]], [row[ycol]],
                           s=marker_size("Confirmed", True, style) * 1.4,
                           marker=STATE_MARKERS["Confirmed"], color=color,
                           edgecolors="black", linewidths=style.output_edge_width * 1.2,
                           alpha=1.0, zorder=4)
                if do_labels:
                    ax.text(row[xcol], row[ycol], f"ID={tid}",
                            fontsize=style.label_fontsize, zorder=5)
        elif label_ends and do_labels:
            last = g.iloc[-1]
            ax.text(last[xcol], last[ycol], f"ID={tid}",
                    fontsize=style.label_fontsize, zorder=5)


def output_title_suffix(args, roi):
    roi_text = roi_title(roi, args.coord)
    return f"\n{roi_text}" if roi_text else ""


def plot_output_global_paths(output_states, args, out_dir, colors, roi, axis_limits, style,
                             dbs_bg=None):
    if output_states.empty:
        print("No final output points available for track_outputs_global_paths.")
        return
    xcol, ycol, xlabel, ylabel = get_plot_coord_columns(args)
    fig, ax = create_spatial_figure(args, dbs_bg=dbs_bg, is_frame=False,
                                    fmt=args.format)
    add_dbs_background(ax, dbs_bg)
    plot_output_track_paths(ax, output_states, colors, args, style)
    ax.set_xlabel(xlabel, fontsize=style.axis_label_fontsize)
    ax.set_ylabel(ylabel, fontsize=style.axis_label_fontsize)
    rid_min = int(output_states["result_id"].min())
    rid_max = int(output_states["result_id"].max())
    ax.set_title(
        f"Final output tracks result_id {rid_min}-{rid_max}{output_title_suffix(args, roi)}",
        fontsize=style.title_fontsize,
    )
    apply_axes(ax, args, style, *axis_limits)
    add_output_legend(ax, style)
    save_named_figure(fig, out_dir, "track_outputs_global_paths", args.format)


def output_frame_stats_text(output_states, rid):
    num_outputs = len(output_states[output_states["result_id"] == rid])
    num_tracks = output_states[output_states["result_id"] <= rid]["track_id"].nunique()
    return f"outputs={num_outputs}, output_tracks_so_far={num_tracks}"


def plot_output_snapshot(output_states, frames, args, out_dir, colors, roi, axis_limits,
                         style, result_id=None, stem="track_outputs_snapshot_last",
                         fmt=None, is_frame=False, dbs_bg=None):
    if output_states.empty:
        print(f"No final output points available for {stem}.")
        return
    xcol, ycol, xlabel, ylabel = get_plot_coord_columns(args)
    rid = int(result_id if result_id is not None else frames["result_id"].max())
    hist = output_states[output_states["result_id"] <= rid]

    fig, ax = create_spatial_figure(args, dbs_bg=dbs_bg, is_frame=is_frame,
                                    fmt=fmt or args.format)
    add_dbs_background(ax, dbs_bg)
    plot_output_track_paths(ax, hist, colors, args, style, current_id=rid,
                            label_ends=False, is_frame=is_frame)
    ax.set_xlabel(xlabel, fontsize=style.axis_label_fontsize)
    ax.set_ylabel(ylabel, fontsize=style.axis_label_fontsize)
    title = f"Final output snapshot, result {rid}\n{output_frame_stats_text(output_states, rid)}"
    if not is_frame:
        title += output_title_suffix(args, roi)
    ax.set_title(title, fontsize=style.title_fontsize)
    apply_axes(ax, args, style, *axis_limits)
    if not is_frame:
        add_output_legend(ax, style)
    save_named_figure(fig, out_dir, stem, fmt or args.format)


def plot_global_paths(states, detections, frames, args, out_dir, colors, roi, axis_limits,
                      style, dbs_bg=None):
    xcol, ycol, xlabel, ylabel = get_plot_coord_columns(args)
    fig, ax = create_spatial_figure(args, dbs_bg=dbs_bg, is_frame=False,
                                    fmt=args.format)
    add_dbs_background(ax, dbs_bg)
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
    apply_axes(ax, args, style, *axis_limits)
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
                  is_frame=False, dbs_bg=None):
    xcol, ycol, xlabel, ylabel = get_plot_coord_columns(args)
    rid = int(result_id if result_id is not None else frames["result_id"].max())
    hist = states[states["result_id"] <= rid]
    current_dets = detections[detections["result_id"] == rid]

    fig, ax = create_spatial_figure(args, dbs_bg=dbs_bg, is_frame=is_frame,
                                    fmt=fmt or args.format)
    add_dbs_background(ax, dbs_bg)
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
    
    apply_axes(ax, args, style, *axis_limits)

    # GIF 帧图不显示图例，普通快照图仍然显示图例
    if not is_frame:
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
                        style, fmt=None, frame_dir=None, dbs_bg=None):
    frame_fmt = fmt or args.format
    frame_dir = frame_dir or os.path.join(out_dir, "frames")
    os.makedirs(frame_dir, exist_ok=True)
    result_ids = sorted(frames["result_id"].unique())
    paths = []
    for rid in result_ids:
        stem = f"frame_{int(rid):06d}"
        plot_snapshot(states, detections, frames, args, frame_dir, colors, roi,
                      axis_limits, style, result_id=int(rid), stem=stem,
                      fmt=frame_fmt, is_frame=True, dbs_bg=dbs_bg)
        paths.append(os.path.join(frame_dir, f"{stem}.{frame_fmt}"))
    return paths


def make_output_frame_sequence(output_states, frames, args, out_dir, colors, roi, axis_limits,
                               style, fmt=None, frame_dir=None, dbs_bg=None):
    frame_fmt = fmt or args.format
    frame_dir = frame_dir or os.path.join(out_dir, "output_frames")
    os.makedirs(frame_dir, exist_ok=True)
    result_ids = sorted(frames["result_id"].unique())
    paths = []
    for rid in result_ids:
        stem = f"output_frame_{int(rid):06d}"
        plot_output_snapshot(output_states, frames, args, frame_dir, colors, roi,
                             axis_limits, style, result_id=int(rid), stem=stem,
                             fmt=frame_fmt, is_frame=True, dbs_bg=dbs_bg)
        path = os.path.join(frame_dir, f"{stem}.{frame_fmt}")
        if os.path.exists(path):
            paths.append(path)
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
    if args.no_dbs_bg:
        args.show_dbs_bg = False
    if args.no_dbs_clip_to_image:
        args.dbs_clip_to_image = False
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

        dbs_bg = build_dbs_background(args)
        if dbs_bg is not None and dbs_bg.get("mode") == "image":
            states, detections = add_dbs_pixel_columns(states, detections, dbs_bg)
            if args.dbs_clip_to_image:
                states = clip_to_dbs_image(states, dbs_bg)
                detections = clip_to_dbs_image(detections, dbs_bg)
            if states.empty:
                print("No tracks remain after DBS image pixel clipping.")
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
        style = resolve_style(args)
        frame_style = resolve_style(args)
        colors = get_track_colors(filtered_states["track_id"].unique())
        axis_limits = compute_axis_limits_with_dbs(
            filtered_states, filtered_detections, roi, args.coord, args.padding_ratio,
            dbs_bg, args.axis_from)
        output_states = filtered_states[filtered_states["is_output"] == 1].copy()
        output_colors = get_track_colors(output_states["track_id"].unique())
        if dbs_bg is not None:
            save_dbs_debug_figures(
                args, dbs_bg["original_image"], dbs_bg["corners"], dbs_bg,
                out_dir, filtered_detections, filtered_states, style)

        plot_global_paths(filtered_states, filtered_detections, frames, args, out_dir,
                          colors, roi, axis_limits, style, dbs_bg=dbs_bg)
        plot_snapshot(filtered_states, filtered_detections, frames, args, out_dir,
                      colors, roi, axis_limits, style, dbs_bg=dbs_bg)
        plot_output_global_paths(output_states, args, out_dir, output_colors, roi,
                                 axis_limits, style, dbs_bg=dbs_bg)
        plot_output_snapshot(output_states, frames, args, out_dir, output_colors, roi,
                             axis_limits, style, dbs_bg=dbs_bg)
        plot_state_timeline(filtered_states, args, out_dir, style)
        plot_speed_timeline(filtered_states, args, out_dir, colors, style)
        plot_short_track_summary(states, args, out_dir, style)

        frame_paths = None
        output_frame_paths = None
        if args.make_frames:
            frame_paths = make_frame_sequence(filtered_states, filtered_detections, frames, args, out_dir,
                                             colors, roi, axis_limits, frame_style,
                                             dbs_bg=dbs_bg)
            if not output_states.empty:
                output_frame_paths = make_output_frame_sequence(
                    output_states, frames, args, out_dir, output_colors, roi,
                    axis_limits, frame_style, dbs_bg=dbs_bg)

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
                                                frame_style, fmt="png", frame_dir=gif_png_dir,
                                                dbs_bg=dbs_bg)
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

        if args.make_gif or args.make_output_gif:
            if not output_states.empty:
                output_temp_dir = None
                if args.format == "png" and output_frame_paths is not None:
                    output_png_paths = output_frame_paths
                else:
                    output_temp_dir = os.path.join(out_dir, "_output_gif_png_frames")
                    output_png_paths = make_output_frame_sequence(
                        output_states, frames, args, out_dir, output_colors, roi,
                        axis_limits, frame_style, fmt="png", frame_dir=output_temp_dir,
                        dbs_bg=dbs_bg)
                output_made = make_gif_from_frames(
                    output_png_paths,
                    os.path.join(out_dir, "track_outputs_growth.gif"),
                    args.gif_duration_ms,
                )
                if output_made and not keep_frame_png and output_temp_dir:
                    shutil.rmtree(output_temp_dir, ignore_errors=True)
                elif output_made and not keep_frame_png and args.format == "png" and output_frame_paths is None:
                    for path in output_png_paths:
                        try:
                            os.remove(path)
                        except OSError:
                            pass
            else:
                print("No final output points available for track_outputs_growth.gif.")
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
