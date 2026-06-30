#!/usr/bin/env python3
"""Generate a human-readable GMTI simulation test report.

The script intentionally treats unavailable algorithm products as "not run" or
"not parsed" instead of inventing metrics. It summarizes the data that is
present under an output directory produced by scripts/run_stage2_one_click_test.sh.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
import struct
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="ignore")


def load_json(path: Path) -> Dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def count_csv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        return max(0, sum(1 for _ in f) - 1)


def read_csv_dicts(path: Path, limit: Optional[int] = None) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        for i, row in enumerate(csv.DictReader(f)):
            rows.append(row)
            if limit is not None and i + 1 >= limit:
                break
    return rows


def parse_stage2_report(path: Path) -> Dict[str, str]:
    text = read_text(path)
    result: Dict[str, str] = {}
    for line in text.splitlines():
        m = re.match(r"-\s+([^:]+):\s*(.+)", line.strip())
        if m:
            result[m.group(1).strip()] = m.group(2).strip()
    return result


def file_size(path: Path) -> int:
    return path.stat().st_size if path.exists() else 0


def human_size(n: int) -> str:
    x = float(n)
    for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
        if x < 1024.0 or unit == "TiB":
            return f"{x:.2f} {unit}"
        x /= 1024.0
    return f"{n} B"


def summarize_target_beams(path: Path) -> Dict[str, object]:
    rows = read_csv_dicts(path)
    visible = [r for r in rows if int(float(r.get("visible_pulse_count", "0") or 0)) > 0]
    injected = [r for r in rows if int(float(r.get("injected_sample_count", "0") or 0)) > 0]
    amps = [float(r.get("mean_target_amplitude", "nan")) for r in rows if r.get("mean_target_amplitude")]
    radial = [float(r.get("mean_radial_velocity_mps", "nan")) for r in rows if r.get("mean_radial_velocity_mps")]
    samples = [float(r.get("mean_range_sample", "nan")) for r in rows if r.get("mean_range_sample")]
    return {
        "beam_rows": len(rows),
        "visible_beams": len(visible),
        "injected_beams": len(injected),
        "visible_beam_ids": [r.get("beam_id", "") for r in visible[:12]],
        "mean_amplitude": statistics.fmean(amps) if amps else None,
        "mean_radial_velocity": statistics.fmean(radial) if radial else None,
        "mean_range_sample": statistics.fmean(samples) if samples else None,
    }


def summarize_algorithm_outputs(result_dir: Path, algorithm_log: Path) -> Dict[str, object]:
    gmti_bins = sorted(result_dir.glob("GMTI[0-9][0-9].bin")) if result_dir.exists() else []
    track_bins = sorted(result_dir.glob("GMTI[0-9][0-9]_track.bin")) if result_dir.exists() else []
    pngs = sorted(result_dir.glob("*.png")) if result_dir.exists() else []
    txts = sorted(result_dir.glob("*.txt")) if result_dir.exists() else []
    log = read_text(algorithm_log)
    normal_exit = bool(log) and ("GMTI 扫描完成" in log or "--- GMTI" in log)
    errors = [line for line in log.splitlines() if "[ERR]" in line or "error" in line.lower()]
    warnings = [line for line in log.splitlines() if "[WARN]" in line or "warn" in line.lower()]
    timing_ms = None
    for m in re.finditer(r"main_total:\s*([0-9.]+)\s*ms", log):
        timing_ms = float(m.group(1))
    return {
        "log_exists": algorithm_log.exists(),
        "normal_exit": normal_exit,
        "error_count": len(errors),
        "warning_count": len(warnings),
        "timing_ms": timing_ms,
        "gmti_bin_count": len(gmti_bins),
        "track_bin_count": len(track_bins),
        "image_count": len(pngs),
        "text_product_count": len(txts),
        "sample_files": [p.name for p in (gmti_bins + track_bins + pngs + txts)[:10]],
    }


def estimate_mdv(system: Dict) -> Dict[str, float]:
    fc = float(system.get("fc_ghz", 16.0)) * 1e9
    prf = float(system.get("prf_hz", 1300.0))
    pulse_num = float(system.get("pulse_num", 130.0))
    c = 299792458.0
    lam = c / fc
    dwell = pulse_num / prf
    doppler_bin_hz = 1.0 / dwell
    velocity_bin = lam * doppler_bin_hz / 2.0
    blind_first = lam * prf / 2.0
    return {
        "lambda_m": lam,
        "dwell_s": dwell,
        "doppler_bin_hz": doppler_bin_hz,
        "velocity_resolution_mps": velocity_bin,
        "first_blind_speed_mps": blind_first,
    }


def sample_ddc_power(data_file: Path, ddc_len: int, max_packets: int = 8) -> Dict[str, float]:
    if not data_file.exists():
        return {}
    packet_bytes = 256 + ddc_len * 16
    powers: List[float] = []
    max_abs = 0.0
    with data_file.open("rb") as f:
        for _ in range(max_packets):
            pkt = f.read(packet_bytes)
            if len(pkt) != packet_bytes:
                break
            payload = memoryview(pkt)[256:]
            for n in range(0, len(payload), 16):
                vals = struct.unpack_from("<ffff", payload, n)
                c1p = vals[0] * vals[0] + vals[1] * vals[1]
                c2p = vals[2] * vals[2] + vals[3] * vals[3]
                powers.append(0.5 * (c1p + c2p))
                max_abs = max(max_abs, *(abs(v) for v in vals))
    if not powers:
        return {}
    return {
        "sampled_complex_power": statistics.fmean(powers),
        "sampled_rms": math.sqrt(statistics.fmean(powers)),
        "sampled_max_component_abs": max_abs,
        "sampled_packet_count": min(max_packets, file_size(data_file) // packet_bytes),
    }


def write_report(args: argparse.Namespace) -> Path:
    out_dir = Path(args.output_dir)
    reports = out_dir / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    report_path = reports / "full_simulation_test_report.md"

    cfg = load_json(Path(args.stage2_config))
    copied_cfg = load_json(out_dir / "config" / "stage2_config.json")
    if copied_cfg:
        cfg = copied_cfg
    system = cfg.get("system", {})
    scene = cfg.get("scene", {})
    simulation = cfg.get("simulation", {})
    targets = cfg.get("targets", {})

    sim_report = parse_stage2_report(reports / "stage2_simulation_report.md")
    data_file = out_dir / "data" / "stage2_statistical_newprotocol.bin"
    scene_truth = out_dir / "truth" / "scene_truth.csv"
    target_pulse = out_dir / "truth" / "target_truth_pulse.csv"
    target_beam = out_dir / "truth" / "target_truth_beam_summary.csv"
    alg_log = out_dir / "logs" / "algorithm_run.log"
    alg_result_dir = out_dir / "algorithm_result"

    target_summary = summarize_target_beams(target_beam)
    alg = summarize_algorithm_outputs(alg_result_dir, alg_log)
    mdv = estimate_mdv(system)
    ddc = sample_ddc_power(data_file, int(system.get("ddc_len", 11820)))

    scene_counts = {
        "scene_truth": count_csv_rows(scene_truth),
        "area": count_csv_rows(out_dir / "truth" / "area_clutter_scatterers.csv"),
        "strong": count_csv_rows(out_dir / "truth" / "strong_scatterers.csv"),
        "line": count_csv_rows(out_dir / "truth" / "line_scatterers.csv"),
        "target_pulse": count_csv_rows(target_pulse),
        "target_beam": count_csv_rows(target_beam),
    }

    lines: List[str] = []
    lines.append("# GMTI 阶段二一键仿真测试报告\n")
    lines.append("## 1. 测试结论\n")
    lines.append(f"- 输出目录：`{out_dir}`")
    lines.append(f"- 回波数据：`{data_file}`，大小 {human_size(file_size(data_file))}")
    lines.append(f"- 仿真 PRT 数：{sim_report.get('packets_written', '未知')}")
    lines.append(f"- 散射点总数：{scene_counts['scene_truth']}，其中面杂波 {scene_counts['area']}，强散射点 {scene_counts['strong']}，线状散射点 {scene_counts['line']}")
    lines.append(f"- 合作目标 truth：pulse {scene_counts['target_pulse']} 行，beam summary {scene_counts['target_beam']} 行")
    lines.append(f"- NaN/Inf：has_nan={sim_report.get('has_nan', '未知')}，has_inf={sim_report.get('has_inf', '未知')}")
    if alg["log_exists"]:
        lines.append(f"- 算法运行：log 已生成，normal_exit={alg['normal_exit']}，错误 {alg['error_count']}，警告 {alg['warning_count']}")
    else:
        lines.append("- 算法运行：未执行或未发现 `logs/algorithm_run.log`，成像/对消/检测/航迹指标仅报告数据侧准备情况。")

    lines.append("\n## 2. 仿真参数设置\n")
    lines.append("| 参数 | 值 |")
    lines.append("|---|---:|")
    for key in [
        "fc_ghz", "bandwidth_mhz", "fs_mhz", "pulse_width_us", "prf_hz",
        "sample_delay_us", "sample_window_us", "ddc_len", "fft_len",
        "pc_crop_start", "pc_crop_len", "scan_min_deg", "scan_step_deg",
        "beam_count", "beam_width_deg", "pulse_num", "platform_height_m",
        "platform_speed_mps", "d_chan_m",
    ]:
        lines.append(f"| {key} | {system.get(key, '未配置')} |")
    lines.append(f"\n- random_seed：`{simulation.get('random_seed', '未配置')}`")
    lines.append(f"- scene_mode：`{args.scene_mode}`")
    lines.append(f"- target_enabled：`{args.target_enabled}`")

    lines.append("\n## 3. 回波数据生成\n")
    lines.append(f"- 数据格式：新协议 PRT，256 字节头 + `{system.get('ddc_len', 11820)} * 16` 字节双通道 float32 IQ。")
    lines.append(f"- `scatterer_echoes`：{sim_report.get('scatterer_echoes', '未知')}")
    lines.append(f"- `scatterer_samples`：{sim_report.get('scatterer_samples', '未知')}")
    lines.append(f"- `max_abs_component`：{sim_report.get('max_abs_component', '未知')}")
    lines.append(f"- `mean_noise_power_per_complex_sample`：{sim_report.get('mean_noise_power_per_complex_sample', '未知')}")
    if ddc:
        lines.append(f"- 抽样 DDC RMS：{ddc['sampled_rms']:.6g}，抽样平均复功率：{ddc['sampled_complex_power']:.6g}，抽样最大分量幅度：{ddc['sampled_max_component_abs']:.6g}")

    lines.append("\n## 4. 合作目标注入\n")
    lines.append(f"- target_config_path：`{targets.get('target_config_path', args.target_config)}`")
    lines.append(f"- target_snr_db：{targets.get('target_snr_db', '未配置')}")
    lines.append(f"- amplitude_mode：`{targets.get('amplitude_mode', '未配置')}`")
    lines.append(f"- target_pulses_injected：{sim_report.get('target_pulses_injected', '未知')}")
    lines.append(f"- target_samples_injected：{sim_report.get('target_samples_injected', '未知')}")
    lines.append(f"- 可见 beam 数：{target_summary.get('visible_beams')}，已注入 beam 数：{target_summary.get('injected_beams')}，示例 beam：{', '.join(target_summary.get('visible_beam_ids', [])) or '无'}")
    if target_summary.get("mean_radial_velocity") is not None:
        lines.append(f"- truth 平均径向速度：{target_summary['mean_radial_velocity']:.6g} m/s")
    if target_summary.get("mean_range_sample") is not None:
        lines.append(f"- truth 平均 DDC range_sample：{target_summary['mean_range_sample']:.6g}")

    lines.append("\n## 5. 成像效果\n")
    if alg["image_count"] or alg["text_product_count"]:
        lines.append(f"- 算法图像/文本产品数量：PNG {alg['image_count']}，TXT {alg['text_product_count']}")
        lines.append(f"- 示例文件：{', '.join(alg['sample_files']) or '无'}")
    else:
        lines.append("- 未发现算法成像产品。当前报告只能确认 DDC 回波和 truth 已生成；请用生成的 XML 跑 `GMTI_core` 后再刷新报告。")

    lines.append("\n## 6. 对消效果\n")
    if alg["log_exists"]:
        log = read_text(alg_log)
        keywords = [line for line in log.splitlines() if any(k in line for k in ["对消", "cancel", "CSI", "clutter", "杂波"])]
        if keywords:
            lines.append("- 算法日志中找到对消/杂波相关记录：")
            for line in keywords[:12]:
                lines.append(f"  - `{line[:180]}`")
        else:
            lines.append("- 算法日志存在，但未解析到明确对消指标。建议后续在算法中输出对消前后功率或 SCR。")
    else:
        lines.append("- 未运行算法，无法评价对消效果。数据侧已生成统计杂波、强散射点、线状散射体和噪声。")

    lines.append("\n## 7. 检测效果\n")
    if alg["gmti_bin_count"]:
        lines.append(f"- 检测结果文件数量：{alg['gmti_bin_count']}")
        lines.append("- 当前脚本尚未解析 `GMTIxx.bin` 二进制字段，检测率/虚警率需下一步接入结果解析。")
    else:
        lines.append("- 未发现 `GMTIxx.bin` 检测结果，检测率、虚警数和漏检数暂不可评估。")

    lines.append("\n## 8. 最小可检测速度评估\n")
    lines.append("- 当前没有自动 SNR/速度梯度实验结果，因此这里给出体制尺度估计，不能替代真实检测门限。")
    lines.append(f"- 波长：{mdv['lambda_m']:.6g} m")
    lines.append(f"- 单波位驻留时间：{mdv['dwell_s']:.6g} s")
    lines.append(f"- 多普勒 bin 宽度：{mdv['doppler_bin_hz']:.6g} Hz")
    lines.append(f"- 对应速度分辨尺度：{mdv['velocity_resolution_mps']:.6g} m/s")
    lines.append(f"- 第一盲速尺度：{mdv['first_blind_speed_mps']:.6g} m/s")
    lines.append("- 建议后续自动生成 0, 1, 2, 3, 5, 10 m/s 目标速度组，并以 90% 检测率定义最小可检测速度。")

    lines.append("\n## 9. 定位精度\n")
    if alg["gmti_bin_count"]:
        lines.append("- 检测结果存在，但当前未解析检测点坐标，尚未与 target truth 匹配。")
    else:
        lines.append("- 未发现检测结果，无法计算距离误差、方位误差和位置误差。")
    lines.append("- truth 文件已经提供 `range_m`、`range_sample_float`、`target_azimuth_deg`、`theta_true_deg`，可作为后续匹配基准。")

    lines.append("\n## 10. 航迹关联\n")
    if alg["track_bin_count"]:
        lines.append(f"- 航迹文件数量：{alg['track_bin_count']}")
        lines.append("- 当前脚本尚未解析 28 字节 track 包字段，TrackID 稳定性需要下一步接入解析。")
    else:
        lines.append("- 未发现 `GMTIxx_track.bin`，confirmed_track_count、id_switch_count、track_break_count 暂不可评估。")

    lines.append("\n## 11. 输出文件索引\n")
    for rel in [
        "config/stage2_config.json",
        "config/temp_config_stage2_newsystem.xml",
        "data/stage2_statistical_newprotocol.bin",
        "truth/scene_truth.csv",
        "truth/target_truth_pulse.csv",
        "truth/target_truth_beam_summary.csv",
        "reports/stage2_simulation_report.md",
        "logs/simulate_stage2.log",
        "logs/algorithm_run.log",
    ]:
        p = out_dir / rel
        lines.append(f"- `{p}`：{'存在' if p.exists() else '未生成'}")

    lines.append("\n## 12. 下一步建议\n")
    lines.append("1. 用单点场景跑 `GMTI_core`，标定 LFM 符号、rect 定义和脉压峰固定偏移。")
    lines.append("2. 增加 `GMTIxx.bin` 和 `GMTIxx_track.bin` 解析器，自动计算检测率、虚警率、定位误差和航迹稳定性。")
    lines.append("3. 增加速度/SNR 梯度批量生成，给出真实最小可检测速度曲线。")
    lines.append("4. 对 30000 点场景增加 beam/range 索引和并行化，降低生成时间。")

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--stage2-config", default="stage2_config.json")
    parser.add_argument("--target-config", default="targets.json")
    parser.add_argument("--scene-mode", default="full")
    parser.add_argument("--target-enabled", default="true")
    args = parser.parse_args()
    path = write_report(args)
    print(path)


if __name__ == "__main__":
    main()

