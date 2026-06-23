#!/usr/bin/env python3
"""
FPM Telemetry Renderer (Diagnostic HUD)
=======================================

One script for turning AxCore JSON telemetry into high-contrast MP4 videos.

Available render modes:
  photon          Causal wave / emergent photon propagation
  torsion         Dynamic torsion ZOMBIE snap
  gravity         Baryonic cluster route-cost / rotation curve probe
  time-dilation   Near/far tick-rate drift probe
  all             Render every available mode

Requires: matplotlib and numpy. FFmpeg is used for MP4 output; the script checks
PATH first, then the winget `Gyan.FFmpeg` package location, then falls back to GIF.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ARTIFACT_DIR = Path("artifacts")
DEFAULT_RESULTS_JSON = ARTIFACT_DIR / "fpm_axcore_results.json"
DEFAULT_PHOTON_JSON = ARTIFACT_DIR / "photon_wave.json"
DEFAULT_TORSION_JSON = ARTIFACT_DIR / "torsion_phase_lock.json"


COLORS = {
    "bg": "#030305",
    "panel": "#070811",
    "grid": "#0a0a14",
    "muted": "#aaaaaa",
    "white": "#ffffff",
    "cyan": "#00ffcc",
    "blue": "#1a3355",
    "green": "#00e676",
    "yellow": "#ffd54f",
    "red": "#ff0033",
    "orange": "#ff8a00",
    "purple": "#9d7cff",
}


def configure_ffmpeg() -> Path | None:
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        mpl.rcParams["animation.ffmpeg_path"] = ffmpeg
        return Path(ffmpeg)

    candidates = []
    local_appdata = os.environ.get("LOCALAPPDATA")
    if local_appdata:
        winget_root = Path(local_appdata) / "Microsoft" / "WinGet" / "Packages"
        candidates.extend(winget_root.glob("Gyan.FFmpeg_*/*/bin/ffmpeg.exe"))
        candidates.extend(winget_root.glob("Gyan.FFmpeg_*/ffmpeg-*/bin/ffmpeg.exe"))

    for candidate in candidates:
        if candidate.exists():
            mpl.rcParams["animation.ffmpeg_path"] = str(candidate)
            return candidate

    return None


def load_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Telemetry payload not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def setup_hud(figsize=(14, 5), dpi=200):
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor(COLORS["bg"])
    ax.set_facecolor(COLORS["bg"])
    for spine in ax.spines.values():
        spine.set_visible(False)
    return fig, ax


def mono(weight="bold") -> dict:
    return {"fontfamily": "monospace", "weight": weight}


def save_animation(fig, ani, output_mp4: Path, fps: int = 20, bitrate: int = 2500) -> bool:
    output_mp4.parent.mkdir(parents=True, exist_ok=True)
    ffmpeg_path = configure_ffmpeg()
    writer = animation.FFMpegWriter(
        fps=fps,
        metadata={"artist": "AxCore FPM"},
        bitrate=bitrate,
    )
    print(f"Encoding AxCore Diagnostic HUD to {output_mp4}...")
    if ffmpeg_path:
        print(f"Using FFmpeg: {ffmpeg_path}")
    try:
        ani.save(output_mp4, writer=writer)
    except FileNotFoundError:
        gif_output = output_mp4.with_suffix(".gif")
        print("WARN: FFmpeg is not installed or not in system PATH. Falling back to GIF.")
        try:
            ani.save(gif_output, writer=animation.PillowWriter(fps=fps))
        except Exception as exc:
            print(f"ERROR: GIF fallback failed: {exc}")
            plt.close(fig)
            return False
        plt.close(fig)
        print(f"[SUCCESS] GIF fallback rendered: {gif_output}")
        return True
    plt.close(fig)
    print(f"[SUCCESS] H.264 video rendered: {output_mp4}")
    return True


def render_photon_wave_mp4(json_path: Path, output_mp4: Path) -> bool:
    data = load_json(json_path)
    ticks = data.get("ticks", [])
    leading_edge = data.get("leading_edge_x", [])
    peak_x = data.get("peak_x", [])
    lattice_x = data.get("lattice", [100, 4, 4])[0]
    if not ticks:
        raise ValueError("Empty photon telemetry.")

    fig, ax = setup_hud(figsize=(14, 5), dpi=200)
    ax.set_xlim(-2, lattice_x + 2)
    ax.set_ylim(-0.3, 1.5)
    ax.set_xticks([])
    ax.set_yticks([])

    for x in range(lattice_x):
        ax.axvline(x, ymin=0.1, ymax=0.9, color=COLORS["grid"], lw=2, zorder=0)

    fig.text(0.05, 0.90, "AXCORE // TOPOLOGICAL CAUSAL WAVE PROBE", color=COLORS["white"], fontsize=14, **mono())
    tick_text = fig.text(0.05, 0.82, "UNIVERSAL TICK : ----", color=COLORS["cyan"], fontsize=12, **mono())
    fig.text(0.05, 0.77, "PHASE VELOCITY : 1.0 CELLS / TICK (c_light)", color=COLORS["muted"], fontsize=10, **mono())
    state_hud = fig.text(0.70, 0.82, "WAKE STATE: ACTIVE", color=COLORS["white"], fontsize=12, **mono())
    fig.text(0.70, 0.77, "PAYLOAD E : 1.000 [MAX]", color=COLORS["muted"], fontsize=10, **mono())

    x_nodes = np.arange(lattice_x)
    noise_scatter = ax.scatter([], [], c=COLORS["blue"], s=20, marker="s", zorder=2)
    zombie_scatter = ax.scatter([], [], c=COLORS["red"], s=30, marker="x", zorder=3)
    photon_core, = ax.plot([], [], "s", color=COLORS["white"], markersize=8, zorder=6)
    glow_lines = []
    for lw, alpha in [(20, 0.05), (15, 0.1), (10, 0.2), (5, 0.5)]:
        glow, = ax.plot([], [], "|", color=COLORS["cyan"], markersize=lw, markeredgewidth=lw / 2, alpha=alpha, zorder=5)
        glow_lines.append(glow)
    scanner = ax.axvline(-10, ymin=0.1, ymax=0.9, color=COLORS["white"], lw=1, alpha=0.5, zorder=4, linestyle=":")

    def init():
        noise_scatter.set_offsets(np.empty((0, 2)))
        zombie_scatter.set_offsets(np.empty((0, 2)))
        photon_core.set_data([], [])
        for glow in glow_lines:
            glow.set_data([], [])
        tick_text.set_text("UNIVERSAL TICK : ----")
        return [noise_scatter, zombie_scatter, photon_core, scanner, tick_text, state_hud] + glow_lines

    def update(frame):
        t = ticks[frame]
        p_x = int(peak_x[frame])
        edge_x = leading_edge[frame]
        noise_x = x_nodes[x_nodes > p_x]
        zombie_x = x_nodes[x_nodes < p_x]
        noise_scatter.set_offsets(
            np.column_stack((noise_x, np.full_like(noise_x, 0.05))) if len(noise_x) else np.empty((0, 2))
        )
        zombie_scatter.set_offsets(
            np.column_stack((zombie_x, np.zeros_like(zombie_x))) if len(zombie_x) else np.empty((0, 2))
        )
        if p_x < lattice_x:
            photon_core.set_data([p_x], [1.0])
            for glow in glow_lines:
                glow.set_data([p_x], [1.0])
        else:
            photon_core.set_data([], [])
            for glow in glow_lines:
                glow.set_data([], [])
        scanner.set_xdata([edge_x, edge_x])
        tick_text.set_text(f"UNIVERSAL TICK : {t:04d}")
        state_hud.set_text("WAKE STATE: ZOMBIE (E=0)" if len(zombie_x) else "WAKE STATE: THERMAL")
        state_hud.set_color(COLORS["red"] if len(zombie_x) else COLORS["muted"])
        return [noise_scatter, zombie_scatter, photon_core, scanner, tick_text, state_hud] + glow_lines

    ani = animation.FuncAnimation(fig, update, frames=len(ticks), init_func=init, blit=True, interval=50)
    return save_animation(fig, ani, output_mp4, fps=20, bitrate=2500)


def render_torsion_snap_mp4(json_path: Path, output_mp4: Path) -> bool:
    data = load_json(json_path)
    samples = data.get("samples", [])
    if not samples:
        raise ValueError("Empty torsion phase-lock telemetry.")

    ticks = np.array([s["tick"] for s in samples], dtype=float)
    energy = np.array([s["mean_E"] for s in samples], dtype=float)
    s_joint = np.array([s["S_joint"] for s in samples], dtype=float)
    modes = [str(s["mode"]).upper() for s in samples]
    snap_tick = data.get("snap_tick", -1)
    tsirelson = 2.0 * np.sqrt(2.0)

    fig, ax = setup_hud(figsize=(14, 6), dpi=200)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])

    fig.text(0.05, 0.92, "AXCORE // DYNAMIC TORSION PHASE-LOCK", color=COLORS["white"], fontsize=14, **mono())
    tick_text = fig.text(0.05, 0.84, "UNIVERSAL TICK : ----", color=COLORS["cyan"], fontsize=12, **mono())
    energy_text = fig.text(0.05, 0.79, "MEAN LEDGER E : ----", color=COLORS["muted"], fontsize=10, **mono())
    mode_text = fig.text(0.74, 0.84, "MODE : FLOW", color=COLORS["green"], fontsize=12, **mono())
    s_text = fig.text(0.74, 0.79, "S_joint : 2.000000", color=COLORS["muted"], fontsize=10, **mono())

    ax.text(0.18, 0.52, "A", color=COLORS["white"], fontsize=18, ha="center", va="center", **mono())
    ax.text(0.82, 0.52, "B", color=COLORS["white"], fontsize=18, ha="center", va="center", **mono())
    link, = ax.plot([0.22, 0.78], [0.52, 0.52], color=COLORS["cyan"], lw=2, alpha=0.25, zorder=1)
    pulse, = ax.plot([], [], color=COLORS["white"], lw=5, alpha=0.7, zorder=2)
    node_a = ax.scatter([0.18], [0.52], s=900, c=[COLORS["green"]], edgecolors=COLORS["white"], linewidths=1.5, zorder=3)
    node_b = ax.scatter([0.82], [0.52], s=900, c=[COLORS["green"]], edgecolors=COLORS["white"], linewidths=1.5, zorder=3)

    meter_ax = fig.add_axes([0.12, 0.12, 0.76, 0.20], facecolor=COLORS["panel"])
    meter_ax.set_xlim(float(ticks[0]), float(ticks[-1]))
    meter_ax.set_ylim(1.85, 2.95)
    meter_ax.tick_params(colors=COLORS["muted"], labelsize=8)
    for spine in meter_ax.spines.values():
        spine.set_color("#1d1f2a")
    meter_ax.axhline(2.0, color=COLORS["yellow"], ls=":", lw=1, label="classical")
    meter_ax.axhline(tsirelson, color=COLORS["cyan"], ls="--", lw=1, label="Tsirelson")
    if snap_tick >= 0:
        meter_ax.axvline(snap_tick, color=COLORS["red"], lw=1, alpha=0.8, label="ZOMBIE snap")
    meter_ax.set_title("CHSH phase-lock trace", color=COLORS["white"], fontsize=9)
    meter_line, = meter_ax.plot([], [], color=COLORS["green"], lw=2)
    cursor = meter_ax.axvline(float(ticks[0]), color=COLORS["white"], lw=1, alpha=0.65)
    meter_ax.legend(loc="lower right", fontsize=7)

    def color_for_mode(mode: str) -> str:
        if mode == "ZOMBIE":
            return COLORS["red"]
        if mode == "FATIGUE":
            return COLORS["yellow"]
        return COLORS["green"]

    def update(frame):
        mode = modes[frame]
        color = color_for_mode(mode)
        alpha = max(0.12, min(0.95, energy[frame] / max(energy[0], 1e-12)))
        fig.patch.set_facecolor((0.03 * alpha, 0.03 * alpha, 0.05 * alpha))
        ax.set_facecolor((0.02 * alpha, 0.02 * alpha, 0.04 * alpha))
        tick_text.set_text(f"UNIVERSAL TICK : {int(ticks[frame]):04d}")
        energy_text.set_text(f"MEAN LEDGER E : {energy[frame]:.6e}")
        mode_text.set_text(f"MODE : {mode}")
        mode_text.set_color(color)
        s_text.set_text(f"S_joint : {s_joint[frame]:.9f}")
        s_text.set_color(COLORS["cyan"] if s_joint[frame] > 2.0 else COLORS["muted"])
        node_a.set_color(color)
        node_b.set_color(color)
        link.set_alpha(0.25 if mode != "ZOMBIE" else 0.85)
        link.set_color(COLORS["cyan"] if mode != "ZOMBIE" else COLORS["red"])

        x0 = 0.22 + 0.56 * ((frame % 20) / 19.0)
        pulse.set_data([x0 - 0.05, x0 + 0.05], [0.52, 0.52])
        pulse.set_color(COLORS["white"] if mode != "ZOMBIE" else COLORS["red"])
        pulse.set_alpha(0.25 if mode == "FLOW" else 0.8)

        meter_line.set_data(ticks[: frame + 1], s_joint[: frame + 1])
        meter_line.set_color(color)
        cursor.set_xdata([ticks[frame], ticks[frame]])
        return [link, pulse, node_a, node_b, meter_line, cursor, tick_text, energy_text, mode_text, s_text]

    ani = animation.FuncAnimation(fig, update, frames=len(samples), blit=False, interval=50)
    return save_animation(fig, ani, output_mp4, fps=20, bitrate=2800)


def render_gravity_mp4(json_path: Path, output_mp4: Path) -> bool:
    data = load_json(json_path)
    gravity = data.get("emergent_probes", {}).get("gravity", {})
    radii = np.array(gravity.get("r", []), dtype=float)
    velocity = np.array(gravity.get("v_emergent_kms", []), dtype=float)
    route_cost = np.array(gravity.get("L_mean", []), dtype=float)
    if len(radii) == 0:
        raise ValueError("Gravity probe telemetry is missing from results JSON.")

    fig, ax = setup_hud(figsize=(14, 6), dpi=200)
    ax.set_aspect("equal")
    max_r = float(max(radii) + 1.0)
    ax.set_xlim(-max_r, max_r)
    ax.set_ylim(-max_r, max_r)
    ax.set_xticks([])
    ax.set_yticks([])

    fig.text(0.05, 0.92, "AXCORE // BARYONIC CLUSTER ROUTE-COST WELL", color=COLORS["white"], fontsize=14, **mono())
    radius_text = fig.text(0.05, 0.84, "SHELL RADIUS : ----", color=COLORS["cyan"], fontsize=12, **mono())
    speed_text = fig.text(0.05, 0.79, "EMERGENT v(r) : ----", color=COLORS["muted"], fontsize=10, **mono())
    cost_text = fig.text(0.68, 0.84, "ROUTE COST L : ----", color=COLORS["orange"], fontsize=12, **mono())

    ax.scatter([0], [0], s=700, c=COLORS["white"], marker="*", zorder=5)
    theta = np.linspace(0, 2 * np.pi, 240)
    shell_line, = ax.plot([], [], color=COLORS["cyan"], lw=2, zorder=4)
    orbit_dot, = ax.plot([], [], "o", color=COLORS["green"], markersize=8, zorder=6)
    arrow = ax.arrow(0, 0, 0, 0, color=COLORS["green"], width=0.02, head_width=0.18, alpha=0.0)

    for r in radii:
        ax.plot(r * np.cos(theta), r * np.sin(theta), color=COLORS["grid"], lw=1, alpha=0.8, zorder=1)

    curve_ax = fig.add_axes([0.62, 0.14, 0.30, 0.24], facecolor=COLORS["panel"])
    curve_ax.tick_params(colors=COLORS["muted"], labelsize=8)
    for spine in curve_ax.spines.values():
        spine.set_color("#1d1f2a")
    curve_ax.set_title("emergent rotation curve", color=COLORS["white"], fontsize=9)
    curve_ax.set_xlim(float(min(radii)) - 0.2, float(max(radii)) + 0.2)
    curve_ax.set_ylim(0, float(max(velocity)) * 1.15)
    curve_line, = curve_ax.plot([], [], color=COLORS["cyan"], lw=2)
    curve_dot, = curve_ax.plot([], [], "o", color=COLORS["white"], markersize=5)

    def update(frame):
        nonlocal arrow
        r = radii[frame]
        v = velocity[frame]
        l_val = route_cost[frame] if frame < len(route_cost) else 0.0
        phase = frame * 0.45
        shell_line.set_data(r * np.cos(theta), r * np.sin(theta))
        orbit_dot.set_data([r * np.cos(phase)], [r * np.sin(phase)])
        arrow.remove()
        arrow = ax.arrow(
            r * np.cos(phase),
            r * np.sin(phase),
            -0.75 * np.sin(phase),
            0.75 * np.cos(phase),
            color=COLORS["green"],
            width=0.02,
            head_width=0.18,
            alpha=0.8,
            zorder=6,
        )
        radius_text.set_text(f"SHELL RADIUS : {r:.2f}")
        speed_text.set_text(f"EMERGENT v(r) : {v:.2f} km/s")
        cost_text.set_text(f"ROUTE COST L : {l_val:.6f}")
        curve_line.set_data(radii[: frame + 1], velocity[: frame + 1])
        curve_dot.set_data([r], [v])
        return [shell_line, orbit_dot, arrow, curve_line, curve_dot, radius_text, speed_text, cost_text]

    frames = max(48, len(radii) * 12)
    frame_map = [min(len(radii) - 1, int(i * len(radii) / frames)) for i in range(frames)]

    def mapped_update(frame):
        return update(frame_map[frame])

    ani = animation.FuncAnimation(fig, mapped_update, frames=frames, blit=False, interval=50)
    return save_animation(fig, ani, output_mp4, fps=20, bitrate=2600)


def render_time_dilation_mp4(json_path: Path, output_mp4: Path) -> bool:
    data = load_json(json_path)
    td = data.get("emergent_probes", {}).get("time_dilation", {})
    ticks_near = int(td.get("ticks_near", 0))
    ticks_far = int(td.get("ticks_far", 0))
    gamma_near = float(td.get("gamma_near", 0.0))
    gamma_far = float(td.get("gamma_far", 0.0))
    ratio = float(td.get("tick_ratio", 0.0))
    if ticks_near <= 0 or ticks_far <= 0:
        raise ValueError("Time-dilation telemetry is missing from results JSON.")

    fig, ax = setup_hud(figsize=(14, 5), dpi=200)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.text(0.05, 0.90, "AXCORE // EMERGENT TIME-DILATION GRADIENT", color=COLORS["white"], fontsize=14, **mono())
    tick_text = fig.text(0.05, 0.82, "UNIVERSAL TICK : ----", color=COLORS["cyan"], fontsize=12, **mono())
    ratio_text = fig.text(0.66, 0.82, f"TICK RATIO : {ratio:.6f}", color=COLORS["muted"], fontsize=11, **mono())

    near_center = (0.28, 0.48)
    far_center = (0.72, 0.48)
    near_clock = plt.Circle(near_center, 0.16, color="#1b1010", ec=COLORS["red"], lw=2, zorder=2)
    far_clock = plt.Circle(far_center, 0.16, color="#101b14", ec=COLORS["green"], lw=2, zorder=2)
    ax.add_patch(near_clock)
    ax.add_patch(far_clock)
    ax.text(near_center[0], 0.20, "NEAR HIGH-LOAD CLUSTER", ha="center", color=COLORS["red"], fontsize=10, **mono())
    ax.text(far_center[0], 0.20, "FAR LOW-COST VACUUM", ha="center", color=COLORS["green"], fontsize=10, **mono())
    near_hand, = ax.plot([], [], color=COLORS["red"], lw=4, zorder=4)
    far_hand, = ax.plot([], [], color=COLORS["green"], lw=4, zorder=4)
    near_count = ax.text(near_center[0], 0.72, "0 ticks", ha="center", color=COLORS["red"], fontsize=12, **mono())
    far_count = ax.text(far_center[0], 0.72, "0 ticks", ha="center", color=COLORS["green"], fontsize=12, **mono())
    ax.plot([0.5, 0.5], [0.18, 0.78], color=COLORS["grid"], lw=2)

    frames = 120

    def hand(center, angle, radius=0.12):
        return [center[0], center[0] + radius * np.cos(angle)], [center[1], center[1] + radius * np.sin(angle)]

    def update(frame):
        frac = frame / (frames - 1)
        near_ticks_now = int(round(frac * ticks_near))
        far_ticks_now = int(round(frac * ticks_far))
        tick_text.set_text(f"UNIVERSAL TICK : {frame:04d}")
        near_hand.set_data(*hand(near_center, -np.pi / 2 + 2 * np.pi * near_ticks_now / max(ticks_near, 1)))
        far_hand.set_data(*hand(far_center, -np.pi / 2 + 2 * np.pi * far_ticks_now / max(ticks_far, 1)))
        near_count.set_text(f"{near_ticks_now} ticks  gamma={gamma_near:.3f}")
        far_count.set_text(f"{far_ticks_now} ticks  gamma={gamma_far:.3f}")
        ratio_text.set_text(f"TICK RATIO : {ratio:.6f}")
        return [near_hand, far_hand, near_count, far_count, tick_text, ratio_text]

    ani = animation.FuncAnimation(fig, update, frames=frames, blit=True, interval=50)
    return save_animation(fig, ani, output_mp4, fps=20, bitrate=2300)


def render_modes(mode: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = {
        "photon": (render_photon_wave_mp4, DEFAULT_PHOTON_JSON, output_dir / "causal_wave_render.mp4"),
        "torsion": (render_torsion_snap_mp4, DEFAULT_TORSION_JSON, output_dir / "torsion_snap_render.mp4"),
        "gravity": (render_gravity_mp4, DEFAULT_RESULTS_JSON, output_dir / "gravity_well_render.mp4"),
        "time-dilation": (render_time_dilation_mp4, DEFAULT_RESULTS_JSON, output_dir / "time_dilation_render.mp4"),
    }

    selected = list(jobs) if mode == "all" else [mode]
    for name in selected:
        renderer, json_path, output_path = jobs[name]
        try:
            renderer(json_path, output_path)
        except Exception as exc:
            print(f"[SKIP] {name}: {exc}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render AxCore telemetry videos.")
    parser.add_argument(
        "--mode",
        choices=["photon", "torsion", "gravity", "time-dilation", "all"],
        default="photon",
        help="Which video to render.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Optional telemetry JSON override for a single mode.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional MP4 output override for a single mode.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ARTIFACT_DIR,
        help="Directory used by --mode all or by default output names.",
    )
    args = parser.parse_args()

    if args.mode == "all" and (args.input or args.output):
        parser.error("--input/--output can only be used with a single render mode")

    if args.input or args.output:
        defaults = {
            "photon": (render_photon_wave_mp4, DEFAULT_PHOTON_JSON, args.output_dir / "causal_wave_render.mp4"),
            "torsion": (render_torsion_snap_mp4, DEFAULT_TORSION_JSON, args.output_dir / "torsion_snap_render.mp4"),
            "gravity": (render_gravity_mp4, DEFAULT_RESULTS_JSON, args.output_dir / "gravity_well_render.mp4"),
            "time-dilation": (render_time_dilation_mp4, DEFAULT_RESULTS_JSON, args.output_dir / "time_dilation_render.mp4"),
        }
        renderer, default_input, default_output = defaults[args.mode]
        renderer(args.input or default_input, args.output or default_output)
        return 0

    render_modes(args.mode, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
