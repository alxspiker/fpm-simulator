#!/usr/bin/env python3
"""
FPM Telemetry Renderer (Diagnostic HUD)
=======================================
Reads C++ FPM substrate telemetry and encodes an H.264 .mp4 video.
Visualizes the discrete memory array, highlighting the causal wavefront,
thermal vacuum noise, and the E=0 ZOMBIE wake.

Requires: matplotlib, numpy, ffmpeg (must be installed on host OS)
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path

def render_photon_wave_mp4(json_path: Path, output_mp4: Path):
    if not json_path.exists():
        print(f"ERROR: Telemetry payload {json_path} not found.")
        return

    with open(json_path, "r") as f:
        data = json.load(f)

    ticks = data.get("ticks", [])
    leading_edge = data.get("leading_edge_x", [])
    peak_x = data.get("peak_x", [])
    lattice_x = data.get("lattice", [100, 4, 4])[0]

    if not ticks:
        print("ERROR: Empty thermodynamic telemetry.")
        return

    # --- CINEMATIC HUD SETUP ---
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(14, 5), dpi=200)
    
    # Deep void background
    fig.patch.set_facecolor('#030305')
    ax.set_facecolor('#030305')
    
    ax.set_xlim(-2, lattice_x + 2)
    ax.set_ylim(-0.3, 1.5)
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Strip standard spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Render discrete memory register grid
    for x in range(lattice_x):
        ax.axvline(x, ymin=0.1, ymax=0.9, color='#0a0a14', lw=2, zorder=0)

    # --- HUD TEXT OVERLAYS ---
    font_mono = {'fontfamily': 'monospace', 'weight': 'bold'}
    
    title_text = fig.text(0.05, 0.90, "AXCORE // TOPOLOGICAL CAUSAL WAVE PROBE", color="#ffffff", fontsize=14, **font_mono)
    tick_text = fig.text(0.05, 0.82, "UNIVERSAL TICK : ----", color="#00ffcc", fontsize=12, **font_mono)
    speed_text = fig.text(0.05, 0.77, "PHASE VELOCITY : 1.0 CELLS / TICK (c_light)", color="#aaaaaa", fontsize=10, **font_mono)
    
    state_hud = fig.text(0.70, 0.82, "WAKE STATE: ACTIVE", color="#ffffff", fontsize=12, **font_mono)
    energy_hud = fig.text(0.70, 0.77, "PAYLOAD E : 1.000 [MAX]", color="#aaaaaa", fontsize=10, **font_mono)

    # --- DISCRETE DATA VISUALIZATION ELEMENTS ---
    x_nodes = np.arange(lattice_x)
    
    # 1. Background Vacuum Noise (E = 0.05)
    noise_scatter = ax.scatter([], [], c='#1a3355', s=20, marker='s', zorder=2)
    
    # 2. ZOMBIE Wake (E = 0.0)
    zombie_scatter = ax.scatter([], [], c='#ff0033', s=30, marker='x', zorder=3)
    
    # 3. Photon Soliton (E = E_max)
    photon_core, = ax.plot([], [], 's', color="#ffffff", markersize=8, zorder=6)
    
    # Fake a neon glow using stacked lines with decreasing alpha
    glow_lines = []
    for lw, alpha in [(20, 0.05), (15, 0.1), (10, 0.2), (5, 0.5)]:
        glow, = ax.plot([], [], '|', color="#00ffcc", markersize=lw, markeredgewidth=lw/2, alpha=alpha, zorder=5)
        glow_lines.append(glow)

    # Leading Edge Scanner
    scanner = ax.axvline(-10, ymin=0.1, ymax=0.9, color='#ffffff', lw=1, alpha=0.5, zorder=4, linestyle=':')

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

        # Determine states
        noise_x = x_nodes[x_nodes > p_x]
        zombie_x = x_nodes[x_nodes < p_x]
        
        # Update offsets (Requires [N, 2] shape for scatter)
        if len(noise_x) > 0:
            noise_scatter.set_offsets(np.column_stack((noise_x, np.full_like(noise_x, 0.05))))
        else:
            noise_scatter.set_offsets(np.empty((0, 2)))
            
        if len(zombie_x) > 0:
            zombie_scatter.set_offsets(np.column_stack((zombie_x, np.zeros_like(zombie_x))))
        else:
            zombie_scatter.set_offsets(np.empty((0, 2)))

        # Update Soliton Core & Glow
        if p_x < lattice_x:
            photon_core.set_data([p_x], [1.0])
            for glow in glow_lines:
                glow.set_data([p_x], [1.0])
        else:
            photon_core.set_data([], [])
            for glow in glow_lines:
                glow.set_data([], [])

        # Move Scanner
        scanner.set_xdata([edge_x])

        # Update HUD Text
        tick_text.set_text(f"UNIVERSAL TICK : {t:04d}")
        
        if len(zombie_x) > 0:
            state_hud.set_text("WAKE STATE: ZOMBIE (E=0)")
            state_hud.set_color("#ff0033")
        else:
            state_hud.set_text("WAKE STATE: THERMAL")
            state_hud.set_color("#aaaaaa")

        return [noise_scatter, zombie_scatter, photon_core, scanner, tick_text, state_hud] + glow_lines

    # Compile Animation
    ani = animation.FuncAnimation(
        fig, update, frames=len(ticks),
        init_func=init, blit=True, interval=50
    )

    writer = animation.FFMpegWriter(fps=20, metadata=dict(artist='AxCore FPM'), bitrate=2500)
    
    print(f"Encoding AxCore Diagnostic HUD to {output_mp4}...")
    try:
        ani.save(output_mp4, writer=writer)
        print(f"[SUCCESS] H.264 High-Fidelity Video rendered: {output_mp4}")
    except FileNotFoundError:
        print("ERROR: FFmpeg is not installed or not in system PATH. Cannot render MP4.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, default="artifacts/photon_wave.json")
    parser.add_argument("--output", type=str, default="artifacts/causal_wave_render.mp4")
    args = parser.parse_args()
    
    render_photon_wave_mp4(Path(args.input), Path(args.output))
