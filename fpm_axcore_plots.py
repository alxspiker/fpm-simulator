#!/usr/bin/env python3
"""
FPM AxCore -- Visualization Companion
Reads fpm_axcore_results.json (produced by the C++ simulator) and generates
all PNG charts. Run after the C++ simulator completes.

Usage:  python fpm_axcore_plots.py [path_to_json]
"""

import json, os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 110

def load_json(path="fpm_axcore_results.json"):
    with open(path) as f:
        return json.load(f)

def plot_all(data, out_dir="fpm_emergent_charts"):
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    d = data.get("derived_constants", {})
    cal = data.get("calibration", {})
    probes = data.get("emergent_probes", {})
    lat = data.get("lattice", {})
    traj = data.get("trajectory", {})

    # 1. Calibration summary
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    ax.axis("off")
    lines = [
        f"N_bit_eq (exact integer) = {d.get('N_bit_eq',0):,}",
        f"G_FPM = {d.get('G_FPM',0):.4e}  (CODATA: 6.6743e-11, rel err {cal.get('G_FPM_rel_err_pct',0):.4f}%)",
        f"n_s = {d.get('n_s',0):.6f}  (Planck: 0.965, rel err {cal.get('n_s_rel_err_pct',0):.4f}%)",
        f"ell_D in range: {cal.get('ell_D_in_range','?')}",
        f"gamma_max > CERN muon: {cal.get('gamma_above_cern_muon','?')}",
        f"alpha_PP = {d.get('alpha_PP',0):.9f}",
        f"Rounding leak: {cal.get('rounding_leak_relative',0):.2e}  (ZERO by exact integer)",
    ]
    for i, line in enumerate(lines):
        ax.text(0.05, 0.90 - 0.12*i, line, fontsize=11, family="monospace",
                transform=ax.transAxes, color="#1a2a4a")
    ax.set_title("FPM AxCore Calibration Audit", fontsize=14, color="#1a2a4a")
    p = os.path.join(out_dir, "fpm_axcore_calibration.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    paths["calibration"] = p

    # 2. Scheduler trajectory
    if traj.get("t"):
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
        t = traj["t"]
        axes[0,0].plot(t, traj["total_E"], color="#1a4a6a")
        axes[0,0].set_title("Total energy (closed ledger)")
        axes[0,0].set_xlabel("universal tick")
        if traj.get("mean_L"):
            axes[0,1].plot(t[:len(traj["mean_L"])], traj["mean_L"], color="#a83232", label="L")
            axes[0,1].axhline(d.get("L_max",0), ls=":", color="grey", lw=0.5, label=f"L_max")
            axes[0,1].axhline(d.get("L_rest",0), ls=":", color="grey", lw=0.5, label=f"L_rest")
            axes[0,1].set_title("Mean route cost L"); axes[0,1].legend(fontsize=7)
        if traj.get("mean_Omega"):
            axes[1,0].plot(t[:len(traj["mean_Omega"])], traj["mean_Omega"], color="#2a5a8a")
            axes[1,0].set_title("Mean viscosity Omega")
        if traj.get("active_fraction"):
            axes[1,1].plot(t[:len(traj["active_fraction"])], traj["active_fraction"], color="#2d7a4a")
            axes[1,1].set_title("Active daemon fraction")
            axes[1,1].set_ylim(0, 1.05)
        p = os.path.join(out_dir, "fpm_axcore_trajectory.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        paths["trajectory"] = p

    # 3. Emergent gravity
    grav = probes.get("gravity", {})
    if grav.get("r"):
        fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
        ax.plot(grav["r"], grav["v_emergent_kms"], "o-", color="#1a4a6a",
                lw=2, markersize=6, label="Emergent v(r) from lattice")
        ax.set_xlabel("Radius (lattice units)")
        ax.set_ylabel("Effective velocity v(r) [km/s]")
        ax.set_title("Emergent gravity: v(r) from viscosity gradient")
        ax.legend(); ax.grid(True, alpha=0.3)
        p = os.path.join(out_dir, "fpm_axcore_gravity.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        paths["gravity"] = p

    # 4. Emergent time dilation
    td = probes.get("time_dilation", {})
    if td:
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        labels = ["Near cluster\n(high L)", "Far from cluster\n(low L)"]
        ticks = [td.get("ticks_near",0), td.get("ticks_far",0)]
        gammas = [td.get("gamma_near",0), td.get("gamma_far",0)]
        bars = ax.bar([0,1], ticks, color=["#a83232","#2d7a4a"], alpha=0.85, width=0.5)
        ax.set_xticks([0,1]); ax.set_xticklabels(labels)
        ax.set_ylabel("Local tick count")
        ax.set_title(f"Emergent time dilation  (tick ratio = {td.get('tick_ratio',0):.3f})")
        for bar, g in zip(bars, gammas):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height(),
                    f"gamma={g:.2f}", ha="center", va="bottom", fontsize=9)
        p = os.path.join(out_dir, "fpm_axcore_time_dilation.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        paths["time_dilation"] = p

    # 5. Emergent Born
    born = probes.get("born", {})
    if born:
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        ax.text(0.1, 0.70, f"Mean TV distance = {born.get('mean_tv',0):.2e}", fontsize=12, transform=ax.transAxes)
        ax.text(0.1, 0.55, f"Max TV distance = {born.get('max_tv',0):.2e}", fontsize=12, transform=ax.transAxes)
        ax.text(0.1, 0.40, f"Verdict: {born.get('verdict','?')}", fontsize=14, fontweight="bold",
                transform=ax.transAxes, color="#2d7a4a" if born.get("verdict")=="PASS" else "#a83232")
        ax.text(0.1, 0.22, "Mechanism: microcell quantization\n(finite substrate, no probability formula)",
                fontsize=10, transform=ax.transAxes, style="italic", color="#555555")
        ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
        ax.set_title("Emergent Born rule: |psi|^2 from lattice counting")
        p = os.path.join(out_dir, "fpm_axcore_born.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        paths["born"] = p

    # 6. Emergent Bell/CHSH
    bell = probes.get("bell", {})
    if bell:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        labels = ["Classical\nbound", "Local\ntorsion", "Joint\ntorsion", "Tsirelson"]
        vals = [2.0, bell.get("S_local",0), bell.get("S_joint",0), bell.get("tsirelson",0)]
        colors = ["#777777","#a83232","#2d7a4a","#1a4a6a"]
        axes[0].bar(labels, vals, color=colors)
        axes[0].axhline(2.0, ls=":", color="#555555", lw=1)
        axes[0].axhline(bell.get("tsirelson",0), ls="--", color="#1a4a6a", lw=1)
        axes[0].set_ylabel("CHSH S"); axes[0].set_ylim(0, 3.1)
        axes[0].set_title(f"S_joint = {bell.get('S_joint',0):.6f}  ({bell.get('verdict','')})")

        # Summary text
        axes[1].axis("off")
        axes[1].text(0.1, 0.80, f"S_local = {bell.get('S_local',0):.6f}", fontsize=12, transform=axes[1].transAxes)
        axes[1].text(0.1, 0.65, f"S_joint = {bell.get('S_joint',0):.6f}", fontsize=12, transform=axes[1].transAxes)
        axes[1].text(0.1, 0.50, f"Tsirelson = {bell.get('tsirelson',0):.6f}", fontsize=12, transform=axes[1].transAxes)
        axes[1].text(0.1, 0.30, f"Verdict: {bell.get('verdict','?')}", fontsize=14, fontweight="bold",
                     transform=axes[1].transAxes, color="#2d7a4a" if bell.get("verdict")=="PASS" else "#a83232")
        axes[1].set_title("Emergent Bell/CHSH from torsion flux")
        p = os.path.join(out_dir, "fpm_axcore_bell.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        paths["bell"] = p

    # 7. Fine structure
    fs = probes.get("fine_structure", {})
    if fs:
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        labels = ["FPM bare\n(dx_univ)", "CODATA macro\n(screened)"]
        vals = [fs.get("one_over_alpha_bare",0), 137.035999084]
        ax.bar(labels, vals, color=["#1a4a6a","#5f6f86"])
        ax.set_ylabel("1/alpha"); ax.set_ylim(135.5, 137.5)
        ax.set_title(f"Bare coupling = {fs.get('one_over_alpha_bare',0):.4f}  "
                     f"({100*fs.get('rel_diff_from_macro',0):.3f}% from macro)")
        p = os.path.join(out_dir, "fpm_axcore_fine_structure.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        paths["fine_structure"] = p

    # 8. Lattice topology
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.axis("off")
    info = [
        f"Lattice: {lat.get('size',[0,0,0])}",
        f"Total daemons: {lat.get('total_daemons',0)}",
        f"Active: {lat.get('active',0)}  Halted: {lat.get('halted',0)}",
        f"Adjacency: {lat.get('adjacency','?')}",
        f"Events processed: {lat.get('scheduler_events_processed',0)}",
        "",
        "Memory: flat std::vector<Daemon>",
        "Neighbor lookup: flat-index modulo (no tuple hashing)",
        "Scheduler: std::priority_queue (thermodynamic gating)",
        "N_bit_eq: exact integer (zero rounding leak)",
    ]
    for i, line in enumerate(info):
        ax.text(0.05, 0.92 - 0.09*i, line, fontsize=11, family="monospace",
                transform=ax.transAxes, color="#1a2a4a")
    ax.set_title("AxCore Z^3 Lattice Architecture", fontsize=14, color="#1a2a4a")
    p = os.path.join(out_dir, "fpm_axcore_architecture.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    paths["architecture"] = p

    return paths

if __name__ == "__main__":
    json_path = sys.argv[1] if len(sys.argv) > 1 else "fpm_axcore_results.json"
    data = load_json(json_path)
    paths = plot_all(data)
    print("Charts generated:")
    for k, p in paths.items():
        print(f"  {k:20s}: {p}")
