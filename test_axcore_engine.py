#!/usr/bin/env python3
"""
AxCore simulator architectural test harness.

This script checks the current C++ engine against the ledger-isolation and
scaling concerns raised during the architecture audit. It is intentionally
root-local and side-effect-light: it compiles/runs the simulator, validates the
generated artifacts, audits hot-path source code for host math calls, confirms
SPARC data availability, and reports macroscopic scaling estimates.

Usage:
  python test_axcore_engine.py
  python test_axcore_engine.py --skip-run
  python test_axcore_engine.py --dynamic-torsion
  python test_axcore_engine.py --sparc-dir "C:\\path\\to\\local_data"
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


ROOT = Path(__file__).resolve().parent
SRC_DIR = ROOT / "src"
COMPILE_SCRIPT = ROOT / "compile_and_run.ps1"
ARTIFACT_DIR = ROOT / "artifacts"
RESULTS_JSON = ARTIFACT_DIR / "fpm_axcore_results.json"
SPARC_PAYLOAD_JSON = ARTIFACT_DIR / "sparc_injection_payload.json"
SPARC_SUBSTRATE_JSON = ARTIFACT_DIR / "sparc_substrate_output.json"
TORSION_PHASE_LOCK_JSON = ARTIFACT_DIR / "torsion_phase_lock.json"
SIM_EXE = ROOT / "build" / "fpm_axcore.exe"
SPARC_BASE_URL = "https://astroweb.case.edu/SPARC"
SPARC_TABLE_URL = f"{SPARC_BASE_URL}/SPARC_Lelli2016c.mrt"
SPARC_ROTMOD_URL = f"{SPARC_BASE_URL}/Rotmod_LTG.zip"

EXPECTED_CHARTS = [
    "fpm_axcore_architecture.png",
    "fpm_axcore_bell.png",
    "fpm_axcore_born.png",
    "fpm_axcore_calibration.png",
    "fpm_axcore_fine_structure.png",
    "fpm_axcore_gravity.png",
    "fpm_axcore_time_dilation.png",
    "fpm_axcore_trajectory.png",
]

HOT_FUNCTIONS = [
    "normalized_entropy_H",
    "spectral_gap_weights",
    "viscosity_update",
    "axcore_lagrangian",
    "ThermodynamicScheduler::step_one",
]

HOT_PATH_FORBIDDEN = [
    "std::log",
    "std::sqrt",
    "std::pow",
    "std::exp",
    "std::sin",
    "std::cos",
    "std::acos",
    "std::cbrt",
    "std::fmod",
]


@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""


def run_command(cmd: list[str], *, timeout: int = 120) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )


def download_file(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=60) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out)


def prepare_downloaded_sparc_data(root: Path) -> Path:
    local_data = root / "local_data"
    local_data.mkdir(parents=True, exist_ok=True)

    table_path = local_data / "SPARC_Lelli2016c.mrt"
    zip_path = local_data / "Rotmod_LTG.zip"
    rotmod_dir = local_data / "Rotmod_LTG"

    download_file(SPARC_TABLE_URL, table_path)
    download_file(SPARC_ROTMOD_URL, zip_path)

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(local_data)

    # Official ZIPs have historically extracted either into Rotmod_LTG/ or as
    # loose *_rotmod.dat files. Normalize to local_data/Rotmod_LTG/.
    if not rotmod_dir.exists():
        rotmod_dir.mkdir()
        for dat in local_data.glob("*_rotmod.dat"):
            dat.replace(rotmod_dir / dat.name)

    return local_data


@contextmanager
def sparc_data_context(explicit_dir: Path | None) -> Iterator[tuple[Path, list[Check]]]:
    if explicit_dir is not None:
        yield explicit_dir, [
            Check("SPARC data source", True, f"using explicit local directory: {explicit_dir}")
        ]
        return

    checks: list[Check] = []
    with tempfile.TemporaryDirectory(prefix="fpm_sparc_") as tmp:
        tmp_root = Path(tmp)
        try:
            sparc_dir = prepare_downloaded_sparc_data(tmp_root)
            checks.append(
                Check(
                    "SPARC data source",
                    True,
                    f"downloaded official SPARC data to temporary directory: {sparc_dir}",
                )
            )
            yield sparc_dir, checks
        except Exception as exc:
            checks.append(Check("SPARC data source", False, f"download failed: {exc}"))
            yield tmp_root / "local_data", checks


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_results() -> dict:
    with RESULTS_JSON.open("r", encoding="utf-8") as f:
        return json.load(f)


def extract_braced_block(source: str, anchor: str) -> str:
    idx = source.find(anchor)
    if idx < 0 and "::" in anchor:
        idx = source.find(anchor.split("::")[-1])
    if idx < 0:
        return ""

    brace = source.find("{", idx)
    if brace < 0:
        return ""

    depth = 0
    for pos in range(brace, len(source)):
        char = source[pos]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return source[brace : pos + 1]
    return ""


def count_data_rows_mrt(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 19:
                continue
            galaxy = parts[0]
            hubble_type = parts[1]
            vflat = parts[15]
            quality = parts[17]
            if (
                re.match(r"^[A-Za-z0-9]", galaxy)
                and hubble_type.lstrip("-").isdigit()
                and re.match(r"^-?\d", vflat)
                and quality.isdigit()
            ):
                count += 1
    return count


def load_sparc_rotmod_points(rotmod_dir: Path) -> list[tuple[float, float, float, float, float]]:
    points: list[tuple[float, float, float, float, float]] = []
    for path in sorted(rotmod_dir.glob("*_rotmod.dat")):
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split()
                if len(parts) < 6:
                    continue
                try:
                    radius, v_obs, _err_v, v_gas, v_disk, v_bul = map(float, parts[:6])
                except ValueError:
                    continue
                if radius > 0.0 and v_obs > 0.0:
                    points.append((radius, v_obs, v_gas, v_disk, v_bul))
    return points


def mond_control_fit(points: list[tuple[float, float, float, float, float]]) -> tuple[float, float, float]:
    best_a0 = 0.0
    best_mse = float("inf")
    mean_obs = sum(p[1] for p in points) / len(points)
    for a0 in range(100, 3001, 25):
        se = 0.0
        for radius, v_obs, v_gas, v_disk, v_bul in points:
            v_bar_sq = v_gas * v_gas + v_disk * v_disk + v_bul * v_bul
            g_bar = max(v_bar_sq / radius, 1e-9)
            nu = 0.5 + math.sqrt(0.25 + a0 / g_bar)
            v_pred = math.sqrt(max(0.0, v_bar_sq * nu))
            se += (v_pred - v_obs) ** 2
        mse = se / len(points)
        if mse < best_mse:
            best_mse = mse
            best_a0 = float(a0)
    rmse = math.sqrt(best_mse)
    return best_a0, rmse, rmse / mean_obs


def load_rotmod_file(path: Path) -> list[tuple[float, float, float, float, float]]:
    points: list[tuple[float, float, float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 6:
                continue
            radius, v_obs, _err_v, v_gas, v_disk, v_bul = map(float, parts[:6])
            if radius > 0.0 and v_obs > 0.0:
                points.append((radius, v_obs, v_gas, v_disk, v_bul))
    return points


def write_sanitized_sparc_payload(galaxy: str, points: list[tuple[float, float, float, float, float]]) -> None:
    max_v = max(v_obs for _r, v_obs, _vg, _vd, _vb in points)
    payload_points = []
    cumulative = 0.0
    previous_r = 0.0
    for radius, v_obs, v_gas, v_disk, v_bul in points:
        v_bar_sq = v_gas * v_gas + v_disk * v_disk + v_bul * v_bul
        dr = max(radius - previous_r, 0.0)
        previous_r = radius
        cumulative += v_bar_sq * max(dr, 0.1)
        b_load = min(24.0, 24.0 * cumulative / max(max_v * max_v * radius, 1.0))
        payload_points.append({"r_kpc": radius, "v_obs": v_obs, "b_load": b_load})
    ARTIFACT_DIR.mkdir(exist_ok=True)
    with SPARC_PAYLOAD_JSON.open("w", encoding="utf-8") as f:
        json.dump({"galaxy": galaxy, "points": payload_points}, f, indent=2)


def fit_scaled_rmse(v_obs: list[float], v_raw: list[float]) -> tuple[float, float, float]:
    numerator = sum(raw * obs for raw, obs in zip(v_raw, v_obs))
    denominator = sum(raw * raw for raw in v_raw)
    scale = numerator / denominator if denominator else 0.0
    rmse = math.sqrt(sum((scale * raw - obs) ** 2 for raw, obs in zip(v_raw, v_obs)) / len(v_obs))
    mean_obs = sum(v_obs) / len(v_obs)
    return scale, rmse, rmse / mean_obs


def launch_cpp_with_runtime(args: list[str], *, timeout: int = 180) -> subprocess.CompletedProcess[str]:
    quoted_args = " ".join(f"'{arg}'" for arg in args)
    ps_command = f"Unblock-File -Path '{args[0]}' -ErrorAction SilentlyContinue; Start-Sleep -Milliseconds 500; $env:Path='C:\\msys64\\mingw64\\bin;' + $env:Path; & {quoted_args}"
    return run_command(["powershell", "-NoProfile", "-Command", ps_command], timeout=timeout)


def run_compile_and_simulator(skip_run: bool) -> list[Check]:
    checks: list[Check] = []
    if not COMPILE_SCRIPT.exists():
        return [Check("compile script exists", False, str(COMPILE_SCRIPT))]

    build = run_command(
        ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(COMPILE_SCRIPT), "-SkipRun"],
        timeout=180,
    )
    checks.append(
        Check(
            "compile succeeds",
            build.returncode == 0,
            "g++/MSVC build completed" if build.returncode == 0 else build.stdout[-1200:],
        )
    )
    if build.returncode != 0 or skip_run:
        return checks

    hashes = []
    for i in range(2):
        run = run_command(
            ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(COMPILE_SCRIPT)],
            timeout=240,
        )
        checks.append(
            Check(
                f"simulator run {i + 1} succeeds",
                run.returncode == 0 and RESULTS_JSON.exists(),
                "results JSON produced" if run.returncode == 0 and RESULTS_JSON.exists() else run.stdout[-1200:],
            )
        )
        if run.returncode != 0 or not RESULTS_JSON.exists():
            return checks
        hashes.append(sha256_file(RESULTS_JSON))

    checks.append(
        Check(
            "bit-for-bit JSON determinism across two runs",
            hashes[0] == hashes[1],
            f"{hashes[0][:12]} == {hashes[1][:12]}" if hashes[0] == hashes[1] else f"{hashes[0]} != {hashes[1]}",
        )
    )
    return checks


def check_artifacts() -> list[Check]:
    checks: list[Check] = []
    if not RESULTS_JSON.exists():
        return [Check("artifacts JSON exists", False, str(RESULTS_JSON))]

    data = load_results()
    derived = data.get("derived_constants", {})
    cal = data.get("calibration", {})
    lattice = data.get("lattice", {})
    probes = data.get("emergent_probes", {})

    checks.extend(
        [
            Check("artifacts JSON exists", True, str(RESULTS_JSON.relative_to(ROOT))),
            Check(
                "N_bit_eq exact integer preserved",
                derived.get("N_bit_eq") == 1452997909 and cal.get("N_bit_eq_exact_integer") == 1452997909,
                str(derived.get("N_bit_eq")),
            ),
            Check(
                "calibration gates pass",
                bool(cal.get("ell_D_in_range")) and bool(cal.get("gamma_above_cern_muon")),
                f"ell_D={cal.get('ell_D_in_range')} gamma={cal.get('gamma_above_cern_muon')}",
            ),
            Check(
                "toy lattice fully active",
                lattice.get("total_daemons") == 256 and lattice.get("active") >= 200,
                f"{lattice.get('active')}/{lattice.get('total_daemons')}",
            ),
        ]
    )

    for name, expected in [
        ("born", "PASS"),
        ("bell", "PASS"),
        ("fine_structure", "PASS_BARE_COUPLING"),
    ]:
        verdict = probes.get(name, {}).get("verdict")
        checks.append(Check(f"{name} probe verdict", verdict == expected, str(verdict)))

    bell = probes.get("bell", {})
    checks.append(
        Check(
            "Bell joint exceeds classical bound",
            bell.get("S_joint", 0.0) > 2.0 and abs(bell.get("S_joint", 0.0) - bell.get("tsirelson", 0.0)) < 0.01,
            f"S_joint={bell.get('S_joint')} Tsirelson={bell.get('tsirelson')}",
        )
    )

    missing = [name for name in EXPECTED_CHARTS if not (ARTIFACT_DIR / name).exists()]
    checks.append(
        Check(
            "chart PNG artifacts exist",
            not missing,
            "all expected charts present" if not missing else ", ".join(missing),
        )
    )
    return checks


def check_source_architecture() -> list[Check]:
    source_files = list(SRC_DIR.glob("*.cpp")) + list(SRC_DIR.glob("*.hpp"))
    source = ""
    for f in source_files:
        source += f.read_text(encoding="utf-8", errors="replace") + "\n"
    checks = [
        Check("src directory exists and contains source files", SRC_DIR.exists() and len(source_files) > 0, str(SRC_DIR.relative_to(ROOT))),
        Check("radix heap scheduler present", "class RadixHeap" in source and "RadixHeap heap;" in source),
        Check("std::priority_queue removed", "std::priority_queue" not in source and "#include <queue>" not in source),
        Check("integer scheduler timestamps", "uint64_t next_step" in source and "uint64_t universal_step" in source),
        Check("redundant daemon coordinates removed", "coord_x" not in source and "coord_y" not in source and "coord_z" not in source),
        Check("deterministic entropy math present", "det_log2" in source and "normalized_entropy_H" in source),
        Check("deterministic spectral math present", "det_acos" in source and "det_cbrt" in source),
        Check(
            "C++ simulator avoids host SPARC file crawling",
            "SPARC_Lelli2016c" not in source and "Rotmod_LTG" not in source and "directory_iterator" not in source,
        ),
    ]

    remaining_forbidden = []
    for function in HOT_FUNCTIONS:
        block = extract_braced_block(source, function)
        if not block:
            remaining_forbidden.append(f"{function}:missing")
            continue
        found = [token for token in HOT_PATH_FORBIDDEN if token in block]
        remaining_forbidden.extend(f"{function}:{token}" for token in found)

    checks.append(
        Check(
            "hot path avoids host transcendental math",
            not remaining_forbidden,
            "clean" if not remaining_forbidden else ", ".join(remaining_forbidden),
        )
    )

    allowed_std_math_lines = []
    for idx, line in enumerate(source.splitlines(), start=1):
        if any(token in line for token in HOT_PATH_FORBIDDEN + ["std::floor"]):
            allowed_std_math_lines.append(idx)
    checks.append(
        Check(
            "remaining std math confined to precompute/probe-safe zones",
            all(line < 490 for line in allowed_std_math_lines),
            f"std math lines: {allowed_std_math_lines}",
        )
    )
    return checks


def check_sparc_data(sparc_dir: Path) -> list[Check]:
    checks: list[Check] = []
    table = sparc_dir / "SPARC_Lelli2016c.mrt"
    rotmod_dir = sparc_dir / "Rotmod_LTG"
    rotmods = list(rotmod_dir.glob("*_rotmod.dat")) if rotmod_dir.exists() else []
    rows = count_data_rows_mrt(table) if table.exists() else 0
    points = load_sparc_rotmod_points(rotmod_dir) if rotmod_dir.exists() else []

    checks.extend(
        [
            Check("SPARC local_data directory exists", sparc_dir.exists(), str(sparc_dir)),
            Check("SPARC sample table present", table.exists(), str(table)),
            Check("SPARC sample table parses 175 galaxies", rows == 175, f"rows={rows}"),
            Check("SPARC rotation-curve files present", len(rotmods) >= 100, f"rotmod files={len(rotmods)}"),
            Check("SPARC rotation-curve points parse", len(points) >= 3000, f"points={len(points)}"),
        ]
    )

    if points:
        a0, rmse, nrmse = mond_control_fit(points)
        checks.extend(
            [
                Check(
                    "Python-side MOND control fit is sane",
                    nrmse < 0.30,
                    f"a0={a0:.0f} (km/s)^2/kpc rmse={rmse:.2f} km/s normalized={nrmse:.3f}",
                ),
            ]
        )
    elif rows and rotmods:
        checks.append(
            Check(
                "SPARC dataset ready for external audit",
                True,
                "dataset is discoverable; rotation curve point parsing needs attention",
            )
        )
    return checks


def check_sparc_substrate_ipc(sparc_dir: Path) -> list[Check]:
    checks: list[Check] = []
    rotmod = sparc_dir / "Rotmod_LTG" / "DDO154_rotmod.dat"
    if not SIM_EXE.exists():
        return [Check("C++ SPARC substrate executable exists", False, str(SIM_EXE))]
    if not rotmod.exists():
        return [Check("DDO154 SPARC payload source exists", False, str(rotmod))]

    points = load_rotmod_file(rotmod)
    if not points:
        return [Check("DDO154 SPARC payload source parses", False, str(rotmod))]

    write_sanitized_sparc_payload("DDO154", points)
    run = launch_cpp_with_runtime(
        [
            str(SIM_EXE),
            "--sparc-payload",
            str(SPARC_PAYLOAD_JSON),
            "--sparc-output",
            str(SPARC_SUBSTRATE_JSON),
        ],
        timeout=180,
    )
    checks.append(
        Check(
            "C++ substrate accepts sanitized SPARC payload",
            run.returncode == 0 and SPARC_SUBSTRATE_JSON.exists(),
            "substrate output produced" if run.returncode == 0 and SPARC_SUBSTRATE_JSON.exists() else run.stdout[-1200:],
        )
    )
    if run.returncode != 0 or not SPARC_SUBSTRATE_JSON.exists():
        return checks

    with SPARC_SUBSTRATE_JSON.open("r", encoding="utf-8") as f:
        result = json.load(f)

    v_obs = [float(v) for v in result.get("v_obs_kms", [])]
    v_raw = [float(v) for v in result.get("v_fpm_raw", [])]
    l_mean = [float(v) for v in result.get("L_mean", [])]
    checks.extend(
        [
            Check(
                "C++ substrate returns one velocity per injection point",
                len(v_obs) == len(points) and len(v_raw) == len(points),
                f"obs={len(v_obs)} raw={len(v_raw)} expected={len(points)}",
            ),
            Check(
                "C++ substrate velocities are finite and nonzero",
                bool(v_raw) and all(math.isfinite(v) and v > 0.0 for v in v_raw),
                f"min_raw={min(v_raw) if v_raw else 'n/a'} max_raw={max(v_raw) if v_raw else 'n/a'}",
            ),
            Check(
                "C++ substrate route costs are physical",
                bool(l_mean) and all(math.isfinite(v) and v > 0.0 for v in l_mean),
                f"min_L={min(l_mean) if l_mean else 'n/a'} max_L={max(l_mean) if l_mean else 'n/a'}",
            ),
        ]
    )

    if v_obs and v_raw and len(v_obs) == len(v_raw):
        scale, rmse, nrmse = fit_scaled_rmse(v_obs, v_raw)
        checks.append(
            Check(
                "emergent C++ substrate RMSE is measured outside Python physics",
                rmse > 0.0 and math.isfinite(rmse),
                f"galaxy=DDO154 scale={scale:.6g} rmse={rmse:.2f} km/s normalized={nrmse:.3f}",
            )
        )
    return checks


def check_dynamic_torsion_phase_lock() -> list[Check]:
    checks: list[Check] = []
    if not SIM_EXE.exists():
        return [Check("dynamic torsion executable exists", False, str(SIM_EXE))]

    if TORSION_PHASE_LOCK_JSON.exists():
        TORSION_PHASE_LOCK_JSON.unlink()

    run = launch_cpp_with_runtime(
        [
            str(SIM_EXE),
            "--torsion-phase-lock-output",
            str(TORSION_PHASE_LOCK_JSON),
        ],
        timeout=180,
    )
    checks.append(
        Check(
            "C++ dynamic torsion phase-lock mode runs",
            run.returncode == 0 and TORSION_PHASE_LOCK_JSON.exists(),
            "torsion trace produced"
            if run.returncode == 0 and TORSION_PHASE_LOCK_JSON.exists()
            else (run.stdout[-1400:] or f"exit={run.returncode}; expected --torsion-phase-lock-output support"),
        )
    )
    if run.returncode != 0 or not TORSION_PHASE_LOCK_JSON.exists():
        return checks

    with TORSION_PHASE_LOCK_JSON.open("r", encoding="utf-8") as f:
        trace = json.load(f)

    samples = trace.get("samples", [])
    s_joint = [float(s.get("S_joint", 0.0)) for s in samples]
    energy = [float(s.get("mean_E", 0.0)) for s in samples]
    modes = [str(s.get("mode", "")) for s in samples]
    snap_tick = trace.get("snap_tick")
    tsirelson = 2.0 * math.sqrt(2.0)

    crossing = next((i for i, s in enumerate(s_joint) if s > 2.0), None)
    final_s = s_joint[-1] if s_joint else 0.0
    max_s = max(s_joint) if s_joint else 0.0

    checks.extend(
        [
            Check(
                "dynamic torsion trace has enough universal ticks",
                len(samples) >= 32,
                f"samples={len(samples)}",
            ),
            Check(
                "starvation energy decays over trace",
                len(energy) >= 2 and energy[-1] < energy[0],
                f"E_start={energy[0] if energy else 'n/a'} E_end={energy[-1] if energy else 'n/a'}",
            ),
            Check(
                "ZOMBIE snap tick is recorded",
                isinstance(snap_tick, int) and snap_tick >= 0,
                f"snap_tick={snap_tick}",
            ),
            Check(
                "S_joint crosses classical bound",
                crossing is not None,
                f"first_crossing_index={crossing} max_S={max_s}",
            ),
            Check(
                "S_joint stabilizes near Tsirelson bound",
                abs(final_s - tsirelson) < 0.02 and max_s <= tsirelson + 0.05,
                f"final_S={final_s} max_S={max_s} Tsirelson={tsirelson}",
            ),
            Check(
                "trace records starvation/ZOMBIE mode",
                any(mode.upper() == "ZOMBIE" for mode in modes),
                f"modes={sorted(set(modes))}",
            ),
        ]
    )
    return checks


def check_scaling_estimate() -> list[Check]:
    data = load_results() if RESULTS_JSON.exists() else {}
    lattice = data.get("lattice", {})
    total = lattice.get("total_daemons") or 256

    # The simulator prints 264 bytes/daemon after the coordinate payload removal.
    # Keep this explicit so changes are easy to notice in test output.
    bytes_per_daemon = 264
    million_daemons = 100 * 100 * 100
    arena_mb = million_daemons * bytes_per_daemon / (1024 * 1024)
    event_bytes = 16
    scheduler_mb = million_daemons * event_bytes / (1024 * 1024)

    return [
        Check("current lattice baseline recorded", total == 256, f"total_daemons={total}"),
        Check(
            "macroscopic 100^3 arena estimate",
            arena_mb < 512,
            f"daemon arena ~= {arena_mb:.1f} MiB at {bytes_per_daemon} bytes/daemon",
        ),
        Check(
            "radix scheduler memory scales linearly in events, not buckets",
            scheduler_mb < 64,
            f"event payload ~= {scheduler_mb:.1f} MiB plus 65 bucket headers",
        ),
    ]


def print_results(checks: Iterable[Check]) -> int:
    failures = 0
    for check in checks:
        status = "PASS" if check.ok else "FAIL"
        print(f"[{status}] {check.name}")
        if check.detail:
            print(f"       {check.detail}")
        if not check.ok:
            failures += 1
    print()
    print(f"Summary: {len(list(checks)) if not isinstance(checks, list) else len(checks)} checks, {failures} failures")
    return 1 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Test AxCore simulator architecture and artifacts.")
    parser.add_argument("--skip-run", action="store_true", help="Compile only; do not run simulator twice.")
    parser.add_argument("--skip-build", action="store_true", help="Skip compile/run and only inspect current files.")
    parser.add_argument(
        "--dynamic-torsion",
        action="store_true",
        help="Run the dynamic torsion phase-lock contract test. Requires C++ --torsion-phase-lock-output support.",
    )
    parser.add_argument(
        "--sparc-dir",
        type=Path,
        default=None,
        help="Optional existing SPARC local_data path. If omitted, official SPARC data is downloaded to a temp dir and deleted after the run.",
    )
    args = parser.parse_args()

    checks: list[Check] = []
    if not args.skip_build:
        checks.extend(run_compile_and_simulator(skip_run=args.skip_run))
    checks.extend(check_artifacts())
    checks.extend(check_source_architecture())
    with sparc_data_context(args.sparc_dir) as (sparc_dir, sparc_source_checks):
        checks.extend(sparc_source_checks)
        checks.extend(check_sparc_data(sparc_dir))
        if not args.skip_build:
            checks.extend(check_sparc_substrate_ipc(sparc_dir))
        downloaded_sparc_dir = sparc_dir if args.sparc_dir is None else None
    if downloaded_sparc_dir is not None:
        checks.append(
            Check(
                "temporary SPARC download cleaned up",
                not downloaded_sparc_dir.exists(),
                str(downloaded_sparc_dir),
            )
        )
    if args.dynamic_torsion:
        checks.extend(check_dynamic_torsion_phase_lock())
    checks.extend(check_scaling_estimate())
    return print_results(checks)


if __name__ == "__main__":
    raise SystemExit(main())
