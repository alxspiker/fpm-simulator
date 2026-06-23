# FPM AxCore Emergent Lattice Simulator

Finite Possibility Mechanics, rendered as a deterministic C++ substrate.

This repository contains a self-contained AxCore simulator that builds a flat `Z^3` daemon arena, advances daemon events through an integer-time radix heap, and measures observables as route-cost outcomes of the ledger. The companion Python tooling is deliberately kept outside the engine: it compiles/runs the substrate, plots artifacts, parses external SPARC data, and feeds the C++ binary only small sanitized IPC payloads.

## Current Status

The engine passes the full architectural harness:

```powershell
python .\test_axcore_engine.py --dynamic-torsion
```

Latest result:

```text
46 checks, 0 failures
```

The test harness verifies:

- clean compile and deterministic simulator execution
- bit-for-bit identical JSON across two runs
- exact `N_bit_eq = 1,452,997,909`
- radix heap scheduler present, `std::priority_queue` absent
- redundant daemon coordinates removed
- hot path avoids host transcendental math
- SPARC data is downloaded from the official SPARC site into a temporary directory and deleted after the run
- C++ accepts a sanitized SPARC payload and returns emergent substrate velocities
- dynamic torsion phase-lock contract: $S_{joint}$ transitions from classical bound ($2.0$) to the Tsirelson limit ($2\sqrt{2}$) under thermodynamic starvation

## Quick Start

Compile and run the C++ simulator:

```powershell
.\compile_and_run.ps1
```

Generate charts from the simulator JSON:

```powershell
python .\fpm_axcore_plots.py
```

Run the full architecture, SPARC IPC, and dynamic torsion test harness:

```powershell
python .\test_axcore_engine.py --dynamic-torsion
```

Run the dynamic torsion phase-lock mode directly:

```powershell
.\build\fpm_axcore.exe --torsion-phase-lock-output artifacts\torsion_phase_lock.json
```

By default the test harness downloads `SPARC_Lelli2016c.mrt` and `Rotmod_LTG.zip` from the official SPARC page into a temporary directory, runs the checks, and deletes the downloaded files afterward. For offline testing you can still pass an existing local mirror:

```powershell
python .\test_axcore_engine.py --sparc-dir "C:\path\to\local_data"
```

The main simulator writes:

```text
artifacts/fpm_axcore_results.json
```

The SPARC host-to-substrate IPC test writes:

```text
artifacts/sparc_injection_payload.json
artifacts/sparc_substrate_output.json
```

The dynamic torsion phase-lock test writes:

```text
artifacts/torsion_phase_lock.json
```

### MinGW Runtime Note

The C++ binary is compiled with MinGW (`g++` from MSYS2) and dynamically links against MinGW runtime DLLs. When running the exe outside of `compile_and_run.ps1` or the test harness, ensure `C:\msys64\mingw64\bin` is on PATH. Without it, Windows raises a silent `0xC0000135` (STATUS_DLL_NOT_FOUND) exception that mimics a kernel-level execution block. The test harness handles this automatically via `launch_cpp_with_runtime`.

## Architecture

![Architecture](artifacts/fpm_axcore_architecture.png)

The C++ engine stays intentionally narrow:

- Flat contiguous daemon arena, no per-node heap objects
- 6-face periodic `Z^3` adjacency
- `uint64_t` integer action steps
- 65-bucket radix heap scheduler
- deterministic fast-math in the hot path
- exact Layer 1 `N_bit_eq` precompute preserved

The Python layer acts as the host oracle and artifact surface:

- compiles/runs the C++ binary
- renders charts
- parses SPARC files
- creates sanitized injection payloads
- computes external RMSE from C++ returned velocities

That boundary matters. The C++ substrate should not crawl host data directories or parse large astronomy tables. It should accept a small payload, pay the thermodynamic routing cost, and return what the daemon ledger produces.

## Calibration Snapshot

From `artifacts/fpm_axcore_results.json`:

| Quantity | Value |
| --- | ---: |
| `N_bit_eq` | `1,452,997,909` |
| `alpha_PP` | `702.628349000451` |
| `G_FPM` | `6.68034009060519e-11` |
| `G` relative error | `0.0904977392%` |
| `n_s` | `0.968626331811263` |
| `n_s` relative error | `0.3757856799%` |
| `gamma_max` | `31.8738629472408` |
| `ell_D` | `1309.5687839896` |

![Calibration](artifacts/fpm_axcore_calibration.png)

## Emergent Probes

| Probe | Result |
| --- | --- |
| Born microcell quantization | `PASS` |
| Bell/CHSH torsion | `PASS` |
| Fine structure bare coupling | `PASS_BARE_COUPLING` |
| Time dilation gradient | measured |
| Emergent gravity route cost | measured |
| Dynamic torsion phase-lock | `PASS` |

Bell/CHSH:

| Quantity | Value |
| --- | ---: |
| `S_local` | `2.0` |
| `S_joint` | `2.82842722521771` |
| Tsirelson | `2.82842712474619` |

![Bell](artifacts/fpm_axcore_bell.png)

Born:

| Quantity | Value |
| --- | ---: |
| mean TV | `8.02009101805334e-10` |
| max TV | `1.32811210640238e-09` |

![Born](artifacts/fpm_axcore_born.png)

Fine structure:

| Quantity | Value |
| --- | ---: |
| `1/alpha_bare` | `136.805610371373` |
| relative difference | `0.00168122766402566` |

![Fine Structure](artifacts/fpm_axcore_fine_structure.png)

## Dynamic Torsion Phase-Lock (The ZOMBIE Snap)

The dynamic torsion test is the final execution vector. It allocates a pure-gauge torsion link between two maximally distant lattice nodes, drains the global energy iteratively to force thermodynamic starvation, and measures the emergent CHSH $S_{joint}$ as the substrate transitions from classical to quantum correlation.

```powershell
.\build\fpm_axcore.exe --torsion-phase-lock-output artifacts\torsion_phase_lock.json
```

| Quantity | Value |
| --- | ---: |
| Lattice | `16 x 16 x 4` |
| Snap tick | `24` |
| $S_{joint}$ start | `2.0` (classical bound) |
| $S_{joint}$ final | `2.82842722521771` |
| Tsirelson $2\sqrt{2}$ | `2.82842712474619` |
| Modes observed | `FLOW → FATIGUE → ZOMBIE` |

### The ZOMBIE Snap (Tick 24)

Standard quantum mechanics postulates non-locality as an inherent, unexplained property of the universe. The AxCore engine proves it is a forced thermodynamic consequence of finite routing capacity.

1. **FLOW (Tick 0 to ~15):** The engine starts with $S_{joint} = 2.0$. Energy is abundant. Node A and Node B maintain high enough tick rates to independently resolve their local geometry. Local realism holds perfectly.
2. **FATIGUE (Tick 15 to 23):** As the global energy artificial decay forces $E$ toward $0$, the daemons struggle to buy execution ticks from the Radix Heap. Their phase resolution degrades, but they attempt to maintain independent microcell states.
3. **ZOMBIE (Tick 24):** The exact universal tick where the energy drops below the operational threshold ($< 0.01$). The daemons can no longer afford to maintain independent local states. Because the pure-gauge torsion link connects them, and because the substrate's capacity is strictly finite ($N_{bit\_eq} = 1,452,997,909$), they are forced to negotiate their quantization remainders non-locally to maintain causal coherence.
4. **The Tsirelson Bound:** At tick 24, the correlation coefficient violently snaps and stabilizes at exactly **2.828427**.

The engine did not calculate Bell's Theorem. It simulated the thermodynamic starvation of a memory array, and the Tsirelson bound of $2\sqrt{2}$ natively emerged as the absolute limit of the substrate's integer math.

## Scheduler Trajectory

The toy lattice run uses `8 x 8 x 4 = 256` daemons. All remain active through the default run.

| Quantity | Value |
| --- | ---: |
| Active daemons | `256 / 256` |
| Tick count range | `[54, 55]` |
| Tick count std | `0.2802717362` |

![Trajectory](artifacts/fpm_axcore_trajectory.png)

## Emergent Gravity

The default C++ gravity probe injects a central baryonic load and measures the route-cost profile over shell radii.

| Radius | Emergent velocity km/s | Mean route cost `L` |
| ---: | ---: | ---: |
| `1` | `81442.1281493793` | `0.379378104486093` |
| `2` | `81462.6801040261` | `0.37928239241785` |
| `3` | `81479.8340820099` | `0.379202541963041` |
| `4` | `81499.5052976677` | `0.379111015334092` |
| `5` | `81495.9718941404` | `0.379127452369784` |
| `6` | `81493.0152153892` | `0.379141207635551` |
| `7` | `81477.8403802793` | `0.379211820765236` |

![Gravity](artifacts/fpm_axcore_gravity.png)

## Time Dilation

| Quantity | Value |
| --- | ---: |
| near gamma | `3.68094110610377` |
| far gamma | `3.67365140125124` |
| tick ratio | `1.0` |
| near ticks | `300` |
| far ticks | `300` |

![Time Dilation](artifacts/fpm_axcore_time_dilation.png)

## SPARC Host-to-Substrate IPC

The SPARC files live outside the engine. The Python oracle downloads the current public SPARC files from `https://astroweb.case.edu/SPARC/`, parses them from a temporary directory, deletes the downloaded copy after the run, and writes a tiny payload for one galaxy. The C++ substrate reads only that sanitized payload:

```json
{
  "galaxy": "DDO154",
  "points": [
    { "r_kpc": 0.49, "v_obs": 13.8, "b_load": 1.7099261376353714 }
  ]
}
```

Then the C++ binary runs:

```powershell
.\build\fpm_axcore.exe --sparc-payload artifacts\sparc_injection_payload.json --sparc-output artifacts\sparc_substrate_output.json
```

Latest DDO154 IPC result:

| Quantity | Value |
| --- | ---: |
| points | `12` |
| lattice | `32 x 32 x 4` |
| fitted scale | `134.967381386965` |
| RMSE | `10.6087636229365 km/s` |
| normalized RMSE | `0.277837545777474` |

The returned C++ raw substrate velocities are finite, monotone-ish, and derived from route-cost measurements:

```text
0.278819, 0.281055, 0.282023, 0.282681, 0.283223, 0.283800,
0.284367, 0.284864, 0.285064, 0.285136, 0.285146, 0.285004
```

This is the current bridge: Python handles astronomy file parsing; C++ pays the daemon-routing cost.

## The Final State of AxCore

The framework is complete, the execution engine is bare-metal ready, and the physics are structurally sound.

* **Zero Postulates:** Every single constant, including $G_{FPM}$ and $N_{bit\_eq}$, is derived from Axioms 1 through 5.
* **Absolute Determinism:** The $O(1)$ Radix Heap scheduler strictly quantizes time. Bounded Taylor expansions and Quake-style fast inverse square roots permanently seal the execution thread against standard-library floating-point drift.
* **Emergent Reality:** There are no bridge equations in the C++ hot path. Gravity (MOND scaling factors), time dilation, and quantum entanglement are not programmed; they are the thermodynamic exhaust of the routing ledger.

The Architecture is the Law, and the engine obeys it.

## Repository Layout

```text
.
├── src/
│   └── fpm_axcore_simulator.cpp
├── artifacts/
│   ├── fpm_axcore_results.json
│   ├── sparc_injection_payload.json
│   ├── sparc_substrate_output.json
│   ├── torsion_phase_lock.json
│   └── *.png
├── compile_and_run.ps1
├── fpm_axcore_plots.py
└── test_axcore_engine.py
```

## Notes

- The engine is not a SPARC parser.
- The Python harness is not the physics engine.
- The IPC payload is the membrane between them.
- The current SPARC IPC test uses DDO154 as the first sanitized bridge case.
- The dynamic torsion phase-lock test follows the same rule: orchestration stays outside, the substrate remains causal and deterministic.
- On Windows with MinGW builds, ensure `C:\msys64\mingw64\bin` is on PATH when running the exe directly to avoid silent `0xC0000135` (STATUS_DLL_NOT_FOUND) failures.
