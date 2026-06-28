#include "modes.hpp"

// =============================================================================
// MAIN
// =============================================================================

int main(int argc, char** argv) {
    if (argc > 1) {
        std::string payload_path;
        std::string output_path = "artifacts/sparc_substrate_output.json";
        std::string torsion_output_path;
        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "--sparc-payload" && i + 1 < argc) {
                payload_path = argv[++i];
            } else if (arg == "--sparc-output" && i + 1 < argc) {
                output_path = argv[++i];
            } else if (arg == "--torsion-phase-lock-output" && i + 1 < argc) {
                torsion_output_path = argv[++i];
            } else if (arg == "--photon-propagation-output" && i + 1 < argc) {
                std::string photon_path = argv[++i];
                std::filesystem::create_directories(std::filesystem::path(photon_path).parent_path());
                return run_photon_propagation_mode(photon_path);
            } else {
                std::cerr << "ERROR: unknown argument: " << arg << "\n";
                return 1;
            }
        }
        if (!torsion_output_path.empty()) {
            std::filesystem::path parent = std::filesystem::path(torsion_output_path).parent_path();
            if (!parent.empty()) std::filesystem::create_directories(parent);
            return run_dynamic_torsion_mode(torsion_output_path);
        }
        if (payload_path.empty()) {
            std::cerr << "ERROR: --sparc-payload requires a sanitized payload path\n";
            return 1;
        }
        std::filesystem::create_directories(std::filesystem::path(output_path).parent_path());
        return run_sparc_payload_mode(payload_path, output_path);
    }

    std::cout << "========================================================================\n";
    std::cout << "FINITE POSSIBILITY MECHANICS v7.0-axcore-cpp -- EMERGENT LATTICE SIMULATOR\n";
    std::cout << "========================================================================\n\n";
    std::cout << "Architecture: flat memory arena + bitwise Z^3 adjacency + radix heap scheduler\n";
    std::cout << "  (No bridge equations. The universe computes itself.)\n\n";

    // Layer 0
    std::cout << "Layer 0: Loading the five axioms...\n";
    Axioms ax;
    std::cout << "  dim_space=" << ax.dim_space << " dim_causal=" << ax.dim_causal
              << " n_directed=" << ax.n_directed() << " n_trace=" << ax.n_trace() << "\n\n";

    // Layer 1
    std::cout << "Layer 1: Deriving all 22 constants (exact integer N_bit_eq in C++)...\n";
    DerivedConstants d = derive_all(ax);
    std::cout << std::setprecision(10);
    std::cout << "  alpha      = " << d.alpha << "     (paper: 0.2)\n";
    std::cout << "  beta       = " << d.beta << "     (paper: 1.8)\n";
    std::cout << "  c0         = " << d.c0 << "     (paper: 0.05)\n";
    std::cout << "  L_max      = " << d.L_max << "   (paper: 3.285)\n";
    std::cout << "  L_rest     = " << d.L_rest << " (paper: 0.1030625)\n";
    std::cout << "  gamma_max  = " << d.gamma_max << "   (paper: 31.8739)\n";
    std::cout << "  alpha_PP   = " << d.alpha_PP << " (paper: 702.628349)\n";
    std::cout << "  N_bit_eq   = " << d.N_bit_eq << " (EXACT INTEGER, zero rounding leak)\n";
    std::cout << "  G_FPM      = " << d.G_FPM << "   (CODATA: 6.6743e-11)\n";
    std::cout << "  dt_univ    = " << d.dt_univ << " s\n";
    std::cout << "  dx_univ    = " << d.dx_univ*1e15 << " fm\n\n";

    // Calibration
    std::cout << "Layer 7: Calibration check...\n";
    CalibrationResult cal = calibrate(d);
    std::cout << "  G_FPM vs CODATA: " << cal.rel_err_G*100 << "%\n";
    std::cout << "  n_s vs Planck: " << cal.rel_err_ns*100 << "%\n";
    std::cout << "  ell_D in range: " << (cal.ell_D_in_range ? "YES" : "NO") << "\n";
    std::cout << "  gamma > CERN muon: " << (cal.gamma_above_cern ? "YES" : "NO") << "\n";
    std::cout << "  N_bit_eq rounding leak: "
              << std::abs((double)cal.N_bit_eq_exact - cal.N_bit_eq_continuous) / cal.N_bit_eq_continuous * 100 << "%\n\n";

    // Layer 2: Z^3 Lattice
    std::cout << "Layer 2: Building Z^3 lattice (flat contiguous arena)...\n";
    Z3Lattice lattice(8, 8, 4, d, 17);
    std::cout << "  Lattice: " << lattice.sx << "x" << lattice.sy << "x" << lattice.sz
              << " = " << lattice.size() << " daemons\n";
    std::cout << "  Memory: " << sizeof(Daemon) * lattice.size() << " bytes ("
              << sizeof(Daemon) << " bytes/daemon)\n";
    std::cout << "  Torsion links: " << lattice.torsion_links.size() << "\n";
    std::cout << "  6-face Z^3 adjacency via flat-index modulo arithmetic\n\n";

    // Layer 4: Thermodynamic Scheduler
    std::cout << "Layer 4: Initializing thermodynamic scheduler (radix heap)...\n";
    ThermodynamicScheduler sched(lattice, d);
    std::cout << "  Radix heap: " << sched.heap.size() << " entries\n";
    std::cout << "  Tick period: dt_local = (L_t / L_rest) * dt_univ\n\n";

    // Run batch
    std::cout << "Running thermodynamic scheduler (200 universal ticks)...\n";
    std::vector<double> t_hist, total_E_hist, mean_L_hist, mean_Omega_hist, active_frac_hist;
    sched.run_batch(200, t_hist, total_E_hist, mean_L_hist, mean_Omega_hist, active_frac_hist);
    int active = 0, halted = 0;
    for (int i = 0; i < lattice.size(); i++) {
        if (lattice.arena[i].is_active) active++; else halted++;
    }
    std::cout << "  Active: " << active << " Halted (E=0): " << halted << "\n";
    std::vector<int> tick_counts(lattice.size());
    for (int i = 0; i < lattice.size(); i++) tick_counts[i] = lattice.arena[i].local_tick_count;
    double mean_ticks = std::accumulate(tick_counts.begin(), tick_counts.end(), 0.0) / lattice.size();
    double std_ticks = 0;
    for (int tc : tick_counts) std_ticks += (tc - mean_ticks) * (tc - mean_ticks);
    std_ticks = det_sqrt(std_ticks / lattice.size());
    std::cout << "  Tick count range: [" << *std::min_element(tick_counts.begin(), tick_counts.end())
              << ", " << *std::max_element(tick_counts.begin(), tick_counts.end()) << "]\n";
    std::cout << "  Tick count std: " << std_ticks << " (thermodynamic spread)\n\n";

    // Layer 5: Emergent Probes
    std::cout << "Layer 5: Running EMERGENT probes (no bridge equations)...\n";

    // Born
    std::cout << "  [Born] Sampling microcell distributions...\n";
    EmergentProbeBorn born_probe(d);
    double born_mean_tv, born_max_tv; std::string born_verdict;
    born_probe.run(1000, born_mean_tv, born_max_tv, born_verdict);
    std::cout << "    Mean TV: " << born_mean_tv << " Max TV: " << born_max_tv
              << " Verdict: " << born_verdict << "\n";

    // Bell
    std::cout << "  [Bell/CHSH] Probing torsion-link correlations...\n";
    EmergentProbeBell bell_probe(d);
    double S_local, S_qm, S_torsion, S_joint; std::string bell_verdict;
    bell_probe.run(S_local, S_qm, S_torsion, S_joint, bell_verdict);
    std::cout << "    S_local: " << S_local << " S_joint: " << S_joint
              << " Tsirelson: " << 2*det_sqrt(2.0) << " Verdict: " << bell_verdict << "\n";

    // Fine structure
    std::cout << "  [Fine structure] Torsion snap at UV cutoff...\n";
    EmergentProbeFineStructure fs_probe(d);
    double fs_one_over, fs_rel_diff; std::string fs_verdict;
    fs_probe.run(fs_one_over, fs_rel_diff, fs_verdict);
    std::cout << "    1/alpha_bare: " << fs_one_over
              << " CODATA: 137.035999084 Verdict: " << fs_verdict << "\n";

    // Gravity
    std::cout << "  [Gravity] Central cluster + test daemons...\n";
    Z3Lattice grav_lattice(10, 10, 4, d, 42);
    ThermodynamicScheduler grav_sched(grav_lattice, d);
    EmergentProbeGravity grav_probe;
    std::vector<double> grav_r, grav_v, grav_L;
    grav_probe.run(grav_lattice, grav_sched, d, grav_r, grav_v, grav_L);
    if (!grav_r.empty()) {
        std::cout << "    Radii: ";
        for (auto r : grav_r) std::cout << r << " ";
        std::cout << "\n    Emergent v(r): ";
        for (auto v : grav_v) std::cout << v << " ";
        std::cout << " km/s\n";
    }

    // Time dilation
    std::cout << "  [Time dilation] Tick rate gradient...\n";
    Z3Lattice td_lattice(8, 8, 4, d, 99);
    ThermodynamicScheduler td_sched(td_lattice, d);
    EmergentProbeTimeDilation td_probe;
    double gamma_near = 0, gamma_far = 0, tick_ratio = 0;
    int ticks_near = 0, ticks_far = 0;
    td_probe.run(td_lattice, td_sched, d, gamma_near, gamma_far, tick_ratio,
                 ticks_near, ticks_far);
    std::cout << "    Near gamma: " << gamma_near << " (" << ticks_near << " ticks)\n";
    std::cout << "    Far gamma:  " << gamma_far << " (" << ticks_far << " ticks)\n";
    std::cout << "    Tick ratio: " << tick_ratio << "\n\n";

    // Write JSON
    std::filesystem::create_directories("artifacts");
    std::string json_path = "artifacts/fpm_axcore_results.json";
    write_json(json_path, ax, d, cal,
               t_hist, total_E_hist, mean_L_hist, mean_Omega_hist, active_frac_hist,
               S_local, S_qm, S_torsion, S_joint, bell_verdict,
               born_mean_tv, born_max_tv, born_verdict,
               fs_one_over, fs_rel_diff, fs_verdict,
               grav_r, grav_v, grav_L,
               gamma_near, gamma_far, tick_ratio, ticks_near, ticks_far,
               lattice.sx, lattice.sy, lattice.sz, active, halted,
               (int)t_hist.size());
    std::cout << "Results JSON saved to: " << json_path << "\n\n";

    std::cout << "FPM v7.0-axcore-cpp emergent lattice simulation complete.\n\n";
    std::cout << "Architectural summary:\n";
    std::cout << "  OLD: Python objects + 1D ring + uniform loop + bridge equation\n";
    std::cout << "  NEW: flat arena + Z^3 modulo + radix heap + emergent probe\n\n";
    std::cout << "  The substrate BUILDS the universe.\n";
    std::cout << "  The scheduler ENFORCES the thermodynamic law.\n";
    std::cout << "  The probes MEASURE what emerges.\n";
    std::cout << "  N_bit_eq is the EXACT INTEGER. Zero rounding leak. Zero patches.\n\n";
    std::cout << "Closure: the universe becomes solid, directional, heavy,\n";
    std::cout << "time-slowed, structured, and stable for one basic reason:\n";
    std::cout << "KEEPING EVERYTHING OPEN IS TOO EXPENSIVE.\n";

    return 0;
}
