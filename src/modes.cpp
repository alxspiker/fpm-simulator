#include "modes.hpp"

bool read_file_string(const std::string& path, std::string& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return true;
}

bool parse_json_number_after_key(const std::string& text, size_t start,
                                        const std::string& key, double& value,
                                        size_t& value_end) {
    std::string needle = "\"" + key + "\"";
    size_t key_pos = text.find(needle, start);
    if (key_pos == std::string::npos) return false;
    size_t colon = text.find(':', key_pos + needle.size());
    if (colon == std::string::npos) return false;
    size_t pos = colon + 1;
    while (pos < text.size() && (text[pos] == ' ' || text[pos] == '\t' || text[pos] == '\n' || text[pos] == '\r')) pos++;
    char* end_ptr = nullptr;
    value = std::strtod(text.c_str() + pos, &end_ptr);
    if (end_ptr == text.c_str() + pos) return false;
    value_end = (size_t)(end_ptr - text.c_str());
    return true;
}

bool parse_sparc_injection_payload(const std::string& path, SparcInjectionPayload& payload) {
    std::string text;
    if (!read_file_string(path, text)) return false;

    size_t pos = 0;
    while (true) {
        double r = 0.0, v = 0.0, b = 0.0;
        size_t end_r = 0, end_v = 0, end_b = 0;
        if (!parse_json_number_after_key(text, pos, "r_kpc", r, end_r)) break;
        if (!parse_json_number_after_key(text, end_r, "v_obs", v, end_v)) return false;
        if (!parse_json_number_after_key(text, end_v, "b_load", b, end_b)) return false;
        if (r > 0.0 && v > 0.0 && b >= 0.0) {
            payload.points.push_back({r, v, b});
        }
        pos = end_b;
    }
    return !payload.points.empty();
}

void inject_core_load(Z3Lattice& lattice, int center, double B_load, double E_inject) {
    double radius = 2.0;
    for (int idx = 0; idx < lattice.size(); idx++) {
        double dist = lattice.euclidean_dist(idx, center);
        if (dist <= radius) {
            Daemon& dm = lattice.arena[idx];
            dm.E = std::min(E_inject, lattice.d.E_max);
            dm.R[0][1] += B_load;
            dm.R[1][0] -= B_load;
            dm.R[0][0] += 0.06 * B_load;
            dm.R[1][1] += 0.06 * B_load;
            dm.R[2][2] += 0.03 * B_load;
        }
    }
}

double measure_shell_L(Z3Lattice& lattice, ThermodynamicScheduler& sched,
                              const DerivedConstants& d, int center, double shell_r) {
    double sum_L = 0.0;
    int n = 0;
    for (int idx = 0; idx < lattice.size(); idx++) {
        double dist = lattice.euclidean_dist(idx, center);
        if (std::abs(dist - shell_r) <= 0.7 && lattice.arena[idx].is_active) {
            Daemon& dm = lattice.arena[idx];
            double O, k, C_N;
            viscosity_update(dm, d, 0.0, O, k, C_N);
            double C_sem, C_geo, sm;
            sum_L += axcore_lagrangian(dm, d, O, sched.cfg, C_sem, C_geo, sm);
            n++;
        }
    }
    return n > 0 ? std::max(d.L_rest, sum_L / n) : d.L_rest;
}

int run_sparc_payload_mode(const std::string& payload_path, const std::string& output_path) {
    SparcInjectionPayload payload;
    if (!parse_sparc_injection_payload(payload_path, payload)) {
        std::cerr << "ERROR: failed to parse sanitized SPARC injection payload: " << payload_path << "\n";
        return 2;
    }

    Axioms ax;
    DerivedConstants d = derive_all(ax);

    double max_r = 0.0;
    for (const auto& p : payload.points) max_r = std::max(max_r, p.r_kpc);
    if (max_r <= 0.0) return 3;

    std::vector<double> raw_v, shell_L, shell_radii;
    raw_v.reserve(payload.points.size());
    shell_L.reserve(payload.points.size());
    shell_radii.reserve(payload.points.size());

    for (size_t pi = 0; pi < payload.points.size(); pi++) {
        const auto& p = payload.points[pi];
        Z3Lattice lattice(32, 32, 4, d, 7200 + (int)pi);
        ThermodynamicScheduler sched(lattice, d);
        int center = lattice.flat(lattice.sx/2, lattice.sy/2, lattice.sz/2);
        double shell_r = 3.0 + 11.0 * p.r_kpc / max_r;
        shell_radii.push_back(shell_r);

        inject_core_load(lattice, center, std::min(24.0, p.b_load), 0.66);

        for (int tick = 0; tick < 70; tick++) {
            for (int i = 0; i < lattice.size(); i++) {
                double L, O, k;
                sched.step_one(L, O, k);
            }
            if (tick % 10 == 0) sched.mean_field_truth_target();
        }

        double L_mean = measure_shell_L(lattice, sched, d, center, shell_r);
        shell_L.push_back(L_mean);

        double gradient = 1.0 + 0.015 * std::min(24.0, p.b_load);
        raw_v.push_back((d.L_rest / L_mean) * gradient);
    }

    std::ofstream f(output_path);
    if (!f) {
        std::cerr << "ERROR: failed to open SPARC substrate output: " << output_path << "\n";
        return 4;
    }
    f << std::setprecision(15);
    f << "{\n";
    f << "  \"mode\": \"sparc_substrate\",\n";
    f << "  \"lattice\": [32,32,4],\n";
    f << "  \"points\": " << payload.points.size() << ",\n";
    f << "  \"r_kpc\": [";
    for (size_t i = 0; i < payload.points.size(); i++) { if (i) f << ","; f << payload.points[i].r_kpc; }
    f << "],\n";
    f << "  \"v_obs_kms\": [";
    for (size_t i = 0; i < payload.points.size(); i++) { if (i) f << ","; f << payload.points[i].v_obs; }
    f << "],\n";
    f << "  \"b_load\": [";
    for (size_t i = 0; i < payload.points.size(); i++) { if (i) f << ","; f << payload.points[i].b_load; }
    f << "],\n";
    f << "  \"shell_radius\": [";
    for (size_t i = 0; i < shell_radii.size(); i++) { if (i) f << ","; f << shell_radii[i]; }
    f << "],\n";
    f << "  \"L_mean\": [";
    for (size_t i = 0; i < shell_L.size(); i++) { if (i) f << ","; f << shell_L[i]; }
    f << "],\n";
    f << "  \"v_fpm_raw\": [";
    for (size_t i = 0; i < raw_v.size(); i++) { if (i) f << ","; f << raw_v[i]; }
    f << "]\n";
    f << "}\n";
    std::cout << "SPARC substrate output saved to: " << output_path << "\n";
    return 0;
}

// =============================================================================
// EMERGENT DYNAMIC TORSION (THE ZOMBIE SNAP)
// =============================================================================

int run_dynamic_torsion_mode(const std::string& output_path) {
    Axioms ax;
    DerivedConstants d = derive_all(ax);
    Z3Lattice lattice(16, 16, 4, d, 42);
    ThermodynamicScheduler sched(lattice, d);

    int idx_A = lattice.flat(2, 8, 2);
    int idx_B = lattice.flat(13, 8, 2);

    double scale = 1.0;
    double A_val[3][3] = {{0,0,0},{0,0,-scale},{0,scale,0}};
    double A_val_T[3][3] = {{0,0,0},{0,0,scale},{0,-scale,0}};

    Daemon& da = lattice.arena[idx_A];
    Daemon& db = lattice.arena[idx_B];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double S_ij_a = 0.5 * (da.R[i][j] + da.R[j][i]);
            da.R[i][j] = S_ij_a + A_val[i][j];
            double S_ij_b = 0.5 * (db.R[i][j] + db.R[j][i]);
            db.R[i][j] = S_ij_b + A_val_T[i][j];
        }
    }

    EmergentProbeBell probe(d);
    std::vector<int> ticks_out;
    std::vector<double> energy_out;
    std::vector<double> s_joint_out;
    std::vector<bool> zombie_out;
    std::vector<std::string> mode_out;
    int snap_tick = -1;

    ticks_out.reserve(150);
    energy_out.reserve(150);
    s_joint_out.reserve(150);
    zombie_out.reserve(150);
    mode_out.reserve(150);

    for (int tick = 0; tick < 150; tick++) {
        double current_E = 0.0;
        for (int i = 0; i < lattice.size(); i++) {
            lattice.arena[i].E *= 0.85;
            current_E += lattice.arena[i].E;
        }
        current_E /= (double)lattice.size();

        bool is_zombie = (da.E < 0.01 && db.E < 0.01);

        for (int i = 0; i < lattice.size(); i++) {
            double L, O, k;
            sched.step_one(L, O, k);
        }

        auto dynamic_correlation = [&](double a, double b) {
            double local_corr = probe.local_torsion_correlation(a, b);
            double lrm_corr = probe.joint_torsion_lrm_correlation(a, b);
            double blend = std::max(0.0, std::min(1.0, current_E / 0.10));
            if (is_zombie) return lrm_corr;
            return blend * local_corr + (1.0 - blend) * lrm_corr;
        };

        double S_joint = std::abs(dynamic_correlation(0.0, PI/4.0) +
                                  dynamic_correlation(0.0, -PI/4.0) +
                                  dynamic_correlation(PI/2.0, PI/4.0) -
                                  dynamic_correlation(PI/2.0, -PI/4.0));

        std::string mode = "FLOW";
        if (is_zombie) mode = "ZOMBIE";
        else if (current_E < 0.10) mode = "FATIGUE";

        ticks_out.push_back(tick);
        energy_out.push_back(current_E);
        s_joint_out.push_back(S_joint);
        zombie_out.push_back(is_zombie);
        mode_out.push_back(mode);

        if (is_zombie && snap_tick == -1) {
            snap_tick = tick;
        }
    }

    std::ofstream f(output_path);
    if (!f) {
        std::cerr << "ERROR: failed to open torsion output: " << output_path << "\n";
        return 5;
    }

    f << std::setprecision(15);
    f << "{\n";
    f << "  \"mode\": \"dynamic_torsion\",\n";
    f << "  \"lattice\": [16,16,4],\n";
    f << "  \"node_A\": " << idx_A << ",\n";
    f << "  \"node_B\": " << idx_B << ",\n";
    f << "  \"snap_tick\": " << snap_tick << ",\n";

    f << "  \"ticks\": [";
    for (size_t i = 0; i < ticks_out.size(); i++) { if (i) f << ","; f << ticks_out[i]; }
    f << "],\n";

    f << "  \"energy\": [";
    for (size_t i = 0; i < energy_out.size(); i++) { if (i) f << ","; f << energy_out[i]; }
    f << "],\n";

    f << "  \"S_joint\": [";
    for (size_t i = 0; i < s_joint_out.size(); i++) { if (i) f << ","; f << s_joint_out[i]; }
    f << "],\n";

    f << "  \"is_zombie\": [";
    for (size_t i = 0; i < zombie_out.size(); i++) { if (i) f << ","; f << (zombie_out[i] ? "true" : "false"); }
    f << "],\n";

    f << "  \"samples\": [\n";
    for (size_t i = 0; i < ticks_out.size(); i++) {
        if (i) f << ",\n";
        f << "    {\"tick\": " << ticks_out[i]
          << ", \"mean_E\": " << energy_out[i]
          << ", \"S_joint\": " << s_joint_out[i]
          << ", \"mode\": \"" << mode_out[i] << "\"}";
    }
    f << "\n  ]\n";
    f << "}\n";

    std::cout << "Dynamic torsion phase-lock output saved to: " << output_path << "\n";
    return 0;
}

// =============================================================================
// EMERGENT PHOTON PROPAGATION (CAUSAL WAVE LIMIT)
// =============================================================================

int run_photon_propagation_mode(const std::string& output_path) {
    Axioms ax;
    DerivedConstants d = derive_all(ax);
    
    // Z^3 Corridor: 100 cells long, 4x4 cross-section
    int len_x = 100;
    Z3Lattice lattice(len_x, 4, 4, d, 1337);
    ThermodynamicScheduler sched(lattice, d);

    // THERMALIZATION: Scatter the background vacuum and misalign timestamps
    for (int t = 0; t < 5; t++) {
        uint64_t thermal_target = sched.universal_step + ACTION_STEPS_PER_UNIVERSAL_TICK;
        while (sched.universal_step < thermal_target && !sched.heap.empty()) {
            double L, O, k;
            sched.step_one(L, O, k);
        }
    }

    // Inject the causal wave (photon) at x = 5
    int inject_x = 5;
    uint64_t sync_time = sched.universal_step + ACTION_STEPS_PER_UNIVERSAL_TICK;
    for (int y = 0; y < 4; y++) {
        for (int z = 0; z < 4; z++) {
            int idx = lattice.flat(inject_x, y, z);
            Daemon& dm = lattice.arena[idx];
            dm.E = d.E_max;
            dm.R[0][0] += 5.0; // Pure longitudinal forward momentum
            dm.is_active = true;
            
            // Synchronize the injection to form a coherent planar wave
            dm.locked_until = sync_time;
            sched.heap.push({sync_time, idx});
        }
    }

    std::vector<int> ticks_out;
    std::vector<double> peak_x_out;
    std::vector<double> leading_edge_out;
    std::vector<double> total_arena_energy;

    uint64_t current_target_step = sched.universal_step;

    // Track the propagation over 80 universal ticks
    for (int tick = 0; tick < 80; tick++) {
        
        // Advance the target ledger by exactly 1 universal tick
        current_target_step += ACTION_STEPS_PER_UNIVERSAL_TICK;
        
        // Execute strictly by the thermodynamic clock, not an arbitrary array count
        while (sched.universal_step < current_target_step && !sched.heap.empty()) {
            double L, O, k;
            sched.step_one(L, O, k);
        }
        
        sched.mean_field_truth_target();

        // Scan the corridor to find the photon's location
        double max_E = 0.0;
        double peak_x = 0.0;
        double leading_edge = 0.0;
        double current_total_E = 0.0;

        // Profile the cross-section of the corridor
        for (int x = 0; x < len_x; x++) {
            double slice_E = 0.0;
            for (int y = 0; y < 4; y++) {
                for (int z = 0; z < 4; z++) {
                    slice_E += lattice.arena[lattice.flat(x, y, z)].E;
                }
            }
            current_total_E += slice_E;
            
            if (slice_E > max_E) {
                max_E = slice_E;
                peak_x = (double)x;
            }
            
            // Define the leading edge as the furthest x-coordinate where the 
            // energy perturbation is strictly above the vacuum baseline (0.5)
            if (slice_E > 16.0 * 0.51) { // 16 nodes in a 4x4 slice, baseline is 0.5
                leading_edge = std::max(leading_edge, (double)x);
            }
        }

        ticks_out.push_back(tick);
        peak_x_out.push_back(peak_x);
        leading_edge_out.push_back(leading_edge);
        total_arena_energy.push_back(current_total_E);
    }

    // Output the telemetry
    std::ofstream f(output_path);
    if (!f) {
        std::cerr << "ERROR: failed to open photon output: " << output_path << "\n";
        return 6;
    }
    
    f << std::setprecision(15);
    f << "{\n";
    f << "  \"mode\": \"photon_propagation\",\n";
    f << "  \"lattice\": [" << len_x << ", 4, 4],\n";
    
    f << "  \"ticks\": [";
    for(size_t i = 0; i < ticks_out.size(); i++) { if(i) f << ","; f << ticks_out[i]; }
    f << "],\n";

    f << "  \"leading_edge_x\": [";
    for(size_t i = 0; i < leading_edge_out.size(); i++) { if(i) f << ","; f << leading_edge_out[i]; }
    f << "],\n";
    
    f << "  \"peak_x\": [";
    for(size_t i = 0; i < peak_x_out.size(); i++) { if(i) f << ","; f << peak_x_out[i]; }
    f << "]\n";
    f << "}\n";

    std::cout << "Photon causal wave output saved to: " << output_path << "\n";
    return 0;
}

// =============================================================================
// CALIBRATION
// =============================================================================



CalibrationResult calibrate(const DerivedConstants& d) {
    CalibrationResult c;
    c.rel_err_G = std::abs(d.G_FPM - G_CODATA) / G_CODATA;
    c.rel_err_ns = std::abs(d.n_s - PLANCK_NS) / PLANCK_NS;
    c.ell_D_in_range = (1100.0 <= d.ell_D && d.ell_D <= 1500.0);
    c.gamma_above_cern = (d.gamma_max > CERN_MUON_GAMMA);
    c.N_bit_eq_exact = d.N_bit_eq;
    c.N_bit_eq_continuous = (4.0 * PI / 3.0) * d.alpha_PP * d.alpha_PP * d.alpha_PP;
    return c;
}

// =============================================================================
// JSON OUTPUT
// =============================================================================

void write_json(const std::string& path,
                const Axioms& ax, const DerivedConstants& d,
                const CalibrationResult& cal,
                const std::vector<double>& t_hist,
                const std::vector<double>& total_E_hist,
                const std::vector<double>& mean_L_hist,
                const std::vector<double>& mean_Omega_hist,
                const std::vector<double>& active_frac_hist,
                double S_local, double S_qm, double S_torsion, double S_joint,
                const std::string& bell_verdict,
                double born_mean_tv, double born_max_tv, const std::string& born_verdict,
                double fs_one_over, double fs_rel_diff, const std::string& fs_verdict,
                const std::vector<double>& grav_r, const std::vector<double>& grav_v,
                const std::vector<double>& grav_L,
                double gamma_near, double gamma_far, double tick_ratio,
                int ticks_near, int ticks_far,
                int lattice_sx, int lattice_sy, int lattice_sz,
                int n_daemons_active, int n_daemons_halted,
                int events_processed) {
    std::ofstream f(path);
    f << std::setprecision(15);
    f << "{\n";
    f << "  \"metadata\": { \"version\": \"v7.0-axcore-cpp\", "
       << "\"architecture\": \"Z^3 flat arena + radix heap scheduler + emergent probes\" },\n";

    // Axioms
    f << "  \"axioms\": { \"dim_space\": " << ax.dim_space
      << ", \"dim_causal\": " << ax.dim_causal
      << ", \"n_directed\": " << ax.n_directed()
      << ", \"n_trace\": " << ax.n_trace() << " },\n";

    // Derived constants
    f << "  \"derived_constants\": {\n";
    f << "    \"alpha\": " << d.alpha << ", \"beta\": " << d.beta
      << ", \"chi_arrow\": " << d.chi_arrow << ",\n";
    f << "    \"Omega_min\": " << d.Omega_min << ", \"Omega_max\": " << d.Omega_max
      << ", \"E_max\": " << d.E_max << ",\n";
    f << "    \"c0\": " << d.c0 << ", \"L_max\": " << d.L_max
      << ", \"L_rest\": " << d.L_rest << ",\n";
    f << "    \"gamma_max\": " << d.gamma_max << ", \"alpha_PP\": " << d.alpha_PP << ",\n";
    f << "    \"N_bit_eq\": " << d.N_bit_eq << ",\n";
    f << "    \"A_FPM\": " << d.A_FPM << ", \"n_s\": " << d.n_s
      << ", \"r_tensor\": " << d.r_tensor << ",\n";
    f << "    \"ell_D\": " << d.ell_D << ", \"dt_univ\": " << d.dt_univ
      << ", \"dx_univ\": " << d.dx_univ << ",\n";
    f << "    \"G_FPM\": " << d.G_FPM << ", \"calib\": " << d.calib << "\n";
    f << "  },\n";

    // Calibration
    f << "  \"calibration\": {\n";
    f << "    \"G_FPM_rel_err_pct\": " << cal.rel_err_G * 100 << ",\n";
    f << "    \"n_s_rel_err_pct\": " << cal.rel_err_ns * 100 << ",\n";
    f << "    \"ell_D_in_range\": " << (cal.ell_D_in_range ? "true" : "false") << ",\n";
    f << "    \"gamma_above_cern_muon\": " << (cal.gamma_above_cern ? "true" : "false") << ",\n";
    f << "    \"N_bit_eq_exact_integer\": " << cal.N_bit_eq_exact << ",\n";
    f << "    \"N_bit_eq_continuous_approx\": " << cal.N_bit_eq_continuous << ",\n";
    f << "    \"rounding_leak_relative\": "
      << std::abs((double)cal.N_bit_eq_exact - cal.N_bit_eq_continuous) / cal.N_bit_eq_continuous << "\n";
    f << "  },\n";

    // Lattice
    f << "  \"lattice\": { \"size\": [" << lattice_sx << "," << lattice_sy << "," << lattice_sz << "],\n";
    f << "    \"total_daemons\": " << lattice_sx * lattice_sy * lattice_sz << ",\n";
    f << "    \"active\": " << n_daemons_active << ", \"halted\": " << n_daemons_halted << ",\n";
    f << "    \"adjacency\": \"6-face Z^3 periodic\",\n";
    f << "    \"scheduler_events_processed\": " << events_processed << "\n";
    f << "  },\n";

    // Emergent probes
    f << "  \"emergent_probes\": {\n";
    f << "    \"born\": { \"mean_tv\": " << born_mean_tv << ", \"max_tv\": " << born_max_tv
      << ", \"verdict\": \"" << born_verdict << "\" },\n";
    f << "    \"bell\": { \"S_local\": " << S_local << ", \"S_qm\": " << S_qm
      << ", \"S_torsion\": " << S_torsion << ", \"S_joint\": " << S_joint
      << ", \"tsirelson\": " << 2.0*det_sqrt(2.0)
      << ", \"verdict\": \"" << bell_verdict << "\" },\n";
    f << "    \"fine_structure\": { \"one_over_alpha_bare\": " << fs_one_over
      << ", \"rel_diff_from_macro\": " << fs_rel_diff
      << ", \"verdict\": \"" << fs_verdict << "\" },\n";
    f << "    \"gravity\": { \"r\": [";
    for (size_t i = 0; i < grav_r.size(); i++) { if(i) f << ","; f << grav_r[i]; }
    f << "], \"v_emergent_kms\": [";
    for (size_t i = 0; i < grav_v.size(); i++) { if(i) f << ","; f << grav_v[i]; }
    f << "], \"L_mean\": [";
    for (size_t i = 0; i < grav_L.size(); i++) { if(i) f << ","; f << grav_L[i]; }
    f << "] },\n";
    f << "    \"time_dilation\": { \"gamma_near\": " << gamma_near
      << ", \"gamma_far\": " << gamma_far
      << ", \"tick_ratio\": " << tick_ratio
      << ", \"ticks_near\": " << ticks_near
      << ", \"ticks_far\": " << ticks_far << " }\n";
    f << "  },\n";

    // Trajectory
    f << "  \"trajectory\": { \"t\": [";
    for (size_t i = 0; i < t_hist.size(); i++) { if(i) f << ","; f << t_hist[i]; }
    f << "], \"total_E\": [";
    for (size_t i = 0; i < total_E_hist.size(); i++) { if(i) f << ","; f << total_E_hist[i]; }
    f << "], \"mean_L\": [";
    for (size_t i = 0; i < mean_L_hist.size(); i++) { if(i) f << ","; f << mean_L_hist[i]; }
    f << "], \"mean_Omega\": [";
    for (size_t i = 0; i < mean_Omega_hist.size(); i++) { if(i) f << ","; f << mean_Omega_hist[i]; }
    f << "], \"active_fraction\": [";
    for (size_t i = 0; i < active_frac_hist.size(); i++) { if(i) f << ","; f << active_frac_hist[i]; }
    f << "] }\n";

    f << "}\n";
    f.close();
}



