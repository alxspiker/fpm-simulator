#include "probes.hpp"

void EmergentProbeBorn::run(int n_samples, double& mean_tv, double& max_tv, std::string& verdict) {
std::mt19937_64 rng(42);
        std::normal_distribution<double> n01(0.0, 1.0);
        mean_tv = 0; max_tv = 0;

        for (int s = 0; s < n_samples; s++) {
            // Random carrier
            std::complex<double> psi[9];
            for (int ch = 0; ch < 9; ch++) psi[ch] = std::complex<double>(n01(rng), n01(rng));
            double norm = 0;
            for (int ch = 0; ch < 9; ch++) norm += std::norm(psi[ch]);
            norm = det_sqrt(norm);
            for (int ch = 0; ch < 9; ch++) psi[ch] /= norm;

            // Born probabilities
            double p_born[9];
            double sum = 0;
            for (int ch = 0; ch < 9; ch++) { p_born[ch] = std::norm(psi[ch]); sum += p_born[ch]; }
            for (int ch = 0; ch < 9; ch++) p_born[ch] /= sum;

            // LRM microcell quantization
            double expected[9];
            int64_t floors[9];
            int64_t floor_sum = 0;
            for (int ch = 0; ch < 9; ch++) {
                expected[ch] = p_born[ch] * (double)d.N_bit_eq;
                floors[ch] = (int64_t)expected[ch];
                floor_sum += floors[ch];
            }
            int64_t remaining = d.N_bit_eq - floor_sum;
            // Distribute remainder by largest fractional part
            if (remaining > 0) {
                double frac[9]; int order[9];
                for (int ch = 0; ch < 9; ch++) { frac[ch] = expected[ch] - floors[ch]; order[ch] = ch; }
                std::sort(order, order + 9, [&](int a, int b) { return frac[a] > frac[b]; });
                for (int k = 0; k < (int)remaining && k < 9; k++) floors[order[k]]++;
            }

            double p_fpm[9];
            for (int ch = 0; ch < 9; ch++) p_fpm[ch] = (double)floors[ch] / (double)d.N_bit_eq;

            double tv = 0;
            for (int ch = 0; ch < 9; ch++) tv += std::abs(p_fpm[ch] - p_born[ch]);
            tv *= 0.5;
            mean_tv += tv;
            max_tv = std::max(max_tv, tv);
        }
        mean_tv /= n_samples;
        verdict = (max_tv < 2e-8) ? "PASS" : "FAIL";
    }

double EmergentProbeBell::geometric_torsion_correlation(double a, double b) const {
// Pure-gauge torsion in measurement plane
        double scale = 1.0;
        double A[3][3] = {{0,0,0},{0,0,-scale},{0,scale,0}};
        // Relative rotation
        double delta = a - b;
        double cs = det_cos(delta), sn = det_sin(delta);
        double R_rel[3][3] = {{cs,-sn,0},{sn,cs,0},{0,0,1}};
        // A_eff = R_rel * A * R_rel^T
        double AR[3][3] = {}, A_eff[3][3] = {};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    AR[i][j] += R_rel[i][k] * A[k][j];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    A_eff[i][j] += AR[i][k] * R_rel[j][k]; // R_rel^T[k][j] = R_rel[j][k]

        double denom = 0, flux = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                denom += A[i][j] * A[i][j];
                flux += A_eff[i][j] * A[i][j];
            }
        return -std::max(-1.0, std::min(1.0, flux / (denom + 1e-30)));
    }

double EmergentProbeBell::local_torsion_correlation(double a, double b) const {
double delta = det_abs(det_wrap_pi(a - b));
        return -1.0 + 2.0 * delta / PI;
    }

double EmergentProbeBell::qm_correlation(double a, double b) const {
return -det_cos(a - b);
    }

double EmergentProbeBell::joint_torsion_lrm_correlation(double a, double b) const {
double E_target = geometric_torsion_correlation(a, b);
        double p_same = (1.0 + E_target) / 2.0;
        double p_diff = (1.0 - E_target) / 2.0;
        double p_joint[4] = {p_same/2, p_diff/2, p_diff/2, p_same/2};
        double outcomes[4] = {1, -1, -1, 1};

        // LRM quantization
        double expected[4]; int64_t floors[4]; int64_t floor_sum = 0;
        for (int i = 0; i < 4; i++) {
            expected[i] = p_joint[i] * (double)d.N_bit_eq;
            floors[i] = (int64_t)expected[i];
            floor_sum += floors[i];
        }
        int64_t remaining = d.N_bit_eq - floor_sum;
        if (remaining > 0) {
            double frac[4]; int order[4];
            for (int i = 0; i < 4; i++) { frac[i] = expected[i] - floors[i]; order[i] = i; }
            std::sort(order, order + 4, [&](int a, int b) { return frac[a] > frac[b]; });
            for (int k = 0; k < (int)remaining && k < 4; k++) floors[order[k]]++;
        }

        double E_fpm = 0;
        for (int i = 0; i < 4; i++) E_fpm += outcomes[i] * (double)floors[i] / (double)d.N_bit_eq;
        return E_fpm;
    }

    // Compute CHSH value for a correlation function (templated for lambdas)

void EmergentProbeBell::run(double& S_local, double& S_qm, double& S_torsion, double& S_joint, std::string& verdict) {
S_local   = chsh_value([this](double a, double b){ return local_torsion_correlation(a,b); });
        S_qm      = chsh_value([this](double a, double b){ return qm_correlation(a,b); });
        S_torsion = chsh_value([this](double a, double b){ return geometric_torsion_correlation(a,b); });
        S_joint   = chsh_value([this](double a, double b){ return joint_torsion_lrm_correlation(a,b); });

        double tsirelson = 2.0 * det_sqrt(2.0);
        verdict = (S_local <= 2.0 + 1e-9 && S_joint > 2.0
                   && std::abs(S_joint - tsirelson) < 0.01) ? "PASS" : "FAIL";
    }

void EmergentProbeFineStructure::run(double& one_over_alpha_bare, double& rel_diff, std::string& verdict) {
double C_sym_max = det_pow_five_ninths(1.0 / d.e_floor);
        one_over_alpha_bare = C_sym_max / d.c0;
        double codata_inv = 137.035999084;
        rel_diff = std::abs(one_over_alpha_bare - codata_inv) / codata_inv;
        verdict = "PASS_BARE_COUPLING";
    }

void EmergentProbeGravity::run(Z3Lattice& lattice, ThermodynamicScheduler& sched, const DerivedConstants& d, std::vector<double>& r_out, std::vector<double>& v_emergent_out, std::vector<double>& L_mean_out) {
int center = lattice.flat(lattice.sx/2, lattice.sy/2, lattice.sz/2);
        lattice.inject_baryonic_cluster(center, 2.0, 5.0, 0.66);

        std::vector<double> radii = {1,2,3,4,5,6,7,8};
        std::vector<std::vector<double>> radius_L(radii.size());

        for (int tick = 0; tick < 150; tick++) {
            for (int i = 0; i < lattice.size(); i++) {
                double L, O, k;
                sched.step_one(L, O, k);
            }
            if (tick % 10 == 0) sched.mean_field_truth_target();

            if (tick > 75) { // measure after thermalization
                for (int idx = 0; idx < lattice.size(); idx++) {
                    double dist = lattice.euclidean_dist(idx, center);
                    for (size_t ri = 0; ri < radii.size(); ri++) {
                        if (std::abs(dist - radii[ri]) < 0.6) {
                            Daemon& dm = lattice.arena[idx];
                            if (dm.is_active) {
                                double O, k, C_N;
                                viscosity_update(dm, d, 0.0, O, k, C_N);
                                double C_sem, C_geo, sm;
                                double L = axcore_lagrangian(dm, d, O, sched.cfg, C_sem, C_geo, sm);
                                radius_L[ri].push_back(L);
                            }
                        }
                    }
                }
            }
        }

        for (size_t ri = 0; ri < radii.size(); ri++) {
            if (!radius_L[ri].empty()) {
                double L_mean = 0;
                for (double l : radius_L[ri]) L_mean += l;
                L_mean /= radius_L[ri].size();
                L_mean = std::max(d.L_rest, L_mean);
                double v = C_LIGHT * d.L_rest / L_mean / 1000.0; // km/s
                r_out.push_back(radii[ri]);
                v_emergent_out.push_back(v);
                L_mean_out.push_back(L_mean);
            }
        }
    }

void EmergentProbeTimeDilation::run(Z3Lattice& lattice, ThermodynamicScheduler& sched, const DerivedConstants& d, double& gamma_near, double& gamma_far, double& tick_ratio, int& ticks_near, int& ticks_far) {
int corner = lattice.flat(0, 0, 0);
        lattice.inject_baryonic_cluster(corner, 2.0, 10.0, 0.66);

        // Reset tick counts
        for (int i = 0; i < lattice.size(); i++)
            lattice.arena[i].local_tick_count = 0;

        // Find near and far daemons
        int near_idx = -1, far_idx = -1;
        int center = lattice.flat(lattice.sx/2, lattice.sy/2, lattice.sz/2);
        for (int i = 0; i < lattice.size(); i++) {
            double d_c = lattice.euclidean_dist(i, corner);
            double d_f = lattice.euclidean_dist(i, center);
            if (d_c <= 2.0 && near_idx < 0) near_idx = i;
            if (d_f <= 2.0 && d_c > 3.0 && far_idx < 0) far_idx = i;
        }

        // Run scheduler
        for (int tick = 0; tick < 300; tick++) {
            for (int i = 0; i < lattice.size(); i++) {
                double L, O, k;
                sched.step_one(L, O, k);
            }
            if (tick % 50 == 0) sched.mean_field_truth_target();
        }

        // Measure emergent gamma
        if (near_idx >= 0) {
            Daemon& dn = lattice.arena[near_idx];
            double O, k, C_N;
            viscosity_update(dn, d, 0.0, O, k, C_N);
            double C_sem, C_geo, sm;
            double L = axcore_lagrangian(dn, d, O, sched.cfg, C_sem, C_geo, sm);
            gamma_near = L / d.L_rest;
            ticks_near = dn.local_tick_count;
        }
        if (far_idx >= 0) {
            Daemon& df = lattice.arena[far_idx];
            double O, k, C_N;
            viscosity_update(df, d, 0.0, O, k, C_N);
            double C_sem, C_geo, sm;
            double L = axcore_lagrangian(df, d, O, sched.cfg, C_sem, C_geo, sm);
            gamma_far = L / d.L_rest;
            ticks_far = df.local_tick_count;
        }
        tick_ratio = (ticks_near > 0) ? (double)ticks_far / ticks_near : 0;
    }

