#include "constants.hpp"

// AxCore cost function
double axcore_cost(double H, double S, double f, double kappa_strat, const Axioms& ax) {
    double Bdt = ax.axcore_dt_offset + ax.axcore_dt_H_coef * H + ax.axcore_dt_S_coef * S;
    Bdt = std::max(ax.axcore_dt_clip_low, std::min(ax.axcore_dt_clip_high, Bdt));
    double f_c = std::max(0.0, std::min(1.0, f));
    double base = ax.axcore_base + ax.axcore_critic_coef * Bdt;
    double critic = ax.axcore_fitness_coef * (1.0 - f_c);
    double total = (base + critic) * kappa_strat;
    return std::max(ax.axcore_min_floor, total);
}

// Exact integer N_bit_eq: count Z^3 lattice points in sphere of radius alpha_PP
int64_t compute_N_bit_eq_exact(double alpha_PP) {
    double R_sq = alpha_PP * alpha_PP;
    int R_int = (int)std::floor(alpha_PP);
    int64_t count = 0;
    for (int i = -R_int; i <= R_int; i++) {
        double i_sq = (double)(i * i);
        if (i_sq > R_sq) continue;
        int j_max = (int)std::floor(std::sqrt(R_sq - i_sq));
        for (int j = -j_max; j <= j_max; j++) {
            double ij_sq = i_sq + (double)(j * j);
            if (ij_sq > R_sq) continue;
            int k_max = (int)std::floor(std::sqrt(R_sq - ij_sq));
            count += 2 * k_max + 1;
        }
    }
    return count;
}

// Alpha_PP counterterm fixed-point iteration
double iterate_alpha_pp_counterterm(double a1, double a2, double L_rest,
                                     double tol, int max_iter) {
    double b1 = 1.5;
    double c2 = 7.5 - L_rest;
    double x = a2;
    for (int iter = 0; iter < max_iter; iter++) {
        double x_new = a1 + b1 / (x + 9.0) + c2 / ((x + 9.0) * (x + 9.0));
        double delta = std::abs(x_new - x);
        if (delta < tol * std::abs(x)) return x_new;
        x = x_new;
    }
    fprintf(stderr, "ERROR: alpha_PP counterterm failed to converge\n");
    exit(1);
}

DerivedConstants derive_all(const Axioms& ax) {
    DerivedConstants d;
    int nd = ax.n_directed(), nt = ax.n_trace();
    d.alpha = 2.0 * nt / (nd + nt);
    d.beta  = 2.0 * nd / (nd + nt);
    d.n_directed = nd;
    d.n_trace = nt;

    double e_floor_paper = 0.0314;
    double pc_iso = 0.2488, pc_dir = 0.50;
    d.chi_arrow = e_floor_paper / ((pc_dir - pc_iso) / 2.0);
    d.e_floor = e_floor_paper;
    d.e_exp = -0.75;
    d.Omega_min = 0.50;
    d.Omega_max = 0.85;

    double L_axcore_mean = (ax.axcore_base + ax.axcore_critic_coef * ax.bench_Bdt
                            + ax.axcore_fitness_coef * (1.0 - ax.bench_fitness))
                           * ax.bench_kappa_strat;
    d.calib = ax.dim_causal * L_axcore_mean; // 80
    d.c0 = 0.05;
    d.E_max = 2.0 * d.c0 / (1.0 - d.Omega_max);
    d.n_blade = 2;
    d.lam = (ax.dim_causal * d.n_directed) / (double)(d.n_directed - d.n_blade);

    double L_axcore_max = (ax.axcore_base + ax.axcore_critic_coef * ax.axcore_dt_clip_high
                           + ax.axcore_fitness_coef * 1.0) * 1.0;
    double C_sem_max = L_axcore_max / d.calib;
    d.L_max = 3.0 * C_sem_max + d.lam * 0.35;

    double factor = d.chi_arrow * (d.n_directed - d.n_blade) / (double)(d.n_directed + d.n_trace);
    double residual = factor * factor / (d.n_directed + d.n_trace);
    d.L_rest = 2.0 * d.c0 + residual;
    d.gamma_max = d.L_max / d.L_rest;
    d.ledger_inertia_ratio = (double)(ax.dim_causal * ax.dim_causal) / ax.dim_space;

    // Theorem 5: alpha_PP 4-step derivation
    int g_PP = 2;
    double C_le9 = 0;
    for (int n = 1; n <= 9; n++) C_le9 += g_PP * n * n; // 570
    double C10 = g_PP * 10 * 10; // 200
    double fT = 2.0 / 3.0;
    double a0 = C_le9 + fT * C10;
    double a1 = a0 - 1.0 / std::sqrt(2.0);
    double b1_val = 1.5;
    double a2 = (-(9.0 - a1) + std::sqrt((9.0 - a1) * (9.0 - a1) + 4.0 * (9.0 * a1 + b1_val))) / 2.0;
    d.alpha_PP = iterate_alpha_pp_counterterm(a1, a2, d.L_rest);

    // EXACT INTEGER N_bit_eq (the reason we are in C++)
    d.N_bit_eq = compute_N_bit_eq_exact(d.alpha_PP);
    d.A_FPM = (2.0 / 3.0) * std::sqrt(d.ledger_inertia_ratio / (double)d.N_bit_eq);
    d.n_s = 1.0 - d.L_rest / d.L_max;
    d.r_tensor = (1.0 / 9.0) * (d.L_rest / d.L_max);
    d.ell_A = 299.82;
    d.ell_freeze = 5720.0;
    d.ell_D = std::sqrt(d.ell_A * d.ell_freeze);
    d.dt_univ = ax.h_planck / (ax.m_e * ax.c_light * ax.c_light * d.alpha_PP);
    d.dx_univ = ax.c_light * d.dt_univ;
    d.f_univ = 1.0 / d.dt_univ;
    d.zeta = d.n_directed / (4.0 * PI * d.L_max);
    d.J_per_bit_eq = (double)d.N_bit_eq * ax.k_B * ax.T_substrate * std::log(2.0);
    double N_bit = (double)d.N_bit_eq;
    double N_bit_2 = N_bit * N_bit;
    double c2 = ax.c_light * ax.c_light;
    d.mu_M_FPM = (2.0 / 3.0) * d.zeta / ((d.alpha_PP + 9.0) * N_bit_2 * N_bit_2);
    d.G_FPM = d.mu_M_FPM * d.zeta * c2 * c2 * d.dx_univ / d.J_per_bit_eq;
    return d;
}

