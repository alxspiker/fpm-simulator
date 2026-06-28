#pragma once
#include "fast_math.hpp"

static const double HBAR       = 1.054571817e-34;
static const double H_PLANCK   = 6.62607015e-34;
static const double K_B        = 1.380649e-23;
static const double C_LIGHT    = 2.99792458e8;
static const double M_ELECTRON = 9.1093837015e-31;
static const double G_CODATA   = 6.67430e-11;
static const double T_SUBSTRATE= 300.0;

static const double PLANCK_NS  = 0.965;
static const double PLANCK_RATIO_LCDM = 5.357;
static const double PLANCK_TT_RMS     = 4.06e-5;
static const double CERN_MUON_GAMMA   = 29.3;

struct Axioms {
    int    dim_space  = 3;
    int    dim_causal = 4;
    int    routed_channels_per_axis = 3;
    double axcore_base           = 4.0;
    double axcore_critic_coef    = 12.0;
    double axcore_fitness_coef   = 8.0;
    double axcore_min_floor      = 0.5;
    double axcore_dt_clip_low    = 0.70;
    double axcore_dt_clip_high   = 2.30;
    double axcore_dt_H_coef      = 0.90;
    double axcore_dt_S_coef      = 0.50;
    double axcore_dt_offset      = 0.80;
    double bench_Bdt             = 1.0;
    double bench_fitness         = 0.5;
    double bench_kappa_strat     = 1.0;
    double c_light    = C_LIGHT;
    double m_e        = M_ELECTRON;
    double h_planck   = H_PLANCK;
    double k_B        = K_B;
    double T_substrate= T_SUBSTRATE;

    int n_directed() const { return routed_channels_per_axis * routed_channels_per_axis; }
    int n_trace()    const { return 1; }
};

struct DerivedConstants {
    double alpha     = 0.0;
    double beta      = 0.0;
    double chi_arrow = 0.0;
    double Omega_min = 0.0;
    double Omega_max = 0.0;
    double E_max     = 0.0;
    double e_exp     = 0.0;
    double e_floor   = 0.0;
    double ledger_inertia_ratio = 0.0;
    double c0        = 0.0;
    double lam       = 0.0;
    double L_max     = 0.0;
    double L_rest    = 0.0;
    double gamma_max = 0.0;
    double alpha_PP  = 0.0;
    int64_t N_bit_eq = 0;
    double A_FPM     = 0.0;
    double n_s       = 0.0;
    double r_tensor  = 0.0;
    double ell_A     = 0.0;
    double ell_freeze= 0.0;
    double ell_D     = 0.0;
    double dt_univ   = 0.0;
    double dx_univ   = 0.0;
    double f_univ    = 0.0;
    double zeta      = 0.0;
    double J_per_bit_eq = 0.0;
    double mu_M_FPM  = 0.0;
    double G_FPM     = 0.0;
    double calib     = 0.0;
    int    n_blade   = 2;
    int    n_directed= 9;
    int    n_trace   = 1;
};

double axcore_cost(double H, double S, double f, double kappa_strat, const Axioms& ax);
int64_t compute_N_bit_eq_exact(double alpha_PP);
double iterate_alpha_pp_counterterm(double a1, double a2, double L_rest, double tol = 1e-15, int max_iter = 200);
DerivedConstants derive_all(const Axioms& ax);
