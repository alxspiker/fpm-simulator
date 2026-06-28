#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include "probes.hpp"

struct SparcInjectionPoint {
    double r_kpc = 0.0;
    double v_obs = 0.0;
    double b_load = 0.0;
};

struct SparcInjectionPayload {
    std::vector<SparcInjectionPoint> points;
};

struct CalibrationResult {
    double rel_err_G;
    double rel_err_ns;
    bool ell_D_in_range;
    bool gamma_above_cern;
    int64_t N_bit_eq_exact;
    double N_bit_eq_continuous;
};

CalibrationResult calibrate(const DerivedConstants& d);
bool read_file_string(const std::string& path, std::string& out);
bool parse_json_number_after_key(const std::string& text, size_t start, const std::string& key, double& value, size_t& value_end);
bool parse_sparc_injection_payload(const std::string& path, SparcInjectionPayload& payload);
void inject_core_load(Z3Lattice& lattice, int center, double B_load, double E_inject);
double measure_shell_L(Z3Lattice& lattice, ThermodynamicScheduler& sched, const DerivedConstants& d, int center, double shell_r);

int run_sparc_payload_mode(const std::string& payload_path, const std::string& output_path);
int run_dynamic_torsion_mode(const std::string& output_path);
int run_photon_propagation_mode(const std::string& output_path);

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
                int events_processed);
