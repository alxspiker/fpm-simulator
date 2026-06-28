#pragma once
#include "viscosity.hpp"

struct LagrangianConfig {
    double w_D = 12.0 / 80.0;
    double w_I = 8.0  / 80.0;
    double w_T = 4.0  / 80.0;
    double w_A = 12.0 / 80.0;
};

double axcore_lagrangian(const Daemon& dm, const DerivedConstants& d, double Omega_new, const LagrangianConfig& cfg, double& C_sem_out, double& C_geo_out, double& smooth_out);

static constexpr uint64_t ACTION_STEPS_PER_UNIVERSAL_TICK = 1000000ULL;

struct SchedulerEvent {
    uint64_t next_step;
    int    flat_idx;
};

class RadixHeap {
public:
    RadixHeap();
    bool empty() const;
    size_t size() const;
    void push(const SchedulerEvent& ev);
    SchedulerEvent pop();
private:
    std::vector<SchedulerEvent> buckets_[65];
    uint64_t last_step_;
    size_t total_size_;
    static int highest_differing_bit(uint64_t value);
    int bucket_index(uint64_t step) const;
};

class ThermodynamicScheduler {
public:
    Z3Lattice& lattice;
    const DerivedConstants& d;
    LagrangianConfig cfg;
    uint64_t universal_step;
    RadixHeap heap;

    ThermodynamicScheduler(Z3Lattice& lat, const DerivedConstants& dc);
    bool step_one(double& L_out, double& O_out, double& k_out);
    void run_batch(int n_universal_ticks, std::vector<double>& t_hist, std::vector<double>& total_E_hist, std::vector<double>& mean_L_hist, std::vector<double>& mean_Omega_hist, std::vector<double>& active_frac_hist);
    void mean_field_truth_target();
private:
    uint64_t period_to_steps(double L_clamped) const;
};
