#pragma once
#include <string>
#include "scheduler.hpp"

struct EmergentProbeBorn {
    const DerivedConstants& d;
    EmergentProbeBorn(const DerivedConstants& dc) : d(dc) {}
    void run(int n_samples, double& mean_tv, double& max_tv, std::string& verdict);
};

struct EmergentProbeBell {
    const DerivedConstants& d;
    EmergentProbeBell(const DerivedConstants& dc) : d(dc) {}
    double geometric_torsion_correlation(double a, double b) const;
    double local_torsion_correlation(double a, double b) const;
    double qm_correlation(double a, double b) const;
    double joint_torsion_lrm_correlation(double a, double b) const;
    template<typename Fn> double chsh_value(Fn corr_fn) const {
        double a = 0, ap = PI/2, b = PI/4, bp = -PI/4;
        return std::abs(corr_fn(a,b) + corr_fn(a,bp) + corr_fn(ap,b) - corr_fn(ap,bp));
    }
    void run(double& S_local, double& S_qm, double& S_torsion, double& S_joint, std::string& verdict);
};

struct EmergentProbeFineStructure {
    const DerivedConstants& d;
    EmergentProbeFineStructure(const DerivedConstants& dc) : d(dc) {}
    void run(double& one_over_alpha_bare, double& rel_diff, std::string& verdict);
};

struct EmergentProbeGravity {
    void run(Z3Lattice& lattice, ThermodynamicScheduler& sched, const DerivedConstants& d,
             std::vector<double>& r_out, std::vector<double>& v_emergent_out,
             std::vector<double>& L_mean_out);
};

struct EmergentProbeTimeDilation {
    void run(Z3Lattice& lattice, ThermodynamicScheduler& sched, const DerivedConstants& d,
             double& gamma_near, double& gamma_far, double& tick_ratio,
             int& ticks_near, int& ticks_far);
};
