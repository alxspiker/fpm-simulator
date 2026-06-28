#pragma once
#include "constants.hpp"

struct Daemon {
    std::complex<double> psi[9];
    double E;
    double b;
    double R[3][3];
    double tau;
    double pi_val;
    double Omega_prev;
    int    local_tick_count;
    bool   is_active;
    uint64_t locked_until;
};

class Z3Lattice {
public:
    int sx, sy, sz;
    std::vector<Daemon> arena;
    const DerivedConstants& d;
    std::mt19937_64 rng;
    std::vector<std::pair<int,int>> torsion_links;
    std::vector<int> torsion_partner;

    Z3Lattice(int sx_, int sy_, int sz_, const DerivedConstants& d_, int seed = 17);

    inline int flat(int x, int y, int z) const {
        return ((x % sx + sx) % sx) * sy * sz + ((y % sy + sy) % sy) * sz + ((z % sz + sz) % sz);
    }
    inline int size() const { return sx * sy * sz; }

    void neighbors6(int idx, int out[6]) const;
    double euclidean_dist(int idx1, int idx2) const;
    void inject_baryonic_cluster(int center_idx, double radius, double B_load, double E_inject);

private:
    void normalize_psi(Daemon& dm);
    void init_torsion_links();
};
