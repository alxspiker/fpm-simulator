#include "lattice.hpp"

Z3Lattice::Z3Lattice(int sx_, int sy_, int sz_, const DerivedConstants& d_, int seed)
        : sx(sx_), sy(sy_), sz(sz_), d(d_), rng(seed)
    {
        int N = sx * sy * sz;
        arena.resize(N);
        torsion_partner.resize(N, -1);

        // Operating-point Omega for initialization
        double e_t_op = 0.75;
        double kappa_op = det_quarter_power(e_t_op);
        double Omega_op = d.Omega_max - (d.Omega_max - d.Omega_min) * kappa_op;

        std::normal_distribution<double> n01(0.0, 1.0);

        for (int x = 0; x < sx; x++) {
            for (int y = 0; y < sy; y++) {
                for (int z = 0; z < sz; z++) {
                    int idx = flat(x, y, z);
                    Daemon& dm = arena[idx];
                    // Initialize psi to uniform normalized carrier
                    for (int ch = 0; ch < 9; ch++) {
                        double p_L = std::max(0.0, std::min(1.0, 0.5 + 0.005 * n01(rng)));
                        double re = det_sqrt(p_L / 5.0);
                        dm.psi[ch] = std::complex<double>(re, 0.0);
                    }
                    // Normalize psi
                    normalize_psi(dm);
                    dm.E = 0.5;
                    dm.b = 0.0;
                    // R = identity * 0.3 + small noise
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            dm.R[i][j] = (i == j ? 0.3 : 0.0) + 0.001 * n01(rng);
                    dm.tau = 0.5;
                    dm.pi_val = std::max(0.0, std::min(1.0, 0.5 + 0.005 * n01(rng)));
                    dm.Omega_prev = Omega_op;
                    dm.local_tick_count = 0;
                    dm.is_active = true;
                    dm.locked_until = 0;
                }
            }
        }
        init_torsion_links();
    }

void Z3Lattice::neighbors6(int idx, int out[6]) const {
int x = idx / (sy * sz);
        int yz = idx % (sy * sz);
        int y = yz / sz;
        int z = yz % sz;
        out[0] = flat(x - 1, y, z);
        out[1] = flat(x + 1, y, z);
        out[2] = flat(x, y - 1, z);
        out[3] = flat(x, y + 1, z);
        out[4] = flat(x, y, z - 1);
        out[5] = flat(x, y, z + 1);
    }

double Z3Lattice::euclidean_dist(int idx1, int idx2) const {
int x1 = idx1 / (sy * sz), yz1 = idx1 % (sy * sz);
        int y1 = yz1 / sz, z1 = yz1 % sz;
        int x2 = idx2 / (sy * sz), yz2 = idx2 % (sy * sz);
        int y2 = yz2 / sz, z2 = yz2 % sz;
        int dx = std::abs(x1 - x2); dx = std::min(dx, sx - dx);
        int dy = std::abs(y1 - y2); dy = std::min(dy, sy - dy);
        int dz = std::abs(z1 - z2); dz = std::min(dz, sz - dz);
        return det_sqrt((double)(dx*dx + dy*dy + dz*dz));
    }

void Z3Lattice::inject_baryonic_cluster(int center_idx, double radius, double B_load, double E_inject) {
for (int i = 0; i < size(); i++) {
            double dist = euclidean_dist(i, center_idx);
            if (dist <= radius) {
                Daemon& dm = arena[i];
                dm.E = std::min(E_inject, d.E_max);
                dm.R[0][1] += B_load;
                dm.R[1][0] -= B_load;
            }
        }
    }

void Z3Lattice::normalize_psi(Daemon& dm) {
double norm = 0.0;
        for (int ch = 0; ch < 9; ch++) norm += std::norm(dm.psi[ch]);
        norm = det_sqrt(norm);
        if (norm <= 0.0) {
            for (int ch = 0; ch < 9; ch++) dm.psi[ch] = std::complex<double>(1.0/3.0, 0.0);
            return;
        }
        for (int ch = 0; ch < 9; ch++) dm.psi[ch] /= norm;
    }

void Z3Lattice::init_torsion_links() {
// Pure-gauge torsion seed
        double scale = 0.015;
        double A_val[3][3] = {{0,0,0},{0,0,-scale},{0,scale,0}};
        double A_val_T[3][3] = {{0,0,0},{0,0,scale},{0,-scale,0}};

        std::vector<bool> visited(size(), false);
        for (int x = 0; x < sx; x++) {
            for (int y = 0; y < sy; y++) {
                for (int z = 0; z < sz; z++) {
                    int idx = flat(x, y, z);
                    if (visited[idx]) continue;
                    // Pair with +x neighbor
                    int partner = flat(x + 1, y, z);
                    if (partner < size() && !visited[partner]) {
                        Daemon& da = arena[idx];
                        Daemon& db = arena[partner];
                        // Symmetric + antisymmetric decomposition
                        for (int i = 0; i < 3; i++)
                            for (int j = 0; j < 3; j++) {
                                double S_ij = 0.5 * (da.R[i][j] + da.R[j][i]);
                                da.R[i][j] = S_ij + A_val[i][j];
                                double S_ij_b = 0.5 * (db.R[i][j] + db.R[j][i]);
                                db.R[i][j] = S_ij_b + A_val_T[i][j];
                            }
                        torsion_links.push_back({idx, partner});
                        torsion_partner[idx] = partner;
                        torsion_partner[partner] = idx;
                        visited[idx] = true;
                        visited[partner] = true;
                    }
                }
            }
        }
    }
