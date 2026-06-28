#include "scheduler.hpp"

double axcore_lagrangian(const Daemon& dm, const DerivedConstants& d,
                         double Omega_new, const LagrangianConfig& cfg,
                         double& C_sem_out, double& C_geo_out, double& smooth_out) {
    double p[2] = {p_L(dm), p_R(dm)};
    double H = normalized_entropy_H(p);
    double S = routing_balance_S(p);
    double Bdt = 0.80 + 0.90 * H + 0.50 * S;
    Bdt = std::max(0.70, std::min(2.30, Bdt));
    double f = 1.0 - std::abs(p[0] - dm.tau);
    f = std::max(0.0, std::min(1.0, f));
    C_sem_out = d.c0 + cfg.w_D * Bdt + cfg.w_I * (1.0 - f);

    double gap = std::abs(p[0] - dm.tau);
    C_geo_out = cfg.w_T * gap + cfg.w_A * dm.b * std::abs(dm.pi_val - dm.tau);

    double dOmega = std::abs(Omega_new - dm.Omega_prev);
    smooth_out = d.lam * dOmega;

    double L = C_sem_out + C_geo_out + smooth_out;
    L = std::max(d.c0, std::min(d.L_max, L));
    return L;
}

// Thermodynamic Scheduler: radix heap, daemons buy quantized ticks with FPM currency
RadixHeap::RadixHeap() : last_step_(0), total_size_(0) {}
bool RadixHeap::empty() const {
return total_size_ == 0;
    }

size_t RadixHeap::size() const {
return total_size_;
    }

void RadixHeap::push(const SchedulerEvent& ev) {
assert(ev.next_step >= last_step_ && "Causal violation: time moved backward");
        buckets_[bucket_index(ev.next_step)].push_back(ev);
        total_size_++;
    }

SchedulerEvent RadixHeap::pop() {
assert(total_size_ > 0 && "Pop from empty scheduler");

        if (buckets_[0].empty()) {
            int source_bucket = 1;
            while (source_bucket < 65 && buckets_[source_bucket].empty()) {
                source_bucket++;
            }
            assert(source_bucket < 65 && "Radix heap size counter drifted");

            uint64_t min_step = std::numeric_limits<uint64_t>::max();
            int min_flat_idx = std::numeric_limits<int>::max();
            for (const auto& ev : buckets_[source_bucket]) {
                if (ev.next_step < min_step ||
                    (ev.next_step == min_step && ev.flat_idx < min_flat_idx)) {
                    min_step = ev.next_step;
                    min_flat_idx = ev.flat_idx;
                }
            }

            last_step_ = min_step;

            std::vector<SchedulerEvent> cracked = std::move(buckets_[source_bucket]);
            buckets_[source_bucket].clear();
            for (const auto& ev : cracked) {
                buckets_[bucket_index(ev.next_step)].push_back(ev);
            }
        }

        auto best_it = buckets_[0].begin();
        for (auto it = buckets_[0].begin() + 1; it != buckets_[0].end(); ++it) {
            if (it->flat_idx < best_it->flat_idx) {
                best_it = it;
            }
        }

        SchedulerEvent result = *best_it;
        *best_it = buckets_[0].back();
        buckets_[0].pop_back();
        total_size_--;
        return result;
    }

int RadixHeap::highest_differing_bit(uint64_t value) {
#if defined(_MSC_VER)
        unsigned long index = 0;
        _BitScanReverse64(&index, value);
        return (int)index + 1;
#else
        return 64 - __builtin_clzll(value);
#endif
    }

int RadixHeap::bucket_index(uint64_t step) const {
uint64_t diff = step ^ last_step_;
        return diff == 0 ? 0 : highest_differing_bit(diff);
    }

ThermodynamicScheduler::ThermodynamicScheduler(Z3Lattice& lat, const DerivedConstants& dc)
    : lattice(lat), d(dc), universal_step(0)
{
    for (int i = 0; i < lattice.size(); i++) {
        heap.push({0, i});
    }
}

bool ThermodynamicScheduler::step_one(double& L_out, double& O_out, double& k_out) {
if (heap.empty()) return false;
        SchedulerEvent ev = heap.pop();
        Daemon& dm = lattice.arena[ev.flat_idx];
        universal_step = std::max(universal_step, ev.next_step);

        if (!dm.is_active) return false;

        // Compute viscosity and Lagrangian
        double C_N_local;
        viscosity_update(dm, d, 0.0, O_out, k_out, C_N_local);
        double C_sem, C_geo, smooth;
        L_out = axcore_lagrangian(dm, d, O_out, cfg, C_sem, C_geo, smooth);

        // Energy check
        if (dm.E <= 0.0) {
            dm.is_active = false;
            return false;
        }

        // ---------------------------------------------------------
        // HARDWARE LOCK INTERCEPTION
        // ---------------------------------------------------------
        if (universal_step < dm.locked_until) {
            // A phantom thermal event popped early. Defer it until the lock expires.
            heap.push({dm.locked_until, ev.flat_idx});
            return false;
        }

        // ---------------------------------------------------------
        // ACTIVE BOUNDARY RESOLVER: Neighbor Spillover & Exhaust
        // ---------------------------------------------------------
        
        double prior_E = dm.E; 
        
        // 1. Resolve spatial flux weights from the routing tensor
        double w[6];
        double base_iso = 0.05; 
        w[0] = std::max(base_iso, -dm.R[0][0]); 
        w[1] = std::max(base_iso,  dm.R[0][0]); 
        w[2] = std::max(base_iso, -dm.R[1][1]); 
        w[3] = std::max(base_iso,  dm.R[1][1]); 
        w[4] = std::max(base_iso, -dm.R[2][2]); 
        w[5] = std::max(base_iso,  dm.R[2][2]); 

        double max_w = w[0];
        int max_n = 0;
        for (int n = 1; n < 6; n++) {
            if (w[n] > max_w) { max_w = w[n]; max_n = n; }
        }

        // THE GAUGE CONSTRAINT (Evaluate BEFORE paying exhaust)
        // Set threshold to 0.995 to distinguish from the Torsion Link which operates at 0.99
        bool is_causal_wave = (prior_E >= 0.995 * d.E_max) && (max_w >= 4.0); 
        
        // 2. Pay the absolute thermodynamic cost 
        // Photons are pure momentum carriers; they do not pay thermal friction.
        double intended_exhaust = is_causal_wave ? 0.0 : L_out;
        double actual_exhaust = std::min(dm.E, intended_exhaust);
        dm.E -= actual_exhaust;

        if (is_causal_wave) {
            for (int n = 0; n < 6; n++) {
                w[n] = (n == max_n) ? 1.0 : 0.0;
            }
        }

        double sum_w = w[0] + w[1] + w[2] + w[3] + w[4] + w[5];
        if (sum_w <= 0.0) sum_w = 1.0; 
        
        int nbrs[6];
        lattice.neighbors6(ev.flat_idx, nbrs);

        // Bulk energy payload for causal waves
        double transferable_E = is_causal_wave ? dm.E : 0.0;

        // 3. Advect Energy and Momentum across the Z^3 boundary
        for (int n = 0; n < 6; n++) {
            if (w[n] <= 0.0) continue; 
            
            double flux_frac = w[n] / sum_w;
            Daemon& target = lattice.arena[nbrs[n]];
            
            // Advect Thermal Exhaust 
            if (actual_exhaust > 0.0) {
                target.E = std::min(d.E_max, target.E + (actual_exhaust * flux_frac));
                if (target.E > 0.0) target.is_active = true;
            }
            
            // Advect Bulk Energy (Soliton Payload)
            if (transferable_E > 0.0) {
                double e_shift = transferable_E * flux_frac;
                target.E = std::min(d.E_max, target.E + e_shift);
                dm.E = std::max(0.0, dm.E - e_shift);
                if (target.E > 0.0) target.is_active = true;
            }
            
            // Momentum Advection 
            double conductivity = is_causal_wave ? 1.0 : (prior_E / d.E_max);
            double transfer_rate = conductivity * flux_frac; 
            
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double delta = dm.R[i][j] * transfer_rate;
                    target.R[i][j] += delta;
                    dm.R[i][j] -= delta;
                }
            }
        }
        
        // ---------------------------------------------------------
        // CAUSAL SCHEDULING DELAY (The Speed of Light Limit)
        // ---------------------------------------------------------
        if (is_causal_wave) {
            for (int n = 0; n < 6; n++) {
                if (w[n] <= 0.0) continue;
                int target_idx = nbrs[n];
                
                // Hardware lock the target node perfectly until the transit time passes
                uint64_t wake_time = universal_step + ACTION_STEPS_PER_UNIVERSAL_TICK;
                lattice.arena[target_idx].locked_until = wake_time;
                
                // Push the execution event so it wakes up
                heap.push({wake_time, target_idx});
            }
        }
        // ---------------------------------------------------------

        // Phase rotation from route cost
        double theta = 0.37;
        for (int ch = 0; ch < 9; ch++) {
            double rc_weight = std::abs(dm.R[ch/3][ch%3]);
            if (rc_weight <= 0) rc_weight = 1.0;
            double L_ch = L_out * rc_weight;
            double angle = theta * L_ch;
            dm.psi[ch] *= std::complex<double>(det_cos(angle), -det_sin(angle));
        }
        // Renormalize psi
        double norm = 0;
        for (int ch = 0; ch < 9; ch++) norm += std::norm(dm.psi[ch]);
        norm = det_sqrt(norm);
        if (norm > 0) for (int ch = 0; ch < 9; ch++) dm.psi[ch] /= norm;

        dm.Omega_prev = O_out;

        // Probability drift toward tau
        double p_l = p_L(dm);
        double new_p_L = p_l + 0.001 * (dm.tau - p_l);
        new_p_L = std::max(0.0, std::min(1.0, new_p_L));
        // Rebuild psi from binary probabilities (preserve phases)
        std::complex<double> phases[9];
        for (int ch = 0; ch < 9; ch++) {
            double amp = det_complex_abs(dm.psi[ch]);
            phases[ch] = (amp > 0) ? dm.psi[ch] / amp : std::complex<double>(1,0);
        }
        double p_r_val = 1.0 - new_p_L;
        for (int ch = 0; ch < 5; ch++) dm.psi[ch] = phases[ch] * det_sqrt(new_p_L / 5.0);
        for (int ch = 5; ch < 9; ch++) dm.psi[ch] = phases[ch] * det_sqrt(p_r_val / 4.0);
        // Renormalize
        norm = 0;
        for (int ch = 0; ch < 9; ch++) norm += std::norm(dm.psi[ch]);
        norm = det_sqrt(norm);
        if (norm > 0) for (int ch = 0; ch < 9; ch++) dm.psi[ch] /= norm;

        // R evolution
        double K1 = trace_curvature(dm.R);
        double S9 = shear_aggregate(dm.R);
        double phi = mobility(K1, S9);
        std::normal_distribution<double> n01(0.0, 1.0);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                dm.R[i][j] += 0.0005 * phi * n01(lattice.rng);

        // Metabolic mode
        double e_frac = dm.E / d.E_max;
        if (e_frac <= 0.20) {
            dm.b = std::min(dm.b + 0.003, 1.0);
            if (dm.E <= 0.0) dm.is_active = false;
        } else if (e_frac <= 0.28) {
            dm.b = std::min(dm.b + 0.001, 1.0);
        }

        dm.local_tick_count++;
        if (dm.E <= 0.0) dm.is_active = false;

        // Re-schedule on the integer action lattice. Route cost still determines
        // the period, but the global event ledger advances only in fixed quanta.
        if (!is_causal_wave) {
            double L_clamped = std::max(d.L_rest, std::min(d.L_max, L_out));
            uint64_t period_steps = period_to_steps(L_clamped);
            heap.push({ev.next_step + period_steps, ev.flat_idx});
        }

        return true;
    }

    // Run batch: process events until the quantized universal ledger reaches target

void ThermodynamicScheduler::run_batch(int n_universal_ticks, std::vector<double>& t_hist, std::vector<double>& total_E_hist, std::vector<double>& mean_L_hist, std::vector<double>& mean_Omega_hist, std::vector<double>& active_frac_hist) {
    uint64_t target_step = (uint64_t)n_universal_ticks * ACTION_STEPS_PER_UNIVERSAL_TICK;
    int events = 0;
    int sample_every = lattice.size();

    while (universal_step < target_step && !heap.empty()) {
        double L, O, k;
        step_one(L, O, k);
        events++;

            // Periodic mean-field update
            if (events % (lattice.size() * 10) == 0) {
                mean_field_truth_target();
            }

            // Record trajectory
            if (events % sample_every == 0) {
                double total_E = 0, sum_L = 0, sum_O = 0;
                int active = 0;
                for (int i = 0; i < lattice.size(); i++) {
                    Daemon& dm = lattice.arena[i];
                    total_E += dm.E;
                    if (dm.is_active) {
                        double Oi, ki, C_Ni;
                        viscosity_update(dm, d, 0.0, Oi, ki, C_Ni);
                        double C_sem, C_geo, sm;
                        sum_L += axcore_lagrangian(dm, d, Oi, cfg, C_sem, C_geo, sm);
                        sum_O += Oi;
                        active++;
                    }
                }
                int n = lattice.size();
                t_hist.push_back((double)universal_step / (double)ACTION_STEPS_PER_UNIVERSAL_TICK);
                total_E_hist.push_back(total_E);
                mean_L_hist.push_back(active > 0 ? sum_L / active : 0);
                mean_Omega_hist.push_back(active > 0 ? sum_O / active : 0);
                active_frac_hist.push_back((double)active / n);
            }
        }
    }

    // Mean-field truth target: Z^3 6-face neighbor average

void ThermodynamicScheduler::mean_field_truth_target() {
double eta_flux = 0.5, eta_geo = 0.5;
        std::vector<double> new_taus(lattice.size());
        double sum_pi = 0, sum_tau = 0;

        for (int i = 0; i < lattice.size(); i++) {
            Daemon& dm = lattice.arena[i];
            int nbrs[6];
            lattice.neighbors6(i, nbrs);
            double w_sum = 0, tau_sum = 0;
            for (int n = 0; n < 6; n++) {
                Daemon& nbr = lattice.arena[nbrs[n]];
                double std_R = 0;
                double mean_R = 0;
                for (int ii = 0; ii < 3; ii++)
                    for (int jj = 0; jj < 3; jj++) mean_R += dm.R[ii][jj];
                mean_R /= 9.0;
                for (int ii = 0; ii < 3; ii++)
                    for (int jj = 0; jj < 3; jj++) std_R += (dm.R[ii][jj] - mean_R) * (dm.R[ii][jj] - mean_R);
                std_R = det_sqrt(std_R / 9.0);
                double w = eta_flux * std_R + eta_geo * std::abs(dm.pi_val - nbr.pi_val);
                w = std::max(1e-9, w);
                w_sum += w;
                tau_sum += w * nbr.pi_val;
            }
            new_taus[i] = (w_sum > 0) ? tau_sum / w_sum : dm.pi_val;
            sum_pi += dm.pi_val;
            sum_tau += new_taus[i];
        }
        // Project to enforce sum(tau) = sum(pi)
        double correction = (sum_pi - sum_tau) / lattice.size();
        for (int i = 0; i < lattice.size(); i++) {
            lattice.arena[i].tau = std::max(0.0, std::min(1.0, new_taus[i] + correction));
        }
    }

uint64_t ThermodynamicScheduler::period_to_steps(double L_clamped) const {
double universal_ticks = L_clamped / d.L_rest;
        double quantized = std::round(universal_ticks * (double)ACTION_STEPS_PER_UNIVERSAL_TICK);
        return std::max<uint64_t>(1, (uint64_t)quantized);
    }

