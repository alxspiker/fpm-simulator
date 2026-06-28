#include "viscosity.hpp"

double shear_aggregate(const double R[3][3]) {
    double sum = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            sum += R[i][j] * R[i][j];
    return det_sqrt(sum / 9.0);
}

double trace_curvature(const double R[3][3]) {
    return std::abs(R[0][0] + R[1][1] + R[2][2]);
}

double mobility(double K1, double S9) {
    return det_mobility(K1, S9);
}

void spectral_gap_weights(const double R[3][3], double& w_H, double& w_S) {
    // Simplified SVD via eigenvalues of R^T R (3x3 symmetric)
    double RtR[3][3] = {};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                RtR[i][j] += R[k][i] * R[k][j];

    // Trace, off-diag sums for characteristic polynomial
    double tr = RtR[0][0] + RtR[1][1] + RtR[2][2];
    double det = RtR[0][0]*(RtR[1][1]*RtR[2][2] - RtR[1][2]*RtR[2][1])
               - RtR[0][1]*(RtR[1][0]*RtR[2][2] - RtR[1][2]*RtR[2][0])
               + RtR[0][2]*(RtR[1][0]*RtR[2][1] - RtR[1][1]*RtR[2][0]);

    // Sum of 2x2 minors (half of (tr^2 - tr(RtR^2)))
    double cof = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            cof += RtR[i][j] * RtR[i][j]; // simplified
        }
    cof = 0.5 * (tr * tr - (RtR[0][0]*RtR[0][0] + RtR[0][1]*RtR[0][1] + RtR[0][2]*RtR[0][2]
                            + RtR[1][0]*RtR[1][0] + RtR[1][1]*RtR[1][1] + RtR[1][2]*RtR[1][2]
                            + RtR[2][0]*RtR[2][0] + RtR[2][1]*RtR[2][1] + RtR[2][2]*RtR[2][2]));

    // Eigenvalues via Cardano (cubic: x^3 - tr*x^2 + cof*x - det = 0)
    // Use Newton's method from tr/3 as starting point
    double sigmas[3];
    double p = cof - tr*tr/3.0;
    double q = -2.0*tr*tr*tr/27.0 + tr*cof/3.0 - det;
    double disc = q*q/4.0 + p*p*p/27.0;
    if (disc > 0) {
        double sq = det_sqrt(disc);
        double u = det_cbrt(-q/2.0 + sq);
        double v = det_cbrt(-q/2.0 - sq);
        sigmas[0] = u + v + tr/3.0;
        // For the remaining eigenvalues, use simplified approximation
        sigmas[1] = std::max(1e-30, (tr - sigmas[0]) / 2.0);
        sigmas[2] = std::max(1e-30, tr - sigmas[0] - sigmas[1]);
    } else {
        // All real roots
        double m = 2.0 * det_sqrt(-p/3.0);
        double acos_arg = (p != 0.0 && m != 0.0) ? 3.0*q/(p*m) : 1.0;
        acos_arg = std::max(-1.0, std::min(1.0, acos_arg));
        double theta = det_acos(acos_arg) / 3.0;
        sigmas[0] = m * det_cos(theta) + tr/3.0;
        sigmas[1] = m * det_cos(theta - 2.0*PI/3.0) + tr/3.0;
        sigmas[2] = m * det_cos(theta - 4.0*PI/3.0) + tr/3.0;
    }
    // Sort descending
    std::sort(sigmas, sigmas + 3, std::greater<double>());
    for (int i = 0; i < 3; i++) sigmas[i] = std::max(1e-30, sigmas[i]);

    double sum_s = sigmas[0] + sigmas[1] + sigmas[2];
    w_H = (sum_s > 0) ? sigmas[0] / sum_s : 1.0/3.0;
    w_S = 1.0 - w_H;
}

double normalized_entropy_H(const double p[2]) {
    double pp[2] = {std::max(p[0], 1e-12), std::max(p[1], 1e-12)};
    double s = pp[0] + pp[1];
    pp[0] /= s; pp[1] /= s;
    double H = -(pp[0] * det_log2(pp[0]) + pp[1] * det_log2(pp[1]));
    return H;
}

double routing_balance_S(const double p[2]) {
    return 1.0 - std::abs(p[0] - p[1]);
}

// Born probabilities from psi
void born_probabilities(const Daemon& dm, double probs[9]) {
    double sum = 0.0;
    for (int ch = 0; ch < 9; ch++) {
        probs[ch] = std::norm(dm.psi[ch]);
        sum += probs[ch];
    }
    if (sum <= 0.0) {
        for (int ch = 0; ch < 9; ch++) probs[ch] = 1.0/9.0;
        return;
    }
    for (int ch = 0; ch < 9; ch++) probs[ch] /= sum;
}

double p_L(const Daemon& dm) {
    double probs[9]; born_probabilities(dm, probs);
    double s = 0; for (int ch = 0; ch < 5; ch++) s += probs[ch]; return s;
}

double p_R(const Daemon& dm) {
    double probs[9]; born_probabilities(dm, probs);
    double s = 0; for (int ch = 5; ch < 9; ch++) s += probs[ch]; return s;
}

double dispersion(const Daemon& dm) {
    // c = sum of psi[:4] . conj(psi[5:9]) / 4
    std::complex<double> c(0,0);
    for (int i = 0; i < 4; i++) c += dm.psi[i] * std::conj(dm.psi[5+i]);
    c /= 4.0;
    return 2.0 * det_complex_abs(c);
}

void viscosity_update(const Daemon& dm, const DerivedConstants& d,
                      double B_load, double& Omega, double& kappa, double& C_N) {
    double p[2] = {p_L(dm), p_R(dm)};
    double H = normalized_entropy_H(p);
    double S = routing_balance_S(p);
    double w_H, w_S;
    spectral_gap_weights(dm.R, w_H, w_S);
    double A_N = w_H * H + w_S * S;
    C_N = std::min(A_N, 1.0);
    double e_B = std::max(1.0 / det_three_quarter_power(1.0 + B_load), d.e_floor);
    double e_t = (dm.E / d.E_max) * e_B;
    e_t = std::max(0.0, std::min(1.0, e_t));
    double g_e = det_quarter_power(e_t);
    kappa = C_N * g_e;
    kappa = std::max(0.0, std::min(1.0, kappa));
    double dOmega = d.Omega_max - d.Omega_min;
    Omega = d.Omega_max - dOmega * kappa;
}

