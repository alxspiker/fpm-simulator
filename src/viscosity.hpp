#pragma once
#include "lattice.hpp"

double shear_aggregate(const double R[3][3]);
double trace_curvature(const double R[3][3]);
double mobility(double K1, double S9);
void spectral_gap_weights(const double R[3][3], double& w_H, double& w_S);
double normalized_entropy_H(const double p[2]);
double routing_balance_S(const double p[2]);
void born_probabilities(const Daemon& dm, double probs[9]);
double p_L(const Daemon& dm);
double p_R(const Daemon& dm);
double dispersion(const Daemon& dm);
void viscosity_update(const Daemon& dm, const DerivedConstants& d, double B_load, double& Omega, double& kappa, double& C_N);
