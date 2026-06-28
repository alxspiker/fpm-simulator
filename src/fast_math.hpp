#pragma once
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <complex>
#include <cstdint>
#include <limits>
#include <functional>
#include <random>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

static const double PI = 3.14159265358979323846;

// DETERMINISTIC FAST-MATH KERNEL
// =============================================================================

inline double det_abs(double x) {
    return x < 0.0 ? -x : x;
}

inline double det_rsqrt(double number) {
    if (number <= 0.0) return 0.0;

    double y = number;
    double x2 = y * 0.5;
    uint64_t bits = 0;
    std::memcpy(&bits, &y, sizeof(bits));
    bits = 0x5fe6eb50c7b537a9ULL - (bits >> 1);
    std::memcpy(&y, &bits, sizeof(y));

    // Four fixed Newton steps tighten the 64-bit seed without libc calls.
    y = y * (1.5 - (x2 * y * y));
    y = y * (1.5 - (x2 * y * y));
    y = y * (1.5 - (x2 * y * y));
    y = y * (1.5 - (x2 * y * y));
    return y;
}

inline double det_sqrt(double number) {
    return number <= 0.0 ? 0.0 : number * det_rsqrt(number);
}

inline double det_pow2_int(int exp) {
    if (exp < -1022) return 0.0;
    if (exp > 1023) exp = 1023;

    uint64_t bits = (uint64_t)(exp + 1023) << 52;
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

inline double det_log2(double x) {
    if (x <= 0.0) return -1024.0;

    uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof(bits));
    int exp = (int)((bits >> 52) & 0x7ff) - 1023;

    // The hot-path probabilities are normal doubles. Keep subnormal handling
    // deterministic anyway by scaling into the normal range before extraction.
    if (exp == -1023) {
        x *= det_pow2_int(54);
        std::memcpy(&bits, &x, sizeof(bits));
        exp = (int)((bits >> 52) & 0x7ff) - 1023 - 54;
    }

    uint64_t mant_bits = (bits & 0x000fffffffffffffULL) | (1023ULL << 52);
    double m = 0.0;
    std::memcpy(&m, &mant_bits, sizeof(m));

    static const double INV_LN2 = 1.44269504088896340736;
    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;
    double term = z;
    double sum = term;
    for (int denom = 3; denom <= 31; denom += 2) {
        term *= z2;
        sum += term / (double)denom;
    }

    return (double)exp + (2.0 * sum * INV_LN2);
}

inline double det_root5(double x) {
    if (x <= 0.0) return 0.0;
    double y = 1.0;
    for (int i = 0; i < 10; i++) {
        double y2 = y * y;
        double y4 = y2 * y2;
        y = 0.2 * (4.0 * y + x / y4);
    }
    return y;
}

inline double det_root9(double x) {
    if (x <= 0.0) return 0.0;
    double y = x > 1.0 ? 1.0 + (x - 1.0) / 9.0 : 1.0;
    for (int i = 0; i < 12; i++) {
        double y2 = y * y;
        double y4 = y2 * y2;
        double y8 = y4 * y4;
        y = (8.0 * y + x / y8) / 9.0;
    }
    return y;
}

inline double det_pow_five_ninths(double x) {
    double root9 = det_root9(x);
    double root9_2 = root9 * root9;
    double root9_4 = root9_2 * root9_2;
    return root9_4 * root9;
}

inline double det_quarter_power(double x) {
    return x <= 0.0 ? 0.0 : det_sqrt(det_sqrt(x));
}

inline double det_three_quarter_power(double x) {
    return x <= 0.0 ? 0.0 : det_sqrt(det_sqrt(x * x * x));
}

inline double det_cbrt(double x) {
    if (x == 0.0) return 0.0;

    double ax = det_abs(x);
    uint64_t bits = 0;
    std::memcpy(&bits, &ax, sizeof(bits));
    int exp = (int)((bits >> 52) & 0x7ff) - 1023;

    if (exp == -1023) {
        ax *= det_pow2_int(54);
        std::memcpy(&bits, &ax, sizeof(bits));
        exp = (int)((bits >> 52) & 0x7ff) - 1023 - 54;
    }

    uint64_t mant_bits = (bits & 0x000fffffffffffffULL) | (1023ULL << 52);
    double m = 0.0;
    std::memcpy(&m, &mant_bits, sizeof(m));

    int q = exp / 3;
    int rem = exp - 3 * q;
    if (rem < 0) {
        q--;
        rem += 3;
    }

    double scaled = m * det_pow2_int(rem);
    double y = 1.0 + (scaled - 1.0) / 3.0;
    for (int i = 0; i < 8; i++) {
        y = (2.0 * y + scaled / (y * y)) / 3.0;
    }

    double result = y * det_pow2_int(q);
    return x < 0.0 ? -result : result;
}

inline double det_mobility(double K1, double S9) {
    double num = det_root5(1.0 + K1);
    double den_base = det_root5(1.0 + S9);
    double den2 = den_base * den_base;
    double den4 = den2 * den2;
    double den8 = den4 * den4;
    double den9 = den8 * den_base;
    return den9 > 0.0 ? num / den9 : 0.0;
}

inline double det_wrap_half_pi(double x, bool& cos_flip) {
    static const double TWO_PI = 6.28318530717958647692;
    static const double HALF_PI = 1.57079632679489661923;

    int64_t turns = (int64_t)(x / TWO_PI + (x >= 0.0 ? 0.5 : -0.5));
    x -= (double)turns * TWO_PI;
    cos_flip = false;

    if (x > HALF_PI) {
        x = PI - x;
        cos_flip = true;
    } else if (x < -HALF_PI) {
        x = -PI - x;
        cos_flip = true;
    }
    return x;
}

inline double det_wrap_pi(double x) {
    static const double TWO_PI = 6.28318530717958647692;
    int64_t turns = (int64_t)(x / TWO_PI + (x >= 0.0 ? 0.5 : -0.5));
    return x - (double)turns * TWO_PI;
}

inline double det_sin(double x) {
    bool cos_flip = false;
    x = det_wrap_half_pi(x, cos_flip);
    double x2 = x * x;
    double x3 = x2 * x;
    double x5 = x3 * x2;
    double x7 = x5 * x2;
    double x9 = x7 * x2;
    return x - (x3 * 0.16666666666666666667)
             + (x5 * 0.00833333333333333333)
             - (x7 * 0.00019841269841269841)
             + (x9 * 0.00000275573192239859);
}

inline double det_cos(double x) {
    bool flip = false;
    x = det_wrap_half_pi(x, flip);
    double x2 = x * x;
    double x4 = x2 * x2;
    double x6 = x4 * x2;
    double x8 = x6 * x2;
    double y = 1.0 - (x2 * 0.5)
             + (x4 * 0.04166666666666666667)
             - (x6 * 0.00138888888888888889)
             + (x8 * 0.00002480158730158730);
    return flip ? -y : y;
}

inline double det_acos(double x) {
    if (x <= -1.0) return PI;
    if (x >= 1.0) return 0.0;

    double ax = det_abs(x);
    double seed = (((-0.0187293 * ax + 0.0742610) * ax - 0.2121144) * ax + 1.5707288)
                * det_sqrt(1.0 - ax);
    double theta = x < 0.0 ? PI - seed : seed;

    for (int i = 0; i < 4; i++) {
        double s = det_sin(theta);
        if (det_abs(s) <= 1e-14) break;
        theta += (det_cos(theta) - x) / s;
        if (theta < 0.0) theta = 0.0;
        if (theta > PI) theta = PI;
    }
    return theta;
}

inline double det_complex_abs(const std::complex<double>& z) {
    return det_sqrt(std::norm(z));
}

