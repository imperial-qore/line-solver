/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_MVN_RECTANGLE_H
#define LINE_SOLVERS_FLUID_FLUID_MVN_RECTANGLE_H

/**
 * Port of `fluid_mvn_rectangle.m`: the rectangle probability
 * P(a <= Y <= b) for Y ~ Normal(m, C), the cell integral behind `getProbAggr`
 * under the moment-closure methods.
 *
 * The integral has no closed form beyond one dimension, so it is evaluated by
 * the separation-of-variables transformation of Genz (1992): the Cholesky factor
 * of C turns the rectangle into an iterated integral over the unit cube whose
 * integrand is a product of normal-CDF differences, and the first coordinate is
 * integrated exactly. The remaining cube is integrated with a DETERMINISTIC
 * Richtmyer lattice rule, frac(k sqrt(p_j)) over the first primes, averaged with
 * its antithetic reflection. Determinism is required here, not merely
 * convenient: the four codebases must return the same number, and a randomized
 * rule would make them agree only in distribution.
 *
 * C may be SINGULAR, which is the common case: a closed population fixes the sum
 * of the station coordinates, so the covariance of a station holding a whole
 * class is rank deficient. A coordinate whose CONDITIONAL variance vanishes is
 * not integrated; it is a hard constraint, contributing 1 when the conditional
 * mean falls inside its interval and 0 otherwise.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {
namespace detail {

/**
 * The first 100 primes, listed rather than sieved so that the MATLAB, Java and
 * Python twins generate the identical lattice.
 */
inline const std::vector<int>& mvn_primes() {
    static const std::vector<int> p = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
        179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
        283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
        419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541};
    return p;
}

/** The standard normal cdf. */
inline double mvn_phi(double x) {
    if (std::isinf(x)) return x > 0.0 ? 1.0 : 0.0;
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

/**
 * The standard normal quantile, MATLAB's `sqrt(2)*erfinv(2u-1)`.
 *
 * Acklam's rational approximation (relative error 1.15e-9) refined by one
 * Halley step on `mvn_phi`, which takes it to double precision. C++ has no
 * `erfinv`, and bisection would be too slow here: this is called once per
 * coordinate per lattice point.
 */
inline double mvn_phi_inv(double u) {
    static const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                                1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                                6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                                -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                                3.754408661907416e+00};
    const double plow = 0.02425, phigh = 1.0 - plow;
    double x;
    if (u < plow) {
        const double q = std::sqrt(-2.0 * std::log(u));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else if (u > phigh) {
        const double q = std::sqrt(-2.0 * std::log(1.0 - u));
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else {
        const double q = u - 0.5, r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }
    // Halley refinement
    const double e = mvn_phi(x) - u;
    const double pdf = std::exp(-0.5 * x * x) / std::sqrt(2.0 * 3.14159265358979323846);
    if (pdf > 0.0) {
        const double t = e / pdf;
        x -= t / (1.0 + 0.5 * x * t);
    }
    return x;
}

/**
 * Cholesky factor of a symmetric positive SEMI-definite matrix. A vanishing
 * pivot leaves a zero row and column, which the caller reads as a deterministic
 * coordinate rather than as a failure.
 */
inline Matrix<double> mvn_chol_psd(const Matrix<double>& C, double dtol) {
    const std::size_t d = C.rows();
    Matrix<double> L(d, d, 0.0);
    for (std::size_t i = 0; i < d; ++i) {
        double v = C(i, i);
        for (std::size_t j = 0; j < i; ++j) v -= L(i, j) * L(i, j);
        if (v > dtol) {
            L(i, i) = std::sqrt(v);
            for (std::size_t r = i + 1; r < d; ++r) {
                double s = C(r, i);
                for (std::size_t j = 0; j < i; ++j) s -= L(r, j) * L(i, j);
                L(r, i) = s / L(i, i);
            }
        } else {
            L(i, i) = 0.0;
            for (std::size_t r = i + 1; r < d; ++r) L(r, i) = 0.0;
        }
    }
    return L;
}

/** One point of the transformed integrand, Genz's recursion over the coordinates. */
inline double mvn_evaluate(const Matrix<double>& L, const std::vector<double>& al,
                           const std::vector<double>& bu, const std::vector<double>& w,
                           std::vector<double>& y, std::size_t last_int, bool has_int, double ctol,
                           bool antithetic) {
    const std::size_t d = al.size();
    double f = 1.0;
    std::size_t kw = 0;
    for (std::size_t i = 0; i < d; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < i; ++j) s += L(i, j) * y[j];
        if (L(i, i) > 0.0) {
            const double dd = mvn_phi((al[i] - s) / L(i, i));
            const double ee = mvn_phi((bu[i] - s) / L(i, i));
            f *= std::max(0.0, ee - dd);
            if (f == 0.0) return 0.0;
            if (!(has_int && i == last_int)) {
                const double wk = antithetic ? 1.0 - w[kw] : w[kw];
                ++kw;
                double u = dd + wk * (ee - dd);
                // the inverse cdf is evaluated strictly inside the unit interval
                u = std::min(std::max(u, 1e-15), 1.0 - 1e-15);
                y[i] = mvn_phi_inv(u);
            }
        } else {
            // zero conditional variance: the coordinate is pinned at s, so the
            // cell is either met or not
            if (s < al[i] - ctol || s > bu[i] + ctol) return 0.0;
            y[i] = 0.0;
        }
    }
    return f;
}

}  // namespace detail

/** Lattice points per antithetic pair. */
const std::size_t FLUID_MVN_POINTS = 4096;

/**
 * P(a <= Y <= b) for Y ~ Normal(m, C).
 *
 * @param m       mean vector
 * @param C       covariance, symmetric positive semi-definite
 * @param a       lower corner, -infinity allowed
 * @param b       upper corner, +infinity allowed
 * @param npoints lattice points per antithetic pair
 */
inline double fluid_mvn_rectangle(const std::vector<double>& m, const Matrix<double>& C,
                                  const std::vector<double>& a, const std::vector<double>& b,
                                  std::size_t npoints = FLUID_MVN_POINTS) {
    const std::size_t d = m.size();
    if (d == 0) return 1.0;
    std::vector<double> al(d), bu(d);
    for (std::size_t i = 0; i < d; ++i) {
        al[i] = a[i] - m[i];
        bu[i] = b[i] - m[i];
        if (bu[i] <= al[i]) return 0.0;
    }

    // scale-relative tolerances: dtol decides which coordinate carries noise,
    // ctol whether a deterministic coordinate satisfies its constraint
    double scale = 1.0;
    for (std::size_t i = 0; i < d; ++i) scale = std::max(scale, std::fabs(C(i, i)));
    const double dtol = 1e-12 * scale;
    const double ctol = 1e-6 * std::sqrt(scale);

    const Matrix<double> L = detail::mvn_chol_psd(C, dtol);
    std::size_t n_int = 0, last_int = 0;
    bool has_int = false;
    for (std::size_t i = 0; i < d; ++i)
        if (L(i, i) > 0.0) {
            ++n_int;
            last_int = i;
            has_int = true;
        }
    const std::size_t nw = (n_int > 0) ? n_int - 1 : 0;
    if (nw > detail::mvn_primes().size())
        throw InputError(
            "fluid_mvn_rectangle: the lattice rule carries generators for at most 100 integration "
            "dimensions. Aggregate classes before evaluating the cell");

    std::vector<double> alpha(nw);
    for (std::size_t j = 0; j < nw; ++j) alpha[j] = std::sqrt(static_cast<double>(detail::mvn_primes()[j]));

    const std::size_t npairs = (nw == 0) ? 1 : npoints;
    const std::size_t n_eval = (nw == 0) ? 1 : 2 * npoints;
    std::vector<double> w(nw), y(d, 0.0);
    double acc = 0.0;
    for (std::size_t k = 1; k <= npairs; ++k) {
        for (std::size_t j = 0; j < nw; ++j) {
            const double v = static_cast<double>(k) * alpha[j];
            w[j] = v - std::floor(v);
        }
        acc += detail::mvn_evaluate(L, al, bu, w, y, last_int, has_int, ctol, false);
        if (nw > 0) acc += detail::mvn_evaluate(L, al, bu, w, y, last_int, has_int, ctol, true);
    }

    const double p = acc / static_cast<double>(n_eval);
    return std::min(std::max(p, 0.0), 1.0);
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_MVN_RECTANGLE_H
