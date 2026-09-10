/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_EXPM_H
#define LINE_UTIL_EXPM_H

/**
 * Matrix exponential by scaling and squaring with a diagonal Pade approximant.
 *
 * This is the numerical primitive behind MATLAB's expm(), which the kpctoolbox
 * MAP counting-process descriptors (map_cdf, map_pdf, map_acfc, map_count_var,
 * map_varcount, map_count_moment) and the LRU(m)-MAP TTL cache approximation
 * (cache_lrum_map_levelstats) all call. It has no MATLAB source file of its
 * own in the tree; the reference is the algorithm of
 *
 *   N. J. Higham, "The scaling and squaring method for the matrix exponential
 *   revisited", SIAM J. Matrix Anal. Appl. 26(4):1179-1193, 2005,
 *
 * which is what MATLAB implements. Given A, pick a scaling s so that
 * ||2^-s A||_1 <= theta_m, evaluate the [m/m] Pade approximant
 *
 *   r_m(X) = D_m(X)^-1 N_m(X),   N_m(X) = sum_k c_k X^k,
 *                                D_m(X) = sum_k (-1)^k c_k X^k,
 *   c_0 = 1,  c_k = c_{k-1} (m-k+1) / (k (2m-k+1)),
 *
 * and square the result s times: exp(A) = r_m(2^-s A)^(2^s).
 *
 * ARITHMETIC: the answer is an approximation controlled by a tolerance -- the
 * Pade truncation error is nonzero for a general A no matter how the arithmetic
 * is carried out -- so this is gated on num_traits<T>::has_transcendental and
 * refuses to instantiate at exact rational arithmetic. It is nevertheless exact
 * (to rounding) on the two cases that matter for testing: A = 0 returns the
 * identity by a shortcut, and a nilpotent A with A^q = 0, q <= 2m+1, is
 * reproduced exactly because the Pade error term is O(X^(2m+1)).
 *
 * PRECISION AWARENESS: Higham's theta table is calibrated for double unit
 * roundoff, and using it unchanged at Real<50> would silently deliver only
 * double accuracy. The bound on the [m/m] Pade error,
 *
 *   |exp(x) - r_m(x)| ~ C_m |x|^(2m+1),   C_m = (m!)^2 / ((2m)! (2m+1)!),
 *
 * is inverted here against the working precision of T, so theta_m shrinks as
 * the precision grows (theta_13 = 5.1 at double, 0.28 at 50 digits, 0.0039 at
 * 100 digits). At double the resulting thresholds agree with Higham's table to
 * within a few percent.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {

namespace detail {

/** Decimal digits carried by the working type, used to size the scaling. */
template <class T>
inline int expm_digits10() {
    const int d = std::numeric_limits<T>::digits10;
    return d > 0 ? d : 15;
}

/**
 * Largest ||X||_1 for which the [m/m] Pade approximant meets the working
 * precision, from C_m theta^(2m+1) <= 10^-digits10. Capped at Higham's 5.37 so
 * the double path never scales less than the published algorithm.
 */
inline double expm_theta(int m, int digits10) {
    const double lgC = 2.0 * std::lgamma(m + 1.0) - std::lgamma(2.0 * m + 1.0) -
                       std::lgamma(2.0 * m + 2.0);
    const double log10C = lgC / std::log(10.0);
    const double t = std::pow(10.0, (-log10C - digits10) / (2.0 * m + 1.0));
    return t < 5.37 ? t : 5.37;
}

/** ||A||_1, the maximum absolute column sum, as a double. */
template <class T>
double expm_norm1(const Matrix<T>& A) {
    double best = 0.0;
    for (std::size_t j = 0; j < A.cols(); ++j) {
        double s = 0.0;
        for (std::size_t i = 0; i < A.rows(); ++i)
            s += num_traits<T>::to_double(num_abs(T(A(i, j))));
        if (s > best) best = s;
    }
    return best;
}

/** Diagonal [m/m] Pade approximant of exp at the matrix X. */
template <class T>
Matrix<T> expm_pade(const Matrix<T>& X, unsigned m) {
    const std::size_t n = X.rows();
    // Pade coefficients c_k, exact ratios of integers evaluated in T.
    std::vector<T> c(m + 1);
    c[0] = num_traits<T>::from_int(1);
    for (unsigned k = 1; k <= m; ++k) {
        const T num = num_traits<T>::from_int(static_cast<long>(m - k + 1));
        const T den = num_traits<T>::from_int(static_cast<long>(k)) *
                      num_traits<T>::from_int(static_cast<long>(2 * m - k + 1));
        c[k] = c[k - 1] * num / den;
    }

    Matrix<T> N = eye<T>(n);   // c_0 I
    Matrix<T> D = eye<T>(n);   // c_0 I
    Matrix<T> P = eye<T>(n);   // X^k
    for (unsigned k = 1; k <= m; ++k) {
        P = matmul(P, X);
        const bool odd = (k % 2u) == 1u;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                const T term = c[k] * P(i, j);
                N(i, j) += term;
                if (odd)
                    D(i, j) -= term;
                else
                    D(i, j) += term;
            }
    }
    return matmul(inverse(D), N);
}

}  // namespace detail

/**
 * Matrix exponential exp(A).
 *
 * @param A square matrix
 * @return exp(A), by scaling and squaring with a Pade approximant whose degree
 *         and scaling are chosen for the working precision of T
 */
template <class T>
Matrix<T> expm(const Matrix<T>& A) {
    static_assert(num_traits<T>::has_transcendental,
                  "expm is a tolerance-controlled approximation and requires "
                  "transcendental (inexact) arithmetic");
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("expm: matrix is not square");
    if (n == 0) throw InputError("expm: empty matrix");

    const T zero = num_traits<T>::from_int(0);
    bool allzero = true;
    for (std::size_t i = 0; i < n && allzero; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (!(A(i, j) == zero)) {
                allzero = false;
                break;
            }
    if (allzero) return eye<T>(n);  // exp(0) = I, exactly

    const double nrm = detail::expm_norm1(A);
    if (!(nrm == nrm) || nrm == std::numeric_limits<double>::infinity())
        throw NumericError("expm: matrix contains a non-finite entry");

    const int d10 = detail::expm_digits10<T>();
    const unsigned degrees[5] = {3u, 5u, 7u, 9u, 13u};
    unsigned m = 13u;
    int s = 0;
    bool picked = false;
    for (int k = 0; k < 5; ++k) {
        if (nrm <= detail::expm_theta(static_cast<int>(degrees[k]), d10)) {
            m = degrees[k];
            picked = true;
            break;
        }
    }
    if (!picked) {
        const double theta13 = detail::expm_theta(13, d10);
        s = static_cast<int>(std::ceil(std::log2(nrm / theta13)));
        if (s < 0) s = 0;
        if (s > 4096) throw NumericError("expm: matrix norm is too large to scale");
    }

    Matrix<T> X = A;
    if (s > 0) {
        const T half = num_traits<T>::from_rational(1, 2);
        T scale = num_traits<T>::from_int(1);
        for (int k = 0; k < s; ++k) scale *= half;  // 2^-s, exact in binary FP
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) X(i, j) *= scale;
    }

    Matrix<T> E = detail::expm_pade(X, m);
    for (int k = 0; k < s; ++k) E = matmul(E, E);
    return E;
}

/** exp(t A), the form every MAP descriptor actually needs. */
template <class T>
Matrix<T> expm(const Matrix<T>& A, const T& t) {
    Matrix<T> B = A;
    for (std::size_t i = 0; i < B.rows(); ++i)
        for (std::size_t j = 0; j < B.cols(); ++j) B(i, j) *= t;
    return expm(B);
}

}  // namespace line

#endif  // LINE_UTIL_EXPM_H
