/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_AOI_DIST2PH_H
#define LINE_API_AOI_AOI_DIST2PH_H

/**
 * MAP to PH conversion for the age-of-information solvers.
 *
 * Templated port of matlab/src/api/aoi/aoi_dist2ph.m. The (D0, D1) pair of a
 * MAP is turned into the (alpha, T) pair the aoi_* algorithms consume:
 *
 *   T     = D0, the sub-generator, unchanged
 *   theta solves theta (D0 + D1) = 0, sum theta = 1
 *   alpha = theta .* (D1 e) normalized, the phase distribution just after a
 *           completion, weighted by where completions actually occur
 *
 * alpha is NOT map_pie. map_pie is theta D1 / (theta D1 e), the embedded
 * arrival chain's stationary vector, which redistributes the mass through the
 * columns of D1; the reference instead keeps the mass in the phase it left
 * from, alpha_i proportional to theta_i (D1 e)_i. For a renewal MAP with
 * D1 = d pi (rank one) the two coincide; for a general MAP they do not, and
 * this port reproduces the reference rather than substituting map_pie.
 *
 * ARITHMETIC. One linear solve, one normalization and sign checks, so this is
 * exact at T = Rational and is instantiated there. MATLAB obtains theta from
 * the overdetermined [Q'; e'] \ [0; 1], which its backslash resolves by QR;
 * on a consistent system that is the exact solution, and the port reaches the
 * same solution through the normal equations, which is a square solve the
 * templated LU can carry exactly.
 *
 * The generator repair of the reference is kept: when the rows of D0 + D1 do
 * not sum to zero to within 1e-10, the diagonal is corrected so that they do.
 * That test is a comparison against a constant and needs no transcendental
 * function, so it carries over verbatim to the exact instantiation, where it
 * is essentially never triggered.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace aoi {

/** The (alpha, T) PH pair produced by aoi_dist2ph. */
template <class T>
struct AoiPh {
    std::vector<T> alpha;  ///< initial probability vector, sums to one
    Matrix<T> Tmat;        ///< sub-generator, MATLAB's T
};

/**
 * Convert a MAP (D0, D1) into the PH pair (alpha, T).
 *
 * @param D0 hidden generator
 * @param D1 arrival matrix
 */
template <class T>
AoiPh<T> aoi_dist2ph(const Matrix<T>& D0, const Matrix<T>& D1) {
    const std::size_t n = D0.rows();
    if (n == 0) throw InputError("aoi_dist2ph: empty process");
    if (D0.cols() != n || D1.rows() != n || D1.cols() != n)
        throw InputError("aoi_dist2ph: D0 and D1 must be square matrices of the same size");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T genTol = num_traits<T>::from_double(1e-10);

    AoiPh<T> out;
    out.Tmat = D0;

    Matrix<T> Q(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Q(i, j) = D0(i, j) + D1(i, j);

    // Repair the generator when the row sums do not vanish (reference behaviour).
    T worst = zero;
    std::vector<T> rowSum(n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) rowSum[i] += Q(i, j);
        const T a = rowSum[i] < zero ? T(-rowSum[i]) : rowSum[i];
        if (a > worst) worst = a;
    }
    if (worst > genTol)
        for (std::size_t i = 0; i < n; ++i) Q(i, i) = Q(i, i) - rowSum[i];

    // theta Q = 0, sum theta = 1, through the normal equations of
    // A theta' = b with A = [Q'; e'] and b = [0; 1].
    Matrix<T> AtA(n, n, zero);
    std::vector<T> Atb(n, one);  // (A' b)_i = 1, the last row of A being e'
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            T s = one;  // contribution of the normalization row
            for (std::size_t k = 0; k < n; ++k) s += Q(i, k) * Q(j, k);
            AtA(i, j) = s;
        }
    std::vector<T> theta = solve(AtA, Atb);

    T mass = zero;
    for (std::size_t i = 0; i < n; ++i) {
        if (theta[i] < zero) theta[i] = zero;
        mass += theta[i];
    }
    if (mass == zero) throw NumericError("aoi_dist2ph: the phase process has no stationary law");
    for (std::size_t i = 0; i < n; ++i) theta[i] = theta[i] / mass;

    std::vector<T> completion(n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) completion[i] += D1(i, j);

    out.alpha.assign(n, zero);
    T num = zero;
    for (std::size_t i = 0; i < n; ++i) {
        out.alpha[i] = theta[i] * completion[i];
        num += out.alpha[i];
    }
    if (num > zero) {
        for (std::size_t i = 0; i < n; ++i) out.alpha[i] = out.alpha[i] / num;
    } else {
        out.alpha = theta;  // no completions anywhere: fall back on theta
        T s = zero;
        for (std::size_t i = 0; i < n; ++i) s += out.alpha[i];
        if (s == zero) throw NumericError("aoi_dist2ph: alpha cannot be normalized");
        for (std::size_t i = 0; i < n; ++i) out.alpha[i] = out.alpha[i] / s;
    }

    for (std::size_t i = 0; i < n; ++i)
        if (out.Tmat(i, i) > zero)
            throw InputError("aoi_dist2ph: T has a positive diagonal entry, which is not a "
                             "sub-generator");
    return out;
}

/** Convenience overload taking the MAP as line::mam::Map. */
template <class T>
AoiPh<T> aoi_dist2ph(const mam::Map<T>& m) {
    return aoi_dist2ph(m.D0, m.D1);
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_AOI_DIST2PH_H
