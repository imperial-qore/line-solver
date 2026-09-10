/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMDP_ISFEASIBLE_H
#define LINE_API_MAM_MMDP_ISFEASIBLE_H

/**
 * Feasibility predicate for a Markov-modulated deterministic process.
 *
 * Templated port of matlab/src/api/mam/mmdp_isfeasible.m. An MMDP is the pair
 * (Q, R) with Q the generator of the modulating chain and R the diagonal
 * matrix of the deterministic service (or inter-arrival) times attached to
 * each modulating state. The pair is feasible when
 *
 *   - Q is square, has non-positive diagonal, non-negative off-diagonal and
 *     zero row sums, i.e. it is a conservative generator;
 *   - R is square of the same order, diagonal, with non-negative diagonal.
 *
 * ARITHMETIC. The predicate is a finite set of sign and sum comparisons, so
 * it instantiates at every arithmetic including Rational. The tolerance is an
 * explicit argument rather than a hard-wired 1e-10, so the exact
 * instantiation can be given tol = 0 and then rejects any deviation however
 * small, which is the right check for an algebraically assembled pair; the
 * MATLAB default of 1e-10 is what the tolerant overload supplies, and it is
 * the right check for the output of a fit.
 */

#include <cstddef>

#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** True when (Q, R) is a valid MMDP pair up to the given tolerance. */
template <class T>
bool mmdp_isfeasible(const Matrix<T>& Q, const Matrix<T>& R, const T& tol) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = Q.rows();
    if (Q.cols() != n) return false;

    for (std::size_t i = 0; i < n; ++i) {
        if (Q(i, i) > tol) return false;
        for (std::size_t j = 0; j < n; ++j)
            if (i != j && Q(i, j) < -tol) return false;
        T rowsum = zero;
        for (std::size_t j = 0; j < n; ++j) rowsum += Q(i, j);
        if (num_abs(T(rowsum)) > tol) return false;
    }

    if (R.rows() != n || R.cols() != n) return false;

    for (std::size_t i = 0; i < n; ++i) {
        if (R(i, i) < -tol) return false;
        for (std::size_t j = 0; j < n; ++j)
            if (i != j && num_abs(T(R(i, j))) > tol) return false;
    }
    return true;
}

/** mmdp_isfeasible with the MATLAB tolerance, 1e-10. */
template <class T>
bool mmdp_isfeasible(const Matrix<T>& Q, const Matrix<T>& R) {
    return mmdp_isfeasible(Q, R, T(num_traits<T>::from_double(1e-10)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMDP_ISFEASIBLE_H
