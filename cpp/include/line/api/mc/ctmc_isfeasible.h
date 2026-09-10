/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_ISFEASIBLE_H
#define LINE_API_MC_CTMC_ISFEASIBLE_H

/**
 * Feasibility predicates for generators and stochastic matrices.
 *
 * Templated port of matlab/src/api/mc/ctmc_isfeasible.m and
 * matlab/lib/kpctoolbox/mc/dtmc_isfeasible.m.
 *
 * The two return different things and that difference is deliberate. The CTMC
 * predicate is a BOOLEAN at a caller-supplied tolerance. The DTMC one returns a
 * PRECISION LEVEL: the largest tol in 1..15 at which the row sums are within
 * 10^-tol of one and the entries are above -10^-tol, and 0 when no level holds.
 * A caller that reads the DTMC result as a boolean is right by accident, since
 * any nonzero level is truthy, but a caller that compares it to 1 is wrong.
 */

#include <cstddef>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** True when Q is square, has nonnegative off-diagonals, nonpositive diagonal and zero row sums. */
template <class T>
bool ctmc_isfeasible(const Matrix<T>& Q, const T& tol) {
    const std::size_t n = Q.rows();
    if (n == 0 || Q.cols() != n) return false;
    for (std::size_t i = 0; i < n; ++i) {
        T rowsum = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j && Q(i, j) < -tol) return false;
            rowsum += Q(i, j);
        }
        if (Q(i, i) > tol) return false;
        if (num_abs(rowsum) > tol) return false;
    }
    return true;
}

/** Default tolerance 1e-10, matching the MATLAB signature. */
template <class T>
bool ctmc_isfeasible(const Matrix<T>& Q) {
    return ctmc_isfeasible(Q, num_traits<T>::from_double(1e-10));
}

/** Largest precision level 1..15 at which P is stochastic, or 0 when none holds. */
template <class T>
int dtmc_isfeasible(const Matrix<T>& P) {
    const std::size_t n = P.rows();
    if (n == 0) return 0;
    T minsum = num_traits<T>::from_int(0), maxsum = num_traits<T>::from_int(0);
    T minel = P(0, 0);
    for (std::size_t i = 0; i < n; ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < P.cols(); ++j) {
            s += P(i, j);
            if (P(i, j) < minel) minel = P(i, j);
        }
        if (i == 0 || s < minsum) minsum = s;
        if (i == 0 || s > maxsum) maxsum = s;
    }
    int res = 0;
    const T one = num_traits<T>::from_int(1);
    T eps = num_traits<T>::from_int(1);
    const T ten = num_traits<T>::from_int(10);
    for (int tol = 1; tol <= 15; ++tol) {
        eps = eps / ten;
        if (minsum > one - eps && maxsum < one + eps && minel > -eps) res = tol;
    }
    return res;
}

}  // namespace mc
}  // namespace line

#endif
