/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_SENS_H
#define LINE_API_MC_CTMC_SENS_H

/**
 * Sensitivity of the steady-state distribution of a CTMC to a scalar
 * parameter.
 *
 * Templated port of matlab/src/api/mc/ctmc_sens.m. Differentiating pi Q = 0 and
 * pi e = 1 with respect to theta gives Trivedi and Bobbio (2017), Eq. (9.81),
 *
 *   (dpi/dtheta) Q = -pi (dQ/dtheta),   sum_i dpi_i/dtheta = 0,
 *
 * whose coefficient matrix is the one the steady-state solve already assembles:
 * a sensitivity costs exactly one extra triangular solve. The normalization
 * replaces the LAST equation of the transposed system, as in ctmc_solve.
 *
 * The whole computation is a linear solve over the field of the rates, so it
 * carries NO gate: at Rational it returns the exact derivative of the exact
 * stationary vector, which is the regime where a finite-difference estimate of
 * the same quantity is worst behaved (it differences two nearly equal vectors)
 * and where this port is therefore most worth having.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * @param Q  generator
 * @param dQ derivative of the generator with respect to theta, same size
 * @param pi steady-state distribution; must sum to one
 * @return dpi/dtheta as a row vector of length n, summing to zero
 */
template <class T>
std::vector<T> ctmc_sens(const Matrix<T>& Q, const Matrix<T>& dQ, const std::vector<T>& pi) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_sens: generator is not square");
    if (dQ.rows() != n || dQ.cols() != n) throw InputError("ctmc_sens: dQ must have the same size as Q");
    if (pi.size() != n) throw InputError("ctmc_sens: pi has the wrong length");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    // b = -pi * dQ, with the last entry replaced by the normalization 0.
    std::vector<T> b(n, zero);
    for (std::size_t j = 0; j < n; ++j) {
        T s = zero;
        for (std::size_t i = 0; i < n; ++i) s += pi[i] * dQ(i, j);
        b[j] = -s;
    }
    b[n - 1] = zero;

    // A = Q' with its last ROW replaced by ones (MATLAB A(n,:) = ones(1,n)).
    Matrix<T> A(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = Q(j, i);
    for (std::size_t j = 0; j < n; ++j) A(n - 1, j) = one;

    return solve(A, b);
}

/** Overload computing the steady-state vector itself, as MATLAB does. */
template <class T>
std::vector<T> ctmc_sens(const Matrix<T>& Q, const Matrix<T>& dQ) {
    return ctmc_sens(Q, dQ, ctmc_solve(Q));
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_SENS_H
