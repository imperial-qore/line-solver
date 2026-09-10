/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_TIMEREVERSE_H
#define LINE_API_MC_CTMC_TIMEREVERSE_H

/**
 * Time-reversed generator and transition matrix.
 *
 * Templated port of matlab/lib/kpctoolbox/mc/ctmc_timereverse.m and
 * dtmc_timereverse.m. The reversed chain is Qrev = diag(pi)^-1 Q' diag(pi),
 * written here as the transpose of the elementwise scaling Q(i,j) pi_i / pi_j,
 * which is the same matrix and is how MATLAB assembles it.
 *
 * The construction needs the stationary law, so it is defined only for an
 * irreducible chain: on a reducible one the stationary vector is not unique and
 * the reversal is not either.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** Generator of the time-reversed CTMC. */
template <class T>
Matrix<T> ctmc_timereverse(const Matrix<T>& Q) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_timereverse: Q is not square");
    const std::vector<T> pie = ctmc_solve(Q);
    Matrix<T> Qrev(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Qrev(j, i) = Q(i, j) * pie[i] / pie[j];
    return Qrev;
}

/** Transition matrix of the time-reversed DTMC. */
template <class T>
Matrix<T> dtmc_timereverse(const Matrix<T>& P) {
    const std::size_t n = P.rows();
    if (P.cols() != n) throw InputError("dtmc_timereverse: P is not square");
    const std::vector<T> pie = dtmc_solve(P);
    Matrix<T> Prev(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Prev(j, i) = P(i, j) * pie[i] / pie[j];
    return Prev;
}

}  // namespace mc
}  // namespace line

#endif
