/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LINEARIZER_H
#define LINE_API_PFQN_LINEARIZER_H

/**
 * Chandy-Neuse Linearizer for single-server stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_linearizer.m, cross-checked
 * against jar/src/main/java/jline/api/pfqn/mva/Pfqn_linearizer.java. Both
 * reference implementations are one line: the original Linearizer is the
 * extended generalized fixed-point Linearizer with every scaling exponent set
 * to one, so this delegates to pfqn_egflinearizer with alpha == 1.
 *
 * Arithmetic: TRANSCENDENTAL-GATED, inherited from pfqn_egflinearizer. At
 * alpha == 1 the power N_r^alpha_r degenerates to a rational operation, but
 * the algorithm still stops the inner Core loop on
 * enorm(Q_{k+1} - Q_k) < tol, so the returned value still depends on the
 * stopping rule and is not the solution of a finite rational problem.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_egflinearizer.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L       (M x R) service demands
 * @param N       (R) population per class
 * @param Z       (K x R) think times, summed over rows; may be empty
 * @param type    (M) scheduling discipline; carried, see pfqn_egflinearizer
 * @param tol     convergence tolerance
 * @param maxiter total inner-iteration budget
 * @param QN0     (M x R) warm start; may be empty
 */
template <class T>
LinearizerResult<T> pfqn_linearizer(const Matrix<T>& L, const std::vector<int>& N,
                                    const Matrix<T>& Z, const std::vector<SchedStrategy>& type,
                                    double tol, int maxiter, const Matrix<T>& QN0) {
    // alpha == 1 makes the real power the identity, so this branch is field
    // arithmetic throughout; see the gate note in pfqn_egflinearizer.
    const std::vector<T> alpha(N.size(), num_traits<T>::from_int(1));
    return pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alpha, QN0);
}

template <class T>
LinearizerResult<T> pfqn_linearizer(const Matrix<T>& L, const std::vector<int>& N,
                                    const Matrix<T>& Z) {
    return pfqn_linearizer(L, N, Z, std::vector<SchedStrategy>(), 1e-8, 1000, Matrix<T>());
}

template <class T>
LinearizerResult<T> pfqn_linearizer(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_linearizer(L, N, Matrix<T>(), std::vector<SchedStrategy>(), 1e-8, 1000,
                           Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LINEARIZER_H
