/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_GFLINEARIZER_H
#define LINE_API_PFQN_GFLINEARIZER_H

/**
 * Generalized fixed-point Linearizer with a single scaling exponent shared by
 * every class (De Souza e Silva and Muntz).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_gflinearizer.m, cross-checked
 * against jar/src/main/java/jline/api/pfqn/mva/Pfqn_gflinearizer.java. Both
 * references simply broadcast the scalar alpha over the R classes and call the
 * extended form.
 *
 * Arithmetic: TRANSCENDENTAL-GATED, inherited from pfqn_egflinearizer. Here
 * the gate is doubly justified: the inner Core loop stops on a tolerance, and
 * a scalar alpha is a genuine real exponent (pfqn_linearizermx uses 2.0, but
 * nothing constrains it to an integer), so N_r^alpha has no meaning in an
 * exact field.
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
 * @param alpha scaling exponent shared by every class
 * @param L (M x R) service demands
 * @param N (R) population per class
 * @param Z (K x R) think times
 * @param type per-station scheduling strategy
 * @param tol convergence tolerance
 * @param maxiter iteration cap
 * @param QN0 queue lengths that warm-start the iteration; empty for a cold start
 * @see pfqn_egflinearizer for the remaining arguments
 */
template <class T>
LinearizerResult<T> pfqn_gflinearizer(const Matrix<T>& L, const std::vector<int>& N,
                                      const Matrix<T>& Z, const std::vector<SchedStrategy>& type,
                                      double tol, int maxiter, const T& alpha,
                                      const Matrix<T>& QN0) {
    // runtime gating rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    const std::vector<T> alphav(N.size(), alpha);
    return pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alphav, QN0);
}

template <class T>
LinearizerResult<T> pfqn_gflinearizer(const Matrix<T>& L, const std::vector<int>& N,
                                      const Matrix<T>& Z, const T& alpha) {
    return pfqn_gflinearizer(L, N, Z, std::vector<SchedStrategy>(), 1e-8, 1000, alpha,
                             Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_GFLINEARIZER_H
