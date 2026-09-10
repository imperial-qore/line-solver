/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_STMONOTONE_H
#define LINE_API_MC_CTMC_STMONOTONE_H

/**
 * Stochastically monotone upper bound of a Markov chain.
 *
 * Templated port of jar/src/main/java/jline/api/mc/Ctmc_stmonotone.java, which
 * has no MATLAB twin. Given a row-stochastic P the algorithm of Abu-Kamel and
 * Stewart builds the smallest st-monotone Q that dominates P in the strong
 * stochastic order, by filling the tail sums from the last column backwards:
 *   Q(0,l..n-1) = P(0,l..n-1),
 *   Q(i,l..n-1) = max( Q(i-1,l..n-1), P(i,l..n-1) ),
 * and recovering Q(i,l) as that tail sum minus the tail already assigned. Every
 * row tail of Q then dominates the corresponding row tail of P and is
 * nondecreasing in i, which is exactly st-monotonicity.
 *
 * The CTMC entry uniformizes Q at max|Q| (the rate the reference uses), repairs
 * the rounding through dtmc_makestochastic, bounds the embedded chain and maps
 * the bound back to a generator with ctmc_makeinfgen. Uniformization is a
 * bijection between the two, so the bound is preserved.
 *
 * ARITHMETIC: field plus comparisons. Exact under Rational.
 */

#include <cstddef>

#include "line/api/mc/ctmc_randomization.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/dtmc_makestochastic.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * @param P (n x n) row-stochastic transition matrix
 * @return the st-monotone upper bound of P, itself row-stochastic
 */
template <class T>
Matrix<T> dtmc_stmonotone(const Matrix<T>& P) {
    const std::size_t n = P.rows();
    if (n == 0) throw InputError("dtmc_stmonotone: the matrix is empty");
    if (P.cols() != n) throw InputError("dtmc_stmonotone: the matrix is not square");
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> Q(n, n, zero);

    Q(0, n - 1) = P(0, n - 1);
    for (std::size_t i = 1; i < n; ++i)
        Q(i, n - 1) = Q(i - 1, n - 1) > P(i, n - 1) ? Q(i - 1, n - 1) : P(i, n - 1);

    for (std::size_t lp = n - 1; lp-- > 0;) {
        Q(0, lp) = P(0, lp);
        for (std::size_t i = 1; i < n; ++i) {
            T tailPrev = zero, tailP = zero, tailRight = zero;
            for (std::size_t k = lp; k < n; ++k) {
                tailPrev += Q(i - 1, k);
                tailP += P(i, k);
                if (k > lp) tailRight += Q(i, k);
            }
            Q(i, lp) = T((tailPrev > tailP ? tailPrev : tailP) - tailRight);
        }
    }
    return Q;
}

/**
 * @param Q (n x n) generator
 * @return the generator of the st-monotone upper bound of Q
 */
template <class T>
Matrix<T> ctmc_stmonotone(const Matrix<T>& Q) {
    const T m = ctmc_maxabs(Q);
    // An all-zero generator has no rate to uniformize at; any positive q gives P = Id.
    const T q = (m == num_traits<T>::from_int(0)) ? num_traits<T>::from_int(1) : m;
    return ctmc_makeinfgen(dtmc_stmonotone(ctmc_randomization(Q, q).P));
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_STMONOTONE_H
