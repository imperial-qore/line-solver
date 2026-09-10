/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MM1_PS_H
#define LINE_API_QSYS_QSYS_MM1_PS_H

/**
 * Exact sojourn-time moments of the multiclass M/M/1-PS queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_mm1_ps.m. No JAR counterpart.
 * Class j arrives Poisson at rate lambda(j) and needs Exp(mu(j)) service on a
 * processor shared equally by every job present, so the class of a job affects
 * its sojourn time both through its own rate and through the MIX of rates it
 * shares the processor with. With alpha = 1 - sum_j lambda_j/mu_j,
 *
 *   E[W_r]   = 1/(alpha mu_r)
 *   E[W_r^2] = 2/(alpha mu_r)^2 * [1 - sum_j lambda_j (mu_j-mu_r)/(mu_j(mu_j+mu_r))]
 *                               / [1 - sum_j lambda_j/(mu_j+mu_r)]
 *
 * which is equation (7) of Mitra and Morrison (1983). Both are EXACT rather
 * than asymptotic: the open system is the N -> infinity limit of the closed
 * terminal-driven system whose moments that paper expands in 1/N, and the
 * leading term of the expansion is exact in the limit. For a single class the
 * second moment reduces to the classical 4/(mu^2 (1-rho)^2 (2-rho)) of
 * Coffman, Muntz and Trotter (1970), which is the identity the test checks.
 *
 * Reference: D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of
 * the Waiting Time in Closed and Open Processor-Sharing Systems with Multiple
 * Job Classes", Adv. Appl. Prob. 15(4), 1983, equation (7).
 *
 * ARITHMETIC: rational in lambda and mu throughout, so the exact
 * instantiation returns both moments with no rounding; there is no
 * transcendental step and no static_assert.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Mirrors MATLAB's [W, W2, alpha] return list. */
template <class T>
struct Mm1PsResult {
    std::vector<T> W;   ///< per-class mean sojourn times
    std::vector<T> W2;  ///< per-class second moments of the sojourn time
    T alpha;            ///< unutilized processor fraction, 1 - sum_j lambda_j/mu_j
};

/**
 * @param lambda per-class Poisson arrival rates, nonnegative
 * @param mu     per-class exponential service rates, positive
 */
template <class T>
Mm1PsResult<T> qsys_mm1_ps(const std::vector<T>& lambda, const std::vector<T>& mu) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const std::size_t R = lambda.size();
    if (mu.size() != R)
        throw InputError("qsys_mm1_ps: lambda and mu must have the same number of classes");
    for (std::size_t j = 0; j < R; ++j) {
        if (lambda[j] < zero) throw InputError("qsys_mm1_ps: lambda must be non-negative");
        if (mu[j] <= zero) throw InputError("qsys_mm1_ps: mu must be positive");
    }

    Mm1PsResult<T> r;
    r.alpha = one;
    for (std::size_t j = 0; j < R; ++j) r.alpha -= lambda[j] / mu[j];
    if (r.alpha <= zero) throw InputError("qsys_mm1_ps: the system is unstable, utilization >= 1");

    r.W.assign(R, zero);
    r.W2.assign(R, zero);
    for (std::size_t rr = 0; rr < R; ++rr) {
        const T mur = mu[rr];
        T num = one, den = one;
        for (std::size_t j = 0; j < R; ++j) {
            num -= lambda[j] * (mu[j] - mur) / (mu[j] * (mu[j] + mur));
            den -= lambda[j] / (mu[j] + mur);
        }
        const T am = r.alpha * mur;
        r.W[rr] = one / am;
        r.W2[rr] = two / (am * am) * num / den;
    }
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MM1_PS_H
