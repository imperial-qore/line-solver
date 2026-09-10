/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_ME_ME_GEGEC_MQL_H
#define LINE_API_ME_ME_GEGEC_MQL_H

/**
 * Mean queue length of a stable infinite-capacity GE/GE/c/FCFS queue.
 *
 * Templated port of `matlab/src/api/me/me_gegec_mql.m`, the exact Maximum
 * Entropy solution of Kouvatsos (1994), equation (3.9).
 *
 * WHY IT EXISTS SEPARATELY FROM `me_oqn`. `me_oqn.h` carries a numerically
 * identical local copy of this formula; the reference does the same and says
 * why -- `me_oqn.m` is compiled to a MEX file by MEXIFY, which cannot call out.
 * This file is the one `me_oqn_blk` uses for the stations whose buffer is
 * INFINITE, so a blocking network with a mix of finite and unbounded queues
 * solves both kinds with the same coefficients.
 *
 * THE GEOMETRIC TAIL IS WHY `x` MATTERS. The state probabilities are a product
 * of `g(1..c)` up to the server count and then a geometric ratio `x`
 * thereafter, so `Z` and the two partial sums below are the closed forms of
 * that split: `S1` covers the states with an idle server, `S2` the saturated
 * tail, whose mean needs both `1/(1-x)` and `x/(1-x)^2`.
 *
 * ARITHMETIC. The tail sums assume `|x| < 1`, i.e. a STABLE queue; the caller
 * is responsible for that, exactly as in the reference. Transcendental only.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace me {

/**
 * Port of `me_gegec_mql`.
 *
 * @param lambda arrival rate
 * @param Ca     squared coefficient of variation of the interarrival times
 * @param mu     service rate of ONE server
 * @param Cs     squared coefficient of variation of the service times
 * @param c      number of servers, finite and at least 1
 * @return the mean number of jobs in the queue
 */
template <class T>
T me_gegec_mql(const T& lambda, const T& Ca, const T& mu, const T& Cs, long c) {
    static_assert(num_traits<T>::has_transcendental,
                  "me_gegec_mql requires transcendental arithmetic");
    if (c < 1) throw InputError("me_gegec_mql: the server count must be at least 1");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    const T alpha2 = T(two / (Cs + one));
    const T alpha1 = T(one - alpha2);
    const T beta2 = T(two / (Ca + one));
    const T beta1 = T(one - beta2);
    const T lambda2 = T(beta2 * lambda);
    const T mu2 = T(alpha2 * mu);

    const std::size_t C = static_cast<std::size_t>(c);
    std::vector<T> g(C, one);
    for (std::size_t j = 1; j + 1 <= C - 1 + 1 && j <= C - 1; ++j) {
        const T jj = num_traits<T>::from_int(static_cast<long>(j));
        const T jm1 = num_traits<T>::from_int(static_cast<long>(j) - 1);
        g[j - 1] = T((lambda2 + jm1 * mu2 * beta1) * alpha2 /
                     (jj * mu2 * (one - alpha1 * beta1)));
    }
    const T cc = num_traits<T>::from_int(c);
    const T cm1 = num_traits<T>::from_int(c - 1);
    g[C - 1] = T((lambda2 + cm1 * mu2 * beta1) * alpha2 / (lambda2 * alpha1 + cc * mu2));
    const T x = T((lambda2 + cc * mu2 * beta1) / (lambda2 * alpha1 + cc * mu2));

    // Gn = cumprod(g)
    std::vector<T> Gn(C, one);
    T acc = one;
    for (std::size_t i = 0; i < C; ++i) {
        acc = T(acc * g[i]);
        Gn[i] = acc;
    }

    const T omx = T(one - x);
    T Z = one;
    for (std::size_t i = 0; i + 1 < C; ++i) Z += Gn[i];
    Z += T(Gn[C - 1] / omx);

    T S1 = num_traits<T>::from_int(0);
    for (std::size_t n = 1; n + 1 <= C; ++n)
        S1 += num_traits<T>::from_int(static_cast<long>(n)) * Gn[n - 1];
    const T S2 = T(Gn[C - 1] * (cc / omx + x / (omx * omx)));
    return T((S1 + S2) / Z);
}

}  // namespace me
}  // namespace line

#endif  // LINE_API_ME_ME_GEGEC_MQL_H
