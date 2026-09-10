/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_ME_ME_GEGECN_H
#define LINE_API_ME_ME_GEGECN_H

/**
 * Censored GE/GE/c/K;N queue by entropy maximisation.
 *
 * Templated port of `matlab/src/api/me/me_gegecn.m`: the single-class censored
 * FCFS queue of Kouvatsos (1994), Section 4.1, equations (4.1)-(4.3).
 *
 * WHAT "CENSORED" MEANS HERE. The queue holds at most N jobs and never fewer
 * than K: arrivals finding N jobs present are turned away, and departures are
 * not allowed from state K. For a queue embedded in an OPEN network K is always
 * 0; a positive K arises in closed networks, where it records the minimum
 * occupancy forced by the remaining stations being full.
 *
 * THE SOLUTION IS CLOSED FORM, NOT AN ITERATION. The ME state probabilities
 * coincide with the global balance solution
 *
 *     p(n) = p(K) G_n x^h(n) y^f(n),      n = K+1,...,N
 *
 * with G_n = prod_{l=K+1}^{m(n)} g(l), J = max(c,K+1), h(n) = max(0,n-J),
 * f(n) = max(0,n-N+1) and m(n) = max{K+1, min(c,n)}. The Lagrangian
 * coefficients g(l), x and y come from raw system data. They are INVARIANT to
 * N and K, which is why letting K -> 0 and N -> infinity recovers the stable
 * GE/GE/c solution `me_oqn` uses -- the same coefficients serve both.
 *
 * WHY THIS PORT WORKS IN LOGS, as the reference does. `x^(N-J)` overflows on a
 * saturated queue with a large buffer, and the log form additionally makes the
 * rho = 1 case (x = 1) fall out of the same expression instead of needing the
 * separate p(K) branch of (4.2).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/me/me_gegecn_pb.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace me {

/** What `me_gegecn` returns: the law and the four means read off it. */
template <class T>
struct GegecnResult {
    std::vector<T> p;  ///< p[idx] = Pr{n = K + idx}, idx = 0..N-K
    T L;               ///< mean number in the queue, sum_n n p(n)
    T U;               ///< utilization, E[min(n,c)]/c
    T PB;              ///< probability an arrival of the aggregate stream is blocked
    T Lq;              ///< mean number WAITING, L - E[min(n,c)]
};

/**
 * Port of `me_gegecn`.
 *
 * @param lambda arrival rate OFFERED to the queue, the arrivals turned away
 *               included
 * @param Ca     squared coefficient of variation of the interarrival times; the
 *               GE distribution needs at least 1
 * @param mu     service rate of ONE server
 * @param Cs     squared coefficient of variation of the service times, >= 1
 * @param c      number of servers, finite and at least 1
 * @param K      minimum number of jobs in the queue, >= 0
 * @param N      buffer capacity in jobs, service included, finite and > K
 */
template <class T>
GegecnResult<T> me_gegecn(const T& lambda, const T& Ca, const T& mu, const T& Cs, long c, long K,
                          long N) {
    static_assert(num_traits<T>::has_transcendental,
                  "me_gegecn requires transcendental arithmetic: the state law is assembled in "
                  "logarithms so that a saturated queue with a large buffer does not overflow");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    if (c < 1) throw InputError("me_gegecn: requires a finite number of servers c >= 1.");
    if (N <= K) throw InputError("me_gegecn: requires N > K.");
    if (num_traits<T>::to_double(Ca) < 1.0 - 1e-12 ||
        num_traits<T>::to_double(Cs) < 1.0 - 1e-12)
        throw InputError(
            "me_gegecn: requires Ca >= 1 and Cs >= 1: the GE distribution is not defined for "
            "scv < 1.");
    if (!(mu > zero)) throw InputError("me_gegecn: requires a positive service rate.");

    const T tau = T(two / (Ca + one));
    const T sigma = T(two / (Cs + one));
    const T cT = num_traits<T>::from_int(c);
    const T rho = T(lambda / (cT * mu));

    const long J = std::max(c, K + 1);
    const T den1 = T(sigma * (one - tau) + tau);
    const T den2 = T(tau * rho * (one - sigma) + sigma);

    // Lagrangian coefficients g(l), l = K+1,...,J, stored 1-based as in the
    // reference so the three-way first-entry branch reads the same.
    std::vector<T> g(static_cast<std::size_t>(J) + 1, one);
    const T Kp1 = num_traits<T>::from_int(K + 1);
    if (K < c - 1) {
        g[static_cast<std::size_t>(K + 1)] = T(tau * cT * rho / (Kp1 * den1));
    } else if (K == c - 1) {
        g[static_cast<std::size_t>(K + 1)] = T(tau * sigma * rho / den2);
    } else {
        g[static_cast<std::size_t>(K + 1)] = T((den1 / den2) * tau * rho);
    }
    for (long l = K + 2; l <= J; ++l) {
        const T lT = num_traits<T>::from_int(l);
        const T lm1 = num_traits<T>::from_int(l - 1);
        if (l < J) {
            g[static_cast<std::size_t>(l)] =
                T((tau * cT * rho + lm1 * sigma * (one - tau)) / (lT * den1));
        } else {
            const T Jm1 = num_traits<T>::from_int(J - 1);
            const T JT = num_traits<T>::from_int(J);
            g[static_cast<std::size_t>(l)] =
                T(sigma * (tau * cT * rho + Jm1 * sigma * (one - tau)) / (JT * den2));
        }
    }

    const T x = T((tau * rho + sigma * (one - tau)) / den2);
    const T y = T(one / (one - (one - sigma) * x));

    // cumlogg(i) = sum_{l=K+1}^{K+i} log g(l), i = 1..J-K
    const std::size_t ng = static_cast<std::size_t>(J - K);
    std::vector<T> cumlogg(ng, zero);
    {
        T acc = zero;
        for (std::size_t i = 0; i < ng; ++i) {
            acc += T(log(g[static_cast<std::size_t>(K + 1) + i]));
            cumlogg[i] = acc;
        }
    }

    const std::size_t nn = static_cast<std::size_t>(N - K + 1);
    const T logx = T(log(x));
    const T logy = T(log(y));
    std::vector<T> logp(nn, zero);
    for (std::size_t idx = 0; idx < nn; ++idx) {
        const long n = K + static_cast<long>(idx);
        if (n > K) {
            const long m = std::max(K + 1, std::min(c, n));
            logp[idx] = cumlogg[static_cast<std::size_t>(m - K) - 1];
        }
        const long h = std::max<long>(0, n - J);
        const long f = std::max<long>(0, n - N + 1);
        if (h > 0) logp[idx] += num_traits<T>::from_int(h) * logx;
        if (f > 0) logp[idx] += num_traits<T>::from_int(f) * logy;
    }
    // Shift by the maximum before exponentiating, then normalize.
    T mx = logp[0];
    for (const T& v : logp)
        if (v > mx) mx = v;
    GegecnResult<T> res;
    res.p.assign(nn, zero);
    T tot = zero;
    for (std::size_t idx = 0; idx < nn; ++idx) {
        res.p[idx] = T(exp(T(logp[idx] - mx)));
        tot += res.p[idx];
    }
    for (T& v : res.p) v = T(v / tot);

    T L = zero, Ebusy = zero;
    for (std::size_t idx = 0; idx < nn; ++idx) {
        const long n = K + static_cast<long>(idx);
        L += num_traits<T>::from_int(n) * res.p[idx];
        Ebusy += num_traits<T>::from_int(std::min<long>(n, c)) * res.p[idx];
    }
    res.L = L;
    res.U = T(Ebusy / cT);
    res.Lq = T(L - Ebusy);
    res.PB = me_gegecn_pb(res.p, K, N, c, Cs, Ca);
    return res;
}

}  // namespace me
}  // namespace line

#endif  // LINE_API_ME_ME_GEGECN_H
