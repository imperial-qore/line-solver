/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_BNDS_EXTREMAL_H
#define LINE_API_QSYS_GIG1_BNDS_EXTREMAL_H

/**
 * Extremal two-moment bounds for the GI/GI/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_bnds_extremal.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_gig1_bnds_extremal.java.
 *
 * Two moments do not determine E[W]; they determine a SET of possible values,
 * and the width of that set is the honest uncertainty in any two-moment
 * approximation. The ends are attained by extremal laws:
 *
 *   lower  D(1) interarrivals and a three-point service law on multiples of it,
 *          E[W] = rho((1+cs^2)rho - 1)^+ / (2(1-rho))                    (2.12)
 *   upper  two-point laws: an interarrival atom at 0, and a service law whose
 *          upper atom runs to infinity as its probability vanishes       (3.2)
 *
 * Making an interarrival time larger only empties the queue once, but making a
 * service time larger delays every customer behind it, which is why the two ends
 * look so different. The upper end reduces to a D(1/p)/RS(D(rho),p)/1 model with
 * p = 1/(1+ca^2), whose mean wait is evaluated by Spitzer's identity
 * sum_n E[Sn^+]/n with Sn = rho(NB(n,1-p)+n) - n/p (Algorithm 1). The closed
 * form (3.4) uses the D/M/1 root delta = exp(-(1-delta)/rho) and is within about
 * 1% of it.
 *
 * ARITHMETIC. The tight bound is a truncated infinite sum and delta comes out of
 * a bisection, so nothing here is exact; the instantiation is restricted to the
 * transcendental types.
 *
 * PARITY. All four codebases form the negative binomial pmf in LOG space,
 *   log P(NB(n,1-p)=k) = lgamma(n+k) - lgamma(k+1) - lgamma(n)
 *                        + n log p + k log(1-p),
 * from one precomputed table of log-gammas, rather than by the ratio recursion
 * Algorithm 1 prints. The two are equivalent, but the recursion accumulates
 * rounding over thousands of multiplications and would leave the ports agreeing
 * only to about 1e-4.
 *
 * Reference: Y. Chen, W. Whitt (2020). Algorithms for the upper bound mean
 * waiting time in the GI/GI/1 queue. Queueing Systems 94, 327-356.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** The bounds, all as TIMES IN QUEUE: add 1/mu for a response time. */
template <class T>
struct Gig1ExtremalResult {
    T trafficIntensity;    ///< rho = lambda/mu
    T lowerBound;          ///< the tight lower bound, eq. (2.12)
    T upperBound;          ///< the conjectured tight upper bound, eq. (3.2)
    T upperBoundClosed;    ///< the closed-form upper bound, eq. (3.4)
    T upperBoundDaley;     ///< Daley's bound, eq. (2.7)
    T upperBoundKingman;   ///< Kingman's bound, eq. (2.6)
    T heavyTraffic;        ///< the heavy-traffic approximation, eq. (2.9)
    T delta;               ///< the D/M/1 root behind upperBoundClosed
    T relativeWidth;       ///< (upper-lower)/upper, what two moments leave undetermined
    bool tightComputed;    ///< whether the O(K*N) bound was evaluated
};

namespace detail {

/**
 * The D/M/1 root of eq. (3.5), delta = exp(-(1-delta)/rho), in (0,1).
 *
 * g(delta) = delta - exp(-(1-delta)/rho) is negative at 0 and positive just
 * below 1, where the second root delta = 1 sits, so bisection on [0,1) finds the
 * wanted root without landing on the trivial one.
 */
template <class T>
T extremal_delta(const T& rho) {
    using std::exp;
    const T two = num_traits<T>::from_int(2);
    T lo = num_traits<T>::from_int(0);
    T hi = num_traits<T>::from_int(1) - num_traits<T>::from_double(1e-15);
    for (int i = 0; i < 200; ++i) {
        const T mid = (lo + hi) / two;
        if (mid - exp(-(num_traits<T>::from_int(1) - mid) / rho) < num_traits<T>::from_int(0)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return (lo + hi) / two;
}

/** Algorithm 1: the mean waiting time of the extremal model, by the log pmf. */
template <class T>
T extremal_tight(const T& rho, const T& ca2, const T& cs2, std::size_t K, std::size_t N) {
    using std::exp;
    using std::lgamma;
    using std::log;
    using std::log1p;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T p = one / (one + ca2);
    const T logp = log(p);
    const T log1mp = log1p(-p);
    std::vector<T> lg(N + K + 2, zero);
    for (std::size_t i = 1; i < lg.size(); ++i)
        lg[i] = lgamma(num_traits<T>::from_int(static_cast<long>(i)));
    T total = rho * ca2 + rho * rho * cs2 / (two * (one - rho));
    for (std::size_t k = 1; k <= K; ++k) {
        T s = zero;
        for (std::size_t n = 1; n <= N; ++n) {
            const T nT = num_traits<T>::from_int(static_cast<long>(n));
            const T step = num_traits<T>::from_int(static_cast<long>(n + k)) * rho - nT / p;
            if (step > zero) {
                const T lpmf = lg[n + k] - lg[k + 1] - lg[n] + nT * logp +
                               num_traits<T>::from_int(static_cast<long>(k)) * log1mp;
                s += exp(lpmf) * step / nT;
            }
        }
        total += s;
    }
    return total;
}

}  // namespace detail

/**
 * @param lambda    arrival rate
 * @param mu        service rate
 * @param ca        coefficient of variation of the interarrival time
 * @param cs        coefficient of variation of the service time
 * @param K         truncation of the negative binomial value
 * @param N         truncation of the random-walk length
 * @param skipTight skip the O(K*N) tight bound and return the closed forms only
 */
template <class T>
Gig1ExtremalResult<T> qsys_gig1_bnds_extremal(const T& lambda, const T& mu, const T& ca,
                                              const T& cs, std::size_t K = 4000,
                                              std::size_t N = 2000, bool skipTight = false) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gig1_bnds_extremal truncates an infinite sum, so it needs inexact arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (lambda <= zero || mu <= zero)
        throw InputError("qsys_gig1_bnds_extremal: the arrival and service rates must be positive");
    const T rho = lambda / mu;
    if (rho >= one)
        throw InputError("qsys_gig1_bnds_extremal: the bounds require a stable queue, rho < 1");
    const T ca2 = ca * ca;
    const T cs2 = cs * cs;
    // The reference sets E[U] = 1, so every waiting time carries the factor
    // 1/lambda, that time unit expressed in the caller's units.
    const T scale = one / lambda;

    Gig1ExtremalResult<T> r;
    r.trafficIntensity = rho;
    const T lowNum = (one + cs2) * rho - one;
    r.lowerBound = scale * rho * (lowNum > zero ? lowNum : zero) / (two * (one - rho));
    r.upperBoundKingman = scale * rho * rho * (ca2 / (rho * rho) + cs2) / (two * (one - rho));
    r.upperBoundDaley =
        scale * rho * rho * ((two - rho) * ca2 / rho + cs2) / (two * (one - rho));
    r.heavyTraffic = scale * rho * rho * (ca2 + cs2) / (two * (one - rho));
    r.delta = detail::extremal_delta(rho);
    r.upperBoundClosed =
        scale * (two * (one - rho) * rho / (one - r.delta) * ca2 + rho * rho * cs2) /
        (two * (one - rho));
    if (skipTight) {
        r.upperBound = r.upperBoundClosed;
        r.tightComputed = false;
    } else {
        r.upperBound = scale * detail::extremal_tight(rho, ca2, cs2, K, N);
        r.tightComputed = true;
    }
    r.relativeWidth = r.upperBound > zero ? T((r.upperBound - r.lowerBound) / r.upperBound) : zero;
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_BNDS_EXTREMAL_H
