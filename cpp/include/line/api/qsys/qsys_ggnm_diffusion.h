/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GGNM_DIFFUSION_H
#define LINE_API_QSYS_GGNM_DIFFUSION_H

/**
 * Diffusion approximation for the G/GI/n/m queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_ggnm_diffusion.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_ggnm_diffusion.java.
 *
 * ONE DIFFUSION WITH TWO REGIONS. Below the staffing level the queue behaves
 * like an infinite-server system, whose limit is NORMAL with variance-to-mean
 * ratio the ASYMPTOTIC PEAKEDNESS
 *
 *   z = 1 + (ca^2 - 1) omega_G,  omega_G = int G^c(x)^2 dx / int G^c(x) dx
 *                                                             (1.6)-(1.7)
 *
 * above it like a single-server queue, whose limit is EXPONENTIAL with
 * variability v = (ca^2 + cs^2)/2 (3.7). The steady-state law is a normal piece
 * spliced to an exponential piece and every measure is an integral of it (3.14):
 *
 *   P(delay) = [1 + b Phi(b)/(phi(b)(1 - e^{-beta gamma/v}))]^-1, b = beta/sqrt(z)
 *   P(block) = f(gamma) v / sqrt(n)                              (7.5)
 *
 * with beta = sqrt(n)(1-rho) and gamma = m/sqrt(n).
 *
 * WHAT z SAYS. The service law enters the delay probability ONLY through
 * omega_G: 1 for deterministic service, 1/2 for exponential, falling toward 0 as
 * service gets more variable. At ca^2 = 1 the delay probability does not depend
 * on the service law at all (z = 1), the long-standing M/GI/n-by-M/M/n
 * approximation; away from ca^2 = 1 it does, and this is how much.
 *
 * At m = Inf this reduces to alpha(beta/sqrt(z)) with the Halfin-Whitt alpha,
 * i.e. to `qsys_mmk_qed` when the service is exponential and ca^2 = 1.
 *
 * ARITHMETIC. erfc, exp and a quadrature: transcendental only.
 *
 * DIVERGENCE. The finite-waiting-room delay function is eq. (3.2) of the paper,
 * whose printed form the available scan does not resolve. What is implemented is
 * the unique form that (i) reduces to eq. (3.10) as gamma -> Inf and (ii)
 * reproduces the exact M/M/n/m delay probability in the QED limit, which was
 * checked numerically against the birth-death chain at n = 100, 400 and 1000.
 *
 * Reference: W. Whitt (2004). A diffusion approximation for the G/GI/n/m queue.
 * Operations Research 52(6), 922-941.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Steady-state measures of the G/GI/n/m diffusion approximation. */
template <class T>
struct QsysGgnmResult {
    T beta;              ///< the QED server slack sqrt(n)(1-rho)
    T gamma;             ///< the scaled waiting room m/sqrt(n)
    T peakedness;        ///< z, the asymptotic peakedness
    T peakednessWeight;  ///< omega_G
    T variability;       ///< v = (ca^2+cs^2)/2
    T probDelay;         ///< P(an arrival waits)
    T probBlock;         ///< P(an arrival is blocked)
    T meanQueueLength;   ///< mean number waiting
    T meanNumber;        ///< mean number in system
    T meanWait;          ///< mean wait of an admitted arrival
    T utilization;       ///< min(rho,1)
    T throughput;        ///< lambda(1-P(block))
    T trafficIntensity;  ///< rho = lambda/(n mu)
};

namespace detail {

/** Standard normal density. */
template <class T>
T ggnm_phi(const T& x) {
    using std::exp;
    using std::sqrt;
    return exp(-x * x / num_traits<T>::from_int(2)) /
           sqrt(num_traits<T>::from_int(2) * num_traits<T>::from_double(M_PI));
}

/** Standard normal cdf, through erfc. */
template <class T>
T ggnm_Phi(const T& x) {
    using std::erfc;
    using std::sqrt;
    return erfc(-x / sqrt(num_traits<T>::from_int(2))) / num_traits<T>::from_int(2);
}

/**
 * omega_G of eq. (1.7) by Simpson on a grid cut where the ccdf is negligible.
 * The denominator is E[S], so only the numerator is integrated.
 */
template <class T, class Ccdf>
T ggnm_omega(Ccdf&& ccdf, const T& ES, double tol, std::size_t panels) {
    T hi = num_traits<T>::from_int(1);
    const T tolT = num_traits<T>::from_double(tol);
    while (ccdf(hi) > tolT) {
        hi *= num_traits<T>::from_int(2);
        if (hi > num_traits<T>::from_double(1e12))
            throw InputError("qsys_ggnm_diffusion: the service ccdf does not decay, so its "
                             "peakedness is undefined");
    }
    const T h = hi / num_traits<T>::from_int(static_cast<long>(panels));
    const T g0 = ccdf(num_traits<T>::from_int(0));
    const T gn = ccdf(hi);
    T sum = g0 * g0 + gn * gn;
    for (std::size_t i = 1; i < panels; ++i) {
        const T g = ccdf(T(num_traits<T>::from_int(static_cast<long>(i)) * h));
        sum += num_traits<T>::from_int(i % 2 == 1 ? 4 : 2) * g * g;
    }
    return (h / num_traits<T>::from_int(3) * sum) / ES;
}

}  // namespace detail

/**
 * @param lambda      arrival rate
 * @param mu          service rate of one server
 * @param n           number of servers, n >= 1
 * @param m           extra waiting spaces; infinity for an unbounded queue
 * @param ca          coefficient of variation of the interarrival time
 * @param cs          coefficient of variation of the service time
 * @param serviceCcdf G^c(x) = P(S > x); empty takes the exponential of rate mu
 * @param tol         service-tail cut for the peakedness integral
 * @param panels      Simpson panels for it
 */
template <class T>
QsysGgnmResult<T> qsys_ggnm_diffusion(
    const T& lambda, const T& mu, unsigned n, double m, const T& ca, const T& cs,
    const std::function<T(const T&)>& serviceCcdf = std::function<T(const T&)>(),
    double tol = 1e-12, std::size_t panels = 4000) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_ggnm_diffusion needs erfc, exp and a quadrature");
    using std::exp;
    using std::expm1;
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (lambda <= zero || mu <= zero)
        throw InputError("qsys_ggnm_diffusion: the arrival and service rates must be positive");
    if (n < 1) throw InputError("qsys_ggnm_diffusion: the number of servers n must be at least 1");
    if (m < 0)
        throw InputError("qsys_ggnm_diffusion: the number of extra waiting spaces m must be "
                         "non-negative");

    const T nT = num_traits<T>::from_int(static_cast<long>(n));
    const T ca2 = ca * ca;
    const T cs2 = cs * cs;
    const T ES = one / mu;
    const T rho = lambda / (nT * mu);
    const T beta = sqrt(nT) * (one - rho);                     // eq. (0.1)
    const bool finiteRoom = std::isfinite(m);
    const T gamma = finiteRoom ? T(num_traits<T>::from_double(m) / sqrt(nT))
                               : num_traits<T>::from_double(std::numeric_limits<double>::infinity());

    const T omega = serviceCcdf ? detail::ggnm_omega<T>(serviceCcdf, ES, tol, panels)
                                : num_traits<T>::from_rational(1, 2);
    const T z = one + (ca2 - one) * omega;                     // eq. (1.6)
    if (z <= zero)
        throw InputError("qsys_ggnm_diffusion: the asymptotic peakedness came out non-positive; "
                         "check ca and the service ccdf");
    const T v = (ca2 + cs2) / two;                             // eq. (3.7), weight w = 1
    const T b = beta / sqrt(z);
    const T r = beta / v;

    // The tail factor is negative together with r when the queue is overloaded,
    // so the ratio below stays positive on both sides of beta = 0.
    const T tail = finiteRoom ? T(-expm1(-r * gamma)) : one;
    T alpha, meanAbove, densityAtTop;
    if (num_abs(r) < num_traits<T>::from_double(1e-14)) {
        // beta = 0: the exponential piece degenerates to a uniform on [0,gamma].
        if (!finiteRoom)
            throw InputError("qsys_ggnm_diffusion: with beta = 0 the queue needs a finite waiting "
                             "room to be stable");
        alpha = one / (one + detail::ggnm_Phi(b) / (detail::ggnm_phi(b) * gamma / sqrt(z)));
        meanAbove = gamma / two;
        densityAtTop = alpha / gamma;
    } else {
        alpha = one / (one + b * detail::ggnm_Phi(b) / (detail::ggnm_phi(b) * tail));
        if (finiteRoom) {
            const T e = exp(-r * gamma);
            meanAbove = (one / r - (gamma + one / r) * e) / tail;
            densityAtTop = alpha * r * e / tail;
        } else {
            meanAbove = one / r;
            densityAtTop = zero;
        }
    }

    // Mean of the normal piece, N(-beta, z) conditioned below 0.
    const T meanBelow = -beta - sqrt(z) * detail::ggnm_phi(b) / detail::ggnm_Phi(b);
    const T meanScaled = (one - alpha) * meanBelow + alpha * meanAbove;

    QsysGgnmResult<T> res;
    res.beta = beta;
    res.gamma = gamma;
    res.peakedness = z;
    res.peakednessWeight = omega;
    res.variability = v;
    res.probDelay = alpha;
    // Eq. (7.5): the loss rate at the upper boundary over the arrival rate is
    // the density there times v / sqrt(n).
    res.probBlock = zero;
    if (finiteRoom) {
        T pb = densityAtTop * v / sqrt(nT);
        if (pb < zero) pb = zero;
        if (pb > one) pb = one;
        res.probBlock = pb;
    }
    res.meanQueueLength = sqrt(nT) * alpha * meanAbove;
    res.meanNumber = nT + sqrt(nT) * meanScaled;
    res.throughput = lambda * (one - res.probBlock);
    res.meanWait = res.throughput > zero ? T(res.meanQueueLength / res.throughput) : zero;
    res.utilization = detail::num_min(rho, one);
    res.trafficIntensity = rho;
    return res;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GGNM_DIFFUSION_H
