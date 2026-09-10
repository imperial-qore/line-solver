/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GGISGI_FLUID_H
#define LINE_API_QSYS_GGISGI_FLUID_H

/**
 * Steady state of the G/GI/s+GI fluid model.
 *
 * Templated port of matlab/src/api/qsys/qsys_ggisgi_fluid.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_ggisgi_fluid.java.
 *
 * Scale the content by s and let s grow. Customers become quanta of fluid but
 * their sojourns do not shrink, so the ages survive the limit: the state is the
 * density b(x) of fluid in service of age x and the density q(x) of fluid
 * waiting of age x. With rho = lambda/(s mu),
 *
 *   rho <= 1   b(x) = rho G^c(x),  q = 0                       (3.2)
 *   rho >  1   b(x) = G^c(x),  q(x) = rho F^c(x) on [0,w]      (3.4)-(3.5)
 *
 * with the queue boundary w solving F^c(w) = 1/rho (3.6): fluid that survives
 * its patience for w enters service, so the surviving fraction must equal the
 * fraction 1/rho the servers can absorb. Then
 *
 *   P(abandon) = 1 - 1/rho,  W = int_0^w F^c = m_a F_e(w),  Q = lambda W.
 *
 * WHAT THE DISTRIBUTIONS CONTRIBUTE (Corollary 3.1): the rates and the number in
 * service depend on G and F only through their means; w, Q and the queue age
 * profile depend on F beyond its mean but on G only through its mean. Neither s
 * nor anything about the arrival process beyond its rate appears.
 *
 * ARITHMETIC. The boundary w comes out of a bisection against a tolerance and
 * the mean wait out of a Simpson quadrature, so this is an approximation of an
 * approximation and there is nothing to gain from exact arithmetic; the
 * instantiation is therefore restricted to the transcendental types.
 *
 * Reference: W. Whitt (2006). Fluid models for multiserver queues with
 * abandonments. Operations Research 54(1), 37-54.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Steady state of the G/GI/s+GI fluid model. */
template <class T>
struct QsysFluidAbandonResult {
    std::string regime;              ///< "underloaded", "balanced" or "overloaded"
    T trafficIntensity;              ///< rho = lambda/(s mu)
    T offeredWait;                   ///< w, the wait of every served customer
    T meanWait;                      ///< E[W] over all customers
    T meanWaitServed;                ///< w again, the fluid wait being deterministic
    T meanWaitAbandon;               ///< E[patience | patience <= w]
    T probAbandon;                   ///< 1 - 1/rho when overloaded
    T meanQueueLength;               ///< Q = lambda W, in customers
    T meanNumberInService;           ///< B = min(lambda/mu, s), in customers
    T meanNumber;                    ///< B + Q
    T utilization;                   ///< min(rho,1)
    T throughput;                    ///< min(lambda, s mu)
    T abandonRate;                   ///< lambda - throughput
    std::vector<T> agePoints;        ///< ages of the densities, empty when not requested
    std::vector<T> serviceAgeDensity;///< b(x) per server
    std::vector<T> queueAgeDensity;  ///< q(x) per server
};

namespace detail {

/**
 * Smallest w with F^c(w) = target, by doubling then bisection. F^c is
 * non-increasing, so the doubling either brackets the crossing or proves that
 * the patience law never decays that far.
 */
template <class T, class Ccdf>
T fluid_inv_ccdf(Ccdf&& ccdf, const T& target, double tol, double maxTime) {
    const T zero = num_traits<T>::from_int(0);
    if (ccdf(zero) < target)
        throw InputError("qsys_ggisgi_fluid: the patience ccdf is below 1/rho at t = 0, so it is "
                         "not a ccdf");
    T lo = zero;
    T hi;
    if (std::isnan(maxTime)) {
        hi = num_traits<T>::from_int(1);
        while (ccdf(hi) > target) {
            hi *= num_traits<T>::from_int(2);
            if (hi > num_traits<T>::from_double(1e12))
                throw InputError("qsys_ggisgi_fluid: the patience ccdf never falls to 1/rho, so "
                                 "the overloaded fluid model has no equilibrium: too little of the "
                                 "fluid is willing to abandon");
        }
    } else {
        hi = num_traits<T>::from_double(maxTime);
        if (ccdf(hi) > target)
            throw InputError("qsys_ggisgi_fluid: the patience ccdf is still above 1/rho at maxTime");
    }
    const T two = num_traits<T>::from_int(2);
    const T tolT = num_traits<T>::from_double(tol);
    const T one = num_traits<T>::from_int(1);
    while (hi - lo > tolT * (hi > one ? hi : one)) {
        const T mid = (lo + hi) / two;
        if (ccdf(mid) > target) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return (lo + hi) / two;
}

/**
 * Composite Simpson rule on a fixed fine grid: the integrand is a ccdf, hence
 * monotone and bounded, so a fixed grid is enough and is reproducible.
 */
template <class T, class Fn>
T fluid_integral(Fn&& f, const T& a, const T& b) {
    const T zero = num_traits<T>::from_int(0);
    if (b <= a) return zero;
    const std::size_t n = 2000;
    const T h = (b - a) / num_traits<T>::from_int(static_cast<long>(n));
    T sum = f(a) + f(b);
    for (std::size_t i = 1; i < n; ++i) {
        const T w = num_traits<T>::from_int((i % 2 == 1) ? 4 : 2);
        sum += w * f(T(a + num_traits<T>::from_int(static_cast<long>(i)) * h));
    }
    return h / num_traits<T>::from_int(3) * sum;
}

}  // namespace detail

/**
 * @param lambda       arrival rate
 * @param mu           service rate of one server
 * @param s            number of servers, s >= 1
 * @param patienceCcdf F^c(t) = P(patience > t)
 * @param servingCcdf  G^c(x) = P(service > x), used only for the in-service age
 *                     density; an empty callable takes the exponential of rate mu
 * @param agePoints    ages at which to return the densities
 * @param tol          bisection tolerance for w
 * @param maxTime      largest age searched for w; NaN grows the search
 */
template <class T>
QsysFluidAbandonResult<T> qsys_ggisgi_fluid(
    const T& lambda, const T& mu, unsigned s, const std::function<T(const T&)>& patienceCcdf,
    const std::function<T(const T&)>& servingCcdf = std::function<T(const T&)>(),
    const std::vector<T>& agePoints = std::vector<T>(), double tol = 1e-12,
    double maxTime = std::numeric_limits<double>::quiet_NaN()) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_ggisgi_fluid bisects against a tolerance, so it needs inexact arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (lambda <= zero) throw InputError("qsys_ggisgi_fluid: the arrival rate lambda must be positive");
    if (mu <= zero) throw InputError("qsys_ggisgi_fluid: the service rate mu must be positive");
    if (s < 1) throw InputError("qsys_ggisgi_fluid: the number of servers s must be at least 1");
    if (!patienceCcdf) throw InputError("qsys_ggisgi_fluid: the patience ccdf must be callable");

    std::function<T(const T&)> gc = servingCcdf;
    if (!gc) {
        const T fmu = mu;
        gc = [fmu](const T& x) {
            using std::exp;
            return exp(-fmu * x);
        };
    }

    const T sT = num_traits<T>::from_int(static_cast<long>(s));
    const T rho = lambda / (sT * mu);
    QsysFluidAbandonResult<T> res;
    res.trafficIntensity = rho;

    T w = zero, meanWait = zero, probAbandon = zero, meanWaitAbandon = zero;
    if (rho <= one) {
        // Underloaded and balanced, eq. (3.2): the queue is empty and the model
        // is the infinite-server fluid model.
        res.regime = (rho == one) ? "balanced" : "underloaded";
    } else {
        res.regime = "overloaded";
        w = detail::fluid_inv_ccdf<T>(patienceCcdf, T(one / rho), tol, maxTime);   // eq. (3.6)
        // Eq. (3.14): W = int_0^w F^c(t) dt = m_a F_e(w), over ALL fluid.
        meanWait = detail::fluid_integral<T>(patienceCcdf, zero, w);
        probAbandon = one - one / rho;
        // E[T | T <= w] = (W - w F^c(w)) / F(w) by parts, F^c(w) = 1/rho.
        meanWaitAbandon = (meanWait - w / rho) / probAbandon;
    }

    res.offeredWait = w;
    res.meanWait = meanWait;
    res.meanWaitServed = w;
    res.meanWaitAbandon = meanWaitAbandon;
    res.probAbandon = probAbandon;
    res.meanQueueLength = lambda * meanWait;                    // eq. (3.11), Little's law
    res.meanNumberInService = detail::num_min(T(lambda / mu), sT);
    res.meanNumber = res.meanNumberInService + res.meanQueueLength;
    res.utilization = detail::num_min(rho, one);
    res.throughput = detail::num_min(lambda, T(sT * mu));
    res.abandonRate = lambda - res.throughput;

    if (!agePoints.empty()) {
        const T sigma = detail::num_min(rho, one);              // rate into service, per server
        res.agePoints = agePoints;
        res.serviceAgeDensity.resize(agePoints.size());
        res.queueAgeDensity.resize(agePoints.size());
        for (std::size_t i = 0; i < agePoints.size(); ++i) {
            res.serviceAgeDensity[i] = sigma * gc(agePoints[i]);
            res.queueAgeDensity[i] = (rho > one && agePoints[i] <= w)
                                         ? T(rho * patienceCcdf(agePoints[i]))
                                         : zero;
        }
    }
    return res;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GGISGI_FLUID_H
