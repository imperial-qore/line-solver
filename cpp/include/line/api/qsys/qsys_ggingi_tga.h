/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GGINGI_TGA_H
#define LINE_API_QSYS_GGINGI_TGA_H

/**
 * Truncated Gaussian approximation (TGA-G) for the G/GI/n+GI queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_ggingi_tga.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_ggingi_tga.java.
 *
 * A FLUID CENTRE PLUS A GAUSSIAN FLUCTUATION, TRUNCATED. In the
 * efficiency-driven regime (rho > 1 fixed as n grows) the fluid limit gives the
 * centre -- every server busy, w = F^-1(1-1/rho), Q = lambda int_0^w F^c -- and
 * the many-server CLT gives a normal fluctuation of order sqrt(n) around it:
 *
 *   sigma_W^2 = [(ca^2-1) + (cs+1)rho] / (2 mu rho^2 f(w))            (24)
 *   sigma_X^2 = mu^2 sigma_W^2
 *               + lambda int_0^w F^c(u)[1 + (ca^2-1)F^c(u)] du        (11)
 *   W = (w + sigma_W Z/sqrt(n))^+,  Q = (nQ + sqrt(n) sigma_X Z)^+    (18)-(20)
 *
 * Adding fluid and fluctuation directly can produce negative queues and waits,
 * so BOTH ARE TRUNCATED at zero; that truncation is what makes the formulas
 * usable down to moderate overload, reportedly rho > 1.02.
 *
 * The three sources of variability enter separately, which is what lets the
 * exponential-service formula be generalized: the service law appears only as
 * the factor (cs+1)rho, which is 2rho at cs = 1.
 *
 * ARITHMETIC. erfc, exp and quadrature: transcendental only.
 *
 * Reference: Y. Liu, W. Whitt, Y. Yu (2016). Approximations for heavily-loaded
 * G/GI/n+GI queues. Naval Research Logistics 63(3), 187-217.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Steady-state measures of the G/GI/n+GI truncated Gaussian approximation. */
template <class T>
struct QsysTgaResult {
    std::string regime;        ///< "underloaded" or "overloaded"
    T trafficIntensity;        ///< rho = lambda/(n mu)
    T fluidWait;               ///< w, the fluid waiting time
    T fluidQueueLength;        ///< the fluid queue content
    T meanWait;                ///< E[W] after truncation
    T varWait;                 ///< Var[W]
    T meanQueueLength;         ///< E[Q] after truncation
    T varQueueLength;          ///< Var[Q]
    T meanNumberInService;     ///< E[B]
    T meanNumber;              ///< E[X] = E[B] + E[Q]
    T probDelay;               ///< P(W > 0)
    T probAbandon;             ///< P(patience < wait)
    T sigmaW;                  ///< the CLT scale of the wait
    T sigmaX;                  ///< the CLT scale of the content
};

namespace detail {

/** Standard normal density. */
template <class T>
T tga_phi(const T& x) {
    using std::exp;
    using std::sqrt;
    return exp(-x * x / num_traits<T>::from_int(2)) /
           sqrt(num_traits<T>::from_int(2) * num_traits<T>::from_double(M_PI));
}

/** Standard normal cdf, through erfc. */
template <class T>
T tga_Phi(const T& x) {
    using std::erfc;
    using std::sqrt;
    return erfc(-x / sqrt(num_traits<T>::from_int(2))) / num_traits<T>::from_int(2);
}

/** Composite Simpson rule on a fixed even panel count. */
template <class T, class Fn>
T tga_simpson(Fn&& f, const T& a, const T& b, std::size_t m = 2000) {
    const T zero = num_traits<T>::from_int(0);
    if (b <= a) return zero;
    const T h = (b - a) / num_traits<T>::from_int(static_cast<long>(m));
    T sum = f(a) + f(b);
    for (std::size_t i = 1; i < m; ++i)
        sum += num_traits<T>::from_int(i % 2 == 1 ? 4 : 2) *
               f(T(a + num_traits<T>::from_int(static_cast<long>(i)) * h));
    return h / num_traits<T>::from_int(3) * sum;
}

/** Smallest w with F^c(w) = target, by doubling then bisection. */
template <class T, class Ccdf>
T tga_inv_ccdf(Ccdf&& ccdf, const T& target) {
    const T two = num_traits<T>::from_int(2);
    T lo = num_traits<T>::from_int(0), hi = num_traits<T>::from_int(1);
    while (ccdf(hi) > target) {
        hi *= two;
        if (hi > num_traits<T>::from_double(1e12))
            throw InputError("qsys_ggingi_tga: the patience ccdf never falls to 1/rho, so the "
                             "overloaded model has no fluid equilibrium");
    }
    const T one = num_traits<T>::from_int(1);
    const T tol = num_traits<T>::from_double(1e-12);
    while (hi - lo > tol * (hi > one ? hi : one)) {
        const T mid = (lo + hi) / two;
        if (ccdf(mid) > target) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return (lo + hi) / two;
}

/** Mean and variance of max(Z,-a) for a standard normal Z. */
template <class T>
void tga_trunc_moments(const T& a, T& m1, T& v) {
    const T one = num_traits<T>::from_int(1);
    const T Pa = tga_Phi(a);
    const T pa = tga_phi(a);
    m1 = pa - a * (one - Pa);
    const T m2 = Pa - a * pa + a * a * (one - Pa);
    v = m2 - m1 * m1;
    if (v < num_traits<T>::from_int(0)) v = num_traits<T>::from_int(0);
}

}  // namespace detail

/**
 * @param lambda       arrival rate
 * @param mu           service rate of one server
 * @param n            number of servers
 * @param ca           coefficient of variation of the interarrival time
 * @param cs           coefficient of variation of the service time
 * @param patienceCcdf F^c(x) = P(patience > x)
 * @param patiencePdf  the patience density; empty differences the ccdf
 * @param serviceCcdf  G^c(x), used only in the underloaded branch
 */
template <class T>
QsysTgaResult<T> qsys_ggingi_tga(
    const T& lambda, const T& mu, unsigned n, const T& ca, const T& cs,
    const std::function<T(const T&)>& patienceCcdf,
    const std::function<T(const T&)>& patiencePdf = std::function<T(const T&)>(),
    const std::function<T(const T&)>& serviceCcdf = std::function<T(const T&)>()) {
    static_assert(num_traits<T>::has_transcendental, "qsys_ggingi_tga needs erfc and exp");
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (lambda <= zero || mu <= zero)
        throw InputError("qsys_ggingi_tga: the arrival and service rates must be positive");
    if (n < 1) throw InputError("qsys_ggingi_tga: the number of servers n must be at least 1");

    const T nT = num_traits<T>::from_int(static_cast<long>(n));
    const T ca2 = ca * ca;
    const T rho = lambda / (nT * mu);
    const T lamPn = lambda / nT;

    std::function<T(const T&)> pdf = patiencePdf;
    if (!pdf) {
        pdf = [&patienceCcdf](const T& x) {
            const T h = num_traits<T>::from_double(1e-6);
            const T lo = x - h < num_traits<T>::from_int(0) ? num_traits<T>::from_int(0) : T(x - h);
            const T d = (patienceCcdf(lo) - patienceCcdf(T(x + h))) / (num_traits<T>::from_int(2) * h);
            return d < num_traits<T>::from_int(0) ? num_traits<T>::from_int(0) : d;
        };
    }

    QsysTgaResult<T> r;
    r.trafficIntensity = rho;
    if (rho <= one) {
        // Underloaded: no queue in the limit; the content is normal with the
        // infinite-server variance (eq. 10).
        T omega = num_traits<T>::from_rational(1, 2);
        if (serviceCcdf) {
            const T hi = detail::tga_inv_ccdf<T>(serviceCcdf, num_traits<T>::from_double(1e-12));
            omega = detail::tga_simpson<T>([&](const T& x) { return T(serviceCcdf(x) * serviceCcdf(x)); },
                                           zero, hi) * mu;
        }
        r.regime = "underloaded";
        r.fluidWait = r.fluidQueueLength = r.meanWait = r.varWait = zero;
        r.meanQueueLength = r.varQueueLength = r.probDelay = r.probAbandon = zero;
        r.meanNumberInService = r.meanNumber = lambda / mu;
        r.sigmaW = zero;
        r.sigmaX = sqrt((lambda / mu) * (one + (ca2 - one) * omega));
        return r;
    }

    // Overloaded: the fluid centre of Theorem 2.1(b).
    const T w = detail::tga_inv_ccdf<T>(patienceCcdf, T(one / rho));
    const T fw = pdf(w);
    if (fw <= zero)
        throw InputError("qsys_ggingi_tga: the patience density vanishes at the fluid waiting "
                         "time, so the Gaussian correction is undefined there");
    const T qPerServer = lamPn * detail::tga_simpson<T>(patienceCcdf, zero, w);

    // Eq. (24): the service law enters only through the (cs+1)rho term.
    const T sigmaW2 = ((ca2 - one) + (cs + one) * rho) / (two * mu * rho * rho * fw);
    const T sigmaX2 =
        mu * mu * sigmaW2 +
        lamPn * detail::tga_simpson<T>(
                    [&](const T& x) { return T(patienceCcdf(x) * (one + (ca2 - one) * patienceCcdf(x))); },
                    zero, w);
    r.sigmaW = sqrt(sigmaW2 < zero ? zero : sigmaW2);
    r.sigmaX = sqrt(sigmaX2 < zero ? zero : sigmaX2);

    const T aW = sqrt(nT) * w / r.sigmaW;              // eq. (21)
    const T aX = sqrt(nT) * qPerServer / r.sigmaX;     // eq. (19)
    T m1W, vW, m1X, vX;
    detail::tga_trunc_moments(aW, m1W, vW);
    detail::tga_trunc_moments(aX, m1X, vX);

    r.regime = "overloaded";
    r.fluidWait = w;
    r.fluidQueueLength = nT * qPerServer;
    r.meanWait = w * (detail::tga_Phi(aW) + detail::tga_phi(aW) / aW);
    r.varWait = (r.sigmaW * r.sigmaW / nT) * vW;
    r.meanQueueLength = nT * qPerServer * (detail::tga_Phi(aX) + detail::tga_phi(aX) / aX);
    r.varQueueLength = nT * r.sigmaX * r.sigmaX * vX;
    // E[B] = E[min(X_n,n)]: every server is busy but for the lower tail.
    r.meanNumberInService =
        nT - sqrt(nT) * r.sigmaX * (detail::tga_phi(aX) - aX * (one - detail::tga_Phi(aX)));
    r.meanNumber = r.meanNumberInService + r.meanQueueLength;
    r.probDelay = detail::tga_Phi(aW);                  // eq. (22)
    // Eq. (23): a customer abandons when its patience falls short of its wait.
    const T hi = w * num_traits<T>::from_int(20) > w + num_traits<T>::from_int(20)
                     ? T(w * num_traits<T>::from_int(20))
                     : T(w + num_traits<T>::from_int(20));
    T pa = detail::tga_simpson<T>(
        [&](const T& x) { return T((one - detail::tga_Phi(T(aW * (x / w - one)))) * pdf(x)); }, zero,
        hi);
    if (pa < zero) pa = zero;
    if (pa > one) pa = one;
    r.probAbandon = pa;
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GGINGI_TGA_H
