/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIGK_APPROX_H
#define LINE_API_QSYS_GIGK_APPROX_H

/**
 * Default G/I/G/k approximation of the mean response time.
 *
 * Templated port of matlab/src/api/qsys/qsys_gigk_approx.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_gigk_approx.java (identical;
 * the JAR is careful to write (k+1)/2.0 so the exponent stays a real, as in
 * MATLAB).
 *
 *   rho   = lambda/(mu k)
 *   alpha = (rho^k+rho)/2        if rho > 0.7
 *   alpha = rho^((k+1)/2)        otherwise
 *   W     = (alpha/mu)(1/(1-rho))(ca^2+cs^2)/(2k) + 1/mu
 *
 * The low-load branch raises rho to a half-integer power whenever k is even,
 * which is a genuine transcendental, so the function requires
 * num_traits<T>::has_transcendental and cannot be instantiated at exact
 * arithmetic. Gating is at function level rather than per branch, since the
 * branch is chosen at run time.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/**
 * @param lambda arrival rate
 * @param mu     service rate of a single server
 * @param ca     coefficient of variation of the interarrival time
 * @param cs     coefficient of variation of the service time
 * @param k      number of servers, k >= 1
 */
template <class T>
QsysResult<T> qsys_gigk_approx(const T& lambda, const T& mu, const T& ca, const T& cs, unsigned k) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gigk_approx requires transcendental arithmetic");
    if (k == 0) throw InputError("qsys_gigk_approx: k must be at least 1");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T kT = num_traits<T>::from_int(static_cast<long>(k));
    const T rho = lambda / (mu * kT);
    detail::require_no_pole(T(one - rho), "qsys_gigk_approx");
    T alpha;
    if (rho > num_traits<T>::from_rational(7, 10)) {
        alpha = (num_pow_int(rho, k) + rho) / two;
    } else {
        alpha = detail::num_pow(rho, num_traits<T>::from_rational(static_cast<long>(k) + 1, 2));
    }
    const T W = (alpha / mu) * (one / (one - rho)) * (num_pow_int(ca, 2) + num_pow_int(cs, 2)) /
                    (two * kT) +
                one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIGK_APPROX_H
