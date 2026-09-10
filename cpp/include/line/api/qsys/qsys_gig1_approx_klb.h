/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_APPROX_KLB_H
#define LINE_API_QSYS_GIG1_APPROX_KLB_H

/**
 * Kraemer and Langenbach-Belz approximation for the G/I/G/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_approx_klb.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_gig1_approx_klb.java
 * (identical).
 *
 *   g = exp(-2(1-rho)(1-ca^2)^2/(3 rho (ca^2+cs^2)))      if ca <= 1
 *   g = exp(-(1-rho)(ca^2-1)/(ca^2+4 cs^2))               otherwise
 *   W = 1/mu * ((rho/(1-rho)) (cs^2+ca^2)/2 g + 1)
 *
 * The correction factor g is an exponential, so the function requires
 * transcendental arithmetic and cannot be instantiated at T = Rational, even
 * though g degenerates to exactly 1 at ca = 1.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"

namespace line {
namespace qsys {

/**
 * @param lambda arrival rate
 * @param mu     service rate
 * @param ca     coefficient of variation of the interarrival time
 * @param cs     coefficient of variation of the service time
 */
template <class T>
QsysResult<T> qsys_gig1_approx_klb(const T& lambda, const T& mu, const T& ca, const T& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gig1_approx_klb requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T four = num_traits<T>::from_int(4);
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_gig1_approx_klb");
    const T ca2 = num_pow_int(ca, 2);
    const T cs2 = num_pow_int(cs, 2);
    T g;
    if (ca <= one) {
        g = detail::num_exp(T(-two * (one - rho) * num_pow_int(one - ca2, 2) / (three * rho * (ca2 + cs2))));
    } else {
        g = detail::num_exp(T(-(one - rho) * (ca2 - one) / (ca2 + four * cs2)));
    }
    const T W = one / mu * ((rho / (one - rho)) * ((cs2 + ca2) / two) * g + one);
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_APPROX_KLB_H
