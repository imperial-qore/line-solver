/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_APPROX_KOBAYASHI_H
#define LINE_API_QSYS_GIG1_APPROX_KOBAYASHI_H

/**
 * Kobayashi diffusion approximation for the G/I/G/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_approx_kobayashi.m,
 * cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gig1_approx_kobayashi.java
 * (identical).
 *
 *   rhohat = exp(-2(1-rho)/(rho(ca^2 + cs^2/rho)))
 *   W      = rhohat/(1-rhohat)/lambda
 *
 * The geometric fit carries an exp, so this requires transcendental
 * arithmetic and cannot be instantiated at T = Rational.
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
QsysResult<T> qsys_gig1_approx_kobayashi(const T& lambda, const T& mu, const T& ca, const T& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gig1_approx_kobayashi requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    const T rhohat =
        detail::num_exp(T(-two * (one - rho) / (rho * (num_pow_int(ca, 2) + num_pow_int(cs, 2) / rho))));
    const T W = rhohat / (one - rhohat) / lambda;
    return {W, rhohat};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_APPROX_KOBAYASHI_H
