/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_APPROX_MARCHAL_H
#define LINE_API_QSYS_GIG1_APPROX_MARCHAL_H

/**
 * Marchal approximation of the mean response time of a G/I/G/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_approx_marchal.m,
 * cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gig1_approx_marchal.java (identical).
 *
 *   W = rho/(1-rho) * (1+cs^2)/2/mu * (ca + rho^2 cs^2)/(1 + rho^2 cs^2) + 1/mu
 *
 * The numerator carries ca and not ca^2. That is what both MATLAB and the JAR
 * compute, so the port reproduces it verbatim rather than "correcting" it to
 * the published form. Pure field arithmetic, exact for T = Rational; it still
 * reduces to M/M/1 at ca = cs = 1 because the correction ratio becomes one.
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
QsysResult<T> qsys_gig1_approx_marchal(const T& lambda, const T& mu, const T& ca, const T& cs) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_gig1_approx_marchal");
    const T Wmm1 = rho / (one - rho);
    const T r2cs2 = num_pow_int(rho, 2) * num_pow_int(cs, 2);
    const T W = Wmm1 * (one + num_pow_int(cs, 2)) / two / mu * (ca + r2cs2) / (one + r2cs2) + one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_APPROX_MARCHAL_H
