/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_APPROX_ALLENCUNNEEN_H
#define LINE_API_QSYS_GIG1_APPROX_ALLENCUNNEEN_H

/**
 * Allen-Cunneen approximation of the mean response time of a G/I/G/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_approx_allencunneen.m,
 * cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gig1_approx_allencunneen.java
 * (identical).
 *
 *   W = (rho/(1-rho))/mu * (cs^2+ca^2)/2 + 1/mu
 *
 * Only squares of the coefficients of variation appear, so this is pure field
 * arithmetic and exact for T = Rational. It reduces to M/M/1 at ca = cs = 1.
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
QsysResult<T> qsys_gig1_approx_allencunneen(const T& lambda, const T& mu, const T& ca,
                                            const T& cs) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_gig1_approx_allencunneen");
    const T W = (rho / (one - rho)) / mu * ((num_pow_int(cs, 2) + num_pow_int(ca, 2)) / two) +
                one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_APPROX_ALLENCUNNEEN_H
