/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_UBND_KINGMAN_H
#define LINE_API_QSYS_GIG1_UBND_KINGMAN_H

/**
 * Kingman upper bound on the mean waiting time of a G/G/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_ubnd_kingman.m,
 * cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gig1_ubnd_kingman.java (identical).
 *
 *   Wq <= lambda (sa^2 + ss^2)/(2(1-rho)),  sa^2 = ca^2/lambda^2, ss^2 = cs^2/mu^2
 *   W  = Wq + 1/mu
 *
 * Reference: Kingman, J.F.C. (1962). Some inequalities for the queue GI/G/1.
 * Biometrika 49(3/4), 315-324.
 *
 * Pure field arithmetic, exact for T = Rational.
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
QsysResult<T> qsys_gig1_ubnd_kingman(const T& lambda, const T& mu, const T& ca, const T& cs) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_gig1_ubnd_kingman");
    const T Wq = lambda *
                 (num_pow_int(ca, 2) / num_pow_int(lambda, 2) + num_pow_int(cs, 2) / num_pow_int(mu, 2)) /
                 (two * (one - rho));
    const T W = Wq + one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_UBND_KINGMAN_H
