/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_GIG1_APPROX_MYSKJA_H
#define LINE_API_QSYS_QSYS_GIG1_APPROX_MYSKJA_H

/**
 * Myskja's third-moment approximation of the mean response time of a G/I/G/1
 * queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_approx_myskja.m.
 *
 *   Wq = rho/(2 mu (1-rho)) [ (1+cs^2) + (q0/qa)^(1/rho-rho) (1/rho) (ca^2-1) ]
 *   W  = Wq + 1/mu
 *
 * exact for M/G/1 (ca = 1), where the bracket collapses to 1 + cs^2. Here qa
 * is the third relative moment E[A^3]/(6 E[A]^3) of the interarrival time and
 * q0 its smallest value for the given mean and SCV, so q0/qa <= 1 and the
 * correction term is a genuine interpolation in the third moment.
 *
 * ARITHMETIC. The exponent 1/rho - rho is a real number, so (q0/qa)^(1/rho-rho)
 * is a transcendental evaluation and the function is gated accordingly.
 *
 * MATLAB-vs-JAR. jline.api.qsys.Qsys_gig1_approx_myskja writes (1+cs) and
 * (ca-1) where MATLAB writes (1+cs^2) and (ca^2-1), i.e. the JAR reads its
 * ca/cs arguments as squared coefficients of variation while MATLAB reads them
 * as coefficients of variation. The two therefore disagree for every input
 * with ca != 1 or cs != 1. This port follows MATLAB, which is consistent with
 * the rest of the qsys_gig1_approx_* family and with the documented signature.
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
 * @param q0     smallest third relative moment for the given mean and SCV
 * @param qa     third relative moment E[A^3]/(6 E[A]^3) of the interarrival time
 */
template <class T>
QsysResult<T> qsys_gig1_approx_myskja(const T& lambda, const T& mu, const T& ca, const T& cs,
                                      const T& q0, const T& qa) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gig1_approx_myskja requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_gig1_approx_myskja");
    const T expo = one / rho - rho;
    const T Wq = rho / (two * mu * (one - rho)) *
                 ((one + cs * cs) +
                  detail::num_pow(T(q0 / qa), expo) * (one / rho) * (ca * ca - one));
    const T W = Wq + one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_GIG1_APPROX_MYSKJA_H
