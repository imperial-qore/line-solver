/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIGK_APPROX_COSMETATOS_H
#define LINE_API_QSYS_GIGK_APPROX_COSMETATOS_H

/**
 * Cosmetatos / Page interpolation approximation for the GI/G/k queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gigk_approx_cosmetatos.m.
 *
 *   gamma = min(0.24, (1-rho)(k-1)(sqrt(4+5k)-2)/(16 k rho))
 *   phi1  = 1 + gamma                                  (M/D/k factor)
 *   phi3  = (1-4 gamma) exp(-2(1-rho)/(3 rho))         (D/M/k factor)
 *   Wq    = [ca^2 cs^2 + ca^2(1-cs^2) phi1/2 + (1-ca^2) cs^2 phi3/2] Wq(M/M/k)
 *
 * for ca^2 <= 1 and cs^2 <= 1; outside the unit box the Lee-Longton scaling
 * Wq = ((ca^2+cs^2)/2) Wq(M/M/k) is used instead. W = Wq + 1/mu.
 *
 * DIVERGENCE: jar/.../Qsys_gigk_approx_cosmetatos.java names its arguments
 * ca2 and cs2 and uses them unsquared, i.e. it expects squared coefficients of
 * variation, whereas MATLAB takes ca, cs and squares them internally. MATLAB
 * is ground truth, so this port takes ca, cs. The JAR also returns
 * {L,W,Q,U} instead of [W,rhohat]; W agrees once the argument convention is
 * matched.
 *
 * References: Cosmetatos (1975) INFOR 13, 328-331; Page (1982) J. Opl. Res.
 * Soc. 33, 453-473; Whitt (1993) eq. (2.17) for the gamma safeguard.
 *
 * Carries sqrt and exp, so it requires transcendental arithmetic and cannot
 * be instantiated at T = Rational.
 */

#include "line/api/qsys/qsys_mmk.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"

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
QsysResult<T> qsys_gigk_approx_cosmetatos(const T& lambda, const T& mu, const T& ca, const T& cs,
                                          unsigned k) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gigk_approx_cosmetatos requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T four = num_traits<T>::from_int(4);
    const T kT = num_traits<T>::from_int(static_cast<long>(k));
    const T ca2 = num_pow_int(ca, 2);
    const T cs2 = num_pow_int(cs, 2);
    const T rho = lambda / (kT * mu);

    const T Wq_mmk = qsys_mmk(lambda, mu, k).W - one / mu;

    T Wq;
    if (ca2 <= one && cs2 <= one) {
        const T gamma = detail::num_min(
            num_traits<T>::from_rational(6, 25),
            (one - rho) * (kT - one) *
                (detail::num_sqrt(T(four + num_traits<T>::from_int(5) * kT)) - two) /
                (num_traits<T>::from_int(16) * kT * rho));
        const T phi1 = one + gamma;
        const T phi3 = (one - four * gamma) * detail::num_exp(T(-two * (one - rho) / (three * rho)));
        Wq = (ca2 * cs2 + ca2 * (one - cs2) * phi1 / two + (one - ca2) * cs2 * phi3 / two) * Wq_mmk;
    } else {
        Wq = ((ca2 + cs2) / two) * Wq_mmk;
    }
    const T W = Wq + one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIGK_APPROX_COSMETATOS_H
