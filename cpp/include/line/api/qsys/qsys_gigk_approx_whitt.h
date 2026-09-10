/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIGK_APPROX_WHITT_H
#define LINE_API_QSYS_GIGK_APPROX_WHITT_H

/**
 * Whitt (1993) approximation for the GI/G/k queue, eqs. (2.16)-(2.25).
 *
 * Templated port of matlab/src/api/qsys/qsys_gigk_approx_whitt.m.
 *
 *   gamma = min(0.24, (1-rho)(k-1)(sqrt(4+5k)-2)/(16 k rho))   (2.17)
 *   phi1  = 1 + gamma                                          (2.16)
 *   phi2  = 1 - 4 gamma                                        (2.18)
 *   phi3  = phi2 exp(-2(1-rho)/(3 rho))                        (2.20)
 *   phi4  = min(1, (phi1+phi3)/2)                              (2.21)
 *   psi   = 1 if c2 >= 1 else phi4^(2(1-c2))                   (2.22)
 *   phi   = psi                                    if ca^2 == cs^2
 *         = 4(ca^2-cs^2)/(4ca^2-3cs^2) phi1 + cs^2/(4ca^2-3cs^2) psi
 *                                                  if ca^2 > cs^2  (2.25)
 *         = (cs^2-ca^2)/(2(ca^2+cs^2)) phi3 + (cs^2+3ca^2)/(2(ca^2+cs^2)) psi
 *                                                  otherwise
 *   Wq    = phi c2 Wq(M/M/k),  c2 = (ca^2+cs^2)/2          (2.24)
 *
 * The equality branch is taken when |ca^2-cs^2| < 1e-12, exactly as in MATLAB
 * and the JAR, so the port keeps that literal tolerance.
 *
 * DIVERGENCE: jar/.../Qsys_gigk_approx_whitt.java names its arguments ca2 and
 * cs2 and uses them unsquared, i.e. it expects squared coefficients of
 * variation, whereas MATLAB takes ca, cs and squares them internally. MATLAB
 * is ground truth, so this port takes ca, cs. The JAR also returns {L,W,Q,U}
 * instead of [W,rhohat].
 *
 * Reference: Whitt, W. (1993). Approximations for the GI/G/m queue.
 * Production and Operations Management 2(2), 114-161.
 *
 * Carries sqrt, exp and a real-valued power, so it requires transcendental
 * arithmetic and cannot be instantiated at T = Rational.
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
QsysResult<T> qsys_gigk_approx_whitt(const T& lambda, const T& mu, const T& ca, const T& cs,
                                     unsigned k) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gigk_approx_whitt requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T four = num_traits<T>::from_int(4);
    const T kT = num_traits<T>::from_int(static_cast<long>(k));
    const T ca2 = num_pow_int(ca, 2);
    const T cs2 = num_pow_int(cs, 2);
    const T rho = lambda / (kT * mu);

    const T Wq_mmk = qsys_mmk(lambda, mu, k).W - one / mu;

    const T gamma = detail::num_min(
        num_traits<T>::from_rational(6, 25),
        (one - rho) * (kT - one) *
            (detail::num_sqrt(T(four + num_traits<T>::from_int(5) * kT)) - two) /
            (num_traits<T>::from_int(16) * kT * rho));
    const T phi1 = one + gamma;
    const T phi2 = one - four * gamma;
    const T phi3 = phi2 * detail::num_exp(T(-two * (one - rho) / (three * rho)));
    const T phi4 = detail::num_min(one, T((phi1 + phi3) / two));

    const T c2 = (ca2 + cs2) / two;
    T psi;
    if (c2 >= one) {
        psi = one;
    } else {
        psi = detail::num_pow(phi4, T(two * (one - c2)));
    }

    T phi;
    if (num_abs(T(ca2 - cs2)) < num_traits<T>::from_double(1e-12)) {
        phi = psi;
    } else if (ca2 > cs2) {
        phi = (four * (ca2 - cs2) / (four * ca2 - three * cs2)) * phi1 +
              (cs2 / (four * ca2 - three * cs2)) * psi;
    } else {
        phi = ((cs2 - ca2) / (two * (ca2 + cs2))) * phi3 +
              ((cs2 + three * ca2) / (two * (ca2 + cs2))) * psi;
    }

    const T Wq = phi * c2 * Wq_mmk;
    const T W = Wq + one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIGK_APPROX_WHITT_H
