/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_APPROX_GELENBE_H
#define LINE_API_QSYS_GIG1_APPROX_GELENBE_H

/**
 * Gelenbe diffusion approximation with instantaneous-return boundary.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_approx_gelenbe.m.
 *
 *   rhat = exp(-2(1-rho)/(rho ca^2 + cs^2))
 *   W    = 1/(mu (1-rhat))
 *
 * DIVERGENCE: jar/src/main/java/jline/api/qsys/Qsys_gig1_approx_gelenbe.java
 * computes rhat = exp(-2(1-rho)/(rho*ca + cs)), treating its arguments as
 * already-squared coefficients of variation. MATLAB is ground truth and
 * squares, so this port squares. The two disagree unless ca = cs = 1.
 *
 * Reference: Gelenbe, E. (1975). On approximate computer system models.
 * Journal of the ACM 22(2), 261-269.
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
QsysResult<T> qsys_gig1_approx_gelenbe(const T& lambda, const T& mu, const T& ca, const T& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gig1_approx_gelenbe requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    const T rhat =
        detail::num_exp(T(-two * (one - rho) / (rho * num_pow_int(ca, 2) + num_pow_int(cs, 2))));
    const T W = one / (mu * (one - rhat));
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_APPROX_GELENBE_H
