/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_GIG1_APPROX_MYSKJA2_H
#define LINE_API_QSYS_QSYS_GIG1_APPROX_MYSKJA2_H

/**
 * Myskja's enhanced third-moment approximation of the mean response time of a
 * G/I/G/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_approx_myskja2.m.
 *
 *   ra = (1+ca^2)/2,  rs = (1+cs^2)/2,  rho = lambda/mu
 *   theta = [ rho(qa-ra) - (qa-ra^2) ] / [ 2 rho (ra-1) ]
 *   d     = (1+1/ra)(1-rs)(1-(q0/qa)^3)(1-rho^3)
 *   D     = (rs-theta)^2 + (2 rs - 1 + d)(ra-1),  clamped at 0
 *   W     = (rho/(1-rho))/lambda [ rs + (1/rho)( sqrt(D) - (rs-theta) ) ]
 *
 * At ca = 1 the interpolation parameter theta is a 0/0 form, so the exact
 * M/G/1 answer is returned instead; that is also the anchor the method
 * interpolates from.
 *
 * ARITHMETIC. sqrt(D) and the cube of q0/qa make this transcendental (the cube
 * alone would not, but the square root does), so the function is gated.
 *
 * MATLAB-vs-JAR. jline.api.qsys.Qsys_gig1_approx_myskja2 sets
 * ra = (1+ca)/2, rs = (1+cs)/2 and branches on |ca-1| < 1e-8, i.e. it reads
 * ca/cs as squared coefficients of variation, and in the M/G/1 branch it calls
 * qsys_mg1 with sqrt(cs) rather than cs. MATLAB reads them as coefficients of
 * variation throughout. The two disagree on every non-Markovian input; this
 * port follows MATLAB.
 */

#include "line/api/qsys/qsys_mg1.h"
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
QsysResult<T> qsys_gig1_approx_myskja2(const T& lambda, const T& mu, const T& ca, const T& cs,
                                       const T& q0, const T& qa) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gig1_approx_myskja2 requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (num_abs(T(ca * ca - one)) < T(num_traits<T>::from_double(1e-8)))
        return qsys_mg1(lambda, mu, cs);  // M/G/1 case: exact

    const T ra = (one + ca * ca) / two;
    const T rs = (one + cs * cs) / two;
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_gig1_approx_myskja2");
    const T theta = (rho * (qa - ra) - (qa - ra * ra)) / (two * rho * (ra - one));
    const T d = (one + one / ra) * (one - rs) * (one - num_pow_int(T(q0 / qa), 3)) *
                (one - num_pow_int(rho, 3));
    T D = num_pow_int(T(rs - theta), 2) + (two * rs - one + d) * (ra - one);
    if (D < zero) D = zero;  // guard small negative values due to round-off
    const T W = (rho / (one - rho)) / lambda *
                (rs + (one / rho) * (detail::num_sqrt(D) - (rs - theta)));
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_GIG1_APPROX_MYSKJA2_H
