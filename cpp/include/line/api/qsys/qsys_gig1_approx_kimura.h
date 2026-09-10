/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_APPROX_KIMURA_H
#define LINE_API_QSYS_GIG1_APPROX_KIMURA_H

/**
 * Kimura diffusion-interpolation approximation for the G/I/G/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_approx_kimura.m.
 *
 *   Wq = rho (ca^2+cs^2) / (mu (1-rho) (1+ca^2)),   W = Wq + 1/mu
 *
 * exact for M/M/1 and M/G/1.
 *
 * DIVERGENCE: jar/src/main/java/jline/api/qsys/Qsys_gig1_approx_kimura.java
 * computes Wq = rho*(ca+cs)/mu/(1-rho)/(1+ca), i.e. it treats its arguments as
 * already-squared coefficients of variation and never squares them. MATLAB is
 * ground truth and squares, so this port squares. The two implementations
 * disagree numerically whenever ca != 1 or cs != 1.
 *
 * Reference: Kimura, T. (1986). A two-moment approximation for the mean
 * waiting time in the GI/G/s queue. Management Science 32(6), 751-763.
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
QsysResult<T> qsys_gig1_approx_kimura(const T& lambda, const T& mu, const T& ca, const T& cs) {
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_gig1_approx_kimura");
    const T ca2 = num_pow_int(ca, 2);
    const T cs2 = num_pow_int(cs, 2);
    const T Wq = rho * (ca2 + cs2) / mu / (one - rho) / (one + ca2);
    const T W = Wq + one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_APPROX_KIMURA_H
