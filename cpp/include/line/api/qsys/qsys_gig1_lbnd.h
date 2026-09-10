/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_LBND_H
#define LINE_API_QSYS_GIG1_LBND_H

/**
 * Fundamental lower bound on the mean response time of a G/G/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_lbnd.m.
 *
 *   W = 1/mu    (the response time is at least the mean service time)
 *
 * DIVERGENCE: jar/src/main/java/jline/api/qsys/Qsys_gig1_lbnd.java returns a
 * map {L=rho, Lq=0, W=1/mu, Wq=0, p0=1-rho} and no rhohat, while MATLAB
 * returns [W,rhohat]. The value of W agrees; the port follows the MATLAB
 * signature. Neither implementation reads ca or cs, which are accepted only
 * for signature uniformity with the rest of the G/I/G/1 family.
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
 * @param ca     coefficient of variation of the interarrival time, unused
 * @param cs     coefficient of variation of the service time, unused
 */
template <class T>
QsysResult<T> qsys_gig1_lbnd(const T& lambda, const T& mu, const T& ca, const T& cs) {
    (void)ca;
    (void)cs;
    const T W = num_traits<T>::from_int(1) / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_LBND_H
