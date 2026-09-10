/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIGK_APPROX_KINGMAN_H
#define LINE_API_QSYS_GIGK_APPROX_KINGMAN_H

/**
 * Kingman (Lee-Longton) scaling of the exact M/M/k waiting time.
 *
 * Templated port of matlab/src/api/qsys/qsys_gigk_approx_kingman.m,
 * cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gigk_approx_kingman.java. The numbers
 * agree; the JAR reads the M/M/k answer back out of the static fields of
 * Ret.qsys, which this port replaces by a plain return value (the port has no
 * global mutable state).
 *
 *   W = (ca^2+cs^2)/2 * (W_MMk - 1/mu) + 1/mu
 *
 * The M/M/k baseline is Erlang-C, all integer powers, so the whole function
 * is pure field arithmetic and exact for T = Rational. At k = 1 it reduces
 * exactly to the Allen-Cunneen / Heyman G/I/G/1 formula.
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
QsysResult<T> qsys_gigk_approx_kingman(const T& lambda, const T& mu, const T& ca, const T& cs,
                                       unsigned k) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T W_mmk = qsys_mmk(lambda, mu, k).W;
    const T W = (num_pow_int(ca, 2) + num_pow_int(cs, 2)) / two * (W_mmk - one / mu) + one / mu;
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIGK_APPROX_KINGMAN_H
