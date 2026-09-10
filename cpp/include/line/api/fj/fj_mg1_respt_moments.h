/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_FJ_MG1_RESPT_MOMENTS_H
#define LINE_API_FJ_FJ_MG1_RESPT_MOMENTS_H

/**
 * Mean and variance of the M/G/1 response time, as ForkTail inputs.
 *
 * Templated port of matlab/src/api/fj/fj_mg1_respt_moments.m. No JAR
 * counterpart. Closes the white-box route of fj_tail_forktail for a fork
 * branch that is an M/G/1 FCFS queue, from the first three moments of its
 * service time:
 *
 *   E[T] = E[S] (1 + rho/(1-rho) (1+SCV_S)/2)
 *   V[T] = E[W]^2 + lambda E[S^3]/(3(1-rho)) + E[S^2] - E[S]^2
 *
 * with rho = lambda E[S] and E[W] = lambda E[S^2]/(2(1-rho)), the
 * Pollaczek-Khinchine mean waiting time. The THIRD moment enters the variance
 * only, so a Markovian branch needs no input beyond what its service law
 * already reports.
 *
 * A service law with no finite third moment (a Pareto branch of shape <= 3,
 * say) leaves the ForkTail variance undefined; that is reported as an input
 * error rather than propagated as an infinity, because every downstream
 * ForkTail fit would then silently return the exponential special case.
 *
 * Reference: M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A
 * Black-Box Fork-Join Tail Latency Prediction Model for User-Facing
 * Datacenter Workloads", ACM HPDC 2018, equations (10) and (11).
 *
 * ARITHMETIC: rational in the four inputs, so the exact instantiation is
 * available and there is no static_assert.
 */

#include <cmath>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** Mirrors MATLAB's [ET, VT] return list. */
template <class T>
struct Mg1ResptMoments {
    T ET;  ///< mean response time
    T VT;  ///< variance of the response time
};

/**
 * @param lambda arrival rate at the branch
 * @param ES     first moment of the service time
 * @param ES2    second moment of the service time
 * @param ES3    third moment of the service time, finite
 */
template <class T>
Mg1ResptMoments<T> fj_mg1_respt_moments(const T& lambda, const T& ES, const T& ES2,
                                        const T& ES3) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T rho = lambda * ES;
    if (rho >= one) throw InputError("fj_mg1_respt_moments: the branch is unstable, rho >= 1");
    if (!std::isfinite(num_traits<T>::to_double(ES3)))
        throw InputError(
            "fj_mg1_respt_moments: the service law has no finite third moment, so the ForkTail "
            "response time variance is undefined");

    const T scvS = (ES2 - ES * ES) / (ES * ES);
    Mg1ResptMoments<T> r;
    r.ET = ES * (one + rho / (one - rho) * (one + scvS) / two);
    const T EW = lambda * ES2 / (two * (one - rho));
    r.VT = EW * EW + lambda * ES3 / (three * (one - rho)) + ES2 - ES * ES;
    return r;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_FJ_MG1_RESPT_MOMENTS_H
