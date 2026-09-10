/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RESPT_NOSPLIT_H
#define LINE_API_FJ_RESPT_NOSPLIT_H

/**
 * Mean response time of the distributed no-splitting parallel system.
 *
 * Templated port of matlab/src/api/fj/fj_respt_nosplit.m.
 *
 * A job of K tasks is routed in one piece to a single server chosen uniformly
 * among the K, so each server is an M/E_K/1 queue of arrival rate lambda/K and
 * service the sum of K exponential stages of rate mu. Pollaczek-Khinchine then
 * reduces to
 *
 *   R = [ K - (K-1) rho/2 ] / (mu - lambda),   rho = lambda/mu,
 *
 * the reference against which the splitting policies are judged. At K = 1 it
 * collapses to the M/M/1 response time.
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [R, rho] of fj_respt_nosplit. */
template <class T>
struct FJResptNosplitResult {
    T R;
    T rho;
};

/**
 * @param K      number of servers, equal to the number of tasks per job
 * @param lambda total job arrival rate
 * @param mu     per-server task service rate, mu > lambda
 * @return       the mean job response time and the per-server utilization
 */
template <class T>
FJResptNosplitResult<T> fj_respt_nosplit(unsigned K, const T& lambda, const T& mu) {
    detail::require_positive_K(K, "fj_respt_nosplit");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1),
            two = num_traits<T>::from_int(2);
    if (!(lambda > zero)) throw InputError("fj_respt_nosplit: the arrival rate must be positive");
    if (!(mu > zero)) throw InputError("fj_respt_nosplit: the service rate must be positive");

    FJResptNosplitResult<T> out;
    out.rho = lambda / mu;
    if (out.rho >= one)
        throw NumericError("fj_respt_nosplit: unstable system, rho = lambda/mu >= 1");
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    // M/M/1 response time at the same utilization
    const T Rmm1 = one / (mu - lambda);
    out.R = (Kt - (Kt - one) * out.rho / two) * Rmm1;
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RESPT_NOSPLIT_H
