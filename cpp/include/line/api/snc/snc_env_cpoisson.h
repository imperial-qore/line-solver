/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_ENV_CPOISSON_H
#define LINE_API_SNC_ENV_CPOISSON_H

/**
 * MGF arrival envelope of a compound Poisson flow with Exp job sizes.
 *
 * Jobs arrive Poisson at rate lambda and each carries an Exp(mu) amount of work,
 * so `rho(theta) = lambda/(mu-theta)` for `0 < theta < mu` and the burst is
 * zero. Fed to a constant-rate server of rate mu (`snc_srv_rate`) this is the
 * network calculus model of the M/M/1 queue in units of WORK, and the delay
 * bound then decays at the exact rate mu-lambda.
 *
 * INFEASIBILITY IS A VALUE, NOT AN ERROR: at `theta >= mu` the job-size MGF
 * diverges and rho is returned as infinity, which the theta search discards.
 *
 * Port of matlab/src/api/snc/snc_env_cpoisson.m.
 */

#include <cmath>
#include <limits>
#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

/**
 * @param lambda job arrival rate, jobs per slot
 * @param mu     rate of the Exp job size, so the mean work per job is 1/mu
 * @param theta  Chernoff parameter, theta > 0
 */
inline Env snc_env_cpoisson(double lambda, double mu, double theta) {
    if (lambda < 0 || mu <= 0)
        throw UnsupportedError("snc_env_cpoisson: lambda must be nonnegative and mu positive");
    if (theta <= 0) throw UnsupportedError("snc_env_cpoisson: theta must be positive");
    if (theta >= mu) return Env{0.0, std::numeric_limits<double>::infinity()};
    return Env{0.0, lambda / (mu - theta)};
}

/** The same envelope as a function of theta. */
inline Envelope snc_env_cpoisson_fn(double lambda, double mu) {
    return [lambda, mu](double theta) { return snc_env_cpoisson(lambda, mu, theta); };
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_ENV_CPOISSON_H
