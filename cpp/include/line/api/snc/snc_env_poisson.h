/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_ENV_POISSON_H
#define LINE_API_SNC_ENV_POISSON_H

/**
 * MGF arrival envelope of a Poisson flow with unit-size jobs.
 *
 * `log E[exp(theta*A(0,t))] = lambda*t*(exp(theta)-1)` exactly, so the envelope
 * is tight with a zero burst term and
 * `rho(theta) = lambda*(exp(theta)-1)/theta`.
 *
 * Port of matlab/src/api/snc/snc_env_poisson.m.
 */

#include <cmath>
#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

/**
 * @param lambda arrival rate, jobs per slot
 * @param theta  Chernoff parameter, theta > 0
 */
inline Env snc_env_poisson(double lambda, double theta) {
    if (lambda < 0)
        throw UnsupportedError("snc_env_poisson: lambda must be nonnegative");
    if (theta <= 0) throw UnsupportedError("snc_env_poisson: theta must be positive");
    return Env{0.0, lambda * (std::exp(theta) - 1.0) / theta};
}

/** The same envelope as a function of theta. */
inline Envelope snc_env_poisson_fn(double lambda) {
    return [lambda](double theta) { return snc_env_poisson(lambda, theta); };
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_ENV_POISSON_H
