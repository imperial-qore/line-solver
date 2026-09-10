/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_SRV_EXP_H
#define LINE_API_SNC_SRV_EXP_H

/**
 * MGF service envelope of an exponential server, in JOB units.
 *
 * A single server with Exp(mu) service times completes jobs at the epochs of a
 * Poisson process of rate mu while it is busy, so its cumulative service counted
 * in JOBS is Poisson with mean mu*(t-s) and
 * `rho(theta) = mu*(1-exp(-theta))/theta` with a zero burst.
 *
 * THIS IS THE SERVICE ELEMENT TO USE WHENEVER THE WORK UNIT IS THE JOB. Pairing
 * `snc_srv_rate` with a job-counting arrival envelope would model a server that
 * completes jobs at deterministic intervals, an M/D/1, and would UNDERSTATE the
 * delay of an exponential server rather than bound it. The M/M/1 read with this
 * element reproduces both exact decay rates: the backlog bound decays as
 * (lambda/mu)^n in jobs and the delay bound as exp(-(mu-lambda)*d) in time,
 * since the optimal theta tends to log(mu/lambda).
 *
 * Job units also compose across hops: a departure envelope from `snc_output` is
 * a job count and is directly the arrival envelope of the next station, whereas
 * service-time work units differ from station to station.
 *
 * Port of matlab/src/api/snc/snc_srv_exp.m.
 */

#include <cmath>
#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

/**
 * @param mu    service rate, jobs per slot
 * @param theta Chernoff parameter, theta > 0
 */
inline Env snc_srv_exp(double mu, double theta) {
    if (mu <= 0) throw UnsupportedError("snc_srv_exp: mu must be positive");
    if (theta <= 0) throw UnsupportedError("snc_srv_exp: theta must be positive");
    return Env{0.0, mu * (1.0 - std::exp(-theta)) / theta};
}

/** The same envelope as a function of theta. */
inline Envelope snc_srv_exp_fn(double mu) {
    return [mu](double theta) { return snc_srv_exp(mu, theta); };
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_SRV_EXP_H
