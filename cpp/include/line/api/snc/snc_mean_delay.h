/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_MEAN_DELAY_H
#define LINE_API_SNC_MEAN_DELAY_H

/**
 * Upper bound on the mean delay, from integrating the delay tail bound.
 *
 * For a nonnegative delay `E[D] = int_0^inf P{D>d} dd`, so integrating the tail
 * bound of `snc_bound_delay` bounds the MEAN. At fixed theta the bound is
 * `K*exp(-a*d)` with `a = theta*rhoS` and
 * `K = exp(theta*(sigmaA+sigmaS))/(1-exp(-theta*(rhoS-rhoA)))`, so, clipping the
 * bound at 1 where it exceeds it, the integral is available in CLOSED FORM:
 * `(log(K)+1)/a` when K >= 1 and `K/a` otherwise. No quadrature is involved, so
 * the result is a bound and not a bound plus a discretization error.
 *
 * IT IS A LOOSE MEAN BOUND AND THAT IS INHERENT: on the M/M/1 read in job units
 * it returns 2.4x the exact 1/(mu-lambda) at rho = 0.1 and 10.4x at rho = 0.95,
 * because the prefactor of the tail bound, not its decay rate, dominates an
 * integral over the whole axis. Use `snc_perc_delay` when the quantile is what
 * matters.
 *
 * Port of matlab/src/api/snc/snc_mean_delay.m. This is what the SolverBA
 * `snc.upper` response-time column calls.
 */

#include <cmath>

#include "line/api/snc/snc_thetaopt.h"
#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

/**
 * @param arv      arrival envelope
 * @param srv      service envelope
 * @param thetamax upper end of the theta search
 */
inline SncResult snc_mean_delay(const Envelope& arv, const Envelope& srv, double thetamax = 1e3) {
    return snc_thetaopt(
        [&](double theta) {
            const detail::SncPair p = detail::snc_pair(arv, srv, theta);
            if (!p.ok) return std::numeric_limits<double>::infinity();
            const double logK = theta * (p.a.sigma + p.s.sigma) -
                                std::log(1.0 - std::exp(-theta * (p.s.rho - p.a.rho)));
            const double a = theta * p.s.rho;
            return logK >= 0.0 ? (logK + 1.0) / a : std::exp(logK) / a;
        },
        thetamax);
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_MEAN_DELAY_H
