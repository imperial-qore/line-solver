/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_MEAN_BACKLOG_H
#define LINE_API_SNC_MEAN_BACKLOG_H

/**
 * Upper bound on the mean backlog, from integrating the backlog tail bound.
 *
 * The counterpart of `snc_mean_delay` with `a = theta`: the clipped integral of
 * `K*exp(-theta*b)` is `(log(K)+1)/theta` when K >= 1 and `K/theta` otherwise.
 * The unit of the answer is the unit of the envelopes: jobs when the pair is
 * `snc_env_poisson` with `snc_srv_exp`, units of work when it is
 * `snc_env_cpoisson` with `snc_srv_rate`.
 *
 * SolverBA does NOT use this for its queue-length column: it applies Little's
 * law to the response-time bound instead, so that Q and R stay consistent with
 * the exact open-network throughput. The two are close but not identical, since
 * each optimizes its own theta.
 *
 * Port of matlab/src/api/snc/snc_mean_backlog.m.
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
inline SncResult snc_mean_backlog(const Envelope& arv, const Envelope& srv, double thetamax = 1e3) {
    return snc_thetaopt(
        [&](double theta) {
            const detail::SncPair p = detail::snc_pair(arv, srv, theta);
            if (!p.ok) return std::numeric_limits<double>::infinity();
            const double logK = theta * (p.a.sigma + p.s.sigma) -
                                std::log(1.0 - std::exp(-theta * (p.s.rho - p.a.rho)));
            return logK >= 0.0 ? (logK + 1.0) / theta : std::exp(logK) / theta;
        },
        thetamax);
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_MEAN_BACKLOG_H
