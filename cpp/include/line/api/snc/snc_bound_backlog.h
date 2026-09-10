/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_BOUND_BACKLOG_H
#define LINE_API_SNC_BOUND_BACKLOG_H

/**
 * Violation probability of a backlog level.
 *
 * For a flow with arrival envelope (sigmaA,rhoA) served by an element with
 * service envelope (sigmaS,rhoS), the backlog of the stable station obeys, for
 * every theta > 0,
 *
 *   P{B(t) > b} <= exp(-theta*(b-sigmaA-sigmaS)) / (1-exp(-theta*(rhoS-rhoA))),
 *
 * the union bound over the start of the backlogged period summed as a geometric
 * series on the unit-slot time axis. The returned value is the infimum over
 * theta, clipped at 1, and is an UPPER BOUND on the tail, never an estimate of
 * it: the decay rate is asymptotically exact and the prefactor is loose.
 *
 * Port of matlab/src/api/snc/snc_bound_backlog.m.
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
 * @param b        backlog level, units of the envelopes
 * @param thetamax upper end of the theta search
 */
inline SncResult snc_bound_backlog(const Envelope& arv, const Envelope& srv, double b,
                                   double thetamax = 1e3) {
    if (b < 0) throw UnsupportedError("snc_bound_backlog: b must be nonnegative");
    const SncResult r = snc_thetaopt(
        [&](double theta) {
            const detail::SncPair p = detail::snc_pair(arv, srv, theta);
            if (!p.ok) return std::numeric_limits<double>::infinity();
            return std::exp(-theta * (b - p.a.sigma - p.s.sigma)) /
                   (1.0 - std::exp(-theta * (p.s.rho - p.a.rho)));
        },
        thetamax);
    // no feasible theta, or the bound is vacuous at this level
    const double eps = (!std::isfinite(r.value) || r.value > 1.0) ? 1.0 : r.value;
    return SncResult{eps, r.theta};
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_BOUND_BACKLOG_H
