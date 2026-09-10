/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_BOUND_DELAY_H
#define LINE_API_SNC_BOUND_DELAY_H

/**
 * Violation probability of a delay target.
 *
 * For a flow with arrival envelope (sigmaA,rhoA) served by an element with
 * service envelope (sigmaS,rhoS), the virtual delay of the stable station obeys,
 * for every theta > 0,
 *
 *   P{D(t) > d} <= exp(-theta*(rhoS*d-sigmaA-sigmaS)) / (1-exp(-theta*(rhoS-rhoA))),
 *
 * the horizontal rather than vertical deviation between the arrival and service
 * envelopes. On the M/M/1 read in job units (`snc_env_poisson` with
 * `snc_srv_exp`) the optimal theta tends to log(mu/lambda), so the bound
 * reproduces the exact asymptotic decay rate exp(-(mu-lambda)*d).
 *
 * Port of matlab/src/api/snc/snc_bound_delay.m.
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
 * @param d        delay target, slots
 * @param thetamax upper end of the theta search
 */
inline SncResult snc_bound_delay(const Envelope& arv, const Envelope& srv, double d,
                                 double thetamax = 1e3) {
    if (d < 0) throw UnsupportedError("snc_bound_delay: d must be nonnegative");
    const SncResult r = snc_thetaopt(
        [&](double theta) {
            const detail::SncPair p = detail::snc_pair(arv, srv, theta);
            if (!p.ok) return std::numeric_limits<double>::infinity();
            return std::exp(-theta * (p.s.rho * d - p.a.sigma - p.s.sigma)) /
                   (1.0 - std::exp(-theta * (p.s.rho - p.a.rho)));
        },
        thetamax);
    const double eps = (!std::isfinite(r.value) || r.value > 1.0) ? 1.0 : r.value;
    return SncResult{eps, r.theta};
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_BOUND_DELAY_H
