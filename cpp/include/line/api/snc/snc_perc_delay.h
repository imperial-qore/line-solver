/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_PERC_DELAY_H
#define LINE_API_SNC_PERC_DELAY_H

/**
 * Delay quantile at a prescribed violation probability.
 *
 * Inverts `snc_bound_delay` in d: at fixed theta,
 *
 *   d(theta) = (sigmaA + sigmaS
 *               - log(eps*(1-exp(-theta*(rhoS-rhoA))))/theta) / rhoS,
 *
 * minimized over the feasible thetas. This is the deliverable of the domain: a
 * statistical delay guarantee, the quantity a service-level objective is written
 * against, as opposed to the mean delay returned by the queueing-theoretic
 * solvers.
 *
 * Port of matlab/src/api/snc/snc_perc_delay.m.
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
 * @param eps      violation probability, 0 < eps < 1
 * @param thetamax upper end of the theta search
 */
inline SncResult snc_perc_delay(const Envelope& arv, const Envelope& srv, double eps,
                                double thetamax = 1e3) {
    if (!(eps > 0.0 && eps < 1.0))
        throw UnsupportedError("snc_perc_delay: eps must lie in (0,1)");
    const SncResult r = snc_thetaopt(
        [&](double theta) {
            const detail::SncPair p = detail::snc_pair(arv, srv, theta);
            if (!p.ok) return std::numeric_limits<double>::infinity();
            return (p.a.sigma + p.s.sigma -
                    std::log(eps * (1.0 - std::exp(-theta * (p.s.rho - p.a.rho)))) / theta) /
                   p.s.rho;
        },
        thetamax);
    return SncResult{std::max(r.value, 0.0), r.theta};
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_PERC_DELAY_H
