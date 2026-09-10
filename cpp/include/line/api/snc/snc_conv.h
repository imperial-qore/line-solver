/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_CONV_H
#define LINE_API_SNC_CONV_H

/**
 * Min-plus convolution of two service envelopes (tandem concatenation).
 *
 * Two stations traversed in series offer the flow their min-plus convolution.
 * For independent servers with exponential-form envelopes, summing the geometric
 * series over the intermediate epoch gives
 *
 *   rho   = min(rho1,rho2),
 *   sigma = sigma1 + sigma2 - log(1-exp(-theta*|rho1-rho2|))/theta.
 *
 * This is the pay-bursts-only-once result: the end-to-end burst term grows
 * additively rather than the per-station delay bounds being summed. The series
 * diverges at equal rates, so those are handled by shifting the slower server
 * down by `delta`, the usual regularization; delta trades rate against burst and
 * is an argument so it can be optimized jointly with theta.
 *
 * Port of matlab/src/api/snc/snc_conv.m. Original: F. Ciucu, A. Burchard,
 * J. Liebeherr, IEEE Trans. Inf. Theory 52(6), 2300-2312, 2006.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

/**
 * @param s1    envelope of the first station
 * @param s2    envelope of the second station
 * @param theta Chernoff parameter, theta > 0
 * @param delta rate separation used when the two rates coincide; a nonpositive
 *              value selects the default 1e-2*min(rho1,rho2)
 */
inline Env snc_conv(const Env& s1, const Env& s2, double theta, double delta = -1.0) {
    if (theta <= 0) throw UnsupportedError("snc_conv: theta must be positive");
    if (delta <= 0) delta = 1e-2 * std::min(s1.rho, s2.rho);
    if (delta <= 0) throw UnsupportedError("snc_conv: delta must be positive");
    if (!std::isfinite(s1.rho) || !std::isfinite(s2.rho))
        return Env{std::numeric_limits<double>::infinity(), std::min(s1.rho, s2.rho)};
    double gap = std::fabs(s1.rho - s2.rho);
    double rho;
    if (gap <= delta) {
        gap = delta;  // equal rates: shift the slower server down to close the series
        rho = std::min(s1.rho, s2.rho) - delta;
    } else {
        rho = std::min(s1.rho, s2.rho);
    }
    if (rho <= 0) return Env{std::numeric_limits<double>::infinity(), rho};
    return Env{s1.sigma + s2.sigma - std::log(1.0 - std::exp(-theta * gap)) / theta, rho};
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_CONV_H
