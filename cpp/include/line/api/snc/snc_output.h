/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_OUTPUT_H
#define LINE_API_SNC_OUTPUT_H

/**
 * Output (departure) arrival envelope of a flow leaving a server.
 *
 * For independent processes and a stable station,
 *
 *   rho   = rhoA,
 *   sigma = sigmaA + sigmaS - log(1-exp(-theta*(rhoS-rhoA)))/theta.
 *
 * The rate is conserved and the server adds burstiness. This is what carries a
 * flow across a feed-forward network one hop at a time; for a tandem traversed
 * by the same flow, `snc_conv` gives the tighter end-to-end result.
 *
 * Port of matlab/src/api/snc/snc_output.m.
 */

#include <cmath>
#include <limits>
#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

/**
 * @param arv   the arrival envelope entering the server
 * @param srv   the service envelope
 * @param theta Chernoff parameter, theta > 0
 */
inline Env snc_output(const Env& arv, const Env& srv, double theta) {
    if (theta <= 0) throw UnsupportedError("snc_output: theta must be positive");
    if (!std::isfinite(arv.rho) || !std::isfinite(srv.rho) || srv.rho <= arv.rho)
        return Env{std::numeric_limits<double>::infinity(), arv.rho};
    return Env{arv.sigma + srv.sigma -
                   std::log(1.0 - std::exp(-theta * (srv.rho - arv.rho))) / theta,
               arv.rho};
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_OUTPUT_H
