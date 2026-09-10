/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_LEFTOVER_H
#define LINE_API_SNC_LEFTOVER_H

/**
 * Leftover service envelope under blind (arbitrary) multiplexing.
 *
 * A server with envelope (sigmaS,rhoS) shared with a cross flow of arrival
 * envelope (sigmaX,rhoX) leaves the flow of interest `S-X`, whose envelope is
 * `rho = rhoS-rhoX`, `sigma = sigmaS+sigmaX`.
 *
 * The subtraction of the two exponential forms is exact when the processes are
 * INDEPENDENT; otherwise the pair must be split by Hoelder's inequality, which
 * this elementary version does not do. A nonpositive rho means the cross traffic
 * can exhaust the server, which the bound functions report as a violation
 * probability of 1.
 *
 * Port of matlab/src/api/snc/snc_leftover.m.
 */

#include "line/api/snc/snc_types.h"

namespace line {
namespace snc {

/**
 * @param srv the service envelope
 * @param cross the cross-flow arrival envelope
 */
inline Env snc_leftover(const Env& srv, const Env& cross) {
    return Env{srv.sigma + cross.sigma, srv.rho - cross.rho};
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_LEFTOVER_H
