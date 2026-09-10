/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_SRV_RATE_H
#define LINE_API_SNC_SRV_RATE_H

/**
 * MGF service envelope of a constant-rate work-conserving server.
 *
 * `S(s,t) = C*(t-s)`, so `(sigma, rho) = (0, C)` for every theta. This is the
 * elementary service element; a station shared by cross traffic follows from it
 * through `snc_leftover` and a tandem through `snc_conv`.
 *
 * USE `snc_srv_exp` INSTEAD WHENEVER THE WORK UNIT IS THE JOB: pairing this
 * element with a job-counting arrival envelope models an M/D/1 and understates
 * the delay of an exponential server.
 *
 * Port of matlab/src/api/snc/snc_srv_rate.m.
 */

#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

/** @param C server capacity, work per slot */
inline Env snc_srv_rate(double C) {
    if (C <= 0) throw UnsupportedError("snc_srv_rate: C must be positive");
    return Env{0.0, C};
}

/** @param theta accepted and ignored, for signature compatibility */
inline Env snc_srv_rate(double C, double theta) {
    if (theta <= 0) throw UnsupportedError("snc_srv_rate: theta must be positive");
    return snc_srv_rate(C);
}

/** The same envelope as a function of theta. */
inline Envelope snc_srv_rate_fn(double C) {
    return [C](double theta) { return snc_srv_rate(C, theta); };
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_SRV_RATE_H
