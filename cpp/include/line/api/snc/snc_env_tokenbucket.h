/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_ENV_TOKENBUCKET_H
#define LINE_API_SNC_ENV_TOKENBUCKET_H

/**
 * Deterministic token-bucket arrival envelope.
 *
 * A flow policed by a (b,r) token bucket satisfies `A(s,t) <= b + r*(t-s)` with
 * probability one, so the envelope is constant in theta. This is the
 * deterministic network calculus arrival curve read as a degenerate MGF
 * envelope, so it mixes freely with the stochastic ones.
 *
 * Port of matlab/src/api/snc/snc_env_tokenbucket.m.
 */

#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

/**
 * @param b bucket depth, units of work
 * @param r token rate, work per slot
 */
inline Env snc_env_tokenbucket(double b, double r) {
    if (b < 0 || r < 0)
        throw UnsupportedError("snc_env_tokenbucket: b and r must be nonnegative");
    return Env{b, r};
}

/** @param theta accepted and ignored, for signature compatibility */
inline Env snc_env_tokenbucket(double b, double r, double theta) {
    if (theta <= 0) throw UnsupportedError("snc_env_tokenbucket: theta must be positive");
    return snc_env_tokenbucket(b, r);
}

/** The same envelope as a function of theta. */
inline Envelope snc_env_tokenbucket_fn(double b, double r) {
    return [b, r](double theta) { return snc_env_tokenbucket(b, r, theta); };
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_ENV_TOKENBUCKET_H
