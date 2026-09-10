/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_TYPES_H
#define LINE_API_SNC_TYPES_H

/**
 * Shared types of the stochastic network calculus domain.
 *
 * An envelope is the pair `(sigma(theta), rho(theta))` in the exponential form
 *
 *   E[exp(theta*A(s,t))] <= exp(theta*(rho*(t-s) + sigma)),   theta > 0,
 *
 * and its service counterpart with the sign of theta reversed. The whole domain
 * passes envelopes as FUNCTIONS of theta rather than as numbers, because every
 * bound is an infimum over theta and a composition (superposition, leftover
 * service, concatenation, departure) has to be re-evaluated at whatever theta
 * the search asks for.
 *
 * ARITHMETIC. This domain is DOUBLE-ONLY, deliberately. Every entry point is an
 * exp/log expression minimized numerically over theta, so there is no exact or
 * multiprecision instantiation to offer: a Rational cannot carry `exp(theta)`
 * and a Real would buy digits the Chernoff search does not have. The templated
 * SolverBA arm converts at the boundary with `num_traits<T>::to_double` and
 * `from_double`, exactly as `sylvester.h` does for the same reason.
 *
 * Port of matlab/src/api/snc, cross-checked against jline.api.snc.
 *
 * Reference: M. Fidler, A. Rizk, "A Guide to the Stochastic Network Calculus",
 * IEEE Communications Surveys and Tutorials 17(1), 92-105, 2015.
 */

#include <functional>

namespace line {
namespace snc {

/** The pair (sigma, rho) of an envelope evaluated at one theta. */
struct Env {
    double sigma = 0.0;
    double rho = 0.0;
};

/** An envelope as a function of the Chernoff parameter. */
using Envelope = std::function<Env(double)>;

/** A bound together with the theta that attains it. */
struct SncResult {
    /** The bound: a violation probability, a quantile or a mean bound. */
    double value = 0.0;
    /** The minimizing theta, NaN when no feasible theta exists. */
    double theta = 0.0;
};

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_TYPES_H
