/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DQSYS_DQSYS_GEOXGEO1_H
#define LINE_API_DQSYS_DQSYS_GEOXGEO1_H

/**
 * Geo^X/Geo/1: the discrete-time single-server queue with batch arrivals.
 *
 * Templated port of matlab/src/api/qsys/dqsys_geoxgeo1.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_geoxgeo1.java.
 *
 * A batch arrives in a slot with probability a and carries X >= 1 jobs; the
 * server completes a job in a slot with probability s. With A(z) the pgf of
 * the number of jobs arriving in one slot the slot-boundary content obeys
 *
 *   P(z) = p_0 s (z-1) A(z) / ( z - A(z)(s + (1-s) z) ),   p_0 = 1 - lambda/s
 *
 * and differentiating at z = 1 gives
 *
 *   E[N] = lambda + ( a E[X(X-1)]/2 + lambda (1-s) ) / (s - lambda)
 *
 * so the batch law enters only through its first two factorial moments. For a
 * geometric batch with parameter beta, E[X] = 1/beta and
 * E[X(X-1)] = 2(1-beta)/beta^2, and at beta = 1 the result collapses onto
 * dqsys_geogeo1.
 *
 * Everything is a rational function of (a, beta, s), so the exact
 * instantiation carries no rounding. That matters in the same regime as for
 * Geo/Geo/1: the interesting case is s - lambda small, and that is exactly
 * where the double evaluation of the E[N] quotient loses digits.
 *
 * The MATLAB entry point returns the pgf as a function handle taking (z, A(z));
 * here that is the free function dqsys_geoxgeo1_pgf, which takes the result
 * struct and the same pair.
 */

#include <string>

#include "line/api/dqsys/dqsys_geogeo1.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace dqsys {

template <class T>
struct GeoXGeo1Result {
    GeoConvention convention = GeoConvention::LAS_DA;
    T batchArrivalProb;               ///< a
    T batchMean;                      ///< E[X]
    T batchSecondFactorialMoment;     ///< E[X(X-1)]
    T serviceProb;                    ///< s
    T arrivalRate;                    ///< lambda = a E[X]
    T throughput;                     ///< lambda
    T utilization;                    ///< lambda/s
    T boundaryEmptyProb;              ///< 1 - lambda/s, at the slot boundary
    T meanQueueLength;
    T meanWaitingQueue;
    T meanSojournTime;
    T meanWaitingTime;
    T meanServiceTime;
};

/**
 * Geo^X/Geo/1 for an arbitrary batch law given by its first two factorial
 * moments. This is MATLAB's local dqsys_geoxgeo1_moments, exposed here because
 * a C++ caller has no other way to reach it.
 */
template <class T>
GeoXGeo1Result<T> dqsys_geoxgeo1_moments(const T& a, const T& batchMean,
                                        const T& batchSecondFactorial, const T& s,
                                        GeoConvention convention = GeoConvention::LAS_DA) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (a <= zero || a > one) throw InputError("dqsys_geoxgeo1: a must lie in (0,1]");
    if (s <= zero || s > one) throw InputError("dqsys_geoxgeo1: s must lie in (0,1]");
    if (batchMean < one)
        throw InputError(
            "dqsys_geoxgeo1: mean batch size must be at least 1: a batch that arrives carries at "
            "least one job");
    if (batchSecondFactorial < zero)
        throw InputError("dqsys_geoxgeo1: E[X(X-1)] must be non-negative");
    // factorial-moment inequality rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
    const T minSecondFactorial = batchMean * batchMean - batchMean;
    if (num_traits<T>::is_exact) {
        if (batchSecondFactorial < minSecondFactorial)
            throw InputError(
                "dqsys_geoxgeo1: E[X(X-1)] is below E[X]^2-E[X], so the batch moments describe no "
                "random variable");
    } else {
        const T slack = T(num_traits<T>::from_double(1e-9)) *
                        (minSecondFactorial > one ? minSecondFactorial : one);
        if (batchSecondFactorial < minSecondFactorial - slack)
            throw InputError(
                "dqsys_geoxgeo1: E[X(X-1)] is below E[X]^2-E[X], so the batch moments describe no "
                "random variable");
    }

    const T lambda = a * batchMean;
    if (lambda >= s) throw InputError("dqsys_geoxgeo1: load lambda/s must be strictly less than 1");

    GeoXGeo1Result<T> r;
    r.convention = convention;
    r.batchArrivalProb = a;
    r.batchMean = batchMean;
    r.batchSecondFactorialMoment = batchSecondFactorial;
    r.serviceProb = s;
    r.arrivalRate = lambda;
    r.throughput = lambda;
    r.utilization = lambda / s;
    r.boundaryEmptyProb = one - r.utilization;

    // P'(1) from the generating function.
    const T meanAtBoundary =
        lambda + (a * batchSecondFactorial / two + lambda * (one - s)) / (s - lambda);
    const T meanSojournAtBoundary = meanAtBoundary / lambda;
    r.meanWaitingTime = meanSojournAtBoundary - one / s;
    r.meanWaitingQueue = lambda * r.meanWaitingTime;

    if (convention == GeoConvention::LAS_DA) {
        r.meanQueueLength = meanAtBoundary;
        r.meanSojournTime = meanSojournAtBoundary;
        r.meanServiceTime = one / s;
    } else {
        // One departure earlier: the epoch drops exactly the departures of the
        // slot, whose rate is lambda, hence one slot of sojourn.
        r.meanQueueLength = meanAtBoundary - lambda;
        r.meanSojournTime = meanSojournAtBoundary - one;
        r.meanServiceTime = (one - s) / s;
    }
    return r;
}

/**
 * @param a    per-slot probability that a batch arrives, in (0,1]
 * @param beta geometric batch-size parameter, in (0,1]; E[X] = 1/beta
 * @param s    per-slot service completion probability, in (0,1]
 * @param convention slot-boundary convention (late arrival, early arrival)
 */
template <class T>
GeoXGeo1Result<T> dqsys_geoxgeo1(const T& a, const T& beta, const T& s,
                                GeoConvention convention = GeoConvention::LAS_DA) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (beta <= zero || beta > one)
        throw InputError("dqsys_geoxgeo1: beta must lie in (0,1]");
    const T batchMean = one / beta;
    const T batchSecondFactorial = two * (one - beta) / (beta * beta);
    return dqsys_geoxgeo1_moments(a, batchMean, batchSecondFactorial, s, convention);
}

/**
 * Probability generating function of the stationary queue length.
 *
 * @param r  a result of dqsys_geoxgeo1
 * @param z  argument in (0,1]
 * @param Az the value A(z) of the slot-arrival pgf at the same z
 *
 * No pmf is offered: for a general batch law the stationary distribution has
 * no elementary closed form, so only the generating function is exact.
 */
template <class T>
T dqsys_geoxgeo1_pgf(const GeoXGeo1Result<T>& r, const T& z, const T& Az) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (z <= zero || z > one) throw InputError("dqsys_geoxgeo1: pgf argument z must lie in (0,1]");
    if (z == one) return one;
    const T s = r.serviceProb;
    const T denom = z - Az * (s + (one - s) * z);
    if (denom == zero) throw NumericError("dqsys_geoxgeo1: pgf denominator vanishes");
    const T boundary = r.boundaryEmptyProb * s * (z - one) * Az / denom;
    if (r.convention == GeoConvention::LAS_DA) return boundary;
    return boundary * (s / z + one - s) + r.boundaryEmptyProb * s * (one - one / z);
}

}  // namespace dqsys
}  // namespace line

#endif  // LINE_API_DQSYS_DQSYS_GEOXGEO1_H
