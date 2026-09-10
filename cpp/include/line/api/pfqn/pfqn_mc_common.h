/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MC_COMMON_H
#define LINE_API_PFQN_MC_COMMON_H

/**
 * Randomness scaffolding shared by the Monte Carlo normalizing-constant
 * estimators (pfqn_mci, pfqn_is, pfqn_ld_is, pfqn_oi_is, pfqn_pas_is, pfqn_ls,
 * pfqn_mmsample2) and by the perfect sampler pfqn_cftp.
 *
 * This header is NOT a port of a MATLAB function. It exists because the MATLAB
 * references draw from the global MATLAB stream (rand / randi / mvnrnd, seeded
 * out of band by `rng(options.seed)`), and a library must not carry a hidden
 * global stream. Every estimator in this tree therefore takes a
 * `std::mt19937_64&` as an explicit argument and consumes it through the
 * helpers below.
 *
 * REPRODUCIBILITY CONTRACT, which every estimator's own header repeats:
 *
 *  - The estimate is comparable to MATLAB only IN DISTRIBUTION, never stream
 *    for stream. MATLAB's Mersenne Twister, its uniform-to-integer mapping and
 *    its normal transform all differ from the ones here, so the same seed does
 *    NOT produce the same sample path and the two estimates agree only up to
 *    Monte Carlo error. What IS comparable is the estimand: both target the
 *    exact constant of pfqn_ca, and both converge to it at the 1/sqrt(n) rate.
 *
 *  - Within this port, reproducibility requires passing a generator in the
 *    same state. Two calls with generators seeded identically, on identical
 *    inputs, in the same build, produce bit-identical output; the generator is
 *    advanced by the call, so a second call on the same generator object does
 *    not repeat the first. No estimator seeds, re-seeds or copies the
 *    generator internally.
 *
 *  - The uniform-to-integer map here is rejection based, not modulo based, so
 *    it is unbiased and identical on every platform and standard-library
 *    version. std::uniform_int_distribution and std::normal_distribution are
 *    deliberately avoided: their output is implementation defined, which would
 *    make a fixed-seed regression test non-portable.
 *
 * Arithmetic: everything here is inherently inexact (uniform deviates are
 * dyadic approximations of a continuous law, the normal transform needs
 * log/cos), so each consumer gates on num_traits<T>::has_transcendental.
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** The generator type every Monte Carlo entry point in this tree accepts. */
using McRng = std::mt19937_64;

/**
 * Uniform deviate on [0,1) with 53 significant bits, as a double. The top 53
 * bits of one 64-bit draw are used, so exactly one generator step is consumed
 * per deviate and the mapping is fully specified.
 */
inline double mc_uniform01(McRng& g) {
    return static_cast<double>(g() >> 11) * (1.0 / 9007199254740992.0);
}

/** The same deviate materialized in the working arithmetic. */
template <class T>
T mc_uniform(McRng& g) {
    return num_traits<T>::from_double(mc_uniform01(g));
}

/**
 * Uniform integer on [0, n), unbiased by rejection. Consumes one generator
 * step per attempt; the rejection probability is below 2^-64 * n, so for the
 * class counts these estimators use it never rejects in practice.
 */
inline std::uint64_t mc_uniform_int(McRng& g, std::uint64_t n) {
    if (n == 0) throw InputError("mc_uniform_int: empty range");
    if (n == 1) return 0;
    const std::uint64_t threshold = (0u - n) % n;  // 2^64 mod n
    std::uint64_t r;
    do {
        r = g();
    } while (r < threshold);
    return r % n;
}

/**
 * Standard normal deviate by the Box-Muller transform. Two uniforms are drawn
 * and only the cosine branch is kept, so the routine holds no state between
 * calls: a generator handed to two different estimators cannot be
 * cross-contaminated by a cached second variate.
 */
inline double mc_normal01(McRng& g) {
    double u1 = mc_uniform01(g);
    const double u2 = mc_uniform01(g);
    // log(0) would be -inf; the smallest representable positive deviate keeps
    // the transform finite without perturbing the distribution measurably.
    if (u1 <= 0.0) u1 = 1.0 / 9007199254740992.0;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586476925286766559 * u2);
}

/**
 * log(mean(exp(v))), computed by factoring out the maximum so that the
 * exponentials stay in range. MATLAB's logmeanexp, which every estimator that
 * averages log-weights calls.
 */
inline double mc_logmeanexp(const std::vector<double>& v) {
    if (v.empty()) return -std::numeric_limits<double>::infinity();
    double m = -std::numeric_limits<double>::infinity();
    for (double x : v)
        if (x > m) m = x;
    if (!std::isfinite(m)) return m;
    double acc = 0.0;
    for (double x : v) acc += std::exp(x - m);
    return m + std::log(acc / static_cast<double>(v.size()));
}

/**
 * exp of a log-domain value, materialized in the working arithmetic. Overflows
 * to infinity in double exactly where the references do, and stays in range
 * for the high-precision backends.
 */
template <class T>
T mc_exp(double lv) {
    using std::exp;
    return exp(num_traits<T>::from_double(lv));
}

/**
 * log(n!) for a non-negative integer n, the factln / gammaln(1+n) of the
 * references, accumulated in the working arithmetic so the high-precision
 * backends do not lose the digits a double lgamma would drop.
 */
template <class T>
double mc_log_factorial(long n) {
    if (n < 0) throw InputError("mc_log_factorial: negative argument");
    if (n < 2) return 0.0;
    return num_traits<T>::log_as_double(num_factorial<T>(static_cast<unsigned>(n)));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MC_COMMON_H
