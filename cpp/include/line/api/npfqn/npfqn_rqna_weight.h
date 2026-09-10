/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_RQNA_WEIGHT_H
#define LINE_API_NPFQN_RQNA_WEIGHT_H

/**
 * Canonical reflected-Brownian-motion correlation weight w*(t) used by the
 * Robust Queueing Network Analyzer (RQNA).
 *
 * Templated port of matlab/src/api/npfqn/npfqn_rqna_weight.m, cross-checked
 * against jar/src/main/java/jline/api/npfqn/Npfqn_rqna_weight.java (identical
 * term for term, including both numerical guards).
 *
 *   w*(t)  = 1 - (1 - c*(t)) / (2 t)
 *   c*(t)  = 2 (1 - 2t - t^2) Phi^c(sqrt(t)) + 2 sqrt(t) phi(sqrt(t)) (1 + t)
 *
 * with Phi^c the standard-normal complementary cdf and phi its density. The
 * weight increases monotonically from w*(0) = 0 to w*(Inf) = 1.
 * Reference: W. Whitt and W. You (2018), "A Robust Queueing Network Analyzer
 * Based on Indices of Dispersion", eqs. (24)-(25).
 *
 * Arithmetic. Phi^c is an erfc and phi is an exp, neither of which exists in
 * the field of the inputs, so this requires transcendental arithmetic and
 * cannot be instantiated at T = Rational.
 */

#include <vector>

#include "line/api/npfqn/npfqn_types.h"
#include "line/num/number.h"

namespace line {
namespace npfqn {

/**
 * @param t nonnegative time argument; t <= 0 gives 0 and t = Inf gives 1
 * @return the weight w*(t), clamped to [0,1] as in MATLAB and the JAR
 */
template <class T>
T npfqn_rqna_weight(const T& t) {
    static_assert(num_traits<T>::has_transcendental,
                  "npfqn_rqna_weight requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    if (t <= zero) return zero;
    if (!detail::num_isfinite(t)) return one;

    const T st = detail::num_sqrt(t);
    // Phi^c(sqrt(t)) = 0.5 erfc(sqrt(t)/sqrt(2))
    const T phic = num_traits<T>::from_rational(1, 2) * detail::num_erfc(T(st / detail::num_sqrt(two)));
    // phi(sqrt(t)) = exp(-t/2)/sqrt(2 pi)
    const T twopi = two * boost::math::constants::pi<T>();
    const T phi = detail::num_exp(T(-t / two)) / detail::num_sqrt(twopi);
    const T cstar = two * (one - two * t - t * t) * phic + two * st * phi * (one + t);

    T w;
    if (t < num_traits<T>::from_double(1e-6)) {
        // limit w*(t) -> 0 as t -> 0; taking the difference below would be
        // catastrophic cancellation there
        w = zero;
    } else {
        w = one - (one - cstar) / (two * t);
    }
    // numerical guard: w* is a weight in [0,1]
    if (w < zero) w = zero;
    if (w > one) w = one;
    return w;
}

/** Elementwise form, mirroring the MATLAB array argument. */
template <class T>
std::vector<T> npfqn_rqna_weight(const std::vector<T>& t) {
    std::vector<T> w(t.size());
    for (std::size_t i = 0; i < t.size(); ++i) w[i] = npfqn_rqna_weight(t[i]);
    return w;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_RQNA_WEIGHT_H
