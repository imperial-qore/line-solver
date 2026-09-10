/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DA_DA_TRAFFIC_SUPERPOS_H
#define LINE_API_DA_DA_TRAFFIC_SUPERPOS_H

/**
 * Superposition of independent renewal flows (Whitt's QNA stationary-interval
 * method).
 *
 * Templated port of matlab/src/api/da/da_traffic_superpos.m. Given flows with
 * rates lambda(i) and squared coefficients of variation a2(i), the merged flow
 * is assigned the rate-weighted mixture
 *
 *   d2 = sum_i lambda(i) a2(i) / sum_i lambda(i),
 *
 * flows with a non-finite rate being dropped first (MATLAB's isfinite mask).
 * The rate of the merged flow is sum_i lambda(i) and is not returned: the
 * decomposition step that calls this already holds it.
 *
 * A weighted mean is a sum and one division, so this is a finite field
 * computation and instantiates at exact arithmetic with no rounding: for
 * rational rates and SCVs the merged SCV is the exact rational mixture. No
 * transcendental gate.
 *
 * The finiteness mask only ever removes anything in an inexact instantiation;
 * an exact rational is finite by construction, so for T = Rational the mask is
 * the identity and every flow is kept.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace da {

namespace detail {

/** MATLAB isfinite: always true in an exact instantiation. */
template <class T>
bool da_is_finite(const T& v) {
    if (num_traits<T>::is_exact) return true;
    return std::isfinite(num_traits<T>::to_double(v));
}

}  // namespace detail

/**
 * @param lambda (m) flow rates; entries that are not finite are ignored
 * @param a2     (m) squared coefficients of variation of the same flows
 * @return the squared coefficient of variation of the superposed flow
 */
template <class T>
T da_traffic_superpos(const std::vector<T>& lambda, const std::vector<T>& a2) {
    if (lambda.size() != a2.size())
        throw InputError("da_traffic_superpos: lambda and a2 have different lengths");
    if (lambda.empty()) throw InputError("da_traffic_superpos: no flows to superpose");

    T num = num_traits<T>::from_int(0);
    T den = num_traits<T>::from_int(0);
    std::size_t kept = 0;
    for (std::size_t i = 0; i < lambda.size(); ++i) {
        if (!detail::da_is_finite(lambda[i])) continue;
        num += a2[i] * lambda[i];
        den += lambda[i];
        ++kept;
    }
    if (kept == 0) throw InputError("da_traffic_superpos: every flow rate is non-finite");
    if (den == num_traits<T>::from_int(0))
        throw NumericError("da_traffic_superpos: the superposed flow has zero rate");
    return num / den;
}

}  // namespace da
}  // namespace line

#endif  // LINE_API_DA_DA_TRAFFIC_SUPERPOS_H
