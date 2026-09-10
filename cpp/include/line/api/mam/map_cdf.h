/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_CDF_H
#define LINE_API_MAM_MAP_CDF_H

/**
 * Cumulative distribution of the inter-arrival time of a MAP.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_cdf.m, cross-checked against
 * jar/src/main/java/jline/api/mam/Map_cdf.java.
 *
 * Conditional on the phase pie seen by an arrival, the next inter-arrival time
 * is phase-type with representation (pie, D0), so
 *
 *   F(t) = 1 - pie exp(D0 t) e,     f(t) = pie exp(D0 t) (-D0) e.
 *
 * ARITHMETIC: exp(D0 t) is a tolerance-controlled approximation, so both
 * functions require transcendental arithmetic; pie itself is exact (a linear
 * solve, see map_moment.h).
 *
 * DIVERGENCE, MATLAB vs JAR: MATLAB evaluates exp(D0 t) directly for every
 * process. The JAR switches to Foxglynn uniformization when D0 has no negative
 * off-diagonal entry and to a direct exponential otherwise, because
 * uniformization is invalid for an ME/RAP representation whose D0 is not a
 * sub-generator. This port follows MATLAB and always uses the exponential,
 * which is correct for both cases. NOTE, MEASURED: the uniformization path is
 * not merely faster, it is also the ACCURATE one at small t. F(t) = 1 - s here
 * loses every significant digit once s is within a few ulp of one, so for a
 * value far below the rounding of s this function returns round-off of either
 * sign; uniformization with an explicit absorbing state sums only nonnegative
 * terms and stays exact. pfqn_stdf therefore evaluates its level CDFs by that
 * construction (detail::stdf_hypoexp_cdf) rather than through this function.
 * Any caller of map_cdf that needs a small CDF accurately has the same problem
 * and no fix. The JAR additionally clamps a NaN CDF value to the previous
 * finite one, which is a workaround, not a definition, and is not reproduced.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Cumulative distribution of the inter-arrival time at the given points.
 *
 * @param m      the MAP (D0, D1)
 * @param points evaluation times, each >= 0
 * @return F(t) = Pr[T <= t] in the order of points
 */
template <class T>
std::vector<T> map_cdf(const Map<T>& m, const std::vector<T>& points) {
    static_assert(num_traits<T>::has_transcendental, "map_cdf requires transcendental arithmetic");
    const std::vector<T> pie = map_pie(m);
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    std::vector<T> out;
    out.reserve(points.size());
    for (std::size_t k = 0; k < points.size(); ++k) {
        if (points[k] < zero) throw InputError("map_cdf: negative evaluation point");
        if (points[k] == zero) {
            out.push_back(zero);
            continue;
        }
        const std::vector<T> v = vecmul(pie, expm(m.D0, points[k]));
        T s = zero;
        for (const T& x : v) s += x;
        out.push_back(one - s);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_CDF_H
