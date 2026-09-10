/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMPP2_FIT2_H
#define LINE_API_MAM_MMPP2_FIT2_H

/**
 * MMPP(2) matching mean, SCV, skewness and the autocorrelation decay rate
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fit2.m).
 *
 * Converts (mean, scv, skew) to raw moments and delegates to mmpp2_fit3. A
 * unit SCV short-circuits to a Poisson process, where the MMPP(2) closed form
 * is degenerate. MATLAB then warns when the result fails map_isfeasible; this
 * port reports it through the result struct instead, since a warning printed
 * from a library is invisible to a caller that only has the return value.
 *
 * Gated on transcendental arithmetic: the third moment needs
 * (E2 - E1^2)^(3/2) and mmpp2_fit3 is itself gated.
 */

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmpp2_fit3.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** Result of the (mean, scv, skew)-parameterized MMPP(2) fits. */
template <class T>
struct Mmpp2FitResult {
    Map<T> map;
    bool feasible;  ///< map_isfeasible of the result at the given tolerance
};

/** MMPP(2) with the given mean, SCV, skewness and decay rate g2. */
template <class T>
Mmpp2FitResult<T> mmpp2_fit2(const T& mean, const T& scv, const T& skew, const T& g2) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmpp2_fit2 requires transcendental arithmetic");
    using fitdetail::num_sqrt;
    using fitdetail::pw;
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);

    Mmpp2FitResult<T> r;
    if (scv == one) {
        r.map = map_exponential_mean(mean);
        r.feasible = true;
        return r;
    }
    const T E1 = mean;
    const T E2 = (one + scv) * E1 * E1;
    const T v = E2 - E1 * E1;
    if (v < num_traits<T>::from_int(0)) throw InputError("mmpp2_fit2: negative variance");
    const T E3 = -(two * pw(E1, 3) - three * E1 * E2 - skew * v * num_sqrt(v));
    r.map = mmpp2_fit3(E1, E2, E3, g2);
    r.feasible = map_isfeasible(r.map, T(num_traits<T>::from_double(1e-10)));
    return r;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMPP2_FIT2_H
