/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMPP2_FIT4_H
#define LINE_API_MAM_MMPP2_FIT4_H

/**
 * MMPP(2) matching mean, SCV, skewness and the lag-1 autocorrelation
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fit4.m).
 *
 * Same conversion as mmpp2_fit2 except that the autocorrelation is given at
 * lag 1 rather than as the decay rate, so it is divided by the lag-0
 * autocorrelation rho0 = (1 - 1/scv)/2 first. skew = -1 is the sentinel that
 * defers the choice of the third moment, which mmpp2_fit3 does not support;
 * unlike MATLAB, which would pass E3 = -1 straight into the closed form and
 * return a meaningless MAP, this port rejects it.
 *
 * Gated on transcendental arithmetic through mmpp2_fit3.
 */

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmpp2_fit2.h"
#include "line/api/mam/mmpp2_fit3.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** MMPP(2) with the given mean, SCV, skewness and lag-1 autocorrelation. */
template <class T>
Mmpp2FitResult<T> mmpp2_fit4(const T& mean, const T& scv, const T& skew, const T& acf1) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmpp2_fit4 requires transcendental arithmetic");
    using fitdetail::num_sqrt;
    using fitdetail::pw;
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);

    if (skew == -one)
        throw InputError("mmpp2_fit4: skew = -1 (automatic third moment) is not supported");
    if (scv == num_traits<T>::from_int(0)) throw InputError("mmpp2_fit4: zero SCV");
    const T rho0 = (one - one / scv) / two;
    if (rho0 == num_traits<T>::from_int(0))
        throw InputError("mmpp2_fit4: unit SCV admits no autocorrelation");

    const T E1 = mean;
    const T E2 = (one + scv) * E1 * E1;
    const T v = E2 - E1 * E1;
    if (v < num_traits<T>::from_int(0)) throw InputError("mmpp2_fit4: negative variance");
    const T E3 = -(two * pw(E1, 3) - three * E1 * E2 - skew * v * num_sqrt(v));

    Mmpp2FitResult<T> r;
    r.map = mmpp2_fit3(E1, E2, E3, T(acf1 / rho0));
    r.feasible = map_isfeasible(r.map, T(num_traits<T>::from_double(1e-10)));
    return r;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMPP2_FIT4_H
