/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMPP2_FIT_H
#define LINE_API_MAM_MMPP2_FIT_H

/**
 * MMPP(2) matching three moments and the lag-1 autocorrelation
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fit.m).
 *
 * The reference duplicates the closed form of mmpp2_fit3 verbatim, differing
 * only in that it converts the lag-1 autocorrelation into the decay rate
 * first, G2 = ACFLAG1 / ((1 - 1/SCV)/2); the denominator is the lag-0
 * autocorrelation of an MMPP(2), so ACFLAG1 is representable only in
 * [0, (1 - 1/SCV)/2]. The port calls mmpp2_fit3 rather than repeating the
 * expression, which keeps the two in step by construction.
 *
 * Gated on transcendental arithmetic through mmpp2_fit3.
 */

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmpp2_fit3.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** MMPP(2) with moments (E1, E2, E3) and lag-1 autocorrelation ACFLAG1. */
template <class T>
Map<T> mmpp2_fit(const T& E1, const T& E2, const T& E3, const T& ACFLAG1) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmpp2_fit requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (E1 == num_traits<T>::from_int(0)) throw InputError("mmpp2_fit: zero first moment");
    const T SCV = (E2 - E1 * E1) / (E1 * E1);
    if (SCV == num_traits<T>::from_int(0)) throw InputError("mmpp2_fit: zero SCV");
    const T rho0 = (one - one / SCV) / two;
    if (rho0 == num_traits<T>::from_int(0))
        throw InputError("mmpp2_fit: unit SCV admits no autocorrelation");
    return mmpp2_fit3(E1, E2, E3, T(ACFLAG1 / rho0));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMPP2_FIT_H
