/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LST_EXP_H
#define LINE_API_AOI_LST_EXP_H

/**
 * Laplace-Stieltjes transform of an exponential distribution.
 *
 * Templated port of matlab/src/api/aoi/aoi_lst_exp.m, cross-checked against
 * Aoi_lst.exp in jar/src/main/java/jline/api/aoi/Aoi_lst.java (identical).
 *
 *   H*(s) = mu / (mu + s)
 *
 * A rational function of s, so it is exact in the field: an exact-arithmetic
 * caller gets the true transform value at any rational s, which is what makes
 * the AoI closed forms checkable term by term.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"

namespace line {
namespace aoi {

/**
 * @param mu rate, > 0 (mean 1/mu)
 * @return   s -> mu/(mu+s)
 */
template <class T>
Lst<T> aoi_lst_exp(const T& mu) {
    detail::require_positive(mu, "aoi_lst_exp", "the rate mu");
    return [mu](const T& s) { return T(mu / (mu + s)); };
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LST_EXP_H
