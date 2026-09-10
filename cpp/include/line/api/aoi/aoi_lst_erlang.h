/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LST_ERLANG_H
#define LINE_API_AOI_LST_ERLANG_H

/**
 * Laplace-Stieltjes transform of an Erlang-k distribution.
 *
 * Templated port of matlab/src/api/aoi/aoi_lst_erlang.m, cross-checked against
 * Aoi_lst.erlang in jar/src/main/java/jline/api/aoi/Aoi_lst.java (identical).
 *
 *   H*(s) = (mu / (mu + s))^k,  mean k/mu
 *
 * The exponent is an integer, so this stays a rational function of s and is
 * exact in the field.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @param k  number of phases, >= 1
 * @param mu per-phase rate, > 0
 * @return   s -> (mu/(mu+s))^k
 */
template <class T>
Lst<T> aoi_lst_erlang(unsigned k, const T& mu) {
    if (k < 1) throw InputError("aoi_lst_erlang: the shape k must be a positive integer");
    detail::require_positive(mu, "aoi_lst_erlang", "the rate mu");
    return [k, mu](const T& s) { return num_pow_int(T(mu / (mu + s)), k); };
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LST_ERLANG_H
