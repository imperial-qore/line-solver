/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LST_DET_H
#define LINE_API_AOI_LST_DET_H

/**
 * Laplace-Stieltjes transform of a deterministic (constant) distribution.
 *
 * Templated port of matlab/src/api/aoi/aoi_lst_det.m, cross-checked against
 * Aoi_lst.det in jar/src/main/java/jline/api/aoi/Aoi_lst.java (identical).
 *
 *   H*(s) = exp(-s d)
 *
 * static_assert(num_traits<T>::has_transcendental) -- this is the one
 * transform in the family that leaves the field. It is also the one that makes
 * an M/D/1 or D/M/1 AoI inexact no matter how the surrounding algebra is done.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"

namespace line {
namespace aoi {

/**
 * @param d constant value, > 0
 * @return  s -> exp(-s d)
 */
template <class T>
Lst<T> aoi_lst_det(const T& d) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lst_det requires transcendental arithmetic");
    detail::require_positive(d, "aoi_lst_det", "the constant d");
    return [d](const T& s) { return detail::num_exp(T(-s * d)); };
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LST_DET_H
