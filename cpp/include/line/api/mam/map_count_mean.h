/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_COUNT_MEAN_H
#define LINE_API_MAM_MAP_COUNT_MEAN_H

/**
 * Mean of the counting process of a MAP at resolution t.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_count_mean.m. The number of
 * arrivals in a window of length t has stationary mean lambda t, with lambda
 * the stationary arrival rate pi D1 e, so the whole function is one linear
 * solve and a product.
 *
 * ARITHMETIC: exact. No exponential appears, and the mean of the counts is a
 * rational function of the entries of (D0, D1) -- which makes it the natural
 * consistency check on the counting-process functions that do need expm
 * (map_count_var, map_count_moment), since M1 from map_count_moment must
 * reproduce this value to the exponential's tolerance.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * @param m the MAP (D0, D1)
 * @param t window lengths
 * @return lambda t for each window length, in the order of t
 */
template <class T>
std::vector<T> map_count_mean(const Map<T>& m, const std::vector<T>& t) {
    const T lam = map_lambda(m);
    std::vector<T> out;
    out.reserve(t.size());
    for (std::size_t k = 0; k < t.size(); ++k) {
        if (t[k] < num_traits<T>::from_int(0))
            throw InputError("map_count_mean: negative window length");
        out.push_back(lam * t[k]);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_COUNT_MEAN_H
