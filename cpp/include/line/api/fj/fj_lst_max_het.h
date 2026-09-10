/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_LST_MAX_HET_H
#define LINE_API_FJ_LST_MAX_HET_H

/**
 * Laplace-Stieltjes transform of the maximum of heterogeneous exponentials.
 *
 * Templated port of matlab/src/api/fj/fj_lst_max_het.m.
 *
 *   ( s + sum_{j=1..m} lambda_j ) L*_m(lambda, s)
 *       = sum_{j=1..m} lambda_j L*_{m-1}(lambda \ j, s)
 *
 * anchored at L*_0 = 1, because the maximum of an empty collection is zero.
 * The recurrence is swept bottom-up over the 2^K sub-collections, each keyed by
 * a bit mask, so every value is computed once.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param lambda the K positive exponential rates
 * @param s      transform argument, s >= 0
 * @return       L*(s) for the maximum of the K variables
 */
template <class T>
T fj_lst_max_het(const std::vector<T>& lambda, const T& s) {
    const std::size_t K = lambda.size();
    if (K < 1) throw InputError("fj_lst_max_het: at least one rate is required");
    if (K > 22) throw InputError("fj_lst_max_het: the recurrence enumerates 2^K sub-collections");
    const T zero = num_traits<T>::from_int(0);
    if (s < zero) throw InputError("fj_lst_max_het: the transform argument must be non-negative");
    for (std::size_t i = 0; i < K; ++i)
        if (!(lambda[i] > zero))
            throw InputError("fj_lst_max_het: all exponential rates must be positive");

    const std::size_t nmask = static_cast<std::size_t>(1) << K;
    std::vector<T> tab(nmask, zero);
    tab[0] = num_traits<T>::from_int(1);
    for (std::size_t mask = 1; mask < nmask; ++mask) {
        T tot = zero, acc = zero;
        for (std::size_t j = 0; j < K; ++j)
            if (mask & (static_cast<std::size_t>(1) << j)) {
                tot += lambda[j];
                acc += lambda[j] * tab[mask ^ (static_cast<std::size_t>(1) << j)];
            }
        tab[mask] = acc / (s + tot);
    }
    return tab[nmask - 1];
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_LST_MAX_HET_H
