/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_MOMENTS_HET_H
#define LINE_API_FJ_XMAX_MOMENTS_HET_H

/**
 * Moments of the maximum of heterogeneous exponentials by recurrence.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_moments_het.m.
 *
 *   M_m(lambda, n) = [ n M_m(lambda, n-1)
 *                      + sum_{j=1..m} lambda_j M_{m-1}(lambda \ j, n) ]
 *                    / sum_{j=1..m} lambda_j
 *
 * with M_m(lambda, 0) = 1 and M_0(., n) = 0 for n >= 1, which is the n-th
 * derivative of the transform recurrence of fj_lst_max_het at the origin.
 *
 * Eq. (30) of the survey prints the second sum WITHOUT the lambda_j weight;
 * that form is not the derivative of Eq. (29) and misses the textbook
 * two-variable answer, so the weight is restored here. fj_xmax_het is the
 * independent inclusion-exclusion check.
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
 * @param n      highest moment order, n >= 1
 * @return       moments of orders 1..n of the maximum
 */
template <class T>
std::vector<T> fj_xmax_moments_het(const std::vector<T>& lambda, unsigned n = 1) {
    const std::size_t K = lambda.size();
    if (K < 1) throw InputError("fj_xmax_moments_het: at least one rate is required");
    if (n < 1) throw InputError("fj_xmax_moments_het: the moment order must be a positive integer");
    if (K > 22)
        throw InputError("fj_xmax_moments_het: the recurrence enumerates 2^K sub-collections");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < K; ++i)
        if (!(lambda[i] > zero))
            throw InputError("fj_xmax_moments_het: all exponential rates must be positive");

    const std::size_t nmask = static_cast<std::size_t>(1) << K;
    // tab[mask][k] is the k-th moment over the sub-collection selected by mask
    std::vector<std::vector<T> > tab(nmask, std::vector<T>(n + 1, zero));
    for (std::size_t mask = 0; mask < nmask; ++mask) tab[mask][0] = one;

    for (unsigned order = 1; order <= n; ++order) {
        // The empty sub-collection has a zero maximum, so all its moments vanish
        tab[0][order] = zero;
        for (std::size_t mask = 1; mask < nmask; ++mask) {
            T tot = zero;
            T acc = num_traits<T>::from_int(static_cast<long>(order)) * tab[mask][order - 1];
            for (std::size_t j = 0; j < K; ++j)
                if (mask & (static_cast<std::size_t>(1) << j)) {
                    tot += lambda[j];
                    acc += lambda[j] * tab[mask ^ (static_cast<std::size_t>(1) << j)][order];
                }
            tab[mask][order] = acc / tot;
        }
    }

    std::vector<T> out(n);
    for (unsigned k = 1; k <= n; ++k) out[k - 1] = tab[nmask - 1][k];
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_MOMENTS_HET_H
