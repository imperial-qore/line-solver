/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_HET_H
#define LINE_API_FJ_XMAX_HET_H

/**
 * Exact moments of the maximum of heterogeneous exponentials.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_het.m.
 *
 *   E[Y^n] = sum over the nonempty subsets S of {1..K} of
 *              (-1)^(|S|+1) n! / ( sum_{i in S} lambda_i )^n
 *
 * Exact, at a cost of 2^K - 1 terms. At n = 1 and K = 2 it collapses to
 * 1/l1 + 1/l2 - 1/(l1+l2), and for equal rates to H_K/lambda.
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
 * @param n      moment order, n >= 1
 * @return       the n-th moment of the maximum
 */
template <class T>
T fj_xmax_het(const std::vector<T>& lambda, unsigned n = 1) {
    const std::size_t K = lambda.size();
    if (K < 1) throw InputError("fj_xmax_het: at least one rate is required");
    if (n < 1) throw InputError("fj_xmax_het: the moment order must be a positive integer");
    if (K > 24)
        throw InputError(
            "fj_xmax_het: inclusion-exclusion needs 2^K terms; use fj_xmax_moments_het instead");
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < K; ++i)
        if (!(lambda[i] > zero))
            throw InputError("fj_xmax_het: all exponential rates must be positive");

    T nfact = num_traits<T>::from_int(1);
    for (unsigned i = 2; i <= n; ++i) nfact *= num_traits<T>::from_int(static_cast<long>(i));

    T acc = zero;
    const std::size_t nmask = static_cast<std::size_t>(1) << K;
    for (std::size_t mask = 1; mask < nmask; ++mask) {
        T rate = zero;
        unsigned card = 0;
        for (std::size_t i = 0; i < K; ++i)
            if (mask & (static_cast<std::size_t>(1) << i)) {
                rate += lambda[i];
                ++card;
            }
        T den = num_traits<T>::from_int(1);
        for (unsigned e = 0; e < n; ++e) den *= rate;
        const T term = nfact / den;
        if (card % 2 == 1) acc += term; else acc -= term;
    }
    return acc;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_HET_H
