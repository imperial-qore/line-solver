/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_BINOMIAL_FROM_NEGBINOMIAL_H
#define LINE_API_MOMENT_MOMENT_BINOMIAL_FROM_NEGBINOMIAL_H

/**
 * Binomial moments from negative-binomial moments.
 *
 * Templated port of matlab/src/api/moment/moment_binomial_from_negbinomial.m. Every operation is integer
 * or rational, so the exact instantiation returns the transform with no
 * rounding at all. That matters more here than almost anywhere else in the API:
 * the alternating binomial sums of the moment conversions cancel
 * catastrophically in double arithmetic once the order grows.
 */

#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace moment {

/** b_0 = 1; b_i = sum_k (-1)^(i-k) C(i-1,k-1) bm_k for i >= 1. */
template <class T>
std::vector<T> moment_binomial_from_negbinomial(const std::vector<T>& bm) {
    const int n = static_cast<int>(bm.size()) - 1;
    std::vector<T> b(bm.size(), num_traits<T>::from_int(0));
    b[0] = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i)
        for (int k = 1; k <= i; ++k) {
            const T term = num_nck<T>(i - 1, k - 1) * bm[static_cast<std::size_t>(k)];
            b[static_cast<std::size_t>(i)] += ((i - k) % 2 == 0) ? term : -term;
        }
    return b;
}

}  // namespace moment
}  // namespace line

#endif
