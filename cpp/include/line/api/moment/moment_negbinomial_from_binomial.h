/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_NEGBINOMIAL_FROM_BINOMIAL_H
#define LINE_API_MOMENT_MOMENT_NEGBINOMIAL_FROM_BINOMIAL_H

/**
 * Negative-binomial moments from binomial moments.
 *
 * Templated port of matlab/src/api/moment/moment_negbinomial_from_binomial.m. Every operation is integer
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

/** bm_0 = 1; bm_i = sum_k C(i-1,k-1) b_k for i >= 1. */
template <class T>
std::vector<T> moment_negbinomial_from_binomial(const std::vector<T>& b) {
    const int n = static_cast<int>(b.size()) - 1;
    std::vector<T> bm(b.size(), num_traits<T>::from_int(0));
    bm[0] = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i)
        for (int k = 1; k <= i; ++k)
            bm[static_cast<std::size_t>(i)] += num_nck<T>(i - 1, k - 1) * b[static_cast<std::size_t>(k)];
    return bm;
}

}  // namespace moment
}  // namespace line

#endif
