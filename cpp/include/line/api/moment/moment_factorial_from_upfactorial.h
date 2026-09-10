/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_FACTORIAL_FROM_UPFACTORIAL_H
#define LINE_API_MOMENT_MOMENT_FACTORIAL_FROM_UPFACTORIAL_H

/**
 * Falling-factorial moments from rising-factorial moments, via the Lah numbers.
 *
 * Templated port of matlab/src/api/moment/moment_factorial_from_upfactorial.m. Every operation is integer
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
#include "line/api/moment/moment_lah.h"

namespace line {
namespace moment {

/** f_0 = 1; f_i = sum_k (-1)^(i-k) L(i,k) fp_k for i >= 1. */
template <class T>
std::vector<T> moment_factorial_from_upfactorial(const std::vector<T>& fp) {
    const int n = static_cast<int>(fp.size()) - 1;
    Matrix<T> L = moment_lah<T>(n);
    std::vector<T> f(fp.size(), num_traits<T>::from_int(0));
    f[0] = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i)
        for (int k = 1; k <= i; ++k) {
            const T term = L(i, k) * fp[static_cast<std::size_t>(k)];
            f[static_cast<std::size_t>(i)] += ((i - k) % 2 == 0) ? term : -term;
        }
    return f;
}

}  // namespace moment
}  // namespace line

#endif
