/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_UPFACTORIAL_FROM_FACTORIAL_H
#define LINE_API_MOMENT_MOMENT_UPFACTORIAL_FROM_FACTORIAL_H

/**
 * Rising-factorial moments from falling-factorial moments, via the Lah numbers.
 *
 * Templated port of matlab/src/api/moment/moment_upfactorial_from_factorial.m. Every operation is integer
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

/** fp_0 = 1; fp_i = sum_k L(i,k) f_k for i >= 1. */
template <class T>
std::vector<T> moment_upfactorial_from_factorial(const std::vector<T>& f) {
    const int n = static_cast<int>(f.size()) - 1;
    Matrix<T> L = moment_lah<T>(n);
    std::vector<T> fp(f.size(), num_traits<T>::from_int(0));
    fp[0] = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i)
        for (int k = 1; k <= i; ++k)
            fp[static_cast<std::size_t>(i)] += L(i, k) * f[static_cast<std::size_t>(k)];
    return fp;
}

}  // namespace moment
}  // namespace line

#endif
