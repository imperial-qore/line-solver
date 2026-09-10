/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_BINOTRANSINV_H
#define LINE_API_MOMENT_MOMENT_BINOTRANSINV_H

/**
 * Inverse binomial transform.
 *
 * Templated port of matlab/src/api/moment/moment_binotransinv.m. Every operation is integer
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

/** x_i = sum_k C(i,k) y_k. */
template <class T>
std::vector<T> moment_binotransinv(const std::vector<T>& y) {
    const std::size_t n = y.size();
    std::vector<T> x(n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k <= i; ++k)
            x[i] += num_nck<T>(static_cast<int>(i), static_cast<int>(k)) * y[k];
    return x;
}

}  // namespace moment
}  // namespace line

#endif
