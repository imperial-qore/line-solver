/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_RAW_FROM_CENTRAL_H
#define LINE_API_MOMENT_MOMENT_RAW_FROM_CENTRAL_H

/**
 * Raw moments from central moments and the mean.
 *
 * Templated port of matlab/src/api/moment/moment_raw_from_central.m. Every operation is integer
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

/** m_i = sum_k C(i,k) mc_k m1^(i-k). */
template <class T>
std::vector<T> moment_raw_from_central(const std::vector<T>& mc, const T& m1) {
    const std::size_t n = mc.size();
    std::vector<T> m(n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k <= i; ++k)
            m[i] += num_nck<T>(static_cast<int>(i), static_cast<int>(k)) * mc[k] *
                    num_pow_int(m1, static_cast<unsigned>(i - k));
    return m;
}

}  // namespace moment
}  // namespace line

#endif
