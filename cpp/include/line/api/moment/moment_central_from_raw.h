/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_CENTRAL_FROM_RAW_H
#define LINE_API_MOMENT_MOMENT_CENTRAL_FROM_RAW_H

/**
 * Central moments from raw moments.
 *
 * Templated port of matlab/src/api/moment/moment_central_from_raw.m. Every operation is integer
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

/** mc_i = sum_k (-1)^(i-k) C(i,k) m_k m1^(i-k), with m1 = m[1]. */
template <class T>
std::vector<T> moment_central_from_raw(const std::vector<T>& m) {
    if (m.size() < 2)
        throw InputError(
            "moment_central_from_raw: the mean m1 is required, so m must have at least 2 entries");
    const std::size_t n = m.size();
    const T m1 = m[1];
    std::vector<T> mc(n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k <= i; ++k) {
            const T term = num_nck<T>(static_cast<int>(i), static_cast<int>(k)) * m[k] *
                           num_pow_int(m1, static_cast<unsigned>(i - k));
            mc[i] += ((i - k) % 2 == 0) ? term : -term;
        }
    return mc;
}

}  // namespace moment
}  // namespace line

#endif
