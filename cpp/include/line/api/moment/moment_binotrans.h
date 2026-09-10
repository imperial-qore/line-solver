/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_BINOTRANS_H
#define LINE_API_MOMENT_MOMENT_BINOTRANS_H

/**
 * Binomial transform with alternating signs.
 *
 * Templated port of matlab/src/api/moment/moment_binotrans.m.
 *
 * Applied to a moment sequence m_i = E[X^i] it returns the moments of the unit
 * DOWNSHIFT, y_i = E[(X-1)^i]. It is NOT an involution -- an easy thing to
 * assume from the alternating signs, and wrong: its inverse is
 * moment_binotransinv, the unsigned transform.
 *
 * Every operation is integer or rational, so the exact instantiation returns
 * the transform with no rounding at all. That matters more here than almost anywhere else in the API:
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

/** y_i = sum_k (-1)^(i-k) C(i,k) x_k. */
template <class T>
std::vector<T> moment_binotrans(const std::vector<T>& x) {
    const std::size_t n = x.size();
    std::vector<T> y(n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k <= i; ++k) {
            const T term = num_nck<T>(static_cast<int>(i), static_cast<int>(k)) * x[k];
            y[i] += ((i - k) % 2 == 0) ? term : -term;
        }
    return y;
}

}  // namespace moment
}  // namespace line

#endif
