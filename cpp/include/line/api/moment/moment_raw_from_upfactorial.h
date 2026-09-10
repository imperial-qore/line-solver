/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_RAW_FROM_UPFACTORIAL_H
#define LINE_API_MOMENT_MOMENT_RAW_FROM_UPFACTORIAL_H

/**
 * Raw moments from rising-factorial moments.
 *
 * Templated port of matlab/src/api/moment/moment_raw_from_upfactorial.m. Every operation is integer
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
#include "line/api/moment/moment_apply.h"
#include "line/api/moment/moment_stirling2.h"

namespace line {
namespace moment {

/** m_i = sum_j (-1)^(i-j) S(i,j) fp_j. */
template <class T>
std::vector<T> moment_raw_from_upfactorial(const std::vector<T>& fp) {
    const int n = static_cast<int>(fp.size()) - 1;
    Matrix<T> S = moment_stirling2<T>(n);
    Matrix<T> Tm(n + 1, n + 1, num_traits<T>::from_int(0));
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= i; ++j) Tm(i, j) = ((i - j) % 2 == 0) ? S(i, j) : -S(i, j);
    return apply_table(Tm, fp);
}

}  // namespace moment
}  // namespace line

#endif
