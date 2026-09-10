/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_FACTORIAL_FROM_RAW_H
#define LINE_API_MOMENT_MOMENT_FACTORIAL_FROM_RAW_H

/**
 * Factorial moments from raw moments, via the signed Stirling table.
 *
 * Templated port of matlab/src/api/moment/moment_factorial_from_raw.m. Every operation is integer
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
#include "line/api/moment/moment_stirling1.h"

namespace line {
namespace moment {

/** f = s * m, with s the signed Stirling numbers of the first kind. */
template <class T>
std::vector<T> moment_factorial_from_raw(const std::vector<T>& m) {
    return apply_table(moment_stirling1<T>(static_cast<int>(m.size()) - 1), m);
}

}  // namespace moment
}  // namespace line

#endif
