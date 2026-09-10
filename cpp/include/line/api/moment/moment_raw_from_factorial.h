/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_RAW_FROM_FACTORIAL_H
#define LINE_API_MOMENT_MOMENT_RAW_FROM_FACTORIAL_H

/**
 * Raw moments from factorial moments, via the Stirling table of the second kind.
 *
 * Templated port of matlab/src/api/moment/moment_raw_from_factorial.m. Every operation is integer
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

/** m = S * f. */
template <class T>
std::vector<T> moment_raw_from_factorial(const std::vector<T>& f) {
    return apply_table(moment_stirling2<T>(static_cast<int>(f.size()) - 1), f);
}

}  // namespace moment
}  // namespace line

#endif
