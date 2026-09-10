/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_FACTORIAL_FROM_BINOMIAL_H
#define LINE_API_MOMENT_MOMENT_FACTORIAL_FROM_BINOMIAL_H

/**
 * Falling-factorial moments from binomial moments.
 *
 * Templated port of matlab/src/api/moment/moment_factorial_from_binomial.m. Every operation is integer
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

/** f_i = i! b_i. */
template <class T>
std::vector<T> moment_factorial_from_binomial(const std::vector<T>& b) {
    std::vector<T> f(b.size());
    for (std::size_t i = 0; i < b.size(); ++i)
        f[i] = num_factorial<T>(static_cast<unsigned>(i)) * b[i];
    return f;
}

}  // namespace moment
}  // namespace line

#endif
