/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_BINOMIAL_FROM_FACTORIAL_H
#define LINE_API_MOMENT_MOMENT_BINOMIAL_FROM_FACTORIAL_H

/**
 * Binomial moments from falling-factorial moments.
 *
 * Templated port of matlab/src/api/moment/moment_binomial_from_factorial.m. Every operation is integer
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

/** b_i = f_i / i!. */
template <class T>
std::vector<T> moment_binomial_from_factorial(const std::vector<T>& f) {
    std::vector<T> b(f.size());
    for (std::size_t i = 0; i < f.size(); ++i)
        b[i] = f[i] / num_factorial<T>(static_cast<unsigned>(i));
    return b;
}

}  // namespace moment
}  // namespace line

#endif
