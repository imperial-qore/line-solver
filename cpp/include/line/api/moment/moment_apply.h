/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_APPLY_H
#define LINE_API_MOMENT_MOMENT_APPLY_H

/**
 * Shared helper applying a triangular moment-transform table.
 *
 * Templated port of matlab/src/api/moment/moment_apply.m. Every operation is integer
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

/** Apply a lower-triangular transform table to a moment vector. */
template <class T>
std::vector<T> apply_table(const Matrix<T>& A, const std::vector<T>& v) {
    std::vector<T> r(v.size(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j <= i; ++j) r[i] += A(i, j) * v[j];
    return r;
}

}  // namespace moment
}  // namespace line

#endif
