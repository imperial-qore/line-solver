/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_STIRLINGCYCLE_H
#define LINE_API_MOMENT_MOMENT_STIRLINGCYCLE_H

/**
 * Unsigned Stirling numbers of the first kind (cycle numbers), orders 0..n.
 *
 * Templated port of matlab/src/api/moment/moment_stirlingcycle.m. Every operation is integer or
 * rational, so the exact instantiation returns the transform with no rounding:
 * moment conversions are exactly where double arithmetic hurts, since the
 * alternating binomial sums cancel catastrophically at high order.
 */

#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace moment {

/** sigma(i,j), (n+1) x (n+1) lower triangular, sigma(0,0) = 1. */
template <class T>
Matrix<T> moment_stirlingcycle(int n) {
    if (n < 0) throw InputError("moment_stirlingcycle: the maximum order n must be nonnegative");
    Matrix<T> s(n + 1, n + 1, num_traits<T>::from_int(0));
    s(0, 0) = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= i; ++j)
            s(i, j) = num_traits<T>::from_int(i - 1) * s(i - 1, j) + s(i - 1, j - 1);
    return s;
}

}  // namespace moment
}  // namespace line

#endif
