/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_STIRLING2_H
#define LINE_API_MOMENT_MOMENT_STIRLING2_H

/**
 * Stirling numbers of the second kind.
 *
 * Templated port of matlab/src/api/moment/moment_stirling2.m. Every operation is integer or
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

/** S(i,j) = j S(i-1,j) + S(i-1,j-1), S(0,0) = 1. */
template <class T>
Matrix<T> moment_stirling2(int n) {
    if (n < 0) throw InputError("moment_stirling2: the maximum order n must be nonnegative");
    Matrix<T> S(n + 1, n + 1, num_traits<T>::from_int(0));
    S(0, 0) = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= i; ++j)
            S(i, j) = num_traits<T>::from_int(j) * S(i - 1, j) + S(i - 1, j - 1);
    return S;
}

}  // namespace moment
}  // namespace line

#endif
