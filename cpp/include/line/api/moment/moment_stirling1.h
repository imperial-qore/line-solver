/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_STIRLING1_H
#define LINE_API_MOMENT_MOMENT_STIRLING1_H

/**
 * Signed Stirling numbers of the first kind.
 *
 * Templated port of matlab/src/api/moment/moment_stirling1.m. Every operation is integer or
 * rational, so the exact instantiation returns the transform with no rounding:
 * moment conversions are exactly where double arithmetic hurts, since the
 * alternating binomial sums cancel catastrophically at high order.
 */

#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"
#include "line/api/moment/moment_stirlingcycle.h"

namespace line {
namespace moment {


/** s(i,j) = (-1)^(i-j) sigma(i,j). */
template <class T>
Matrix<T> moment_stirling1(int n) {
    Matrix<T> sigma = moment_stirlingcycle<T>(n);
    Matrix<T> s(n + 1, n + 1, num_traits<T>::from_int(0));
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= i; ++j)
            s(i, j) = ((i - j) % 2 == 0) ? sigma(i, j) : -sigma(i, j);
    return s;
}

}  // namespace moment
}  // namespace line

#endif
