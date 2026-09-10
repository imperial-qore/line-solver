/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_LAH_H
#define LINE_API_MOMENT_MOMENT_LAH_H

/**
 * Unsigned Lah numbers.
 *
 * Templated port of matlab/src/api/moment/moment_lah.m. Every operation is integer or
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

/** L(i,j) = L(i-1,j-1) + (i+j-1) L(i-1,j), L(0,0) = 1. */
template <class T>
Matrix<T> moment_lah(int n) {
    if (n < 0) throw InputError("moment_lah: the maximum order n must be nonnegative");
    Matrix<T> L(n + 1, n + 1, num_traits<T>::from_int(0));
    L(0, 0) = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= i; ++j)
            L(i, j) = L(i - 1, j - 1) + num_traits<T>::from_int(i + j - 1) * L(i - 1, j);
    return L;
}

}  // namespace moment
}  // namespace line

#endif
