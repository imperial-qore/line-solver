/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_DTMC_MAKESTOCHASTIC_H
#define LINE_API_MC_DTMC_MAKESTOCHASTIC_H

/**
 * Normalize a non-negative matrix into a stochastic transition matrix.
 *
 * Templated port of matlab/lib/kpctoolbox/mc/dtmc_makestochastic.m and
 * jar/src/main/java/jline/api/mc/Dtmc_makestochastic.java. Each row with a
 * positive sum is divided by that sum, and the diagonal entry then absorbs
 * whatever deficit is left, clipped into [0,1]; a row that sums to zero is
 * replaced by the unit vector on its own state, which turns a dead state into
 * an absorbing one rather than leaving a substochastic row behind.
 *
 * After the division the row already sums to one, so the diagonal update is an
 * identity in exact arithmetic; it is kept because in floating point it is the
 * step that removes the accumulated rounding of the division, and dropping it
 * would make the double and exact paths disagree in the last bit.
 *
 * Every operation is a field operation, so this is exact at Rational.
 */

#include <cstddef>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * @param Pin matrix with non-negative entries
 * @return row-stochastic matrix of the same size
 */
template <class T>
Matrix<T> dtmc_makestochastic(const Matrix<T>& Pin) {
    const std::size_t n = Pin.rows();
    if (Pin.cols() != n) throw InputError("dtmc_makestochastic: matrix is not square");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    Matrix<T> P = Pin;
    for (std::size_t i = 0; i < n; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += P(i, j);
        if (s > zero) {
            for (std::size_t j = 0; j < n; ++j) P(i, j) /= s;
            T off = zero;
            for (std::size_t j = 0; j < n; ++j)
                if (j != i) off += P(i, j);
            T d = one - off;
            if (d < zero) d = zero;
            if (d > one) d = one;
            P(i, i) = d;
        } else {
            for (std::size_t j = 0; j < n; ++j) P(i, j) = zero;
            P(i, i) = one;
        }
    }
    return P;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_DTMC_MAKESTOCHASTIC_H
