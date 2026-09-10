/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_COV_H
#define LINE_API_TRACE_MTRACE_COV_H

/**
 * Class-pair covariance matrices of a marked trace.
 *
 * For every ordered pair (c1,c2) the two masked series
 *
 *   X0(i) = T_i     if A_i = c1     else 0,
 *   X1(i) = T_{i+1} if A_{i+1} = c2 else 0,   i = 1..N-1,
 *
 * are formed and their 2x2 sample covariance matrix (denominator N-2, i.e.
 * MATLAB's cov on N-1 observations) is returned.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_cov.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_cov.java (identical,
 * including the unbiased denominator and the indexing of classes by the raw
 * label 1..max(A) rather than by position in unique(A)).
 *
 * ARITHMETIC: sample second moments, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/**
 * @param Tv inter-arrival times
 * @param A  class labels, positive
 * @return cov[c1-1][c2-1], each a 2x2 matrix
 */
template <class T>
std::vector<std::vector<Matrix<T>>> mtrace_cov(const std::vector<T>& Tv,
                                               const std::vector<int>& A) {
    detail::require_marked(Tv, A, "mtrace_cov");
    const std::size_t N = A.size();
    if (N < 3) throw InputError("mtrace_cov: at least three events are required");
    int C = 0;
    for (std::size_t k = 0; k < N; ++k)
        if (A[k] > C) C = A[k];
    if (C <= 0) throw InputError("mtrace_cov: class labels must be positive");

    const T zero = num_traits<T>::from_int(0);
    const T den = num_traits<T>::from_int(static_cast<long>(N - 2));
    std::vector<std::vector<Matrix<T>>> COV(
        static_cast<std::size_t>(C),
        std::vector<Matrix<T>>(static_cast<std::size_t>(C), Matrix<T>(2, 2, zero)));

    std::vector<T> x(N - 1), y(N - 1);
    for (int c1 = 1; c1 <= C; ++c1) {
        for (int c2 = 1; c2 <= C; ++c2) {
            for (std::size_t i = 0; i + 1 < N; ++i) {
                x[i] = (A[i] == c1) ? Tv[i] : zero;
                y[i] = (A[i + 1] == c2) ? Tv[i + 1] : zero;
            }
            T mx = zero, my = zero;
            for (std::size_t i = 0; i + 1 < N; ++i) {
                mx += x[i];
                my += y[i];
            }
            const T n1 = num_traits<T>::from_int(static_cast<long>(N - 1));
            mx /= n1;
            my /= n1;
            T sxx = zero, sxy = zero, syy = zero;
            for (std::size_t i = 0; i + 1 < N; ++i) {
                const T dx = x[i] - mx, dy = y[i] - my;
                sxx += dx * dx;
                sxy += dx * dy;
                syy += dy * dy;
            }
            Matrix<T> c(2, 2, zero);
            c(0, 0) = sxx / den;
            c(0, 1) = sxy / den;
            c(1, 0) = c(0, 1);
            c(1, 1) = syy / den;
            COV[static_cast<std::size_t>(c1 - 1)][static_cast<std::size_t>(c2 - 1)] = c;
        }
    }
    return COV;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_COV_H
