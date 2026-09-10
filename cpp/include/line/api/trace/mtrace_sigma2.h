/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_SIGMA2_H
#define LINE_API_TRACE_MTRACE_SIGMA2_H

/**
 * Two-step class transition frequencies of a marked trace,
 *
 *   sigma(i,j,h) = #{t : A_t = i, A_{t+1} = j, A_{t+2} = h} / (N-2).
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_sigma2.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_sigma2.java (identical).
 *
 * The three-index array is returned flattened into a C x C^2 matrix, row i
 * and column j*C + h, so that it needs no tensor type; summing it gives 1 and
 * summing over h reproduces mtrace_sigma up to the different denominator.
 *
 * ARITHMETIC: counts over N-2, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/** @param L class labels; @return (C x C*C) matrix, entry (i, j*C+h). */
template <class T>
Matrix<T> mtrace_sigma2(const std::vector<int>& L) {
    if (L.size() < 3) throw InputError("mtrace_sigma2: at least three events are required");
    const std::vector<int> marks = detail::unique_labels(L);
    const std::size_t C = marks.size();
    Matrix<T> sigma(C, C * C, num_traits<T>::from_int(0));
    const T den = num_traits<T>::from_int(static_cast<long>(L.size() - 2));
    for (std::size_t i = 0; i < C; ++i) {
        for (std::size_t j = 0; j < C; ++j) {
            for (std::size_t h = 0; h < C; ++h) {
                long count = 0;
                for (std::size_t t = 0; t + 2 < L.size(); ++t)
                    if (L[t] == marks[i] && L[t + 1] == marks[j] && L[t + 2] == marks[h]) ++count;
                sigma(i, j * C + h) = num_traits<T>::from_int(count) / den;
            }
        }
    }
    return sigma;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_SIGMA2_H
