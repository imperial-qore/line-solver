/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_SIGMA_H
#define LINE_API_TRACE_MTRACE_SIGMA_H

/**
 * One-step class transition frequencies of a marked trace,
 *
 *   sigma(i,j) = #{t : A_t = i, A_{t+1} = j} / (N-1).
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_sigma.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_sigma.java (identical).
 *
 * Note that this is the JOINT frequency of the pair, not the conditional
 * transition probability: the whole matrix sums to 1, its rows do not.
 *
 * ARITHMETIC: counts over N-1, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/** @param L class labels; @return (C x C) matrix over the labels of unique(L). */
template <class T>
Matrix<T> mtrace_sigma(const std::vector<int>& L) {
    if (L.size() < 2) throw InputError("mtrace_sigma: at least two events are required");
    const std::vector<int> marks = detail::unique_labels(L);
    const std::size_t C = marks.size();
    Matrix<T> sigma(C, C, num_traits<T>::from_int(0));
    const T den = num_traits<T>::from_int(static_cast<long>(L.size() - 1));
    for (std::size_t i = 0; i < C; ++i) {
        for (std::size_t j = 0; j < C; ++j) {
            long count = 0;
            for (std::size_t t = 0; t + 1 < L.size(); ++t)
                if (L[t] == marks[i] && L[t + 1] == marks[j]) ++count;
            sigma(i, j) = num_traits<T>::from_int(count) / den;
        }
    }
    return sigma;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_SIGMA_H
