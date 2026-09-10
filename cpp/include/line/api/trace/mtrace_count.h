/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_COUNT_H
#define LINE_API_TRACE_MTRACE_COUNT_H

/**
 * Per-class count process of a marked trace on a fixed resolution: the number
 * of events of each class in the successive windows of length `t` that start
 * at the first arrival epoch.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_count.m.
 *
 * REFERENCE DEFECT (MATLAB): the class test is `A == c` with c the LOOP INDEX
 * 1..length(unique(A)), not the label unique(A)(c). Every trace whose labels
 * are not exactly 1..C is therefore counted against the wrong classes, and
 * labels outside that range are never counted at all -- with 0/1 labels, for
 * instance, all the zeros are dropped and the class-1 column is reported as
 * class 1 while column 2 stays empty. This port compares against the label,
 * which is what the header comment of the function describes.
 *
 * DIVERGENCE, MATLAB vs JAR: Mtrace_count.java is not this function. It
 * computes windowed count STATISTICS (mean, variance, index of dispersion and
 * skewness of the count process, plus a multiscale sweep) over
 * floor(total/window) windows, and its generateCountProcess increments every
 * window from the previous index to the current one, so its counts are
 * partial cumulative sums rather than per-window counts. It has no MATLAB
 * counterpart and is not ported; the count-statistics wrapper it exists for
 * is a composition of this function with trace_var / trace_skew.
 *
 * ARITHMETIC: comparisons of cumulative sums, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/** Return value of mtrace_count. */
template <class T>
struct MtraceCountResult {
    std::vector<int> labels;  ///< the distinct labels, increasing
    Matrix<long> counts;      ///< (periods x C) events per window and class
};

/**
 * @param Tv inter-arrival times
 * @param A  class labels
 * @param t  window length (the resolution)
 */
template <class T>
MtraceCountResult<T> mtrace_count(const std::vector<T>& Tv, const std::vector<int>& A,
                                  const T& t) {
    detail::require_marked(Tv, A, "mtrace_count");
    if (t <= num_traits<T>::from_int(0))
        throw InputError("mtrace_count: the resolution must be positive");
    const std::size_t n = Tv.size();

    T span = num_traits<T>::from_int(0);
    for (std::size_t i = 1; i < n; ++i) span += Tv[i];
    long periods = 0;
    while (num_traits<T>::from_int(periods) * t < span) ++periods;  // ceil(span/t)

    MtraceCountResult<T> out;
    out.labels = detail::unique_labels(A);
    const std::size_t C = out.labels.size();
    out.counts = Matrix<long>(static_cast<std::size_t>(periods), C, 0);
    const std::vector<T> cs = detail::cumsum0(Tv);

    for (long i = 1; i <= periods; ++i) {
        const T tstart = Tv[0] + num_traits<T>::from_int(i - 1) * t;
        const T tend = Tv[0] + num_traits<T>::from_int(i) * t;
        for (std::size_t k = 0; k < n; ++k) {
            const T& epoch = cs[k + 1];
            if (!(epoch > tstart) || !(epoch <= tend)) continue;
            for (std::size_t c = 0; c < C; ++c)
                if (A[k] == out.labels[c]) out.counts(static_cast<std::size_t>(i - 1), c) += 1;
        }
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_COUNT_H
