/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_CROSS_MOMENT_H
#define LINE_API_TRACE_MTRACE_CROSS_MOMENT_H

/**
 * Class-pair cross moments of a marked trace: the k-th moment of the interval
 * that separates an event of class i from the next event, of class j,
 *
 *   MC(i,j) = mean{ T_t^k : A_{t-1} = i, A_t = j }.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_cross_moment.m,
 * cross-checked against
 * jar/src/main/java/jline/api/trace/Mtrace_cross_moment.java. The two agree
 * on every observed pair; they differ only in what they report for a pair
 * that never occurs, where MATLAB returns 0/0 = NaN and the JAR sets NaN
 * explicitly. Rational arithmetic has no NaN, so this port reports the
 * observed counts alongside the moments and leaves the unobserved entries at
 * zero: count(i,j) == 0 marks them.
 *
 * ARITHMETIC: sums of integer powers and a division, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/** Return value of mtrace_cross_moment. */
template <class T>
struct MtraceCrossMomentResult {
    Matrix<T> mc;         ///< (C x C) moments; entry undefined where count = 0
    Matrix<long> count;   ///< (C x C) number of observed transitions
};

/**
 * @param Tv inter-arrival times
 * @param L  class labels
 * @param k  moment order
 */
template <class T>
MtraceCrossMomentResult<T> mtrace_cross_moment(const std::vector<T>& Tv,
                                               const std::vector<int>& L, unsigned k) {
    detail::require_marked(Tv, L, "mtrace_cross_moment");
    const std::vector<int> marks = detail::unique_labels(L);
    const std::size_t C = marks.size();
    MtraceCrossMomentResult<T> out;
    out.mc = Matrix<T>(C, C, num_traits<T>::from_int(0));
    out.count = Matrix<long>(C, C, 0);

    std::vector<std::size_t> idx(L.size());
    for (std::size_t t = 0; t < L.size(); ++t) {
        std::size_t c = 0;
        while (c < C && marks[c] != L[t]) ++c;
        idx[t] = c;
    }
    for (std::size_t t = 1; t < Tv.size(); ++t) {
        const std::size_t i = idx[t - 1], j = idx[t];
        out.mc(i, j) += num_pow_int(Tv[t], k);
        out.count(i, j) += 1;
    }
    for (std::size_t i = 0; i < C; ++i)
        for (std::size_t j = 0; j < C; ++j)
            if (out.count(i, j) > 0) out.mc(i, j) /= num_traits<T>::from_int(out.count(i, j));
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_CROSS_MOMENT_H
