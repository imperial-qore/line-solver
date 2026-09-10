/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_SPLIT_H
#define LINE_API_TRACE_MTRACE_SPLIT_H

/**
 * Splits a marked trace into its per-class traces: for each class, the
 * inter-arrival times BETWEEN CONSECUTIVE EVENTS OF THAT CLASS, with the
 * first interval measured from the origin.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_split.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_split.java (identical:
 * both prepend a zero epoch before differencing).
 *
 * The per-class traces partition the arrivals, and the sum of every class
 * trace equals the epoch of that class's last event, so the sum over classes
 * of the class sums is generally NOT the trace length.
 *
 * ARITHMETIC: cumulative sums and differences, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of mtrace_split. */
template <class T>
struct MtraceSplitResult {
    std::vector<int> labels;               ///< the distinct labels, increasing
    std::vector<std::vector<T>> traces;    ///< per-class inter-arrival times
};

/**
 * @param Tv inter-arrival times of the marked process
 * @param L  class labels
 */
template <class T>
MtraceSplitResult<T> mtrace_split(const std::vector<T>& Tv, const std::vector<int>& L) {
    detail::require_marked(Tv, L, "mtrace_split");
    MtraceSplitResult<T> out;
    out.labels = detail::unique_labels(L);
    const std::vector<T> cs = detail::cumsum0(Tv);  // cs[k] = epoch of event k
    out.traces.resize(out.labels.size());
    for (std::size_t c = 0; c < out.labels.size(); ++c) {
        T prev = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < L.size(); ++i) {
            if (L[i] != out.labels[c]) continue;
            const T epoch = cs[i + 1];
            out.traces[c].push_back(epoch - prev);
            prev = epoch;
        }
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_SPLIT_H
