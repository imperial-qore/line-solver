/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_MEAN_H
#define LINE_API_TRACE_MTRACE_MEAN_H

/**
 * Per-type sample mean of a trace.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/mtrace_mean.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_mean.java (identical,
 * including the 0-indexed type vector: type c ranges over 0..ntypes-1, unlike
 * the m3a mtrace_* family which uses the values of unique(A)).
 *
 * Both references return NaN for a type that never occurs. Rational
 * arithmetic has no NaN, so the per-type counts are returned alongside the
 * means and count[c] == 0 marks an undefined entry.
 *
 * ARITHMETIC: a sum and one division per type, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of mtrace_mean. */
template <class T>
struct MtraceMeanResult {
    std::vector<T> mean;      ///< one mean per type; undefined where count = 0
    std::vector<long> count;  ///< samples of each type
};

/**
 * @param Tv     the trace
 * @param ntypes number of types
 * @param type   type of each sample, in 0..ntypes-1
 */
template <class T>
MtraceMeanResult<T> mtrace_mean(const std::vector<T>& Tv, long ntypes,
                                const std::vector<int>& type) {
    detail::require_marked(Tv, type, "mtrace_mean");
    if (ntypes <= 0) throw InputError("mtrace_mean: ntypes must be positive");
    MtraceMeanResult<T> out;
    out.mean.assign(static_cast<std::size_t>(ntypes), num_traits<T>::from_int(0));
    out.count.assign(static_cast<std::size_t>(ntypes), 0);
    for (long c = 0; c < ntypes; ++c) {
        T sum = num_traits<T>::from_int(0);
        long ctr = 0;
        for (std::size_t i = 0; i < Tv.size(); ++i) {
            if (type[i] != static_cast<int>(c)) continue;
            sum += Tv[i];
            ++ctr;
        }
        out.count[static_cast<std::size_t>(c)] = ctr;
        if (ctr > 0)
            out.mean[static_cast<std::size_t>(c)] = sum / num_traits<T>::from_int(ctr);
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_MEAN_H
