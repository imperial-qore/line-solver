/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_IAT2COUNTS_H
#define LINE_API_TRACE_MTRACE_IAT2COUNTS_H

/**
 * Per-class counting processes of a marked trace: for each arrival, how many
 * events of each class fall in the window of length `scale` that starts at
 * that arrival.
 *
 * Templated port of the pure-MATLAB branch of
 * matlab/lib/m3a/m3a/mtrace/mtrace_iat2counts.m (the function prefers a MEX
 * implementation, mtrace_iat2counts_native, when it is compiled; the MATLAB
 * fallback is the specification and is what is ported), cross-checked against
 * jar/src/main/java/jline/api/trace/Mtrace_iat2counts.java.
 *
 * DIVERGENCE, MATLAB vs JAR: the JAR compares CT[cur+1] against CT[i+1] where
 * MATLAB compares CT(cur+1) against CT(i), so its window ends one arrival
 * early; it also seeds the search with max(i, previousCur) instead of
 * MATLAB's (i-1)+C(i-1). The seeds are only a speed-up (the scan is monotone,
 * so any lower start converges to the same window), but the endpoint is not,
 * and the JAR counts are systematically one arrival short. MATLAB is the
 * reference.
 *
 * Note also that MATLAB's own speed-up seed `cur = (i-1) + C(i-1)` uses
 * LINEAR indexing into the count matrix, so it reads the class-1 count rather
 * than the total; it is an underestimate of the previous window end and
 * therefore harmless, but it makes the loop slower, not faster, on
 * multi-class traces.
 *
 * The series is truncated at the first window that reaches the end of the
 * trace, because from there on the counts are censored.
 *
 * ARITHMETIC: additions and comparisons, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/** Return value of mtrace_iat2counts. */
template <class T>
struct MtraceCountsResult {
    std::vector<int> labels;  ///< the distinct labels, increasing
    Matrix<long> counts;      ///< (rows x C), column c for labels[c]
};

/**
 * @param Tv    inter-arrival times
 * @param A     class labels
 * @param scale window length
 */
template <class T>
MtraceCountsResult<T> mtrace_iat2counts(const std::vector<T>& Tv, const std::vector<int>& A,
                                        const T& scale) {
    detail::require_marked(Tv, A, "mtrace_iat2counts");
    if (scale <= num_traits<T>::from_int(0))
        throw InputError("mtrace_iat2counts: the time scale must be positive");
    const long n = static_cast<long>(Tv.size());
    if (n < 2) throw InputError("mtrace_iat2counts: at least two events are required");
    MtraceCountsResult<T> out;
    out.labels = detail::unique_labels(A);
    const std::size_t K = out.labels.size();
    const std::vector<T> cs = detail::cumsum0(Tv);  // cs[k] = MATLAB CT(k)

    std::vector<std::vector<long>> rows;
    for (long i = 1; i <= n - 1; ++i) {
        long cur = i;
        bool censored = false;
        while (cs[static_cast<std::size_t>(cur + 1)] - cs[static_cast<std::size_t>(i)] <= scale) {
            ++cur;
            if (cur == n) {
                censored = true;
                break;
            }
        }
        std::vector<long> row(K, 0);
        for (long t = i + 1; t <= cur; ++t)
            for (std::size_t j = 0; j < K; ++j)
                if (A[static_cast<std::size_t>(t - 1)] == out.labels[j]) ++row[j];
        rows.push_back(row);
        if (censored) break;
    }
    out.counts = Matrix<long>(rows.size(), K, 0);
    for (std::size_t r = 0; r < rows.size(); ++r)
        for (std::size_t j = 0; j < K; ++j) out.counts(r, j) = rows[r][j];
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_IAT2COUNTS_H
