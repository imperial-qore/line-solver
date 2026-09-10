/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_IDI_H
#define LINE_API_TRACE_TRACE_IDI_H

/**
 * Index of dispersion for intervals,
 *
 *   IDI(k) = k * var(S_t + ... + S_{t+k-1}) / mean(S_t + ... + S_{t+k-1})^2,
 *
 * the standard burstiness descriptor of Sriram and Whitt (JSAC 6, 1986).
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_idi.m, cross-checked
 * against `jar/src/main/java/jline/api/trace/Trace_var.java#trace_idi`.
 *
 * REFERENCE DEFECT (both codebases, identical): the vector of aggregated
 * samples is allocated with length(S)-k entries but the loop fills only
 * length(S)-k-1 of them, so a spurious ZERO sample is always included in the
 * variance and the mean. The bias is O(1/(n-k)) and vanishes for long traces,
 * but it is real and it is reproduced here deliberately: removing it would
 * silently change every IDI/IDC value LINE reports. The `drop_trailing_zero`
 * flag computes the intended statistic instead.
 *
 * DIVERGENCE, MATLAB vs JAR: MATLAB's var uses the n-1 denominator, the JAR's
 * trace_var uses n; MATLAB is the reference. The 'aggregate-mix' branch of
 * trace_idi.m is not ported -- it partitions the trace by a per-sample
 * aggregation-count vector that no caller in LINE supplies, and its MATLAB
 * implementation overwrites Sk on every outer iteration, so what it returns
 * is the statistic of the last partition only.
 *
 * ARITHMETIC: sums, a variance and a division, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_mean.h"
#include "line/api/trace/trace_types.h"
#include "line/api/trace/trace_var.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of trace_idi, mirroring [IDIk, support]. */
template <class T>
struct TraceIdiResult {
    std::vector<T> idi;       ///< one value per requested k
    std::vector<long> support;  ///< the number of points each value rests on
};

/**
 * @param k_set               aggregation levels
 * @param aggregate_n         0 for raw samples; n > 0 when S is already the
 *                            sum of n inter-arrivals, which uses k/n windows
 *                            (MATLAB's 'aggregate' option)
 * @param drop_trailing_zero  true removes the spurious zero sample described
 *                            above; false (default) reproduces the references
 * @param S the interarrival-time trace
 */
template <class T>
TraceIdiResult<T> trace_idi(const std::vector<T>& S, const std::vector<long>& k_set,
                            long aggregate_n = 0, bool drop_trailing_zero = false) {
    detail::require_nonempty(S, "trace_idi");
    const long n = static_cast<long>(S.size());
    TraceIdiResult<T> out;
    for (std::size_t a = 0; a < k_set.size(); ++a) {
        const long k = k_set[a];
        if (k <= 0) throw InputError("trace_idi: the aggregation level must be positive");
        long keff = k;
        long support = n - k - 1;
        if (aggregate_n > 0) {
            keff = k / aggregate_n;
            if (keff <= 0) throw InputError("trace_idi: the aggregation level is below n");
            support = n / keff;
        }
        const long len = n - keff;  // allocated length in both references
        const long filled = n - keff - 1;
        if (filled < 1) throw InputError("trace_idi: the aggregation level exceeds the trace length");
        const long used = drop_trailing_zero ? filled : len;
        if (used < 2) throw InputError("trace_idi: too few aggregated samples for a variance");

        std::vector<T> Sk(static_cast<std::size_t>(used), num_traits<T>::from_int(0));
        for (long t = 0; t < filled; ++t) {
            T s = num_traits<T>::from_int(0);
            for (long j = t; j < t + keff; ++j) s += S[static_cast<std::size_t>(j)];
            Sk[static_cast<std::size_t>(t)] = s;
        }
        const T mu = trace_mean(Sk);
        if (mu == num_traits<T>::from_int(0))
            throw NumericError("trace_idi: the aggregated samples have zero mean");
        out.idi.push_back(num_traits<T>::from_int(k) * trace_var(Sk) / (mu * mu));
        out.support.push_back(support);
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_IDI_H
