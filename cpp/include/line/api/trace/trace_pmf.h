/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_PMF_H
#define LINE_API_TRACE_TRACE_PMF_H

/**
 * Empirical probability mass function of a discrete trace (counts, batch
 * sizes, queue-length samples).
 *
 * Templated port of `jar/src/main/java/jline/api/trace/Trace_var.java#trace_pmf`,
 * cross-checked against matlab/lib/kpctoolbox/trace/trace_pmf.m.
 *
 * DIVERGENCE, MATLAB vs JAR: MATLAB computes `hist(X, max(X))' ./ numel(X)`,
 * i.e. it spreads max(X) equally spaced BINS over the range of the data and
 * returns their relative frequencies, while separately returning unique(X) as
 * the support. The two outputs then have different lengths and are not
 * aligned whenever the observed values are not exactly 1..max(X) -- the
 * MATLAB pmf is a histogram, not a pmf on the returned support. The JAR
 * counts the distinct observed values, which is the documented intent and is
 * what this port implements. The pmf sums to 1 by construction, so its
 * cumulative sum is a proper empirical CDF.
 *
 * ARITHMETIC: counts divided by the sample size, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of trace_pmf, mirroring [pmf, px]. */
template <class T>
struct TracePmfResult {
    std::vector<T> pmf;      ///< relative frequency of each distinct value
    std::vector<int> values;  ///< the distinct values, in increasing order
};

/** @param X the discrete trace. */
template <class T>
TracePmfResult<T> trace_pmf(const std::vector<int>& X) {
    if (X.empty()) throw InputError("trace_pmf: the trace is empty");
    TracePmfResult<T> out;
    out.values = detail::unique_labels(X);
    out.pmf.reserve(out.values.size());
    const T n = num_traits<T>::from_int(static_cast<long>(X.size()));
    for (std::size_t i = 0; i < out.values.size(); ++i) {
        long count = 0;
        for (std::size_t k = 0; k < X.size(); ++k)
            if (X[k] == out.values[i]) ++count;
        out.pmf.push_back(num_traits<T>::from_int(count) / n);
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_PMF_H
