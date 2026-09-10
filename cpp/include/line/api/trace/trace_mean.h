/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_MEAN_H
#define LINE_API_TRACE_TRACE_MEAN_H

/**
 * Sample mean of a trace.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_mean.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Trace_mean.java (identical: both
 * are sum/n).
 *
 * ARITHMETIC: a sum and one division, exact in Rational.
 */

#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** (1/n) sum_i S(i). */
template <class T>
T trace_mean(const std::vector<T>& S) {
    detail::require_nonempty(S, "trace_mean");
    T s = num_traits<T>::from_int(0);
    for (const T& v : S) s += v;
    return s / num_traits<T>::from_int(static_cast<long>(S.size()));
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_MEAN_H
