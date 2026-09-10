/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_SCV_H
#define LINE_API_TRACE_TRACE_SCV_H

/**
 * Squared coefficient of variation of a trace, var/mean^2.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_scv.m, cross-checked
 * against `jar/src/main/java/jline/api/trace/Trace_var.java#trace_scv`. The
 * expression is the same; the two differ only through the variance
 * denominator, see trace_var.h.
 *
 * ARITHMETIC: a ratio of sample moments, exact in Rational.
 */

#include <vector>

#include "line/api/trace/trace_mean.h"
#include "line/api/trace/trace_types.h"
#include "line/api/trace/trace_var.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/**
 * @param S        the trace
 * @param unbiased variance denominator, see trace_var
 */
template <class T>
T trace_scv(const std::vector<T>& S, bool unbiased = true) {
    const T mu = trace_mean(S);
    if (mu == num_traits<T>::from_int(0)) throw NumericError("trace_scv: the trace has zero mean");
    return trace_var(S, unbiased) / (mu * mu);
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_SCV_H
