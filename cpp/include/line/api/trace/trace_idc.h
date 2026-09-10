/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_IDC_H
#define LINE_API_TRACE_TRACE_IDC_H

/**
 * Index of dispersion for counts, estimated by its asymptotic equality with
 * the index of dispersion for intervals at a large aggregation level.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_idc.m, cross-checked
 * against `jar/src/main/java/jline/api/trace/Trace_var.java#trace_idc`.
 *
 * DIVERGENCE, MATLAB vs JAR: the aggregation level is
 * min(1000, ceil(n/30)) in MATLAB but min(1000, n/30) with integer division
 * in the JAR. They differ on every trace whose length is not a multiple of
 * 30, and for n < 30 the JAR asks for k = 0, which its trace_idi turns into a
 * division by zero (NaN). MATLAB is the reference.
 *
 * ARITHMETIC: as trace_idi, exact in Rational.
 */

#include <vector>

#include "line/api/trace/trace_idi.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** @param S the trace; the aggregation level is min(1000, ceil(n/30)). */
template <class T>
T trace_idc(const std::vector<T>& S) {
    detail::require_nonempty(S, "trace_idc");
    const long n = static_cast<long>(S.size());
    const long k = std::min<long>(1000, (n + 29) / 30);
    const TraceIdiResult<T> r = trace_idi(S, std::vector<long>(1, k));
    return r.idi[0];
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_IDC_H
