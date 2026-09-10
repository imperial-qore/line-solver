/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_BACKWARD_MOMENT_H
#define LINE_API_TRACE_MTRACE_BACKWARD_MOMENT_H

/**
 * Backward moments of a marked trace: the moments of the inter-arrival time
 * that PRECEDES an event of each class,
 *
 *   B(c,k) = (1/N) sum_{i: A_i = c} T_i^k,
 *
 * normalized by N/count_c when norm is set, so that M_k = sum_c B(c,k) p_c.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_backward_moment.m. In
 * MATLAB this is mtrace_moment(T,A,orders,0,NORM) with NORM defaulting to on,
 * since T_i is by construction the interval ending at event i.
 *
 * DIVERGENCE, MATLAB vs JAR: Mtrace_backward_moment.java computes a
 * DIFFERENT statistic. It accumulates absolute arrival epochs and, for each
 * class, forms the intervals between CONSECUTIVE EVENTS OF THAT SAME CLASS
 * (the class-c recurrence times, with the first one measured from the origin)
 * and returns their raw moments, indexed by the raw label 0..max(A) rather
 * than by the position in unique(A). That is the per-class inter-event time,
 * i.e. what mtrace_split feeds to a single-class estimator, not the backward
 * moment of the marked process, and it agrees with MATLAB only when every
 * event belongs to the same class. Its second entry point,
 * mtrace_backward_moment_conditional, has no MATLAB counterpart at all and is
 * not ported. MATLAB is the reference here.
 *
 * ARITHMETIC: sums of integer powers and a division, exact in Rational.
 */

#include <vector>

#include "line/api/trace/mtrace_moment.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/**
 * @param Tv     inter-arrival times
 * @param A      class labels
 * @param orders moment orders
 * @param norm   normalize by N/count_c (the MATLAB default)
 */
template <class T>
Matrix<T> mtrace_backward_moment(const std::vector<T>& Tv, const std::vector<int>& A,
                                 const std::vector<unsigned>& orders, bool norm = true) {
    return mtrace_moment(Tv, A, orders, false, norm);
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_BACKWARD_MOMENT_H
