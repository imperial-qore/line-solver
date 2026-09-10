/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_FORWARD_MOMENT_H
#define LINE_API_TRACE_MTRACE_FORWARD_MOMENT_H

/**
 * Forward moments of a marked trace: the moments of the inter-arrival time
 * that FOLLOWS an event of each class,
 *
 *   F(c,k) = (1/(N-1)) sum_{i<N: A_i = c} T_{i+1}^k,
 *
 * normalized by N/count_c when norm is set, so that M_k = sum_c F(c,k) p_c.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_forward_moment.m,
 * cross-checked against
 * jar/src/main/java/jline/api/trace/Mtrace_forward_moment.java.
 *
 * In MATLAB this function is literally mtrace_moment(T,A,orders,1,NORM) with
 * NORM defaulting to on, and it is implemented here the same way. The two JAR
 * divergences described in mtrace_moment.h (the sum divided by count_c
 * already in the unnormalized branch, and the (N-1)/count_c normalization
 * factor) apply verbatim to Mtrace_forward_moment.java as well.
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
Matrix<T> mtrace_forward_moment(const std::vector<T>& Tv, const std::vector<int>& A,
                                const std::vector<unsigned>& orders, bool norm = true) {
    return mtrace_moment(Tv, A, orders, true, norm);
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_FORWARD_MOMENT_H
