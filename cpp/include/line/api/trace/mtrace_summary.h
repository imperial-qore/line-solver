/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_SUMMARY_H
#define LINE_API_TRACE_MTRACE_SUMMARY_H

/**
 * Descriptor set of a marked trace: the first five raw moments of the
 * inter-arrival times, the autocorrelation function, the first two forward
 * and backward moments, the first two class-pair cross moments, the class
 * probabilities and the one-step class transition frequencies.
 *
 * This is exactly the input a marked MAP fitting procedure consumes.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_summary.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_summary.java. The two
 * assemble the same list; their entries differ wherever the underlying
 * function does (see mtrace_moment.h for the normalization defect that
 * affects F1, F2 and, in the JAR, B1 and B2 as well, and
 * mtrace_backward_moment.h for the JAR's different backward statistic).
 *
 * The acf lag set is 1..100 in both references; it is a parameter here
 * because trace_acf drops the lags a short trace cannot support.
 *
 * ARITHMETIC: all components are sample moments and frequencies, exact in
 * Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/mtrace_backward_moment.h"
#include "line/api/trace/mtrace_cross_moment.h"
#include "line/api/trace/mtrace_forward_moment.h"
#include "line/api/trace/mtrace_pc.h"
#include "line/api/trace/mtrace_sigma.h"
#include "line/api/trace/trace_acf.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/** Return value of mtrace_summary. */
template <class T>
struct MtraceSummary {
    std::vector<T> M;    ///< raw moments of order 1..5
    std::vector<T> acf;  ///< autocorrelation at the requested lags
    Matrix<T> F1, F2;    ///< forward moments of order 1 and 2
    Matrix<T> B1, B2;    ///< backward moments of order 1 and 2
    Matrix<T> C1, C2;    ///< class-pair cross moments of order 1 and 2
    std::vector<T> Pc;   ///< class probabilities
    Matrix<T> Pab;       ///< one-step class transition frequencies
};

/**
 * @param Tv       inter-arrival times
 * @param A        class labels
 * @param max_lag  largest acf lag (100 in both references)
 */
template <class T>
MtraceSummary<T> mtrace_summary(const std::vector<T>& Tv, const std::vector<int>& A,
                                int max_lag = 100) {
    detail::require_marked(Tv, A, "mtrace_summary");
    MtraceSummary<T> out;
    for (unsigned k = 1; k <= 5; ++k) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < Tv.size(); ++i) s += num_pow_int(Tv[i], k);
        out.M.push_back(s / num_traits<T>::from_int(static_cast<long>(Tv.size())));
    }
    std::vector<int> lags;
    for (int l = 1; l <= max_lag; ++l) lags.push_back(l);
    out.acf = trace_acf(Tv, lags);

    const std::vector<unsigned> o1(1, 1u), o2(1, 2u);
    out.F1 = mtrace_forward_moment(Tv, A, o1);
    out.F2 = mtrace_forward_moment(Tv, A, o2);
    out.B1 = mtrace_backward_moment(Tv, A, o1);
    out.B2 = mtrace_backward_moment(Tv, A, o2);
    out.C1 = mtrace_cross_moment(Tv, A, 1u).mc;
    out.C2 = mtrace_cross_moment(Tv, A, 2u).mc;
    out.Pc = mtrace_pc<T>(A);
    out.Pab = mtrace_sigma<T>(A);
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_SUMMARY_H
