/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_MOMENT_H
#define LINE_API_TRACE_MTRACE_MOMENT_H

/**
 * Empirical class-dependent moments of a marked trace.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_moment.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_moment.java.
 *
 *   after = false (Horvath variables):  M(c,j) = (1/N)     sum_{i: A_i=c} T_i^k
 *   after = true  (Buchholz variables): M(c,j) = (1/(N-1)) sum_{i<N: A_i=c} T_{i+1}^k
 *
 * and with norm = true each entry is multiplied by N/count_c, which turns the
 * contribution into the class-conditional moment normalized so that
 * M_k = sum_c M(c,k) * p_c with p_c = count_c/N the class probabilities of
 * mtrace_pc.
 *
 * REFERENCE DEFECT (JAR): Mtrace_moment.java divides the class sum by
 * count_c ALREADY in the unnormalized branch, so its norm = 0 output is
 * MATLAB's norm = 1 output, and its norm = 1 output is that value multiplied
 * by N/count_c a SECOND time -- a quantity that is neither of the two
 * documented normalizations and that diverges as the class becomes rare.
 * Consequently sum_c M_java(c,k) is not the class-independent moment for
 * either flag value. MATLAB is the reference and is implemented here.
 * Measured on T = 1,2,3,4,5 with A = 1,2,1,2,1 and order 1: MATLAB returns
 * (9/5, 6/5) unnormalized and (3, 3) normalized, the JAR returns (3, 3) and
 * (5, 15/2). Only the MATLAB pair sums to the trace mean 3.
 *
 * A second, smaller divergence: for after = true the JAR normalizes by
 * (N-1)/count_c whereas MATLAB uses N/count_c (its `length(T-1)` is the
 * length of the elementwise T-1, i.e. N). MATLAB's factor is the one
 * consistent with p_c = count_c/N, so it is kept.
 *
 * ARITHMETIC: sums of integer powers and a division, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/**
 * @param Tv     inter-arrival times
 * @param A      class labels
 * @param orders moment orders
 * @param after  false for Horvath variables, true for Buchholz variables
 * @param norm   true to normalize by N/count_c
 * @return (C x |orders|) matrix, row c for the c-th smallest label
 */
template <class T>
Matrix<T> mtrace_moment(const std::vector<T>& Tv, const std::vector<int>& A,
                        const std::vector<unsigned>& orders, bool after = false,
                        bool norm = false) {
    detail::require_marked(Tv, A, "mtrace_moment");
    const std::vector<int> marks = detail::unique_labels(A);
    const std::size_t C = marks.size();
    const std::size_t N = Tv.size();
    if (after && N < 2) throw InputError("mtrace_moment: the Buchholz form needs N >= 2");

    Matrix<T> M(C, orders.size(), num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < orders.size(); ++j) {
        const unsigned k = orders[j];
        for (std::size_t c = 0; c < C; ++c) {
            T sum = num_traits<T>::from_int(0);
            long count = 0;
            const std::size_t last = after ? N - 1 : N;
            for (std::size_t i = 0; i < last; ++i) {
                if (A[i] != marks[c]) continue;
                sum += num_pow_int(Tv[after ? i + 1 : i], k);
                ++count;
            }
            // unnormalized-entry mean rationale: see _kb/03-api-layer.md (cpp port notes: trace)
            T val = sum / num_traits<T>::from_int(static_cast<long>(last));
            if (norm && count > 0)
                val *= num_traits<T>::from_int(static_cast<long>(N)) /
                       num_traits<T>::from_int(count);
            M(c, j) = val;
        }
    }
    return M;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_MOMENT_H
