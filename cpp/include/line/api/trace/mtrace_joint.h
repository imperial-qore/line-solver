/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_JOINT_H
#define LINE_API_TRACE_MTRACE_JOINT_H

/**
 * Class-dependent joint moments of a marked trace,
 *
 *   JM(a) = (1/N_a) sum_{j : A_{j+1} = a} T_j^{i1} T_{j+1}^{i2},
 *
 * the empirical estimate of E[(X_j)^{i1} (X_{j+1})^{i2}] conditioned on the
 * middle event being of class a; the sum runs over the interior events, so
 * the first and last event of the trace are excluded.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_joint.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_joint.java (identical,
 * including the index range and the normalization by the interior class
 * count).
 *
 * Both references index classes by the RAW label 1..max(A) rather than by
 * position in unique(A), so a label alphabet with gaps yields zero rows. That
 * convention is preserved here: entry a-1 of the result refers to label a.
 *
 * ARITHMETIC: sums of products of integer powers, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of mtrace_joint. */
template <class T>
struct MtraceJointResult {
    std::vector<T> jm;       ///< one moment per label 1..max(A)
    std::vector<long> count;  ///< interior events of each label; 0 = undefined
};

/**
 * @param Tv inter-event times
 * @param A  class labels, positive
 * @param i1 exponent of the first interval
 * @param i2 exponent of the second interval
 */
template <class T>
MtraceJointResult<T> mtrace_joint(const std::vector<T>& Tv, const std::vector<int>& A,
                                  unsigned i1, unsigned i2) {
    detail::require_marked(Tv, A, "mtrace_joint");
    const std::size_t N = A.size();
    if (N < 3) throw InputError("mtrace_joint: at least three events are required");
    int maxlab = 0;
    for (std::size_t k = 0; k < N; ++k) {
        if (A[k] < 0) throw InputError("mtrace_joint: class labels must be nonnegative");
        if (A[k] > maxlab) maxlab = A[k];
    }
    MtraceJointResult<T> out;
    out.jm.assign(static_cast<std::size_t>(maxlab), num_traits<T>::from_int(0));
    out.count.assign(static_cast<std::size_t>(maxlab), 0);

    for (int a = 1; a <= maxlab; ++a) {
        const std::size_t ai = static_cast<std::size_t>(a - 1);
        for (std::size_t j = 1; j + 1 < N; ++j)
            if (A[j] == a) out.count[ai] += 1;
        T tmp = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j + 2 < N; ++j)
            if (A[j + 1] == a) tmp += num_pow_int(Tv[j], i1) * num_pow_int(Tv[j + 1], i2);
        if (out.count[ai] > 0) out.jm[ai] = tmp / num_traits<T>::from_int(out.count[ai]);
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_JOINT_H
