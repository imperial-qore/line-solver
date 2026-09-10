/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_PC_H
#define LINE_API_TRACE_MTRACE_PC_H

/**
 * Class probabilities of a marked trace, p_c = count_c / N.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_pc.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_pc.java (identical).
 *
 * The inter-arrival times are not used; both references take them only for
 * signature orthogonality with the rest of the mtrace family, so they are not
 * a parameter here.
 *
 * ARITHMETIC: counts over the sample size, exact in Rational. The entries sum
 * to exactly 1 in the exact backends.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** @param A class labels; @return one probability per label of unique(A). */
template <class T>
std::vector<T> mtrace_pc(const std::vector<int>& A) {
    if (A.empty()) throw InputError("mtrace_pc: the trace is empty");
    const std::vector<int> labels = detail::unique_labels(A);
    std::vector<T> pc;
    pc.reserve(labels.size());
    const T N = num_traits<T>::from_int(static_cast<long>(A.size()));
    for (std::size_t i = 0; i < labels.size(); ++i) {
        long count = 0;
        for (std::size_t k = 0; k < A.size(); ++k)
            if (A[k] == labels[i]) ++count;
        pc.push_back(num_traits<T>::from_int(count) / N);
    }
    return pc;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_PC_H
