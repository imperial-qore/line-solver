/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_VAR_H
#define LINE_API_TRACE_MTRACE_VAR_H

/**
 * Per-class variance of a marked trace.
 *
 * Port of `mtrace_var` in python/line_solver/api/trace/trace_analysis.py.
 * PYTHON-ONLY: no MATLAB or JAR twin. The unmarked twin is `trace_var.h`.
 *
 * IT IS THE POPULATION VARIANCE, `E[X^2] - E[X]^2`, not the sample one: the
 * reference divides by the class count, not by count minus one. That matters
 * when a class is short, and the two differ by the factor n/(n-1).
 *
 * A CLASS WITH ONE OR NO SAMPLE GETS NaN, NOT ZERO. A single observation has no
 * dispersion to report, and zero would be indistinguishable from a class whose
 * samples happen to be identical -- which is a real and different finding. The
 * reference returns NaN and so does this; a caller aggregating these must test.
 *
 * ARITHMETIC: field.
 */

#include <cstddef>
#include <limits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/**
 * @param tv     the trace values
 * @param ntypes number of classes
 * @param types  0-based class label of each value
 * @return       (ntypes) per-class variance, NaN where a class has under two
 */
template <class T>
std::vector<T> mtrace_var(const std::vector<T>& tv, std::size_t ntypes,
                          const std::vector<int>& types) {
    if (tv.size() != types.size())
        throw InputError("mtrace_var: the trace and its labels must agree in length");
    const T zero = num_traits<T>::from_int(0);
    const T nan = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());

    std::vector<T> out(ntypes, nan);
    for (std::size_t c = 0; c < ntypes; ++c) {
        std::size_t cnt = 0;
        T s1 = zero, s2 = zero;
        for (std::size_t i = 0; i < tv.size(); ++i) {
            if (types[i] != static_cast<int>(c)) continue;
            ++cnt;
            s1 += tv[i];
            s2 += tv[i] * tv[i];
        }
        if (cnt <= 1) continue;  // no dispersion to report; NaN, not zero
        const T nd = num_traits<T>::from_int(static_cast<long>(cnt));
        const T e1 = T(s1 / nd), e2 = T(s2 / nd);
        out[c] = T(e2 - e1 * e1);
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_VAR_H
