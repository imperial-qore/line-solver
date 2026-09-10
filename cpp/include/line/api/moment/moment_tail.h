/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_TAIL_H
#define LINE_API_MOMENT_MOMENT_TAIL_H

/**
 * Binomial moments from tail moments and the inverse.
 *
 * Templated port of matlab/src/api/moment/moment_binomial_from_tail.m and
 * moment_tail_from_binomial.m. These are the only edges of the house of moments
 * whose table is UPPER triangular, so both directions read the WHOLE remaining
 * sequence: b_i depends on t_i,...,t_n. Truncating the input therefore
 * truncates the information, not just the output order.
 *
 * Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 */

#include <vector>

#include "line/api/moment/moment_housematrix.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace moment {

/** b_i = sum_{k>=i} C(k-1,i-1) t_k. */
template <class T>
std::vector<T> moment_binomial_from_tail(const std::vector<T>& t) {
    if (t.empty()) throw InputError("moment_binomial_from_tail: t must be nonempty");
    return moment_apply_full<T>(
        moment_housematrix<T>(MomentEdge::BinomialFromTail, static_cast<int>(t.size()) - 1), t);
}

/** t_i = sum_{k>=i} (-1)^(k-i) C(k-1,i-1) b_k. */
template <class T>
std::vector<T> moment_tail_from_binomial(const std::vector<T>& b) {
    if (b.empty()) throw InputError("moment_tail_from_binomial: b must be nonempty");
    return moment_apply_full<T>(
        moment_housematrix<T>(MomentEdge::TailFromBinomial, static_cast<int>(b.size()) - 1), b);
}

}  // namespace moment
}  // namespace line

#endif
