/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_TYPES_H
#define LINE_API_TRACE_TRACE_TYPES_H

/**
 * Shared declarations for the empirical trace statistics domain.
 *
 * Templated port of the kpctoolbox trace primitives
 * (matlab/lib/kpctoolbox/trace/, matlab/lib/m3a/m3a/mtrace/) and of
 * jar/src/main/java/jline/api/trace/.
 *
 * CONVENTIONS
 *  - A single-class trace is a vector of inter-arrival (or service) times,
 *    passed as std::vector<T>.
 *  - A marked (multi-class) trace is the pair (T, A) with A a vector of
 *    integer class labels of the same length. The label alphabet is whatever
 *    distinct values occur in A, taken in increasing order; row c of every
 *    returned matrix refers to the c-th smallest label, exactly as MATLAB's
 *    `unique(A)` and the JAR's TreeSet do. Two JAR functions instead index by
 *    the raw label value (Mtrace_cov, Mtrace_joint use 1..max(A)); those
 *    divergences are documented on the individual headers.
 *
 * ARITHMETIC
 *  Sample moments, autocovariances, class-transition frequencies, counting
 *  processes and order statistics are sums, products and divisions of the
 *  data, so the whole of this domain except trace_skew, trace_summary and
 *  trace_gamma's residual comparison is a field computation and is
 *  instantiated for double, Rational and Real50. Only the functions that take
 *  a square root (skewness normalization, standard deviations in the summary)
 *  carry a has_transcendental static_assert.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

namespace detail {

/** exp(v), resolved by ADL. */
template <class T>
inline T num_exp(const T& v) {
    using std::exp;
    return exp(v);
}

/** log(v), resolved by ADL. */
template <class T>
inline T num_log(const T& v) {
    using std::log;
    return log(v);
}

/** sqrt(v), resolved by ADL. */
template <class T>
inline T num_sqrt(const T& v) {
    using std::sqrt;
    return sqrt(v);
}

/** base^exponent for a real-valued exponent, resolved by ADL. */
template <class T>
inline T num_pow(const T& base, const T& exponent) {
    using std::pow;
    return pow(base, exponent);
}

/** The distinct class labels in increasing order, i.e. MATLAB's unique(A). */
inline std::vector<int> unique_labels(const std::vector<int>& A) {
    const std::set<int> s(A.begin(), A.end());
    return std::vector<int>(s.begin(), s.end());
}

/** Throws unless the trace is non-empty. */
template <class T>
inline void require_nonempty(const std::vector<T>& S, const char* who) {
    if (S.empty()) throw InputError(std::string(who) + ": the trace is empty");
}

/** Throws unless times and labels have the same length. */
template <class T>
inline void require_marked(const std::vector<T>& S, const std::vector<int>& A, const char* who) {
    require_nonempty(S, who);
    if (S.size() != A.size())
        throw InputError(std::string(who) + ": times and class labels have different lengths");
}

/**
 * Cumulative sums with a leading zero: cs[k] is the sum of the first k
 * entries, so cs has size n+1 and cs[0] = 0. This is MATLAB's [0; cumsum(S)]
 * and lets the 1-based index arithmetic of the reference sources be
 * transcribed literally.
 */
template <class T>
inline std::vector<T> cumsum0(const std::vector<T>& S) {
    std::vector<T> cs(S.size() + 1, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < S.size(); ++i) cs[i + 1] = cs[i] + S[i];
    return cs;
}

}  // namespace detail
}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_TYPES_H
