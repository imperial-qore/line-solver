/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_ORDSTAT_EXP_H
#define LINE_API_FJ_ORDSTAT_EXP_H

/**
 * Mean of the k-th smallest of n independent EXPONENTIAL branch completion
 * times, i.e. the instant a k-of-n (quorum) join fires. k = n is the ordinary
 * AND-join, the maximum, and k = 1 the minimum.
 *
 * Templated port of matlab/src/api/fj/fj_ordstat_exp.m, cross-checked against
 * jar/src/main/java/jline/api/fj/FJ_ordstat_exp.java.
 *
 * With lambda_i = 1/ri(i) and m = n-k stragglers allowed,
 *
 *   E[X_(k)] = sum_{j=m+1..n} (-1)^(j-m-1) C(j-1,m) e_j,
 *   e_j      = sum_{|S|=j} 1 / sum_{i in S} lambda_i
 *
 * the inclusion-exclusion identity for the order statistics of independent
 * exponentials. At m = 0 it collapses to sum_j (-1)^(j-1) e_j, the classical
 * expression for the maximum, TERM BY TERM: a full join therefore evaluates
 * exactly as it did before this header existed, which is what keeps the MMT
 * fork-join fixed point bit-identical on a standard join.
 *
 * The exact path is RATIONAL -- reciprocals, sums and integer binomials only --
 * so it holds in an exact field and is not gated on has_transcendental. Only
 * the large-n fallback is, since fj_quorum_moments fits a standard deviation.
 *
 * The sum has 2^n terms and its signs alternate, so it is evaluated exactly
 * only while the branch count is small. Beyond FJ_ORDSTAT_MAX_EXACT branches a
 * genuine quorum (k < n) is evaluated by fj_quorum_moments instead, whose
 * Poisson-binomial recurrence adds no cancellation; a full join keeps the exact
 * path at every n so that no existing result moves.
 *
 * Reference: A. Thomasian, "Analysis of Fork/Join and Related Queueing
 * Systems", ACM Computing Surveys 47(2), Article 17, 2014, Sec. 3 (Eq. 18-19).
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/fj/fj_quorum_moments.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** Branch count above which a genuine quorum leaves the exact alternating sum. */
static const std::size_t FJ_ORDSTAT_MAX_EXACT = 15;

namespace detail {

/**
 * The large-n quorum fallback, compiled only where the arithmetic supports it.
 * An exact field has no sqrt, so there it refuses rather than silently
 * degrading to the alternating sum it was reached to avoid.
 */
template <class T, bool HasTranscendental = num_traits<T>::has_transcendental>
struct OrdstatExpFallback {
    static T eval(const std::vector<T>& ri, std::size_t k) {
        // Branch times are taken as exponential, so the variance is the square
        // of the mean.
        std::vector<T> var(ri.size());
        for (std::size_t i = 0; i < ri.size(); ++i) var[i] = ri[i] * ri[i];
        return fj_quorum_moments(ri, var, k).m;
    }
};

template <class T>
struct OrdstatExpFallback<T, false> {
    static T eval(const std::vector<T>&, std::size_t) {
        throw InputError(
            "fj_ordstat_exp: a quorum over more than FJ_ORDSTAT_MAX_EXACT branches needs the "
            "two-moment fallback, which requires transcendental arithmetic");
    }
};

}  // namespace detail

/**
 * @param ri branch completion time means
 * @param k  quorum, 1 <= k <= ri.size()
 * @return   the mean instant the k-of-n join fires
 */
template <class T>
T fj_ordstat_exp(const std::vector<T>& ri, std::size_t k) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<T> r;
    r.reserve(ri.size());
    for (std::size_t i = 0; i < ri.size(); ++i)
        // an exact field carries no NaN or infinity, so the test goes through
        // the double image, which is where a branch time can be either
        if (std::isfinite(num_traits<T>::to_double(ri[i]))) r.push_back(ri[i]);
    std::size_t n = r.size();
    if (n == 0) return zero;
    if (k < 1 || k > n) throw InputError("fj_ordstat_exp: k must satisfy 1 <= k <= n");
    // A branch of zero mean completes instantly: it never delays the join and it
    // counts toward the quorum at once. Removing it here keeps the reciprocal
    // below finite, which the exact field requires and IEEE only tolerates.
    std::size_t nzero = 0;
    for (std::size_t i = 0; i < n; ++i)
        if (!(r[i] > zero)) ++nzero;
    if (nzero > 0) {
        if (k <= nzero) return zero;
        k -= nzero;
        std::vector<T> pos;
        pos.reserve(n - nzero);
        for (std::size_t i = 0; i < n; ++i)
            if (r[i] > zero) pos.push_back(r[i]);
        r.swap(pos);
        n = r.size();
    }
    if (n == 1) return r[0];

    if (k < n && n > FJ_ORDSTAT_MAX_EXACT) return detail::OrdstatExpFallback<T>::eval(r, k);

    std::vector<T> lambdai(n);
    for (std::size_t i = 0; i < n; ++i) lambdai[i] = one / r[i];

    const std::size_t nstrag = n - k;
    T total = zero;
    // Walk the subsets of each size j by an index vector, so no 2^n table is
    // materialised and the exact field never leaves the loop.
    for (std::size_t j = nstrag + 1; j <= n; ++j) {
        T ej = zero;
        std::vector<std::size_t> idx(j);
        for (std::size_t i = 0; i < j; ++i) idx[i] = i;
        while (true) {
            T s = zero;
            for (std::size_t i = 0; i < j; ++i) s = s + lambdai[idx[i]];
            ej = ej + one / s;
            // next combination in lexicographic order
            std::size_t p = j;
            while (p > 0 && idx[p - 1] == n - j + (p - 1)) --p;
            if (p == 0) break;
            ++idx[p - 1];
            for (std::size_t i = p; i < j; ++i) idx[i] = idx[i - 1] + 1;
        }
        // C(j-1, nstrag) is an integer and stays one in the field
        long long binom = 1;
        for (std::size_t i = 0; i < nstrag; ++i)
            binom = binom * static_cast<long long>(j - 1 - i) / static_cast<long long>(i + 1);
        const T term = num_traits<T>::from_int(static_cast<int>(binom)) * ej;
        if (((j - nstrag - 1) % 2) == 0)
            total = total + term;
        else
            total = total - term;
    }
    return total;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_ORDSTAT_EXP_H
