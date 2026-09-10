/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_CHAR_MAX_DISCRETE_H
#define LINE_API_FJ_CHAR_MAX_DISCRETE_H

/**
 * Characteristic maximum of a lattice random variable.
 *
 * Templated port of matlab/src/api/fj/fj_char_max_discrete.m.
 *
 * With m_K the smallest integer at which P(X > m_K) <= 1/K,
 *
 *   M_K = m_K + K sum_{k >= m_K} P(X > k),
 *
 * which upper bounds the expected maximum of K i.i.d. copies at O(1) instead of
 * the alternating binomial sum. Two lattice laws close the tail sum:
 *
 *   geometric, P(X = k) = (1-p) p^k:
 *     m_K = ceil(-ln K / ln p),  M_K = m_K + K p^(m_K+1)/(1-p),
 *     exact E[Y_K] = sum_k C(K,k) (-1)^(k+1) p^k/(1-p^k);
 *   Poisson:
 *     M_K = m_K (1 - K P(X > m_K)) + K theta P(X > m_K - 1).
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** The lattice laws for which the characteristic maximum is closed. */
enum class FJDiscreteDist { Geometric, Poisson };

/** [MK, mK, exact] of fj_char_max_discrete. */
template <class T>
struct FJCharMaxDiscreteResult {
    T MK;
    unsigned mK;
    T exact;
};

/**
 * @param K    number of i.i.d. copies, K >= 1
 * @param dist the lattice law
 * @param par  p in (0,1) for the geometric, theta > 0 for the Poisson
 * @return     characteristic maximum, its threshold, and the exact maximum
 */
template <class T>
FJCharMaxDiscreteResult<T> fj_char_max_discrete(unsigned K, FJDiscreteDist dist, const T& par) {
    detail::require_positive_K(K, "fj_char_max_discrete");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    FJCharMaxDiscreteResult<T> out;

    if (dist == FJDiscreteDist::Geometric) {
        const T p = par;
        if (!(p > zero) || !(p < one))
            throw InputError("fj_char_max_discrete: the geometric parameter must lie in (0,1)");
        // Smallest integer k with p^k <= 1/K
        const double mkd = -std::log(static_cast<double>(K)) /
                           std::log(num_traits<T>::to_double(p));
        long mk = static_cast<long>(std::ceil(mkd - 1e-12));
        if (mk < 0) mk = 0;
        out.mK = static_cast<unsigned>(mk);
        T pw = one;
        for (unsigned e = 0; e <= out.mK; ++e) pw *= p;
        out.MK = num_traits<T>::from_int(mk) + num_traits<T>::from_int(static_cast<long>(K)) * pw /
                                                   (one - p);
        // Exact maximum by inclusion-exclusion on the geometric tail
        T acc = zero;
        T pk = one;
        for (unsigned k = 1; k <= K; ++k) {
            pk *= p;
            const T term = detail::fj_binom<T>(K, k) * pk / (one - pk);
            if (k % 2 == 1) acc += term; else acc -= term;
        }
        out.exact = acc;
        return out;
    }

    const T theta = par;
    if (!(theta > zero))
        throw InputError("fj_char_max_discrete: the Poisson mean must be positive");
    const double th = num_traits<T>::to_double(theta);
    const std::size_t kmax = static_cast<std::size_t>(std::ceil(th + 12 * std::sqrt(th) + 40));
    std::vector<T> pmf(kmax + 1), cdf(kmax + 1), tail(kmax + 1);
    T term = detail::num_exp<T>(-theta);
    T acc = zero;
    for (std::size_t k = 0; k <= kmax; ++k) {
        if (k > 0) term = term * theta / num_traits<T>::from_int(static_cast<long>(k));
        pmf[k] = term;
        acc += term;
        cdf[k] = (acc > one) ? one : acc;
        tail[k] = one - cdf[k];
    }
    std::size_t mk = kmax + 1;
    const T thr = one / num_traits<T>::from_int(static_cast<long>(K));
    for (std::size_t k = 0; k <= kmax; ++k)
        if (tail[k] <= thr) { mk = k; break; }
    if (mk > kmax)
        throw NumericError("fj_char_max_discrete: the Poisson lattice truncation never reached 1/K");
    out.mK = static_cast<unsigned>(mk);
    const T tail_prev = (mk == 0) ? one : tail[mk - 1];
    out.MK = num_traits<T>::from_int(static_cast<long>(mk)) *
                 (one - num_traits<T>::from_int(static_cast<long>(K)) * tail[mk]) +
             num_traits<T>::from_int(static_cast<long>(K)) * theta * tail_prev;
    // Exact maximum as the sum over the lattice of 1 - F(k)^K
    T ex = zero;
    for (std::size_t k = 0; k <= kmax; ++k) {
        T pw = one;
        for (unsigned e = 0; e < K; ++e) pw *= cdf[k];
        ex += one - pw;
    }
    out.exact = ex;
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_CHAR_MAX_DISCRETE_H
