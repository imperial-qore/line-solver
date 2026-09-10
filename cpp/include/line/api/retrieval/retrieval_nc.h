/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_RETRIEVAL_RETRIEVAL_NC_H
#define LINE_API_RETRIEVAL_RETRIEVAL_NC_H

/**
 * Exact normalizing constant E(v,m) of a delayed-hit (list-based) cache.
 *
 * Templated port of matlab/src/api/retrieval/retrieval_nc.m, cross-checked
 * against jar/src/main/java/jline/api/retrieval/Retrieval_nc.java.
 *
 * The cache holds h lists of capacities m(1..h) filled from n items, and a
 * missed item is fetched by a retrieval system with one infinite-server
 * station (index s=0) and r processor-sharing stations (s=1..r). The exact
 * recurrence peels off item k:
 *
 *   E(v,m) = (1 + lambda_k eta_{0,k}) E_k(v,m)
 *          + sum_{s=1}^r lambda_k eta_{s,k} (v_s+1) E_k(v+1_s, m)
 *          + sum_{j=1}^h m_j gamma_{k,j} E_k(v, m-1_j)
 *
 * with E_k the constant of the system without item k, E = 1 once no item is
 * left, and E = 0 as soon as sum_j m_j exceeds the number of remaining items
 * or some m_j is negative. The plain constant is E(m) = E(0,m); the vector v
 * carries the moment order at each PS station, so E(1_s,m) is what the
 * delayed-hit probability at station s needs.
 *
 * ARITHMETIC: every step is a multiplication and an addition of the inputs,
 * so this is a finite field computation and the exact instantiation returns
 * the true constant with no rounding at all. That is what makes it the oracle
 * for retrieval_fpi, which only approximates the same quantities.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace retrieval {

namespace detail {

/** E over items 1..k of the n-item system. */
template <class T>
T retrieval_nc_aux(std::vector<int>& v, std::vector<int>& m, const std::vector<T>& lambda,
                   const Matrix<T>& eta, const Matrix<T>& gamma, int k) {
    const std::size_t r = v.size();
    const std::size_t h = m.size();
    long msum = 0;
    int mmin = 0;
    bool first = true;
    for (int x : m) {
        msum += x;
        if (first || x < mmin) mmin = x;
        first = false;
    }
    if (msum > k || mmin < 0) return num_traits<T>::from_int(0);
    if (k == 0) return num_traits<T>::from_int(1);

    const std::size_t ki = static_cast<std::size_t>(k - 1);

    // item k out of the system, or fetched at the IS station s = 0
    T E = (num_traits<T>::from_int(1) + lambda[ki] * eta(ki, 0)) *
          retrieval_nc_aux(v, m, lambda, eta, gamma, k - 1);

    // item k fetched at PS station s = 1..r
    for (std::size_t s = 0; s < r; ++s) {
        const int vs = v[s];
        v[s] = vs + 1;
        const T sub = retrieval_nc_aux(v, m, lambda, eta, gamma, k - 1);
        v[s] = vs;
        E += lambda[ki] * eta(ki, s + 1) * num_traits<T>::from_int(static_cast<long>(vs) + 1) * sub;
    }

    // item k stored in cache list j = 1..h
    for (std::size_t j = 0; j < h; ++j) {
        if (m[j] > 0) {
            const int mj = m[j];
            m[j] = mj - 1;
            const T sub = retrieval_nc_aux(v, m, lambda, eta, gamma, k - 1);
            m[j] = mj;
            E += gamma(ki, j) * num_traits<T>::from_int(static_cast<long>(mj)) * sub;
        }
    }
    return E;
}

}  // namespace detail

/**
 * @param v      (r) moment order at each PS station; all zeros for the plain constant
 * @param m      (h) cache list capacities
 * @param lambda (n) per-item arrival rates
 * @param eta    (n x (r+1)) fetching demands, column 0 = IS station, columns 1..r = PS
 * @param gamma  (n x h) access factors
 * @return E(v,m)
 */
template <class T>
T retrieval_nc(const std::vector<int>& v, const std::vector<int>& m, const std::vector<T>& lambda,
               const Matrix<T>& eta, const Matrix<T>& gamma) {
    const std::size_t n = lambda.size();
    if (eta.rows() != n) throw InputError("retrieval_nc: eta and lambda disagree on the item count");
    if (gamma.rows() != n)
        throw InputError("retrieval_nc: gamma and lambda disagree on the item count");
    if (eta.cols() != v.size() + 1)
        throw InputError("retrieval_nc: eta and v disagree on the number of PS stations");
    if (gamma.cols() != m.size())
        throw InputError("retrieval_nc: gamma and m disagree on the number of lists");
    std::vector<int> vv = v;
    std::vector<int> mm = m;
    return detail::retrieval_nc_aux(vv, mm, lambda, eta, gamma, static_cast<int>(n));
}

}  // namespace retrieval
}  // namespace line

#endif  // LINE_API_RETRIEVAL_RETRIEVAL_NC_H
