/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_EREC_H
#define LINE_API_CACHE_EREC_H

/**
 * Exact recursive normalizing constant of a multi-list cache model.
 *
 * Templated port of matlab/src/api/cache/cache_erec.m, cross-checked against
 * jar/src/main/java/jline/api/cache/Cache_erec.java.
 *
 * The cache holds h lists of capacities m(1..h) filled from n items; the
 * steady-state distribution of the list-based replacement model is
 * proportional to prod over placed items of gamma(item,list), and E(gamma,m)
 * is the sum of that product over every admissible placement, weighted by the
 * list multiplicities. The recursion peels off item k:
 *
 *   E_k(m) = E_{k-1}(m) + sum_j gamma(k,j) m(j) E_{k-1}(m - e_j),
 *
 * with E(0) = 1, E(m) = 0 whenever sum(m) exceeds the number of remaining
 * items or any capacity is negative, and E_1(e_j) = gamma(1,j).
 *
 * Every step is a multiplication and an addition, so the exact instantiation
 * returns E as a rational with no rounding at all: for rational access factors
 * that is the true normalizing constant, which is what makes it usable as the
 * oracle for the approximate members of the family (cache_spm, cache_xi_fp).
 *
 * DIVERGENCE, MATLAB vs JAR: the MATLAB entry point recurses from
 * k = length(gamma), and MATLAB's length() on an (n x h) matrix is max(n,h).
 * That is the item count n only while n >= h; for a cache with more lists than
 * items MATLAB starts the recursion at h and reads gamma rows that do not
 * exist. The JAR uses gamma.getNumRows() and is right. This port follows the
 * JAR (n = number of rows).
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

namespace detail {

/**
 * gamma with row i deleted, i.e. MATLAB's gamma(setdiff(1:n,i),:). Shared by
 * every member of the family that conditions on "item i is absent".
 */
template <class T>
Matrix<T> gamma_without_row(const Matrix<T>& gamma, std::size_t i) {
    if (gamma.rows() == 0) return gamma;
    Matrix<T> g(gamma.rows() - 1, gamma.cols());
    std::size_t r = 0;
    for (std::size_t a = 0; a < gamma.rows(); ++a) {
        if (a == i) continue;
        for (std::size_t b = 0; b < gamma.cols(); ++b) g(r, b) = gamma(a, b);
        ++r;
    }
    return g;
}

/** E over items 1..k of the (n x h) access-factor matrix. */
template <class T>
T cache_erec_aux(const Matrix<T>& gamma, const std::vector<int>& m, int k) {
    long mt = 0;
    int mmin = 0;
    bool first = true;
    for (int v : m) {
        mt += v;
        if (first || v < mmin) mmin = v;
        first = false;
    }
    if (mt == 0) return num_traits<T>::from_int(1);
    if (mt > k || mmin < 0) return num_traits<T>::from_int(0);

    const std::size_t h = m.size();
    if (k == 1 && mt == 1) {
        for (std::size_t j = 0; j < h; ++j)
            if (m[j] != 0) return gamma(0, j);
        return num_traits<T>::from_int(0);
    }

    T E = cache_erec_aux(gamma, m, k - 1);
    for (std::size_t j = 0; j < h; ++j) {
        if (m[j] > 0) {
            std::vector<int> mj = m;
            mj[j] -= 1;
            E += gamma(static_cast<std::size_t>(k - 1), j) *
                 num_traits<T>::from_int(static_cast<long>(m[j])) *
                 cache_erec_aux(gamma, mj, k - 1);
        }
    }
    return E;
}

/**
 * E(m,k) over the (residual capacity, residual cost cap) lattice, the
 * cost-capped recursion of Casale-Gast, IEEE/ACM Trans. Networking 29(2),
 * 2021, Sec. IX:
 *
 *   E(m,k) = E_i(m,k) + sum_j m_j gamma(i,j) E_i(m - e_j, k - sigma_i e_j),
 *
 * with the extra boundary E(m,k) = 0 whenever a residual cap is negative. The
 * plain recursion cannot carry the cap, because the k - sigma_i e_j argument
 * couples item sizes into the recursion graph, so this is evaluated as a
 * dynamic program over the product lattice of residual capacities and residual
 * caps: O(n prod_j (m_j+1)(k_j+1)) time and space.
 */
template <class T>
T cache_erec_cost(const Matrix<T>& gamma, const std::vector<int>& m,
                  const std::vector<int>& sigma, const std::vector<int>& k) {
    const std::size_t n = gamma.rows();
    const std::size_t h = m.size();
    if (sigma.size() != n)
        throw InputError("cache_erec: the item size vector must have one entry per item");
    if (k.size() != h)
        throw InputError("cache_erec: the cost cap vector must have one entry per cache list");
    long mt = 0;
    for (std::size_t j = 0; j < h; ++j) {
        if (m[j] < 0 || k[j] < 0) return num_traits<T>::from_int(0);
        mt += m[j];
    }
    for (std::size_t i = 0; i < n; ++i)
        if (sigma[i] <= 0)
            throw InputError("cache_erec: item sizes must be positive integers");
    if (mt > static_cast<long>(n)) return num_traits<T>::from_int(0);
    if (mt == 0) return num_traits<T>::from_int(1);

    std::vector<std::size_t> dims(2 * h), stride(2 * h);
    for (std::size_t j = 0; j < h; ++j) {
        dims[j] = static_cast<std::size_t>(m[j]) + 1;
        dims[h + j] = static_cast<std::size_t>(k[j]) + 1;
    }
    std::size_t size = 1;
    for (std::size_t d = 0; d < 2 * h; ++d) {
        stride[d] = size;
        size *= dims[d];
        if (size > 10000000u)
            throw InputError(
                "cache_erec: the cost-constrained lattice exceeds the exact method "
                "limit; use the sampling method");
    }
    std::vector<T> F(size, num_traits<T>::from_int(0));
    std::vector<T> Fprev(size, num_traits<T>::from_int(0));
    std::vector<std::size_t> sub(2 * h, 0);
    for (std::size_t idx = 0; idx < size; ++idx) {
        std::size_t rem = idx, mc = 0;
        for (std::size_t d = 2 * h; d-- > 0;) {
            sub[d] = rem / stride[d];
            rem -= sub[d] * stride[d];
        }
        for (std::size_t j = 0; j < h; ++j) mc += sub[j];
        if (mc == 0) F[idx] = num_traits<T>::from_int(1);
    }
    for (std::size_t t = 0; t < n; ++t) {
        Fprev = F;
        for (std::size_t idx = 0; idx < size; ++idx) {
            std::size_t rem = idx, mc = 0;
            for (std::size_t d = 2 * h; d-- > 0;) {
                sub[d] = rem / stride[d];
                rem -= sub[d] * stride[d];
            }
            for (std::size_t j = 0; j < h; ++j) mc += sub[j];
            if (mc > t + 1) {
                F[idx] = num_traits<T>::from_int(0);
                continue;
            }
            T val = Fprev[idx];
            for (std::size_t j = 0; j < h; ++j) {
                const std::size_t mj = sub[j];
                const std::size_t kj = sub[h + j];
                const std::size_t si = static_cast<std::size_t>(sigma[t]);
                if (mj > 0 && kj >= si && !(gamma(t, j) == num_traits<T>::from_int(0))) {
                    val += gamma(t, j) * num_traits<T>::from_int(static_cast<long>(mj)) *
                           Fprev[idx - stride[j] - si * stride[h + j]];
                }
            }
            F[idx] = val;
        }
    }
    return F[size - 1];
}

}  // namespace detail

/**
 * @param gamma (n x h) access factors, item by list
 * @param m     (h) list capacities
 * @return the normalizing constant E
 */
template <class T>
T cache_erec(const Matrix<T>& gamma, const std::vector<int>& m) {
    if (gamma.empty()) {
        // empty-cache base case rationale: see _kb/09-ldes-and-cache.md (cpp port notes)
        return detail::cache_erec_aux(gamma, m, 0);
    }
    if (gamma.cols() != m.size())
        throw InputError("cache_erec: gamma and m disagree on the number of lists");
    return detail::cache_erec_aux(gamma, m, static_cast<int>(gamma.rows()));
}

/**
 * The normalizing constant under per-list storage cost caps. An EMPTY sigma or
 * k selects the unconstrained recursion, so a caller can pass the model's
 * (possibly absent) sizes through unconditionally.
 *
 * @param gamma (n x h) access factors, item by list
 * @param m     (h) list capacities
 * @param sigma (n) per-item storage costs, positive integers
 * @param k     (h) per-list storage cost caps, non-negative integers
 * @return the normalizing constant E(m,k)
 */
template <class T>
T cache_erec(const Matrix<T>& gamma, const std::vector<int>& m,
             const std::vector<int>& sigma, const std::vector<int>& k) {
    if (sigma.empty() || k.empty()) return cache_erec(gamma, m);
    if (gamma.cols() != m.size())
        throw InputError("cache_erec: gamma and m disagree on the number of lists");
    return detail::cache_erec_cost(gamma, m, sigma, k);
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_EREC_H
