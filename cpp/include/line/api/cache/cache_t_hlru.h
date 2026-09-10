/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_T_HLRU_H
#define LINE_API_CACHE_T_HLRU_H

/**
 * Characteristic times of the h-LRU / LRU(m) TTL approximation.
 *
 * Templated port of matlab/src/api/cache/cache_t_hlru.m, cross-checked against
 * jar/src/main/java/jline/api/cache/Cache_t_hlru.java.
 *
 * Under the characteristic-time (Che) approximation each list l of an h-LRU
 * cache has a time T(l) such that an item of request rate lam is evicted from
 * list l if it is not requested within T(l). The level of an item is then a
 * birth-death chain with e(l) = exp(-lam T(l)), giving unnormalized weights
 *
 *   w(0) = 1,   w(l) = w(l-1) (1 - e(l)) / e(l),
 *
 * and pi_l = w(l)/sum(w). The times are fixed by sum_k pi_l(k;T) = m(l), one
 * equation per list. sum_k pi_l is increasing in T(l), so each equation is
 * solved by bisection (bracketed by doubling the upper end until the occupancy
 * reaches the capacity, then 100 halvings), swept Gauss-Seidel over the lists
 * until the times move by less than 1e-8 relative.
 *
 * ARITHMETIC: exp is required and the answer is defined by a bisection
 * tolerance, so this needs transcendental arithmetic.
 *
 * Reference: Gast and Van Houdt, SIGMETRICS 2015. For h = 1 the fixed point
 * reduces exactly to the Che approximation of LRU.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

namespace detail {

/** LINE's GlobalConstants.Zero and FineTol, as set by lineStart. */
template <class T>
T hlru_zero() {
    return num_traits<T>::from_double(1e-14);
}
template <class T>
T hlru_finetol() {
    return num_traits<T>::from_double(1e-8);
}

/**
 * Birth-death level probabilities pi(k,l), l = 0..h, for request rates lam and
 * characteristic times T.
 */
template <class T>
Matrix<T> hlru_levelprobs(const std::vector<T>& lam, const std::vector<T>& T_) {
    using std::exp;
    const std::size_t n = lam.size();
    const std::size_t h = T_.size();
    const T one = num_traits<T>::from_int(1);
    const T zeroTol = hlru_zero<T>();
    Matrix<T> P(n, h + 1, num_traits<T>::from_int(0));
    for (std::size_t k = 0; k < n; ++k) {
        std::vector<T> w(h + 1, one);
        for (std::size_t l = 0; l < h; ++l) {
            const T e = exp(T(-lam[k] * T_[l]));
            const T den = e > zeroTol ? e : zeroTol;
            w[l + 1] = w[l] * (one - e) / den;
        }
        T s = num_traits<T>::from_int(0);
        for (std::size_t l = 0; l <= h; ++l) s += w[l];
        for (std::size_t l = 0; l <= h; ++l) P(k, l) = w[l] / s;
    }
    return P;
}

/** Total occupancy of list l when its characteristic time is set to Tl. */
template <class T>
T hlru_occ(const std::vector<T>& lam, std::vector<T> T_, std::size_t l, const T& Tl) {
    T_[l] = Tl;
    const Matrix<T> P = hlru_levelprobs(lam, T_);
    T occ = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < lam.size(); ++k) occ += P(k, l + 1);
    return occ;
}

/** The shared Gauss-Seidel + bisection fixed point for the times. */
template <class T>
std::vector<T> hlru_solve_times(const std::vector<T>& lam, const std::vector<int>& m) {
    const std::size_t n = lam.size();
    const std::size_t h = m.size();
    if (n == 0) throw InputError("cache_t_hlru: no items");
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);

    T mean = zero;
    for (std::size_t k = 0; k < n; ++k) mean += lam[k];
    mean /= num_traits<T>::from_int(static_cast<long>(n));
    const T scale = mean > hlru_finetol<T>() ? mean : hlru_finetol<T>();
    const T T0 = num_traits<T>::from_int(1) / scale;
    const T hicap = num_traits<T>::from_double(1e12);

    std::vector<T> Tv(h, T0);
    for (int sweep = 0; sweep < 200; ++sweep) {
        const std::vector<T> Told = Tv;
        for (std::size_t l = 0; l < h; ++l) {
            const T ml = num_traits<T>::from_int(static_cast<long>(m[l]));
            T lo = zero;
            T hi = Tv[l] > T0 ? Tv[l] : T0;
            while (hlru_occ(lam, Tv, l, hi) < ml && hi < hicap) hi = two * hi;
            for (int b = 0; b < 100; ++b) {
                const T mid = (lo + hi) / two;
                if (hlru_occ(lam, Tv, l, mid) < ml)
                    lo = mid;
                else
                    hi = mid;
            }
            Tv[l] = (lo + hi) / two;
        }
        T rel = zero;
        for (std::size_t l = 0; l < h; ++l) {
            const T den = Told[l] > hlru_zero<T>() ? Told[l] : hlru_zero<T>();
            const T d = num_abs(T(Tv[l] - Told[l])) / den;
            if (d > rel) rel = d;
        }
        if (rel < hlru_finetol<T>()) break;
    }
    return Tv;
}

}  // namespace detail

/**
 * @param gamma (n x 1) per-item request rates; an (n x h) matrix is accepted
 *              for backward compatibility and its first column is used, as in
 *              the MATLAB reference
 * @param m     (h) list capacities
 * @return (h) characteristic times
 */
template <class T>
std::vector<T> cache_t_hlru(const Matrix<T>& gamma, const std::vector<int>& m) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_t_hlru requires transcendental arithmetic");
    if (gamma.cols() == 0) throw InputError("cache_t_hlru: empty rate matrix");
    std::vector<T> lam(gamma.rows());
    for (std::size_t k = 0; k < gamma.rows(); ++k) lam[k] = gamma(k, 0);
    return detail::hlru_solve_times(lam, m);
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_T_HLRU_H
