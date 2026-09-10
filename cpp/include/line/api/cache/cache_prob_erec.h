/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_PROB_EREC_H
#define LINE_API_CACHE_PROB_EREC_H

/**
 * Exact per-item hit and miss probabilities of a multi-list cache.
 *
 * Templated port of matlab/src/api/cache/cache_prob_erec.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_prob_erec.java.
 *
 * The probability that item i sits in list j is the ratio of normalizing
 * constants
 *
 *   prob(i,1+j) = m(j) gamma(i,j) E(gamma without item i, m - e_j) / E(gamma,m)
 *
 * and prob(i,1) = 1 - sum_j prob(i,1+j) is the miss probability. Both
 * constants come from cache_erec, so every operation is a field operation and
 * the exact instantiation returns the probabilities as rationals; in that
 * arithmetic the row sums are exactly one, which is the conservation law the
 * test suite asserts with == rather than a tolerance.
 *
 * The MATLAB reference wraps the miss probability in abs(), which only matters
 * when the hit probabilities sum above one -- impossible in exact arithmetic
 * and a rounding artefact in double. The abs() is kept for bit-compatibility
 * with MATLAB and the JAR.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_erec.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

template <class T>
Matrix<T> cache_prob_erec(const Matrix<T>& gamma, const std::vector<int>& m,
                          const std::vector<int>& sigma, const std::vector<int>& k);

/**
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities
 * @return (n x (h+1)) matrix; column 0 is the miss probability, column 1+j the
 *         probability of a hit in list j
 */
template <class T>
Matrix<T> cache_prob_erec(const Matrix<T>& gamma, const std::vector<int>& m) {
    return cache_prob_erec(gamma, m, std::vector<int>(), std::vector<int>());
}

/**
 * Per-item hit and miss probabilities under per-list storage cost caps,
 * pi_ij = m_j gamma(i,j) E_i(m - e_j, k - sigma_i e_j) / E(m,k). An EMPTY
 * sigma or k selects the unconstrained expression.
 *
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities
 * @param sigma (n) per-item storage costs; empty for none
 * @param k     (h) per-list storage cost caps; empty for none
 * @return (n x (h+1)) matrix; column 0 is the miss probability, column 1+j the
 *         probability of a hit in list j
 */
template <class T>
Matrix<T> cache_prob_erec(const Matrix<T>& gamma, const std::vector<int>& m,
                          const std::vector<int>& sigma, const std::vector<int>& k) {
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    if (h != m.size())
        throw InputError("cache_prob_erec: gamma and m disagree on the number of lists");
    const bool capped = !sigma.empty() && !k.empty();

    const T E = cache_erec(gamma, m, sigma, k);
    if (E == num_traits<T>::from_int(0))
        throw NumericError("cache_prob_erec: the normalizing constant is zero");

    Matrix<T> prob(n, h + 1, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        T hits = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < h; ++j) {
            std::vector<int> mj = m;
            mj[j] -= 1;
            T Ei = num_traits<T>::from_int(0);
            if (!capped) {
                Ei = cache_erec(detail::gamma_without_row(gamma, i), mj);
            } else {
                std::vector<int> kij = k;
                kij[j] -= sigma[i];
                if (kij[j] >= 0) {
                    std::vector<int> si;
                    si.reserve(n - 1);
                    for (std::size_t a = 0; a < n; ++a)
                        if (a != i) si.push_back(sigma[a]);
                    Ei = cache_erec(detail::gamma_without_row(gamma, i), mj, si, kij);
                }
            }
            const T p = num_traits<T>::from_int(static_cast<long>(m[j])) * gamma(i, j) * Ei / E;
            prob(i, j + 1) = p;
            hits += p;
        }
        prob(i, 0) = num_abs(T(num_traits<T>::from_int(1) - hits));
    }
    return prob;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_PROB_EREC_H
