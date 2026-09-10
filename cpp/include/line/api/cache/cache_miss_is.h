/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_MISS_IS_H
#define LINE_API_CACHE_MISS_IS_H

/**
 * Cache miss rates from the importance-sampling hit probabilities.
 *
 * Templated port of matlab/src/api/cache/cache_miss_is.m: the per-item miss
 * probabilities come from cache_prob_is, and the miss rates are the request
 * rates weighted by them. The normalizing constant is estimated alongside by
 * cache_is and returned, as the reference does.
 *
 * With no request rates the reference returns the mean miss probability as the
 * global rate and leaves the per-user and per-item vectors empty; that
 * contract is kept.
 *
 * Arithmetic. static_assert(has_transcendental) -- it is the importance
 * sampler of cache_prob_is with a linear map on top.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_is.h"
#include "line/api/cache/cache_prob_is.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

template <class T>
struct CacheMissIsResult {
    T M;                 ///< global miss rate
    std::vector<T> MU;   ///< (u) per-user miss rate; empty when no lambda given
    std::vector<T> MI;   ///< (n) per-item miss rate; empty when no lambda given
    std::vector<T> pi0;  ///< (n) per-item miss probability
    T lE;                ///< log of the normalizing constant estimate
};

/**
 * @param gamma   (n x h) access factors
 * @param m       (h) list capacities
 * @param lambda  (u x n) per-user per-item request rates, MATLAB's
 *                lambda(:,:,1); empty for the mean miss probability alone
 * @param samples number of Monte Carlo samples
 * @param seed    seed of the sampling stream
 * @param sigma (n) per-item storage cost; empty for uncapped lists
 * @param cap (h) per-list cost cap; empty for uncapped lists
 */
template <class T>
CacheMissIsResult<T> cache_miss_is(const Matrix<T>& gamma, const std::vector<int>& m,
                                   const Matrix<T>& lambda, std::size_t samples,
                                   std::uint64_t seed, const std::vector<int>& sigma,
                                   const std::vector<int>& cap) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_miss_is requires transcendental arithmetic: it is the cache_prob_is "
                  "importance sampler with a linear map on top");
    if (gamma.cols() != m.size())
        throw InputError("cache_miss_is: gamma and m disagree on the number of lists");

    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = gamma.rows();

    CacheMissIsResult<T> res;
    res.M = zero;
    res.lE = cache_is(gamma, m, samples, seed, sigma, cap).lE;
    const Matrix<T> pij = cache_prob_is(gamma, m, samples, seed, sigma, cap);
    res.pi0.assign(n, zero);
    for (std::size_t k = 0; k < n; ++k) res.pi0[k] = pij(k, 0);

    if (lambda.empty()) {
        for (std::size_t k = 0; k < n; ++k) res.M += res.pi0[k];
        if (n > 0) res.M /= num_traits<T>::from_int(static_cast<long>(n));
        return res;
    }
    if (lambda.cols() != n)
        throw InputError("cache_miss_is: lambda and gamma disagree on the number of items");

    const std::size_t u = lambda.rows();
    res.MU.assign(u, zero);
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t k = 0; k < n; ++k) res.MU[v] += lambda(v, k) * res.pi0[k];

    res.MI.assign(n, zero);
    for (std::size_t k = 0; k < n; ++k) {
        T s = zero;
        for (std::size_t v = 0; v < u; ++v) s += lambda(v, k);
        res.MI[k] = s * res.pi0[k];
        res.M += res.MI[k];
    }
    return res;
}

/** cache_miss_is without storage cost caps. */
template <class T>
CacheMissIsResult<T> cache_miss_is(const Matrix<T>& gamma, const std::vector<int>& m,
                                   const Matrix<T>& lambda, std::size_t samples,
                                   std::uint64_t seed) {
    return cache_miss_is(gamma, m, lambda, samples, seed, std::vector<int>(), std::vector<int>());
}

/** cache_miss_is with the MATLAB default of 1e5 samples. */
template <class T>
CacheMissIsResult<T> cache_miss_is(const Matrix<T>& gamma, const std::vector<int>& m,
                                   const Matrix<T>& lambda) {
    return cache_miss_is(gamma, m, lambda, static_cast<std::size_t>(100000),
                         static_cast<std::uint64_t>(0));
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_MISS_IS_H
