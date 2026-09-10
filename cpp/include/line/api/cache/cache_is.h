/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_IS_H
#define LINE_API_CACHE_IS_H

/**
 * Importance-sampling estimate of the cache normalizing constant.
 *
 * Templated port of matlab/src/api/cache/cache_is.m. A configuration is drawn
 * by choosing mt = sum(m) of the n items uniformly without replacement and
 * splitting them across the lists in the given order, a proposal of constant
 * density 1/(C(n,mt) multinomial(mt;m)); the estimator is the sample mean of
 * the product of the access factors divided by that density, so it is
 * unbiased. It scales to item counts at which the exact enumeration of
 * cache_erec is out of reach, at the price of a Monte Carlo error.
 *
 * MATLAB draws the sample with randperm(n,mt) and then shuffles it again
 * before splitting. This port draws it with a partial Fisher-Yates pass, which
 * produces a uniformly random ordered sequence of mt distinct items directly;
 * the second shuffle of the reference is a permutation of an already exchangeable
 * sequence and changes nothing in distribution. The stream is therefore NOT
 * reproducible against MATLAB run for run: only the estimate is comparable,
 * within its own Monte Carlo error.
 *
 * Arithmetic. static_assert(has_transcendental) -- the estimator is formed in
 * the log domain (log-factorials, a log-sum-exp average) because the raw
 * weights overflow, and it is a Monte Carlo average in any case, so exact
 * arithmetic is meaningless here.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

#include "line/api/cache/cache_erec.h"
#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

namespace detail {

/** gamma with the all-zero rows deleted, i.e. gamma(sum(gamma,2)>0,:). */
template <class T>
Matrix<T> gamma_drop_zero_rows(const Matrix<T>& gamma) {
    const T zero = num_traits<T>::from_int(0);
    std::size_t keep = 0;
    for (std::size_t i = 0; i < gamma.rows(); ++i) {
        T s = zero;
        for (std::size_t j = 0; j < gamma.cols(); ++j) s += gamma(i, j);
        if (s > zero) ++keep;
    }
    Matrix<T> g(keep, gamma.cols());
    std::size_t r = 0;
    for (std::size_t i = 0; i < gamma.rows(); ++i) {
        T s = zero;
        for (std::size_t j = 0; j < gamma.cols(); ++j) s += gamma(i, j);
        if (s <= zero) continue;
        for (std::size_t j = 0; j < gamma.cols(); ++j) g(r, j) = gamma(i, j);
        ++r;
    }
    return g;
}

/**
 * mt distinct items out of n, uniformly, as an ordered sequence: a partial
 * Fisher-Yates pass over a lazily materialized identity permutation.
 */
inline void sample_without_replacement(std::size_t n, std::size_t mt, std::mt19937_64& rng,
                                       std::vector<std::size_t>& perm,
                                       std::vector<std::size_t>& out) {
    perm.resize(n);
    for (std::size_t i = 0; i < n; ++i) perm[i] = i;
    out.resize(mt);
    for (std::size_t k = 0; k < mt; ++k) {
        std::uniform_int_distribution<std::size_t> d(k, n - 1);
        const std::size_t j = d(rng);
        const std::size_t t = perm[k];
        perm[k] = perm[j];
        perm[j] = t;
        out[k] = perm[k];
    }
}

}  // namespace detail

template <class T>
struct CacheIsResult {
    T E;   ///< normalizing constant estimate
    T lE;  ///< its logarithm
};

/**
 * @param gamma_in (n x h) access factors
 * @param m       (h) list capacities
 * @param samples number of Monte Carlo samples (MATLAB default 1e5)
 * @param seed    seed of the sampling stream
 * @param sigma_in (n) per-item storage cost; empty for uncapped lists
 * @param k (h) per-list cost cap; empty for uncapped lists
 */
template <class T>
CacheIsResult<T> cache_is(const Matrix<T>& gamma_in, const std::vector<int>& m,
                          std::size_t samples, std::uint64_t seed,
                          const std::vector<int>& sigma_in, const std::vector<int>& k) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_is requires transcendental arithmetic: the estimator is a Monte Carlo "
                  "average formed in the log domain");
    if (gamma_in.cols() != m.size())
        throw InputError("cache_is: gamma and m disagree on the number of lists");
    if (samples == 0) throw InputError("cache_is: at least one sample is required");

    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    const bool capped = !sigma_in.empty() && !k.empty();
    const Matrix<T> gamma = detail::gamma_drop_zero_rows(gamma_in);
    std::vector<int> sigma;
    if (capped) {
        // gamma_drop_zero_rows keeps only the rows with a positive row sum
        for (std::size_t i = 0; i < gamma_in.rows(); ++i) {
            T rs = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < gamma_in.cols(); ++j) rs += gamma_in(i, j);
            if (rs > num_traits<T>::from_int(0)) sigma.push_back(sigma_in[i]);
        }
    }
    const std::size_t n = gamma.rows();
    const std::size_t h = m.size();
    long mt = 0;
    for (int v : m) {
        if (v < 0) throw InputError("cache_is: negative list capacity");
        mt += v;
    }

    CacheIsResult<T> res;
    if (n == 0 || mt == 0) {
        res.E = one;
        res.lE = zero;
        return res;
    }
    if (static_cast<long>(n) < mt) {
        // Fewer items than cache slots: no valid configuration exists, exactly
        // as the reference reports (it warns and returns zero).
        res.E = zero;
        res.lE = -T(std::numeric_limits<T>::infinity());
        return res;
    }
    if (static_cast<long>(n) == mt) {
        // Every item must be cached: a single configuration, taken exactly.
        res.E = cache_erec(gamma, m, sigma, k);
        res.lE = log(res.E);
        return res;
    }

    Matrix<T> lgam(n, h, zero);
    const T floorv = num_traits<T>::from_double(1e-300);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < h; ++j) lgam(i, j) = log(T(gamma(i, j) + floorv));

    T logMFact = zero;
    for (std::size_t j = 0; j < h; ++j)
        logMFact += pfqn::detail::num_logfact_int<T>(static_cast<long>(m[j]));
    const T logComb = pfqn::detail::num_logfact_int<T>(static_cast<long>(n)) -
                      pfqn::detail::num_logfact_int<T>(mt) -
                      pfqn::detail::num_logfact_int<T>(static_cast<long>(n) - mt);
    const T logMultinom = pfqn::detail::num_logfact_int<T>(mt) - logMFact;
    const T logProposal = -logComb - logMultinom;

    std::mt19937_64 rng(seed);
    std::vector<std::size_t> perm, sel;
    std::vector<T> lZ(samples, zero);
    for (std::size_t s = 0; s < samples; ++s) {
        detail::sample_without_replacement(n, static_cast<std::size_t>(mt), rng, perm, sel);
        T logState = logMFact;
        std::size_t idx = 0;
        bool feasible = true;
        for (std::size_t j = 0; j < h && feasible; ++j) {
            if (capped) {
                long listCost = 0;
                for (int c = 0; c < m[j]; ++c) listCost += sigma[sel[idx + c]];
                if (listCost > k[j]) {
                    feasible = false;
                    break;
                }
            }
            for (int c = 0; c < m[j]; ++c) logState += lgam(sel[idx++], j);
        }
        if (!feasible) {
            lZ[s] = -T(std::numeric_limits<T>::infinity());  // I{S_v in O} = 0
            continue;
        }
        lZ[s] = logState - logProposal;
    }
    res.lE = pfqn::detail::logsumexp(lZ) - log(num_traits<T>::from_int(static_cast<long>(samples)));
    using std::exp;
    res.E = exp(res.lE);
    return res;
}

/** cache_is without storage cost caps. */
template <class T>
CacheIsResult<T> cache_is(const Matrix<T>& gamma_in, const std::vector<int>& m,
                          std::size_t samples, std::uint64_t seed) {
    return cache_is(gamma_in, m, samples, seed, std::vector<int>(), std::vector<int>());
}

/** cache_is with the MATLAB default of 1e5 samples. */
template <class T>
CacheIsResult<T> cache_is(const Matrix<T>& gamma, const std::vector<int>& m) {
    return cache_is(gamma, m, static_cast<std::size_t>(100000), static_cast<std::uint64_t>(0));
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_IS_H
