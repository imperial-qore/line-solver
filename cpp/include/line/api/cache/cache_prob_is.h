/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_PROB_IS_H
#define LINE_API_CACHE_PROB_IS_H

/**
 * Importance-sampling estimate of the cache hit-probability distribution.
 *
 * Templated port of matlab/src/api/cache/cache_prob_is.m. The same proposal as
 * cache_is is used; prob(i,1+j) is the ratio of the accumulated importance
 * weight of the configurations that place item i on list j to the total
 * weight, and prob(i,0) is the residual miss probability. The self-normalized
 * form makes the constant proposal density cancel, so the estimator depends on
 * the weights only through their ratios.
 *
 * The reference scales the weights by exp(-50) before accumulating them, to
 * keep the sum inside the double range; the constant cancels in the ratio and
 * is kept here for entry-by-entry agreement with MATLAB.
 *
 * Arithmetic. static_assert(has_transcendental) -- Monte Carlo weights formed
 * through logs and an exponential, as in cache_is.
 */

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "line/api/cache/cache_is.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/**
 * @param gamma   (n x h) access factors
 * @param m       (h) list capacities
 * @param samples number of Monte Carlo samples (MATLAB default 1e5)
 * @param seed    seed of the sampling stream
 * @param sigma (n) per-item storage cost; empty for uncapped lists
 * @param k (h) per-list cost cap; empty for uncapped lists
 * @return (n x h+1); column 0 is the miss probability, column 1+j the
 *         probability that the item sits on list j
 */
template <class T>
Matrix<T> cache_prob_is(const Matrix<T>& gamma, const std::vector<int>& m, std::size_t samples,
                        std::uint64_t seed, const std::vector<int>& sigma,
                        const std::vector<int>& k) {
    const bool capped = !sigma.empty() && !k.empty();
    static_assert(num_traits<T>::has_transcendental,
                  "cache_prob_is requires transcendental arithmetic: the weights are Monte Carlo "
                  "quantities formed through logs and an exponential");
    if (gamma.cols() != m.size())
        throw InputError("cache_prob_is: gamma and m disagree on the number of lists");
    if (samples == 0) throw InputError("cache_prob_is: at least one sample is required");

    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t n = gamma.rows();
    const std::size_t h = m.size();
    long mt = 0;
    for (int v : m) {
        if (v < 0) throw InputError("cache_prob_is: negative list capacity");
        mt += v;
    }

    Matrix<T> prob(n, h + 1, zero);
    if (n == 0) return prob;
    if (mt == 0 || static_cast<long>(n) < mt) {
        // No slot, or fewer items than slots: every item misses, as the
        // reference reports.
        for (std::size_t i = 0; i < n; ++i) prob(i, 0) = one;
        return prob;
    }
    if (static_cast<long>(n) == mt) return cache_prob_erec(gamma, m, sigma, k);

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
    const T shift = num_traits<T>::from_int(50);  // the reference's overflow guard

    Matrix<T> acc(n, h, zero);
    T total = zero;
    std::mt19937_64 rng(seed);
    std::vector<std::size_t> perm, sel;
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
        if (!feasible) continue;  // I{S_v in O} = 0
        const T wgt = exp(T(logState - logProposal - shift));
        total += wgt;
        idx = 0;
        for (std::size_t j = 0; j < h; ++j)
            for (int c = 0; c < m[j]; ++c) acc(sel[idx++], j) += wgt;
    }

    if (total <= zero) {
        for (std::size_t i = 0; i < n; ++i) prob(i, 0) = one;
        return prob;
    }
    for (std::size_t i = 0; i < n; ++i) {
        T hit = zero;
        for (std::size_t j = 0; j < h; ++j) {
            prob(i, 1 + j) = acc(i, j) / total;
            hit += prob(i, 1 + j);
        }
        const T miss = one - hit;
        prob(i, 0) = miss > zero ? miss : zero;
    }
    return prob;
}

/** cache_prob_is without storage cost caps. */
template <class T>
Matrix<T> cache_prob_is(const Matrix<T>& gamma, const std::vector<int>& m, std::size_t samples,
                        std::uint64_t seed) {
    return cache_prob_is(gamma, m, samples, seed, std::vector<int>(), std::vector<int>());
}

/** cache_prob_is with the MATLAB default of 1e5 samples. */
template <class T>
Matrix<T> cache_prob_is(const Matrix<T>& gamma, const std::vector<int>& m) {
    return cache_prob_is(gamma, m, static_cast<std::size_t>(100000), static_cast<std::uint64_t>(0));
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_PROB_IS_H
