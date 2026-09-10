/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SPN_SPN_CONV_H
#define LINE_API_SPN_SPN_CONV_H

/**
 * Convolution algorithm for the normalising constant of an S-invariant
 * reachable product-form stochastic Petri net.
 *
 * J. Coleman, W. Henderson, P. Taylor, "Product form equilibrium distributions
 * and a convolution algorithm for stochastic Petri nets", Performance
 * Evaluation 26(3), 1996, 159-180; presented as the point of comparison for
 * MDD-rec in S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020), Sec. 5.1.
 *
 * With S the minimal-support S-invariant matrix and V = S m0 the load vector,
 * the reachability set of an S-INVARIANT REACHABLE net is exactly
 * {m >= 0 : S m = V}, and conditioning on the marking of one place partitions it
 * (Lemma 5.1). Writing G_j(W) for the mass of the markings supported on the
 * first j places with S m = W,
 *
 *     G_0(W) = [W == 0],     G_j(W) = sum_i g_j(i) G_{j-1}(W - i S_j),
 *
 * and G = G_n(V). On a net whose only invariant is "the tokens are conserved"
 * this is Buzen's convolution for a closed queueing network, one place per
 * station.
 *
 * NO ILP IS SOLVED. The paper obtains the marking set M_p(P',W) from the
 * feasibility of an integer program (Prop. 5.2) so that the sum skips the terms
 * that contribute nothing. Here the sum simply runs over i whose residual
 * W - i S_j stays non-negative and the recursion returns zero on an infeasible
 * residual, which gives the same value: the ILP is an optimisation of the
 * enumeration, not part of the definition. Memoising on (j, W) keeps the walk
 * over the reachable residuals rather than over all of them.
 *
 * S-INVARIANT REACHABILITY IS NOT CHECKED, and cannot be cheaply: no algorithm
 * is known that decides it without generating the reachability set (FGCS, Sec.
 * 5.1). On a net that fails it, {m : S m = V} is strictly larger than the
 * reachable set and this returns a normalising constant over unreachable
 * markings too, which is why `mdd::mdd_rec` -- which walks the reachable set
 * itself -- is the general algorithm and this one the special case. Compare the
 * two on a new net before trusting this one on it.
 */

#include <cstddef>
#include <map>
#include <vector>

#include "line/api/spn/spn_sinvariants.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace spn {

namespace detail {

template <class T>
T spn_conv_rec(std::size_t j, const std::vector<long long>& W,
               const std::vector<std::vector<long long>>& S,
               const std::vector<std::vector<T>>& g,
               std::map<std::pair<std::size_t, std::vector<long long>>, T>& memo) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (j == 0) {
        for (std::size_t r = 0; r < W.size(); ++r)
            if (W[r] != 0) return zero;
        return one;
    }
    const std::pair<std::size_t, std::vector<long long>> key(j, W);
    const typename std::map<std::pair<std::size_t, std::vector<long long>>, T>::const_iterator it =
        memo.find(key);
    if (it != memo.end()) return it->second;

    const std::size_t p = j - 1;
    T acc = zero;
    for (std::size_t i = 0; i < g[p].size(); ++i) {
        std::vector<long long> rem = W;
        bool feasible = true;
        for (std::size_t r = 0; r < W.size() && feasible; ++r) {
            rem[r] -= static_cast<long long>(i) * S[r][p];
            if (rem[r] < 0) feasible = false;
        }
        if (!feasible) break;  // S is non-negative, so larger i only gets worse
        if (g[p][i] == zero) continue;
        acc += T(g[p][i] * spn_conv_rec(j - 1, rem, S, g, memo));
    }
    memo[key] = acc;
    return acc;
}

}  // namespace detail

/**
 * The normalising constant by convolution over the invariant load vector.
 *
 * @param S S[i][p], the minimal-support S-invariants, one row per invariant
 * @param V the load vector S m0, one entry per invariant
 * @param g g[p][i] is g_p(i), the product-form factor of i tokens in place p;
 *        its length bounds the marking of place p
 */
template <class T>
T spn_conv(const std::vector<std::vector<long long>>& S, const std::vector<long long>& V,
           const std::vector<std::vector<T>>& g) {
    if (S.empty()) throw InputError("spn_conv: the net has no S-invariant to convolve over");
    if (S.size() != V.size())
        throw InputError("spn_conv: one load-vector entry per invariant is required");
    const std::size_t n = g.size();
    for (std::size_t r = 0; r < S.size(); ++r) {
        if (S[r].size() != n)
            throw InputError("spn_conv: the invariant matrix and g must agree on the place count");
        for (std::size_t p = 0; p < n; ++p)
            if (S[r][p] < 0)
                throw InputError("spn_conv: an S-invariant has a negative weight, so the residual "
                                 "recursion has no monotone bound on the marking");
    }
    std::map<std::pair<std::size_t, std::vector<long long>>, T> memo;
    return detail::spn_conv_rec(n, V, S, g, memo);
}

/** Convolve straight off the invariant basis `spn_sinvariants` returned. */
template <class T>
T spn_conv(const SpnInvariants& inv, const std::vector<std::vector<T>>& g) {
    return spn_conv(inv.S, inv.V, g);
}

}  // namespace spn
}  // namespace line

#endif  // LINE_API_SPN_SPN_CONV_H
