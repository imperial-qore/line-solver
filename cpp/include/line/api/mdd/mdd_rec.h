/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MDD_MDD_REC_H
#define LINE_API_MDD_MDD_REC_H

/**
 * MDD-rec: the normalising constant of a product-form model whose reachable set
 * is held in a decision diagram.
 *
 * S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
 * product-form models of distributed systems with synchronisation", Future
 * Generation Computer Systems 111 (2020) 475-490, Sec. 4.
 *
 * A product-form model has P(s) = (1/G) prod_k g_k(s_k) over its levels, and
 *
 *     G = sum_{s in S} prod_k g_k(s_k).
 *
 * Summing state by state is exponential and numerically unstable. MDD-rec
 * instead walks the diagram that already encodes S, accumulating the
 * unnormalised mass of each node ONCE (Def. 4.4, Algorithm 1):
 *
 *     M(<l.p>) = sum_{v in S_l} g_l(v) * M(<l.p>[v]),   M(TRUE) = 1, M(FALSE) = 0
 *
 * so the cost is O(sum_l |nodes_l| * |S_l|) rather than O(|S|), and G = M(root).
 *
 * FORMALISM-AGNOSTIC. Nothing here knows what a level is: the paper's Appendix B
 * shows that on the lattice sum_k s_k = n of a closed queueing network this
 * collapses to Buzen's convolution, and Sec. 5 that on an S-invariant reachable
 * Petri net it collapses to the Coleman-Henderson-Taylor convolution
 * (`spn::spn_conv`). Unlike either, it needs only that the reachable set be
 * finite and encoded -- no lattice, no S-invariant reachability.
 *
 * THE MASK is how Sec. 5.3 computes measures. Restricting the sum at level l to
 * a subset of its local values gives the unnormalised mass of the corresponding
 * subset of S, so P(m_l = k) and P(e_j >= k) are the same recursion under a
 * different mask rather than three separate algorithms.
 *
 * WHAT IS NOT HERE. The g_l themselves, and the test that the model has a
 * product form at all, are the caller's: the paper declares that out of scope
 * (Sec. 3.2) and refers to the per-formalism conditions instead. Passing g_l
 * that do not describe a product-form model returns a number that is not the
 * normalising constant of anything, and nothing here can detect it.
 */

#include <cstddef>
#include <vector>

#include "line/api/mdd/mdd.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mdd {

/**
 * Per-level admissible local values, the restriction of Sec. 5.3.
 *
 * mask[l][v] false drops local value v of level l from the sum. An EMPTY mask
 * admits everything, which is the plain MDD-rec of Algorithm 1.
 */
typedef std::vector<std::vector<bool>> MddMask;

namespace detail {

/** Unnormalised mass below one node, memoised per (level, node). */
template <class T>
T mdd_rec_node(const MddStruct& mdds, const std::vector<std::vector<T>>& g, const MddMask& mask,
               std::size_t l, int id, std::vector<std::vector<T>>& memo,
               std::vector<std::vector<bool>>& done) {
    if (done[l][id - 1]) return memo[l][id - 1];
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T acc = zero;
    const std::vector<int>& a = mdds.node[l][id - 1];
    for (int v = 0; v < mdds.domain[l]; ++v) {
        if (!mask.empty() && !mask[l][static_cast<std::size_t>(v)]) continue;
        const T gv = g[l][static_cast<std::size_t>(v)];
        if (gv == zero) continue;
        const int ch = a[v];
        T below;
        if (l + 1 == mdds.K) {
            if (ch != TERM_TRUE) continue;
            below = one;
        } else {
            if (ch == TERM_FALSE) continue;
            below = mdd_rec_node(mdds, g, mask, l + 1, ch, memo, done);
        }
        acc += T(gv * below);
    }
    memo[l][id - 1] = acc;
    done[l][id - 1] = true;
    return acc;
}

template <class T>
void mdd_rec_check(const MddStruct& mdds, const std::vector<std::vector<T>>& g,
                   const MddMask& mask, const char* caller) {
    if (g.size() != mdds.K)
        throw InputError(std::string(caller) + ": one g_l per level is required");
    for (std::size_t l = 0; l < mdds.K; ++l) {
        if (g[l].size() != static_cast<std::size_t>(mdds.domain[l]))
            throw InputError(std::string(caller) + ": g_l must have one entry per local state");
        if (!mask.empty() && mask[l].size() != static_cast<std::size_t>(mdds.domain[l]))
            throw InputError(std::string(caller) + ": the mask must have one entry per local state");
    }
    if (!mask.empty() && mask.size() != mdds.K)
        throw InputError(std::string(caller) + ": the mask must have one row per level");
}

}  // namespace detail

/**
 * Unnormalised mass of the masked subset of the reachable set (Algorithm 1).
 *
 * @param mdds the reachable set, in MDD orientation (level 0 is the root)
 * @param g g[l][v] is g_l(v), the per-level factor of the product form
 * @param mask per-level admissible values; empty admits everything and returns G
 */
template <class T>
T mdd_rec_masked(const MddStruct& mdds, const std::vector<std::vector<T>>& g,
                 const MddMask& mask) {
    detail::mdd_rec_check(mdds, g, mask, "mdd_rec");
    const T zero = num_traits<T>::from_int(0);
    if (mdds.root == TERM_FALSE) return zero;
    std::vector<std::vector<T>> memo(mdds.K);
    std::vector<std::vector<bool>> done(mdds.K);
    for (std::size_t l = 0; l < mdds.K; ++l) {
        memo[l].assign(static_cast<std::size_t>(mdds.nnodes[l]), zero);
        done[l].assign(static_cast<std::size_t>(mdds.nnodes[l]), false);
    }
    return detail::mdd_rec_node(mdds, g, mask, 0, mdds.root, memo, done);
}

/** The normalising constant G = sum_{s in S} prod_l g_l(s_l). */
template <class T>
T mdd_rec(const MddStruct& mdds, const std::vector<std::vector<T>>& g) {
    return mdd_rec_masked(mdds, g, MddMask());
}

/**
 * Unnormalised masses of {s in S : s_l = k}, one per local value k of level l.
 *
 * Divided by G these are P(m_l = k) of Sec. 5.3: the mean occupancy of a level
 * is sum_k k * P(m_l = k), and its utilization 1 - P(m_l = 0).
 */
template <class T>
std::vector<T> mdd_rec_marginal(const MddStruct& mdds, const std::vector<std::vector<T>>& g,
                                std::size_t l) {
    if (l >= mdds.K) throw InputError("mdd_rec_marginal: level index is out of range");
    const std::size_t d = static_cast<std::size_t>(mdds.domain[l]);
    std::vector<T> out(d, num_traits<T>::from_int(0));
    for (std::size_t k = 0; k < d; ++k) {
        MddMask mask(mdds.K);
        for (std::size_t j = 0; j < mdds.K; ++j)
            mask[j].assign(static_cast<std::size_t>(mdds.domain[j]), true);
        for (std::size_t v = 0; v < d; ++v) mask[l][v] = (v == k);
        out[k] = mdd_rec_masked(mdds, g, mask);
    }
    return out;
}

}  // namespace mdd
}  // namespace line

#endif  // LINE_API_MDD_MDD_REC_H
