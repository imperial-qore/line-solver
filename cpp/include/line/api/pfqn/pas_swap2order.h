/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PAS_SWAP2ORDER_H
#define LINE_API_PFQN_PAS_SWAP2ORDER_H

/**
 * Global placement-order DAG of a closed two-station pass-and-swap tandem.
 *
 * Templated port of matlab/src/api/pfqn/pas_swap2order.m. With a non-empty
 * swap graph the ordered-state chain is reducible (Comte and Dorsman, 2021,
 * arXiv:2009.12299) and its recurrent class consists of the splits of the
 * orderings that are the linear extensions of a single placement partial
 * order. The port enumerates that recurrent class from the all-in-queue-1
 * state of the minimal single-job-per-class instance, reads the full ordering
 * c = [l1, reverse(l2)] off every reachable state, and sets H(i, j) = 1 iff i
 * precedes j in EVERY such ordering. The placement order is a class-level
 * property, so the result is valid for any population.
 *
 * PARTNER FUNCTION. pas_placement takes the H produced here and returns its
 * transitive closure together with the placeable-next test; the two are meant
 * to be used together and are tested against each other.
 *
 * SERVICE RATES. listRate[m](c) is the total service rate of queue m on the
 * ordered prefix c, and the marginal rate of the p-th customer is the forward
 * difference. The reference guards p == 1 explicitly because some rate handles
 * return a nonzero constant on the empty prefix, and that guard is reproduced:
 * the empty-prefix rate is zero by definition.
 *
 * ARITHMETIC. Only the rate comparison touches the number type, against the
 * reference's fixed 1e-12 threshold, which is a representable rational. No
 * transcendental function is used and the header is instantiated at Rational.
 * The threshold is NOT relaxed to an exact "> 0" test at exact arithmetic:
 * that would be a different algorithm, admitting transitions the reference
 * prunes.
 *
 * INDEXING. Classes are 1-based, as in MATLAB, everywhere a class appears in
 * an ordering or in N0. H is (R x R) and 0-based in the C++ sense, so
 * H(a - 1, b - 1) is MATLAB's H(a, b).
 */

#include <algorithm>
#include <cstddef>
#include <functional>
#include <set>
#include <utility>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Total service rate of a queue on an ordered prefix of classes (1-based). */
template <class T>
using PasRateFun = std::function<T(const std::vector<int>&)>;

namespace detail {

/**
 * One pass-and-swap completion at position p of the ordered queue c
 * (pas_swap_local in the reference).
 *
 * The job at p chases the first later job it may swap with, that one chases
 * the next, and so on; the last job of the chain departs, and every earlier
 * job of the chain moves up into the slot of the next one, position p being
 * removed.
 *
 * @return the new order and the departing class
 */
template <class T>
std::pair<std::vector<int>, int> pas_swap_local(const std::vector<int>& c, std::size_t p,
                                                const Matrix<T>& G) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = c.size();
    std::vector<std::size_t> chain;
    chain.push_back(p);
    int moving = c[p - 1];
    std::size_t cur = p;
    while (true) {
        std::size_t q = 0;
        for (std::size_t j = cur + 1; j <= n; ++j)
            if (G.rows() != 0 && G(static_cast<std::size_t>(moving) - 1,
                                   static_cast<std::size_t>(c[j - 1]) - 1) != zero) {
                q = j;
                break;
            }
        if (q == 0) break;
        chain.push_back(q);
        moving = c[q - 1];
        cur = q;
    }
    const int dep = c[chain.back() - 1];
    std::vector<int> tmp(c);
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) tmp[chain[i + 1] - 1] = c[chain[i] - 1];
    tmp.erase(tmp.begin() + static_cast<std::ptrdiff_t>(chain[0] - 1));
    return std::make_pair(tmp, dep);
}

}  // namespace detail

/**
 * Placement-order DAG of a two-station pass-and-swap tandem.
 *
 * @param swap     one graph, applied to both queues, or one graph per queue;
 *                 G(a, b) nonzero means class a chases class b
 * @param listRate the two ordered service-rate functions, queue 1 then queue 2
 * @param N0       (R) minimal probing population; empty means one job per class
 * @return (R x R) H, with H(i, j) = 1 iff class i+1 must precede class j+1
 */
template <class T>
Matrix<T> pas_swap2order(const std::vector<Matrix<T>>& swap,
                         const std::vector<PasRateFun<T>>& listRate,
                         const std::vector<int>& N0 = std::vector<int>()) {
    const std::size_t M = 2;
    if (swap.empty()) throw InputError("pas_swap2order: no swap graph");
    if (listRate.size() != M) throw InputError("pas_swap2order: two rate functions are required");
    std::vector<Matrix<T>> G(M);
    if (swap.size() == 1) {
        G[0] = swap[0];
        G[1] = swap[0];
    } else if (swap.size() == M) {
        G[0] = swap[0];
        G[1] = swap[1];
    } else {
        throw InputError("pas_swap2order: expected one swap graph or one per queue");
    }

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<int> pop(N0);
    if (pop.empty()) {
        if (G[0].rows() == 0) throw InputError("pas_swap2order: cannot infer the class count");
        pop.assign(G[0].rows(), 1);
    }
    const std::size_t R = pop.size();

    bool anySwap = false;
    for (std::size_t m = 0; m < M; ++m)
        for (std::size_t i = 0; i < G[m].rows(); ++i)
            for (std::size_t j = 0; j < G[m].cols(); ++j)
                if (G[m](i, j) != zero) anySwap = true;
    if (!anySwap) return Matrix<T>(R, R, zero);  // pure OI: no placement constraint

    // depth-first enumeration of the reachable (communicating) class
    typedef std::pair<std::vector<int>, std::vector<int>> State;
    State init;
    for (std::size_t r = 0; r < R; ++r)
        for (int k = 0; k < pop[r]; ++k) init.first.push_back(static_cast<int>(r) + 1);

    std::set<State> seen;
    std::vector<State> frontier, classStates;
    seen.insert(init);
    frontier.push_back(init);
    const T rateTol = num_traits<T>::from_double(1e-12);
    while (!frontier.empty()) {
        const State st = frontier.back();
        frontier.pop_back();
        classStates.push_back(st);
        for (std::size_t m = 0; m < M; ++m) {
            const std::vector<int>& c = m == 0 ? st.first : st.second;
            for (std::size_t p = 1; p <= c.size(); ++p) {
                std::vector<int> prefix(c.begin(), c.begin() + static_cast<std::ptrdiff_t>(p));
                T prevRate = zero;
                if (p > 1) {
                    std::vector<int> shorter(c.begin(),
                                             c.begin() + static_cast<std::ptrdiff_t>(p - 1));
                    prevRate = listRate[m](shorter);
                }
                const T rate = T(listRate[m](prefix) - prevRate);
                if (!(rate > rateTol)) continue;
                const std::pair<std::vector<int>, int> mv = detail::pas_swap_local(c, p, G[m]);
                State stn = st;
                (m == 0 ? stn.first : stn.second) = mv.first;
                (m == 0 ? stn.second : stn.first).push_back(mv.second);
                if (seen.insert(stn).second) frontier.push_back(stn);
            }
        }
    }

    // the orderings exposed by the reachable states
    std::set<std::vector<int>> D;
    for (std::size_t s = 0; s < classStates.size(); ++s) {
        std::vector<int> c(classStates[s].first);
        for (std::size_t k = classStates[s].second.size(); k > 0; --k)
            c.push_back(classStates[s].second[k - 1]);
        D.insert(c);
    }

    Matrix<T> H(R, R, zero);
    for (std::size_t a = 1; a <= R; ++a)
        for (std::size_t b = 1; b <= R; ++b) {
            if (a == b) continue;
            bool both = false, forced = true;
            for (std::set<std::vector<int>>::const_iterator it = D.begin(); it != D.end(); ++it) {
                const std::vector<int>& c = *it;
                std::size_t maxa = 0, minb = 0;
                for (std::size_t k = 0; k < c.size(); ++k) {
                    if (c[k] == static_cast<int>(a)) maxa = k + 1;
                    if (c[k] == static_cast<int>(b) && minb == 0) minb = k + 1;
                }
                if (maxa == 0 || minb == 0) continue;
                both = true;
                if (!(maxa < minb)) {
                    forced = false;
                    break;
                }
            }
            if (both && forced) H(a - 1, b - 1) = one;
        }
    return H;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PAS_SWAP2ORDER_H
