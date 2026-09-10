/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_TRAFFIC_SPLIT_RR_H
#define LINE_API_NPFQN_TRAFFIC_SPLIT_RR_H

/**
 * Deterministic (round-robin) split degrees of every station-class departure
 * stream.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_traffic_split_rr.m,
 * cross-checked against
 * jar/src/main/java/jline/api/npfqn/Npfqn_traffic_split_rr.java and
 * python/line_solver/api/npfqn/split_rr.py (identical).
 *
 * kRR(i,r) = k > 1 says the class-r departures of station i are dispatched
 * one-in-k by a round-robin node, so a downstream flow carrying a fraction p of
 * them is the k-fold convolution thinned at q = k p and has SCV 1 + p (d2 - k),
 * against the Markovian 1 + p (d2 - 1). kRR(i,r) = 1 is an ordinary
 * probabilistic split, which is what every entry stays at on a model with no
 * round-robin dispatcher; the QNA/MNA equations then reduce exactly to the
 * Bernoulli-thinning form they had before this factor existed.
 *
 * Only two topologies admit the deterministic rule: the station dispatches
 * round-robin itself, or it feeds WITH PROBABILITY ONE a non-station node that
 * does and whose pointer no other flow advances. A router shared by two
 * upstream streams interleaves them, so neither sees a clean one-in-k split and
 * both fall back to k = 1.
 *
 * THE DEGREE IS READ OFF `rtnodes`, NOT off `connmatrix` and not off an
 * `outlinks` list. The reference prefers `sn.nodeparam{ind}{r}.outlinks`, the
 * per-class destination list of the dispatcher, and falls back to the row of
 * `sn.connmatrix`. This port's NetworkStruct carries neither field: `rtnodes`
 * is the same graph once refresh_routing has expanded the strategies (the same
 * substitution downstream_stations documents in network_struct.h). Counting the
 * DISTINCT destination NODES of row (ind,r) reproduces the per-class `outlinks`
 * semantics, which is the branch the reference takes whenever a dispatcher was
 * actually built, and differs from the connmatrix fallback only on a link the
 * model declares and then routes no mass over.
 *
 * ARITHMETIC: field, and integral in fact. The result is a matrix of counts;
 * only comparisons against zero and against 1 - FineTol are taken, so this is
 * exact at T = Rational.
 */

#include <cstddef>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

namespace detail {

/** The routing strategy of (1-based node `ind`, 1-based class `r`). */
template <class T>
qn::RoutingStrategy routing_of(const qn::NetworkStruct<T>& sn, std::size_t ind, std::size_t r) {
    const std::vector<qn::RoutingStrategy>& rs = sn.nodes[ind - 1].routing;
    return rs.size() >= r ? rs[r - 1] : qn::RoutingStrategy::PROB;
}

/** Number of destinations the round-robin pointer of (`ind`, `r`) cycles through. */
template <class T>
std::size_t rr_degree(const qn::NetworkStruct<T>& sn, std::size_t ind, std::size_t r) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = sn.nclasses, I = sn.nodes.size();
    std::size_t k = 0;
    for (std::size_t j = 1; j <= I; ++j) {
        bool linked = false;
        for (std::size_t s = 1; s <= K && !linked; ++s)
            if (sn.rtnodes((ind - 1) * K + (r - 1), (j - 1) * K + (s - 1)) > zero) linked = true;
        if (linked) ++k;
    }
    return k < 1 ? 1 : k;
}

}  // namespace detail

/**
 * Port of `npfqn_traffic_split_rr.m`.
 *
 * @param sn the model, after refresh(): `rtnodes` carries the expanded routing
 * @return an (nstations x nclasses) matrix of split degrees, 1 where the split
 *         is Markovian
 */
template <class T>
Matrix<T> npfqn_traffic_split_rr(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = sn.nstations, K = sn.nclasses, I = sn.nodes.size();
    Matrix<T> kRR(M, K, one);

    bool any_rr = false;
    for (const qn::NodeDef& nd : sn.nodes)
        for (qn::RoutingStrategy rs : nd.routing)
            if (rs == qn::RoutingStrategy::RROBIN) any_rr = true;
    if (!any_rr) return kRR;
    if (sn.rtnodes.rows() < I * K) return kRR;

    const T sure = T(one - num_traits<T>::from_double(lang::GlobalConstants::FineTol));
    for (std::size_t ist = 1; ist <= M; ++ist) {
        const std::size_t ind = sn.node_of_station(ist);
        if (ind == 0) continue;
        for (std::size_t r = 1; r <= K; ++r) {
            if (detail::routing_of(sn, ind, r) == qn::RoutingStrategy::RROBIN) {
                kRR(ist - 1, r - 1) =
                    num_traits<T>::from_int(static_cast<long>(detail::rr_degree(sn, ind, r)));
                continue;
            }
            // a sure transition into a node that dispatches round-robin
            const std::size_t row = (ind - 1) * K + (r - 1);
            std::size_t dest = 0, ndest = 0;
            for (std::size_t col = 0; col < sn.rtnodes.cols(); ++col)
                if (sn.rtnodes(row, col) > zero) {
                    dest = col;
                    ++ndest;
                }
            if (ndest != 1) continue;
            const std::size_t jnd = dest / K + 1, s = dest % K + 1;
            if (sn.nodes[jnd - 1].station != 0) continue;
            if (detail::routing_of(sn, jnd, s) != qn::RoutingStrategy::RROBIN) continue;
            if (sn.rtnodes(row, dest) < sure) continue;
            // the round-robin pointer must be advanced by this stream alone
            std::size_t nfeed = 0;
            for (std::size_t a = 0; a < sn.rtnodes.rows(); ++a)
                if (sn.rtnodes(a, dest) > zero) ++nfeed;
            if (nfeed != 1) continue;
            kRR(ist - 1, r - 1) =
                num_traits<T>::from_int(static_cast<long>(detail::rr_degree(sn, jnd, s)));
        }
    }
    return kRR;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_TRAFFIC_SPLIT_RR_H
