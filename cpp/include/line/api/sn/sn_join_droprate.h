/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_JOIN_DROPRATE_H
#define LINE_API_SN_SN_JOIN_DROPRATE_H

#include <cmath>
#include <cstddef>
#include <map>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace sn {

/**
 * The number of siblings the Join node `joinNode` (1-based) fires on, out of
 * `nbranches` forked. Port of `matlab/src/api/fj/sn_join_quorum.m`.
 *
 * A standard join, an absent declaration, a non-positive quorum and a quorum
 * that is not smaller than the sibling count all return `nbranches`, the
 * ordinary AND-join.
 *
 * NOTE the index space: `joindecl` is keyed by the BASE node index and carries
 * ONE quorum per join rather than one per class, so a multiclass model with a
 * different quorum per class cannot be expressed here. The reference and the
 * JSON interchange have the same limitation.
 */
template <class T>
std::size_t sn_join_quorum(const qn::NetworkStruct<T>& sn, std::size_t joinNode,
                           std::size_t nbranches) {
    if (joinNode == 0) return nbranches;
    typename std::map<std::size_t, typename qn::NetworkStruct<T>::JoinDecl>::const_iterator it =
        sn.joindecl.find(joinNode);
    if (it == sn.joindecl.end()) return nbranches;
    if (it->second.strategy == lang::JoinStrategy::STD) return nbranches;
    const double q = it->second.quorum;
    if (q > 0.0) {
        const std::size_t k = static_cast<std::size_t>(q + 0.5);
        if (k > 0 && k < nbranches) return k;
    }
    return nbranches;
}

/**
 * Rate at which sibling tasks are discarded at each Join, as an
 * (nstations x nclasses) matrix that is zero away from the Join rows.
 * Port of `matlab/src/api/fj/sn_join_droprate.m`.
 *
 * A Join is the one station where the loss identity ArvR - Tput does NOT hold,
 * because the two rates are in different units: AN counts the SIBLINGS offered
 * to the join (N per parent job) while TN counts the PARENT jobs released by it
 * (one per synchronisation). Reading ArvR - Tput there reports (N-1)/N of the
 * offered traffic as lost at every join, standard joins included, when a
 * standard join loses nothing at all.
 *
 * The siblings a join actually consumes are K per synchronisation, where K is
 * the quorum (K = N on a standard join), so
 *
 *     DropRateJoin = max(0, AN - K*TN)
 *
 * which is 0 for a standard join and (N-K)*TN for a quorum, the rate at which
 * the stragglers of an already-fired parent are discarded on arrival.
 *
 * This is the DERIVED value, exact given TN and AN. A solver that MEASURES the
 * discards on its own sample path reports its own: the LDES engine fills
 * `LdesResult::DropRateJoin` from the sibling it actually threw away, and a
 * measured zero is a better answer than this one, which carries the
 * finite-sample gap between the two rates.
 */
template <class T>
Matrix<T> sn_join_droprate(const qn::NetworkStruct<T>& sn, const Matrix<T>& TN,
                           const Matrix<T>& AN) {
    const std::size_t M = sn.stations.size(), K = sn.classes.size();
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> out(M, K, zero);
    if (sn.fj.empty() || TN.rows() == 0 || AN.rows() == 0) return out;
    for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
        if (sn.nodes[ind - 1].nodetype != qn::NodeType::Join) continue;
        const std::size_t ist = sn.nodes[ind - 1].station;
        if (ist == 0 || ist > M) continue;
        for (std::size_t r = 0; r < K; ++r) {
            if (ist - 1 >= AN.rows() || r >= AN.cols()) continue;
            const T a = AN(ist - 1, r);
            const T t = (ist - 1 < TN.rows() && r < TN.cols()) ? TN(ist - 1, r) : zero;
            if (!std::isfinite(num_traits<T>::to_double(a))) continue;
            if (!std::isfinite(num_traits<T>::to_double(t))) continue;
            if (!(a > zero)) continue;
            // PER CLASS: a variable forking level makes the sibling count differ
            // between classes, so it cannot be hoisted out of this loop.
            const std::size_t nsib = sn.join_siblings(ind, r + 1);
            if (nsib == 0) continue;
            const std::size_t kreq = sn_join_quorum(sn, ind, nsib);
            const T d = T(a - num_traits<T>::from_int(static_cast<long>(kreq)) * t);
            out(ist - 1, r) = (d > zero) ? d : zero;
        }
    }
    return out;
}

}  // namespace sn
}  // namespace line

#endif  // LINE_API_SN_SN_JOIN_DROPRATE_H
