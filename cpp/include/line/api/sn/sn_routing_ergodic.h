/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_ROUTING_ERGODIC_H
#define LINE_API_SN_SN_ROUTING_ERGODIC_H

/**
 * Reducibility of a network's ROUTING, and the repair that makes it ergodic.
 *
 * Port of the MATLAB `@MNetwork` methods `isRoutingErgodic`,
 * `getReducibilityInfo`, `getAbsorbingStations` and `makeErgodic`, and twin of
 * the JAR `jline.lang.RoutingErgodicity` and the native Python
 * `Network.get_reducibility_info` / `make_ergodic`.
 *
 * This inspects the ROUTING STRUCTURE ONLY, not the state space. Ergodic
 * routing is NECESSARY but not sufficient for the chain to be ergodic:
 * state-dependent routing, finite buffers and blocking can still make the CTMC
 * reducible. The converse is the useful direction -- reducible routing is a
 * defect the modeller can see and fix before any solver runs, and it is what
 * otherwise surfaces as "the generator has no recurrent state" inside
 * ctmc_solve.
 *
 * The adjacency is built over the NODE-level blocks `sn.P[(r,s)]`, folding
 * every class pair into one graph: a job that leaves node i as another class
 * has still left i. A station is ABSORBING when it has routing defined and no
 * edge to any other node; the Sink is removed from that list, being
 * legitimately absorbing in an open network.
 *
 * ARITHMETIC: comparisons against zero only, so this is exact in every
 * instantiation and needs no transcendental gate.
 */

#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/api/mc/stronglyconncomp.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace sn {

/** The reducibility structure of a routing matrix, the MATLAB info struct. */
struct RoutingErgodicityInfo {
    bool isRoutingErgodic = true;             ///< the routing is irreducible
    bool isReducible = false;                 ///< the routing is reducible
    std::vector<std::string> absorbingStations;  ///< names, the Sink excluded
    std::vector<std::string> transientStations;  ///< names in a transient component
    std::size_t numSCCs = 1;                  ///< strongly connected components
    std::vector<std::string> suggestedFixes;  ///< one repair per absorbing station
};

namespace detail {

/** Node-level adjacency of the routing, every class pair folded into one graph. */
template <class T>
Matrix<T> routing_adjacency(const qn::NetworkStruct<T>& sn) {
    const std::size_t I = sn.nodes.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> adj(I, I, zero);
    for (typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T>>::const_iterator it =
             sn.P.begin();
         it != sn.P.end(); ++it) {
        const Matrix<T>& blk = it->second;
        const std::size_t nr = std::min(blk.rows(), I);
        const std::size_t nc = std::min(blk.cols(), I);
        for (std::size_t i = 0; i < nr; ++i)
            for (std::size_t j = 0; j < nc; ++j)
                if (blk(i, j) > zero) adj(i, j) = one;
    }
    return adj;
}

/** True when node IDX has any outgoing routing probability at all. */
template <class T>
bool has_any_routing(const qn::NetworkStruct<T>& sn, std::size_t idx) {
    const T zero = num_traits<T>::from_int(0);
    for (typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T>>::const_iterator it =
             sn.P.begin();
         it != sn.P.end(); ++it) {
        const Matrix<T>& blk = it->second;
        if (idx >= blk.rows()) continue;
        for (std::size_t j = 0; j < blk.cols(); ++j)
            if (blk(idx, j) > zero) return true;
    }
    return false;
}

/** First Delay, else the first station that is not itself absorbing. */
template <class T>
std::string default_target(const qn::NetworkStruct<T>& sn,
                           const std::vector<std::string>& absorbing) {
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].nodetype == lang::NodeType::Delay) return sn.nodes[i].name;
    std::string first;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        if (sn.nodes[i].station == 0) continue;
        if (first.empty()) first = sn.nodes[i].name;
        if (std::find(absorbing.begin(), absorbing.end(), sn.nodes[i].name) == absorbing.end())
            return sn.nodes[i].name;
    }
    return first;
}

}  // namespace detail

/** Ergodicity of the routing, with the structure behind the verdict. */
template <class T>
RoutingErgodicityInfo sn_is_routing_ergodic(const qn::NetworkStruct<T>& sn) {
    RoutingErgodicityInfo info;
    if (sn.P.empty()) return info;  // no routing defined yet, nothing to refute

    const Matrix<T> adj = detail::routing_adjacency(sn);
    const mc::SccResult scc = mc::stronglyconncomp(adj);
    info.numSCCs = scc.numSCC();

    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        if (sn.nodes[i].station == 0) continue;
        bool outgoing = false;
        for (std::size_t j = 0; j < sn.nodes.size(); ++j)
            if (j != i && adj(i, j) > zero) {
                outgoing = true;
                break;
            }
        if (!outgoing && detail::has_any_routing(sn, i))
            info.absorbingStations.push_back(sn.nodes[i].name);
    }
    // The Sink is legitimately absorbing in an open network
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].nodetype == lang::NodeType::Sink) {
            info.absorbingStations.erase(std::remove(info.absorbingStations.begin(),
                                                     info.absorbingStations.end(),
                                                     sn.nodes[i].name),
                                         info.absorbingStations.end());
        }

    std::size_t recurrentCount = 0;
    for (std::size_t c = 1; c <= info.numSCCs; ++c) {
        if (scc.recurrent[c - 1]) {
            ++recurrentCount;
            continue;
        }
        for (std::size_t v = 0; v < scc.scc.size(); ++v) {
            if (scc.scc[v] != c || sn.nodes[v].station == 0) continue;
            const std::string& nm = sn.nodes[v].name;
            if (std::find(info.absorbingStations.begin(), info.absorbingStations.end(), nm) ==
                info.absorbingStations.end())
                info.transientStations.push_back(nm);
        }
    }

    info.isReducible = !info.absorbingStations.empty() || recurrentCount > 1;
    info.isRoutingErgodic = !info.isReducible;
    return info;
}

/** The reducibility structure plus a suggested repair per absorbing station. */
template <class T>
RoutingErgodicityInfo sn_reducibility_info(const qn::NetworkStruct<T>& sn) {
    RoutingErgodicityInfo info = sn_is_routing_ergodic(sn);
    if (info.isRoutingErgodic) return info;
    const std::string target = detail::default_target(sn, info.absorbingStations);
    for (std::size_t i = 0; i < info.absorbingStations.size(); ++i) {
        const std::string& abs = info.absorbingStations[i];
        if (!target.empty() && target != abs)
            info.suggestedFixes.push_back("Route jobs from " + abs + " back to " + target +
                                          " (e.g., P{class}(" + abs + ", " + target + ") = 1.0)");
    }
    return info;
}

/** NODE indices, 1-based, of the absorbing stations. */
template <class T>
std::vector<std::size_t> sn_absorbing_stations(const qn::NetworkStruct<T>& sn) {
    const RoutingErgodicityInfo info = sn_is_routing_ergodic(sn);
    std::vector<std::size_t> out;
    for (std::size_t k = 0; k < info.absorbingStations.size(); ++k)
        for (std::size_t i = 0; i < sn.nodes.size(); ++i)
            if (sn.nodes[i].name == info.absorbingStations[k]) out.push_back(i + 1);
    return out;
}

/**
 * A routing map that makes the network ergodic, by redirecting every absorbing
 * station to TARGETNAME (empty for the default: the first Delay, else the first
 * non-absorbing station).
 *
 * Returns the modified blocks; the caller assigns them to `sn.P` and refreshes.
 * Nothing is mutated here, which is the C++ reading of the reference's "does
 * NOT relink": the modeller sees the repair before it is applied.
 */
template <class T>
std::map<std::pair<std::size_t, std::size_t>, Matrix<T>> sn_make_ergodic(
    const qn::NetworkStruct<T>& sn, const std::string& targetName = std::string()) {
    std::map<std::pair<std::size_t, std::size_t>, Matrix<T>> P = sn.P;
    const RoutingErgodicityInfo info = sn_is_routing_ergodic(sn);
    if (info.isRoutingErgodic || info.absorbingStations.empty()) return P;

    const std::string target =
        targetName.empty() ? detail::default_target(sn, info.absorbingStations) : targetName;
    if (target.empty()) throw InputError("sn_make_ergodic: no suitable target node for routing");
    std::size_t targetIdx = sn.nodes.size();
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].name == target) targetIdx = i;
    if (targetIdx == sn.nodes.size())
        throw InputError("sn_make_ergodic: target node \"" + target + "\" not found");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    for (std::size_t k = 0; k < info.absorbingStations.size(); ++k) {
        std::size_t absIdx = sn.nodes.size();
        for (std::size_t i = 0; i < sn.nodes.size(); ++i)
            if (sn.nodes[i].name == info.absorbingStations[k]) absIdx = i;
        if (absIdx == sn.nodes.size() || absIdx == targetIdx) continue;
        for (typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T>>::iterator it =
                 P.begin();
             it != P.end(); ++it) {
            if (absIdx >= it->second.rows()) continue;
            for (std::size_t j = 0; j < it->second.cols(); ++j) it->second(absIdx, j) = zero;
        }
        for (std::size_t r = 1; r <= sn.nclasses; ++r) {
            typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T>>::iterator it =
                P.find(std::make_pair(r, r));
            if (it != P.end() && absIdx < it->second.rows() && targetIdx < it->second.cols())
                it->second(absIdx, targetIdx) = one;
        }
    }
    return P;
}

}  // namespace sn
}  // namespace line

#endif  // LINE_API_SN_SN_ROUTING_ERGODIC_H
