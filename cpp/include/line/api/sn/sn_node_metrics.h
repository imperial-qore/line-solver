/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_NODE_METRICS_H
#define LINE_API_SN_SN_NODE_METRICS_H

/**
 * Ports of matlab/src/api/sn/sn_get_node_arvr_from_tput.m and
 * sn_get_node_tput_from_tput.m: the NODE-level arrival rate and throughput
 * tables behind `getAvgNode`, `getAvgNodeTable` and the node chain getters.
 *
 * Station metrics are what a solver returns; node metrics are what the user
 * asks for when the model has nodes that are not stations -- a Router, a
 * ClassSwitch, a Cache, a Fork or a Join. Those nodes have no service and no
 * queue, so their only meaningful averages are the rates flowing through them,
 * and those come from the station throughputs plus the NODE-level routing
 * `sn.rtnodes` and the per-chain node visits.
 *
 * REFERENCE STATION INDEX. Both functions read `sn.refstat(c)` with c a CHAIN
 * index, while `sn.refstat` is indexed by CLASS. The two agree whenever the
 * chain's classes are numbered from its own index, which every single-chain
 * and every non-switching model satisfies, and the refresh already refuses a
 * chain whose classes disagree on their reference station. Reproduced as the
 * reference has it: correcting the index here alone would make these tables
 * disagree with MATLAB, the JAR and native Python at once.
 *
 * CACHE HIT AND MISS RATES. The reference splits a Cache node's flow with the
 * ACTUAL hit and miss probabilities the cache fixed point converged to, which
 * this struct does not store on the node -- `refresh_cacheqn_actual_visits`
 * folds them into the visits instead. They are therefore optional arguments:
 * a caller holding the converged probabilities passes them and gets the
 * reference's first branch, and a caller that does not gets its second, the
 * visit-ratio split, which is what the folded visits already encode.
 *
 * ARITHMETIC: field. Sums and quotients only.
 */

#include <cstddef>
#include <map>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

/**
 * The converged cache split, per ORIGINAL class of one Cache node.
 *
 * Empty vectors mean "not available", which selects the reference's
 * visit-ratio branch. `delayed_hit` may stay empty even when the other two are
 * given; it is zero then, as the reference defaults it.
 */
template <class T>
struct CacheActualProb {
    std::vector<T> hit, miss, delayed_hit;
};

namespace detail {

/** True when class `r` (1-based) is a hit or miss class of ANY Cache node. */
template <class T>
bool sn_is_any_cache_class(const qn::NetworkStruct<T>& sn, std::size_t r) {
    for (const auto& kv : sn.nodeparam) {
        const qn::CacheParam<T>& cp = kv.second;
        for (std::size_t x : cp.hitclass)
            if (x == r) return true;
        for (std::size_t x : cp.missclass)
            if (x == r) return true;
    }
    return false;
}

/**
 * `sum(TN(refstat, inchain))`, the chain's total throughput at its reference
 * station, over the classes that HAVE a reference visit.
 *
 * THIS SUM AND `sn_chain_refvisits` ARE THE TWO HALVES OF ONE RATIO -- a rate
 * per reference visit -- so they have to run over the same classes or the
 * quotient is not that rate. The restriction is inert on every model but a
 * cache one, because a class with no visit at a station has no throughput there
 * either. A CACHE IS THE ONE EXCEPTION: its over-route sends a hit or miss
 * class THROUGH a station that never serves it, `getAvgHandles` keeps the Tput
 * handle of a cache class so that pass-through flow IS reported, and the
 * refreshed struct gives the same pair no visit. The reference never has the
 * term at all -- its implicit ClassSwitch absorbs the switch before the station
 * -- so counting it here alone scaled every cache node rate by the ratio of the
 * two sums: exactly 2x on cache_replc_fifo, whose Delay carries the whole
 * hit+miss flow a second time.
 */
template <class T>
T sn_chain_tput(const qn::NetworkStruct<T>& sn, const Matrix<T>& TN, std::size_t c,
                std::size_t refstat, const std::vector<std::size_t>& inchain) {
    T acc = num_traits<T>::from_int(0);
    if (refstat == 0 || refstat > sn.nstations) return acc;
    const T zero = num_traits<T>::from_int(0);
    const bool have_visits = c < sn.visits.size();
    const std::size_t isf = have_visits ? sn.stateful_of_station(refstat) - 1 : 0;
    for (std::size_t r : inchain) {
        if (have_visits && sn.visits[c](isf, r - 1) == zero && sn_is_any_cache_class(sn, r))
            continue;
        acc = T(acc + TN(refstat - 1, r - 1));
    }
    return acc;
}

/** `sum(sn.visits{c}(stationToStateful(refstat), inchain))`. */
template <class T>
T sn_chain_refvisits(const qn::NetworkStruct<T>& sn, std::size_t c, std::size_t refstat,
                     const std::vector<std::size_t>& inchain) {
    T acc = num_traits<T>::from_int(0);
    if (refstat == 0 || refstat > sn.nstations || c >= sn.visits.size()) return acc;
    const std::size_t isf = sn.stateful_of_station(refstat) - 1;
    for (std::size_t r : inchain) acc = T(acc + sn.visits[c](isf, r - 1));
    return acc;
}

/**
 * `sn.refstat(c)` as the reference spells it: the CLASS-indexed reference
 * station read at the chain's index. See the file header.
 */
template <class T>
std::size_t sn_chain_refstat(const qn::NetworkStruct<T>& sn, std::size_t c) {
    if (c >= sn.classes.size()) return 0;
    return sn.classes[c].refstat;
}

/** True when class `r` (1-based) is one of the node's hit or miss classes. */
inline bool sn_is_cache_class(const std::vector<std::size_t>& hit,
                              const std::vector<std::size_t>& miss, std::size_t r) {
    for (std::size_t x : hit)
        if (x == r) return true;
    for (std::size_t x : miss)
        if (x == r) return true;
    return false;
}

}  // namespace detail

/**
 * Port of sn_get_node_arvr_from_tput.
 *
 * `AN` is the STATION arrival-rate table (sn_get_arvr_from_tput). Stations
 * copy their own row; every other node takes the chain's reference throughput
 * scaled by its node visits.
 */
template <class T>
Matrix<T> sn_get_node_arvr_from_tput(const qn::NetworkStruct<T>& sn, const Matrix<T>& TN,
                                     const Matrix<T>& AN) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t I = sn.nodes.size(), C = sn.nchains, M = sn.nstations, R = sn.nclasses;
    Matrix<T> ANn(I, R, zero);
    if (TN.rows() == 0 || AN.rows() == 0) return ANn;

    for (std::size_t ist = 0; ist < M; ++ist) {
        const std::size_t ind = sn.station_to_node[ist];
        if (ind == 0) continue;
        for (std::size_t r = 0; r < R; ++r) ANn(ind - 1, r) = AN(ist, r);
    }
    for (std::size_t ind = 1; ind <= I; ++ind) {
        if (sn.nodes[ind - 1].nodetype == qn::NodeType::Source) continue;
        const bool is_cache = sn.nodes[ind - 1].nodetype == qn::NodeType::Cache;
        const bool is_station = sn.nodes[ind - 1].station != 0;
        for (std::size_t c = 0; c < C; ++c) {
            const std::vector<std::size_t>& inchain = sn.inchain[c];
            const std::size_t refstat = detail::sn_chain_refstat(sn, c);
            const T den = detail::sn_chain_refvisits(sn, c, refstat, inchain);
            if (!(den > zero)) continue;
            const T tot = detail::sn_chain_tput(sn, TN, c, refstat, inchain);
            for (std::size_t r : inchain) {
                if (is_cache) {
                    typename std::map<std::size_t, qn::CacheParam<T>>::const_iterator it =
                        sn.nodeparam.find(ind);
                    if (it != sn.nodeparam.end() &&
                        detail::sn_is_cache_class(it->second.hitclass, it->second.missclass, r))
                        continue;
                } else if (is_station) {
                    continue;
                }
                ANn(ind - 1, r - 1) = T(sn.nodevisits[c](ind - 1, r - 1) / den * tot);
            }
        }
    }
    return ANn;
}

/**
 * Port of sn_get_node_tput_from_tput.
 *
 * Three passes, in the reference's order. A Cache node first splits the
 * arriving flow into its hit and miss classes; every station then overwrites
 * its own node row with its station throughput; finally every non-station node
 * that is not a Sink accumulates what `sn.rtnodes` carries out of it.
 *
 * `cache_prob` maps a 1-based Cache NODE index to its converged split; an
 * absent or empty entry selects the visit-ratio branch.
 */
template <class T>
Matrix<T> sn_get_node_tput_from_tput(
    const qn::NetworkStruct<T>& sn, const Matrix<T>& TN, const Matrix<T>& ANn,
    const std::map<std::size_t, CacheActualProb<T>>& cache_prob =
        std::map<std::size_t, CacheActualProb<T>>()) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t I = sn.nodes.size(), C = sn.nchains, M = sn.nstations, R = sn.nclasses;
    Matrix<T> TNn(I, R, zero);
    if (TN.rows() == 0) return TNn;

    // pass 1: the Cache hit and miss split
    for (std::size_t ind = 1; ind <= I; ++ind) {
        if (sn.nodes[ind - 1].nodetype != qn::NodeType::Cache) continue;
        typename std::map<std::size_t, qn::CacheParam<T>>::const_iterator np =
            sn.nodeparam.find(ind);
        if (np == sn.nodeparam.end()) continue;
        const std::vector<std::size_t>& hitclass = np->second.hitclass;
        const std::vector<std::size_t>& missclass = np->second.missclass;
        typename std::map<std::size_t, CacheActualProb<T>>::const_iterator cp =
            cache_prob.find(ind);
        const bool have_actual = cp != cache_prob.end() && !cp->second.hit.empty();
        for (std::size_t c = 0; c < C; ++c) {
            const std::vector<std::size_t>& inchain = sn.inchain[c];
            const std::size_t refstat = detail::sn_chain_refstat(sn, c);
            const T tot = detail::sn_chain_tput(sn, TN, c, refstat, inchain);
            for (std::size_t r : inchain) {
                if (have_actual) {
                    for (std::size_t o = 0; o < hitclass.size(); ++o) {
                        T arv = zero;
                        if (ANn.rows() == I && o < R) {
                            arv = ANn(ind - 1, o);
                        } else {
                            bool in = false;
                            for (std::size_t q : inchain)
                                if (q == o + 1) in = true;
                            if (in && refstat >= 1 && refstat <= M) arv = TN(refstat - 1, o);
                        }
                        if (hitclass[o] == r && o < cp->second.hit.size()) {
                            T p = cp->second.hit[o];
                            if (o < cp->second.delayed_hit.size())
                                p = T(p + cp->second.delayed_hit[o]);
                            TNn(ind - 1, r - 1) = T(TNn(ind - 1, r - 1) + arv * p);
                        } else if (o < missclass.size() && missclass[o] == r &&
                                   o < cp->second.miss.size()) {
                            TNn(ind - 1, r - 1) = T(TNn(ind - 1, r - 1) + arv * cp->second.miss[o]);
                        }
                    }
                } else if (detail::sn_is_cache_class(hitclass, missclass, r)) {
                    const T den = detail::sn_chain_refvisits(sn, c, refstat, inchain);
                    if (den > zero)
                        TNn(ind - 1, r - 1) =
                            T(sn.nodevisits[c](ind - 1, r - 1) / den * tot);
                }
            }
        }
    }

    // pass 2: every station's own throughput
    for (std::size_t ist = 0; ist < M; ++ist) {
        const std::size_t ind = sn.station_to_node[ist];
        if (ind == 0) continue;
        for (std::size_t r = 0; r < R; ++r) TNn(ind - 1, r) = TN(ist, r);
    }

    // pass 3: what the node routing carries out of a non-station node
    if (ANn.rows() != I || sn.rtnodes.rows() != I * R) return TNn;
    for (std::size_t ind = 1; ind <= I; ++ind) {
        const qn::NodeType nt = sn.nodes[ind - 1].nodetype;
        if (nt == qn::NodeType::Sink) continue;
        if (nt == qn::NodeType::Source) {
            const std::size_t ist = sn.nodes[ind - 1].station;
            if (ist != 0)
                for (std::size_t s = 0; s < R; ++s) TNn(ind - 1, s) = TN(ist - 1, s);
            continue;
        }
        // a Join is accumulated even though it IS a station: the reference's
        // Join arm has no isstation guard, and a Join's throughput is the rate
        // its outgoing edges carry, not the rate its (auxiliary) station serves
        if (nt != qn::NodeType::Join && sn.nodes[ind - 1].station != 0) continue;
        for (std::size_t c = 0; c < C; ++c) {
            const std::vector<std::size_t>& inchain = sn.inchain[c];
            for (std::size_t r : inchain) {
                bool anystateful = false;
                if (c < sn.visits.size())
                    for (std::size_t a = 0; a < sn.visits[c].rows(); ++a)
                        if (sn.visits[c](a, r - 1) != zero) anystateful = true;
                if (!anystateful) continue;
                for (std::size_t s : inchain)
                    for (std::size_t jnd = 1; jnd <= I; ++jnd) {
                        if (nt == qn::NodeType::Cache && ind == jnd) continue;
                        TNn(ind - 1, s - 1) =
                            T(TNn(ind - 1, s - 1) +
                              ANn(ind - 1, r - 1) *
                                  sn.rtnodes((ind - 1) * R + (r - 1), (jnd - 1) * R + (s - 1)));
                    }
            }
        }
    }
    return TNn;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_NODE_METRICS_H
