/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_FJ_HT_H
#define LINE_SOLVERS_MVA_FJ_HT_H

/**
 * The Heidelberger-Trivedi fork-join transform, `options.config.fork_join='ht'`.
 *
 * Port of matlab/src/io/@@ModelAdapter/ht.m (Heidelberger and Trivedi, "Queueing
 * network models for parallel processing with asynchronous tasks", IEEE TC
 * C-31(11), 1982). It is the second arm of the fork-join fixed point that
 * `fj_driver.h` drives, beside the MMT transform of `fj_mmt.h`, and it answers a
 * different question: where MMT keeps the circulating job on ONE branch and
 * carries the remaining branches by auxiliary OPEN classes, H-T sends the
 * circulating job STRAIGHT PAST the branches and gives every branch its own
 * auxiliary CLOSED class, one per (forked class, branch), whose population
 * matches the original's. The four moves are:
 *
 *   1. THE FORK BECOMES A ROUTER, and the ORIGINAL classes are routed from it
 *      directly to the join. The original job therefore spends no time on the
 *      branches; the whole fork-join span is charged to it as one delay at the
 *      join.
 *   2. THE JOIN BECOMES A DELAY. Its service for the original class is the
 *      instant the join fires, `E[X_(k)] * fanOut`, and for an auxiliary class
 *      the residual `E[X_(k)] - R_branch` that branch still waits at the
 *      synchronisation point.
 *   3. AN AUXILIARY DELAY IS ADDED PER JOIN, "Auxiliary Delay - <join>". The
 *      join routes into it and it routes back to the fork, so an auxiliary token
 *      cycles fork -> its own branch -> join -> auxiliary delay -> fork. Its
 *      service for an auxiliary class is the response time the ORIGINAL class
 *      accumulates OUTSIDE the span, which is what keeps the auxiliary token's
 *      cycle time equal to the original's.
 *   4. ONE CLOSED AUXILIARY CLASS PER BRANCH AND PER FORKED CLASS, of the
 *      original's population, referencing the auxiliary delay. It carries the
 *      original's service demand at every station of its branch, so each branch
 *      is loaded by a population equal to the original's, which is what a fork
 *      actually generates.
 *
 * WHAT IT REFUSES, by name and for the same reasons the reference does:
 *   - `tasksPerLink > 1`. The transform has no way to send w identical tasks
 *     down a link, since a branch carries exactly one auxiliary class.
 *   - an OPEN class through the fork. The auxiliary class is a ClosedClass of
 *     the original's population, and an open class has none.
 *   - a fork with no join. H-T charges the whole span at the synchronisation
 *     point, so without one there is nothing to charge.
 *
 * The synchronisation delays themselves, and the merge-back of the auxiliary
 * columns, live in `fj_driver.h`: the reference keeps both arms in one
 * `fjFixedPoint`, and so does this port.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/fj_mmt.h"
#include "line/util/error.h"

namespace line {
namespace mva {

namespace detail {

/**
 * `P{a,b} = P{r,s}` of the reference: copy a whole routing block onto another
 * class pair. An absent source block is an all-zero one, in which case nothing
 * is copied and the destination keeps whatever the caller writes into it -- the
 * same state the reference reaches, its cell array holding a zero matrix there.
 */
template <class T>
void fj_copy_block(qn::NetworkStruct<T>& V, std::size_t r, std::size_t s, std::size_t a,
                   std::size_t b) {
    const typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T>>::const_iterator it =
        V.P.find(std::make_pair(r, s));
    if (it == V.P.end()) return;
    const Matrix<T> src = it->second;  // by value: set_route below may rehash V.P
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < src.rows(); ++i)
        for (std::size_t j = 0; j < src.cols(); ++j)
            if (src(i, j) > zero || src(i, j) < zero) V.set_route(a, b, i + 1, j + 1, src(i, j));
}

/** `P{r,s}(dst,:) = P{r,s}(src,:)` for every class pair the routing holds. */
template <class T>
void fj_copy_row_all_blocks(qn::NetworkStruct<T>& V, std::size_t src, std::size_t dst) {
    std::vector<std::pair<std::size_t, std::size_t>> keys;
    for (typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T>>::const_iterator it =
             V.P.begin();
         it != V.P.end(); ++it)
        keys.push_back(it->first);
    for (std::size_t k = 0; k < keys.size(); ++k) {
        typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T>>::iterator it =
            V.P.find(keys[k]);
        if (it == V.P.end() || it->second.rows() < src || it->second.rows() < dst) continue;
        for (std::size_t j = 0; j < it->second.cols(); ++j)
            it->second(dst - 1, j) = it->second(src - 1, j);
    }
}

}  // namespace detail

/**
 * Build the H-T transform of `L`. Returns an inactive record when the layer
 * holds no fork, in which case the caller solves the layer directly.
 */
template <class T>
FjMmt<T> fj_ht(const qn::NetworkStruct<T>& L) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    FjMmt<T> tr;
    tr.heidelberger_trivedi = true;
    tr.norig = L.nclasses;

    std::vector<std::size_t> forkNodes;
    for (std::size_t i = 0; i < L.nodes.size(); ++i)
        if (L.nodes[i].nodetype == NodeType::Fork) forkNodes.push_back(i + 1);
    if (forkNodes.empty()) return tr;

    // ---- 1. one record per fork, with its join, its fan-out and the ORDERED
    //         head node of each of its branches, read off the BASE routing.
    // The reference reads the heads from `outputStrategy{r}{3}`, which is the
    // order the links were declared in; this port has only the routing matrix,
    // so it takes them in ascending node index. The branches are interchangeable
    // in everything that follows -- each auxiliary class walks one of them and
    // the order statistic is over the whole set -- so the ordering names the
    // classes, it does not change the answer.
    std::vector<std::vector<std::vector<std::size_t>>> heads(forkNodes.size());
    for (std::size_t a = 0; a < forkNodes.size(); ++a) {
        typename FjMmt<T>::ForkRec F;
        F.node = forkNodes[a];
        F.fanOut = L.nodes[F.node - 1].tasks_per_link;
        if (F.fanOut > 1.0)
            throw UnsupportedError(
                "fj_ht: the fork node '" + L.nodes[F.node - 1].name + "' of model '" + L.name +
                "' sends more than one task per link (tasksPerLink=" +
                std::to_string(static_cast<long>(F.fanOut + 0.5)) +
                "); multiple tasks per link are not supported in H-T, use fork_join='mmt'");
        for (std::size_t p = 0; p < L.fj.size(); ++p) {
            if (L.fj[p].first != F.node) continue;
            if (F.joinNode != 0)
                throw UnsupportedError("fj_ht: model '" + L.name + "' pairs fork node " +
                                       std::to_string(F.node) +
                                       " with more than one join station; the reference supports "
                                       "one join per fork");
            F.joinNode = L.fj[p].second;
        }
        if (F.joinNode == 0)
            throw UnsupportedError(
                "fj_ht: the fork node '" + L.nodes[F.node - 1].name + "' of model '" + L.name +
                "' has no join; the H-T method charges the whole fork-join span at the "
                "synchronisation point, so it needs one. Use fork_join='mmt'");
        F.origfanout.assign(L.nclasses + 1, 0);
        heads[a].assign(L.nclasses + 1, std::vector<std::size_t>());
        for (std::size_t r = 1; r <= L.nclasses; ++r) {
            for (std::size_t j = 1; j <= L.nodes.size(); ++j) {
                bool linked = false;
                for (std::size_t s = 1; s <= L.nclasses && !linked; ++s)
                    if (L.get_route(r, s, F.node, j) > zero) linked = true;
                if (linked) heads[a][r].push_back(j);
            }
            F.origfanout[r] = heads[a][r].size();
        }
        tr.forks.push_back(F);
    }

    tr.V = L;
    qn::NetworkStruct<T>& V = tr.V;

    // ---- 2. every fork becomes a Router. Its out-edges are NOT divided as the
    // MMT transform divides them: the original class is rerouted straight to the
    // join below and never takes a branch at all.
    for (std::size_t a = 0; a < tr.forks.size(); ++a)
        V.nodes[tr.forks[a].node - 1].nodetype = NodeType::Router;
    V.fj.clear();

    // ---- 3. every Join becomes a zero-service Delay and gains its auxiliary
    // Delay. As in the MMT transform, EVERY Join of the model is converted, not
    // only the ones a fork claims, so an unpaired join does not survive as a
    // Join the inner solver would then reject.
    for (std::size_t j = 1; j <= L.nodes.size(); ++j) {
        if (L.nodes[j - 1].nodetype != NodeType::Join) continue;
        const std::size_t st = V.nodes[j - 1].station;
        if (st == 0)
            throw InputError("fj_ht: the join node '" + V.nodes[j - 1].name + "' of model '" +
                             L.name + "' is not a station");
        V.stations[st - 1].nodetype = NodeType::Delay;
        V.stations[st - 1].sched = SchedStrategy::INF;
        V.stations[st - 1].nservers = std::numeric_limits<double>::infinity();
        V.nodes[j - 1].nodetype = NodeType::Delay;
        for (std::size_t k = 1; k <= V.classes.size(); ++k)
            V.set_service(st, k, Distrib<T>::immediate());
        tr.joinStations.push_back(st);

        qn::Station<T> ad;
        ad.name = "Auxiliary Delay - " + V.nodes[j - 1].name;
        ad.nodetype = NodeType::Delay;
        ad.sched = SchedStrategy::INF;
        ad.nservers = std::numeric_limits<double>::infinity();
        const std::size_t adstat = V.add_station(ad);
        const std::size_t adnode = V.station_to_node[adstat - 1];
        tr.auxDelayStation[j] = adstat;
        tr.auxDelayNode[j] = adnode;
        for (std::size_t k = 1; k <= V.classes.size(); ++k)
            V.set_service(adstat, k, Distrib<T>::immediate());

        // The auxiliary delay INHERITS the join's out-edges and the join is then
        // routed into it, so every class leaves the span through the new delay.
        // The reference copies the row in every block and rewrites only the
        // DIAGONAL one, which is what keeps a class switch declared at the join
        // pointing where it pointed.
        detail::fj_copy_row_all_blocks(V, j, adnode);
        for (std::size_t r = 1; r <= V.classes.size(); ++r) {
            detail::fj_clear_row(V, r, r, j);
            V.set_route(r, r, j, adnode, one);
        }
    }
    for (std::size_t a = 0; a < tr.forks.size(); ++a)
        tr.forks[a].joinStation = V.nodes[tr.forks[a].joinNode - 1].station;

    // ---- 4. one auxiliary CLOSED class per forked class and per branch -------
    tr.fjclassmap.assign(V.classes.size() + 1, 0);
    tr.fjforkmap.assign(V.classes.size() + 1, 0);
    tr.fanout.assign(V.classes.size() + 1, 0.0);
    tr.auxdisabled.assign(V.classes.size() + 1, false);
    tr.auxbranch.assign(V.classes.size() + 1, 0);

    for (std::size_t fa = 0; fa < tr.forks.size(); ++fa) {
        const typename FjMmt<T>::ForkRec& F = tr.forks[fa];
        const std::size_t adstat = tr.auxDelayStation[F.joinNode];
        const std::size_t adnode = tr.auxDelayNode[F.joinNode];

        // the classes that reach this fork, i.e. `Vnodes(f,:) > 0`
        std::vector<bool> forked(L.nclasses + 1, false);
        for (std::size_t c = 0; c < L.nchains; ++c)
            for (std::size_t r = 1; r <= L.nclasses; ++r)
                if (L.nodevisits[c](F.node - 1, r - 1) > zero) forked[r] = true;

        for (std::size_t c = 0; c < L.nchains; ++c) {
            const std::vector<std::size_t>& ic = L.inchain[c];
            bool any = false;
            for (std::size_t x = 0; x < ic.size(); ++x)
                if (forked[ic[x]]) any = true;
            if (!any) continue;

            // aux[r] holds this chain's auxiliary classes for original class r,
            // in branch order; empty for a class that does not reach the fork.
            std::vector<std::vector<std::size_t>> aux(L.nclasses + 1);
            for (std::size_t x = 0; x < ic.size(); ++x) {
                const std::size_t r = ic[x];
                if (!(L.nodevisits[c](F.node - 1, r - 1) > zero)) continue;
                if (std::isinf(L.classes[r - 1].population))
                    throw UnsupportedError(
                        "fj_ht: class '" + L.classes[r - 1].name + "' of model '" + L.name +
                        "' is open and reaches a fork; the H-T method can be used only on closed "
                        "models, use fork_join='mmt'");
                for (std::size_t par = 1; par <= F.origfanout[r]; ++par) {
                    qn::JobClass xc;
                    xc.name = V.classes[r - 1].name + "." + V.nodes[F.node - 1].name + ".B" +
                              std::to_string(par);
                    xc.type = JobClassType::CLOSED;
                    // tasksPerLink * population, and tasksPerLink is 1 here: a
                    // fork emitting more was refused above.
                    xc.population = F.fanOut * L.classes[r - 1].population;
                    xc.refstat = adstat;
                    xc.completes = true;
                    xc.is_ref_class = false;
                    xc.attr_kind = -1;
                    xc.attr_idx = 0;
                    xc.prio = 0;  // `ClosedClass(..., 0)` of the reference
                    const std::size_t s = V.add_class(xc);

                    // The service demand of the original class, at the BASE
                    // stations only: a Join is Immediate (the driver overwrites
                    // it with the residual synchronisation delay), a Source or a
                    // Fork is a no-op, and every OTHER auxiliary delay is left
                    // disabled, which is where the reference's `1:sn.nnodes`
                    // over the BASE struct leaves it.
                    for (std::size_t i = 1; i <= L.stations.size(); ++i) {
                        const NodeType nt = L.stations[i - 1].nodetype;
                        if (nt == NodeType::Join) {
                            V.set_service(i, s, Distrib<T>::immediate());
                        } else if (nt == NodeType::Source || nt == NodeType::Fork) {
                            // no-op
                        } else {
                            V.set_service(i, s, L.service[i - 1][r - 1]);
                        }
                    }
                    V.set_service(adstat, s, Distrib<T>::immediate());

                    tr.fjclassmap.resize(V.classes.size() + 1, 0);
                    tr.fjforkmap.resize(V.classes.size() + 1, 0);
                    tr.fanout.resize(V.classes.size() + 1, 0.0);
                    tr.auxdisabled.resize(V.classes.size() + 1, false);
                    tr.auxbranch.resize(V.classes.size() + 1, 0);
                    tr.fjclassmap[s] = r;
                    tr.fjforkmap[s] = fa;
                    tr.fanout[s] = static_cast<double>(F.origfanout[r]) * F.fanOut;
                    tr.auxbranch[s] = par;
                    tr.auxclasses.push_back(s);
                    aux[r].push_back(s);
                }
            }

            // ---- the routing of the auxiliary classes, and of the originals --
            for (std::size_t x = 0; x < ic.size(); ++x) {
                const std::size_t r = ic[x];
                if (aux[r].empty()) continue;
                for (std::size_t y = 0; y < ic.size(); ++y) {
                    const std::size_t s = ic[y];
                    if (aux[s].empty()) continue;
                    for (std::size_t par = 1; par <= aux[r].size() && par <= aux[s].size(); ++par) {
                        const std::size_t a = aux[r][par - 1], b = aux[s][par - 1];
                        detail::fj_copy_block(V, r, s, a, b);
                        // the auxiliary token takes ITS OWN branch, waits at the
                        // join, and returns to the fork through the auxiliary
                        // delay that carries the rest of the original's cycle
                        detail::fj_clear_row(V, a, b, F.node);
                        V.set_route(a, b, F.node, heads[fa][r][par - 1], one);
                        V.set_route(a, b, F.joinNode, adnode, one);
                        detail::fj_clear_row(V, a, b, adnode);
                        V.set_route(a, b, adnode, F.node, one);
                    }
                    // The original class is routed straight to the join, so it
                    // does not interfere with the auxiliary tokens on the
                    // branches; the whole span is charged to it as the join's
                    // own service time.
                    detail::fj_clear_row(V, r, s, F.node);
                    V.set_route(r, s, F.node, F.joinNode, one);
                }
            }
        }
    }

    tr.fjclassmap.resize(V.classes.size() + 1, 0);
    tr.fjforkmap.resize(V.classes.size() + 1, 0);
    tr.fanout.resize(V.classes.size() + 1, 0.0);
    tr.auxdisabled.resize(V.classes.size() + 1, false);
    tr.auxbranch.resize(V.classes.size() + 1, 0);

    // `outer` and `parent` are read by the MMT arm of the driver only, but they
    // are filled here too so that a record built by either transform answers the
    // same questions. A fork whose span holds another fork is inner on that
    // class, exactly as in `fj_mmt`.
    fj_sort_forks(tr);
    fj_widen_csmatrix(tr);

    // The effective routing has to be re-derived, not just the chains: every
    // edit above is written through `set_route`, i.e. into `P`, while every
    // consumer reads `route_eff`, which returns `Peff` whenever `Peff` is
    // non-empty. See the same note in `fj_mmt`.
    V.refresh_routing();
    V.refresh_chains();
    // And the capacities, for the reason spelled out at the tail of `fj_mmt`:
    // this transform adds a delay station and a block of auxiliary classes, so
    // the copied `cap`/`classcap` no longer span the stations `buffer_size`
    // indexes.
    V.refresh_capacity();
    return tr;
}

/**
 * `options.config.fork_join` -> the transform it names.
 *
 * The reference spells the same switch twice (`fjFixedPoint.m:47` and `:130`);
 * this port resolves it once, at the point the transform is built. An unknown
 * name is refused rather than silently taking the default, since answering with
 * the MMT transform under another method's name is the one failure mode a
 * method switch must not have.
 */
template <class T>
FjMmt<T> fj_fork_join_transform(const qn::NetworkStruct<T>& L, const std::string& method) {
    if (method.empty() || method == "default" || method == "mmt" || method == "fjt")
        return fj_mmt(L);
    if (method == "ht" || method == "heidelberger-trivedi") return fj_ht(L);
    throw InputError("fork-join: '" + method +
                     "' is not a fork-join method; options.config.fork_join is one of "
                     "'default', 'mmt', 'fjt' (the MMT transform) or 'ht', "
                     "'heidelberger-trivedi'");
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_FJ_HT_H
