/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_FJ_MMT_H
#define LINE_SOLVERS_MVA_FJ_MMT_H

/**
 * The fork-join transform SolverMVA applies before solving a layer that
 * contains a Fork.
 *
 * Port of matlab/src/io/@@ModelAdapter/mmt.m. The transform turns a layer with a
 * fork into a plain queueing network with no fork at all, in three moves:
 *
 *   1. THE FORK BECOMES A ROUTER. Its outgoing routing is divided by its
 *      fan-out, so a circulating job takes ONE branch chosen at random rather
 *      than all of them, and the node stops being a Fork -- which also stops
 *      the fork corrections in `chain_visits` from firing, since the reference
 *      transform leaves no Fork in the model for them to key on.
 *   2. THE JOIN BECOMES A DELAY whose per-class service time is the
 *      synchronisation delay `E[X_(k) over branches] * fanOut - mean(branch)`,
 *      re-set on every pass of the fixed point from the branch response times
 *      the previous pass measured. `k` is the branch count on a standard join,
 *      which makes X_(k) the maximum, and the declared quorum on a PARTIAL one;
 *      the floor at zero the quorum needs is in `fj_driver.h`.
 *   3. THE BRANCHES THE JOB DID NOT TAKE ARE CARRIED BY AUXILIARY OPEN CLASSES.
 *      A Source and a Sink are added, and one auxiliary open class is minted
 *      per class of every chain that reaches the fork. Auxiliary tokens arrive
 *      at rate `(fanout - 1) * forkLambda`, enter at the fork, take one branch
 *      each by the same split, and leave at the Sink after the Join. With the
 *      circulating class contributing 1/fanout of a visit to each branch and
 *      the auxiliary stream the remaining (fanout-1)/fanout, every branch
 *      station sees the load of one branch traversal per fork event, which is
 *      what the fork actually generates.
 *
 * Moves 1 and 2 alone are exact only when every station of the layer is an
 * infinite server, because auxiliary tokens then add no waiting. Move 3 is what
 * makes the transform correct at a station that can QUEUE, and it is also what
 * makes the layer MIXED open-closed -- so solver_amvald and solver_mva both
 * have to accept open chains before any of this can be solved.
 *
 * WHAT IS DELIBERATELY REPRODUCED RATHER THAN IMPROVED:
 *
 *   - An auxiliary class is minted for EVERY class of the forked chain, not
 *     only for the classes that reach the fork. The ones that do not get a
 *     disabled arrival and an inert Source -> Sink route; they exist so that
 *     the class indexing of the auxiliary block mirrors the original block
 *     one-for-one, which is what the merge-back keys on.
 *   - The auxiliary routing is CONFINED to the fork-join scope by a class-aware
 *     BFS from the fork that stops at the join. Without it the copied
 *     return-path cycles form recurrent components disconnected from the
 *     Source, which capture the whole stationary mass of the auxiliary chain
 *     and destroy the visit ratios (on lqn_bpmn the reference class switch is
 *     amplified to 27 instead of 1, saturating the layer beyond solvability).
 *
 * SEVERAL FORKS, AND NESTED ONES (2026-07-30). Every Fork of the layer gets its
 * own record: its own join, its own per-class fan-out, and its own block of
 * auxiliary classes, so `fjforkmap[s]` is needed alongside `fjclassmap[s]` to say
 * WHICH fork an auxiliary class stands in for. Two properties then have to be
 * derived rather than assumed, and `fj_sort_forks` derives them (the port of
 * `ModelAdapter.sortForks`):
 *
 *   outer_forks(f, r)  fork f is the OUTERMOST fork on class r's path, i.e. no
 *                      other fork encloses it. Only an outer fork writes a
 *                      synchronisation delay onto the ORIGINAL class r: an inner
 *                      fork's delay is charged while the outer fork's branch is
 *                      being walked, by `fj_find_paths` itself, and writing it
 *                      again on the original class would count it twice.
 *   parent_forks(f)    the fork whose node visits measure how often f fires. For
 *                      an outer fork that is f; for a nested one it is the
 *                      enclosing fork, because a nested fork fires once per
 *                      traversal of the enclosing branch, not once per class
 *                      completion.
 *
 * A fork with NO join is admitted (the reference allows it): there is no
 * synchronisation point, so no delay is computed and `forkLambda` is driven from
 * the fork's own throughput instead of a join's.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/api/fj/fj_ordstat_exp.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"

namespace line {
namespace mva {

using lang::Distrib;
using lang::GlobalConstants;
using lang::JobClassType;
using lang::NodeType;
using lang::SchedStrategy;

/**
 * The transformed layer and the bookkeeping the fixed point needs to drive it
 * and to merge its results back.
 *
 * `fjclassmap[s]` is the ORIGINAL class an auxiliary class s stands for and 0
 * when s is not auxiliary; `fjforkmap[s]` is the index into `forks` of the fork
 * it stands in for; `fanout[s]` is that fork's fan-out for that original class,
 * which is what scales the auxiliary arrival rate.
 */
template <class T>
struct FjMmt {
    /** One record per Fork of the base layer, ascending in node index. */
    struct ForkRec {
        std::size_t node = 0;         ///< 1-based node index, shared with the base layer
        std::size_t joinNode = 0;     ///< 0 when the fork has no join
        std::size_t joinStation = 0;  ///< 0 when the fork has no join
        /**
         * `sn.nodeparam{f}.fanOut`, which MATLAB sets to the fork's tasksPerLink
         * and NOT to its number of output links. A fork built from an LQN
         * POST_AND emits one task per link, so this is 1 and the synchronisation
         * delay is `E[max] - mean`. It is kept as a named quantity because the
         * reference multiplies E[max] by it.
         */
        double fanOut = 1.0;
        /** (nclasses+1) the fork's number of output links for each class. */
        std::vector<std::size_t> origfanout;
        /**
         * (nclasses+1) `outer_forks(f, r)`: this fork is the outermost one on
         * class r's path. See the header. Only then does the synchronisation
         * delay reach the original class, and only then is one computed at all.
         */
        std::vector<bool> outer;
        /**
         * `parent_forks(f)` as an index into `forks`: the fork whose node visits
         * measure how often this one fires. Equal to this fork's own index when
         * it is not nested.
         */
        std::size_t parent = 0;
    };

    qn::NetworkStruct<T> V;  ///< the transformed ("nonfj") layer
    std::vector<ForkRec> forks;
    std::size_t sourceStation = 0;
    std::size_t sourceNode = 0;
    std::size_t sinkNode = 0;
    std::size_t norig = 0;                ///< class count of the base layer
    std::vector<std::size_t> fjclassmap;  ///< (nclasses+1) auxiliary -> original, 0 if not auxiliary
    std::vector<std::size_t> fjforkmap;   ///< (nclasses+1) auxiliary -> index into `forks`
    std::vector<double> fanout;           ///< (nclasses+1)
    std::vector<bool> auxdisabled;        ///< (nclasses+1)
    std::vector<std::size_t> auxclasses;  ///< auxiliary class indices, ascending

    /**
     * Every station that was a Join and is now a zero-service Delay, whether or
     * not a fork claims it. The passes that must hold Immediate at a join, and
     * the merge-back that must keep a join's original-class throughput, iterate
     * this rather than the fork records.
     */
    std::vector<std::size_t> joinStations;

    /**
     * The transform is Heidelberger-Trivedi (`options.config.fork_join='ht'`)
     * rather than MMT. The two share this record because the fixed point that
     * drives them is one function in the reference too (`fjFixedPoint.m`); what
     * differs is the transform (fj_ht.h against fj_mmt.h), the per-pass update
     * -- H-T has no auxiliary open stream whose arrival rate to re-set -- and
     * the synchronisation delays and merge-back the driver applies afterwards.
     */
    bool heidelberger_trivedi = false;
    /**
     * H-T only: the join node -> the auxiliary Delay NODE that carries the time
     * the ORIGINAL class spends outside the fork-join span. `ht.m` adds one such
     * delay per Join, routes the Join into it and it back to the fork, so that an
     * auxiliary token cycles over its own branch alone.
     */
    std::map<std::size_t, std::size_t> auxDelayNode;
    /** H-T only: the same delay as a STATION index. */
    std::map<std::size_t, std::size_t> auxDelayStation;
    /** H-T only: the 1-based branch a given auxiliary class walks; 0 elsewhere. */
    std::vector<std::size_t> auxbranch;

    bool active() const { return !auxclasses.empty(); }

    bool is_join_station(std::size_t st) const {
        return std::find(joinStations.begin(), joinStations.end(), st) != joinStations.end();
    }
};

/**
 * The instant the join fires: E[X_(k)] of independent exponentials with the
 * given means, by inclusion-exclusion. This is the `d0` of the reference,
 * generalised from its maximum to the k-th order statistic so that a quorum
 * join is charged the completion it actually waits for. `k = means.size()` is
 * the ordinary AND-join and reproduces the reference term for term.
 *
 * A zero-length branch no longer degenerates the whole expression to zero, as
 * this port used to have it: it completes instantly, so it neither delays the
 * join nor suppresses the others, which is what the reference computes through
 * its 1/Inf terms.
 */
template <class T>
T fj_expected_ordstat(const std::vector<T>& means, std::size_t k) {
    return fj::fj_ordstat_exp(means, k);
}

/** The AND-join case, i.e. E[max]. */
template <class T>
T fj_expected_max(const std::vector<T>& means) {
    if (means.empty()) return num_traits<T>::from_int(0);
    return fj_expected_ordstat(means, means.size());
}

/**
 * The number of siblings the Join node fires on, out of `nbranches` forked.
 * Port of matlab/src/api/fj/sn_join_quorum.m; a standard join, an absent
 * declaration, a non-positive quorum and a quorum that is not smaller than the
 * sibling count all return `nbranches`, the ordinary AND-join.
 *
 * NOTE the index space: `joindecl` is keyed by the BASE node index, and this
 * port's `JoinDecl` carries ONE quorum per join rather than one per class, so a
 * multiclass model with a different quorum per class cannot be expressed here
 * (the JSON interchange has the same limitation, see network_writer.h).
 */
template <class T>
std::size_t fj_join_quorum(const qn::NetworkStruct<T>& sn, std::size_t joinNode,
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
 * `Exp.fitMean(m)`, including its clamp.
 *
 * The reference builds every synchronisation delay through this factory, and
 * the clamp is load-bearing at both ends: a zero delay would otherwise be an
 * exponential of infinite rate, and it comes out as GlobalConstants.Immediate
 * instead, which is exactly the rate the join carried before the fork loop
 * touched it.
 */
template <class T>
Distrib<T> fj_exp_fit_mean(const T& mean) {
    const T zero = num_traits<T>::from_int(0);
    const T hi = num_traits<T>::from_double(GlobalConstants::Immediate);
    const T lo = num_traits<T>::from_double(GlobalConstants::Zero);
    T rate = hi;
    if (mean > zero) rate = T(num_traits<T>::from_int(1) / mean);
    if (rate < lo) rate = lo;
    if (rate > hi) rate = hi;
    return Distrib<T>::exp_rate(rate);
}

namespace detail {

/**
 * The class-switch-aware node routing r -> s over the arc (i, j).
 *
 * WHY NOT `get_route`. `link()` SYNTHESIZES a ClassSwitch node for every
 * class-switching arc, and the synthesized node's `P` block is class-PRESERVING:
 * the switch itself lives in `csmatrix`, exactly as `refresh_routing` assembles
 * `Peff` (`rtnodes(r, (j-1)K+s) = Pcs(r,s) * P(s,s,i,j)`). A walk over `P` alone
 * therefore stops dead at the first switch -- which is what made the fork
 * reachability BFS below report `reach[fork][r] = false` on every fj_cs_*
 * model, disabling the auxiliary stream and halving the branch load.
 *
 * `route_eff` is not usable here: the transform mutates `P` (the fan-out split,
 * the auxiliary blocks) and `Peff` is only re-derived at the end, so mid-transform
 * it still describes the untransformed layer. This reproduces the one rule that
 * matters on the CURRENT `P`. The reference walks `sn.rtnodes`, which is the same
 * product.
 */
template <class T>
T fj_route_cs(const qn::NetworkStruct<T>& V, std::size_t r, std::size_t s, std::size_t i,
              std::size_t j) {
    const T zero = num_traits<T>::from_int(0);
    const typename std::map<std::size_t, Matrix<T>>::const_iterator it = V.csmatrix.find(i);
    if (it == V.csmatrix.end()) return V.get_route(r, s, i, j);
    const Matrix<T>& C = it->second;
    if (r > C.rows() || s > C.cols()) return zero;
    if (!(C(r - 1, s - 1) > zero)) return zero;
    return T(C(r - 1, s - 1) * V.get_route(s, s, i, j));
}

/** Zero the whole outgoing row of node `nd` in the (r,s) routing block. */
template <class T>
void fj_clear_row(qn::NetworkStruct<T>& V, std::size_t r, std::size_t s, std::size_t nd) {
    auto it = V.P.find(std::make_pair(r, s));
    if (it == V.P.end() || it->second.rows() < nd) return;
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < it->second.cols(); ++j) it->second(nd - 1, j) = zero;
}

/**
 * Port of the local `nestedForks` of `ModelAdapter.sortForks`: walk the class-r
 * routing of the transformed layer from `startNode` to `endNode`, clearing the
 * flag of every Fork node met strictly along the way.
 *
 * `flags[nd]` enters true for a Fork and leaves false when this walk PASSED
 * THROUGH it, which is what makes it an inner fork of the (startNode, endNode)
 * span. The reference intersects the result of each branch (`forks & nested`), so
 * a fork counts as inner as soon as ONE path reaches it.
 *
 * The `seen` set is this port's, not the reference's: `nestedForks` recurses on a
 * routing matrix that may hold a cycle (a call whose mean exceeds one leaves a
 * self-loop) and would not terminate on it, where MATLAB's recursion limit turns
 * that into an error. Marking (node) visited per walk is sound here because the
 * result is a pure reachability question -- revisiting a node cannot clear a flag
 * that the first visit left set.
 */
template <class T>
void fj_nested_forks(const qn::NetworkStruct<T>& V, const std::vector<bool>& isFork,
                     std::size_t startNode, std::size_t endNode, std::size_t cls,
                     std::vector<bool>& flags, std::vector<bool>& seen) {
    if (startNode == endNode) return;
    if (seen[startNode - 1]) return;
    seen[startNode - 1] = true;
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t j = 1; j <= V.nodes.size(); ++j) {
        if (!(V.get_route(cls, cls, startNode, j) > zero)) continue;
        // `isFork` and NOT `V.nodes[j].nodetype`: the reference asks the BASE
        // struct for the node type while walking the TRANSFORMED routing, and by
        // this point every Fork of the transformed struct is a Router, so a type
        // test here would never fire and no fork would ever be found nested.
        if (isFork[j - 1]) flags[j - 1] = false;
        fj_nested_forks(V, isFork, j, endNode, cls, flags, seen);
    }
}

}  // namespace detail

/**
 * Port of `ModelAdapter.sortForks`: fill in `outer` and `parent` on every fork
 * record.
 *
 * Runs on the TRANSFORMED layer, as the reference does
 * (`nonfjmodel.getLinkedRoutingMatrix{r,r}`), and on the ORIGINAL class r, whose
 * routing support the transform leaves alone -- the fork is a Router by now, but
 * with the same out-edges.
 *
 * A fork is outer on class r until some other fork's span is found to pass
 * through it. The parent of every fork found inside f's span becomes f's OWN
 * parent, not f, so a three-deep nest collapses onto the outermost fork exactly
 * as the reference's single assignment does. The reference also writes `parents`
 * at non-fork nodes, where it is never read; this port keeps it per fork record.
 */
template <class T>
void fj_sort_forks(FjMmt<T>& tr) {
    const qn::NetworkStruct<T>& V = tr.V;
    const std::size_t nf = tr.forks.size();
    std::vector<bool> isFork(V.nodes.size(), false);
    for (std::size_t a = 0; a < nf; ++a) {
        tr.forks[a].outer.assign(tr.norig + 1, true);
        tr.forks[a].parent = a;
        isFork[tr.forks[a].node - 1] = true;
    }
    for (std::size_t a = 0; a < nf; ++a) {
        typename FjMmt<T>::ForkRec& F = tr.forks[a];
        if (F.joinNode == 0) continue;
        // The reference walks once per AUXILIARY class of this fork, which
        // includes the classes whose fan-out is zero (they carry a disabled
        // arrival but still have a span).
        for (std::size_t s : tr.auxclasses) {
            if (tr.fjforkmap[s] != a) continue;
            const std::size_t r = tr.fjclassmap[s];
            std::vector<bool> flags(V.nodes.size(), false);
            for (std::size_t b = 0; b < nf; ++b) flags[tr.forks[b].node - 1] = true;
            std::vector<bool> seen(V.nodes.size(), false);
            detail::fj_nested_forks(V, isFork, F.node, F.joinNode, r, flags, seen);
            for (std::size_t b = 0; b < nf; ++b) {
                if (flags[tr.forks[b].node - 1]) continue;
                tr.forks[b].outer[r] = false;
                tr.forks[b].parent = F.parent;
            }
        }
    }
}

/**
 * Widen every explicit ClassSwitch matrix to the auxiliary-expanded class set.
 *
 * AN EXPLICIT ClassSwitch MATRIX IS (nclasses x nclasses), and a fork-join
 * transform has just widened nclasses. `refresh_routing` checks that width and
 * refuses the model outright when it no longer holds, which is what killed
 * SolverMVA on every fork-join model carrying a switch -- and since `link()` now
 * SYNTHESIZES a ClassSwitch node for any class-switching arc, that is every
 * fj_cs_* example ("the class-switch matrix of node 'CS_Fork1_to_Queue1' is not
 * (nclasses x nclasses)"). `fj_tag` and `tag_chain` already widen it at their
 * own augmentation points.
 *
 * An auxiliary class stands for the original it was split from, so it must
 * switch the way that original does, among auxiliaries: C(a,b) = old(x,y) where
 * a,b are the auxiliaries of x,y. Any auxiliary block the map does not pair
 * keeps the identity, which is the only admissible row for a class that never
 * switches.
 */
template <class T>
void fj_widen_csmatrix(FjMmt<T>& tr) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    qn::NetworkStruct<T>& V = tr.V;
    const std::size_t Kaug = V.classes.size();
    for (typename std::map<std::size_t, Matrix<T>>::iterator it = V.csmatrix.begin();
         it != V.csmatrix.end(); ++it) {
        if (it->second.rows() == Kaug && it->second.cols() == Kaug) continue;
        const Matrix<T> old = it->second;
        Matrix<T> C(Kaug, Kaug, zero);
        for (std::size_t x = 0; x < old.rows() && x < Kaug; ++x)
            for (std::size_t y = 0; y < old.cols() && y < Kaug; ++y) C(x, y) = old(x, y);
        for (std::size_t a = old.rows(); a < Kaug; ++a) C(a, a) = one;
        for (std::size_t ai = 0; ai < tr.auxclasses.size(); ++ai)
            for (std::size_t bi = 0; bi < tr.auxclasses.size(); ++bi) {
                const std::size_t a = tr.auxclasses[ai], b = tr.auxclasses[bi];
                if (a > Kaug || b > Kaug || a >= tr.fjclassmap.size() || b >= tr.fjclassmap.size())
                    continue;
                const std::size_t x = tr.fjclassmap[a], y = tr.fjclassmap[b];
                if (x == 0 || y == 0 || x > old.rows() || y > old.cols()) continue;
                if (a == b || num_traits<T>::to_double(old(x - 1, y - 1)) != 0.0)
                    C(a - 1, b - 1) = old(x - 1, y - 1);
            }
        it->second = C;
    }
}

/**
 * Build the transformed layer. Returns an inactive FjMmt when the layer has no
 * fork, in which case the caller solves the layer directly.
 *
 * The auxiliary arrivals are left at GlobalConstants.FineTol, which is what a
 * cold `ModelAdapter.mmt` call sets and what `refreshServicesFromBase` restores
 * at the start of every outer iteration; the fixed point overwrites them with
 * `(fanout - 1) * forkLambda` from its second pass onwards.
 */
template <class T>
FjMmt<T> fj_mmt(const qn::NetworkStruct<T>& L) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    FjMmt<T> tr;
    tr.norig = L.nclasses;

    std::vector<std::size_t> forkNodes;
    for (std::size_t i = 0; i < L.nodes.size(); ++i)
        if (L.nodes[i].nodetype == NodeType::Fork) forkNodes.push_back(i + 1);
    if (forkNodes.empty()) return tr;

    // ---- 1. one record per fork, each with its own join and fan-out --------
    for (std::size_t node : forkNodes) {
        typename FjMmt<T>::ForkRec F;
        F.node = node;
        // tasksPerLink (== sn.nodeparam{f}.fanOut) scales the auxiliary arrival
        // rate and the sync delay; without this it stays 1 and any fork emitting
        // more than one task per link is under-loaded (halved throughput at
        // tasksPerLink=2).
        F.fanOut = L.nodes[node - 1].tasks_per_link;
        for (const std::pair<std::size_t, std::size_t>& p : L.fj) {
            if (p.first != node) continue;
            if (F.joinNode != 0)
                throw UnsupportedError("fj_mmt: layer '" + L.name + "' pairs fork node " +
                                       std::to_string(node) +
                                       " with more than one join station; the reference supports "
                                       "one join per fork");
            F.joinNode = p.second;
        }
        // A fork with NO join is admitted, as in the reference: there is simply
        // no synchronisation point, so no delay is charged and forkLambda is
        // driven from the fork's own throughput instead of a join's.
        // origfanout rationale: see _kb/06-solver-catalog.md (cpp port notes: fj_mmt.h)
        F.origfanout.assign(L.nclasses + 1, 0);
        for (std::size_t r = 1; r <= L.nclasses; ++r) {
            std::vector<bool> dest(L.nodes.size(), false);
            for (std::size_t s = 1; s <= L.nclasses; ++s)
                for (std::size_t j = 1; j <= L.nodes.size(); ++j)
                    if (L.get_route(r, s, node, j) > zero) dest[j - 1] = true;
            for (bool b : dest)
                if (b) ++F.origfanout[r];
        }
        tr.forks.push_back(F);
    }

    tr.V = L;
    qn::NetworkStruct<T>& V = tr.V;

    // Every fork's out-edges are divided by ITS OWN fan-out before any auxiliary
    // class is copied, so a second fork's split is already in place when the
    // first fork's auxiliary block copies the routing. Forks are distinct nodes,
    // so the rows they rewrite are disjoint and the order does not matter.
    for (const typename FjMmt<T>::ForkRec& F : tr.forks) {
        for (std::size_t r = 1; r <= L.nclasses; ++r) {
            if (F.origfanout[r] == 0) continue;
            const T f = num_traits<T>::from_int(static_cast<long>(F.origfanout[r]));
            for (std::size_t s = 1; s <= L.nclasses; ++s)
                for (std::size_t j = 1; j <= L.nodes.size(); ++j) {
                    const T p = V.get_route(r, s, F.node, j);
                    if (p > zero) V.set_route(r, s, F.node, j, T(p / f));
                }
        }
        V.nodes[F.node - 1].nodetype = NodeType::Router;
    }
    V.fj.clear();

    // ---- 2. EVERY Join becomes a zero-service Delay ------------------------
    // The reference converts every Join node of the model, not only the ones
    // paired with a fork, so an unpaired join does not survive as a Join that
    // the inner solver would then reject.
    for (std::size_t j = 1; j <= L.nodes.size(); ++j) {
        if (L.nodes[j - 1].nodetype != NodeType::Join) continue;
        const std::size_t st = V.nodes[j - 1].station;
        if (st == 0)
            throw InputError("fj_mmt: the join node '" + V.nodes[j - 1].name + "' of layer '" +
                             L.name + "' is not a station");
        V.stations[st - 1].nodetype = NodeType::Delay;
        V.stations[st - 1].sched = SchedStrategy::INF;
        V.stations[st - 1].nservers = std::numeric_limits<double>::infinity();
        V.nodes[j - 1].nodetype = NodeType::Delay;
        for (std::size_t k = 1; k <= V.classes.size(); ++k)
            V.set_service(st, k, Distrib<T>::immediate());
        tr.joinStations.push_back(st);
    }
    for (typename FjMmt<T>::ForkRec& F : tr.forks)
        if (F.joinNode != 0) F.joinStation = V.nodes[F.joinNode - 1].station;

    // ---- 3. Source and Sink ----------------------------------------------
    // THE MODEL'S OWN SOURCE IS REUSED WHEN IT HAS ONE. `mmt.m` adds a pair only
    // in the closed case (`if nonfjmodel.hasOpenClasses, source = getSource ...
    // else Source(...) / Sink(...) end`). Adding a SECOND Source to an already
    // open layer left the original one with no `sourceIdx`, so the base chain
    // lost its arrivals altogether: on fj_cs_multi_visits every original-class
    // throughput came back 0 while the auxiliary chain carried the whole flow.
    std::size_t existing_src = 0, existing_sink = 0;
    for (std::size_t i = 1; i <= V.stations.size(); ++i)
        if (V.stations[i - 1].nodetype == NodeType::Source) existing_src = i;
    for (std::size_t j = 1; j <= V.nodes.size(); ++j)
        if (V.nodes[j - 1].nodetype == NodeType::Sink) existing_sink = j;
    if (existing_src != 0 && existing_sink != 0) {
        tr.sourceStation = existing_src;
        tr.sourceNode = V.station_to_node[existing_src - 1];
        tr.sinkNode = existing_sink;
    } else {
        qn::Station<T> src;
        src.name = "Source";
        src.nodetype = NodeType::Source;
        src.sched = SchedStrategy::EXT;
        src.nservers = 1.0;
        tr.sourceStation = V.add_station(src);
        tr.sourceNode = V.station_to_node[tr.sourceStation - 1];
        tr.sinkNode = V.add_node("Sink", NodeType::Sink, false);
    }
    V.sourceIdx = tr.sourceStation;
    V.sinkNode = tr.sinkNode;

    // ---- 4. one auxiliary open class per class of every forked chain, PER
    //         FORK. `fjforkmap` is what tells the two blocks apart.
    const std::size_t nnodes = V.nodes.size();
    for (std::size_t fa = 0; fa < tr.forks.size(); ++fa) {
      const typename FjMmt<T>::ForkRec& F = tr.forks[fa];
      std::vector<bool> forked(L.nclasses + 1, false);
      for (std::size_t c = 0; c < L.nchains; ++c)
        for (std::size_t r = 1; r <= L.nclasses; ++r)
            if (L.nodevisits[c](F.node - 1, r - 1) > zero) forked[r] = true;

      for (std::size_t c = 0; c < L.nchains; ++c) {
        bool any = false;
        for (std::size_t r : L.inchain[c])
            if (forked[r]) any = true;
        if (!any) continue;
        const std::vector<std::size_t>& ic = L.inchain[c];

        // reachability BFS rationale: see _kb/06-solver-catalog.md (cpp port notes: fj_mmt.h)
        const std::size_t refstat = L.classes[ic[0] - 1].refstat;
        const std::size_t refnode = L.station_to_node[refstat - 1];
        std::vector<std::vector<bool>> reach(L.nodes.size(),
                                             std::vector<bool>(L.nclasses + 1, false));
        std::vector<std::pair<std::size_t, std::size_t>> q;
        for (std::size_t r : ic) {
            reach[refnode - 1][r] = true;
            q.emplace_back(refnode, r);
        }
        for (std::size_t h = 0; h < q.size(); ++h) {
            const std::size_t nd = q[h].first, cl = q[h].second;
            for (std::size_t j = 1; j <= L.nodes.size(); ++j)
                for (std::size_t s : ic)
                    if (!reach[j - 1][s] && detail::fj_route_cs(L, cl, s, nd, j) > zero) {
                        reach[j - 1][s] = true;
                        q.emplace_back(j, s);
                    }
        }
        std::vector<std::size_t> aux(ic.size(), 0);
        for (std::size_t a = 0; a < ic.size(); ++a) {
            const std::size_t r = ic[a];
            qn::JobClass xc;
            xc.name = V.classes[r - 1].name + "." + V.nodes[F.node - 1].name;
            xc.type = JobClassType::OPEN;
            xc.population = std::numeric_limits<double>::infinity();
            xc.refstat = tr.sourceStation;
            xc.completes = false;
            xc.is_ref_class = false;
            xc.attr_kind = -1;
            xc.attr_idx = 0;
            xc.prio = V.classes[r - 1].prio;
            aux[a] = V.add_class(xc);

            for (std::size_t i = 1; i <= V.stations.size(); ++i) {
                if (i == tr.sourceStation) continue;
                // EVERY join station takes Immediate for the new class, not only
                // this fork's: the reference's per-auxiliary-class loop keys on
                // the node TYPE, so a second fork's join is Immediate here too.
                if (tr.is_join_station(i)) {
                    V.set_service(i, aux[a], Distrib<T>::immediate());
                    continue;
                }
                V.set_service(i, aux[a], L.service[i - 1][r - 1]);
            }
            const bool disable = (F.origfanout[r] == 0) || !reach[F.node - 1][r];
            V.set_service(tr.sourceStation, aux[a],
                          disable ? Distrib<T>::disabled_dist()
                                  : Distrib<T>::exp_rate(
                                        num_traits<T>::from_double(GlobalConstants::FineTol)));

            tr.fjclassmap.resize(V.classes.size() + 1, 0);
            tr.fjforkmap.resize(V.classes.size() + 1, 0);
            tr.fanout.resize(V.classes.size() + 1, 0.0);
            tr.auxdisabled.resize(V.classes.size() + 1, false);
            tr.fjclassmap[aux[a]] = r;
            tr.fjforkmap[aux[a]] = fa;
            tr.fanout[aux[a]] = static_cast<double>(F.origfanout[r]) * F.fanOut;
            tr.auxdisabled[aux[a]] = disable;
            tr.auxclasses.push_back(aux[a]);
        }

        // copy the (already fork-split) routing onto the auxiliary block
        for (std::size_t x = 0; x < ic.size(); ++x)
            for (std::size_t y = 0; y < ic.size(); ++y)
                for (std::size_t i = 1; i <= nnodes; ++i)
                    for (std::size_t j = 1; j <= nnodes; ++j) {
                        const T p = V.get_route(ic[x], ic[y], i, j);
                        if (p > zero) V.set_route(aux[x], aux[y], i, j, p);
                    }

        // The Source is the auxiliary entry and the Sink its exit. Only THIS
        // fork's join row is cleared; another fork's join is an ordinary node of
        // this span and keeps the routing that was copied onto it.
        for (std::size_t x = 0; x < ic.size(); ++x) {
            for (std::size_t y = 0; y < ic.size(); ++y) {
                detail::fj_clear_row(V, aux[x], aux[y], tr.sourceNode);
                if (F.joinNode != 0) detail::fj_clear_row(V, aux[x], aux[y], F.joinNode);
            }
            if (F.origfanout[ic[x]] > 0) {
                V.set_route(aux[x], aux[x], tr.sourceNode, F.node, one);
                if (F.joinNode != 0)
                    V.set_route(aux[x], aux[x], F.joinNode, tr.sinkNode, one);
            }
        }

        // auxiliary-class scope BFS rationale: see _kb/06-solver-catalog.md (cpp port notes: fj_mmt.h)
        std::vector<std::vector<bool>> vis(nnodes, std::vector<bool>(L.nclasses + 1, false));
        std::vector<std::pair<std::size_t, std::size_t>> bq;
        for (std::size_t x = 0; x < ic.size(); ++x) {
            const std::size_t r = ic[x];
            bool out = false;
            for (std::size_t y = 0; y < ic.size() && !out; ++y)
                for (std::size_t j = 1; j <= nnodes; ++j)
                    if (V.get_route(r, ic[y], F.node, j) > zero) {
                        out = true;
                        break;
                    }
            if (out && !vis[F.node - 1][r]) {
                vis[F.node - 1][r] = true;
                bq.emplace_back(F.node, r);
            }
        }
        for (std::size_t h = 0; h < bq.size(); ++h) {
            const std::size_t cn = bq[h].first, cc = bq[h].second;
            for (std::size_t y = 0; y < ic.size(); ++y) {
                const std::size_t s = ic[y];
                for (std::size_t nd = 1; nd <= nnodes; ++nd)
                    if (!vis[nd - 1][s] && detail::fj_route_cs(V, cc, s, cn, nd) > zero) {
                        vis[nd - 1][s] = true;
                        if (nd != F.joinNode) bq.emplace_back(nd, s);
                    }
            }
        }

        // This fork, ITS join, the Source and the Sink are the transform's own
        // infrastructure for this span: their auxiliary rows were set explicitly
        // above and are never cleared. Another fork is NOT infrastructure here --
        // it is an ordinary node that survives only if this span reaches it,
        // which is what the reference's per-fork `infra_mask` says.
        std::vector<bool> infra(nnodes, false);
        infra[F.node - 1] = true;
        if (F.joinNode != 0) infra[F.joinNode - 1] = true;
        infra[tr.sourceNode - 1] = true;
        infra[tr.sinkNode - 1] = true;
        // per-(node,class) scope clearing rationale: see _kb/06-solver-catalog.md (cpp port notes: fj_mmt.h)
        for (std::size_t nd = 1; nd <= nnodes; ++nd) {
            if (infra[nd - 1]) continue;
            for (std::size_t x = 0; x < ic.size(); ++x) {
                if (vis[nd - 1][ic[x]]) continue;
                for (std::size_t y = 0; y < ic.size(); ++y)
                    detail::fj_clear_row(V, aux[x], aux[y], nd);
            }
        }
        // unreached-auxiliary-class rationale: see _kb/06-solver-catalog.md (cpp port notes: fj_mmt.h)
        for (std::size_t x = 0; x < ic.size(); ++x) {
            bool anyvis = false;
            for (std::size_t nd = 1; nd <= nnodes && !anyvis; ++nd)
                if (vis[nd - 1][ic[x]]) anyvis = true;
            if (!anyvis) V.set_route(aux[x], aux[x], tr.sourceNode, tr.sinkNode, one);
        }
      }
    }

    tr.fjclassmap.resize(V.classes.size() + 1, 0);
    tr.fjforkmap.resize(V.classes.size() + 1, 0);
    tr.fanout.resize(V.classes.size() + 1, 0.0);
    tr.auxdisabled.resize(V.classes.size() + 1, false);
    // The nesting bookkeeping reads the finished auxiliary maps and the
    // transformed routing, so it runs last, exactly where the reference calls
    // sortForks (immediately after mmt returns).
    fj_sort_forks(tr);

    fj_widen_csmatrix(tr);

    // THE EFFECTIVE ROUTING MUST BE RE-DERIVED, not just the chains. Every edit
    // above -- the fan-out division, the auxiliary blocks, the source and sink
    // arcs -- is written through `set_route`, i.e. into `P`; but every consumer
    // reads `route_eff`, which returns `Peff` whenever `Peff` is non-empty. This
    // used to be moot: `refresh_routing` leaves `Peff` EMPTY on a model whose
    // routing is all-PROB and which holds no ClassSwitch node, so `route_eff`
    // fell through to `P` and the edits were visible by accident. A model that
    // DOES hold a ClassSwitch -- now including every model where `link()`
    // synthesizes one -- populates `Peff`, and the whole transform became
    // invisible: on fj_cs_multi_visits the fork's rows kept their undivided
    // fan-out of 2, and the Join's class-2 visit came out 2 instead of 1.
    // `refresh_struct` runs refresh_routing immediately BEFORE refresh_chains;
    // this restores that order for the transformed struct.
    V.refresh_routing();
    V.refresh_chains();
    // AND THE CAPACITIES, which `refresh_struct` derives immediately after the
    // chains for the same reason it derives them at all: a station's buffer is
    // bounded by the population of the chains that reach it, and this transform
    // has just added a Source station and a block of auxiliary OPEN classes. A
    // copy of the layer's `cap`/`classcap` is then both too short and wrong, and
    // on an LN layer it is EMPTY, because `build_layer` stops at refresh_chains
    // and `solve_layer` refreshes the ensemble member rather than this copy.
    // `buffer_size` indexes both vectors by station with no size test -- it is
    // the implementation of Kendall's K, not a permissive default -- so
    // `has_blocking`, and through it the BCMP gate `has_product_form`, read out
    // of range. That was latent until solver_amva began asking a mixed model for
    // `has_product_form` (2026-09-04): lqn_workflows' AND-join layer, whose
    // auxiliary classes are what make it mixed, segfaulted on an empty `cap`.
    V.refresh_capacity();
    // `rt` MUST be re-derived too, and not only the chains. MVA and NC read the
    // per-chain `visits` that refresh_chains rebuilds, so a stale `rt` went
    // unnoticed; the fluid drift reads `sn.rt` directly (stochastic complement
    // over the stateful nodes, which is what absorbs the Router this transform
    // puts where the fork was). Left stale it keeps the ORIGINAL class count, so
    // `rt.rows() != nstateful*nclasses`, every route reads zero and the ODE
    // returns its initial condition with every job parked at the reference
    // station. SolverLN calls refresh_rt for the same reason after it edits a
    // layer.
    V.refresh_rt();
    return tr;
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_FJ_MMT_H
