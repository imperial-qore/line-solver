/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_FJ_VISITS_SPN_H
#define LINE_API_SN_SN_FJ_VISITS_SPN_H

/**
 * Fork-join node visit ratios, via the auxiliary closed SPN.
 *
 * Port of matlab/src/api/sn/sn_fj_visits_spn.m (twins: `SnFjVisitsSpn.java`,
 * `python/line_solver/api/sn/sn_fj_visits_spn.py`).
 *
 * WHAT IT COMPUTES. For each class that passes through a fork-join pair, the
 * reference builds an auxiliary closed stochastic Petri net that carries the
 * fork/join synchronization exactly -- one Place per station, a "done" Place
 * per branch feeding a Join, an immediate Join transition consuming one token
 * from each branch, and B tokens circulating, B being the largest leaf count
 * over the outermost forks -- and reads the per-Place throughputs as the visit
 * ratios, normalized so the chain's reference station is one.
 *
 * THE SPN SOLVE IS NOT WHAT PRODUCES THE ANSWER, but the rule it produces is
 * NOT the uniform one `sn_fj_visits_spn.m:113` claims. That comment says the net
 * is population preserving, so every station Place fires at the same rate; the
 * reference's own CTMC solve says otherwise, and it is the solve that defines
 * the function. The pre-fork transition consumes all B tokens at once and the
 * Join returns them, so a station INSIDE a fork-join region fires once per B
 * firings of the cycle: normalized on the reference station, a station outside
 * the region is 1, a station inside it is 1/B, and a Fork or a Join, which
 * holds no Place, is 0. Measured against MATLAB and the JAR on a two-branch, a
 * three-branch, a nested, a chained-branch, a pre/post-fork-station and a
 * two-class model: 1, 0.5 and 1/3 exactly where predicted.
 *
 * This port evaluates that rule instead of enumerating a state space whose size
 * is exponential in B. WHAT IS LOST, stated plainly: MATLAB's solve would DETECT
 * a construction that violates the assumptions (an unbalanced Join, a fork whose
 * leaves are not conserved), where this port assumes them. That is why the
 * structural preconditions are CHECKED here -- an unresolvable fork, or a Join
 * with no branch feeding it, is refused by name.
 *
 * ARITHMETIC: field. No transcendental function is involved; before
 * normalization the answer is exactly 0, 1 or the unit fraction 1/B.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

namespace fjvdetail {

/**
 * Recursively resolve a Fork's destinations down to station nodes.
 *
 * A fork branch may open another fork (nested fork-join), so the leaves are
 * the stations reachable without crossing a Join.
 */
template <class T>
void resolve_fork_dests(const qn::NetworkStruct<T>& sn, const Matrix<T>& P_r,
                        const std::vector<bool>& visited, std::size_t fork_nd, std::size_t cls,
                        double weight, std::vector<std::size_t>* out, std::vector<double>* wout,
                        std::vector<bool>* seen) {
    const std::size_t I = sn.nodes.size();
    if ((*seen)[fork_nd]) return;  // a routing cycle through forks would not terminate
    (*seen)[fork_nd] = true;
    const qn::ForkParam<T>* fk = sn.fork_param_of(fork_nd + 1);
    const bool variable = fk != 0;
    for (std::size_t j = 0; j < I; ++j) {
        if (!(num_traits<T>::to_double(P_r(fork_nd, j)) > 0.0) || !visited[j]) continue;
        // EXPECTED tasks this link carries. A plain fork gives every link 1, so
        // the weighted leaf count collapses to the leaf count it always was.
        double w = weight;
        if (variable && j < fk->fan_out_link.rows() && cls < fk->fan_out_link.cols())
            w *= num_traits<T>::to_double(fk->fan_out_prob(j, cls)) *
                 num_traits<T>::to_double(fk->fan_out_link(j, cls));
        if (sn.nodes[j].nodetype == qn::NodeType::Fork) {
            resolve_fork_dests(sn, P_r, visited, j, cls, w, out, wout, seen);
        } else if (sn.nodes[j].station != 0) {
            out->push_back(j);
            wout->push_back(w);
        }
    }
}

/** The station nodes a fork ultimately feeds, with the expected tasks each receives. */
template <class T>
std::vector<std::size_t> fork_leaves(const qn::NetworkStruct<T>& sn, const Matrix<T>& P_r,
                                     const std::vector<bool>& visited, std::size_t fork_nd,
                                     std::size_t cls = 0, std::vector<double>* weights = 0) {
    std::vector<std::size_t> out;
    std::vector<double> w;
    std::vector<bool> seen(sn.nodes.size(), false);
    resolve_fork_dests(sn, P_r, visited, fork_nd, cls, 1.0, &out, &w, &seen);
    if (weights) *weights = w;
    return out;
}

}  // namespace fjvdetail

/** What the auxiliary construction found, for one class. */
struct FjSpnStructure {
    std::vector<std::size_t> stationNodes, forkNodes, joinNodes;  ///< 0-based node indices
    /**
     * Circulating population: the largest EXPECTED outermost leaf count.
     *
     * On a plain fork every link carries one certain task, so this is the leaf
     * count it has always been and the value is integral. Under a variable
     * forking level it is sum over links of P(branch fires) * E[tasks on it],
     * which is generally FRACTIONAL -- and a fractional token population is not
     * a net anyone can enumerate, which is why the reference solve is skipped
     * exactly there and the closed form is the answer everywhere.
     */
    double B = 1.0;
    std::vector<std::size_t> joinLeaves;  ///< per node, the leaves a Join synchronizes
    /**
     * Per node, whether it sits INSIDE a fork-join region: reachable from an
     * outermost Fork without crossing a Join. These are the stations the net
     * runs at 1/B of the reference station's rate; everything else on the cycle
     * runs at the reference's own rate. A branch is followed to its end, not
     * only to its first station, because the reference's net gives a chained
     * branch station the same rate as the leaf it feeds.
     */
    std::vector<bool> inRegion;
};

/**
 * Classify one class's visited subgraph and size the auxiliary net.
 *
 * Exposed because it is what carries the model content: B and the per-Join leaf
 * counts are the construction, and a test that only checked the ones and zeros
 * of the visit vector would not be testing anything.
 */
template <class T>
FjSpnStructure fj_spn_structure(const qn::NetworkStruct<T>& sn, const Matrix<T>& P_r,
                                const std::vector<bool>& visited, std::size_t cls = 0) {
    const std::size_t I = sn.nodes.size();
    FjSpnStructure out;
    out.joinLeaves.assign(I, 0);
    out.inRegion.assign(I, false);
    for (std::size_t nd = 0; nd < I; ++nd) {
        if (!visited[nd]) continue;
        const qn::NodeType t = sn.nodes[nd].nodetype;
        if (t == qn::NodeType::Fork) out.forkNodes.push_back(nd);
        else if (t == qn::NodeType::Join) out.joinNodes.push_back(nd);
        else if (sn.nodes[nd].station != 0 && t != qn::NodeType::Source &&
                 t != qn::NodeType::Sink)
            out.stationNodes.push_back(nd);
    }

    // Leaf count of each Join, bottom-up. The reference sweeps the Join list
    // once per Join, which is enough to resolve any nesting depth because each
    // pass resolves at least the innermost unresolved level.
    for (std::size_t pass = 0; pass < out.joinNodes.size(); ++pass)
        for (std::size_t ji = 0; ji < out.joinNodes.size(); ++ji) {
            const std::size_t jnd = out.joinNodes[ji];
            std::size_t lc = 0;
            for (std::size_t src = 0; src < I; ++src) {
                if (!(num_traits<T>::to_double(P_r(src, jnd)) > 0.0) || !visited[src]) continue;
                if (sn.nodes[src].nodetype == qn::NodeType::Join && out.joinLeaves[src] > 0)
                    lc += out.joinLeaves[src];
                else
                    lc += 1;  // a station, or a Join not yet resolved
            }
            out.joinLeaves[jnd] = lc;
        }

    // B is the largest leaf count over the OUTERMOST forks. A fork fed by a
    // Join is a serial fork-join stage, not an outer one, and its branches
    // circulate inside the population the outer stage already fixed.
    for (std::size_t fi = 0; fi < out.forkNodes.size(); ++fi) {
        const std::size_t fnd = out.forkNodes[fi];
        bool outermost = true;
        for (std::size_t src = 0; src < I; ++src)
            if (num_traits<T>::to_double(P_r(src, fnd)) > 0.0 && visited[src] &&
                sn.nodes[src].nodetype == qn::NodeType::Join)
                outermost = false;
        if (!outermost) continue;
        // Everything this fork opens, down to its Join, runs at the branch rate.
        std::vector<std::size_t> frontier(1, fnd);
        for (std::size_t h = 0; h < frontier.size(); ++h) {
            const std::size_t nd = frontier[h];
            for (std::size_t j = 0; j < I; ++j) {
                if (!(num_traits<T>::to_double(P_r(nd, j)) > 0.0) || !visited[j]) continue;
                if (sn.nodes[j].nodetype == qn::NodeType::Join || out.inRegion[j]) continue;
                out.inRegion[j] = true;
                frontier.push_back(j);
            }
        }
        std::vector<double> leafw;
        const std::vector<std::size_t> leaves =
            fjvdetail::fork_leaves(sn, P_r, visited, fnd, cls, &leafw);
        if (leaves.empty())
            throw InputError(
                "sn_fj_visits_spn: a Fork reaches no station on any branch, so the auxiliary net "
                "has nothing to synchronize; the routing of this class is malformed");
        // The EXPECTED sibling count, which is the leaf count exactly when every
        // link is certain and carries one task.
        double expected = 0.0;
        for (std::size_t li = 0; li < leaves.size(); ++li) expected += leafw[li];
        if (expected <= 0.0)
            throw InputError(
                "sn_fj_visits_spn: a Fork emits no task in expectation, so its Join can never "
                "fire; at least one branch must be certain to emit at least one task");
        if (expected > out.B) out.B = expected;
    }
    for (std::size_t ji = 0; ji < out.joinNodes.size(); ++ji)
        if (out.joinLeaves[out.joinNodes[ji]] == 0)
            throw InputError(
                "sn_fj_visits_spn: a Join has no branch feeding it, so the auxiliary net cannot "
                "be balanced and its throughputs would not be uniform");
    return out;
}

/**
 * Per-chain fork-join node visit ratios.
 *
 * @return one (nnodes x nclasses) matrix per chain, normalized so the chain's
 *         reference station carries one
 */
template <class T>
std::vector<Matrix<T>> sn_fj_visits_spn(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t I = sn.nodes.size(), K = sn.nclasses, C = sn.nchains;
    std::vector<Matrix<T>> nodevisits(C, Matrix<T>(I, K, zero));
    if (sn.fj.empty()) return nodevisits;  // no fork-join structure at all

    const double fineTol = 1e-8;
    for (std::size_t c = 0; c < C; ++c) {
        if (c >= sn.inchain.size() || sn.inchain[c].empty()) continue;
        const std::vector<std::size_t>& classes = sn.inchain[c];

        for (std::size_t ci = 0; ci < classes.size(); ++ci) {
            const std::size_t r = classes[ci] - 1;  // inchain is 1-based

            // The single-class node routing, read out of the (I*K) block form.
            Matrix<T> P_r(I, I, zero);
            if (sn.rtnodes.rows() == I * K && sn.rtnodes.cols() == I * K)
                for (std::size_t i = 0; i < I; ++i)
                    for (std::size_t j = 0; j < I; ++j) P_r(i, j) = sn.rtnodes(i * K + r, j * K + r);

            // The nodes this class reaches from its reference station.
            const std::size_t refstat = sn.classes[r].refstat;
            if (refstat == 0 || refstat > sn.station_to_node.size()) continue;
            const std::size_t refnode = sn.station_to_node[refstat - 1] - 1;
            std::vector<bool> visited(I, false);
            visited[refnode] = true;
            bool changed = true;
            while (changed) {
                changed = false;
                for (std::size_t i = 0; i < I; ++i) {
                    if (!visited[i]) continue;
                    for (std::size_t j = 0; j < I; ++j)
                        if (num_traits<T>::to_double(P_r(i, j)) > 0.0 && !visited[j]) {
                            visited[j] = true;
                            changed = true;
                        }
                }
            }

            bool hasFork = false;
            for (std::size_t i = 0; i < I; ++i)
                if (visited[i] && sn.nodes[i].nodetype == qn::NodeType::Fork) hasFork = true;

            if (!hasFork) {
                // No fork on this class's path: every visited node is on the
                // ordinary cycle and is visited once.
                for (std::size_t i = 0; i < I; ++i)
                    if (visited[i]) nodevisits[c](i, r) = one;
                continue;
            }

            // The construction is checked, then its rates are used: a station
            // inside a fork-join region fires once per B firings of the cycle,
            // because the pre-fork transition consumes all B tokens and the Join
            // returns them, so it carries 1/B where a station outside carries 1.
            // A Fork and a Join hold no Place at all, hence zero.
            const FjSpnStructure st = fj_spn_structure(sn, P_r, visited, r);
            const T invB = one / num_traits<T>::from_double(st.B);
            for (std::size_t si = 0; si < st.stationNodes.size(); ++si) {
                const std::size_t nd = st.stationNodes[si];
                nodevisits[c](nd, r) = st.inRegion[nd] ? invB : one;
            }
        }

        // Normalize on the chain's reference station, as the reference does.
        const std::size_t r0 = classes[0] - 1;
        const std::size_t refstat0 = sn.classes[r0].refstat;
        if (refstat0 == 0 || refstat0 > sn.station_to_node.size()) continue;
        const std::size_t refnode_c = sn.station_to_node[refstat0 - 1] - 1;
        for (std::size_t ci = 0; ci < classes.size(); ++ci) {
            const std::size_t r = classes[ci] - 1;
            const T nv = nodevisits[c](refnode_c, r);
            if (num_traits<T>::to_double(nv) > fineTol)
                for (std::size_t i = 0; i < I; ++i) nodevisits[c](i, r) = nodevisits[c](i, r) / nv;
        }
    }
    return nodevisits;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_FJ_VISITS_SPN_H
