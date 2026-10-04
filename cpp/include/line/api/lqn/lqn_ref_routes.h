/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LQN_LQN_REF_ROUTES_H
#define LINE_API_LQN_LQN_REF_ROUTES_H

/**
 * @file
 * @ingroup api_lqn
 * Synchronous call DAG carrying reference-task customers into a layer.
 *
 * Port of `matlab/src/api/lqn/lqn_ref_routes.m`. Resolves, for the caller set
 * of one layer, the synchronous call graph along which reference (REF) task
 * customers descend to those callers, and the mean number of times each entry
 * and each call is invoked per REF cycle. SolverLN's `interlock_method='refpath'`
 * reads it to merge the callers that are the same REF customers arriving by
 * different routes into ONE client chain.
 *
 * Nodes are ENTRIES. Visits are computed TOPOLOGICALLY, v(u) = sum over parents
 * of v(p)*w(a)*callmean, and routes are COUNTED in the same pass rather than
 * enumerated. Only SYNC calls are followed: an ASYNC call is send-no-reply, and
 * forwarding is already flattened into pseudo-SYNC arcs by lqn_fwd_rendezvous.
 *
 * INDEXING. Element and call indices are the 1-based ones of LqnStruct. The
 * POSITIONS into `entries` (`LqnRefCall::from`, `::to`, `prefix_pos`) are
 * 0-based, where the reference's are 1-based.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/num/number.h"

namespace line {
namespace lqn {

/** One row of `calls`: the reference's [cidx, fromPos, toPos, aidx, vCall]. */
template <class T>
struct LqnRefCall {
    std::size_t cidx = 0;  ///< call index
    std::size_t from = 0;  ///< 0-based position of the calling entry in `entries`
    std::size_t to = 0;    ///< 0-based position of the called entry in `entries`
    std::size_t aidx = 0;  ///< activity issuing the call
    T vcall{};             ///< mean invocations of the call per REF cycle
};

/** One group, the reference's R(g): the DAG below one REF task. */
template <class T>
struct LqnRefGroup {
    std::size_t reftask = 0;             ///< task index of the REF task at the root
    bool head_is_caller = false;         ///< the REF task is itself a caller of the layer
    std::vector<std::size_t> members;    ///< callers of the layer on the DAG, descent order
    std::vector<std::size_t> entries;    ///< every DAG entry, topological order, roots first
    std::vector<std::size_t> etask;      ///< lqn.parent of each entry
    std::vector<bool> ismember;          ///< that entry's task is a caller of the layer
    std::vector<T> ventry;               ///< mean invocations of each entry per REF cycle
    /// per entry, (aidx, executions per invocation of the entry)
    std::vector<std::vector<std::pair<std::size_t, T>>> actweight;
    std::vector<LqnRefCall<T>> calls;
    std::vector<std::size_t> prefix_pos;  ///< 0-based positions forming the prefix, topological order
    std::vector<bool> prefix_term;        ///< that prefix position is a first caller
    double npaths = 0.0;                  ///< distinct REF-to-caller routes, counted
    /// min of lqn.maxmult over the DAG tasks. A DIAGNOSTIC: the chain population is never capped by it.
    double poolmin = std::numeric_limits<double>::infinity();
};

/** What lqn_ref_routes returns: the groups, or a non-empty `why` and no groups. */
template <class T>
struct LqnRefRoutes {
    std::vector<LqnRefGroup<T>> groups;
    std::string why;  ///< non-empty when the layer must fall back to another interlock method
};

namespace detail {

/** One synchronous successor of an entry: [cidx, called entry, calling activity, callmean]. */
template <class T>
struct LqnRefSucc {
    std::size_t cidx, to, aidx;
    T mean;
};

/** Printable name of an LQN element, nameOf of the reference. */
template <class T>
std::string lqn_ref_name(const LqnStruct<T>& lqn, std::size_t idx) {
    if (idx < lqn.hashnames.size() && !lqn.hashnames[idx].empty()) return lqn.hashnames[idx];
    return "#" + std::to_string(idx);
}

/** MATLAB's %g. */
inline std::string lqn_ref_fmt_g(double x) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%g", x);
    return buf;
}

/** Mean number of invocations carried by call CIDX; 1 when it is not finite. */
template <class T>
T lqn_ref_call_mean(const LqnStruct<T>& lqn, std::size_t cidx) {
    if (cidx < lqn.callproc_mean.size()) {
        const T m = lqn.callproc_mean[cidx];
        if (std::isfinite(num_traits<T>::to_double(m))) return m;
    }
    return num_traits<T>::from_int(1);
}

/**
 * Per calling ENTRY, its synchronous calls. The calling entry is NOT lqn.parent
 * of the activity, which is its TASK: actsof is inverted over the entry range
 * instead, the last entry listing an activity winning, as in the reference.
 */
template <class T>
std::vector<std::vector<LqnRefSucc<T>>> lqn_ref_sync_successors(const LqnStruct<T>& lqn) {
    std::vector<std::vector<LqnRefSucc<T>>> succ(lqn.nidx + 1);
    std::vector<std::size_t> entry_of_act(lqn.nidx + 1, 0);
    for (std::size_t e = 1; e <= lqn.nentries; ++e) {
        const std::size_t eidx = lqn.eshift + e;
        for (std::size_t a : lqn.actsof[eidx])
            if (a <= lqn.nidx) entry_of_act[a] = eidx;
    }
    for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
        if (lqn.calltype[cidx] != lang::CallType::SYNC) continue;
        const std::size_t aidx = lqn.callpair_src[cidx];
        const std::size_t eto = lqn.callpair_dst[cidx];
        if (aidx < 1 || eto < 1 || aidx > lqn.nidx) continue;
        const std::size_t efrom = entry_of_act[aidx];
        if (efrom < 1) continue;
        succ[efrom].push_back({cidx, eto, aidx, lqn_ref_call_mean(lqn, cidx)});
    }
    return succ;
}

/** Every entry reachable from task TIDX over SYNC calls, without pruning. */
template <class T>
std::vector<bool> lqn_ref_reachable(const LqnStruct<T>& lqn,
                                    const std::vector<std::vector<LqnRefSucc<T>>>& succ,
                                    std::size_t tidx) {
    std::vector<bool> seen(lqn.nidx + 1, false);
    std::vector<std::size_t> stack(lqn.entriesof[tidx].begin(), lqn.entriesof[tidx].end());
    while (!stack.empty()) {
        const std::size_t e = stack.back();
        stack.pop_back();
        if (seen[e]) continue;
        seen[e] = true;
        for (const LqnRefSucc<T>& s : succ[e]) stack.push_back(s.to);
    }
    return seen;
}

/**
 * Mean executions of each activity of entry EIDX per invocation of that entry,
 * propagated over lqn.graph so an OR-branch splits by its probabilities. An
 * activity-graph loop is refused rather than truncated.
 */
template <class T>
std::vector<std::pair<std::size_t, T>> lqn_ref_act_weights(const LqnStruct<T>& lqn,
                                                            std::size_t eidx, std::string& why) {
    std::vector<std::pair<std::size_t, T>> aw;
    const std::vector<std::size_t>& acts = lqn.actsof[eidx];
    if (acts.empty()) return aw;
    const T zero = num_traits<T>::from_int(0);
    std::vector<std::size_t> nodeset;
    nodeset.push_back(eidx);
    nodeset.insert(nodeset.end(), acts.begin(), acts.end());
    const std::size_t n = nodeset.size();
    std::vector<std::size_t> pos(lqn.nidx + 1, 0);  // 1-based, 0 = absent
    for (std::size_t i = 0; i < n; ++i) pos[nodeset[i]] = i + 1;
    std::vector<std::vector<T>> A(n, std::vector<T>(n, zero));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t v : lqn.graph.succ(nodeset[i]))
            if (v <= lqn.nidx && pos[v] > 0) A[i][pos[v] - 1] = lqn.graph.get(nodeset[i], v);
    std::vector<long> remaining(n, 0);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i)
            if (A[i][j] > zero) ++remaining[j];
    remaining[0] = 0;  // the entry is the source
    std::vector<T> w(n, zero);
    w[0] = num_traits<T>::from_int(1);
    std::vector<std::size_t> queue;
    for (std::size_t j = 0; j < n; ++j)
        if (remaining[j] == 0) queue.push_back(j);
    std::vector<bool> done(n, false);
    std::size_t ndone = 0, head = 0;
    while (head < queue.size()) {
        const std::size_t i = queue[head++];
        if (done[i]) continue;
        done[i] = true;
        ++ndone;
        for (std::size_t j = 0; j < n; ++j) {
            if (!(A[i][j] > zero)) continue;
            w[j] = T(w[j] + w[i] * A[i][j]);
            --remaining[j];
            if (remaining[j] <= 0 && !done[j]) queue.push_back(j);
        }
    }
    if (ndone < n) {
        why = "the activity graph of entry '" + lqn_ref_name(lqn, eidx) + "' contains a loop";
        return {};
    }
    for (std::size_t i = 1; i < n; ++i) aw.emplace_back(nodeset[i], w[i]);
    return aw;
}

/** One group: the DAG below REF task R, its visits, and the prefix above the first callers. */
template <class T>
bool lqn_ref_build_group(const LqnStruct<T>& lqn,
                         const std::vector<std::vector<LqnRefSucc<T>>>& succ, std::size_t r,
                         const std::vector<bool>& is_caller, double maxpaths,
                         const std::vector<std::size_t>& server_set, LqnRefGroup<T>& g,
                         std::string& why) {
    const T zero = num_traits<T>::from_int(0);
    const std::vector<std::size_t>& roots = lqn.entriesof[r];
    if (roots.empty()) return false;

    // Depth-first sweep with an on-stack marker: a back edge is REFUSED, since a
    // cycle makes v(u) a geometric series the chain would have to express as a
    // self-loop through the layer's own server.
    enum : int { WHITE = 0, GREY = 1, BLACK = 2 };
    std::vector<int> color(lqn.nidx + 1, WHITE);
    std::vector<std::size_t> post;
    for (std::size_t e0 : roots) {
        if (color[e0] != WHITE) continue;
        std::vector<std::pair<std::size_t, std::size_t>> stack{{e0, 0}};
        while (!stack.empty()) {
            const std::size_t u = stack.back().first;
            const std::size_t ci = stack.back().second;
            if (ci == 0) color[u] = GREY;
            if (ci < succ[u].size()) {
                stack.back().second = ci + 1;
                const std::size_t v = succ[u][ci].to;
                if (color[v] == GREY) {
                    why = "the synchronous call graph below '" + lqn_ref_name(lqn, v) +
                          "' is cyclic";
                    return false;
                } else if (color[v] == WHITE) {
                    stack.emplace_back(v, 0);
                }
            } else {
                color[u] = BLACK;
                post.push_back(u);
                stack.pop_back();
            }
        }
    }
    const std::vector<std::size_t> entries(post.rbegin(), post.rend());  // topological order
    if (entries.empty()) return false;
    const std::size_t nE = entries.size();
    std::vector<std::size_t> pos(lqn.nidx + 1, 0);  // 1-based, 0 = absent
    for (std::size_t i = 0; i < nE; ++i) pos[entries[i]] = i + 1;
    std::vector<std::size_t> etask(nE);
    std::vector<bool> ismem(nE);
    bool anymem = false;
    for (std::size_t i = 0; i < nE; ++i) {
        etask[i] = lqn.parent[entries[i]];
        ismem[i] = is_caller[etask[i]];
        anymem = anymem || ismem[i];
    }
    if (!anymem) return false;  // this REF reaches none of the callers

    // Per-entry activity weights, then the call list; calls carry the
    // per-invocation weight w*callmean until the visits below rescale them.
    std::vector<std::vector<std::pair<std::size_t, T>>> actweight(nE);
    std::vector<LqnRefCall<T>> calls;
    for (std::size_t i = 0; i < nE; ++i) {
        const std::size_t u = entries[i];
        std::string awwhy;
        actweight[i] = lqn_ref_act_weights(lqn, u, awwhy);
        if (!awwhy.empty()) {
            why = awwhy;
            return false;
        }
        for (const LqnRefSucc<T>& s : succ[u]) {
            if (pos[s.to] == 0) continue;
            T w = zero;
            for (const std::pair<std::size_t, T>& aw : actweight[i])
                if (aw.first == s.aidx) {
                    w = aw.second;
                    break;
                }
            calls.push_back({s.cidx, i, pos[s.to] - 1, s.aidx, T(w * s.mean)});
        }
    }

    // Topological visits: a root entry is visited 1/nentries times per REF cycle
    std::vector<T> ventry(nE, zero);
    std::vector<std::size_t> rootpos;
    for (std::size_t e : roots)
        if (pos[e] > 0) rootpos.push_back(pos[e] - 1);
    const T rootshare =
        T(num_traits<T>::from_int(1) / num_traits<T>::from_int(int(roots.size())));
    for (std::size_t p : rootpos) ventry[p] = rootshare;
    for (std::size_t i = 0; i < nE; ++i)
        for (const LqnRefCall<T>& c : calls)
            if (c.from == i) ventry[c.to] = T(ventry[c.to] + ventry[i] * c.vcall);
    for (LqnRefCall<T>& c : calls) c.vcall = T(ventry[c.from] * c.vcall);

    // The prefix: entries reachable from a root without passing THROUGH a caller
    std::vector<bool> in_prefix(nE, false), pref_term(nE, false);
    std::vector<double> npath(nE, 0.0);
    for (std::size_t p : rootpos) {
        in_prefix[p] = true;
        npath[p] = 1.0;
    }
    for (std::size_t i = 0; i < nE; ++i) {
        if (!in_prefix[i]) continue;
        if (ismem[i]) {
            pref_term[i] = true;
            continue;  // do not descend past a caller
        }
        // a HOP whose task is a server of this layer would be in two places at once
        for (std::size_t s : server_set)
            if (s == etask[i]) {
                why = "task '" + lqn_ref_name(lqn, etask[i]) +
                      "' is both an intermediate on the reference path and a server of this "
                      "layer";
                return false;
            }
        for (const LqnRefCall<T>& c : calls)
            if (c.from == i) {
                in_prefix[c.to] = true;
                npath[c.to] += npath[i];
            }
    }
    double np = 0.0;
    for (std::size_t i = 0; i < nE; ++i)
        if (pref_term[i]) np += npath[i];
    if (np > maxpaths) {
        why = "the reference path into this layer carries " + lqn_ref_fmt_g(np) +
              " distinct routes, above config.interlock_maxpaths = " + lqn_ref_fmt_g(maxpaths);
        return false;
    }

    g = LqnRefGroup<T>();
    g.reftask = r;
    g.head_is_caller = is_caller[r];
    for (std::size_t i = 0; i < nE; ++i) {
        if (!ismem[i]) continue;
        bool seen = false;
        for (std::size_t t : g.members) seen = seen || t == etask[i];
        if (!seen) g.members.push_back(etask[i]);
    }
    g.entries = entries;
    g.etask = etask;
    g.ismember = ismem;
    g.ventry = ventry;
    g.actweight = actweight;
    g.calls = calls;
    for (std::size_t i = 0; i < nE; ++i)
        if (in_prefix[i]) {
            g.prefix_pos.push_back(i);
            g.prefix_term.push_back(pref_term[i]);
        }
    g.npaths = np;
    // DIAGNOSTIC ONLY: capping the chain at this would delete customers from the REF think stage
    for (std::size_t t : etask)
        if (t > 0 && t < lqn.maxmult.size()) g.poolmin = std::min(g.poolmin, lqn.maxmult[t]);
    return true;
}

}  // namespace detail

/**
 * Resolve the reference routes into the layer whose callers are CALLERS.
 *
 * @param lqn       the layered struct (after lqn_fwd_rendezvous, as SolverLN holds it)
 * @param callers   task indices that call the layer's server
 * @param maxpaths  refuse the layer above this many REF routes into it (default 32)
 * @param server_set server elements of the layer; a prefix node whose task is one of them refuses
 *
 * A caller reachable from two REF tasks is two INDEPENDENT customer pools, and
 * the whole LAYER falls back rather than that caller alone: a refused caller
 * may lie on another group's path, which would count its threads twice.
 */
template <class T>
LqnRefRoutes<T> lqn_ref_routes(const LqnStruct<T>& lqn, const std::vector<std::size_t>& callers,
                               double maxpaths = 32.0,
                               const std::vector<std::size_t>& server_set = {}) {
    LqnRefRoutes<T> out;
    if (lqn.ncalls == 0 || callers.empty()) return out;
    std::vector<bool> is_caller(lqn.nidx + 1, false);
    for (std::size_t c : callers)
        if (c <= lqn.nidx) is_caller[c] = true;
    const std::vector<std::vector<detail::LqnRefSucc<T>>> succ =
        detail::lqn_ref_sync_successors(lqn);

    std::vector<std::size_t> reftasks;
    for (std::size_t t = 1; t <= lqn.ntasks; ++t)
        if (lqn.isref[lqn.tshift + t]) reftasks.push_back(lqn.tshift + t);
    if (reftasks.empty()) return out;

    std::vector<int> nref_of(lqn.nidx + 1, 0);
    for (std::size_t r : reftasks) {
        const std::vector<bool> seen = detail::lqn_ref_reachable(lqn, succ, r);
        std::vector<bool> mem(lqn.nidx + 1, false);
        for (std::size_t e = 1; e <= lqn.nidx; ++e)
            if (seen[e] && is_caller[lqn.parent[e]]) mem[lqn.parent[e]] = true;
        if (is_caller[r]) mem[r] = true;
        for (std::size_t t = 1; t <= lqn.nidx; ++t)
            if (mem[t]) ++nref_of[t];
    }
    for (std::size_t t = 1; t <= lqn.nidx; ++t)
        if (nref_of[t] > 1) {
            out.why = "task '" + detail::lqn_ref_name(lqn, t) + "' is reachable from " +
                      std::to_string(nref_of[t]) +
                      " reference tasks, whose customer pools are independent";
            return out;
        }

    for (std::size_t r : reftasks) {
        LqnRefGroup<T> g;
        std::string gwhy;
        const bool ok =
            detail::lqn_ref_build_group(lqn, succ, r, is_caller, maxpaths, server_set, g, gwhy);
        if (!gwhy.empty()) {
            out.why = gwhy;
            out.groups.clear();
            return out;
        }
        if (ok) out.groups.push_back(std::move(g));
    }
    return out;
}

}  // namespace lqn
}  // namespace line

#endif  // LINE_API_LQN_LQN_REF_ROUTES_H
