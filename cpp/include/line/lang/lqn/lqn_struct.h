/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_LQN_LQN_STRUCT_H
#define LINE_LANG_LQN_LQN_STRUCT_H

/**
 * LayeredNetworkStruct, the flattened description of a layered queueing network.
 *
 * Port of matlab/src/lang/layered/LayeredNetworkStruct.m and of the fields that
 * matlab/src/lang/layered/@@LayeredNetwork/getStruct.m populates. Only the
 * fields SolverLN and SolverMVA read are carried; the process-descriptor
 * families (hostdem_proc, itemproc, setuptime, delayofftime) exist in MATLAB to
 * serve solvers this port does not have and are omitted rather than filled with
 * placeholders.
 *
 * INDEXING. Element indices are 1-based and live in one flat space shared by
 * the four element kinds, exactly as in MATLAB:
 *
 *   hosts       1                 .. nhosts          (hshift = 0)
 *   tasks       tshift+1          .. tshift+ntasks   (tshift = nhosts)
 *   entries     eshift+1          .. eshift+nentries (eshift = nhosts+ntasks)
 *   activities  ashift+1          .. ashift+nacts    (ashift = eshift+nentries)
 *
 * with nidx = ashift + nacts and cshift = nidx, so a call index cidx is
 * addressed as nidx+cidx inside the entry-service matrix. Calls are numbered
 * separately, 1..ncalls. Vectors are sized nidx+1 (or ncalls+1) and slot 0 is
 * unused; this keeps every index expression identical to the reference, which
 * is worth more here than the one wasted slot, because the arithmetic on these
 * indices (parent-of-parent to reach a host, aidx-ashift to reach a phase) is
 * dense and a systematic off-by-one would be silent.
 *
 * ARITHMETIC. Means, think times and call multiplicities are T. Element
 * multiplicities, replication counts and job populations are `double`, matching
 * MATLAB: they are counts, they can be infinite (an inf-scheduled task), and
 * they enter the solvers as integer populations, so nothing is gained by
 * carrying them in the exact type and the infinity would have to be emulated.
 */

#include <cmath>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lqn {

using lang::CallType;
using lang::CdScaling;
using lang::Distrib;
using lang::LqnElement;
using lang::PrecedenceType;
using lang::ReplacementStrategy;
using lang::SchedStrategy;

/**
 * Heterogeneous server pools declared on a layer server, the twin of the
 * `nservertypes` / `servertypenames` / `serverspertype` / `servercompat` /
 * `heterorates` block a Network carries in `sn.nodeparam{i}`.
 *
 * `compat(t, j)` is nonzero when pool t may serve OPERAND j -- task j of a host,
 * entry j of a task, in declaration order, the same operand space the
 * class- and joint-dependence handles read. SolverLN lowers the pools to the
 * activated-server rate mu(n) = sum_t counts(t)*rates(t)*[pool t compatible with
 * some operand present], which is order-independent, and hands it to the layer
 * station as a joint dependence; see _kb/06-solver-catalog.md (LN section).
 */
template <class T>
struct ServerPools {
    std::vector<std::string> names;  ///< (npools) declared pool name
    std::vector<double> counts;      ///< (npools) servers held by each pool
    std::vector<T> rates;            ///< (npools) per-pool rate multiplier
    Matrix<T> compat;                ///< (npools x noperands), nonzero = eligible

    bool empty() const { return names.empty(); }
    std::size_t npools() const { return names.size(); }
};

/**
 * A sparse square matrix over element indices, held as a dense vector of rows
 * with an explicit nonzero list per row.
 *
 * The LQN graph is very sparse (each activity has one or two successors) and
 * the algorithms walk it by "successors of i", never by column, except for
 * `find(graph(:,j))` in the reply search and `find(taskgraph(:,t))` in the
 * multiplicity correction. Both column queries are rare, so they scan.
 */
template <class T>
struct SparseGraph {
    std::size_t n = 0;
    std::vector<std::vector<std::pair<std::size_t, T>>> row;  ///< 1-based, row[0] unused

    void resize(std::size_t nn) {
        n = nn;
        row.assign(nn + 1, {});
    }
    void set(std::size_t i, std::size_t j, const T& v) {
        for (auto& e : row[i])
            if (e.first == j) {
                e.second = v;
                return;
            }
        row[i].emplace_back(j, v);
    }
    T get(std::size_t i, std::size_t j) const {
        for (const auto& e : row[i])
            if (e.first == j) return e.second;
        return num_traits<T>::from_int(0);
    }
    /** Successors of i in ascending index order, as MATLAB's find() returns them. */
    std::vector<std::size_t> succ(std::size_t i) const {
        std::vector<std::size_t> out;
        const T zero = num_traits<T>::from_int(0);
        for (const auto& e : row[i])
            if (e.second != zero) out.push_back(e.first);
        std::sort(out.begin(), out.end());
        return out;
    }
    /** Predecessors of j in ascending index order. */
    std::vector<std::size_t> pred(std::size_t j) const {
        std::vector<std::size_t> out;
        const T zero = num_traits<T>::from_int(0);
        for (std::size_t i = 1; i <= n; ++i)
            for (const auto& e : row[i])
                if (e.first == j && e.second != zero) {
                    out.push_back(i);
                    break;
                }
        return out;
    }
    void erase(std::size_t i, std::size_t j) {
        for (std::size_t k = 0; k < row[i].size(); ++k)
            if (row[i][k].first == j) {
                row[i].erase(row[i].begin() + k);
                return;
            }
    }
};

/** A boolean sparse relation over element indices, e.g. iscaller. */
struct BoolGraph {
    std::size_t n = 0;
    std::vector<std::vector<std::size_t>> row;

    void resize(std::size_t nn) {
        n = nn;
        row.assign(nn + 1, {});
    }
    void set(std::size_t i, std::size_t j) {
        for (std::size_t v : row[i])
            if (v == j) return;
        row[i].push_back(j);
    }
    bool get(std::size_t i, std::size_t j) const {
        for (std::size_t v : row[i])
            if (v == j) return true;
        return false;
    }
    bool any_row(std::size_t i) const { return !row[i].empty(); }
    bool any_col(std::size_t j) const {
        for (std::size_t i = 1; i <= n; ++i)
            if (get(i, j)) return true;
        return false;
    }
    std::vector<std::size_t> col(std::size_t j) const {
        std::vector<std::size_t> out;
        for (std::size_t i = 1; i <= n; ++i)
            if (get(i, j)) out.push_back(i);
        return out;
    }
};

/** One activity precedence of a task, with its activities resolved to indices. */
template <class T>
struct LqnPrecedence {
    PrecedenceType pretype = PrecedenceType::PRE_SEQ;
    PrecedenceType posttype = PrecedenceType::POST_SEQ;
    std::vector<std::size_t> preacts;   ///< absolute activity indices
    std::vector<std::size_t> postacts;  ///< absolute activity indices
    std::vector<T> preparams;           ///< PRE_OR shares, or a PRE_AND quorum
    std::vector<T> postparams;          ///< POST_OR probabilities or the POST_LOOP count
};

/**
 * One routed call group: an activity, the strategy that picks among its
 * targets, and the target ENTRIES in declaration order.
 *
 * Not a template: it holds no numeric parameter, the per-target call means
 * living on the member calls in `callproc` as usual.
 */
struct LqnCallGroup {
    std::size_t caller = 0;  ///< absolute index of the dispatching activity
    lang::RoutingStrategy strategy = lang::RoutingStrategy::RROBIN;
    std::vector<std::size_t> targets;  ///< absolute entry indices, in declaration order
};

template <class T>
struct LqnStruct {
    std::size_t nidx = 0, nhosts = 0, ntasks = 0, nentries = 0, nacts = 0, ncalls = 0;
    std::size_t hshift = 0, tshift = 0, eshift = 0, ashift = 0, cshift = 0;

    std::vector<std::string> names;      ///< (nidx+1) declared name
    std::vector<std::string> hashnames;  ///< (nidx+1) name prefixed by kind: P:/T:/R:/E:/A:
    std::vector<LqnElement> type;        ///< (nidx+1)
    std::vector<std::size_t> parent;     ///< (nidx+1) host of a task, task of an entry/activity
    std::vector<SchedStrategy> sched;    ///< (tshift+ntasks+1)
    std::vector<double> mult;            ///< (tshift+ntasks+1) declared multiplicity, may be Inf
    std::vector<double> maxmult;         ///< (tshift+ntasks+1) sustainable multiplicity
    std::vector<double> repl;            ///< (tshift+ntasks+1) replication

    /**
     * Queue-dependent service rates declared on a layer server (a host or a
     * task), by element index, empty where absent.
     *
     * `lldscaling[i]` is the vector alpha(n) applied at total population n;
     * `cdscaling[i]` and `jdscaling[i]` are the per-OPERAND handles beta(n) and
     * eta(n), whose argument counts the jobs the layer station holds on behalf
     * of task j of a host or entry j of a task, in declaration order. The peak
     * vectors are required beside the handles, since utilization at such a
     * station is reported as U = T*S/peak. `pools` carries a compatibility
     * declaration, which SolverLN lowers to a jdscaling of its own.
     *
     * Only the class-switching layer builders emit these; the composed
     * phase-type law replaces the station by an entry law and so refuses them
     * by name -- see _kb/04-networkstruct.md and _kb/06-solver-catalog.md.
     */
    std::vector<std::vector<T>> lldscaling;     ///< (tshift+ntasks+1)
    std::vector<CdScaling<T>> cdscaling;        ///< (tshift+ntasks+1)
    std::vector<std::vector<T>> cdscalingpeak;  ///< (tshift+ntasks+1)
    std::vector<CdScaling<T>> jdscaling;        ///< (tshift+ntasks+1)
    std::vector<std::vector<T>> jdscalingpeak;  ///< (tshift+ntasks+1)
    std::vector<ServerPools<T>> pools;          ///< (tshift+ntasks+1)

    /**
     * Fan-out and fan-in, keyed by task element index, absent = 0.
     *
     * `fanout[{caller, callee}]` is how many callee replicas one caller replica
     * addresses, and is read by SolverLN to decide whether a replicated task
     * layer can be pooled into a single station instead of materialised once
     * per replica (`getStruct.m` builds the same matrix as `lsn.fanout`).
     * `fanin` is the mirror declaration and is carried for round-trip fidelity
     * only: no analyzer in ANY codebase reads it, MATLAB and the JAR likewise
     * park it on the Task and never consult it.
     */
    std::map<std::pair<std::size_t, std::size_t>, double> fanout;
    std::map<std::pair<std::size_t, std::size_t>, double> fanin;

    /** fan-out from caller task `i` to callee task `j`; 0 when undeclared. */
    double fanout_at(std::size_t i, std::size_t j) const {
        const std::map<std::pair<std::size_t, std::size_t>, double>::const_iterator it =
            fanout.find(std::make_pair(i, j));
        return it == fanout.end() ? 0.0 : it->second;
    }

    std::vector<bool> isref;             ///< (tshift+ntasks+1)
    std::vector<bool> iscache;           ///< (tshift+ntasks+1)
    std::vector<bool> hassetup;        ///< (tshift+ntasks+1)

    /**
     * Cache tasks and item entries.
     *
     * `nitems` is indexed by ELEMENT and carries the item population on BOTH
     * the cache task and each of its item entries, which is how getStruct.m
     * writes it (`iscache` is the nitems>0 test, over hosts+tasks only).
     * `itemcap` is the capacity of each cache list, so a plain single-level
     * cache has one entry. `itemproc` is the item POPULARITY of an item entry,
     * as an explicit pmf over its `nitems`: the reference stores a discrete
     * Distribution (a Zipf, typically) and only ever reads its pmf, and this
     * port has no discrete-distribution type to put there.
     */
    std::vector<std::size_t> nitems;                    ///< (nidx+1)
    std::vector<std::vector<int>> itemcap;              ///< (tshift+ntasks+1)
    std::vector<ReplacementStrategy> replacestrat;      ///< (tshift+ntasks+1)
    std::vector<std::vector<T>> itemproc;               ///< (nidx+1) popularity pmf

    /**
     * Setup tasks: the server powers down when idle and pays to restart.
     *
     * (tshift+ntasks+1); `hassetup` is the "a setup time is declared and is
     * neither Immediate nor sub-tolerance" test, matching how the reference
     * gates the feature in LQN2QN's functionTimesOf rather than the bare
     * ~isempty of getStruct.
     */
    std::vector<Distrib<T>> setuptime;
    std::vector<Distrib<T>> delayofftime;

    std::vector<Distrib<T>> hostdem;   ///< (nidx+1) host demand per activity (Immediate elsewhere)
    std::vector<Distrib<T>> think;     ///< (nidx+1) task think time
    std::vector<Distrib<T>> actthink;  ///< (nidx+1) activity think time
    std::vector<bool> has_arrival;     ///< (nidx+1) entry with an open arrival
    std::vector<Distrib<T>> arrival;   ///< (nidx+1) open arrival process of an entry

    std::vector<std::vector<std::size_t>> tasksof;    ///< (nhosts+1)
    std::vector<std::vector<std::size_t>> entriesof;  ///< (tshift+ntasks+1)
    std::vector<std::vector<std::size_t>> actsof;     ///< (ashift+1) by task and by entry
    std::vector<std::vector<std::size_t>> callsof;    ///< (nidx+1) call indices issued by an activity

    std::vector<std::size_t> callpair_src;  ///< (ncalls+1) calling activity (entry for FWD)
    std::vector<std::size_t> callpair_dst;  ///< (ncalls+1) called entry
    std::vector<CallType> calltype;         ///< (ncalls+1)
    std::vector<T> callproc_mean;           ///< (ncalls+1) mean number of calls
    std::vector<std::string> callnames;     ///< (ncalls+1)
    std::vector<std::string> callhashnames; ///< (ncalls+1)

    SparseGraph<T> graph;      ///< element call/precedence graph, edge weights are branch shares
    SparseGraph<T> dag;        ///< graph with entry-task edges reversed and loop back-edges removed
    SparseGraph<T> taskgraph;  ///< task-to-task calls
    BoolGraph iscaller, issynccaller, isasynccaller;

    /**
     * Admission constraint `A n <= b` on the layer station of a host or task.
     *
     * (tshift+ntasks+1); an empty A means unconstrained. The COLUMNS are that
     * host's tasks (`tasksof`) or that task's entries (`entriesof`), in that
     * order -- element space, not class space. SolverLN::build_layer expands
     * them into the layer's own classes and emits a Region on the server, an
     * entry column becoming the CALL classes that target it and a task column
     * the ACTIVITY classes of that task. Carried by the JSON interchange as
     * `admissionConstraints`, NOT by .lqnx. See _kb/04-networkstruct.md.
     */
    std::vector<Matrix<T>> lincon_A;
    std::vector<std::vector<T>> lincon_b;

    /**
     * Activity precedences of each task, as DECLARED, indexed by the task's
     * absolute index.
     *
     * `graph` is the same information after the reader has expanded it into
     * arcs, and that expansion is lossy for a loop: the back edge carries the
     * branch PROBABILITY 1 - 1/count, so recovering the count from it inverts a
     * division and loses the distinction between a loop and an ordinary cycle.
     * SolverLN method 'srvn.ph' composes the activity graph into a phase-type law
     * instead of routing it, and needs the count, so it reads this. Every other
     * consumer reads `graph`.
     */
    std::vector<std::vector<LqnPrecedence<T>>> precedences;  ///< (tshift+ntasks+1)

    /**
     * Synchronous calls DISPATCHED AS A GROUP, `lsn.callgroups`.
     *
     * `synchCallRoundRobin` / `synchCallJSQ` issue one call per invocation whose
     * destination cycles over, or is chosen among, several target entries. The
     * members are ordinary SYNC calls of mean `total/n` and are already in
     * `callpair`; what this adds is that they are ONE dispatch decision rather
     * than n independent Bernoulli draws, which is the whole point -- the mean
     * call rate is the same and the variance is not.
     *
     * Representable only under the squashed layering, since the targets must
     * share a submodel for a dispatch among them to mean anything, and only
     * under a layer solver that resolves the strategy from the state.
     */
    std::vector<LqnCallGroup> callgroups;

    std::vector<PrecedenceType> actpretype;   ///< (nidx+1)
    std::vector<PrecedenceType> actposttype;  ///< (nidx+1)
    std::vector<std::size_t> actquorum;       ///< (nidx+1) AND-join quorum, on the join target
    std::vector<int> actphase;                ///< (nacts+1) phase of each activity, 1-based by act

    /** Index of the host of the element, for a task or anything owned by one. */
    std::size_t host_of(std::size_t idx) const {
        if (type[idx] == LqnElement::HOST) return idx;
        if (type[idx] == LqnElement::TASK) return parent[idx];
        return parent[parent[idx]];
    }
};

}  // namespace lqn
}  // namespace line

#endif  // LINE_LANG_LQN_LQN_STRUCT_H
