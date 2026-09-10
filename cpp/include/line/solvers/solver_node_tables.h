/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SOLVER_NODE_TABLES_H
#define LINE_SOLVERS_SOLVER_NODE_TABLES_H

/**
 * The NODE-indexed view of a station result, behind `getAvgNodeTable`.
 *
 * ONE DEFINITION, THREE CALLERS. `line-cli -a node` and `-a nodechain` were the
 * first two, and the comment on `node_metrics` already said why they must not
 * each build the scatter: "recomputing them in the second arm would be a second
 * definition of the same scatter, free to disagree with the first". The third
 * is the example corpus. Four `cache_replc_*` references call
 * `getAvgNodeTable()` rather than `getAvgTable()` -- a cache model's answer
 * lives at the Cache and the ClassSwitch, neither of which is a station -- and
 * their C++ twins printed the station table instead, so every Cache and
 * ClassSwitch row was absent and the goldens' 71 cells went unmeasured. This
 * header is what lets a twin print the table its reference printed without a
 * fourth copy of the rule.
 *
 * ARITHMETIC: field. Sums, quotients and the two `sn_*` ports it delegates to.
 */

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "line/api/sn/sn_node_metrics.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/matrix.h"

namespace line {
namespace solvers {

/**
 * The station AvgResult of a solver whose runner returns its own solution type,
 * i.e. SSA and Fluid: their QN/UN/RN/TN are the same six columns the AvgTable
 * arm already prints, so a VIEW of them is a view of the same numbers.
 *
 * The two columns those runners do not carry are filled the way the reference
 * fills them, and the SAME way their own `-a avg` arms do:
 * - ResidT is `sn_get_residt_from_respt`, the per-JOB residence time. It is a
 *   pure function of `sn` and RN, which is why SolverMVA's runner calls it
 *   rather than deriving it inside the analyzer, and why a solver that reports
 *   no residence time of its own is not thereby excused from reporting one.
 *   ResidT = RespT holds only when every station is visited once per cycle;
 *   assuming it cost `sdroute_closed` a factor of 3 on the two Queues and
 *   `init_state_ps` a factor of 17 on Queue1.
 * - ArvR is the throughput except at a Source, which has no arrivals TO ITSELF.
 * Restating either rule here rather than sharing it would let `-a avg` and
 * `-a node` disagree on one column of one model, which is the divergence
 * `-s mva` vs `-s fluid` already produced once on gallery_mm1.
 */
template <class T>
line::mva::AvgResult<T> avg_result_from_sim(const line::qn::NetworkStruct<T>& sn,
                                            const line::Matrix<double>& QN,
                                            const line::Matrix<double>& UN,
                                            const line::Matrix<double>& RN,
                                            const line::Matrix<double>& TN,
                                            const std::vector<double>& CN,
                                            const std::vector<double>& XN,
                                            const std::string& method) {
    const std::size_t M = sn.nstations, R = sn.nclasses;
    const T zero = line::num_traits<T>::from_int(0);
    line::mva::AvgResult<T> r;
    r.QN = line::Matrix<T>(M, R, zero);
    r.UN = line::Matrix<T>(M, R, zero);
    r.RN = line::Matrix<T>(M, R, zero);
    r.WN = line::Matrix<T>(M, R, zero);
    r.AN = line::Matrix<T>(M, R, zero);
    r.TN = line::Matrix<T>(M, R, zero);
    std::vector<std::vector<bool>> srcmask(M, std::vector<bool>(R, false));
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].sched == line::lang::SchedStrategy::EXT)
            for (std::size_t c = 0; c < R; ++c) srcmask[i][c] = true;
        for (std::size_t c = 0; c < R; ++c) {
            r.QN(i, c) = line::num_traits<T>::from_double(QN(i, c));
            r.UN(i, c) = line::num_traits<T>::from_double(UN(i, c));
            r.RN(i, c) = line::num_traits<T>::from_double(RN(i, c));
            r.TN(i, c) = line::num_traits<T>::from_double(TN(i, c));
        }
    }
    r.WN = line::mva::sn_get_residt_from_respt<T>(sn, r.RN);
    // ARRIVAL RATE IS AN INFLOW, NOT A THROUGHPUT. `AN = TN` holds only where
    // the station's own flow balances, and a fluid fixed point need not: on
    // cache_replc_routing MATLAB's rmf reports Tput 0 at Delay1 and Delay2 while
    // 0.4 and 0.6 arrive there, so the shortcut dropped those four rows out of
    // the node table entirely. `SolverFLD/runAnalyzer.m:372` and
    // `SolverSSA/runAnalyzer.m:199` both derive AN through sn_get_arvr_from_tput
    // from the class-expanded routing, exactly as CTMC, NC, MAM, BA and QNS
    // already do in this port; the Source mask is `getAvg`'s own zeroMask.
    r.AN = line::mva::filter_metric(sn, line::mva::sn_get_arvr_from_tput(sn, r.TN),
                                    line::mva::MetricKind::ArvR, &srcmask);
    for (std::size_t c = 0; c < CN.size(); ++c)
        r.CN.push_back(line::num_traits<T>::from_double(CN[c]));
    for (std::size_t c = 0; c < XN.size(); ++c)
        r.XN.push_back(line::num_traits<T>::from_double(XN[c]));
    r.method = method;
    r.actualmethod = method;
    return r;
}

/**
 * The station table scattered to the NODE index space, plus the two flow
 * columns the reference recomputes there.
 *
 * Factored out of `-a node` because `-a nodechain` aggregates exactly these six
 * matrices by chain: recomputing them in the second arm would be a second
 * definition of the same scatter, free to disagree with the first.
 */
template <class T>
struct NodeMetrics {
    line::Matrix<T> QN, UN, RN, WN, AN, TN;  ///< (nnodes x nclasses)
};

template <class T>
NodeMetrics<T> node_metrics(const line::qn::NetworkStruct<T>& sn,
                            const line::mva::AvgResult<T>& r) {
    const std::size_t I = sn.nodes.size(), R = sn.nclasses;
    const T zero = line::num_traits<T>::from_int(0);
    NodeMetrics<T> m;
    m.QN = line::Matrix<T>(I, R, zero);
    m.UN = line::Matrix<T>(I, R, zero);
    m.RN = line::Matrix<T>(I, R, zero);
    m.WN = line::Matrix<T>(I, R, zero);
    for (std::size_t ist = 0; ist < sn.nstations; ++ist) {
        const std::size_t ind = sn.station_to_node[ist];
        if (!ind) continue;
        for (std::size_t c = 0; c < R; ++c) {
            m.QN(ind - 1, c) = r.QN(ist, c);
            m.UN(ind - 1, c) = r.UN(ist, c);
            m.RN(ind - 1, c) = r.RN(ist, c);
            m.WN(ind - 1, c) = r.WN(ist, c);
        }
    }
    // THE TWO RECOMPUTED COLUMNS READ THE REFRESHED STRUCT, exactly as the
    // station table already does (`solver_mva_runner.h`, `refL`). A cache's
    // hit/miss split is a SOLVER RESULT: `link()` leaves the self-switch at the
    // offered 1/2-1/2, and both of these columns are visit-weighted, so reading
    // the base struct here reports the split nobody computed. Measured on
    // cache_replc_routing, Cache/Router/Sink hit and miss came out 1/1 against a
    // golden of 0.8/1.2, and on cache_replc_fifo the Cache read-class ArvR came
    // out 2 against 1 -- the two halves of one guess, summed.
    //
    // NULL for every model with no cache branch, which is the base struct and
    // the behaviour every other model already had.
    const line::qn::NetworkStruct<T>& refsn = r.refreshed_struct ? *r.refreshed_struct : sn;
    // AND THE SPLIT ITSELF IS PASSED, not left to be inferred from the visits.
    // A refreshed struct is only ONE of the two ways a solver states what it
    // measured, and the arms that do not build one -- CTMC, SSA, and the
    // Source-Cache-Sink branch of MVA and NC -- still report a hit probability
    // in `AvgResult::cache`. Without it every one of them fell through to the
    // visit-ratio branch of sn_get_node_tput_from_tput, i.e. to the 1/2-1/2
    // `link()` offers before anything is solved: 0.5/0.5 against a golden of
    // 0.4/0.6 on cache_replc_fifo and cache_replc_lru, 1/1 against 1.1127/
    // 0.8873 on cache_replc_rr. Empty on every model with no Cache node.
    std::map<std::size_t, line::api::CacheActualProb<T>> cacheprob;
    for (std::size_t ci = 0; ci < r.cache.caches.size(); ++ci) {
        const line::solvers::CacheNodeMetrics<T>& cm = r.cache.caches[ci];
        if (cm.node == 0 || cm.hitprob.empty()) continue;
        line::api::CacheActualProb<T> p;
        p.hit = cm.hitprob;
        p.miss = cm.missprob;
        p.delayed_hit = cm.delayedprob;
        cacheprob[cm.node] = p;
    }
    m.AN = line::api::sn_get_node_arvr_from_tput(refsn, r.TN, r.AN);
    m.TN = line::api::sn_get_node_tput_from_tput(refsn, r.TN, m.AN, cacheprob);
    // WHAT A CACHE SENDS ON IS WHAT THE NEXT NODE RECEIVES, and the visit ratios
    // cannot say so: `sn.nodevisits` still carries the 1/2-1/2 `link()` offered,
    // so a ClassSwitch or a Sink downstream of a Cache reported that guess as its
    // hit and miss ARRIVAL rate even once the departure rate above it was right.
    // `getAvgNode.m` lines 73-92 close the same gap with the same assignment, and
    // AFTER both tables are formed rather than before: the fix is to ANn alone
    // and must not feed back into the TNn that produced it.
    for (std::size_t cind = 1; cind <= I; ++cind) {
        if (refsn.nodes[cind - 1].nodetype != line::qn::NodeType::Cache) continue;
        typename std::map<std::size_t, line::qn::CacheParam<T>>::const_iterator np =
            refsn.nodeparam.find(cind);
        if (np == refsn.nodeparam.end()) continue;
        for (std::size_t ind = 1; ind <= I; ++ind) {
            const line::qn::NodeType nt = refsn.nodes[ind - 1].nodetype;
            if (nt != line::qn::NodeType::ClassSwitch && nt != line::qn::NodeType::Sink) continue;
            for (std::size_t rr = 1; rr <= R; ++rr)
                if (line::api::detail::sn_is_cache_class(np->second.hitclass,
                                                         np->second.missclass, rr))
                    m.AN(ind - 1, rr - 1) = m.TN(cind - 1, rr - 1);
        }
    }
    return m;
}
}  // namespace solvers
}  // namespace line

#endif  // LINE_SOLVERS_SOLVER_NODE_TABLES_H
