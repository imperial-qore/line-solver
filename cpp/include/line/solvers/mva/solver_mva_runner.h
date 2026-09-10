/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_RUNNER_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_RUNNER_H

/**
 * The SolverMVA class surface: `@@SolverMVA/runAnalyzer.m` and the gates around
 * it.
 *
 * What sits here rather than in the dispatch is everything that happens BEFORE
 * and AFTER one inner solve: the method whitelist, the structural gates, the
 * conversions from response time to residence time and from throughput to
 * arrival rate, and the metric filter. `mvaDispatch` is the callback; on a
 * model without forks it runs exactly once.
 *
 * THE FILTER IS LOAD BEARING, not cosmetic. `@@NetworkSolver/getAvg` is not a
 * getter: between the analyzer and the caller it zeroes a metric wherever the
 * class has no visit, snaps anything below FineTol to zero, and masks the
 * queue length and the utilization where the response time is below
 * 10 FineTol. In a layered model every entry, task and call class is served by
 * an Immediate distribution somewhere, whose response time is exactly 1e-8, so
 * without the mask a layer reports the immediate classes' share of the
 * population as real queue length.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/mva/fj_driver.h"
#include "line/solvers/mva/fj_ht.h"
#include "line/solvers/mva/fj_mmt.h"
#include "line/solvers/mva/solver_mva_sjn.h"
#include "line/api/sn/sn_has_bursty_arrival.h"
#include "line/api/sn/sn_patience_handles.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/mva/mva_dispatch.h"

namespace line {
namespace mva {

/** The metrics `getAvg` returns, after filtering. */
template <class T>
struct AvgResult {
    Matrix<T> QN;   ///< queue length
    Matrix<T> UN;   ///< utilization
    Matrix<T> RN;   ///< response time, per visit
    Matrix<T> TN;   ///< throughput
    Matrix<T> AN;   ///< arrival rate
    Matrix<T> WN;   ///< residence time, per job
    std::vector<T> CN;  ///< system response time per class
    std::vector<T> XN;  ///< system throughput per class
    std::string method;      ///< the method asked for
    std::string actualmethod;///< the algorithm that ran
    /**
     * The reference's own warning text, verbatim, empty when it did not warn.
     *
     * It belongs beside `actualmethod` because it is the same kind of thing: a
     * statement about HOW the numbers were produced rather than a metric. It is
     * carried, not thrown, because the answer is usable -- the SJN starvation
     * cap leaves the population law exact and only the ACCURACY unwarranted,
     * which is solver_nc_cdf.h:83's precedent, not solver_mam_retrial.h:437's.
     * Dropping it silently is worse than a wrong number: the answer then looks
     * authoritative exactly where the reference declines to stand behind it.
     */
    std::string warning;
    /**
     * (h) mean storage cost held by each cache list, K_j = sum_i sigma_i pi_ij,
     * filled only by the NC cache branch on a model with item sizes; EMPTY
     * otherwise. See ton21cache Sec. IX and [[09-ldes-and-cache]].
     */
    std::vector<T> listcost;
    int iter = 0;
    /**
     * Whether the fixed point met its tolerance, empty when the handler reports
     * none. Carried beside `iter` because the count alone cannot decide it on
     * the AMVA route, where it aggregates the nested inner sweeps and reaches
     * the cap on solves whose outer residual is exactly zero. A caller that
     * reads the count as a convergence test raises false warnings; that is why
     * the flag exists and why it has to reach the CLI envelope.
     */
    std::optional<bool> converged;
    /**
     * `@@SolverNC/getProbNormConstAggr`, i.e. the reference's
     * `result.Prob.logNormConstAggr`, filled by the NC runner alone and EMPTY
     * for every other solver.
     *
     * It sits here for `listcost`'s reason: it is a statement about the solve
     * that only one solver can make, and the alternative -- calling
     * `solver_nc_solve` a second time to recover it -- would solve the model
     * twice for one number. Absent is not zero: log G = 0 is the constant of an
     * empty network, so a caller must test the option rather than the value.
     */
    std::optional<double> lognormconst;
    /**
     * What the cache branches observed, EMPTY on a model with no Cache node and
     * on every solver that does not analyze one. `getAvgCacheTable` and
     * `getAvgItemTable` are built from it; see `solvers::CacheMetrics`.
     */
    solvers::CacheMetrics<T> cache;
    /**
     * The struct whose cache self-switch carries the CONVERGED hit/miss split,
     * filled by the cacheqn branches and NULL for every other model.
     *
     * It has to travel with the result because the split is a SOLVER OUTPUT, not
     * model structure: `link()` leaves the self-switch at the offered 1/2-1/2,
     * and every quantity derived from a visit ratio -- ArvR, ResidT, and the
     * node-level hit and miss throughputs -- is wrong by that ratio until the
     * solve replaces it. The runner already uses this struct for the STATION
     * table (see `refL` below); handing the caller the base struct instead makes
     * `-a avg` and `-a node` disagree on one column of one model, which is the
     * divergence `node_metrics` exists to avoid. Measured on cache_replc_routing:
     * Cache/Router/Sink hit and miss read 1/1 off the base struct against 0.8/1.2.
     */
    std::shared_ptr<qn::NetworkStruct<T>> refreshed_struct;
};

/**
 * Port of `SolverMVA.listValidMethods`.
 *
 * A LISTED NAME MUST ACTUALLY RUN. The reference records why: the bound family
 * was delisted when it moved to SolverBA, because listing names that error made
 * the list a false claim. `qdaql` is implemented in `solver_amva` and dispatched
 * by no analyzer, so it stays off the list here too and the gate below refuses
 * it.
 */
template <class T>
std::vector<std::string> list_valid_methods(const qn::NetworkStruct<T>& L) {
    std::vector<std::string> m{"default", "mva",   "exact",       "amva",
                               "sum",     "esum",  "qdlin",       "amva.qdlin",
                               "bs",      "amva.bs", "sqni",      "qd",
                               "amva.qd", "qli",   "amva.qli",    "fli",
                               "amva.fli", "ab",   "amva.ab",     "schmidt",
                               "amva.schmidt", "schmidt-ext", "amva.schmidt-ext",
                               "scat", "amva.scat",
                               "lcp", "amva.lcp", "chow", "amva.chow",
                               "pamb", "amva.pamb", "pami", "amva.pami",
                               "pamt", "amva.pamt", "clust", "amva.clust",
                               "dmlin", "amva.dmlin"};

    // AQL (pfqn_aql), QSA (pfqn_qsa) and Tay (pfqn_tay) reject multiserver
    // stations in solver_amva, exactly as MATLAB SolverMVA.listValidMethods
    // does, so all three are advertised for single-server models only.
    if (!L.has_multi_server()) {
        m.push_back("aql");
        m.push_back("amva.aql");
        m.push_back("qsa");
        m.push_back("amva.qsa");
        m.push_back("tay");
        m.push_back("amva.tay");
    }

    bool open = true, closed = false;
    for (const qn::JobClass& c : L.classes) {
        if (std::isfinite(c.population) && c.population > 0.0) {
            closed = true;
            open = false;
        }
    }
    (void)closed;
    // SQD is the Smith queue-decomposition handler for a closed single-class
    // Blocking-After-Service model; mva_is_bas_model advertises it only there --
    // a closed single class with a BAS drop rule -- not for every single-chain
    // closed model, which would let check_method accept 'sqd' where MATLAB
    // rejects it by name.
    bool has_bas = false;
    for (const qn::Station<T>& st : L.stations)
        for (int dr : st.droprule)
            if (dr == static_cast<int>(lang::DropStrategy::BAS)) has_bas = true;
    if (L.nclasses == 1 && !open && has_bas) m.push_back("sqd");

    // MVAC recurs on the chains and has no population to remove from an open
    // one, so SolverMVA.m:80 advertises it on ~any(isinf(njobs)) -- FULLY
    // closed, not merely "has a closed class": a mixed model would otherwise
    // pass the gate here and be refused again by solver_mvac_analyzer.
    bool any_open_class = false;
    for (const qn::JobClass& c : L.classes)
        if (std::isinf(c.population)) any_open_class = true;
    if (!any_open_class) m.push_back("mvac");

    // SJN (shortest-job-next, pfqn_mvasjn / pfqn_amvasjn): the conditional
    // waiting time equation is a population recursion, so the family runs on a
    // CLOSED model with an SJF station and nowhere else -- mva_dispatch rejects
    // an open one by name. Advertised only there, for the reason 'sqni' is
    // gated: a listed name must actually run, and an unlisted one that IS
    // dispatched is refused by check_method before the analyzer sees it.
    if (!any_open_class && sn_has_sjn(L)) {
        m.push_back("sjn.mva");
        m.push_back("sjn.amva");
    }

    if (open) {
        // QNA and RQNA are open-network decomposition analyzers; on a model with
        // closed chains they have no fixed point and the reference advertises
        // them only for a fully open model.
        m.insert(m.begin() + 4, "rqt");
        m.insert(m.begin() + 4, "rqna");
        m.insert(m.begin() + 4, "qna");
    } else {
        // Marie is a closed-network method, and is withheld when the classes
        // are not routed alike: measured against CTMC it is exact for a single
        // class and for classes routed identically, and 6.2 off once they visit
        // different queues, so the reference withholds it rather than let it
        // answer wrongly.
        bool classdep_routing = false;
        for (std::size_t r = 0; r + 1 < L.nclasses; ++r)
            for (std::size_t i = 1; i <= L.nof_nodes() && !classdep_routing; ++i)
                for (std::size_t j = 1; j <= L.nof_nodes(); ++j)
                    if (L.route_eff(r + 1, r + 1, i, j) != L.route_eff(r + 2, r + 2, i, j)) {
                        classdep_routing = true;
                        break;
                    }
        if (!classdep_routing) {
            m.push_back("marie");
            m.push_back("amva.marie");
        }
    }
    // priomva: preemptive-resume priority arm (Chandy-Lakshmi [ChaL83]), offered
    // only when a station actually uses FCFSPRPRIO. The arm itself lives in
    // solver_amvald's forward step. Mirrors SolverMVA.m.
    for (const qn::Station<T>& st : L.stations)
        if (st.sched == qn::SchedStrategy::FCFSPRPRIO) {
            m.push_back("priomva");
            m.push_back("amva.priomva");
            break;
        }
        // amva.mapqn: the horizontal-cut MVA for one exponential delay and one FCFS MAP
        // queue (mapqn_amva); offered only on that shape, which mva_mapqn_reason
        // judges for the list and the run alike.
        if (mva_mapqn_reason(L, "mapqn").empty()) {
            m.push_back("amva.mapqn");
        }

    m.push_back("lin");
    m.push_back("egflin");
    m.push_back("gflin");
    m.push_back("amva.lin");

    // The queueing-system closed forms, for an open single-class two-station
    // model, which is exactly the shape the qsys analyzer serves.
    if (open && L.nstations == 2 && L.nclasses == 1) {
        const char* qsys[] = {"mm1",        "mmk",           "mg1",
                              "mgi1",       "gm1",           "gig1",
                              "gim1",       "gig1.kingman",  "gigk",
                              "gigk.kingman_approx",         "gig1.gelenbe",
                              "gig1.heyman", "gig1.kimura",  "gig1.allen",
                              "gig1.kobayashi", "gig1.klb",  "gig1.marchal",
                              "gigk.whitt", "qed", "gig1.extremal", "gigk.diffusion"};
        for (const char* q : qsys) m.push_back(q);
        // The two abandonment methods are listed only when the station actually
        // reneges: they have nothing to say about a queue nobody walks away
        // from, and listing them there would name a method that cannot run.
        const std::size_t qi = detail::station_of_type(L, qn::NodeType::Queue);
        if (qi != 0 && api::sn_patience_handles(L, qi - 1, 0).present) {
            m.push_back("erlanga");
            m.push_back("mgisrgi");
        }
    }

    // A LISTED NAME MUST ACTUALLY RUN, and the rules the feature registry cannot
    // express have to be applied HERE: a product form, a class count and a server
    // count have no registry name, and `auto_family_refusal` reaches them only
    // through this list. Each predicate is the one the analyzer itself raises on,
    // so a listed row is a row that runs and a withheld one is one that would have
    // errored -- or, for the closed-population family on an open model, answered
    // with a table of zeros.
    std::size_t n_inf = 0;
    for (const qn::Station<T>& st : L.stations)
        if (st.sched == qn::SchedStrategy::INF) ++n_inf;
    // The schmidt-ext arm recurs on the CHAIN populations, so the gate asks about
    // those and not about the class ones; see mva_schmidt_ext_reason. Built lazily
    // and once: it is the only rule here that costs a chain aggregation, and every
    // other name in the list would pay for it.
    std::vector<double> sx_n;
    std::vector<bool> sx_fcfs;
    bool sx_ready = false;
    std::vector<std::string> keep;
    keep.reserve(m.size());
    for (const std::string& name : m) {
        if (!mva_closed_population_reason(L, name).empty()) continue;
        if (!mva_single_class_open_reason(L, name).empty()) continue;
        if (!mva_mvac_reason(L, name).empty()) continue;
        if (!mva_mapqn_reason(L, name).empty()) continue;
        if (qn::mva_base_method(name) == "schmidt-ext") {
            if (!sx_ready) {
                const ChainDemands<T> cd = sn_get_demands_chain(L);
                for (std::size_t c = 0; c < L.nchains; ++c) sx_n.push_back(cd.Nchain[c]);
                for (const qn::Station<T>& st : L.stations)
                    if (st.sched != qn::SchedStrategy::INF && st.sched != qn::SchedStrategy::EXT)
                        sx_fcfs.push_back(st.sched == qn::SchedStrategy::FCFS);
                sx_ready = true;
            }
            if (!mva_schmidt_ext_reason(sx_n, sx_fcfs, name).empty()) continue;
        }
        // pfqn_sqni is a closed form for ONE queueing station with a delay; listing
        // it elsewhere named a method solver_amva refuses by name.
        if (qn::mva_base_method(name) == "sqni" && !(L.nstations == 2 && n_inf == 1)) continue;
        keep.push_back(name);
    }
    m.swap(keep);
    return m;
}

/**
 * Port of `runAnalyzerChecks`' method gate: a method the solver does not list
 * is refused before any analyzer sees it.
 */
template <class T>
void check_method(const qn::NetworkStruct<T>& L, const std::string& method) {
    const std::vector<std::string> valid = list_valid_methods(L);
    if (std::find(valid.begin(), valid.end(), method) != valid.end()) return;
    throw UnsupportedError("SolverMVA: the '" + method + "' method is unsupported by this solver");
}

/**
 * Port of `SolverMVA.resolveMethod`: the feature-driven `default` -> `rqna`
 * upgrade for a bursty single-class open network.
 *
 * The reference computes this ONLY to pick the method feature set
 * (`getMethodFeatureSet`), never writing it back to options; `mva_dispatch`
 * makes the same decision again on its own terminal branch. Without it the gate
 * would read the base envelope, which declares neither MAP nor MMPP2, and would
 * refuse the very models `solver_rqna` exists to solve.
 */
template <class T>
std::string resolve_method(const qn::NetworkStruct<T>& L, const std::string& method) {
    if (method != "default" || L.nclasses != 1) return method;
    // JobClass::population is a plain double, as mva_dispatch.h:854 reads it; a
    // num_traits<T> round trip would throw on an infinite one at T = Rational.
    for (std::size_t r = 0; r < L.nclasses; ++r)
        if (std::isfinite(L.classes[r].population)) return method;
    if (api::sn_has_bursty_arrival(L)) return std::string("rqna");
    if (L.nstations == 2) {
        // A single-class open station customers ABANDON is a different model,
        // not a correction to a G/G/k one: the resolution has to happen here as
        // well as in mva_dispatch, because the feature gate runs on the
        // resolved name and Reneging is admitted for these two methods only.
        const std::size_t qi = detail::station_of_type(L, qn::NodeType::Queue);
        if (qi != 0) {
            const api::PatienceHandles<T> h = api::sn_patience_handles(L, qi - 1, 0);
            if (h.present) return h.isExponential ? std::string("erlanga") : std::string("mgisrgi");
        }
    }
    return method;
}

namespace detail {

/**
 * The visit ratios summed over chains, at STATION level.
 *
 * `sn.visits` is indexed by stateful node and by chain; both conversions below
 * need one (nstations x nclasses) matrix, which is what `cellsum` followed by
 * the statefulToStation mapping produces in the reference.
 */
template <class T>
Matrix<T> station_visits(const qn::NetworkStruct<T>& L) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> V(L.nstations, L.nclasses, zero);
    for (std::size_t c = 0; c < L.nchains; ++c)
        for (std::size_t i = 0; i < L.nstations; ++i) {
            const std::size_t sf = L.stateful_of_station(i + 1);
            for (std::size_t k = 0; k < L.nclasses; ++k)
                V(i, k) = T(V(i, k) + L.visits[c](sf - 1, k));
        }
    return V;
}

}  // namespace detail

/**
 * Port of `sn_get_residt_from_respt`: the per-JOB residence time.
 *
 * `RN` is per visit; multiplying by this station's visits and dividing by the
 * reference station's visits of the chain's reference class turns it into the
 * time a job spends here per completion. A response time already below FineTol
 * (an Immediate class) is passed through unscaled, since scaling noise is
 * still noise.
 */
template <class T>
Matrix<T> sn_get_residt_from_respt(const qn::NetworkStruct<T>& L, const Matrix<T>& RN) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, K = L.nclasses;
    const Matrix<T> V = detail::station_visits(L);
    Matrix<T> WN(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            // A PLACE HAS NO SERVICE PROCESS, SO `disabled` IS NOT ABSENCE --
            // the same caveat `refresh_capacity` already carries. The reference
            // gates this loop on the metric HANDLE (`WH{ist,k}.disabled`), which
            // a Place holding a class does have; `sn.disabled` is only a
            // stand-in for it, and on an SPN it is true at every Place, which
            // zeroed the whole ResidT column (`-s jmt` reported ResidT 0 beside
            // RespT 0.35308 on spn_basic_closed while every other codebase
            // reported both).
            const bool is_place = L.stations[i].nodetype == qn::NodeType::Place;
            if ((L.disabled[i][k] && !is_place) || !(RN(i, k) > zero)) continue;
            if (num_traits<T>::to_double(RN(i, k)) < GlobalConstants::FineTol) {
                WN(i, k) = RN(i, k);
                continue;
            }
            std::size_t c = L.nchains;
            for (std::size_t cc = 0; cc < L.nchains; ++cc)
                if (L.chains[cc][k]) c = cc;
            if (c == L.nchains) continue;
            const std::size_t rs = L.classes[k].refstat;
            T den = zero;
            if (L.refclass[c] > 0) {
                den = V(rs - 1, L.refclass[c] - 1);
            } else {
                for (std::size_t r : L.inchain[c]) den += V(rs - 1, r - 1);
            }
            if (den > zero) WN(i, k) = T(RN(i, k) * V(i, k) / den);
        }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            if (num_traits<T>::to_double(WN(i, k)) < 10.0 * GlobalConstants::FineTol)
                WN(i, k) = zero;
    return WN;
}

/**
 * Port of `sn_get_arvr_from_tput`: the arrival rate each station sees, from the
 * throughputs and the class-expanded routing `sn.rt`.
 *
 * A Source's row is left as the routing computes it and zeroed by the caller's
 * mask, which is where `getAvg` does it.
 *
 * A CACHE IS PROPAGATED FIRST, exactly as the reference propagates it in a
 * loop of its own before the generic one. The hit/miss split rides in the
 * effective routing as a SELF-LOOP at the cache node (`set_route_effective(r,
 * hitclass[r], ci, ci, hitprob)`), so the accumulation below fills the cache's
 * hit and miss rows from its own arrival row -- but only if the cache is
 * reached before whatever it feeds. In node order a Router usually precedes
 * the Cache, and taking the nodes in that order left the Router at zero and,
 * through it, every station downstream: on cache_replc_routing the two Delay
 * stations reported ArvR 0 against the reference's 0.4 and 0.6, which then
 * dropped both rows from the table as all-zero.
 */
/**
 * A Join reports the PER-SIBLING waiting time, QLen over the SIBLING ARRIVAL
 * RATE, and not QLen over its own firing rate.
 *
 * The two differ by exactly the fork degree: a Join fires once per parent job
 * while it takes in one sibling per branch, so Little's law applied with the
 * firing rate answers a question about parents with a queue length measured in
 * siblings. On the closed two-branch fork-join with N=2 that is 1.44444 where
 * the settled JMT convention is 0.72222.
 *
 * Called by every engine that folds sibling classes back (CTMC and SSA), after
 * their own `Q/T` pass, so the override lands on the same table the caller
 * reads. A model with no Join leaves the matrix untouched.
 */
template <class T>
void sn_apply_join_respt(const qn::NetworkStruct<T>& L, const Matrix<T>& QN, const Matrix<T>& AN,
                         Matrix<T>& RN) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t nd = 0; nd < L.nodes.size(); ++nd) {
        if (L.nodes[nd].nodetype != qn::NodeType::Join) continue;
        const std::size_t ist = L.nodes[nd].station;
        if (ist == 0 || ist > RN.rows()) continue;
        for (std::size_t r = 0; r < RN.cols(); ++r)
            if (AN(ist - 1, r) > zero) RN(ist - 1, r) = QN(ist - 1, r) / AN(ist - 1, r);
    }
}

template <class T>
Matrix<T> sn_get_arvr_from_tput(const qn::NetworkStruct<T>& L, const Matrix<T>& TN) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, R = L.nclasses;
    const std::size_t S = L.nof_stateful();
    Matrix<T> AN(M, R, zero);
    if (L.rt.rows() != S * R) return AN;  // no routing was refreshed

    // Throughput per stateful node: a station contributes its own, and any
    // other stateful node inherits what the routing carries into it.
    Matrix<T> TS(S, R, zero);
    for (std::size_t sf = 0; sf < S; ++sf) {
        const std::size_t nd = L.stateful_nodes[sf];
        const std::size_t ist = L.nodes[nd - 1].station;
        if (ist == 0) continue;
        for (std::size_t r = 0; r < R; ++r) TS(sf, r) = TN(ist - 1, r);
    }
    for (int pass = 0; pass < 2; ++pass) {
        for (std::size_t sf = 0; sf < S; ++sf) {
            const std::size_t nd = L.stateful_nodes[sf];
            if (L.nodes[nd - 1].station != 0) continue;
            const bool is_cache = L.nodes[nd - 1].nodetype == qn::NodeType::Cache;
            if (is_cache != (pass == 0)) continue;
            for (std::size_t k = 0; k < R; ++k) {
                T acc = zero;
                for (std::size_t sf2 = 0; sf2 < S; ++sf2)
                    for (std::size_t r = 0; r < R; ++r)
                        acc += T(TS(sf2, r) * L.rt(sf2 * R + r, sf * R + k));
                TS(sf, k) = acc;
            }
        }
    }
    for (std::size_t ist = 0; ist < M; ++ist) {
        const std::size_t sf_i = L.stateful_of_station(ist + 1) - 1;
        for (std::size_t k = 0; k < R; ++k) {
            T acc = zero;
            for (std::size_t sf2 = 0; sf2 < S; ++sf2)
                for (std::size_t r = 0; r < R; ++r)
                    acc += T(TS(sf2, r) * L.rt(sf2 * R + r, sf_i * R + k));
            AN(ist, k) = acc;
        }
    }
    return AN;
}

/**
 * Which metric is being filtered. `@@MNetwork/getAvgHandles.m` disables a metric
 * handle per KIND, not per station, and `filterMetric` turns a disabled handle
 * into a zero -- so the kind is what decides which rows survive:
 *
 *   Q, R, W   disabled at a Source and at a Sink
 *   U         disabled at a Source, a Sink, a Fork AND A JOIN -- a Join is an
 *             infinite server whose "utilization" would just restate its queue
 *             length, and the reference declines to report it
 *   T, A      never disabled by station kind, only by an undefined service
 *
 * Every kind is additionally disabled where the class has no service defined.
 */
enum class MetricKind { QLen, Util, RespT, ResidT, Tput, ArvR };

/**
 * Port of `filterMetric`: what `@@NetworkSolver/getAvg` does between the
 * analyzer and the caller.
 *
 * `zero_mask` carries the RN < 10 FineTol test the queue length and the
 * utilization take and the other metrics do not.
 */
template <class T>
Matrix<T> filter_metric(const qn::NetworkStruct<T>& L, const Matrix<T>& metric, MetricKind kind,
                        const std::vector<std::vector<bool>>* zero_mask) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, K = L.nclasses;
    // the fork-join, cache and spawn-fed exemptions from the reachability test
    // below: those three constructs are exactly the ones whose visit equations
    // do not describe how a job actually reaches a station (a Join is entered by
    // the parent class, which has no routing visit there), so a metric the
    // analyzer reports as positive is trusted over the visit ratio
    // ... and the STOCHASTIC PETRI NET, the fourth such construct and for the
    // same reason: where a token sits is the MARKING, and a Place is reached by
    // a transition firing rather than by a routing visit, so `sn.visits` is
    // identically zero over every Place. Without the exemption this loop erased
    // the whole table an analyzer had just computed -- `spnlp.upper` on a
    // fork-join net reported four zeros in place of the token bounds, and
    // `r.QN >= tokens` then failed against a row the bound had got right.
    bool has_fj = false, has_cache = false, has_spn = false;
    for (const qn::NodeDef& nd : L.nodes) {
        if (nd.nodetype == qn::NodeType::Fork || nd.nodetype == qn::NodeType::Join) has_fj = true;
        if (nd.nodetype == qn::NodeType::Cache) has_cache = true;
        if (nd.nodetype == qn::NodeType::Transition) has_spn = true;
    }
    // The Tput/ArvR HANDLE is disabled iff `~isServiceDefined && ~isCacheClass`
    // (getAvgHandles.m): a service-undefined pair keeps its throughput only when
    // the class is a cache hit/miss class -- NOT for fork-join. The fork-join
    // exemption belongs to the reachability test (part 2 below), not here.
    std::vector<bool> is_cache_class(K, false);
    for (const auto& kv : L.nodeparam) {
        const qn::CacheParam<T>& cp = kv.second;
        for (std::size_t r = 0; r < cp.hitclass.size(); ++r) {
            if (cp.hitclass[r] != 0 && cp.hitclass[r] <= K) is_cache_class[cp.hitclass[r] - 1] = true;
            if (r < cp.missclass.size() && cp.missclass[r] != 0 && cp.missclass[r] <= K)
                is_cache_class[cp.missclass[r] - 1] = true;
        }
    }

    Matrix<T> out(M, K, zero);
    for (std::size_t i = 0; i < M; ++i) {
        const qn::NodeType nt = L.stations[i].nodetype;
        const bool src_or_sink = (nt == qn::NodeType::Source || nt == qn::NodeType::Sink);
        bool kind_disabled = false;
        switch (kind) {
            case MetricKind::QLen:
            case MetricKind::RespT:
            case MetricKind::ResidT:
                kind_disabled = src_or_sink;
                break;
            case MetricKind::Util:
                kind_disabled = src_or_sink || nt == qn::NodeType::Fork || nt == qn::NodeType::Join;
                break;
            case MetricKind::Tput:
            case MetricKind::ArvR:
                kind_disabled = false;
                break;
        }
        if (kind_disabled) continue;
        // `hasServiceTunnel`, the wrapper this port was missing. Every metric in
        // `getAvgHandles.m` gates its service test on `if ~hasServiceTunnel(ist)`
        // -- a station whose SERVER is a ServiceTunnel is never disabled for
        // having no service law, only by the kind-specific rules just above
        // (which already match, arm for arm).
        //
        // Those stations are Source, Fork, Join and an ORDINARY Place. A Place is
        // the one that has to be decided rather than looked up: `installQueueServer`
        // REPLACES the tunnel on a queueing Place, so `hasServiceTunnel` is false
        // there, and this struct carries no queueing flag. It is recovered the way
        // MATLAB defines it, PER STATION: an ordinary Place has a law for no class
        // at all, a queueing one has a real server. Deciding it per class instead
        // would keep a queueing Place's explicitly Disabled class, which MATLAB
        // disables.
        //
        // The Place arm is what made a pure Petri net print an EMPTY AvgTable:
        // a non-queueing Place has NaN rates, which is `disabled` here, so every
        // metric was zeroed and the printer drops an all-zero row -- while MATLAB
        // reports e.g. P1 QLen 0.25966 on spn_basic_open. Source and Join change
        // little in practice (a Source's ungenerated class and a Join both compute
        // 0 or are already kept), but they are MATLAB's rule and are encoded here
        // rather than left to coincidence.
        bool service_tunnel = (nt == qn::NodeType::Source || nt == qn::NodeType::Fork ||
                               nt == qn::NodeType::Join);
        if (nt == qn::NodeType::Place) {
            service_tunnel = true;
            for (std::size_t k = 0; k < K && k < L.service[i].size(); ++k)
                if (!L.service[i][k].disabled) { service_tunnel = false; break; }
        }
        for (std::size_t k = 0; k < K; ++k) {
            if (!L.disabled[i][k] || service_tunnel) {
                out(i, k) = metric(i, k);
            } else if (is_cache_class[k] &&
                       (kind == MetricKind::Tput || kind == MetricKind::ArvR)) {
                // a cache over-routing sends a hit/miss class through a station
                // that never serves it (a pass-through); its Tput/ArvR handle is
                // not disabled, so the reported flow is kept
                out(i, k) = metric(i, k);
            }
        }
    }
    if (zero_mask)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k)
                if ((*zero_mask)[i][k]) out(i, k) = zero;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            if (num_traits<T>::to_double(out(i, k)) < GlobalConstants::FineTol) out(i, k) = zero;
    // a class with no visit at a station holds nothing there
    for (std::size_t k = 0; k < K; ++k) {
        std::size_t c = L.nchains;
        for (std::size_t cc = 0; cc < L.nchains; ++cc)
            if (L.chains[cc][k]) c = cc;
        if (c == L.nchains) continue;
        for (std::size_t i = 0; i < M; ++i) {
            if (L.visits[c](L.stateful_of_station(i + 1) - 1, k) != zero) continue;
            if ((has_fj || has_cache || has_spn) &&
                num_traits<T>::to_double(metric(i, k)) > GlobalConstants::FineTol)
                continue;
            out(i, k) = zero;
        }
    }
    return out;
}

/**
 * Port of `SolverMVA.supportsFiniteCapacity` (SolverMVA.m:158-186), the
 * structural capacity gate SolverMVA.supportsModelMethod (SolverMVA.m:143-153)
 * lays on top of the universal feature gate.
 *
 * WHY IT IS SEPARATE FROM `feature_gate`. Finite capacity has no
 * LINE_QN_FEATURE_LIST enumerator: it is a NUMBER on a station, not a
 * construct, so no declared feature set can express it. Without this check a
 * model built with setCapacity is solved as if uncapacitated and returns a
 * WRONG NUMBER instead of a refusal (BUG-39).
 *
 * WHY IT RUNS AFTER THE FEATURE GATE. The reference calls
 * supportsModelMethod@@NetworkSolver first and only tests capacity `if bool`, so
 * a model that fails both is refused by the FEATURE, whose message names the
 * construct to remove. The order is therefore load-bearing, not incidental.
 *
 * THE TWO EXEMPTIONS.
 *  - A BAS model, because solver_mva_analyzer routes it to solver_sqd under
 *    `default` too, so the finite buffers ARE honoured on every MVA path. The
 *    predicate is the api/sn one (single CLASS), not the analyzer's dispatch
 *    test mva_is_bas_model (single CHAIN); see solver_mva.h.
 *  - A single-station M/M/1/K with tail drop under any method but `exact`,
 *    which the mg1k.mgs branch of solver_mva_qsys_analyzer solves. That branch
 *    is an approximation away from scv=1, so `exact` is NOT exempted and must
 *    refuse here. `method` is the RESOLVED method, which is what
 *    runAnalyzerChecks (NetworkSolver.m:184) passes on.
 */
template <class T>
void mva_check_finite_capacity(const qn::NetworkStruct<T>& L, const std::string& method) {
    if (sn_is_bas_model(L)) return;
    if (method != "exact" && sn_is_mm1k_loss(L)) return;
    qn::check_binding_capacity("SolverMVA", L);
}

/**
 * Port of `@@SolverMVA/runAnalyzer.m` for the `lang='matlab'` path: gate, solve,
 * convert, filter.
 *
 * A model with a Fork is solved through the shared fork-join fixed point
 * (fj_driver.h), which drives `mva_dispatch` as its inner solve on the
 * transformed (auxiliary-expanded) model; a model without one runs the dispatch
 * exactly once. The filtering below is common to both paths.
 */
template <class T>
AvgResult<T> solver_mva_run_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt_in,
                            const Matrix<T>& init_sol) {
    check_method(L, opt_in.method);

    // Finite Capacity Region: MVA does not enforce the aggregate per-region job
    // limit and would otherwise silently return the unconstrained product-form
    // answer. Port of `runAnalyzer.m:21-25`.
    if (!L.regions.empty())
        throw UnsupportedError(
            "SolverMVA: this model uses a Finite Capacity Region (addRegion), which is not "
            "supported by SolverMVA. Use SolverJMT, or setCapacity for a single-station limit.");

    // resolveMethod (SolverMVA.resolveMethod) computes the bursty default->rqna
    // upgrade only for the feature gate; the reference never writes it back to
    // options. The substitution itself lives in mva_dispatch's terminal
    // ordinary-QN branch, so the structural special cases (e.g. a 3-node
    // MAP/M/1, a G/M/1 the qsys analyzer solves) keep method='default' rather
    // than being forced into rqna and refused.
    // runAnalyzerChecks' universal feature gate, AFTER the method gate so the
    // specific message wins wherever one exists. Gated on `L`, before fj_mmt
    // rewrites the Fork into a Router that MVA does not declare.
    // The resolved method is passed so that a `default` -> `rqna` refusal names
    // rqna, which is NetworkSolver.m:184-190's second message form: rqna's set
    // withdraws ClosedClass, and nothing the user typed mentions rqna.
    const std::string rm = resolve_method(L, opt_in.method);
    // The MODEL-aware overload, because one grant in getMethodFeatureSet reads
    // the model: every method but 'exact' carries FiniteCapacity on a
    // single-station M/M/1/K with tail drop, the shape mva_dispatch answers
    // under any name. `mva_check_finite_capacity` below is the same rule as a
    // structural refusal, and the two must agree.
    qn::feature_gate("SolverMVA", qn::mva_feature_set(rm, L), L, opt_in.method, rm);
    mva_check_finite_capacity(L, rm);
    const MvaOptions& opt = opt_in;

    MvaSolution<T> s;
    std::string actualmethod, warning;
    // The integrated cacheqn analyzer returns a refreshed struct whose routing
    // carries the actual hit/miss probabilities; ArvR and ResidT are read from
    // it. Empty for every other model, in which case the base struct is used.
    std::shared_ptr<qn::NetworkStruct<T>> refreshed;
    // What a cache branch measured; see AvgResult::cache. Left empty on the
    // fork-join path, whose inner solves are of a TRANSFORMED network.
    solvers::CacheMetrics<T> cache;
    if (L.has_fork()) {
        // The transform turns the fork into a router, the join into a delay and
        // the branches into auxiliary open classes; the driver iterates the
        // auxiliary arrival rates to their fixed point and merges the auxiliary
        // columns back before returning. `mva_dispatch` sees a plain mixed
        // network, so its own fork checks never fire.
        FjMmt<T> tr = fj_fork_join_transform(L, opt.fork_join);
        std::vector<T> lam(tr.V.classes.size() + 1,
                           num_traits<T>::from_double(GlobalConstants::FineTol));
        MvaOptions inner = opt;
        inner.base_has_fork = true;  // the transform removes the fork; see MvaOptions
        std::string am, wn;
        s = fj_fixed_point(L, tr, lam, opt, [&inner, &am, &wn](qn::NetworkStruct<T>& V) {
            const DispatchResult<T> d = mva_dispatch(V, inner, Matrix<T>());
            am = d.actualmethod;
            // the last inner solve is the one the returned metrics came from
            wn = d.warning;
            return d.sol;
        });
        actualmethod = am;
        warning = wn;
    } else {
        const DispatchResult<T> dr = mva_dispatch(L, opt, init_sol);
        s = dr.sol;
        actualmethod = dr.actualmethod;
        warning = dr.warning;
        refreshed = dr.refreshed_struct;
        cache = dr.cache;
    }
    // Values (not topology) for ArvR/ResidT come from the refreshed struct when
    // the cacheqn analyzer supplied one; filter_metric still keys on the base L
    // so its cache pass-through masking is unaffected by the ClassSwitch relabel.
    const qn::NetworkStruct<T>& refL = refreshed ? *refreshed : L;

    const std::size_t M = L.nstations, K = L.nclasses;
    std::vector<std::vector<bool>> mask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            mask[i][k] = num_traits<T>::to_double(s.R(i, k)) < 10.0 * GlobalConstants::FineTol;

    // `getAvg` zeroes the arrival rate at a Source explicitly: a Source is not
    // fed by anything, and the routing would otherwise report whatever the
    // stochastic complement carries back into it.
    std::vector<std::vector<bool>> srcmask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Source)
            for (std::size_t k = 0; k < K; ++k) srcmask[i][k] = true;

    AvgResult<T> out;
    // Carried, not consumed and dropped: the node table needs the same visits
    // the station table just used, or the two disagree on the cache columns.
    out.refreshed_struct = refreshed;
    // The cache answer a cache branch measured. Empty on every other model, and
    // on the fork-join path, whose inner solves are of a TRANSFORMED network.
    out.cache = cache;
    out.QN = filter_metric(L, s.Q, MetricKind::QLen, &mask);
    out.UN = filter_metric(L, s.U, MetricKind::Util, &mask);
    out.RN = filter_metric(L, s.R, MetricKind::RespT, nullptr);
    out.TN = filter_metric(L, s.Tp, MetricKind::Tput, nullptr);
    // Both ArvR and ResidT come from the refreshed struct when the cacheqn
    // analyzer supplied one: its cache self-switch is normalized at the ACTUAL
    // hit/miss split, so the visits back the carried rate without the over-route
    // inflation. filter_metric still keys on the base L for its topology masking.
    out.WN = filter_metric(L, sn_get_residt_from_respt(refL, out.RN), MetricKind::ResidT, nullptr);
    out.AN = filter_metric(L, sn_get_arvr_from_tput(refL, out.TN), MetricKind::ArvR, &srcmask);

    // Unstable open-station cap (@@NetworkSolver/getAvg.m:198-235): a finite-server
    // queueing station whose open-class offered load sum_r TN/(c*rate) >= 1 is
    // saturated. MATLAB normalizes its utilization to sum 1.0 and reports queue
    // length and response time as Inf; it does NOT rescale throughput (Python
    // native does -- a MATLAB-vs-Python divergence -- and MATLAB is ground truth).
    // Skipped under exact arithmetic, where infinity is not representable.
    if constexpr (num_traits<T>::has_transcendental) {
        const double dinf = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < M; ++i) {
            const double c = num_traits<T>::to_double(L.stations[i].nservers);
            const qn::NodeType nt = L.stations[i].nodetype;
            const qn::SchedStrategy sch = L.stations[i].sched;
            if (!std::isfinite(c) || c <= 0.0) continue;  // delay / infinite server
            if (nt == qn::NodeType::Source || nt == qn::NodeType::Sink) continue;
            if (sch == qn::SchedStrategy::INF || sch == qn::SchedStrategy::EXT) continue;
            std::vector<double> rho(K, 0.0);
            double rho_open = 0.0, rho_tot = 0.0;
            for (std::size_t r = 0; r < K; ++r) {
                const double rate = num_traits<T>::to_double(L.rates(i, r));
                const double tp = num_traits<T>::to_double(out.TN(i, r));
                if (rate > 0.0 && tp > 0.0) rho[r] = tp / (c * rate);
                rho_tot += rho[r];
                if (std::isinf(L.classes[r].population)) rho_open += rho[r];
            }
            if (rho_open < 1.0 || rho_tot <= 0.0) continue;
            for (std::size_t r = 0; r < K; ++r) {
                out.UN(i, r) = num_traits<T>::from_double(rho[r] / rho_tot);
                if (std::isinf(L.classes[r].population) && rho[r] > 0.0) {
                    out.QN(i, r) = num_traits<T>::from_double(dinf);
                    out.RN(i, r) = num_traits<T>::from_double(dinf);
                }
            }
        }
    }
    out.CN = s.C;
    out.XN = s.X;
    out.method = opt.method;
    out.actualmethod = actualmethod;
    // Immediate feedback is APPROXIMATED, not refused: a fed-back job holds its
    // server, and mean-value analysis has no way to express that, so the visit
    // it makes is counted as an ordinary re-entry. The reference warns and
    // solves (`runAnalyzer.m:26-27`); refusing here would reject a model MATLAB
    // and the JAR both answer. Appended rather than assigned, so a warning the
    // dispatch already raised is not swallowed by this one.
    if (L.has_immediate_feedback()) {
        if (!warning.empty()) warning += " ";
        warning +=
            "SolverMVA does not handle immediate feedback (immfeed); the solver will treat "
            "self-loops as class-switching with re-queueing.";
    }
    // Verbatim: a user must be able to match it against the reference output.
    out.warning = warning;
    out.iter = s.iter;
    out.converged = s.converged;
    return out;
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_RUNNER_H
