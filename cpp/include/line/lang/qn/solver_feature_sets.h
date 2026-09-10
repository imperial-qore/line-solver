/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_SOLVER_FEATURE_SETS_H
#define LINE_LANG_QN_SOLVER_FEATURE_SETS_H

/**
 * The DECLARED side of the gate: one feature set per solver.
 *
 * SOURCE. Each function transcribes that solver's MATLAB declaration:
 *   MVA    matlab/src/solvers/MVA/@@SolverMVA/SolverMVA.m:188-217
 *          per-method deltas SolverMVA.m:126-141
 *   NC     matlab/src/solvers/NC/@@SolverNC/SolverNC.m:180-209
 *   CTMC   matlab/src/solvers/CTMC/@@SolverCTMC/SolverCTMC.m:126-168
 *   SSA    matlab/src/solvers/SSA/@@SolverSSA/SolverSSA.m:67-110
 *   FLUID  matlab/src/solvers/FLD/@@SolverFLD/SolverFLD.m:88-119
 *   MAM    matlab/src/solvers/MAM/@@SolverMAM/SolverMAM.m:141-184 (four setTrue calls, unioned)
 *   ENV    matlab/src/solvers/@@SolverENV/SolverENV.m:1399-1444 (inline, inside supports)
 *   LDES   matlab/src/solvers/wrappers/LDES/@@SolverLDES/SolverLDES.m:108-186
 *   BA     matlab/src/solvers/BA/@@SolverBA/SolverBA.m (getFeatureSet)
 *   QNS    matlab/src/solvers/wrappers/QNS/@@SolverQNS/SolverQNS.m (getFeatureSet)
 *
 * WHY `method`. MATLAB's gate is METHOD-aware, not solver-aware:
 * NetworkSolver.runAnalyzerChecks resolves the method first and then reads
 * getMethodFeatureSet(method), which only SolverMVA specialises in the
 * reference (SolverMVA.m:126-141). The other five take the argument and ignore
 * it. MAM's C++ port is the one exception: mam_feature_set reads `method`
 * because its LoadDependence declaration is a WIDENING over MATLAB (see below)
 * and only holds for the methods that actually consume st.lldscaling; MATLAB
 * itself has no such delta to transcribe.
 *
 * WHY THE SETS ARE NOT ALWAYS THE MATLAB ONES. feature_set.h:30-33: a feature
 * DECLARED but not implemented yields a wrong number, an undeclared one yields
 * a clean refusal. So each divergence below narrows the MATLAB set to what the
 * C++ code actually handles, and carries the C++ refusal site that proves it.
 * The one WIDENING left is MAM's LoadDependence, method-gated as above, which is
 * a false refusal otherwise. MVA's five size-based disciplines were a second
 * widening until they were removed; see mva_feature_set below. Verified by
 * diffing each of the nine narrowed sets against its MATLAB declaration name by
 * name, not by inspection.
 *
 * LDES AND QNS ARE THE EXCEPTIONS TO THAT PARAGRAPH, and only because neither is
 * an engine here: SolverLDES forwards the model document to the SSJ engine, and
 * SolverQNS writes the JMVA document `qnsolver` reads or converts the model for
 * `lqns`. In both the code that implements a declared feature is the external
 * binary the reference ALSO drives, so `ldes_feature_set` and `qns_feature_set`
 * are the MATLAB declarations transcribed WHOLE; narrowing either would refuse
 * models the reference answers with the very process this port is about to
 * start.
 *
 * MATLAB FEATURE NAMES WITH NO LINE_QN_FEATURE_LIST ENUMERATOR: none. Every
 * name in all twelve declarations, and in the MVA per-method deltas, has an
 * exact 1:1 registry entry, so the mapping is the identity and needs no alias.
 * `Normal` is the one name dropped on purpose rather than for want of an
 * enumerator; `qns_feature_set` records why.
 */

#include <string>

#include "line/lang/qn/feature_set.h"

namespace line {
namespace qn {

/**
 * SolverMVA.getFeatureSet, plus getMethodFeatureSet's per-method deltas.
 *
 * SolverMVA.m:189-213 declares 59 names; this set is those 59, name for name.
 * RQNA natively consumes non-renewal MAP/MMPP/RAP processes but is open-network
 * only; QNA is a two-moment open-network method.
 *
 * THE SIZE-BASED NAMES WERE REMOVED ONCE AND ARE BACK, AND THE REFERENCE IS WHY.
 * SchedStrategy_SRPT, _PSJF, _FB, _LRPT and _SETF were dropped on 2026-07-29 on
 * the ground that MATLAB's getFeatureSet declared none of them, making them a
 * WIDENING of the reference envelope rather than a transcription of it.
 *
 * That premise expired. MATLAB's SolverMVA.m now declares all five ("Listed
 * only now that the analyzer runs"), and the JAR and native Python gained
 * solver_mva_qsys_sizebased_analyzer and the same five names on 2026-08-16.
 * Keeping them out made C++ NARROWER than the reference: branch 3 of
 * mva_dispatch.h and the qsys_mg1_srpt/_psjf/_fb/_lrpt/_setf closed forms
 * behind it were reachable code held behind a closed gate, and a size-based
 * model the other three solve was refused here.
 *
 * The lesson is about the shape of the note, not the names: a featset comment
 * that records what the reference declared ON A DATE goes stale silently. When
 * the reference moves, this file has to move with it.
 *
 * Every OTHER declared name was checked against the code that consumes it, not
 * assumed:
 *   Cache, CacheClassSwitcher          mva_dispatch.h:754, :765
 *   CacheRetrieval                     mva_dispatch.h:741 (open); the CLOSED variant
 *                                      is refused by name at :734-740, a shape
 *                                      restriction the boolean cannot express
 *   ReplacementStrategy_LRU, _HLRU     solver_mva_cache.h:144-165
 *   SchedStrategy_OI, _PAS             mva_dispatch.h, branch 1, solver_mva_oi.h
 *   SchedStrategy_POLLING              mva_dispatch.h, branch 5, solver_mva_polling.h
 *   SchedStrategy_SJF                  mva_dispatch.h, branch 0, solver_mva_sjn.h; the
 *                                      OPEN case is refused by name there, a shape
 *                                      restriction the boolean cannot express
 *   SchedStrategy_SRPT, _PSJF, _FB,    mva_dispatch.h, branch 3,
 *   _LRPT, _SETF                       solver_mva_qsys_sizebased
 *   Fork, Forker, Join, Joiner         solver_mva_runner.h:414-426, fj_mmt + fj_fixed_point
 *   LoadDependence                     solver_mva.h:543-561, :655-662, :1486-1487
 *   ClassDependence                    solver_mva.h:553-560, :671-684 (pfqn_cdfun)
 * Pareto, Weibull, Lognormal, Uniform, Det, Replayer, BMAP and MMAP carry no
 * consumer of their own by design: there is no `ProcessType::` reference
 * anywhere under cpp/include/line/solvers/mva/, because MVA is a two-moment
 * method that reads only rates and SCVs. That is exactly how MATLAB's MVA
 * treats them, which is why it declares them. BMAP is the weakest of the eight:
 * the batch structure is dropped to two moments rather than refused, in both
 * codebases.
 */
/**
 * The method name with any leading "amva." alias stripped.
 *
 * The mva registry strips the prefix before it selects an algorithm, so every
 * rule keyed on a method name has to see the canonical spelling or it would
 * hold for "bs" and not for "amva.bs".
 */
inline std::string mva_base_method(const std::string& method) {
    return (method.compare(0, 5, "amva.") == 0) ? method.substr(5) : method;
}

/**
 * The AMVA algorithms whose recursion is over a CLOSED population vector.
 *
 * Each approximates the arrival-instant queue length E[Q(N-1_r)] from E[Q(N)]
 * and is handed (L, N, Z) alone, with no arrival rate and no rate function, so
 * an open chain gives it nothing to recur on and the recursion falls out with
 * every metric at zero. ONE list: `mva_feature_set` drops OpenClass for these,
 * `mva::mva_closed_population_reason` refuses them by name at solve time, and
 * `mva::list_valid_methods` withholds them where they cannot run.
 */
inline bool mva_is_closed_population_method(const std::string& method) {
    static const char* kNames[] = {"bs", "aql", "qsa", "sqni", "tay", "scat", "lcp", "chow",
                                   "pamb", "pami", "pamt", "clust", "dmlin", "ab",
                                   "schmidt", "schmidt-ext"};
    const std::string base = mva_base_method(method);
    for (const char* n : kNames)
        if (base == n) return true;
    return false;
}

/**
 * Drops every scheduling name OUTSIDE the BCMP set {INF, PS, FCFS, SIRO,
 * LCFS-PR} that the base MVA envelope declares.
 *
 * The chain algorithms that walk the stations one by one -- the summation
 * method, MVAC and QNA -- accept the BCMP set and refuse the rest, so each drops
 * these from its own envelope rather than answering a station it has no arm for.
 */
inline void mva_unset_non_bcmp_sched(FeatureSet& f) {
    f.unset(Feature::SchedStrategy_HOL);
    f.unset(Feature::SchedStrategy_DPS);
    f.unset(Feature::SchedStrategy_FCFSPRPRIO);
    f.unset(Feature::SchedStrategy_LCFS);
    f.unset(Feature::SchedStrategy_POLLING);
    f.unset(Feature::SchedStrategy_SJF);
    f.unset(Feature::SchedStrategy_SRPT);
    f.unset(Feature::SchedStrategy_PSJF);
    f.unset(Feature::SchedStrategy_FB);
    f.unset(Feature::SchedStrategy_LRPT);
    f.unset(Feature::SchedStrategy_SETF);
    f.unset(Feature::SchedStrategy_OI);
    f.unset(Feature::SchedStrategy_PAS);
}

inline FeatureSet mva_feature_set(const std::string& raw_method) {
    const std::string method = mva_base_method(raw_method);
    FeatureSet f;
    f.set({Feature::Sink, Feature::Source});
    f.set({Feature::ClassSwitch, Feature::Delay, Feature::DelayStation, Feature::Queue});
    f.set({Feature::APH, Feature::Coxian, Feature::Cox2, Feature::Erlang, Feature::Exp,
           Feature::HyperExp, Feature::BMAP});
    f.set({Feature::Pareto, Feature::Weibull, Feature::Lognormal, Feature::Uniform, Feature::Det});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer, Feature::SharedServer,
           Feature::Buffer, Feature::Dispatcher});
    f.set({Feature::CacheClassSwitcher, Feature::Cache});
    f.set(Feature::CacheRetrieval);
    f.set({Feature::Server, Feature::JobSink, Feature::RandomSource, Feature::ServiceTunnel});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS});
    f.set({Feature::SchedStrategy_DPS, Feature::SchedStrategy_FCFS, Feature::SchedStrategy_SIRO,
           Feature::SchedStrategy_HOL, Feature::SchedStrategy_FCFSPRPRIO});
    f.set({Feature::SchedStrategy_LCFS, Feature::SchedStrategy_LCFSPR,
           Feature::SchedStrategy_POLLING});
    f.set({Feature::SchedStrategy_OI, Feature::SchedStrategy_PAS});
    // Closed models only, solver_mva_sjn.h; the open case is refused by name.
    f.set(Feature::SchedStrategy_SJF);
    // Size-based M/G/1, mva_dispatch.h branch 3, solver_mva_qsys_sizebased
    // (Wierman and Harchol-Balter, SIGMETRICS 2003)
    f.set({Feature::SchedStrategy_SRPT, Feature::SchedStrategy_PSJF,
           Feature::SchedStrategy_FB, Feature::SchedStrategy_LRPT,
           Feature::SchedStrategy_SETF});
    // JoinPartial: the MMT fixed point charges the k-th branch completion
    // (fj_ordstat_exp), not the maximum
    f.set({Feature::Fork, Feature::Forker, Feature::Join, Feature::Joiner,
           Feature::JoinPartial});
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND});
    f.set({Feature::ReplacementStrategy_RR, Feature::ReplacementStrategy_FIFO,
           Feature::ReplacementStrategy_LRU, Feature::ReplacementStrategy_HLRU});
    f.set(Feature::MMAP);
    f.set({Feature::ClosedClass, Feature::SelfLoopingClass, Feature::OpenClass, Feature::Replayer});
    // lldscaling: solver_mva.h:543-561 and :1486-1487; cdscaling: :553-560 via pfqn_cdfun.
    f.set({Feature::LoadDependence, Feature::ClassDependence});
    f.set(Feature::JointDependence);
    // c-server stations, SolverMVA.m:1214-1221: the exact recursion, every AMVA
    // kernel, qna/rqt and the M/M/k and G/G/k closed forms carry the count; the
    // per-method block below withdraws it from the single-server names.
    // FiniteCapacity is deliberately NOT here, as in the reference: it is
    // granted per method, and on an M/M/1/K loss model per MODEL.
    f.set(Feature::MultiServer);
    if (mva_is_closed_population_method(method)) {
        // The closed-population AMVA family estimates the arrival-instant queue
        // length as a function of the population vector N and is handed (L, N, Z)
        // alone, so an open chain gives it nothing to recur on: solver_amva has no
        // arm for any of these outside its closed product-form branch, and falling
        // through returned the qd-family answer under their name. Strict product
        // form has no registry name and is applied by
        // mva::mva_closed_population_reason instead.
        f.unset(Feature::OpenClass);
    }
    // The load-dependent analyzer serves a load-, class- or joint-dependent model
    // through "exact"/"mva" (load dependence only, it has no class- or
    // joint-dependent recursion) and through the default/amva/qd/lin/qdlin arms,
    // and refuses every other name by name. The queueing-system closed forms are
    // intercepted by mva_dispatch upstream of that analyzer and keep the base
    // envelope.
    {
        static const char* kLdRefused[] = {"sum", "esum", "mvac", "qli", "fli",
                                           "gflin", "egflin", "qna", "rqna", "rqt"};
        bool ld_refused = mva_is_closed_population_method(method);
        for (const char* n : kLdRefused)
            if (method == n) ld_refused = true;
        if (ld_refused) {
            f.unset(Feature::LoadDependence);
            f.unset(Feature::ClassDependence);
            f.unset(Feature::JointDependence);
        } else if (method == "mva" || method == "exact") {
            f.unset(Feature::ClassDependence);
            f.unset(Feature::JointDependence);
        }
    }
    if (method != "default" && method != "exact") {
        // An order-independent or pass-and-swap station is served by solver_mva_oi
        // alone, which mva_dispatch reaches only under "default" or "exact"; every
        // other name is refused there by name, so it must not be advertised here.
        f.unset(Feature::SchedStrategy_OI);
        f.unset(Feature::SchedStrategy_PAS);
    }
    if (method == "sum" || method == "esum") {
        // solver_mva_sum passes each station to the summation kernel as an INF, PS,
        // LCFS-PR, FCFS or SIRO centre and refuses every other discipline by name.
        mva_unset_non_bcmp_sched(f);
    } else if (method == "mvac") {
        // pfqn_mvac recurs on the closed chains over single-server fixed-rate
        // (SSFR) and infinite-server centres; solver_mvac refuses the rest by name.
        f.unset(Feature::OpenClass);
        mva_unset_non_bcmp_sched(f);
    }
    if (method == "mapqn") {
        // The horizontal-cut MVA consumes a MAP service natively (a closed delay +
        // FCFS queue model, see mva::mva_mapqn_reason); MAP is granted here alone.
        f.set({Feature::MAP, Feature::MMPP2});
        for (Feature x : {Feature::OpenClass, Feature::Source, Feature::Sink, Feature::Fork, Feature::Forker,
                          Feature::Join, Feature::Joiner, Feature::JoinPartial, Feature::ClassSwitch,
                          Feature::StatelessClassSwitcher, Feature::Cache, Feature::CacheClassSwitcher,
                          Feature::CacheRetrieval, Feature::LoadDependence, Feature::ClassDependence,
                          Feature::JointDependence, Feature::SchedStrategy_PS, Feature::SchedStrategy_SIRO,
                          Feature::SchedStrategy_LCFSPR, Feature::SchedStrategy_SRPT, Feature::SchedStrategy_PSJF,
                          Feature::SchedStrategy_FB, Feature::SchedStrategy_LRPT, Feature::SchedStrategy_SETF,
                          Feature::SchedStrategy_OI, Feature::SchedStrategy_PAS})
            f.unset(x);
        mva_unset_non_bcmp_sched(f);
    }
    if (method == "rqna") {
        // RAP is INERT: ProcessType carries no RAP, so it can never be emitted.
        f.set({Feature::MAP, Feature::MMPP2, Feature::MMAP, Feature::RAP});
        f.unset(Feature::ClosedClass);
        f.unset(Feature::SelfLoopingClass);
    } else if (method == "rqt") {
        // Robust queueing theory: single-class open networks, the primitives
        // entering the uncertainty sets are two moments.
        f.unset(Feature::ClosedClass);
        f.unset(Feature::SelfLoopingClass);
    }
    if (method == "qna") {
        // Round-robin dispatching enters as a deterministic traffic split
        // (npfqn_traffic_split_rr), which SolverMVA.m:149 sets true for this
        // method and the exact-MVA paths have no counterpart for; the C++ set
        // omitted this and only carried the two setFalse calls below.
        f.set(Feature::RoutingStrategy_RROBIN);
        f.unset(Feature::ClosedClass);
        f.unset(Feature::SelfLoopingClass);
        // solver_qna's station loop has an arm for INF, PS and FCFS and none for
        // anything else, so on a SIRO, LCFS-PR, HOL or priority station it left
        // that row of Q, U, R and T at zero and reported the table as a solution.
        f.unset(Feature::SchedStrategy_SIRO);
        f.unset(Feature::SchedStrategy_LCFSPR);
        mva_unset_non_bcmp_sched(f);
    } else if (method == "erlanga" || method == "mgisrgi") {
        // The only analytical methods in LINE that carry an abandonment rate
        // (qsys_erlanga, qsys_mgisrgi_whitt). Reneging stays OUT of the base
        // MVA envelope: every other method here would silently ignore the
        // patience law and report the no-abandonment answer.
        f.set(Feature::Reneging);
        f.unset(Feature::ClosedClass);
        f.unset(Feature::SelfLoopingClass);
    }
    if (method == "rqna" || method == "rqt") {
        // A Join is a synchronisation node, not a queue: it carries no service
        // process, so the index-of-dispersion curve these two read off every
        // station does not exist for it, and neither analyzer has a
        // synchronisation term to put in its place. QNA keeps Fork/Join -- its
        // station loop has an explicit Join arm.
        f.unset(Feature::Fork);
        f.unset(Feature::Forker);
        f.unset(Feature::Join);
        f.unset(Feature::Joiner);
        f.unset(Feature::JoinPartial);
    }
    // THE SINGLE-SERVER RECURSIONS, SolverMVA.m:378-382. AQL, QSA and Tay
    // (mva_closed_population_reason), MVAC's SSFR chain recursion, RQNA's
    // GI/G/1 workload, Kant's SJN recursion and the single-server closed forms
    // of the queueing-system analyzer, i.e. every M/G/1, G/M/1 and G/G/1 name.
    // Each structural predicate stays and words the refusal for the run; the
    // delta is what makes it nameable. RQT, QNA, M/M/k, G/G/k and the rest of
    // the envelope carry a server count.
    {
        static const char* kSingleServer[] = {"aql", "qsa", "tay", "mvac", "rqna", "sjn.mva",
                                              "sjn.amva", "mm1", "mg1", "mgi1", "gm1", "gim1"};
        bool single = method.compare(0, 4, "gig1") == 0;
        for (const char* n : kSingleServer)
            if (method == n) single = true;
        if (single) f.unset(Feature::MultiServer);
    }
    // FINITECAPACITY, SolverMVA.m:384-399, is NOT in the base envelope: the
    // product-form recursions solve a buffer away, which is what
    // `mva_check_finite_capacity` refuses. 'default' and 'sqd' reach solver_sqd,
    // the one Blocking-After-Service arm. The single-station M/M/1/K with tail
    // drop is served by the qsys_mg1k_loss_mgs branch of the queueing-system
    // analyzer, which mva_dispatch reaches under EVERY name on that shape,
    // 'exact' excepted since the branch is exact at scv=1 only; that grant is
    // judged on the MODEL and lives in the overload below, because no name
    // carries it.
    if (method == "default" || method == "sqd") f.set(Feature::FiniteCapacity);
    return f;
}

/**
 * `mva_feature_set` on a MODEL, i.e. the whole of SolverMVA.getMethodFeatureSet.
 *
 * The reference's last arm reads `self.model` (SolverMVA.m:396-399): every
 * method but 'exact' is granted FiniteCapacity on a single-station M/M/1/K with
 * tail drop, because mva_dispatch ORs `sn_is_mm1k_loss` into its single-station
 * test and answers that shape whatever the caller named. A name cannot carry
 * the grant, so a caller that HAS the struct asks this overload and one that
 * does not gets the static table above.
 */
template <class T>
FeatureSet mva_feature_set(const std::string& raw_method, const NetworkStruct<T>& sn) {
    FeatureSet f = mva_feature_set(raw_method);
    if (mva_base_method(raw_method) != "exact" && sn.is_mm1k_loss())
        f.set(Feature::FiniteCapacity);
    return f;
}

/**
 * SolverNC.getFeatureSet, 48 names, transcribed unchanged.
 *
 * `Region` is declared DELIBERATELY (SolverNC.m:185-190): NC solves the open
 * single-Delay loss network exactly and a boolean cannot express that split, so
 * the queueing-station case stays as the imperative refusal at
 * solver_nc_runner.h:315-332. DPS, GPS and HOL appear under cpp/.../nc/ but are
 * absent here because MATLAB's NC refuses them and MATLAB is ground truth.
 */
inline FeatureSet nc_feature_set(const std::string& method) {
    FeatureSet f;
    f.set({Feature::Sink, Feature::Source});
    f.set({Feature::ClassSwitch, Feature::Delay, Feature::DelayStation, Feature::Queue});
    f.set({Feature::APH, Feature::Coxian, Feature::Cox2, Feature::Erlang, Feature::Det,
           Feature::Exp, Feature::HyperExp});
    // Geometric is admitted for the discrete-time route only (opt.slotted,
    // solver_nc_dt); on the continuous-time routes it is treated by its mean
    // and SCV like any other renewal law.
    f.set(Feature::Geometric);
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer});
    f.set({Feature::SharedServer, Feature::Buffer, Feature::Dispatcher});
    f.set(Feature::Region);
    f.set({Feature::Server, Feature::JobSink, Feature::RandomSource, Feature::ServiceTunnel});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS, Feature::SchedStrategy_SIRO});
    f.set({Feature::SchedStrategy_LCFS, Feature::SchedStrategy_LCFSPR});
    // DPS is served ONLY in Morrison's closed think+DPS shape (nc_is_dps_model
    // -> solver_nc_dps). A boolean feature cannot express that restriction, so
    // the runner keeps an imperative check for every other DPS model, the same
    // pattern as Region above.
    f.set(Feature::SchedStrategy_DPS);
    // SDR is state dependent yet keeps a product form of its own (eq. 16),
    // which solver_nc_sdr evaluates exactly; see _kb/16-state-dependent-routing.md
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND,
           Feature::RoutingStrategy_SDR});
    f.set({Feature::SchedStrategy_FCFS, Feature::SchedStrategy_OI, Feature::SchedStrategy_PAS});
    // JoinPartial: the MMT fixed point charges the k-th branch completion
    // (fj_ordstat_exp), not the maximum
    f.set({Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner,
           Feature::JoinPartial});
    f.set({Feature::ClosedClass, Feature::SelfLoopingClass});
    f.set({Feature::Cache, Feature::CacheClassSwitcher, Feature::OpenClass});
    f.set(Feature::CacheRetrieval);
    f.set(Feature::CacheItemSize);  // exact and sampling methods only
    f.set({Feature::ReplacementStrategy_RR, Feature::ReplacementStrategy_FIFO});
    f.set(Feature::ReplacementStrategy_HLRU);
    // Petri nets: the 'rec' route (solver_nc_spn) walks the reachable set in a
    // decision diagram, so a Place is a token container rather than a station
    // with a service process. A queueing Place is refused by spn_pf, which is
    // where the product-form class is decided.
    f.set({Feature::Place, Feature::Transition, Feature::Linkage, Feature::Enabling,
           Feature::Inhibiting, Feature::Timing, Feature::Firing, Feature::Storage});
    f.set({Feature::LoadDependence, Feature::ClassDependence});
    f.set(Feature::JointDependence);
    // c-server stations, SolverNC.m:454-460: every route but 'divdiff' folds the
    // count into its own kernel. FiniteCapacity is deliberately NOT here, as in
    // the reference: 'mem', 'default' and 'exact' are granted it below and the
    // rest solve a buffer away.
    f.set(Feature::MultiServer);

    // PER-METHOD DELTAS. A feature set declares what the method ACCEPTS, so it
    // can refuse a model for HAVING a construct and never for lacking one:
    // "closed population only" and "no think time" are sayable here, while
    // "requires a cache" or "requires a loss network" are not and live in
    // `nc::nc_method_refusal`, which the runner and `auto_family_refusal` ask.
    if (method == "divdiff") {
        // The divided-difference closed form of Casale (SIGMETRICS 2017), Eqs.
        // (15)-(16), covers load-independent queues; a think time would ask for
        // the integral form of Cor. 3.4, which is not implemented, so pfqn_nc and
        // pfqn_ncld both refuse one by name. An infinite server is where a think
        // time comes from, so the envelope drops it.
        f.unset(Feature::SchedStrategy_INF);
    } else if (method == "rd" || method == "nrp" || method == "nrl" || method == "nre" ||
               method == "comomld" || method == "panald" || method == "is") {
        // The load-dependent normalizing-constant evaluators are reached by
        // `solver_ncld` only on its CLOSED branch, where pfqn_ncld reads the
        // method name; an open chain sends the model to the mixed route, which
        // never reads it, so every one of these names silently became 'ncldmx'.
        // 'is' is here for the same reason on the load-independent side: the
        // sample-an-ordering estimator integrates over a closed population
        // simplex and has no open-class form.
        f.unset(Feature::OpenClass);
    }
    // MULTISERVER, SolverNC.m:236-243: the divided-difference closed form of
    // 'divdiff' covers load-independent SINGLE-server queues, and
    // `nc::nc_method_refusal` keeps wording why (a c-server station enters the
    // constant as Seidmann's surrogate delay).
    if (method == "divdiff") f.unset(Feature::MultiServer);
    // FINITECAPACITY, SolverNC.m:244-254: 'mem' represents the buffer as a
    // GE/GE/c/0;N queue (`solver_nc_mem_supports`.blocking) and the single-station
    // M/M/1/K with tail drop is solved in closed form (qsys_mm1k_loss) under
    // 'default' and 'exact'. The shape half of each rule stays structural, in
    // `solver_nc_runner`'s own `check_binding_capacity` call.
    if (method == "mem" || method == "default" || method == "exact")
        f.set(Feature::FiniteCapacity);
    return f;
}

/**
 * SolverCTMC.getFeatureSet, the reference's 104 MATLAB names in full.
 *
 * Balking, Reneging and Breakdown stay declared and are INERT: they are call
 * parameters of state_events.h rather than struct fields, so
 * used_lang_features cannot emit them and the declaration is unenforced.
 */
inline FeatureSet ctmc_feature_set(const std::string& method) {
    FeatureSet f;
    f.set({Feature::Source, Feature::Sink});
    // Router is INERT: used_lang_features emits a Router's routing strategy and
    // not the node name, as getUsedLangFeatures.m:77-78 does, so nothing can
    // emit it. Kept declared rather than deleted -- it is what SolverCTMC.m
    // declares, and the AUTO ranking reads these sets by name.
    f.set({Feature::ClassSwitch, Feature::Delay, Feature::DelayStation, Feature::Queue,
           Feature::Router});
    // ME is INERT: ProcessType carries no ME, so it can never be emitted.
    f.set({Feature::MAP, Feature::APH, Feature::MMPP2, Feature::MMAP, Feature::PH,
           Feature::Coxian, Feature::Cox2, Feature::Erlang, Feature::Exp, Feature::HyperExp,
           Feature::ME});
    f.set({Feature::Det, Feature::Gamma, Feature::Lognormal, Feature::Pareto, Feature::Uniform,
           Feature::Weibull});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer, Feature::SharedServer,
           Feature::Buffer, Feature::Dispatcher});
    // The blocked-job overflow buffer is part of the chain, so the FCR is exact.
    f.set(Feature::Region);
    f.set({Feature::Cache, Feature::CacheClassSwitcher});
    f.set(Feature::CacheRetrieval);
    f.set({Feature::Server, Feature::JobSink, Feature::RandomSource, Feature::ServiceTunnel});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS});
    f.set({Feature::SchedStrategy_DPS, Feature::SchedStrategy_GPS});
    f.set({Feature::SchedStrategy_SIRO, Feature::SchedStrategy_SEPT});
    f.set({Feature::SchedStrategy_LEPT, Feature::SchedStrategy_FCFS});
    f.set({Feature::SchedStrategy_HOL, Feature::SchedStrategy_LCFS});
    f.set({Feature::SchedStrategy_LCFSPR, Feature::SchedStrategy_LCFSPRPRIO,
           Feature::SchedStrategy_FCFSPRPRIO});
    // 2026-07-29: the rest of the preempt family, added to SolverCTMC.m:159 in
    // the same change. after_event_station_arv and _dep both carry a
    // buffer_is_tag_phase_pairs arm covering all eight, mirroring
    // afterEventStation.m:1043-1330, so declaring only three gated off five
    // chains the generator builds correctly.
    f.set({Feature::SchedStrategy_LCFSPI, Feature::SchedStrategy_FCFSPR,
           Feature::SchedStrategy_FCFSPI});
    f.set({Feature::SchedStrategy_LCFSPIPRIO, Feature::SchedStrategy_FCFSPIPRIO});
    f.set({Feature::SchedStrategy_PSPRIO, Feature::SchedStrategy_DPSPRIO,
           Feature::SchedStrategy_GPSPRIO});
    f.set({Feature::SchedStrategy_LPS, Feature::SchedStrategy_PAS, Feature::SchedStrategy_OI,
           Feature::SchedStrategy_POLLING});
    // WRROBIN, JSQ and SQ are still refused by NetworkStruct::refresh_routing,
    // so declaring them cannot let a model through. RROBIN is NOT: the refresh
    // now expands it uniformly for QNA/MNA, and this port's generator carries no
    // round-robin pointer (nvars reserves the column, nothing writes it), so a
    // declared RROBIN would answer a RANDOM-routing model under a round-robin
    // name. Undeclared until the pointer is resolved in after_event_router.
    f.set({Feature::RoutingStrategy_WRROBIN, Feature::RoutingStrategy_JSQ,
           Feature::RoutingStrategy_SQ});
    // SDR is the ONE state-dependent strategy this generator evaluates for real:
    // `rt_state` rebuilds the routing at every state from eq. (10) of Krzesinski
    // (1987), so the chain carries the routing the model declares rather than the
    // uniform placeholder `refresh_routing` leaves in `sn.rt`.
    f.set(Feature::RoutingStrategy_SDR);
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND});
    // ROUND-ROBIN DISPATCH, now that the pointer is a coordinate of the state:
    // `refresh_local_vars` allocates it, `append_local_vars` enumerates it,
    // `after_event_station`/`after_event_router` advance it on every departure,
    // and the generator reads the destination out of the ACTIVE node's
    // post-departure row. Before that the uniform expansion was all there was,
    // and declaring it would have answered a random-routing model.
    f.set({Feature::RoutingStrategy_RROBIN, Feature::RoutingStrategy_WRROBIN});
    f.set({Feature::ReplacementStrategy_RR, Feature::ReplacementStrategy_FIFO,
           Feature::ReplacementStrategy_SFIFO, Feature::ReplacementStrategy_LRU});
    f.set({Feature::ReplacementStrategy_HLRU, Feature::ReplacementStrategy_CLIMB,
           Feature::ReplacementStrategy_QLRU});
    f.set({Feature::ClosedClass, Feature::SelfLoopingClass, Feature::OpenClass, Feature::Replayer});
    f.set({Feature::OpenSignal, Feature::ClosedSignal});
    f.set({Feature::SignalType_NEGATIVE, Feature::SignalType_CATASTROPHE,
           Feature::SignalType_REPLY});
    f.set({Feature::SignalBatchRemoval, Feature::SignalRemovalPolicy});
    f.set({Feature::Place, Feature::Transition, Feature::Linkage, Feature::Enabling,
           Feature::Inhibiting, Feature::Timing, Feature::Firing, Feature::Storage});
    f.set({Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner});
    f.set(Feature::Retrial);
    // INERT: balking and reneging are call parameters of state_events.h, not
    // struct fields, so used_lang_features cannot emit them and declaring them
    // costs nothing.
    //
    // BREAKDOWN IS DIFFERENT SINCE 2026-08-15, when `sn.breakdownparam` gave it
    // a field and `used_lang_features` started emitting it. The declaration
    // stays, because MATLAB's SolverCTMC genuinely implements breakdowns
    // (`State.afterEventStation`'s FAILURE and REPAIR branches, on the shared
    // trailing marker column) and this set is that declaration transcribed. What
    // keeps the declaration from becoming a silent wrong answer is a REFUSAL BY
    // NAME in `solver_ctmc_analyzer.h`: this port's CTMC does not read the
    // field, so it names the gap instead of solving a model whose server never
    // fails. Delete that refusal in the same change that teaches solver_ctmc.h
    // the marker column, not before.
    f.set({Feature::Balking, Feature::Reneging, Feature::Breakdown});
    f.set({Feature::LoadDependence, Feature::ClassDependence, Feature::JointDependence});
    // phi(n) over the FULL population matrix; solver_ctmc.h tabulates it per state
    f.set(Feature::GlobalDependence);
    // c-server stations and binding buffers are both State constructs
    // (`from_marginal` / `after_event_station`), served by the explicit
    // generator and withdrawn from `cftp` and `mdd` below. SolverCTMC.m:443-446.
    f.set({Feature::MultiServer, Feature::FiniteCapacity});

    // PER-METHOD DELTAS. Four of the six methods enumerate the generator above
    // and share the whole envelope; `cftp` and `mdd` build no generator at all,
    // so what the rest of the envelope describes does not apply to them. The
    // rules the registry has no name for -- the class count that both need, the
    // station count `cftp` needs and the phase count it needs -- are structural
    // and live in `ctmc::solver_ctmc_cftp_supports` / `ctmc::solver_ctmc_mdd_supports`,
    // which the analyzers themselves call.
    if (method == "cftp" || method == "cftp.approx" || method == "cftp.exact") {
        // PERFECT SAMPLING FROM A BALANCE FUNCTION, not from a generator: the
        // sampler encodes the closed single-class product form of Gordon-Newell
        // and nothing else, so every construct outside it leaves the envelope.
        //
        // The one-phase-per-station rule is deliberately NOT spelled as a list
        // of distribution names: the sampler refuses phases(i,0) > 1, and a name
        // is not a phase count -- a one-phase Coxian passes and a HyperExp does
        // not, while Det/Gamma/Pareto only acquire their phases in
        // sn_nonmarkov_toph.
        f.unset(Feature::OpenClass);
        // Queue, Delay and Router are the only node kinds the sampler walks.
        f.unset(Feature::Source);
        f.unset(Feature::Sink);
        f.unset(Feature::RandomSource);
        f.unset(Feature::JobSink);
        f.unset(Feature::ClassSwitch);
        f.unset(Feature::StatelessClassSwitcher);
        f.unset(Feature::Cache);
        f.unset(Feature::CacheClassSwitcher);
        f.unset(Feature::CacheRetrieval);
        f.unset(Feature::ReplacementStrategy_RR);
        f.unset(Feature::ReplacementStrategy_FIFO);
        f.unset(Feature::ReplacementStrategy_SFIFO);
        f.unset(Feature::ReplacementStrategy_LRU);
        f.unset(Feature::ReplacementStrategy_HLRU);
        f.unset(Feature::ReplacementStrategy_CLIMB);
        f.unset(Feature::ReplacementStrategy_QLRU);
        f.unset(Feature::Fork);
        f.unset(Feature::Join);
        f.unset(Feature::Forker);
        f.unset(Feature::Joiner);
        f.unset(Feature::Place);
        f.unset(Feature::Transition);
        f.unset(Feature::Linkage);
        f.unset(Feature::Enabling);
        f.unset(Feature::Inhibiting);
        f.unset(Feature::Timing);
        f.unset(Feature::Firing);
        f.unset(Feature::Storage);
        // Disciplines outside INF/PS/FCFS/SIRO/LCFSPR have no product form.
        f.unset(Feature::SchedStrategy_DPS);
        f.unset(Feature::SchedStrategy_GPS);
        f.unset(Feature::SchedStrategy_SEPT);
        f.unset(Feature::SchedStrategy_LEPT);
        f.unset(Feature::SchedStrategy_HOL);
        f.unset(Feature::SchedStrategy_LCFS);
        f.unset(Feature::SchedStrategy_LCFSPRPRIO);
        f.unset(Feature::SchedStrategy_FCFSPRPRIO);
        f.unset(Feature::SchedStrategy_FCFSPR);
        f.unset(Feature::SchedStrategy_LCFSPI);
        f.unset(Feature::SchedStrategy_FCFSPI);
        f.unset(Feature::SchedStrategy_LCFSPIPRIO);
        f.unset(Feature::SchedStrategy_FCFSPIPRIO);
        f.unset(Feature::SchedStrategy_PSPRIO);
        f.unset(Feature::SchedStrategy_DPSPRIO);
        f.unset(Feature::SchedStrategy_GPSPRIO);
        f.unset(Feature::SchedStrategy_LPS);
        f.unset(Feature::SchedStrategy_PAS);
        f.unset(Feature::SchedStrategy_OI);
        f.unset(Feature::SchedStrategy_POLLING);
        f.unset(Feature::Region);
        f.unset(Feature::LoadDependence);
        f.unset(Feature::ClassDependence);
        f.unset(Feature::JointDependence);
        f.unset(Feature::GlobalDependence);
        // A state-dependent decision is not Markovian routing.
        f.unset(Feature::RoutingStrategy_RROBIN);
        f.unset(Feature::RoutingStrategy_WRROBIN);
        f.unset(Feature::RoutingStrategy_JSQ);
        f.unset(Feature::RoutingStrategy_SQ);
        f.unset(Feature::RoutingStrategy_SDR);
        // The Gordon-Newell balance function has no buffer; `solver_ctmc_cftp_supports`
        // refuses a finite one by name. The server count is NOT withdrawn: the
        // reference keeps it, the balance function carrying min(n,c).
        f.unset(Feature::FiniteCapacity);
    } else if (method == "mdd") {
        // The decision diagram holds the MARKING of a closed network; an open
        // stream makes it unbounded, so there is no finite diagram to hold. A
        // stochastic Petri net keeps the Place/Transition names: `spn_mdd` reads
        // the marking directly and is exempt from the closed-network rule.
        //
        // A FORK-JOIN MODEL IS NEITHER of the two shapes it serves: the tag
        // augmentation a fork needs adds one auxiliary class per branch, so the
        // struct that reaches the analyzer is never single-class however the
        // model was written, and the level decomposition has no meaning for a
        // firing that does not conserve the per-chain population.
        f.unset(Feature::OpenClass);
        f.unset(Feature::Source);
        f.unset(Feature::Sink);
        f.unset(Feature::RandomSource);
        f.unset(Feature::JobSink);
        f.unset(Feature::Fork);
        f.unset(Feature::Join);
        f.unset(Feature::Forker);
        f.unset(Feature::Joiner);
        f.unset(Feature::JoinPartial);
        // The level decomposition reads rates, servers and phases and no
        // sn.cap/classcap, so a buffer would be dropped. SolverCTMC.m:307-309.
        f.unset(Feature::FiniteCapacity);
    }
    return f;
}

/**
 * SolverSSA.getFeatureSet, 98 MATLAB names.
 *
 * IT IS THE SOLVER'S REACH, NOT ONE ENGINE'S, and that is what makes it usable:
 * `solver_ssa` under `default` runs the NRM when the NRM can run the model and
 * the event-driven serial engine when it cannot, so a construct EITHER engine
 * covers is a construct SolverSSA covers. AUTO chooses this solver on the
 * strength of this set (`solver_auto.h:453`), so under-declaring costs a model
 * the only solver that can answer it and over-declaring hands AUTO a refusal.
 *
 * The serial engine applies the SAME event handlers `state_events.h` gives the
 * CTMC generator, which is why the two sets are now close: a cache access, an
 * SPN firing, a fork firing, a polling controller and the priority shares are
 * all reached through handlers both solvers share. What is DECLARED here is
 * what a test in `test_ssa_fj_fcr.cpp` or `test_ssa_serial.cpp` exercises end
 * to end, never what a handler merely appears to support.
 *
 * Still dropped, with the C++ refusal that proves each:
 *   RoutingStrategy_RROBIN, _WRROBIN,
 *     _JSQ, _SQ, _RL                            ssa_dispatch.h refuses them in the NRM and
 *                                               `NetworkStruct::refresh_routing` refuses the
 *                                               state-dependent ones when the struct is built,
 *                                               so no engine here can receive one
 *   SignalType_REPLY                            the REPLY signal is a LAYERED construct whose
 *                                               guard sits in the model layer; SolverCTMC
 *                                               declares it, this solver has no test for it
 * `ssa_check_phases` (ssa_dispatch.h) stays imperative: non-exponential service
 * AT a given discipline is a per-(station,class) join of two features and the
 * boolean registry cannot state it. It gates the NRM only -- the serial engine
 * expands a phase-type process wherever the CTMC generator does -- so a model it
 * refuses is answered by the fallback rather than refused outright.
 */
inline FeatureSet ssa_feature_set(const std::string& /*method*/) {
    FeatureSet f;
    // Router is INERT here for the same reason as in ctmc_feature_set: nothing
    // emits it. SolverSSA.m:71 declares it all the same, so it stays.
    f.set({Feature::Sink, Feature::Source, Feature::Router});
    f.set({Feature::ClassSwitch, Feature::Delay, Feature::DelayStation, Feature::Queue});
    f.set({Feature::MAP, Feature::MMPP2, Feature::MMAP, Feature::APH, Feature::PH,
           Feature::Replayer});
    f.set({Feature::Coxian, Feature::Cox2, Feature::Erlang, Feature::Exp, Feature::HyperExp});
    f.set({Feature::Det, Feature::Gamma, Feature::Lognormal, Feature::Pareto, Feature::Uniform,
           Feature::Weibull});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer});
    f.set({Feature::SharedServer, Feature::Buffer, Feature::Dispatcher});
    f.set({Feature::Server, Feature::JobSink, Feature::RandomSource, Feature::ServiceTunnel});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS});
    f.set({Feature::SchedStrategy_DPS, Feature::SchedStrategy_FCFS});
    f.set({Feature::SchedStrategy_GPS, Feature::SchedStrategy_LPS, Feature::SchedStrategy_SIRO});
    f.set({Feature::SchedStrategy_HOL, Feature::SchedStrategy_LCFS});
    f.set({Feature::SchedStrategy_SEPT, Feature::SchedStrategy_LEPT});
    f.set(Feature::SchedStrategy_LCFSPR);
    // The priority variants the reference declares, reached through the serial
    // engine: `after_event_station_arv` and `_dep` carry the [class, phase]
    // buffer arm, and `test_ssa_serial.cpp` runs them against the exact CTMC.
    //
    // DELIBERATELY NOT THE REST OF THE PREEMPT FAMILY (LCFSPI, FCFSPR, FCFSPI,
    // LCFSPRIO, LCFSPIPRIO, FCFSPIPRIO). The C++ serial engine runs all of them
    // -- the same test file proves it -- but `SolverSSA.getFeatureSet` does not
    // declare them in MATLAB, and this set is what AUTO selects on. Declaring
    // more here would make the C++ AUTO pick SolverSSA for models the MATLAB
    // AUTO sends elsewhere, which is a cross-codebase divergence introduced by a
    // DECLARATION rather than by an algorithm. The reference's own list is
    // arguably short -- its serial engine goes through the same
    // `State.afterEvent` handlers -- and that question belongs on the MATLAB
    // side, recorded in `_kb/06-solver-catalog.md`, not settled by diverging.
    f.set({Feature::SchedStrategy_LCFSPRPRIO, Feature::SchedStrategy_FCFSPRPRIO});
    f.set({Feature::SchedStrategy_PSPRIO, Feature::SchedStrategy_DPSPRIO,
           Feature::SchedStrategy_GPSPRIO});
    // PAS / OI and the polling controller, reached through the serial engine:
    // the NRM refuses all three by name and `default` falls back. PAS and OI
    // needed `to_marginal`'s ordered-list arm (ported 2026-07-31) before the
    // queue length at such a station was a number at all.
    f.set({Feature::SchedStrategy_PAS, Feature::SchedStrategy_OI,
           Feature::SchedStrategy_POLLING});
    // Fork-join, on the TAG-AUGMENTED copy the analyzer builds; the sibling
    // classes are folded back before the table is returned.
    f.set({Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner});
    // Finite capacity regions: DROP censors the path exactly as it censors the
    // CTMC's space, and WAITQ runs the shared token-FIFO relation.
    f.set(Feature::Region);
    // A cache access is a state-dependent class switch plus a replacement, both
    // in `after_event_cache`; the analyzer writes the realized hit ratio back.
    f.set({Feature::Cache, Feature::CacheClassSwitcher, Feature::CacheRetrieval});
    // The stochastic Petri net path: a firing is a GLOBAL synchronization, which
    // the serial engine walks alongside the ordinary ones.
    f.set({Feature::Place, Feature::Transition, Feature::Linkage, Feature::Enabling,
           Feature::Inhibiting, Feature::Timing, Feature::Firing, Feature::Storage});
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND});
    // SDR, through the serial engine: the NRM refuses every state-dependent
    // strategy by name and `default` falls back, and the serial engine reads the
    // same per-state table (`rt_state`) the CTMC generator does.
    f.set(Feature::RoutingStrategy_SDR);
    // THE DISPATCHERS, both engines. The NRM walks the outgoing arcs itself and
    // keeps its own cursor (`resolve_state_dependent_dest`); the SERIAL engine
    // reads the pointer the state now carries, the same one the CTMC generator
    // reads. JSQ and SQ are NRM-only in their exact form -- the serial engine
    // consumes the uniform expansion for them -- which is the reference's own
    // position: `getRoutingMatrix.m:117` spreads both uniformly.
    f.set({Feature::RoutingStrategy_RROBIN, Feature::RoutingStrategy_WRROBIN,
           Feature::RoutingStrategy_JSQ, Feature::RoutingStrategy_SQ});
    f.set({Feature::ReplacementStrategy_RR, Feature::ReplacementStrategy_FIFO,
           Feature::ReplacementStrategy_SFIFO, Feature::ReplacementStrategy_LRU});
    f.set({Feature::ReplacementStrategy_HLRU, Feature::ReplacementStrategy_CLIMB,
           Feature::ReplacementStrategy_QLRU});
    f.set({Feature::SchedStrategy_EXT, Feature::ClosedClass, Feature::SelfLoopingClass,
           Feature::OpenClass});
    f.set({Feature::OpenSignal, Feature::ClosedSignal});
    f.set({Feature::SignalType_NEGATIVE, Feature::SignalType_CATASTROPHE});
    f.set({Feature::SignalBatchRemoval, Feature::SignalRemovalPolicy});
    f.set(Feature::Retrial);
    // INERT: balking and reneging are call parameters of state_events.h, not struct fields.
    f.set({Feature::Balking, Feature::Reneging});
    f.set(Feature::LoadDependence);
    // Class- and joint-dependent scaling: the NRM refuses both because it builds
    // its rates without evaluating the handle, and the serial engine evaluates
    // it through `state_events.h` exactly as the CTMC does.
    f.set({Feature::ClassDependence, Feature::JointDependence});
    // A global (Whittle) dependence is likewise refused by the NRM and carried by
    // the serial engine, which evaluates phi(n) once per sample-path state.
    f.set(Feature::GlobalDependence);
    // c-server stations and binding buffers: the serial engine walks the same
    // State arms ctmc_feature_set declares and the NRM honours both.
    // SolverSSA.m:178-181.
    f.set({Feature::MultiServer, Feature::FiniteCapacity});
    return f;
}

/**
 * SolverFLD.getFeatureSet, transcribed, MINUS what the requested method cannot
 * evaluate -- the port of `@@SolverFLD/getMethodFeatureSet`.
 *
 * THE PER-METHOD SUBTRACTION IS NOT COSMETIC. A station whose discipline has no
 * branch in `ode_rates_closing_factors` keeps g = x, i.e. it is integrated as an
 * INFINITE SERVER, and the answer is wrong with no warning at all: on
 * Delay(Z=1) -> Queue(c=1), N=4, whose exact Q2 is 3.0154, the fall-through
 * returns 2.0000. The closing family therefore has to reject the disciplines it
 * lacks a branch for rather than accept them into the wrong drift. Only the
 * closing family is affected: matrix/pnorm build a PS drift, which is the right
 * aggregate for any work-conserving discipline.
 *
 * `CacheRetrieval` is DELIBERATELY absent upstream (SolverFLD.m:95-101 records
 * the flow-conservation failure on examples/basic/cacheModel/retrieval_simple).
 * `NHPP`, `MAPt` and `PHt` are the time-inhomogeneous families. Every method
 * accepts them, as the reference does: the first-order methods integrate the
 * time-averaged NOMINAL pair, `kp` integrates the schedule itself, and a
 * time-varying rate multiplier is refused on the options rather than here.
 * What the C++ fluid
 * solver ignores silently (Region, Fork/Join, Router, Retrial, Place/Transition,
 * the signals, ClassDependence) is already absent from the MATLAB set.
 */
inline FeatureSet fluid_feature_set(const std::string& method) {
    FeatureSet f;
    f.set({Feature::ClassSwitch, Feature::Delay, Feature::DelayStation, Feature::Queue});
    f.set({Feature::Cache, Feature::CacheClassSwitcher});
    f.set({Feature::Cox2, Feature::Coxian, Feature::Erlang, Feature::Exp, Feature::HyperExp});
    f.set({Feature::APH, Feature::Det, Feature::MAP, Feature::MMPP2, Feature::NHPP});
    // The time-inhomogeneous families, which only `kp` integrates; the gate is
    // per-method below, so declaring them here is not a claim that every fluid
    // method reads a schedule.
    f.set({Feature::MAPt, Feature::PHt});
    // Non-Markovian renewal laws reach the ODE as acyclic PH via sn_nonmarkov_toph.
    f.set({Feature::Gamma, Feature::Lognormal, Feature::Pareto, Feature::Uniform,
           Feature::Weibull});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer, Feature::SharedServer,
           Feature::Buffer, Feature::Dispatcher});
    f.set({Feature::Server, Feature::ServiceTunnel});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS});
    f.set({Feature::SchedStrategy_DPS, Feature::SchedStrategy_FCFS, Feature::SchedStrategy_HOL});
    f.set({Feature::SchedStrategy_GPS});
    f.set({Feature::SchedStrategy_SIRO, Feature::SchedStrategy_LCFS,
           Feature::SchedStrategy_LCFSPR});
    f.set(Feature::LoadDependence);
    f.set({Feature::ReplacementStrategy_RR, Feature::ReplacementStrategy_FIFO});
    f.set(Feature::ReplacementStrategy_SFIFO);
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND});
    f.set({Feature::ClosedClass, Feature::SelfLoopingClass, Feature::Replayer});
    f.set({Feature::RandomSource, Feature::Sink, Feature::Source, Feature::OpenClass,
           Feature::JobSink});
    // Fork-join through the MMT (or Heidelberger-Trivedi) transform, driven by
    // `mva::fj_fixed_point` with a fluid inner solve; see fluid_runner.h. The
    // transform emits only Source, Delay, Queue, Router and ClassSwitch, all of
    // which the drift already carries, so nothing fork-join specific enters the
    // ODE. The per-method gate below removes it again where the fixed point
    // cannot run.
    f.set({Feature::Fork, Feature::Forker, Feature::Join, Feature::Joiner,
           Feature::JoinPartial});
    // c-server stations: the drifts carry min(n,c); withdrawn from 'diffusion'
    // and 'mfq' below. SolverFLD.m:605-607.
    f.set(Feature::MultiServer);
    // A binding buffer, SolverFLD.m:608-614: 'dae' carries it as an algebraic
    // constraint, 'mol' IS the Mt/G/s/0 loss system and the AoI arm of 'mfq' is
    // a bufferless or single-buffer queue. WHICH method serves one is the
    // structural rule `fluid_runner`'s own `check_binding_capacity` call asks,
    // so no per-method delta duplicates it here; the reference withdraws it from
    // no method for exactly that reason.
    f.set(Feature::FiniteCapacity);

    // THE CANONICAL NAME, which every test below is on. Stripping only the
    // `fluid.` prefix here left two holes: `butools` and `aoi` are aliases of
    // `mfq` (fluid_runner.h `fluid_unqualify`) and lost the HOL grant `mfq`
    // keeps, reachable through AUTO, which passes the token unmapped; and the
    // short spellings `ggisgi`/`tga` are aliases of the two single-station
    // limits (fluid_qsys.h `fluid_qsys_canonical`). Those two functions cannot be
    // called from here without a cycle (both headers read this one), so the three
    // rules are mirrored inline and must move together with them.
    std::string m =
        (method.size() > 6 && method.compare(0, 6, "fluid.") == 0) ? method.substr(6) : method;
    if (m == "butools" || m == "aoi") m = "mfq";
    else if (m == "ggisgi") m = "ggisgi.fluid";
    else if (m == "tga") m = "ggingi.tga";
    const bool closing_family =
        m == "closing" || m == "statedep" || m == "softmin" || m == "tbi" || m == "minnormal" ||
        m == "refined" || m == "dae";
    // Limited load dependence composes with the closure as a rate multiplier
    // alpha(n_i) on the scheduling share, which only the closing family evaluates.
    // The matrix, pnorm, softmin, statedep, tbi, diffusion, mfq, kp and rmf paths
    // build their drift independently and would silently ignore alpha.
    if (!(m == "closing" || m == "minnormal" || m == "refined" || m == "dae"))
        f.unset(Feature::LoadDependence);
    // SIRO, LCFS and LCFSPR have no branch in the closing drift. The metric reader
    // of `solver_fluid_closing` caught LCFS/LCFSPR/HOL, but `solver_fluid_moments`
    // reads the event representation directly and never reached that guard, and
    // SIRO is worse still: the reader accepts it AS FCFS, so the ODE integrated it
    // as INF while the metrics were read as if it shared the server.
    if (closing_family) {
        f.unset(Feature::SchedStrategy_SIRO);
        f.unset(Feature::SchedStrategy_LCFS);
        f.unset(Feature::SchedStrategy_LCFSPR);
    }
    // GPS splits the server by weight among the BACKLOGGED classes, so its share
    // is a function of the backlog INDICATOR. A first-order closure cannot express
    // it at all: with continuous x_k > 0 every class is always backlogged and the
    // share collapses to the constant w_k/sum_j w_j, the heavy-traffic limit,
    // regardless of load. Only `minnormal` supplies the P(X_k >= 1) it needs --
    // `refined` reads its rates at the MEAN-FIELD variance, which is exactly the
    // degenerate case.
    if (m != "minnormal") f.unset(Feature::SchedStrategy_GPS);
    // HOL allocates capacity in PRIORITY order, not in proportion to population,
    // and no fluid drift reads the class priorities except the single-queue `mfq`
    // priority branch.
    if (m != "mfq") f.unset(Feature::SchedStrategy_HOL);
    // MULTISERVER, SolverFLD.m:387-395: the drifts carry min(n,c) except
    // `fluid_diffusion`, whose SDE is written for one or infinitely many
    // servers, and `solver_mfq`, a single-queue model on the same server counts
    // (off them `mfq` resolves to `matrix`, so the delta binds only where it
    // runs as itself). fluid_diffusion.h:89-93 keeps wording that refusal.
    if (m == "diffusion" || m == "mfq") f.unset(Feature::MultiServer);
    // `dae` is the one method that WIDENS this set, and the widening is why the
    // gate has to be per-method: a Finite Capacity Region is a linear inequality
    // on the state, which the DAE form carries as an algebraic equation beside
    // the drift and no ODE method can carry at all. Every other fluid method has
    // to go on refusing it -- an ODE integrated through a cap it cannot see
    // returns the UNCONSTRAINED answer with no warning -- so the widening is
    // stated here rather than in the static set. `fluid_dae_constraints` still
    // refuses by name the region forms that are not constraints on this drift
    // (drop, BAS/BBS/RSRD, retrial, per-class admission weights), and
    // `solver_fluid_dae` refuses the FCR transient, which needs event restarts.
    //
    // DPS goes the other way: it closes on the covariance BETWEEN a station's
    // class coordinates, where the DAE form carries one scalar variance per
    // station. `minnormal` carries those matrix blocks through its outer
    // iteration; adding an unknown block per station here would restore the
    // quartic cost that keeping Sigma out of the Newton vector avoids. GPS is
    // already unset above for everything but `minnormal`, for the same reason.
    if (m == "dae") {
        f.set(Feature::Region);
        f.unset(Feature::SchedStrategy_DPS);
    }
    // A STOCHASTIC PETRI NET HAS NO DRIFT OUTSIDE THE DAE FORM: its conserved
    // quantities are P-invariants rather than chain populations, an immediate
    // transition is an algebraic FLOW rather than an event with a rate, and a
    // bounded place is a linear inequality on the marking (fluid_petri.h). Every
    // other fluid method builds its drift from the station/class/phase encoding,
    // where a Place contributes no coordinate at all, so it would integrate the
    // net as an empty model and report zeros without a warning -- which is why
    // the names are declared for `dae` alone rather than in the static set.
    // `dae` ALONE, not `default`: the reference's own gate removes them for every
    // other method including the unresolved default (SolverFLD.m:240), and the
    // Petri route is reached by name or by the runner's own branch, never by a
    // default that has not yet resolved.
    if (m == "dae")
        f.set({Feature::Place, Feature::Transition, Feature::Linkage, Feature::Enabling,
               Feature::Inhibiting, Feature::Timing, Feature::Firing, Feature::Storage});
    // THE MMT FIXED POINT, and which methods can actually run it. This list used
    // to hold seven names on the argument that the transform hands the inner
    // solve a MIXED network -- the parallelism rides on auxiliary OPEN classes
    // even when the model is closed -- so a method that takes only closed
    // models, only open ones, a single queue, or a decomposition could not
    // drive it. MEASURED on a SYMMETRIC closed fork-join (Delay -> Fork -> two
    // identical FCFS queues -> Join, N = 2), against the exact chain
    // Q1 = Q2 = 0.664, J = 0.624, D = 1.024, five of the seven answer it and
    // answer it symmetrically, which is the one property no approximation of a
    // symmetric model may lose:
    //
    //   statedep, mfq, rmf   0.5714 0.5714 0.5714 1.1429  (the `matrix` drift)
    //   refined              0.2715 0.2715 0.9713 0.4857
    //   tbi                  0.2500 0.2500 1.0000 0.5000
    //
    // So the argument does not hold for those five and the names come back;
    // MATLAB and native python have always run them there. It DOES hold for the
    // other two, which is why they stay:
    //
    //   diffusion   0 0 2 0        the whole population on ONE station and zero
    //               (and 0 1.997 0 0.003 on a rerun -- a DIFFERENT station), so
    //               the SDE is not integrating this model at all
    //   kp          0 0 0 0        on a symmetric OPEN fork-join fed at rate
    //               0.5, i.e. an empty network where jobs are arriving
    //
    // Both are silent wrong answers rather than refusals, so keeping the names
    // off is what makes the model refused instead of mis-answered. MATLAB and
    // native python are over-permissive here and now drop them too.
    if (m == "diffusion" || m == "kp") {
        f.unset(Feature::Fork);
        f.unset(Feature::Forker);
        f.unset(Feature::Join);
        f.unset(Feature::Joiner);
        f.unset(Feature::JoinPartial);
    }
    // MAPt/PHt/NHPP are NOT removed per method, and that is the reference's own
    // position rather than an oversight. `solver_fluid.m:20-35` and
    // `solver_fluid_matrix.m:18-27` install the width-weighted NOMINAL pair for a
    // schedule and integrate it deliberately: the nominal is the stationary
    // carrier of the phase structure the schedule modulates, so a steady-state
    // answer at the time-averaged rate is the intended approximation, not a
    // silent substitution. `kp` is the method that integrates the SCHEDULE, and
    // `getMethodFeatureSet` in the reference restricts none of the three.
    //
    // The refusal that does exist lives on the OPTIONS, not here: a time-varying
    // rate multiplier (`nhpp_sched` / `rate_traj` / `rate_sched`) makes the drift
    // non-autonomous, which the moment closures cannot take, and
    // `fluid_minnormal_applicable` declines `minnormal` for it while
    // `fluid_moment_terms` raises. Those are conditions on the OPTIONS a caller
    // passes and cannot be decided from the model alone, which is why a featset
    // is the wrong place for them.
    // ON THE CANONICAL NAME, not the raw one. These two branches used to compare
    // `method`, so a direct call with a `fluid.`-qualified spelling silently got
    // a DIFFERENT envelope than the bare one; it was masked on the solver path
    // only because `fluid_resolve_method` unqualifies on entry. test_feature_set
    // asserts the invariant for `minnormal` alone, and the same assertion written
    // for `fluid.tvms` would have failed.
    if (m == "ggisgi.fluid" || m == "ggingi.tga" || m == "tvms") {
        // The only fluid methods in LINE stated for a queue customers ABANDON.
        // Reneging stays out of the base FLD envelope: the network drift carries
        // no abandonment flow, so every other method would integrate the model
        // as if nobody left.
        f.set(Feature::Reneging);
    }
    if (m == "ggisgi.fluid" || m == "ggingi.tga" || m == "tvms" || m == "mtginf" ||
        m == "mol") {
        // Every one of them is stated for a single open station; the base
        // envelope's closed classes have no meaning there.
        f.unset(Feature::ClosedClass);
        f.unset(Feature::SelfLoopingClass);
    }
    // CLOSED MODELS ONLY for `refined`, which the reference's runAnalyzer has
    // always enforced by name and this set never stated: the 1/N correction is
    // solved on orth(D) over the FULL state, so on an open model it adds a
    // perturbation to the SOURCE POOL mass, a normalisation constant rather than
    // a population. Only `minnormal` was validated open. Stating it here is what
    // lets a report withdraw the pair instead of offering a run that stops.
    if (m == "refined") {
        f.unset(Feature::OpenClass);
        f.unset(Feature::Source);
        f.unset(Feature::Sink);
        f.unset(Feature::RandomSource);
        f.unset(Feature::JobSink);
    }
    // The diffusion SDE PROJECTS each class back onto its own fixed population
    // at every step, which is the closed-network constraint itself: an open
    // class has no population to project onto, and a Source is not a station the
    // SDE has a coordinate for. `fluid_diffusion` refuses both by name, and
    // stating them here is what lets a CALLER see the refusal before running.
    if (m == "diffusion") {
        f.unset(Feature::OpenClass);
        f.unset(Feature::Source);
        f.unset(Feature::Sink);
        f.unset(Feature::RandomSource);
        f.unset(Feature::JobSink);
    }
    // The Ko-Pender limits are proved for an OPEN network of stations fed by
    // external arrival processes: a closed class has no arrival process to
    // modulate and no source phase to carry, and the cache and class-switch
    // machinery has no counterpart in the paper's event set. Only the fork-join
    // names were withheld here, five features against the nine MATLAB, the JAR and
    // native python all unset, so this port alone offered `kp` on closed, cache
    // and class-switching models the other three refuse -- and `doc/solver-methods.md`
    // states the narrowing as the documented behaviour.
    if (m == "kp") {
        f.unset(Feature::ClosedClass);
        f.unset(Feature::SelfLoopingClass);
        f.unset(Feature::Cache);
        f.unset(Feature::CacheClassSwitcher);
        f.unset(Feature::ClassSwitch);
        f.unset(Feature::StatelessClassSwitcher);
        f.unset(Feature::ReplacementStrategy_RR);
        f.unset(Feature::ReplacementStrategy_FIFO);
        f.unset(Feature::ReplacementStrategy_SFIFO);
    }
    // Trajectory-based iteration decomposes the CLOSED population into cells and
    // relaxes the waveforms between them; there is no cell for an unbounded open
    // stream. A cache model is solved by decomposition rather than by one drift,
    // so the cell partition has nothing to partition -- use `rmf`.
    if (m == "tbi") {
        f.unset(Feature::OpenClass);
        f.unset(Feature::Source);
        f.unset(Feature::Sink);
        f.unset(Feature::RandomSource);
        f.unset(Feature::JobSink);
        f.unset(Feature::Cache);
        f.unset(Feature::CacheClassSwitcher);
        f.unset(Feature::ReplacementStrategy_RR);
        f.unset(Feature::ReplacementStrategy_FIFO);
        f.unset(Feature::ReplacementStrategy_SFIFO);
    }
    return f;
}

/**
 * SolverMAM.getFeatureSet, the union of its four setTrue calls: 55 MATLAB
 * names, WIDENED for 'default'/'ldqbd' only, and NARROWED by one name no MAM
 * path ever serves.
 *
 * THE METHOD PARAMETER IS NOW READ, where it used to be discarded. Two
 * consequences of discarding it were both the over-declaration
 * feature_set.h:30-33 warns against:
 *
 *   +LoadDependence, method in {default, ldqbd} ONLY
 *                     solver_mam_ldqbd.h:152 reads st.lldscaling, but only the
 *                     'ldqbd' method and 'default''s single-class closed
 *                     Delay+Queue route (mam_dispatch.h's branch 2e) ever call
 *                     it. Declaring it for every method let 'dec.source',
 *                     'dec.poisson' and 'mna' pass a load-dependent model
 *                     straight to solver_mam_basic / solver_mna, neither of
 *                     which reads the field, so those methods silently solved
 *                     at nominal rates. The featset cannot see topology, so
 *                     'default' still needs the structural check in
 *                     check_model_method (solver_mam_runner.h) that refuses a
 *                     load-dependent model whose shape is not the one ldqbd
 *                     can solve; 'ldqbd' itself is already self-checking
 *                     (solver_mam_ldqbd.h:109-131).
 *   -SchedStrategy_FCFSPRPRIO, every method
 *                     No MAM analyzer serves it: solve_fcfs_station in
 *                     solver_mam_basic.h dispatches only FCFS and HOL
 *                     (solver_mam_basic.h:1144), so an FCFSPRPRIO station
 *                     always reaches the station ladder's unconditional throw
 *                     at solver_mam_basic.h:1147-1154 regardless of this
 *                     declaration. The enumerator itself DOES exist in C++
 *                     (lang_types.h:144); what is absent is a serving path, so
 *                     the fix is here, not in SchedStrategy. Undeclaring lets
 *                     the gate refuse first, by feature name, instead of
 *                     reaching the deeper and less specific station-ladder
 *                     message.
 * DMAP, ME, RAP and SetupDelayOff stay declared and are inert: no ProcessType
 * or NetworkStruct field can make used_lang_features emit them.
 */
inline FeatureSet mam_feature_set(const std::string& method) {
    FeatureSet f;
    f.set({Feature::Sink, Feature::Source});
    f.set({Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner});
    f.set({Feature::Delay, Feature::DelayStation, Feature::Queue});
    // ME and RAP are INERT: ProcessType carries no such enumerator.
    f.set({Feature::APH, Feature::Coxian, Feature::Cox2, Feature::Erlang, Feature::Exp,
           Feature::HyperExp, Feature::MMPP2, Feature::MAP, Feature::MMAP, Feature::DMAP,
           Feature::ME, Feature::RAP});
    // The slotted path (mam_dispatch.h branch -1). A Geometric or a
    // DiscreteUniform reaches solver_mam_dt on the slot lattice; off it,
    // sn_is_discrete_time reports the model continuous and the continuous
    // analyzers below take the same laws through their moment fits.
    f.set({Feature::Geometric, Feature::DiscreteUniform});
    f.set({Feature::Det, Feature::Gamma, Feature::Lognormal, Feature::Pareto, Feature::Uniform,
           Feature::Weibull});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer});
    f.set(Feature::ClassSwitch);
    f.set({Feature::SharedServer, Feature::Buffer, Feature::Dispatcher});
    f.set({Feature::Server, Feature::JobSink, Feature::RandomSource, Feature::ServiceTunnel});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS, Feature::SchedStrategy_HOL});
    // FCFSPRPRIO is deliberately absent; see the header note above.
    f.set(Feature::SchedStrategy_FCFS);
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND});
    f.set({Feature::ClosedClass, Feature::SelfLoopingClass});
    f.set(Feature::OpenClass);
    // The G-network names are NOT here. They moved to ag_feature_set with the
    // RCAT methods that read `issignal`: no MAM algorithm reads the field, so
    // any of them would solve the model with the signals turned into ordinary
    // customers and report that as the answer.
    // Call 3, the BMAP/PH/N/N retrial queue.
    f.set({Feature::Retrial, Feature::BMAP, Feature::PH});
    // Call 4, qbd_setupdelayoff, reached through NetworkStruct::setupparam.
    f.set(Feature::SetupDelayOff);
    // c-server stations (every analyzer reads sn.nservers; the slotted path's
    // single-server rule stays structural, solver_mam_dt.h:73) and finite
    // buffers as LOSS buffers: solver_mam_basic's M/M/c/K and MMAP[K]/G/1/K,
    // solver_mam and solver_mna_open's truncate-and-renormalize, and
    // solver_mam_retrial's bufferless N/N station. SolverMAM.m:363-368.
    f.set({Feature::MultiServer, Feature::FiniteCapacity});
    // solver_mam_decmmap is an OPEN-network departure-process fixed point: it
    // iterates on arrival streams a closed population does not have, and its
    // opening station loop serves EXT, FCFS, HOL, FCFSPRPRIO and PS only (INF is
    // NOT in that list, so a model carrying a Delay belongs to 'dec.source').
    // Both restrictions are things the model HAS, so both belong in the featset;
    // solver_mam_decmmap.h raises the matching message when the method is named
    // by hand. Until this delta existed the gate offered dec.mmap on every closed
    // model, which the reference and the JAR answered with a table of zeros.
    // Fork-join goes with them: the sweep uses solver_mam_traffic, the PLAIN
    // traffic step, so a Fork or a Join is refused by the traffic step itself
    // rather than synchronized (solver_mam_decmmap.h says so in its header) --
    // which is why mam_dispatch routes a fork-join topology from
    // 'default'/'dec.source' to solver_mam_fj and never here.
    if (method == "dec.mmap")
        for (Feature ff : {Feature::ClosedClass, Feature::SelfLoopingClass,
                           Feature::SchedStrategy_INF, Feature::Fork, Feature::Join,
                           Feature::Forker, Feature::Joiner})
            f.unset(ff);
    // Widening over MATLAB, and only for the methods that actually read it.
    if (method == "default" || method == "ldqbd") f.set(Feature::LoadDependence);
    // Round-robin dispatching enters as a deterministic traffic split
    // (npfqn_traffic_split_rr), which only solver_mna_open resolves; SolverMAM.m
    // :60 sets it true for this method alone. The closed branch has no
    // counterpart and is refused structurally in check_model_method.
    if (method == "mna") f.set(Feature::RoutingStrategy_RROBIN);
    // FINITECAPACITY, SolverMAM.m:143-147: the two chains that read no sn.cap
    // withdraw it. The LD-QBD levels run to the population or the cutoff, and
    // the background chain to the state cap. A buffer a CLOSED class can fill is
    // refused structurally by `runner_detail::check_model_method` instead (a
    // closed job blocks, a loss formula does not).
    if (method == "ldqbd" || method == "bgchain") f.unset(Feature::FiniteCapacity);
    return f;
}

/**
 * SolverAG.getFeatureSet: what the RCAT decomposition can represent.
 *
 * THE G-NETWORK NAMES LIVE HERE AND NOWHERE ELSE. solver_ag.h's build_rcat is
 * the only code in LINE that reads `issignal`, so declaring the signal features
 * on any other solver would promise what that solver cannot deliver -- it would
 * answer with every signal turned into an ordinary customer.
 *
 * Every AG method is the same decomposition, differing only in how the reversed
 * rate is read off an agent, so the envelope does not vary by method; the
 * genuine restrictions (Markovian service law, single server) are structural and
 * are applied in ag::runner_detail::check_model_method.
 */
inline FeatureSet ag_feature_set(const std::string& /*method*/) {
    FeatureSet f;
    f.set({Feature::Sink, Feature::Source});
    f.set({Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner});
    f.set({Feature::Delay, Feature::DelayStation, Feature::Queue});
    f.set({Feature::APH, Feature::Coxian, Feature::Erlang, Feature::Exp, Feature::HyperExp,
           Feature::MAP, Feature::MMPP2});
    f.set({Feature::Det, Feature::Gamma, Feature::Lognormal, Feature::Pareto, Feature::Uniform,
           Feature::Weibull});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer});
    f.set({Feature::SharedServer, Feature::Buffer, Feature::Dispatcher});
    f.set({Feature::Server, Feature::JobSink, Feature::RandomSource, Feature::ServiceTunnel});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS, Feature::SchedStrategy_FCFS});
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND});
    // A self-looping class is its own single-station component (`solver_ag`
    // reads the flag), so it is declared beside the two ordinary class types.
    // MultiServer and FiniteCapacity are deliberately absent, as in
    // SolverAG.m:211-216: `check_model_method` words both refusals, RCAT driving
    // rho = lambda/mu with no buffer.
    f.set({Feature::ClosedClass, Feature::SelfLoopingClass, Feature::OpenClass});
    f.set({Feature::OpenSignal, Feature::ClosedSignal});
    f.set({Feature::SignalType_NEGATIVE, Feature::SignalType_CATASTROPHE});
    f.set(Feature::SignalBatchRemoval);
    return f;
}

/**
 * SolverENV.supports's inline set, 28 names, transcribed unchanged.
 *
 * The narrowest declared set in the codebase. `RoutingStrategy_RROBIN` is
 * marked "with SolverJMT" upstream and was kept here while the CTMC stage
 * solver declared it; it is now DROPPED, because ctmc_feature_set no longer
 * does. env_dispatch.h:164 routes the statevec coupling to that same CTMC, so
 * declaring round robin here would promise what the stage cannot deliver.
 */
inline FeatureSet env_feature_set(const std::string& /*method*/) {
    FeatureSet f;
    f.set({Feature::ClassSwitch, Feature::Delay, Feature::DelayStation, Feature::Queue,
           Feature::Sink, Feature::Source});
    f.set({Feature::Coxian, Feature::Cox2, Feature::Erlang, Feature::Exp, Feature::HyperExp});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer, Feature::SharedServer,
           Feature::Buffer, Feature::Dispatcher, Feature::Server, Feature::JobSink,
           Feature::RandomSource, Feature::ServiceTunnel});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS, Feature::SchedStrategy_FCFS});
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND});
    f.set({Feature::ClosedClass, Feature::OpenClass});
    return f;
}

/**
 * SolverLDES.getFeatureSet, transcribed WHOLE.
 *
 * THE ONE SET IN THIS FILE THAT IS NOT NARROWED, and for a reason that does not
 * generalise: every other function here declares what the C++ code implements,
 * because a feature declared and not implemented yields a wrong number. The
 * LDES client implements NONE of these -- it forwards the model.json document to
 * the SSJ engine byte for byte and reads the result back -- so the code that
 * has to support a feature is the engine, the same engine MATLAB and native
 * Python drive. Narrowing here would refuse models that the reference simulates
 * with the very binary this port is about to run, which is the opposite of the
 * usual danger.
 *
 * SOURCE: `matlab/src/solvers/wrappers/LDES/@@SolverLDES/SolverLDES.m:108-230`,
 * name for name, and every name there has a registry enumerator. `BatchArrival`
 * IS one of them: the note that used to sit here called it a name no solver
 * declares and no recorder sets, and both halves were false -- SolverLDES.m:210
 * declares it and MNetwork.getUsedLangFeatures marks it, as `used_lang_features`
 * has here since the Source arm learned `arrival_batch`. Left undeclared, this
 * set refused every model built with `Source.setArrivalBatch`, which is the one
 * engine that simulates the batch.
 *
 * WHAT THE GATE STILL CANNOT SEE. The set is consulted against
 * `used_lang_features`, which is built from a PARSED `NetworkStruct`; the LDES
 * arm of the CLI never parses the document, so a model outside this port's
 * reader reaches the engine without passing here at all. That is intended --
 * forwarding losslessly is the whole point -- and it means this set governs
 * SolverAUTO's choice rather than the client's own admission.
 */
inline FeatureSet ldes_feature_set(const std::string& /*method*/) {
    FeatureSet f;
    f.set({Feature::Sink, Feature::Source, Feature::Queue, Feature::Delay});
    // JoinPartial: the engine fires the join at the k-th sibling and discards
    // the stragglers on arrival, so the quorum is an exact sample-path event
    f.set({Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner,
           Feature::JoinPartial});
    f.set({Feature::Place, Feature::Transition, Feature::QueueingPlace});
    f.set({Feature::Linkage, Feature::Enabling, Feature::Inhibiting, Feature::Timing,
           Feature::Firing, Feature::Storage});
    f.set({Feature::Logger, Feature::LogTunnel, Feature::Buffer, Feature::Region});
    f.set({Feature::Exp, Feature::Erlang, Feature::HyperExp, Feature::PH, Feature::APH,
           Feature::Coxian, Feature::Cox2});
    f.set({Feature::MAP, Feature::DMAP, Feature::MMAP, Feature::BMAP, Feature::MMPP2,
           Feature::ME, Feature::RAP});
    f.set({Feature::Immediate, Feature::Disabled, Feature::Replayer, Feature::Trace});
    f.set({Feature::Det, Feature::Uniform, Feature::Gamma, Feature::Pareto, Feature::Weibull,
           Feature::Lognormal});
    // Lattice-valued and counting laws: Geometric backs Geo/Geo/1 and the
    // slotted mode, and the zero atom of a counting law becomes an immediate
    // interval in continuous mode.
    f.set({Feature::Geometric, Feature::Bernoulli, Feature::Binomial, Feature::Poisson});
    // Source.setArrivalBatch: the model document carries `arrivalBatch` and the
    // engine releases the whole batch at one arrival epoch (SolverLDES.m:210).
    f.set(Feature::BatchArrival);
    // Time-inhomogeneous processes: the piecewise-constant schedule is simulated
    // exactly by carrying the phase across a breakpoint.
    f.set({Feature::NHPP, Feature::MAPt, Feature::PHt});
    f.set({Feature::Server, Feature::JobSink, Feature::RandomSource, Feature::InfiniteServer,
           Feature::SharedServer, Feature::ServiceTunnel, Feature::DelayStation});
    f.set({Feature::SchedStrategy_FCFS, Feature::SchedStrategy_INF, Feature::SchedStrategy_HOL,
           Feature::SchedStrategy_FCFSPRIO});
    f.set({Feature::SchedStrategy_PS, Feature::SchedStrategy_DPS, Feature::SchedStrategy_GPS});
    f.set({Feature::SchedStrategy_LCFS, Feature::SchedStrategy_LCFSPR,
           Feature::SchedStrategy_LCFSPI});
    f.set({Feature::SchedStrategy_FCFSPR, Feature::SchedStrategy_FCFSPI,
           Feature::SchedStrategy_LPS, Feature::SchedStrategy_SIRO});
    f.set({Feature::SchedStrategy_SJF, Feature::SchedStrategy_LJF, Feature::SchedStrategy_LEPT,
           Feature::SchedStrategy_SEPT});
    f.set({Feature::SchedStrategy_SRPT, Feature::SchedStrategy_SRPTPRIO,
           Feature::SchedStrategy_PSJF, Feature::SchedStrategy_FB, Feature::SchedStrategy_LRPT});
    f.set({Feature::SchedStrategy_EXT, Feature::SchedStrategy_POLLING});
    f.set({Feature::SchedStrategy_PSPRIO, Feature::SchedStrategy_DPSPRIO,
           Feature::SchedStrategy_GPSPRIO});
    f.set({Feature::SchedStrategy_LCFSPRIO, Feature::SchedStrategy_LCFSPRPRIO,
           Feature::SchedStrategy_LCFSPIPRIO});
    f.set({Feature::SchedStrategy_FCFSPRPRIO, Feature::SchedStrategy_FCFSPIPRIO});
    f.set({Feature::SchedStrategy_FSP, Feature::SchedStrategy_PAS, Feature::SchedStrategy_OI});
    f.set({Feature::SchedStrategy_EDD, Feature::SchedStrategy_EDF, Feature::SchedStrategy_SETF});
    f.set({Feature::Router, Feature::Dispatcher, Feature::ClassSwitch,
           Feature::StatelessClassSwitcher});
    f.set({Feature::Cache, Feature::CacheClassSwitcher, Feature::CacheRetrieval,
           Feature::CacheItemSize});
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND,
           Feature::RoutingStrategy_RROBIN, Feature::RoutingStrategy_WRROBIN,
           Feature::RoutingStrategy_JSQ, Feature::RoutingStrategy_SQ});
    // HeteroServers: each job occupies one server of one compatible pool and is
    // served at that pool's own law (`free_slot_for` / `start_service`). The
    // sharing disciplines are refused by name in the engine rather than
    // flattened, so the declaration does not over-promise.
    f.set(Feature::HeteroServers);
    // SDR: the engine re-evaluates Krzesinski eq. (10) at the live marking
    // (`draw_sdr` in ldes_engine.h), so the sample path has an exact
    // counterpart in solver_nc_sdr rather than the uniform split the static
    // routing row carries. Declared here as well as in the reference because
    // the native engine now honours it; before that it was the ONE name whose
    // absence mattered, since the engine answered the uniform model in
    // silence instead of refusing.
    f.set(Feature::RoutingStrategy_SDR);
    f.set({Feature::OpenClass, Feature::ClosedClass, Feature::SelfLoopingClass});
    f.set({Feature::OpenSignal, Feature::ClosedSignal, Feature::SignalType_NEGATIVE,
           Feature::SignalType_REPLY, Feature::SignalType_CATASTROPHE,
           Feature::SignalBatchRemoval, Feature::SignalRemovalPolicy});
    f.set({Feature::LoadDependence, Feature::ClassDependence, Feature::JointDependence,
           Feature::SetupDelayOff});
    f.set({Feature::Balking, Feature::Reneging, Feature::Retrial});
    // Breakdown is NOT inert here, unlike in the CTMC set above: the LDES
    // engine reads `sn.breakdownparam` and simulates the outage, so the
    // declaration is enforced by a struct field a caller can actually set.
    f.set(Feature::Breakdown);
    f.set({Feature::ReplacementStrategy_RR, Feature::ReplacementStrategy_FIFO,
           Feature::ReplacementStrategy_SFIFO, Feature::ReplacementStrategy_LRU,
           Feature::ReplacementStrategy_HLRU, Feature::ReplacementStrategy_CLIMB,
           Feature::ReplacementStrategy_QLRU});
    // c-server stations and a finite buffer with its drop rule: the engine holds
    // both on the sample path, so the simulated station is the one the model
    // declares. SolverLDES.m:230.
    f.set({Feature::MultiServer, Feature::FiniteCapacity});
    return f;
}

/**
 * Port of `SolverJMT.getFeatureSet` (`@@SolverJMT/SolverJMT.m:180-290`).
 *
 * JMT IS THE REFERENCE FCR ENGINE -- every other solver's finite-capacity
 * refusal points users here -- so `Region` is declared, and so is the whole SPN
 * section family, which JMT simulates natively.
 *
 * WHAT IS NARROWED AGAINST THE REFERENCE, and why each is a transcription of a
 * refusal this port already makes rather than a feature JMT lacks:
 *   Normal            the reference lists it, and `saveServiceStrategy` has no
 *                     branch for it: a Normal service would reach the analytic
 *                     switch with no `javaClass` and error. It is not declared.
 *   HLRU/CLIMB/QLRU   JMT has cache objects for RR, FIFO, SFIFO and LRU only;
 *                     `save_cache_strategy` refuses the other three by name
 *                     instead of falling back to LRU, so declaring them would
 *                     promise an admission the writer then denies.
 *   SchedStrategy_LPS is declared: the limit is exported as a single-node
 *                     blocking region, not dropped.
 */
/**
 * True for the JMVA algorithms that solve a CLOSED product-form network only.
 *
 * RECAL, CoMoM, Chow, Bard-Schweitzer (both spellings), AQL, Linearizer and De
 * Souza-Muntz Linearizer. Measured against JMT 1.2.x: each answers an open or a
 * mixed model with `jmt.common.exception.UnsupportedModelException: The selected
 * solver cannot handle open classes, please choose another.` and a
 * load-dependent one with the same exception naming load-dependent stations,
 * while the exact MVA engine behind 'jmva' and 'jmva.mva' serves both. The
 * single-server half of the same restriction is a station COUNT and has no
 * feature name, so it is refused structurally by `jmt::jmt_method_refusal`.
 */
inline bool jmva_is_closed_only(const std::string& method) {
    return method == "jmva.amva" || method == "jmva.recal" || method == "jmva.comom" ||
           method == "jmva.chow" || method == "jmva.bs" || method == "jmva.aql" ||
           method == "jmva.lin" || method == "jmva.dmlin";
}

inline FeatureSet jmt_feature_set(const std::string& method) {
    FeatureSet f;
    f.set({Feature::Sink, Feature::Source, Feature::Queue, Feature::Delay, Feature::DelayStation});
    f.set({Feature::Router, Feature::ClassSwitch, Feature::StatelessClassSwitcher,
           Feature::Dispatcher});
    // JoinPartial: written out as a jmt PartialJoin with numRequired = k
    f.set({Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner,
           Feature::JoinPartial});
    // A variable forking level: the writer turns isSimplifiedFork off and emits the
    // per-branch entries, so jmt reads the counts, the probabilities and the degree
    // distribution rather than sending one job down every link.
    f.set({Feature::ForkFanoutVector, Feature::ForkFanoutRandom,
           Feature::ForkBranchProbability});
    f.set({Feature::Logger, Feature::LogTunnel, Feature::Buffer, Feature::Region,
           Feature::Linkage});
    f.set({Feature::Enabling, Feature::Inhibiting, Feature::Timing, Feature::Firing,
           Feature::Storage, Feature::Place, Feature::Transition});
    f.set({Feature::Server, Feature::JobSink, Feature::RandomSource, Feature::InfiniteServer,
           Feature::SharedServer, Feature::ServiceTunnel});
    // Heterogeneous server pools: `jmt_writer` emits serverTypesNames,
    // serverTypesNumOfServers and serverTypesCompatibilities (jmt_writer.h:1023
    // and :1138), so jsim simulates the pools rather than a station of the same
    // total size. Undeclared until 2026-08-22, which made this port refuse a
    // model its own writer serialises and the MATLAB, JAR and python SolverJMT
    // all accept.
    f.set(Feature::HeteroServers);
    f.set({Feature::Exp, Feature::Erlang, Feature::HyperExp, Feature::PH, Feature::APH,
           Feature::Coxian, Feature::Cox2});
    f.set({Feature::Det, Feature::Gamma, Feature::Lognormal, Feature::Pareto, Feature::Weibull,
           Feature::Uniform});
    f.set({Feature::MAP, Feature::MMPP2, Feature::Replayer, Feature::Trace, Feature::Immediate,
           Feature::Disabled});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_FCFS, Feature::SchedStrategy_PS,
           Feature::SchedStrategy_DPS, Feature::SchedStrategy_GPS, Feature::SchedStrategy_LPS});
    f.set({Feature::SchedStrategy_SIRO, Feature::SchedStrategy_HOL,
           Feature::SchedStrategy_FCFSPRIO});
    f.set({Feature::SchedStrategy_PSPRIO, Feature::SchedStrategy_DPSPRIO,
           Feature::SchedStrategy_GPSPRIO});
    f.set({Feature::SchedStrategy_LCFS, Feature::SchedStrategy_LCFSPR,
           Feature::SchedStrategy_LCFSPI, Feature::SchedStrategy_LCFSPRIO,
           Feature::SchedStrategy_LCFSPRPRIO, Feature::SchedStrategy_LCFSPIPRIO});
    f.set({Feature::SchedStrategy_FCFSPR, Feature::SchedStrategy_FCFSPI,
           Feature::SchedStrategy_FCFSPRPRIO, Feature::SchedStrategy_FCFSPIPRIO});
    // save_class_parallelism writes it and server_pools carries it, so the gate
    // was refusing a model this writer emits. The reference, the JAR and native
    // Python all declare it.
    f.set(Feature::ServerParallelism);
    f.set({Feature::SchedStrategy_SEPT, Feature::SchedStrategy_LEPT, Feature::SchedStrategy_SJF,
           Feature::SchedStrategy_LJF});
    f.set({Feature::SchedStrategy_SRPT, Feature::SchedStrategy_SRPTPRIO,
           Feature::SchedStrategy_POLLING, Feature::SchedStrategy_EXT});
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND,
           Feature::RoutingStrategy_RROBIN, Feature::RoutingStrategy_WRROBIN,
           Feature::RoutingStrategy_JSQ, Feature::RoutingStrategy_SQ});
    f.set({Feature::OpenClass, Feature::ClosedClass, Feature::SelfLoopingClass});
    f.set({Feature::Cache, Feature::CacheClassSwitcher});
    f.set({Feature::ReplacementStrategy_RR, Feature::ReplacementStrategy_FIFO,
           Feature::ReplacementStrategy_SFIFO, Feature::ReplacementStrategy_LRU});
    f.set({Feature::SetupDelayOff, Feature::Reneging, Feature::Balking});
    // Queue.setRetrial: the JSIM writer picks the retrial Queue constructor and
    // writes the per-class orbit delay, so jsim simulates the orbit rather than
    // losing the refused job. Batch arrivals are NOT declared, in the reference
    // either: the arrival-strategy writer has no batch element, so the stream
    // would be written single. SolverJMT.m:474-478.
    f.set(Feature::Retrial);
    // c-server stations (`save_number_of_servers`) and finite buffers
    // (`save_buffer_capacity` with the drop rule): the JSIM writer exports both,
    // and `io::jmt_buffer_capacity_refusal` keeps refusing the buffers JMT would
    // answer unconstrained (closed WAITQ, BBS, RSRD, RETRIAL_WITH_LIMIT). The
    // JMVA block below withdraws the buffer. SolverJMT.m:479-484.
    f.set({Feature::MultiServer, Feature::FiniteCapacity});
    // Limited load dependence reaches JMT only as a SERVER COUNT:
    // `save_number_of_servers` exports max(nservers, max alpha) and `write_jmva`
    // writes the matching <ldstation>. That is exact for alpha(n) = min(n,c) and
    // for nothing else, so `solver_jmt_run_analyzer` refuses any other scaling by name --
    // a featset cannot inspect the vector.
    f.set(Feature::LoadDependence);

    // SOLVERJMT DRIVES TWO ENGINES, AND THEY ACCEPT DIFFERENT MODELS. Everything
    // above is the JSIM SIMULATOR's envelope, which is what 'default', 'jsim' and
    // 'replication' run. The 'jmva.*' names run the JMVA ANALYTICAL engine, and
    // `io::write_jmva` emits for it only a station type (delay / load-independent
    // / load-dependent), a per-chain service demand, a per-chain visit count, the
    // class populations or arrival rates and a reference station. NOTHING ELSE IN
    // THE MODEL REACHES JMVA, so declaring the JSIM envelope for jmva was a
    // promise the writer could not keep: on a three-class LRU cache model all
    // eight closed-form jmva methods returned an entirely zero table with no
    // error, jmva.mva labelled 'exact' among them.
    //
    // The DISTRIBUTIONS stay: JMVA consumes a mean service demand, so any renewal
    // law with a finite mean is admissible, exactly as it is for mva_feature_set
    // and nc_feature_set. What goes is every construct whose effect is not
    // carried by (station type, demand, visits, population).
    if (method.compare(0, 4, "jmva") == 0) {
        for (Feature ff :
             {Feature::Cache, Feature::CacheClassSwitcher, Feature::ReplacementStrategy_RR,
              Feature::ReplacementStrategy_FIFO, Feature::ReplacementStrategy_SFIFO,
              Feature::ReplacementStrategy_LRU,
              // no fork element exists, and a visit ratio cannot express the join
              Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner,
              Feature::JoinPartial, Feature::ForkFanoutVector, Feature::ForkFanoutRandom,
              Feature::ForkBranchProbability,
              // no Petri-net counterpart, and no finite capacity region
              Feature::Place, Feature::Transition, Feature::Enabling, Feature::Inhibiting,
              Feature::Timing, Feature::Firing, Feature::Storage, Feature::Region,
              // no impatience element: the abandonment would simply not happen,
              // and no retrial element either: the orbit would be lost
              Feature::Reneging, Feature::Balking, Feature::Retrial,
              // the JMVA document has no capacity element at all
              Feature::FiniteCapacity,
              // server-side attributes the JMVA document has no slot for
              Feature::SetupDelayOff, Feature::ServerParallelism, Feature::HeteroServers,
              // the writer emits NO discipline, so a priority, weighted, size-based
              // or limited-sharing station would be solved as an ordinary
              // load-independent one. Only the four BCMP station types survive the
              // encoding, the same line nc_feature_set draws.
              Feature::SchedStrategy_DPS, Feature::SchedStrategy_GPS, Feature::SchedStrategy_HOL,
              Feature::SchedStrategy_FCFSPRIO, Feature::SchedStrategy_PSPRIO,
              Feature::SchedStrategy_DPSPRIO, Feature::SchedStrategy_GPSPRIO,
              Feature::SchedStrategy_LCFSPI, Feature::SchedStrategy_LCFSPIPRIO,
              Feature::SchedStrategy_LCFSPRIO, Feature::SchedStrategy_LCFSPRPRIO,
              Feature::SchedStrategy_FCFSPR, Feature::SchedStrategy_FCFSPI,
              Feature::SchedStrategy_FCFSPRPRIO, Feature::SchedStrategy_FCFSPIPRIO,
              Feature::SchedStrategy_SEPT, Feature::SchedStrategy_LEPT,
              Feature::SchedStrategy_SJF, Feature::SchedStrategy_LJF,
              Feature::SchedStrategy_SRPT, Feature::SchedStrategy_SRPTPRIO,
              Feature::SchedStrategy_LPS, Feature::SchedStrategy_POLLING,
              // the document carries MEAN visit counts, which is not what makes a
              // join-the-shortest-queue model behave as it does
              Feature::RoutingStrategy_RROBIN, Feature::RoutingStrategy_WRROBIN,
              Feature::RoutingStrategy_JSQ, Feature::RoutingStrategy_SQ})
            f.unset(ff);
        // The eight closed-form algorithms are single-server ones
        // (`jmt::jmt_method_refusal` words it); exact MVA carries the count.
        if (jmva_is_closed_only(method))
            for (Feature ff :
                 {Feature::OpenClass, Feature::LoadDependence, Feature::MultiServer})
                f.unset(ff);
    }
    return f;
}

/**
 * SolverBA.getFeatureSet, transcribed name for name.
 *
 * WHAT MAKES A DISTRIBUTION ADMISSIBLE HERE IS ITS MEAN, which is why the list
 * of laws is longer than any product-form solver's and shorter than a
 * simulator's. `solver_ba_analyzer` reads `sn.rates` and `sn.visits` and nothing
 * else: every bound in the ABA/BJB/PB/GB/SB/Harel/MWBA families is a function of
 * the demands D = V./rates and the think time, so any renewal law with a finite
 * mean is admissible whatever its higher moments. `solver_ba_qrf_analyzer` is
 * the one that needs more, and what it needs is the {D0,D1} pair, which the
 * phase-type laws below carry; it refuses the rest through
 * `sn_to_qrf_blocking`'s own message rather than through this set.
 *
 * THE MODULATED PROCESSES (MAP, MMPP2, MMAP, BMAP) ARE DELIBERATELY OUT. Their
 * mean rate exists, so the utilization law still holds, but the bounds are
 * derived for a product-form network in which the correlation between successive
 * services does not exist: bounding such a model would report a bracket for a
 * DIFFERENT system rather than refuse it. The same reasoning keeps Cache,
 * Fork/Join and the Petri-net constructs out, none of which the analyzer
 * represents at all.
 *
 * THE STRUCTURAL NARROWING IS NOT HERE. `ba::list_valid_methods(L)` already
 * drops the reduction bounds off a model that is not a single-class closed
 * network of single servers, and the three open-network families off a closed
 * one; a flat feature set cannot express either, and duplicating the shape test
 * here is how the two drift apart.
 *
 * THE ONE ADDITION TO THE MATLAB LIST IS THE SUB-NODE SECTIONS, and it is not a
 * widening. `MNetwork.getUsedLangFeatures` emits no Buffer, Dispatcher, Server,
 * InfiniteServer, SharedServer, JobSink, RandomSource or ServiceTunnel -- its
 * node loop marks the distribution, the discipline and the routing and stops --
 * while `used_lang_features` here marks all of them for every Queue, Delay,
 * Source and Sink. Transcribing the MATLAB list literally therefore refused
 * EVERY model, an ordinary closed exponential network included, on `Buffer`;
 * `env_feature_set` and `ag_feature_set` carry the same eight names for the same
 * reason. The JAR has MATLAB's list and a JAR-side emitter that behaves like
 * this one, which is why `SolverAUTO.listValidMethods` there offers no 'ba.*'
 * method name on a model whose bounds SolverBA computes: the same defect, unfixed.
 */
/**
 * The family prefix of a bound method name: everything before the first dot.
 *
 * Local to this header so `ba_feature_set` can apply its per-method deltas
 * without reaching into `line/solvers/ba`, which includes this file.
 */
inline std::string ba_family_of(const std::string& method) {
    const std::string::size_type dot = method.find('.');
    return (dot == std::string::npos) ? method : method.substr(0, dot);
}

/**
 * The 'default'/'auto'/'qr'/'lr' aliases, duplicated from `ba::resolve_method`
 * for the same reason `ba_family_of` sits here: this header is included BY the
 * solver and cannot include it back. The two are one line apiece and are
 * asserted equal in `cpp/tests/test_gate_ba.cpp`.
 */
inline std::string ba_resolve_method_name(const std::string& method) {
    if (method == "default") return "gb.upper";
    if (method == "auto") return "auto.upper";
    if (method == "lr") return "lr.upper";
    if (method == "qr") return "qrf.mmi";
    return method;
}

inline FeatureSet ba_feature_set(const std::string& method) {
    FeatureSet f;
    f.set({Feature::ClassSwitch, Feature::Delay, Feature::DelayStation, Feature::Queue,
           Feature::Sink, Feature::Source, Feature::Router});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer, Feature::SharedServer,
           Feature::Buffer, Feature::Dispatcher, Feature::Server, Feature::JobSink,
           Feature::RandomSource, Feature::ServiceTunnel});
    // A self-looping class is a closed chain of one station, which the
    // demand-parameterized bounds read as any other chain (SolverBA.m:826-829).
    f.set({Feature::ClosedClass, Feature::SelfLoopingClass, Feature::OpenClass});
    f.set({Feature::APH, Feature::Coxian, Feature::Cox2, Feature::Erlang, Feature::Exp,
           Feature::HyperExp, Feature::PH});
    f.set({Feature::Det, Feature::Lognormal, Feature::Pareto, Feature::Uniform,
           Feature::Weibull});
    // MAP and MMPP2 are declared HERE and withdrawn below from every family but
    // `mapamva`, which is the one bound derived for a CORRELATED arrival stream.
    // Neither name was in this set at all, so `mapamva.upper`/`.lower` -- both
    // advertised by `solver_ba_runner` and sieved by `ba::method_refusal` -- were
    // refused on the only models they exist for, and this port offered no BA
    // method whatever on a MAP model where the other three offer `mapamva` alone.
    f.set({Feature::MAP, Feature::MMPP2});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS, Feature::SchedStrategy_FCFS,
           Feature::SchedStrategy_LCFSPR});
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND});
    // Petri-net constructs, for the spnlp family (solver_ba_spnlp). They were
    // out on the "no representation" ground the rest of this set rests on, and
    // that ground is gone: the spnlp relaxation is indexed by the MARKING,
    // reads the enabling, inhibiting and firing arcs out of the node
    // parameters, and refuses by name the modes it cannot carry (immediate,
    // multi-server, marking-dependent, and phase-type on its Markovian side).
    // QueueingPlace is deliberately absent: a place with an embedded queue has
    // local state the relaxation has no variable for. Same division nc draws.
    f.set({Feature::Place, Feature::Transition, Feature::Linkage, Feature::Enabling,
           Feature::Inhibiting, Feature::Timing, Feature::Firing, Feature::Storage});

    // THE PER-METHOD DELTAS, i.e. the premises a feature name CAN state. A
    // feature set says "I accept this construct", so it can refuse a model for
    // HAVING one and never for lacking one; that is exactly the shape of the
    // delay-station and closed-class premises below, and exactly not the shape
    // of "one class" or "one server", which have no feature name and are
    // refused structurally by `ba::method_refusal`. Judged on the RESOLVED name
    // so that 'default' carries the envelope of the gb.upper it runs as.
    //
    // DELAY STATIONS. 'sb' and 'lr' reject an infinite-server station outright,
    // and 'harel', 'sib' and 'scb' reject a nonzero think time, which on these
    // models is the same station: harel extrapolates the exact normalizing
    // constant of a delay-free network, SIB Section 3.2 is the extension that
    // would carry Z and is not implemented, and SCB Theorem 3 rests on the
    // delay-free balanced-network throughput. The three OPEN families reject
    // one too, each being derived for one server per station.
    //
    // CLASS TYPES. The three OPEN families drop ClosedClass, which is the whole
    // of their class premise. The MIRROR delta -- dropping OpenClass from every
    // demand-parameterized family -- is deliberately NOT applied: "supports
    // single-class closed networks only" is one rule, its single-class half has
    // no feature name, and splitting it across the two mechanisms would report
    // the closed half here and the single-class half in `ba::method_refusal`
    // for the same model. It is stated once, structurally. 'spnlp' takes no
    // delta at all: it is indexed by the marking, and whether that marking is
    // bounded is a question about the P-invariants of the net, which
    // `spn_lpbnd` answers.
    const std::string resolved = ba_resolve_method_name(method);
    const std::string fam = ba_family_of(resolved);
    if (fam == "bpt" || fam == "bgt" || fam == "snc") {
        f.unset(Feature::ClosedClass);
        f.unset(Feature::Delay);
        f.unset(Feature::DelayStation);
        f.unset(Feature::InfiniteServer);
        f.unset(Feature::SchedStrategy_INF);
        if (fam != "snc") {
            // SERVICE AND ARRIVAL LAWS. 'bpt' and 'bgt' are derived for a
            // MARKOVIAN open network and read the mean alone, so a
            // non-exponential law anywhere is not something they refuse at run
            // time -- it is something they silently bound as if it were
            // Poisson. Measured on the M/M/1 shape: replacing the Exp(1) source
            // by an Erlang of the same mean leaves bgt.upper at QLen 32.6667
            // and bpt.lower at 1.0, digit for digit. That is a bound on a
            // DIFFERENT system, so the laws are dropped here rather than left
            // to a run-time check the analyzers do not make: their own procid
            // test covers the queueing stations only and would miss exactly the
            // source case.
            //
            // 'snc' is excluded: it CONSUMES the arrival law (the same
            // substitution moves it from 3.8244 to 3.0092) and its analyzer
            // branches on a non-exponential source deliberately. Its rule is
            // about the SERVICE only, which no feature name can say, so it
            // lives in `ba::method_refusal` instead.
            f.unset(Feature::APH);
            f.unset(Feature::Coxian);
            f.unset(Feature::Cox2);
            f.unset(Feature::Erlang);
            f.unset(Feature::HyperExp);
            f.unset(Feature::PH);
            f.unset(Feature::Det);
            f.unset(Feature::Lognormal);
            f.unset(Feature::Pareto);
            f.unset(Feature::Uniform);
            f.unset(Feature::Weibull);
        }
    } else if (fam == "sb" || fam == "harel" || fam == "sib" || fam == "scb" || fam == "lr") {
        f.unset(Feature::Delay);
        f.unset(Feature::DelayStation);
        f.unset(Feature::InfiniteServer);
        f.unset(Feature::SchedStrategy_INF);
    }
    // THE CORRELATED-ARRIVAL BOUND, which is the one family whose premise a
    // feature name states positively. Every other bound reads the MEAN demand
    // alone, so a MAP or MMPP2 anywhere would be bounded as if it were Poisson --
    // a bracket for a DIFFERENT system, and silent. Written as a delta on the
    // complement rather than as a grant because a feature set can refuse a model
    // for HAVING a construct and never for lacking one.
    //
    // 'mapamva' takes the delay delta instead: its LP is a network of queues, and
    // Casale-Smirni name the delay extension as open work.
    if (fam == "mapamva") {
        f.unset(Feature::Delay);
        f.unset(Feature::DelayStation);
        f.unset(Feature::SchedStrategy_INF);
    } else {
        f.unset(Feature::MAP);
        f.unset(Feature::MMPP2);
    }
    // MULTISERVER, SolverBA.m:574-580, is out of the base envelope: every
    // demand-parameterized family reads one server per station (`ba::method_refusal`
    // names 'ssd' as the alternative), the alpha-free QRF arms refuse a c-server
    // station through `sn_to_qrf_alpha` and the open families through their own
    // refusal. What carries the count is granted here: 'ssd' (the multiserver
    // bound), 'ldbcmp' (its fixed-rate form runs on the c-server rate law),
    // 'auto' (which picks among them) and the two load-dependent QRF arms, whose
    // alpha(i,n) IS min(n,c).
    if (fam == "ssd" || fam == "ldbcmp" || fam == "auto" || resolved == "qrf.mmi.ld" ||
        resolved == "qrf.mmi.linear")
        f.set(Feature::MultiServer);
    // FINITECAPACITY, SolverBA.m:581-591: only the QRF blocking bounds carry the
    // buffer (the MM, MM1, ZZ, ZM, BB, F tables), which is the same split
    // `ba::ignores_blocking` makes; the structural refusal keeps naming them.
    // 'spnlp' is NOT granted, though `ignores_blocking` exempts it: its polytope
    // reads no Place capacity, so a capped place would be relaxed away.
    if (resolved.compare(0, 7, "qrf.bas") == 0 || resolved.compare(0, 8, "qrf.rsrd") == 0)
        f.set(Feature::FiniteCapacity);
    return f;
}

/**
 * SolverQNS.getFeatureSet, transcribed WHOLE.
 *
 * THE SECOND SET IN THIS FILE THAT IS NOT NARROWED, for `ldes_feature_set`'s
 * reason and no other: SolverQNS implements no numerics. It writes the JMVA
 * document `qnsolver` reads, or converts the model through `qn2lqn` and hands it
 * to `lqns` -- the same two routes, through the same two external binaries, that
 * `solver_qns.m` drives. The code that has to support a declared feature is
 * therefore the binary the reference also runs, so narrowing here would refuse
 * models the reference answers with the very process this port is about to
 * start.
 *
 * WHERE THE REAL REFUSALS LIVE, all of them at solve time and none expressible
 * as a boolean over feature names:
 *   qns::detail::check_supported   a station that is not a Queue, a Delay or a
 *                                  Source, and class priorities, on the JMVA
 *                                  route only -- the layered route legitimately
 *                                  carries Join nodes the document cannot encode
 *   qn::check_binding_capacity     a finite buffer, which neither binary models
 *   the multiserver gate           'suri' and 'schmidt' reach LQNS and not
 *                                  `qnsolver -m`
 *
 * SOURCE: `matlab/src/solvers/wrappers/QNS/@@SolverQNS/SolverQNS.m`, name for
 * name. `Normal` is the one name dropped: it has a registry enumerator here, and
 * `write_jmva` reduces every law to its mean, so a Normal service is written as
 * its mean like any other -- but a Normal law admits negative samples and the
 * reference's own JMT writer refuses it (see `jmt_feature_set`), so declaring it
 * on this route alone would be the only place in the port that accepts it.
 */
inline FeatureSet qns_feature_set(const std::string& /*method*/) {
    FeatureSet f;
    f.set({Feature::Sink, Feature::Source, Feature::Router, Feature::ClassSwitch,
           Feature::Delay, Feature::DelayStation, Feature::Queue});
    f.set({Feature::Fork, Feature::Join, Feature::Forker, Feature::Joiner});
    f.set({Feature::Logger, Feature::LogTunnel});
    f.set({Feature::Coxian, Feature::Cox2, Feature::APH, Feature::Erlang, Feature::Exp,
           Feature::HyperExp, Feature::PH});
    f.set({Feature::Det, Feature::Gamma, Feature::Lognormal, Feature::Pareto, Feature::Uniform,
           Feature::Weibull});
    f.set({Feature::MAP, Feature::MMPP2, Feature::Replayer, Feature::Trace});
    f.set({Feature::StatelessClassSwitcher, Feature::InfiniteServer, Feature::SharedServer,
           Feature::Buffer, Feature::Dispatcher, Feature::Server, Feature::JobSink,
           Feature::RandomSource, Feature::ServiceTunnel, Feature::Linkage});
    f.set({Feature::Enabling, Feature::Timing, Feature::Firing, Feature::Storage,
           Feature::Place, Feature::Transition});
    f.set({Feature::SchedStrategy_INF, Feature::SchedStrategy_PS, Feature::SchedStrategy_DPS,
           Feature::SchedStrategy_FCFS, Feature::SchedStrategy_GPS, Feature::SchedStrategy_SIRO,
           Feature::SchedStrategy_HOL, Feature::SchedStrategy_LCFS,
           Feature::SchedStrategy_LCFSPR});
    f.set({Feature::SchedStrategy_SEPT, Feature::SchedStrategy_LEPT, Feature::SchedStrategy_SJF,
           Feature::SchedStrategy_LJF, Feature::SchedStrategy_EXT});
    f.set({Feature::RoutingStrategy_PROB, Feature::RoutingStrategy_RAND,
           Feature::RoutingStrategy_RROBIN, Feature::RoutingStrategy_WRROBIN,
           Feature::RoutingStrategy_SQ});
    f.set({Feature::ClosedClass, Feature::OpenClass});
    // c-server stations: the JMVA document carries the count as an <ldstation>
    // and the LQN as a host multiplicity; 'suri' and 'schmidt' refuse one on the
    // qnsolver path, which stays structural (the multiserver gate in solver_qns.h).
    // FiniteCapacity is NOT declared, as in SolverQNS.m:183-188: neither document
    // has a buffer, and `qn::check_binding_capacity` words the refusal.
    f.set(Feature::MultiServer);
    return f;
}

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_SOLVER_FEATURE_SETS_H
