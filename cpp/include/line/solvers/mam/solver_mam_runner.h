/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_RUNNER_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_RUNNER_H

/**
 * The SolverMAM class surface: `@@SolverMAM/runAnalyzer.m` and the gates around
 * it.
 *
 * Everything that happens BEFORE and AFTER one inner solve lives here: the
 * method whitelist (`SolverMAM.listValidMethods`), the per-method structural
 * gates (`SolverMAM.supportsModelMethod`), the process gate, and the metric
 * filter. `runAnalyzer`'s Finite Capacity Region rejection needs no counterpart:
 * the C++ `NetworkStruct` has no region field, so such a model cannot be built.
 *
 * THE METRIC FILTER IS SHARED WITH SolverMVA AND SolverNC VERBATIM
 * (`filter_metric`, `sn_get_residt_from_respt`, `sn_get_arvr_from_tput` in
 * `solver_mva_runner.h`), because `@@NetworkSolver/getAvg` is
 * solver-independent: it is literally the same code path for all three in the
 * reference.
 *
 * NON-MARKOVIAN SERVICE IS CONVERTED, NOT REFUSED (since 2026-08-01).
 * `solver_mam_analyzer.m:40` runs `sn_nonmarkov_toph` on its own copy of the
 * struct, and so does `solver_mam_solve` here: a Gamma, Weibull, Lognormal,
 * Pareto or Uniform is replaced by the two-moment concentrated-ME fit, or by the
 * 20-phase Bernstein density fit under `phfit_cme = false`, and a Det by a
 * 20-phase Erlang when `preserve_det` is off. The conversion runs AFTER the
 * slotted test, which would otherwise see a continuous fit where a Geometric was
 * declared. What survives the conversion and still reaches `check_processes` is
 * NHPP / MAPt / PHt, whose content is a schedule that no homogeneous surrogate
 * can carry, and those are refused by name.
 *
 * TWO GATES REMAIN STRICTER HERE THAN IN THE REFERENCE, and both are the honest
 * outcome rather than a limitation:
 *
 *  - A CLASS-SWITCHING CHAIN AT AN FCFS STATION. The reference's own chain-1
 *    mark collapse makes the marking handed to the queue solver disagree with
 *    the class count, and MATLAB then fails on a dimension mismatch. The
 *    analyzer refuses by name with both counts; see solver_mam_basic.h.
 *  - LOAD DEPENDENCE UNDER 'default' OUTSIDE THE LD-QBD SHAPE. MATLAB has no
 *    counterpart because it has no LoadDependence declaration for MAM at all
 *    (mam_feature_set's declaration is a C++-only widening, solver_feature_sets.h).
 *    `check_model_method` refuses a load-dependent model here when its shape is
 *    not the one `solver_mam_ldqbd` requires, because the featset gate cannot
 *    see topology and 'default' would otherwise fall through to
 *    `solver_mam_basic`, which never reads `st.lldscaling`.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/api/sn/sn_is_discrete_time.h"
#include "line/api/sn/sn_nonmarkov_toph.h"
#include "line/solvers/mam/mam_dispatch.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mam/solver_mam_fj.h"
#include "line/solvers/mam/solver_mam_ldqbd_transient.h"
#include "line/solvers/mam/solver_mam_passage_time.h"
#include "line/solvers/mam/solver_mam_prob.h"
#include "line/solvers/mam/solver_mam_transient_qbd.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/error.h"

namespace line {
namespace mam {

using lang::GlobalConstants;

/**
 * Port of `SolverMAM.listValidMethods`.
 *
 * The reference's list verbatim, in its order (the comment there records that
 * test files index into it, so nothing may be inserted in the middle). A listed
 * name that has no ported analyzer passes this gate and is then refused BY NAME
 * by the dispatch, which is the honest outcome and different from silently
 * solving a different model.
 *
 * `exact` at the end is the RCAT ALIAS -- `mam_dispatch.h` groups it with
 * `inap`/`inapplus`/`inapinf`, exactly as `solver_mam_analyzer.m` does -- and
 * NOT the retired autocat of `solver_mam_autocat.h`, whose own refusal message
 * predates the reference re-advertising the name. `retrial` names the
 * BMAP/PH/N/N analyzer that `default` also resolves to on a retrial topology.
 */
inline std::vector<std::string> list_valid_methods() {
    // SolverMAM.m verbatim, in its order, MINUS the RCAT names
    // ('inap', 'inapplus', 'inapinf', 'exact') which moved to SolverAG;
    // mam_dispatch refuses those by name and says where they went.
    return {"default", "dec.source", "dec.mmap", "dec.poisson", "mna", "ldqbd", "dec.source.mmap", "bgchain", "retrial"};
}

/**
 * True for the four RCAT names that moved to SolverAG, with the redirect message.
 *
 * It lives here, and is called from `check_method` BEFORE the unlisted-method
 * throw, because that ordering is the whole point: the names are no longer in
 * `list_valid_methods`, so a caller carrying an old `options.method` would
 * otherwise be told only that the method is "unsupported by this solver" and
 * left to guess. `check_model_method` calls the same helper so the two paths
 * cannot drift into two different messages.
 */
inline bool rcat_moved_to_ag(const std::string& method, std::string& why) {
    if (method != "inap" && method != "inapplus" && method != "inapinf"
        && method != "exact") {
        return false;
    }
    why = "SolverMAM: the '" + method + "' method moved to SolverAG: RCAT decomposes the "
          "model into cooperating agents rather than decomposing traffic, and no MAM "
          "algorithm shares its machinery. Solve it with -s ag";
    return true;
}

/** An unlisted method is refused; one that MOVED is redirected by name. */
inline void check_method(const std::string& method) {
    const std::vector<std::string> valid = list_valid_methods();
    if (std::find(valid.begin(), valid.end(), method) != valid.end()) return;
    std::string why;
    if (rcat_moved_to_ag(method, why)) throw UnsupportedError(why);
    throw UnsupportedError("SolverMAM: the '" + method + "' method is unsupported by this solver");
}

namespace runner_detail {

/**
 * The distributions `getFeatureSet` lists, plus the Markovian ones.
 *
 * Runs AFTER `sn_nonmarkov_toph`, so the families that have a Markovian
 * surrogate are already gone; what can still arrive is a schedule family, whose
 * whole content is the time dependence a single generator cannot carry.
 *
 * `is_discrete` is set when sn_is_discrete_time claimed the model: the lattice
 * laws are then solved on their own time scale by solver_mam_dt and must not be
 * refused for lacking a continuous phase-type fit, which is the only reason the
 * refusal below exists.
 */
template <class T>
void check_processes(const qn::NetworkStruct<T>& L, bool is_discrete) {
    using lang::ProcessType;
    if (is_discrete) return;
    for (std::size_t i = 0; i < L.nstations; ++i)
        for (std::size_t r = 0; r < L.nclasses; ++r) {
            if (L.disabled[i][r]) continue;
            const ProcessType p = L.service[i][r].type;
            if (basic_detail::is_markovian_type(p) || p == ProcessType::DET) continue;
            throw UnsupportedError(
                std::string("SolverMAM: station '") + L.stations[i].name + "' class '" +
                L.classes[r].name + "' uses a " + lang::process_to_text(p) +
                " service process, which sn_nonmarkov_toph has no Markovian surrogate for: its "
                "content is a piecewise-constant SCHEDULE, and collapsing it to one homogeneous "
                "generator would answer for a stationary model the caller did not describe. Solve "
                "it with SolverFLD or SolverLDES, which integrate the schedule");
        }
}

/**
 * True when any station declares a non-trivial load-dependent scaling, the
 * same test `used_lang_features` (feature_set.h) applies to derive
 * Feature::LoadDependence.
 */
template <class T>
bool has_load_dependence(const qn::NetworkStruct<T>& L) {
    for (const qn::Station<T>& st : L.stations)
        if (!st.lldscaling.empty()) return true;
    return false;
}

/** True when the model declares at least one G-network signal class. */
template <class T>
bool has_signal_class(const qn::NetworkStruct<T>& L) {
    for (std::size_t r = 0; r < L.issignal.size(); ++r)
        if (L.issignal[r]) return true;
    return false;
}

/** `SolverMAM.supportsModelMethod`, the per-method structural gates. */
template <class T>
void check_model_method(const qn::NetworkStruct<T>& L, const std::string& method) {
    using lang::ProcessType;
    // The RCAT methods moved to SolverAG. `check_method` normally reports this
    // first; the same helper is used here so a direct caller of
    // check_model_method (the twin of SolverMAM.supportsModelMethod, which the
    // reference also answers false for) gets the identical message.
    {
        std::string why;
        if (rcat_moved_to_ag(method, why)) throw UnsupportedError(why);
    }
    // G-network signals belong to the RCAT analyzer alone: no MAM algorithm
    // reads issignal, so every one of them would solve the model with the
    // signals turned into ordinary customers and report that as the answer.
    if (has_signal_class(L)) {
        throw UnsupportedError(
            "SolverMAM: the " + method + " method does not support G-network signals: no MAM "
            "algorithm reads issignal, so this method would solve the model with every signal "
            "turned into an ordinary customer. Use SolverAG (-s ag), whose RCAT methods are "
            "the only ones that model signals");
    }
    if (method == "mna") {
        bool open = false, closed = false;
        for (const qn::JobClass& c : L.classes) {
            if (std::isinf(c.population)) open = true;
            else closed = true;
        }
        if (open && closed)
            throw UnsupportedError(
                "SolverMAM: the mna method does not support mixed open/closed models");
        // mam_feature_set declares RoutingStrategy_RROBIN for 'mna' because
        // solver_mna_open resolves the deterministic split through
        // npfqn_traffic_split_rr. solver_mna_closed has no counterpart -- the
        // reference corrects the open traffic equations only -- so a closed
        // round-robin model would be solved as if the dispatcher were random.
        if (!open) {
            for (const qn::NodeDef& nd : L.nodes)
                for (qn::RoutingStrategy rs : nd.routing)
                    if (rs == qn::RoutingStrategy::RROBIN)
                        throw UnsupportedError(
                            "SolverMAM: the mna method supports round-robin routing in open "
                            "models only");
        }
    } else if (method == "default") {
        // mam_feature_set declares LoadDependence for 'default' because the
        // single-class closed Delay+Queue shape routes here to solver_mam_ldqbd
        // (mam_dispatch.h branch 2e), which reads st.lldscaling. The featset
        // cannot see topology, so a load-dependent model outside that shape --
        // or one carrying a setup/delay-off pair, which steals the same shape
        // ahead of ldqbd -- must be refused here rather than falling through
        // to solver_mam_basic, which never reads the field and would solve
        // every level at the nominal rate.
        if (has_load_dependence(L) &&
            !(L.setupparam.empty() && dispatch_detail::is_closed_delay_queue(L)))
            throw UnsupportedError(
                "SolverMAM: this model uses load-dependent service rates outside the "
                "single-class closed Delay+Queue shape that solver_mam_analyzer.m routes to the "
                "level-dependent QBD (solver_mam_ldqbd); the 'default' method would otherwise "
                "fall through to solver_mam_basic's dec.source decomposition, which does not "
                "read the load-dependent scaling and would solve every level at the nominal "
                "rate. Call method 'ldqbd' directly, which refuses by name if the shape still "
                "does not fit, or restructure the model to the required shape");
    } else if (method == "ldqbd") {
        if (L.nclasses != 1)
            throw UnsupportedError("SolverMAM: the ldqbd method requires a single-class model");
    } else if (method == "retrial") {
        // A "MUST BE PRESENT" RULE, which is why it cannot live in a feature
        // set: a FeatureSet says "I accept this construct", so it can refuse a
        // model for HAVING something and never for LACKING it.
        // solver_mam_retrial needs the BMAP/PH/N/N bufferless retrial topology
        // of Dudin et al. to analyze, and a model without one is not a smaller
        // retrial model, it is a different one. mam_retrial_refusal is the same
        // predicate mam_dispatch's 'retrial' arm and auto_family_refusal ask,
        // and it reports WHICH requirement the model missed.
        const std::string retrialWhy = mam_retrial_refusal(L);
        if (!retrialWhy.empty()) throw UnsupportedError(retrialWhy);
    } else if (method == "bgchain") {
        // The closed classes ARE the background chain, so a purely open model has
        // nothing to build it from. A purely CLOSED one is accepted: it is the
        // degenerate case where the chain answers alone, with no open work to
        // modulate it.
        bool hasClosed = false;
        for (std::size_t k = 0; k < L.nclasses; ++k) {
            if (!std::isinf(L.classes[k].population)) hasClosed = true;
        }
        if (!hasClosed)
            throw UnsupportedError(
                "SolverMAM: the bgchain method requires at least one closed class: the background "
                "chain IS the closed population vector, which a purely open model does not have. "
                "Use the dec.source method");
        // Priority disciplines need the per-class QBD of MMAPPH1PRPR, which has
        // no counterpart in the modulated level-dependent QBD this method builds.
        bool prioSched = false;
        for (std::size_t i = 0; i < L.nstations; ++i) {
            const lang::SchedStrategy sc = L.stations[i].sched;
            if (sc == lang::SchedStrategy::HOL || sc == lang::SchedStrategy::FCFSPRPRIO)
                prioSched = true;
        }
        if (prioSched && L.has_distinct_priorities())
            throw UnsupportedError(
                "SolverMAM: the bgchain method does not support class priorities: it aggregates "
                "the open classes into one phase-type mixture per station, which cannot express a "
                "priority order. Use the dec.source method");
        if (L.has_fork())
            throw UnsupportedError(
                "SolverMAM: the bgchain method does not support fork-join: the background chain "
                "conserves the closed population per station, which a fork violates. Use the "
                "dec.source method");
    }
}

}  // namespace runner_detail

/**
 * `check_model_method` asked WITHOUT raising: the same verdict as a sentence.
 *
 * `autosolver::auto_family_refusal` needs a reason rather than an exception --
 * the report answers yes or no per (family, method) pair and prints why -- and
 * asking the gate itself is what keeps the offered pairs and the runnable ones
 * the same set. Duplicating the rules there is how a method comes to be offered
 * on a model its own runner refuses, which is exactly what happened to
 * 'dec.mmap' and 'retrial'.
 */
template <class T>
std::string mam_model_method_refusal(const qn::NetworkStruct<T>& L, const std::string& method) {
    try {
        runner_detail::check_model_method(L, method);
    } catch (const std::exception& e) {
        return std::string(e.what());
    }
    return std::string();
}

/**
 * The gates and the dispatch of `@@SolverMAM/runAnalyzer.m`, without the metric
 * filter.
 */
template <class T>
MamSolution<T> solver_mam_solve(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    check_method(opt.method);
    runner_detail::check_model_method(L, opt.method);
    // Whether the slotted path will claim the model, which decides both the
    // process gate below and branch -1 of the dispatch. Asking here keeps the
    // two answers from diverging.
    bool is_discrete = false;
    if constexpr (num_traits<T>::has_transcendental) {
        api::DiscreteTimeOptions dtopt;
        dtopt.timescale = opt.timescale;
        dtopt.slotlength = opt.slotlength;
        double slot = 1.0;
        api::DiscreteTimeInfo dtinfo;
        is_discrete = api::sn_is_discrete_time(L, dtopt, &slot, &dtinfo);
    }

    // `solver_mam_analyzer.m:40` runs sn_nonmarkov_toph on its own copy, AFTER
    // the slotted test (which would otherwise see a continuous CME where a
    // Geometric was declared) and with preserveDet on, so the exact MAP/D/c
    // branch still sees its Det. The struct is copied only when something in it
    // would actually be replaced.
    // The RCAT methods build a CTMC per component out of (D0,D1), so a preserved
    // Det would reach them with no matrix at all and be read back as its mean
    // rate, and a concentrated matrix exponential is not a generator at all.
    // They need a genuine phase-type, as SSA, Fluid and JMT do.
    const bool preserve_det = opt.preserve_det;
    qn::NetworkStruct<T> converted;
    const qn::NetworkStruct<T>* Lp = &L;
    if constexpr (num_traits<T>::has_transcendental) {
        if (!is_discrete && api::sn_has_nonmarkov(L, preserve_det)) {
            converted = L;
            api::NonmarkovOptions no;
            no.order = opt.nonmkv_order;
            no.phfit = opt.phfit_cme ? api::PhFit::Cme : api::PhFit::Ph;
            no.preserve_det = preserve_det;
            api::sn_nonmarkov_toph(converted, no);
            Lp = &converted;
        }
    }
    const qn::NetworkStruct<T>& Lc = *Lp;

    runner_detail::check_processes(Lc, is_discrete);
    // runAnalyzerChecks' universal feature gate, AFTER the three structural
    // checks so their specific messages (map_bernstein, single-class ldqbd,
    // single-server inap, and the load-dependence shape check just above)
    // still win. The declared set carries LoadDependence over MATLAB's for
    // 'default' and 'ldqbd' only, without which every ldqbd model would
    // refuse here; every other method leaves it undeclared, so a
    // load-dependent model under 'dec.source', 'dec.poisson' or 'mna' is
    // refused cleanly by the gate instead of silently solved at nominal rates.
    qn::feature_gate("SolverMAM", qn::mam_feature_set(opt.method), Lc);
    return mam_dispatch(Lc, opt);
}

/**
 * Port of `@@SolverMAM/runAnalyzer.m` for the `lang='matlab'` path: solve, then
 * apply the metric filter `@@NetworkSolver/getAvg` puts between the analyzer and
 * the caller.
 */
template <class T>
mva::AvgResult<T> solver_mam_run_analyzer(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    const std::string origmethod = opt.method;
    const MamSolution<T> d = solver_mam_solve(L, opt);
    const mva::MvaSolution<T>& s = d.sol;

    const std::size_t M = L.nstations, K = L.nclasses;
    std::vector<std::vector<bool>> mask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            mask[i][k] = num_traits<T>::to_double(s.R(i, k)) < 10.0 * GlobalConstants::FineTol;
    std::vector<std::vector<bool>> srcmask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Source)
            for (std::size_t k = 0; k < K; ++k) srcmask[i][k] = true;

    mva::AvgResult<T> out;
    out.QN = mva::filter_metric(L, s.Q, mva::MetricKind::QLen, &mask);
    out.UN = mva::filter_metric(L, s.U, mva::MetricKind::Util, &mask);
    out.RN = mva::filter_metric(L, s.R, mva::MetricKind::RespT, nullptr);
    out.TN = mva::filter_metric(L, s.Tp, mva::MetricKind::Tput, nullptr);
    out.WN = mva::filter_metric(L, mva::sn_get_residt_from_respt(L, out.RN),
                                mva::MetricKind::ResidT, nullptr);
    out.AN = mva::filter_metric(L, mva::sn_get_arvr_from_tput(L, out.TN), mva::MetricKind::ArvR,
                                &srcmask);
    out.CN = s.C;
    out.XN = s.X;
    out.method = origmethod;
    // `runAnalyzer` reports 'default/<algorithm>' when the caller asked for the
    // default, so the algorithm that produced the numbers is never lost.
    out.actualmethod = (origmethod == "default" && !d.actualmethod.empty())
                           ? "default/" + d.actualmethod
                           : d.actualmethod;
    out.iter = s.iter;
    return out;
}

// ---------------------------------------------------------------------------
// The rest of the @@SolverMAM class surface.
//
// Each of these mirrors one accessor and does its own gating; they take the
// converged `AvgResult` where the reference calls `self.getAvg()` first, so a
// caller cannot reach them on an unsolved model by accident.
// ---------------------------------------------------------------------------

/** `@@SolverMAM/getProb.m`: the joint (level, phase) table at a node. */
template <class T>
ProbTable<T> solver_mam_get_prob(const qn::NetworkStruct<T>& L, const MamOptions& opt,
                                 std::size_t node, const mva::AvgResult<T>& avg) {
    return solver_mam_getprob(L, opt, node, avg);
}

/** `@@SolverMAM/getProbMarg.m`: P(n jobs of one class) at a station. */
template <class T>
std::vector<T> solver_mam_get_prob_marg(const qn::NetworkStruct<T>& L, const MamOptions& opt,
                                        std::size_t ist, std::size_t jobclass,
                                        const mva::AvgResult<T>& avg) {
    return solver_mam_getprobmarg(L, opt, ist, jobclass, avg);
}

/** `@@SolverMAM/getMAMResult.m`: the M/G/1-type internals of a single queue. */
template <class T>
qsys::BmapM1Result<T> solver_mam_get_mam_result(const qn::NetworkStruct<T>& L) {
    return solver_mam_getmamresult(L);
}

/**
 * `@@SolverMAM/getCdfRespT.m`: the response-time CDF per class.
 *
 * `getSjrnT` and `sjrnT` are aliases of this in the reference and are not given
 * separate entry points here; there is nothing for them to do that this does
 * not already do, and a second name for one function is a maintenance cost
 * rather than a feature.
 */
template <class T>
std::vector<RespTCdf<T>> solver_mam_get_cdf_respt(const qn::NetworkStruct<T>& L,
                                                  const MamOptions& opt) {
    check_method(opt.method);
    runner_detail::check_model_method(L, opt.method);
    // The same non-Markovian conversion the steady-state entry runs. NOT
    // exempted for a slotted model: solver_mam_passage_time reads the CONTINUOUS
    // phase-type fit, so a lattice law has to be refused here even though the
    // steady-state entry solves it on its own time scale.
    qn::NetworkStruct<T> converted;
    const qn::NetworkStruct<T>* Lp = &L;
    if constexpr (num_traits<T>::has_transcendental) {
        if (api::sn_has_nonmarkov(L, opt.preserve_det)) {
            converted = L;
            api::NonmarkovOptions no;
            no.order = opt.nonmkv_order;
            no.phfit = opt.phfit_cme ? api::PhFit::Cme : api::PhFit::Ph;
            no.preserve_det = opt.preserve_det;
            api::sn_nonmarkov_toph(converted, no);
            Lp = &converted;
        }
    }
    runner_detail::check_processes(*Lp, false);
    return solver_mam_passage_time(*Lp, opt);
}

/** `@@SolverMAM/getSjrnT.m` and `sjrnT.m`, both aliases of getCdfRespT. */
template <class T>
std::vector<RespTCdf<T>> solver_mam_get_sjrn_t(const qn::NetworkStruct<T>& L,
                                               const MamOptions& opt) {
    return solver_mam_get_cdf_respt(L, opt);
}

/**
 * Whether `getPerctRespT` reads the FJ_codes table rather than inverting a CDF.
 *
 * The reference decides by whether the last solve left a `result.Percentile`
 * behind, which happens exactly when the analyzer took branch 2a. That is a
 * state test on a stateless interface here, so the condition is recomputed from
 * the model. Exposed rather than inlined because the CLI needs the SAME answer:
 * a fork-join model has no response-time CDF at all, so a caller that asked for
 * one must be told which of the two it is getting.
 *
 * Compiled out below double: `mam_fj_is_homogeneous` fits a phase-type per
 * branch and is not instantiable in an exact field.
 */
template <class T>
bool mam_has_fj_percentiles(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L;
        (void)opt;
        return false;
    } else {
        // Branch 2a is reachable from 'default' and 'dec.source' only; any other
        // method names a different analyzer, which stores no percentile table.
        if (opt.method != "default" && opt.method != "dec.source") return false;
        if (L.fj.empty()) return false;
        return mam_fj_is_homogeneous(L).ok;
    }
}

/**
 * `@@SolverMAM/getPerctRespT.m`: response-time percentiles per class.
 *
 * TWO PATHS, split by `mam_has_fj_percentiles`: a homogeneous fork-join model
 * reads the table FJ_codes stored, and everything else inverts the
 * response-time CDF.
 */
template <class T>
std::vector<std::vector<T>> solver_mam_get_perct_respt(const qn::NetworkStruct<T>& L,
                                                       const MamOptions& opt,
                                                       const std::vector<double>& percentiles) {
    if (mam_has_fj_percentiles(L, opt))
        return solver_mam_fj_percentiles(L, opt, percentiles);
    const std::vector<RespTCdf<T>> rd = solver_mam_get_cdf_respt(L, opt);
    std::vector<std::vector<T>> out;
    out.reserve(rd.size());
    for (const RespTCdf<T>& c : rd) out.push_back(mam_percentiles_from_cdf(c, percentiles));
    return out;
}

/**
 * Port of `@@SolverMAM/getTranAvg.m`: transient queue length, utilization and
 * throughput.
 *
 * The reference forces `options.method = 'ldqbd'` and `options.lang = 'matlab'`
 * before delegating, so the transient path is NOT selected by the caller's
 * method; it is fixed. Reproduced: this ignores `opt.method` entirely.
 *
 * Which engine runs is decided by `mam_transient_qbd_applicable`: a correlated
 * or non-Poisson arrival, or a non-renewal service, goes to the Laplace-domain
 * transient QBD, and everything else to the QBD fast path.
 */
template <class T>
TranResult<T> solver_mam_get_tran_avg(const qn::NetworkStruct<T>& L, const MamOptions& opt_in) {
    MamOptions opt = opt_in;
    opt.method = "ldqbd";
    runner_detail::check_model_method(L, opt.method);
    // As in the passage-time entry: the transient engines are continuous-time,
    // so a lattice law is refused rather than read through its continuous fit.
    runner_detail::check_processes(L, false);
    // The transient entry never reaches solver_mam_solve, so it needs its own
    // gate, and that gate is NOT the steady-state ldqbd envelope even though the
    // method name says ldqbd. `ldqbd` withdraws FiniteCapacity because the
    // steady-state LD-QBD solves the buffer away, but BOTH transient engines
    // honour it: solver_mam_transient_qbd builds the level-dependent generator up
    // to the cap, and solver_mam_ldqbd_transient integrates that same generator.
    // The M/M/1/5 case in test_mam.cpp checks the exact stationary law the
    // construction never uses, so it is a genuine oracle for the limit. Gating
    // this entry on the steady-state set refused a model these engines solve.
    qn::FeatureSet tranSet = qn::mam_feature_set(opt.method);
    tranSet.set(qn::Feature::FiniteCapacity);
    qn::feature_gate("SolverMAM", tranSet, L);
    if (mam_transient_qbd_applicable(L)) return solver_mam_transient_qbd(L, opt);
    return solver_mam_ldqbd_transient(L, opt);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_RUNNER_H
