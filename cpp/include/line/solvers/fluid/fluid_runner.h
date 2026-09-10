/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_RUNNER_H
#define LINE_SOLVERS_FLUID_FLUID_RUNNER_H

/**
 * The fluid solver's outermost entry point: `@@SolverFLD/runAnalyzer.m`'s method
 * resolution over `solver_fluid_analyzer.m`'s dispatch.
 *
 * WHY THIS IS A SEPARATE HEADER AND NOT ANOTHER BRANCH IN `fluid_dispatch`.
 * Three of the reference's branches are ports that sit ABOVE `solver_fluid.h`
 * and include it -- `solver_fluid_closing.m` (fluid_closing.h),
 * `solver_fld_cacheqn_analyzer.m` (fluid_cacheqn.h) and
 * `@@SolverFLD/exportODEs.m` (fluid_export_odes.h). A dispatcher inside
 * `solver_fluid.h` therefore cannot call them without a cyclic include, and
 * lowering them into it would put three unrelated solvers in one file. The
 * routing goes here instead, in the one place that may include all four, which
 * is also the pattern the CTMC, MAM, MVA, NC and SSA runners follow.
 *
 * WHAT `fluid_dispatch` STILL REFUSES, and why that is not a duplicate gate.
 * It cannot see the cache analyzer, so a Cache model reaching it directly is
 * refused BY NAME rather than solved by the wrong method; the refusal names
 * this header as the entry point that does route it. Both are true statements
 * about the function that carries them.
 *
 * THE CORRECTION IS APPLIED ONCE, HERE. `solver_fluid_analyzer.m` runs its
 * utilization and response-time correction after the method switch and on every
 * branch alike, and the per-method functions it calls return the uncorrected
 * table for exactly that reason. `solver_fluid` is itself that switch plus the
 * correction, so the branches routed around it are corrected here and the ones
 * routed through it are not corrected twice.
 */

#include <cstddef>
#include <algorithm>
#include <vector>
#include <string>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/fluid/fluid_cacheqn.h"
#include "line/solvers/fluid/fluid_closing.h"
#include "line/solvers/fluid/fluid_dae.h"
#include "line/solvers/fluid/fluid_export_odes.h"
#include "line/solvers/fluid/fluid_kp.h"
#include "line/solvers/fluid/fluid_qsys.h"
#include "line/solvers/fluid/fluid_moments.h"
#include "line/solvers/fluid/fluid_nonhyperbolic.h"
#include "line/solvers/fluid/fluid_passage.h"
#include "line/solvers/fluid/fluid_stiff.h"
#include "line/solvers/fluid/fluid_petri.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/solvers/mva/fj_driver.h"
#include "line/solvers/mva/fj_ht.h"
#include "line/util/error.h"

namespace line {
namespace fluid {

/**
 * Port of `SolverFLD.listValidMethods`.
 *
 * Every method the dispatch below accepts, INCLUDING the `fluid.`-qualified
 * spelling of each: the reference's own dispatch took `fluid.tbi`, `fluid.dae`
 * and six more while its list named none of them, and native python listed the
 * whole set, so the three lists disagreed. `butools` (the MFQ backend) and
 * `aoi` (its age-of-information reading) are aliases of `mfq`.
 */
inline std::vector<std::string> fluid_list_valid_methods() {
    return {"default",
            "matrix",    "fluid.matrix",    "pnorm",  "fluid.pnorm",
            "softmin",   "fluid.softmin",
            "statedep",  "fluid.statedep",
            "closing",   "fluid.closing",
            "minnormal", "fluid.minnormal",
            "refined",   "fluid.refined",
            "tbi",       "fluid.tbi",
            "diffusion", "fluid.diffusion",
            "mfq",       "fluid.mfq",       "butools",
            "rmf",       "fluid.rmf",
            "aoi",       "fluid.aoi",
            "kp",        "fluid.kp",
            "dae",       "fluid.dae",
            // The single-station fluid limits (Source -> Queue -> Sink, one
            // class), which `fluid_qsys_handles` routes to `solver_fluid_qsys`.
            // `ggisgi` and `tga` are the short spellings `fluid_qsys_canonical`
            // maps onto the two primary names.
            "ggisgi.fluid", "fluid.ggisgi", "ggisgi",
            "ggingi.tga",   "fluid.tga",    "tga",
            "tvms",         "fluid.tvms",
            "mtginf",       "fluid.mtginf",
            "mol",          "fluid.mol"};
}

/** Port of `runAnalyzerChecks`' method gate: an unlisted method is refused. */
inline void fluid_check_method(const std::string& method) {
    const std::vector<std::string> valid = fluid_list_valid_methods();
    if (std::find(valid.begin(), valid.end(), method) != valid.end()) return;
    throw UnsupportedError("SolverFLD: the '" + method + "' method is unsupported by this solver");
}


namespace detail {

/** `sn.nodetype == NodeType.Cache`, which is what selects `rmf`. */
template <class T>
bool fluid_has_cache(const qn::NetworkStruct<T>& sn) {
    for (const qn::NodeDef& nd : sn.nodes)
        if (nd.nodetype == qn::NodeType::Cache) return true;
    return false;
}

/** A DPS station, which the matrix method cannot express. */
template <class T>
bool fluid_has_dps(const qn::NetworkStruct<T>& sn) {
    for (const auto& st : sn.stations)
        if (st.sched == lang::SchedStrategy::DPS) return true;
    return false;
}

/**
 * The method name without its `fluid` prefix, which names the same method.
 *
 * `butools` and `aoi` are folded onto `mfq` here for the same reason: they name
 * the MFQ branch's BACKEND and its age-of-information reading rather than a
 * method of their own, which is how native python spells them and how
 * `@@SolverFLD/getAvgAoI` reaches it (it requires method='mfq').
 */
inline std::string fluid_unqualify(const std::string& method) {
    std::string m = method;
    if (m.size() > 6 && m.compare(0, 6, "fluid.") == 0) m = m.substr(6);
    if (m == "butools" || m == "aoi") return "mfq";
    return m;
}

/**
 * Port of `fluid_minnormal_applicable.m`: whether the moment closure can answer
 * this model, which is what `default` consults before preferring it.
 *
 * THE TEST IS STATIC. It inspects the model and the options, never the solution,
 * so a feature the model DECLARES is decided here and only here. The one
 * condition that cannot be static is a NON-HYPERBOLIC fixed point (balanced
 * bottlenecks, a saturated multiclass station, an overloaded open station): it
 * exists only once the mean is solved, `fluid_lyapunov` raises
 * FluidNonHyperbolicError on it, and the runner switches `minnormal` -- resolved
 * or requested -- to a first-order method on that exception alone, since the
 * closure has no stationary covariance either way and the mean is still
 * available.
 *
 * Every condition below mirrors a refusal that `fluid_moment_terms` or
 * `solver_fluid_moments` would otherwise raise, so this function and those
 * refusals must move together.
 */
template <class T>
bool fluid_minnormal_applicable(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                                std::string& reason) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const FluidLayout L = fluid_layout(sn);

    // Open and mixed models are supported: the covariance projects the EXT source
    // pool out. What it cannot take is a NON-POISSON arrival stream, whose source
    // coordinates track the phase of a single arrival process, not a population.
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].sched != lang::SchedStrategy::EXT) continue;
        for (std::size_t r = 0; r < K; ++r)
            if (L.kic[i][r] > 1) {
                reason = "class " + std::to_string(r + 1) + " has a " +
                         std::to_string(L.kic[i][r]) + "-phase (non-Poisson) arrival process";
                return false;
            }
    }
    // A cache model is answered through the decomposition analyzer, with the
    // closure in its network step, so cache nodes no longer decline the method.
    // What still declines is a replacement strategy with no drift-based fluid
    // model, mirroring the runtime gate in `fluid_cacheqn_gate`: LRU, HLRU,
    // CLIMB and QLRU are answered by a characteristic-time fixed point, which is
    // not a fluid method and has no covariance.
    for (std::size_t nd = 0; nd < sn.nodes.size(); ++nd) {
        if (sn.nodes[nd].nodetype != qn::NodeType::Cache) continue;
        const std::size_t key = nd + 1;  // nodeparam is keyed 1-based
        if (sn.nodeparam.count(key) == 0) {
            reason = "a cache uses a replacement strategy with no drift-based fluid model";
            return false;
        }
        const lang::ReplacementStrategy rs = sn.nodeparam.at(key).replacestrat;
        if (rs != lang::ReplacementStrategy::RR && rs != lang::ReplacementStrategy::FIFO &&
            rs != lang::ReplacementStrategy::SFIFO) {
            reason = "a cache uses a replacement strategy with no drift-based fluid model";
            return false;
        }
    }
    // The disciplines with a branch in `ode_rates_closing_factors`; anything else
    // falls through to g = x, i.e. an infinite server.
    for (std::size_t i = 0; i < M; ++i) {
        const lang::SchedStrategy s = sn.stations[i].sched;
        const bool ok = s == lang::SchedStrategy::INF || s == lang::SchedStrategy::EXT ||
                        s == lang::SchedStrategy::PS || s == lang::SchedStrategy::FCFS ||
                        s == lang::SchedStrategy::DPS || s == lang::SchedStrategy::GPS;
        if (!ok) {
            reason = "station " + std::to_string(i + 1) + " uses " + lang::sched_to_text(s) +
                     ", which has no fluid drift branch";
            return false;
        }
    }
    // The Lyapunov solve is cubic in the phase-resolved state, so the same cap
    // `solver_fluid_moments` enforces decides selection rather than being hit
    // later.
    if (L.nstates > opt.moment_maxstate) {
        reason = "the phase-resolved state has " + std::to_string(L.nstates) +
                 " coordinates, above the moment_maxstate limit of " +
                 std::to_string(opt.moment_maxstate);
        return false;
    }
    // The moment methods need an autonomous drift.
    // `fluid_minnormal_applicable.m:99`: a time-varying rate multiplier makes the
    // drift non-autonomous, so there is no stationary covariance to solve for and
    // `default` must resolve to a first-order method instead.
    if (fluid_has_time_varying_rates(opt)) {
        reason = "a rate schedule (nhpp_sched / rate_traj / rate_sched) makes the drift "
                 "time-varying";
        return false;
    }
    reason.clear();
    return true;
}

/**
 * Port of `fluid_dae_applicable.m`: whether the differential-algebraic route can
 * answer this model, which is what the fallback ladder consults before trying
 * `dae` on a declined `minnormal`.
 *
 * WHY THIS EXISTS SEPARATELY FROM `fluid_minnormal_applicable`. The two methods
 * state the SAME closure and differ only in how the coupled equations are
 * discharged, so a model `minnormal` accepts is almost always one `dae` accepts
 * too. Almost: `dae` carries a finite-difference Jacobian over the whole unknown
 * vector rather than one Lyapunov solve, so its state cap is lower; it closes on
 * the per-station variance only, so DPS and GPS are out; and it has no
 * decomposition route, so a cache model is out. Those three are exactly the
 * difference set, and naming them here keeps the ladder from entering a rung
 * that would refuse the model a moment later.
 *
 * The test is STATIC, for the same reason the min-normal one is. The one
 * condition that cannot be static is the NON-HYPERBOLIC fixed point the ladder
 * exists to route around -- it exists only once the mean is solved -- and `dae`
 * fails on it loudly, which is what moves the ladder to its last rung.
 *
 * Every condition below mirrors a refusal `solver_fluid_dae` would otherwise
 * raise, so this function and those refusals must move together.
 */
template <class T>
bool fluid_dae_applicable(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                          std::string& reason,
                          const FluidDaeOptions& dopt_in = FluidDaeOptions()) {
    const FluidDaeOptions dopt = fluid_dae_options(opt, dopt_in);
    const std::size_t M = sn.nstations;
    const FluidLayout L = fluid_layout(sn);

    // A cache model is a DECOMPOSITION, and `dae` has no arm for it: `minnormal`
    // on a cache model routes through `solver_fld_cacheqn_analyzer` with the
    // closure in its network step, while this route would read the cache nodes as
    // ordinary stations -- a decomposition has no single drift for the algebraic
    // constraint to attach to.
    if (fluid_has_cache(sn)) {
        reason = "a cache model is answered by the decomposition analyzer, which has no dae route";
        return false;
    }
    // The DPS and GPS shares close on the covariance BETWEEN station coordinates,
    // not on the station total, so their closure state is a matrix block rather
    // than the scalar the Newton vector carries. Mirrors the `solver_fluid_dae`
    // refusal.
    for (std::size_t i = 0; i < M; ++i) {
        const lang::SchedStrategy s = sn.stations[i].sched;
        if (s == lang::SchedStrategy::DPS || s == lang::SchedStrategy::GPS) {
            reason = "station " + std::to_string(i + 1) + " uses " + lang::sched_to_text(s) +
                     ", whose share closes on the covariance between its class coordinates "
                     "rather than on the station variance";
            return false;
        }
    }
    // The simultaneous solve is quartic overall, against the cubic of one Lyapunov
    // solve, so its crossover is lower than the 200 `minnormal` permits and it
    // carries its own cap. Same count the min-normal test forms, so the two limits
    // are read on the same scale.
    if (L.nstates > dopt.maxstate) {
        reason = "the phase-resolved state has " + std::to_string(L.nstates) +
                 " coordinates, above the dae_maxstate limit of " +
                 std::to_string(dopt.maxstate);
        return false;
    }
    reason.clear();
    return true;
}

/**
 * Port of `fluid_resolve_default_method.m`: what `default` stands for on THIS
 * model.
 *
 * Preference order: `rmf` for a cache model, then `minnormal` wherever
 * `fluid_minnormal_applicable` accepts the model, then the historical choice of
 * `closing` for DPS and `matrix` otherwise. The second-order closure dominates
 * the first-order methods on every family measured against the exact CTMC and is
 * the only method that can represent GPS at all, so it is preferred wherever it
 * applies.
 *
 * THIS RUNS BEFORE THE FEATURE GATE, NOT INSIDE THE DISPATCH SWITCH. The gate
 * validates `fluid_feature_set(method)`, and `default` is not `minnormal`, so
 * resolving later would let the gate reject a GPS model that the resolved method
 * supports. Keeping the decision in one function is what stops the gate and the
 * dispatch from disagreeing.
 */
template <class T>
std::string fluid_resolve_method(const qn::NetworkStruct<T>& sn, const std::string& method,
                                 const FluidOptions& opt, std::string& reason) {
    const std::string m = fluid_unqualify(method);
    reason.clear();
    // `rmf` NAMES THE CACHE DECOMPOSITION, and a model with no Cache node has
    // nothing to decompose: the reference's solver_fld_cacheqn_analyzer has no
    // cache-node test and its netsolve loop degenerates to the plain matrix
    // solve, so `rmf` is an ordinary fluid solve there. This port's analyzer
    // refuses an empty cache set -- rightly, it is the transient entry point's
    // precondition -- so resolve the name here instead of erroring on a model
    // the other three codebases answer. Resolving (rather than branching at the
    // dispatch) also keeps solver_fluid_export_odes coherent: without a cache
    // there IS one ODE system to export.
    if (m == "rmf" && !fluid_has_cache(sn)) {
        reason = "the model has no cache node, so the refined mean field has nothing to decompose";
        return "matrix";
    }
    if (m != "default") return m;
    if (fluid_has_cache(sn)) {
        reason = "the model has cache nodes";
        return "rmf";
    }
    // A BINDING BUFFER OR A CAPACITY REGION ALSO HAS ONE FLUID ROUTE, for the
    // same reason a Petri net does: nothing else in the fluid tree reads the cap,
    // the per-class cap or the region limit, so every other method integrates the
    // capped station as an unbounded one -- which is why
    // `fluid_check_finite_capacity` refuses them. Resolving `default` to one of
    // those turned a model this solver CAN answer into an error whose advice was
    // to type the very method the resolution should have picked. The capacity
    // test is the gate's own (`binding_capacity_reason`), so the two cannot
    // disagree, and where the DAE route declines the fall-through leaves the
    // refusal to the gate, which names the blocking feature.
    if (!sn.regions.empty() || qn::has_binding_capacity(sn)) {
        std::string dae_why;
        if (fluid_dae_applicable(sn, opt, dae_why)) {
            reason = sn.regions.empty() ? "the model has a binding finite buffer"
                                        : "the model has a finite capacity region";
            return "dae";
        }
    }
    if (fluid_minnormal_applicable(sn, opt, reason)) return "minnormal";
    return fluid_has_dps(sn) ? "closing" : "matrix";
}

/** The same, for a caller with no interest in why `minnormal` was declined. */
template <class T>
std::string fluid_resolve_method(const qn::NetworkStruct<T>& sn, const std::string& method,
                                 const FluidOptions& opt = FluidOptions()) {
    std::string reason;
    return fluid_resolve_method(sn, method, opt, reason);
}

/**
 * `@@SolverFLD/runAnalyzer.m:169-175`: limited load dependence is honoured only by
 * the closing family, whose drift multiplies the scheduling share by alpha(n_i).
 *
 * The C++ featset DOES emit LoadDependence from the model (unlike MATLAB's
 * getUsedLangFeatures, which is why the reference has to refuse here explicitly),
 * so this is the same refusal reached one step earlier -- kept as its own check
 * because `solver_fluid` is callable without the gate.
 */
template <class T>
void fluid_check_load_dependence(const qn::NetworkStruct<T>& sn, const std::string& m) {
    if (m == "closing" || m == "minnormal" || m == "refined" || m == "dae") return;
    for (std::size_t i = 0; i < sn.stations.size(); ++i)
        for (std::size_t k = 0; k < sn.stations[i].lldscaling.size(); ++k)
            if (std::fabs(num_traits<T>::to_double(sn.stations[i].lldscaling[k]) - 1.0) >
                lang::GlobalConstants::Zero)
                throw UnsupportedError(
                    "solver_fluid_run_analyzer: this model uses load dependence (setLoadDependence), which "
                    "the '" +
                    m +
                    "' method does not evaluate. Use method 'closing' for the mean-field answer, or "
                    "'minnormal'/'refined' for the moment-closure correction");
}

/**
 * Can `m` run the fluid fork-join fixed point on this model?
 *
 * A fork-join model is not integrated as one drift: `fj_fork_join_transform`
 * replaces the fork by auxiliary classes and the answer is the fixed point of
 * solving that transformed model repeatedly. On a CLOSED model the transform
 * stays closed and every fluid method takes it. On an OPEN one the auxiliary
 * classes arrive at a Source, and the DAE form has no unknowns for them: the
 * inner solve fails on the class count rather than returning a drift, so the
 * method is refused by name instead.
 *
 * `refined` is NOT listed here even though it fails the same way, because it is
 * already refused on every open model, fork-join or not, by its own closed-model
 * restriction (`qn::fluid_feature_set("refined")`).
 *
 * ASKED AS A PREDICATE so that a REPORT can reach the rule too: Fork and
 * OpenClass are both declared names, so a feature set cannot state a rule that
 * is their CONJUNCTION.
 *
 * @param sn the refreshed struct of the model
 * @param m the resolved method name
 * @return an empty string when the method may run the fixed point here
 */
template <class T>
std::string fluid_forkjoin_supports(const qn::NetworkStruct<T>& sn, const std::string& m) {
    if (m != "dae" && m != "fluid.dae") return std::string();
    if (!sn.has_fork()) return std::string();
    bool any_open = false;
    for (std::size_t r = 0; r < sn.classes.size(); ++r)
        if (!std::isfinite(sn.classes[r].population)) any_open = true;
    if (!any_open) return std::string();
    return "solver_fluid_run_analyzer: the dae method has no route through the fork-join fixed "
           "point on an OPEN model: the transform hands the inner solve a mixed network whose "
           "auxiliary open classes the DAE form carries no unknowns for. Use method 'minnormal', "
           "which is the same closure and does run that fixed point";
}

/**
 * `@@SolverFLD/runAnalyzer.m:159-183`: the two limits of the DAE route that the
 * feature gate cannot express, refused before the Newton solve so the message
 * names the model feature rather than surfacing from inside it.
 *
 * A CACHE MODEL is a decomposition, not one system: the caches are solved in
 * isolation and the network with them relabeled as class switches, so there is
 * no single drift for the constraint to be attached to. `minnormal` and `rmf`
 * reach `solver_fld_cacheqn_analyzer`, which carries the closure inside its
 * network step; no such route exists for the DAE form. Cache IS declared in the
 * fluid featset -- for those two methods -- so the gate passes it and this is
 * the only place the distinction can be drawn.
 *
 * DPS and GPS are NOT repeated here: the per-method featset already unsets both
 * for `dae` (solver_feature_sets.h), so the gate refuses them one step earlier
 * and by feature name, which is the finer message of the two.
 */
template <class T>
void fluid_check_dae(const qn::NetworkStruct<T>& sn, const std::string& m) {
    if (m != "dae") return;
    if (fluid_has_cache(sn))
        throw UnsupportedError(
            "solver_fluid_run_analyzer: the dae method does not support caching stations: a cache model is "
            "solved by decomposition, so it has no single drift to constrain. Use "
            "method 'minnormal' for the same closure, or 'rmf'");
    const std::string fj = fluid_forkjoin_supports(sn, m);
    if (!fj.empty()) throw UnsupportedError(fj);
}

/**
 * The structural finite-buffer gate `@@SolverFLD/runAnalyzer.m` applies, which
 * this port did not carry: nothing in the fluid tree reads `sn.cap` or
 * `sn.classcap`, so a capped station was integrated as an unbounded one and the
 * table reported more jobs in the buffer than the buffer holds.
 *
 * `dae` carries the buffer as an algebraic constraint on the drift, and `mol` is
 * stated for the Mt/G/s/0 LOSS system, where the server count IS the buffer: for
 * those two a finite capacity is the model rather than something ignored. Every
 * other method keeps the guard. There is no feature-registry name for plain
 * capacity, hence the structural test -- `SolverMVA`, `SolverNC` and `SolverAG`
 * gate the same way, through the same `check_binding_capacity`.
 *
 * `m` is the RESOLVED method, as the reference's gate receives it.
 */
template <class T>
void fluid_check_finite_capacity(const qn::NetworkStruct<T>& sn, const std::string& m) {
    if (m == "dae" || m == "fluid.dae" || m == "mol" || m == "fluid.mol") return;
    qn::check_binding_capacity("SolverFLD", sn);
}

}  // namespace detail

/**
 * Port of `@@SolverFLD/runAnalyzer.m`: resolve the method, route to the function
 * the reference routes to, correct once.
 *
 * This is what a host should call. `solver_fluid` remains the port of
 * `solver_fluid_analyzer.m` alone and stays reachable for a caller that wants
 * one named method and no model-dependent resolution.
 *
 * `sn_out`, when given, receives `result.solverSpecific.sn`: the struct the
 * returned table was integrated on, which the FCFS refit may have re-fitted to
 * a different phase count. Anything that reads `FluidSolution::xvec` afterwards
 * must be handed THAT struct, because the vector is laid out by its phases.
 *
 * `refreshed_out`, when given, receives the cache decomposition's own struct:
 * the base one with the cache self-switch renormalized at the CONVERGED hit and
 * miss split. It is DELIBERATELY NOT `sn_out`: that one is the inner, relabeled
 * struct the ODE ran on and it lays out `xvec`, whereas this one is the caller's
 * own topology carrying a solved quantity. Left untouched on every model with no
 * cache, and on the cache branches that report no split.
 *
 * `cache_out`, when given, receives the SAME split as a `CacheMetrics`, which is
 * the form every other solver states it in (`mva::AvgResult::cache`). The
 * refreshed struct alone is not enough for a caller: a host that solves through
 * the CLI reads the hit and miss fractions back onto its own Cache node from
 * that block, and MATLAB's `CPPLINE.restoreCacheResults` CLEARS the node when it
 * is absent -- which is how `-s fluid` came to report link()'s offered 1/2-1/2
 * on cache_replc_rr for a split it had converged to 0.5623/0.4377.
 */
template <class T>
FluidSolution solver_fluid_run_analyzer(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                               qn::NetworkStruct<T>* sn_out = nullptr,
                               qn::NetworkStruct<T>* refreshed_out = nullptr,
                               solvers::CacheMetrics<T>* cache_out = nullptr);

namespace detail {

/**
 * `@@SolverFLD/fldDispatch.m`: one inner solve of the fixed point, as an
 * `mva::MvaSolution`.
 *
 * The driver is `mva::fj_fixed_point`, the one MVA and NC drive; it speaks
 * MvaSolution, so the fluid table is carried into that shape here rather than
 * a second fork-join loop being written for the fluid solver. Nothing else is
 * fluid-specific: the transform emits Source, Delay, Queue, Router and
 * ClassSwitch, all of which the drift already carries.
 */
template <class T>
mva::MvaSolution<T> fluid_as_mva_solution(const FluidSolution& f) {
    mva::MvaSolution<T> s;
    const std::size_t M = f.QN.rows(), K = f.QN.cols();
    const T zero = num_traits<T>::from_int(0);
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            s.Q(i, k) = num_traits<T>::from_double(f.QN(i, k));
            s.U(i, k) = num_traits<T>::from_double(f.UN(i, k));
            s.R(i, k) = num_traits<T>::from_double(f.RN(i, k));
            s.Tp(i, k) = num_traits<T>::from_double(f.TN(i, k));
        }
    s.C.reserve(f.CN.size());
    for (double v : f.CN) s.C.push_back(num_traits<T>::from_double(v));
    s.X.reserve(f.XN.size());
    for (double v : f.XN) s.X.push_back(num_traits<T>::from_double(v));
    s.method = f.method;
    s.iter = static_cast<int>(f.iters);
    // The fluid solver produces no normalizing constant; the reference leaves it
    // at zero rather than fabricating one.
    s.lG = 0.0;
    return s;
}

/** The merged fork-join table, back in the fluid solver's own result shape. */
template <class T>
FluidSolution fluid_from_mva_solution(const mva::MvaSolution<T>& s, const std::string& method) {
    FluidSolution f;
    const std::size_t M = s.Q.rows(), K = s.Q.cols();
    f.QN = Matrix<double>(M, K, 0.0);
    f.UN = Matrix<double>(M, K, 0.0);
    f.RN = Matrix<double>(M, K, 0.0);
    f.TN = Matrix<double>(M, K, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            f.QN(i, k) = num_traits<T>::to_double(s.Q(i, k));
            f.UN(i, k) = num_traits<T>::to_double(s.U(i, k));
            f.RN(i, k) = num_traits<T>::to_double(s.R(i, k));
            f.TN(i, k) = num_traits<T>::to_double(s.Tp(i, k));
        }
    f.CN.reserve(s.C.size());
    for (const T& v : s.C) f.CN.push_back(num_traits<T>::to_double(v));
    f.XN.reserve(s.X.size());
    for (const T& v : s.X) f.XN.push_back(num_traits<T>::to_double(v));
    f.iters = static_cast<std::size_t>(s.iter < 0 ? 0 : s.iter);
    f.method = method;
    return f;
}

/** The leading `M` station rows of a metric matrix. */
inline Matrix<double> fluid_leading_rows(const Matrix<double>& m, std::size_t M) {
    if (m.rows() <= M) return m;
    Matrix<double> out(M, m.cols(), 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < m.cols(); ++k) out(i, k) = m(i, k);
    return out;
}

/**
 * The fork-join arm of `solver_fluid_run_analyzer`.
 *
 * Only the steady-state means are produced. The transient tables are indexed by
 * the ORIGINAL stations and classes, whereas each pass integrates a DIFFERENT
 * transformed network, so a trajectory read off the last pass would not be the
 * trajectory of the model the caller built; `xvec` is left empty for the same
 * reason, which is what makes the passage-time and transient entry points
 * decline a fork-join model rather than answer it from the wrong state.
 */
template <class T>
FluidSolution fluid_fork_join_run(const qn::NetworkStruct<T>& sn, const FluidOptions& opt) {
    mva::FjMmt<T> tr = mva::fj_fork_join_transform(sn, opt.fork_join);
    std::vector<T> lam(tr.V.classes.size() + 1,
                       num_traits<T>::from_double(lang::GlobalConstants::FineTol));
    mva::MvaOptions mopt;
    mopt.method = opt.method;
    mopt.tol = opt.tol;
    mopt.iter_tol = opt.iter_tol;
    mopt.iter_max = opt.iter_max;
    mopt.fork_join = opt.fork_join;
    mopt.base_has_fork = true;
    FluidOptions inner = opt;
    const mva::MvaSolution<T> merged =
        mva::fj_fixed_point(sn, tr, lam, mopt, [&inner](qn::NetworkStruct<T>& V) {
            return fluid_as_mva_solution<T>(solver_fluid_run_analyzer(V, inner));
        });
    FluidSolution out = fluid_from_mva_solution<T>(merged, opt.method);
    // The driver drops the auxiliary CLASS columns but not the transform's own
    // STATION row: it appends a Source to carry the auxiliary arrivals, and the
    // base stations are the prefix of that list. `solver_mva_run_analyzer` trims it in
    // `filter_metric` against the base struct; this entry point returns the
    // table directly, so it trims here.
    out.QN = detail::fluid_leading_rows(out.QN, sn.nstations);
    out.UN = detail::fluid_leading_rows(out.UN, sn.nstations);
    out.RN = detail::fluid_leading_rows(out.RN, sn.nstations);
    out.TN = detail::fluid_leading_rows(out.TN, sn.nstations);
    return out;
}

}  // namespace detail

namespace detail {

/** Whether the model holds a Transition node, i.e. is a stochastic Petri net. */
template <class T>
inline bool fluid_is_petri_net(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].nodetype == lang::NodeType::Transition) return true;
    return false;
}

/**
 * The 'dae' method's Petri arm.
 *
 * 'dae' is the ONLY fluid method that can carry a net: the P-invariants, the
 * firing flow of an immediate transition and a bounded place are all EQUATIONS,
 * and the other methods have nowhere to put them. A method named explicitly is
 * therefore checked rather than silently redirected, so a caller who asked for
 * 'closing' on a net is told why it cannot answer.
 */
template <class T>
inline FluidSolution fluid_petri_run(const qn::NetworkStruct<T>& sn, const FluidOptions& opt) {
    const std::string m = (opt.method.rfind("fluid.", 0) == 0) ? opt.method.substr(6) : opt.method;
    if (!(m == "default" || m == "dae"))
        throw UnsupportedError(
            "SolverFluid: method '" + opt.method +
            "' cannot solve a Petri net: its conserved quantities are P-invariants rather than "
            "chain populations, an immediate transition is an algebraic FLOW rather than an "
            "event with a rate, and a bounded place is a linear inequality on the marking. Only "
            "'dae' states those as equations; every other fluid method builds its drift from the "
            "station/class/phase encoding, where a Place contributes no coordinate at all, and "
            "would integrate the net as an empty model and report zeros without a warning");
    const petri::PetriSolution ps = petri::solver_fluid_petri(sn, petri::PetriOptions());
    FluidSolution out;
    out.QN = ps.QN;
    out.UN = ps.UN;
    out.RN = ps.RN;
    out.TN = ps.TN;
    out.xvec = ps.x;
    out.iters = ps.iters;
    out.method = "dae";
    const std::size_t M = ps.QN.rows(), K = ps.QN.cols();
    out.CN.assign(K, 0.0);
    out.XN.assign(K, 0.0);
    for (std::size_t k = 0; k < K; ++k) {
        double q = 0.0, x = 0.0;
        for (std::size_t i = 0; i < M; ++i) {
            q += ps.QN(i, k);
            x = std::max(x, ps.TN(i, k));
        }
        out.XN[k] = x;
        out.CN[k] = (x > 1e-14) ? q / x : 0.0;
    }
    // The marking covariance of the linear noise approximation: the second
    // moment no other solver of a Petri net in LINE returns.
    out.has_moments = true;
    out.moments.Sigma = ps.Sigma;
    out.moments.QVar = ps.QVar;
    out.moments.QStd = ps.QStd;
    return out;
}

}  // namespace detail

template <class T>
FluidSolution solver_fluid_run_analyzer(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                               qn::NetworkStruct<T>* sn_out,
                               qn::NetworkStruct<T>* refreshed_out,
                               solvers::CacheMetrics<T>* cache_out) {
    std::string why;
    if (sn_out) *sn_out = sn;
    // A MODEL HOLDING ANY Transition NODE IS A DIFFERENT FORMALISM and goes to a
    // dedicated runner rather than through the queueing dispatch. The queueing
    // path rewrites UN from sn.rates (NaN at a Place) and RN from the station
    // scheduling, and would overwrite the Petri conventions the Petri route
    // reports -- a place's utilization IS its token count, and its throughput is
    // the token departure rate, not a service completion rate. The branch is
    // taken before the feature gate, which is written about queueing stations.
    if (detail::fluid_is_petri_net(sn)) return detail::fluid_petri_run(sn, opt);
    // runAnalyzerChecks' METHOD gate, before the feature gate and before the
    // resolution, exactly as the reference orders them: an unknown name is
    // reported as unknown rather than as whatever the feature gate finds.
    fluid_check_method(opt.method);
    const std::string m = detail::fluid_resolve_method(sn, opt.method, opt, why);
    // runAnalyzerChecks' universal feature gate, on the RESOLVED method. Only
    // the runner is gated: `solver_fluid` stays the ungated port of
    // solver_fluid_analyzer.m for a caller that wants one named method.
    qn::feature_gate("SolverFluid", qn::fluid_feature_set(m), sn);
    detail::fluid_check_load_dependence(sn, m);
    detail::fluid_check_dae(sn, m);
    detail::fluid_check_finite_capacity(sn, m);
    // Fork-join: the same fixed point MVA and NC drive, with a fluid inner
    // solve. Taken after the gate so an unsupported feature is still named by
    // its own message rather than by a failure inside the transformed model.
    if (sn.has_fork()) {
        FluidOptions fo = opt;
        fo.method = m;
        return detail::fluid_fork_join_run(sn, fo);
    }
    FluidOptions o = opt;
    o.method = m;

    // The reference's FCFS non-exponential refit loop, on the branches that are
    // routed AROUND `solver_fluid` (which runs it itself). `solve` is the same
    // per-method function the branch used for its first integration, so the
    // sweeps stay on one drift; the loop is a no-op for a method the reference
    // does not refit and for a model with no FCFS station.
    const auto closing_solve = [](const qn::NetworkStruct<T>& s, const FluidOptions& op) {
        return solver_fluid_closing(s, op);
    };
    const auto moments_solve = [](const qn::NetworkStruct<T>& s, const FluidOptions& op) {
        return solver_fluid_moments(s, op);
    };
    const auto dae_solve = [](const qn::NetworkStruct<T>& s, const FluidOptions& op) {
        return solver_fluid_dae(s, op);
    };

    FluidSolution out;
    if (m == "rmf" || (m == "minnormal" && detail::fluid_has_cache(sn))) {
        // A cache model is a DECOMPOSITION, not one ODE: the caches are solved
        // in isolation and the network with them relabeled as class switches.
        // "minnormal" takes the same route, with the closure inside `netfun`.
        // The analyzer's refreshed struct is CARRIED, not dropped: the split it
        // converged on is what every visit-weighted quantity downstream needs.
        const FluidCacheqnSolution<T> cq = solver_fld_cacheqn_analyzer(sn, o);
        out = cq.sol;
        if (refreshed_out) *refreshed_out = cq.refreshed;
        if (cache_out) *cache_out = solvers::cache_metrics_of_matrix(sn, cq.hitprob, cq.missprob);
    } else if (m == "closing" || m == "statedep" || m == "softmin" || m == "tbi") {
        out = detail::fluid_fcfs_nonexp_refit(sn, o, solver_fluid_closing(sn, o), closing_solve,
                                              sn_out);
    } else if (m == "minnormal" || m == "refined") {
        // A non-hyperbolic fixed point cannot be seen before the mean is solved,
        // so it cannot be declined in advance.
        //
        // THE LADDER HAS TWO RUNGS, AND THE FIRST ONE KEEPS THE CLOSURE. Most of
        // these failures are not a property of the model at all:
        // `solver_fluid_moments` must start its alternation at sigma2 = 0, where
        // min(n,c) has no derivative, so a saturated or balanced model's
        // first-order fixed point lands on the kink, sits on a continuum of
        // equilibria, and the Jacobian there is neutral. `solver_fluid_dae` seeds
        // the variance POSITIVE and never adopts sigma2 = 0 as an iterate, so the
        // smoothed E[min(X,c)] breaks the degeneracy and the fixed point is
        // isolated and hyperbolic -- it answers the same closure, with a
        // covariance, where the alternation cannot. Dropping straight to first
        // order is not merely a lost second moment: on a balanced two-station PS
        // cycle at N=10 it returns [9 1] against the exact [5 5], because a
        // first-order method has no reason to prefer one point of the continuum
        // over another.
        //
        // The second rung is the first-order method, taken when `dae` declines the
        // model in advance (`fluid_dae_applicable`) or fails on the same error,
        // which is the genuinely non-hyperbolic case: an unstable open station has
        // no stationary distribution to approximate under any closure. Fall back
        // whether the closure was RESOLVED from `default` or REQUESTED outright:
        // it has no stationary covariance either way, so refusing an explicit
        // request would only deny the caller the mean that is still available.
        try {
            out = detail::fluid_fcfs_nonexp_refit(sn, o, solver_fluid_moments(sn, o),
                                                  moments_solve, sn_out);
        } catch (const FluidNonHyperbolicError&) {
            std::string dae_reason;
            bool solved = false;
            if (detail::fluid_dae_applicable(sn, o, dae_reason)) {
                FluidOptions daeopt = o;
                daeopt.method = "dae";
                try {
                    out = detail::fluid_fcfs_nonexp_refit(sn, daeopt, solver_fluid_dae(sn, daeopt),
                                                          dae_solve, sn_out);
                    o = daeopt;
                    solved = true;
                } catch (const FluidNonHyperbolicError&) {
                    // genuinely neutral or unstable: the last rung answers it
                }
            }
            if (!solved) {
                o.method = detail::fluid_has_dps(sn) ? "closing" : "matrix";
                qn::feature_gate("SolverFluid", qn::fluid_feature_set(o.method), sn);
                detail::fluid_check_load_dependence(sn, o.method);
                // `solver_fluid` is the analyzer, refit included, so the matrix arm
                // is complete; the closing arm is routed around it and runs the loop
                // here.
                if (o.method != "closing") return solver_fluid(sn, o, sn_out);
                out = detail::fluid_fcfs_nonexp_refit(sn, o, solver_fluid_closing(sn, o),
                                                      closing_solve, sn_out);
                detail::fluid_analyzer_correct(sn, out.QN, out.UN, out.RN, out.TN);
                detail::fluid_snap_all(out.QN, out.UN, out.RN, out.TN);
                return out;
            }
        }
    } else if (m == "dae") {
        // The same closure `minnormal` computes, stated and solved as ONE
        // differential-algebraic system: population conservation is an equation
        // rather than a consequence of the drift, and a region cap is an
        // algebraic inequality beside it. No non-hyperbolic fallback, unlike the
        // branch above: the constraint rows supply exactly the rank the drift
        // Jacobian is missing along the conserved directions, so the Lyapunov
        // step that fails there is not what this route solves through. Refusing
        // rather than falling back is also the only honest answer once a REGION
        // is present -- the first-order methods would return the unconstrained
        // population, which is a different model, not a coarser answer to this
        // one.
        out = detail::fluid_fcfs_nonexp_refit(sn, o, solver_fluid_dae(sn, o), dae_solve, sn_out);
    } else if (m == "kp") {
        out = solver_fluid_kp(sn, o);
    } else if (detail::fluid_qsys_handles(m)) {
        // The single-station fluid limits are closed forms, not integrations of
        // the network drift: they take the whole model in one call, so they are
        // returned unmodified -- neither the FCFS refit nor the analyzer
        // correction below applies to a closed form.
        return solver_fluid_qsys(sn, o);
    } else {
        // `solver_fluid` is the switch, the refit AND the correction, so its
        // branches are complete already.
        return solver_fluid(sn, o, sn_out);
    }
    detail::fluid_analyzer_correct(sn, out.QN, out.UN, out.RN, out.TN);
    detail::fluid_snap_all(out.QN, out.UN, out.RN, out.TN);
    return out;
}

/**
 * `-a tran` / `@@SolverFLD/getTranAvg` with the method HONOURED, which is the
 * one place the reference does not force `closing`.
 *
 * `solver_fluid_tran_avg` integrates the first-order closing drift whatever
 * `opt.method` says, and that IS the reference's rule for every steady-state
 * device -- matrix, pnorm, softmin and the rest have no trajectory of their own,
 * so the reference substitutes one and warns. `dae` is the exception
 * (`getTranAvg.m:155-166`): it integrates the closure itself, with conservation
 * as an algebraic equation, so substituting the first-order drift for it would
 * hand back a DIFFERENT method's trajectory under this method's name.
 *
 * Every other method keeps the byte-identical path it had, gate included -- this
 * routes, it does not re-gate.
 */
template <class T>
std::vector<FluidTranPoint> solver_fluid_run_transient(const qn::NetworkStruct<T>& sn,
                                                       const FluidOptions& opt,
                                                       std::size_t points = 101) {
    // The single-station fluid limits produce their OWN trajectory over the
    // horizon -- that is what they are -- so substituting the first-order
    // closing drift for them would answer a different model with a
    // time-averaged rate.
    if (detail::fluid_qsys_handles(detail::fluid_unqualify(opt.method))) {
        FluidOptions o = opt;
        o.method = detail::fluid_unqualify(opt.method);
        qn::feature_gate("SolverFluid", qn::fluid_feature_set(o.method), sn);
        std::vector<FluidTranPoint> traj;
        solver_fluid_qsys(sn, o, &traj);
        return traj;
    }
    if (detail::fluid_unqualify(opt.method) != "dae")
        return solver_fluid_tran_avg(sn, opt, points);
    qn::feature_gate("SolverFluid", qn::fluid_feature_set("dae"), sn);
    detail::fluid_check_dae(sn, std::string("dae"));
    return solver_fluid_dae_transient(sn, opt, fluid_default_horizon(sn, opt), points);
}

/**
 * Port of `@@SolverFLD/getCdfRespT`: the response-time law of every (station,
 * class) pair, read off the marked-fluid trajectory started from the steady
 * state.
 *
 * THE SOLVE COMES FIRST AND IS NOT THE CALLER'S. The reference clears its
 * cached result, re-runs `getAvg` to fill `odeStateVec`, and hands the passage
 * time BOTH that vector and `result.solverSpecific.sn` -- the struct the FCFS
 * refit may have re-fitted -- because the vector is laid out by the phase counts
 * of that struct and not of the model's. Bundling the two here is what keeps a
 * caller from pairing a state vector with the wrong struct, which would be a
 * length error at best and a law of the wrong model at worst.
 *
 * A Source has no response time and is left empty, as is a pair whose class the
 * station does not serve; an empty entry is an ABSENT law, never a degenerate
 * one.
 */
template <class T>
std::vector<std::vector<FluidPassage> > solver_fluid_cdf_respt(const qn::NetworkStruct<T>& sn,
                                                               const FluidOptions& opt,
                                                               std::size_t points = 201) {
    qn::NetworkStruct<T> snr = sn;
    const FluidSolution r = solver_fluid_run_analyzer(sn, opt, &snr);
    if (r.xvec.empty())
        throw UnsupportedError(
            "solver_fluid_cdf_respt: the '" + r.method +
            "' method reports no fluid state vector to mark a job in, so it has no passage time; "
            "use 'closing', 'matrix' or a moment closure");
    std::vector<std::vector<FluidPassage> > RD(
        snr.nstations, std::vector<FluidPassage>(snr.nclasses));
    for (std::size_t i = 0; i < snr.nstations; ++i) {
        if (snr.stations[i].nodetype == qn::NodeType::Source) continue;
        for (std::size_t c = 0; c < snr.nclasses; ++c) {
            if (snr.disabled[i][c]) continue;
            RD[i][c] = fluid_passage_time(snr, r.xvec, i + 1, c + 1, opt.tol, points, r.closure);
        }
    }
    return RD;
}

/**
 * Port of `@@SolverFLD/exportODEs.m` at the runner's own method resolution, so
 * the exported system is the one `solver_fluid_run_analyzer` would integrate.
 *
 * `rmf` has no single exported drift -- the cache decomposition rewrites the
 * routing between passes, so there is one fluid system per sweep and none of
 * them is "the" model's -- and is refused by name rather than exported at the
 * routing of whichever sweep happened to be last.
 */
template <class T>
std::string solver_fluid_export_odes(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                                     const std::string& notation = "scalar",
                                     const std::string& model_name = "model") {
    const std::string m = detail::fluid_resolve_method(sn, opt.method);
    if (m == "rmf")
        throw UnsupportedError(
            "solver_fluid_export_odes: a Cache model is solved by alternating a cache drift with a "
            "queueing drift whose routing is rewritten at every sweep, so it has no one ODE system "
            "to export; export the model without its Cache node, or read the sweeps through "
            "solver_fld_cacheqn_tran");
    FluidOptions o = opt;
    o.method = m;
    return fluid_export_odes(sn, o, notation, model_name);
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_RUNNER_H
