/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_LN_SOLVER_LN_H
#define LINE_SOLVERS_LN_SOLVER_LN_H

/**
 * SolverLN: layered decomposition of a layered queueing network.
 *
 * Port of the matlab/src/solvers/LN SolverLN class folder -- construct, buildLayers,
 * buildLayersRecursive, init, initInterlock, converged, analyze, post,
 * updateMetrics, updatePopulations, updateThinkTimes, updateLayers,
 * updateRoutingProbabilities, getEntryServiceMatrix and getEnsembleAvg --
 * driven by the EnsembleSolver iteration in @@EnsembleSolver/iterate.m.
 *
 * The method is a Picard iteration over a fixed point. Each processor and each
 * software task becomes a closed queueing network (a "layer") in which the
 * element is the server and its callers are the customers; solving all layers
 * gives new residence times, which become the service times of the calls one
 * level up and the think times of the callers one level down, and the layers
 * are re-solved until the queue lengths stop moving.
 *
 * WHAT IS IMPLEMENTED. Every path a plain layered model takes is here:
 * synchronous calls, entry selection by throughput ratio, multi-entry tasks,
 * reference tasks with think time, infinite-server tasks and processors,
 * server replication, AND-forks and AND-joins, and the LQNS V5 interlock
 * analysis. Forwarding is solved too: construct() runs the
 * lqn_fwd_rendezvous rewrite of lqn_helpers.h (:306), which flattens every
 * forwarding chain into caller-side pseudo rendezvous arcs, so nothing below
 * this point needs a forwarding case. Asynchronous calls and entry open
 * arrivals are solved by giving the layer its own Source/Sink pair and an
 * open chain (build_layer's async_here/open_entries construction, roughly
 * :622-668). Cache tasks and item entries are solved by a Cache node placed
 * in the host layer, with the hit and miss branches switching class at the
 * server (build_layer's iscachelayer path, :591-616 and the read walk at
 * :883-908). Setup tasks with a setup time are solved by charging the cold
 * start to the ENTRY, at the probability the thread was found powered down
 * (setup_charge, added to entry_servt and to the caller's think time); no layer
 * station carries a SetupDelayOffParam. Admission constraints are solved by a Region on
 * the server station, which forces that layer's fixed point through
 * SolverCTMC rather than MVA or fluid (:1051-1096, dispatch at :1753-1760).
 * reject_unsupported() (:433) is kept but empty: every construct the
 * reference LQN model can carry is now accounted for here. A few narrower
 * combinations remain unsupported and are refused where they are discovered
 * instead of up front: a fork sharing a layer with an open stream (:492-496),
 * a fork's fixed point outside layer_solver 'mva' (:1723-1728), and a fluid
 * layer or an AND-join's order-statistic fit under a non-double/non-real
 * arithmetic backend (:120-126, :2198-2201).
 *
 * WHAT ELSE THIS CLASS ANSWERS, beyond the mean-value table. `options.method`
 * selects between three DIFFERENT questions rather than three routes to one:
 * `default` is the mean-based update; `moment3` is updateMetricsMomentBased,
 * which fits an APH to each layer's response-time CDF and convolves along the
 * entry's activity sequence, so it additionally reports a per-entry response
 * time DISTRIBUTION through get_cdf_respt(); `mwba.upper` / `mwba.lower` report
 * Majumdar-Woodside robust box bounds and never solve a layer at all.
 * get_tran_avg() is the layered transient, decoupled (demands frozen at the
 * fixed point) or coupled (waveform relaxation over time-varying demands
 * injected through the fluid rate schedule). get_sensitivity_table() delegates
 * to each layer's own sensitivity table; every derivative in it is a partial
 * WITHIN its layer, not a total derivative of the layered model.
 *
 * THE LAYER ENGINE is `layer_solver`: `mva`, `nc`, `fluid` or `ssa`, the C++
 * spelling of the reference's solver FACTORY. They converge to DIFFERENT fixed
 * points, because each feeds different demands back into the next outer sweep.
 * `ssa` is noisy, so the deterministic convergence test is replaced by
 * LnStochController (Robbins-Monro relaxation, Polyak-Ruppert reporting);
 * running it under the deterministic test would simply reach iter_max.
 *
 * PHASE-2 ACTIVITIES ARE SOLVED, not refused. `servt` keeps both phases, so the
 * server's utilization is unchanged, while `residt` carries the CALLER's view:
 * phase 1 in full plus the phase-2 time the caller is actually overtaken by,
 * with the overtaking probability from lqn_analyzers.h's
 * lqn_overtake_prob_markov. Every phase-2 branch below is gated on `has_phase2`
 * and is inert on a single-phase model.
 *
 * ARITHMETIC. All model quantities are T. Iteration counters, populations,
 * multiplicities and the convergence tolerances are double, matching the
 * reference: they are properties of the algorithm and of the model's integer
 * structure, not quantities whose precision is under study.
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "line/api/fj/fj_branch_members.h"
#include "line/api/fj/fj_quorum_moments.h"
#include "line/api/lqn/lqn_boxbounds.h"
#include "line/api/mam/aph_convseq.h"
#include "line/api/mam/aph_fit.h"
#include "line/api/mam/aph_simplify.h"
#include "line/api/sn/sn_fj_nodevisits_mmt.h"
#include "line/api/sn/sn_compat_rate.h"
#include "line/lang/distribution.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/lang/qn/qn_layer.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/fluid/fluid_passage.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/solvers/mam/mam_dispatch.h"
#include "line/api/lqn/lqn_ph.h"
#include "line/solvers/ln/lqn_helpers.h"
#include "line/solvers/mva/fj_driver.h"
#include "line/solvers/mva/fj_mmt.h"
#include "line/solvers/mva/solver_mva.h"
#include "line/solvers/mva/mva_dispatch.h"
#include "line/solvers/nc/nc_dispatch.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/sens/solver_sens_table.h"
#include "line/solvers/ssa/ssa_dispatch.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace ln {

using lang::CallType;
using lang::Distrib;
using lang::GlobalConstants;
using lang::JobClassType;
using lang::LqnElement;
using lang::NodeType;
using lang::PrecedenceType;
using lang::RoutingStrategy;
using lang::SchedStrategy;
using lqn::LqnCallGroup;
using lqn::LqnStruct;

/** One (station, class) rate trajectory injected into a layer's closing ODE. */
using FluidRateSched = fluid::FluidOptions::RateSched;

namespace detail {

/**
 * Run the fluid analyzer on one layer and write its metrics into the shape the
 * MVA layer path returns.
 *
 * Overloaded rather than gated: LSODA is double-only, so the template below is
 * what any other backend resolves to, and it refuses by name instead of failing
 * to compile inside a branch that could never run.
 */
inline void ln_fluid_solve(const qn::NetworkStruct<double>& L, const fluid::FluidOptions& fo,
                           mva::MvaSolution<double>& out) {
    // THE RUNNER, NOT THE SWITCH. `solver_fluid` is the ungated port of
    // solver_fluid_analyzer.m and cannot reach `minnormal`, `rmf`, `refined` or
    // `kp` at all; `solver_fluid_run_analyzer` is the port of runAnalyzer's resolution,
    // which is what `LN(model, @(m) Fluid(m))` runs in the reference. Calling
    // the switch here silently ran the FIRST-ORDER `matrix` method in every
    // layer: on lqn_basic that reported T3's residence as its bare demand
    // (0.02 against MATLAB's 0.0206615) because the hard min charges no
    // queueing below the server count, and the LN table then disagreed with
    // MATLAB by 3.3% on a metric neither codebase flagged.
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(L, fo);
    for (std::size_t i = 0; i < L.nstations; ++i)
        for (std::size_t r = 0; r < L.nclasses; ++r) {
            out.Q(i, r) = s.QN(i, r);
            out.U(i, r) = s.UN(i, r);
            out.R(i, r) = s.RN(i, r);
            out.Tp(i, r) = s.TN(i, r);
        }
    for (std::size_t r = 0; r < L.nclasses && r < s.XN.size(); ++r) {
        out.X[r] = s.XN[r];
        out.C[r] = s.CN[r];
    }
    out.method = s.method;
    out.iter = static_cast<int>(s.iters);
}
template <class T>
void ln_fluid_solve(const qn::NetworkStruct<T>&, const fluid::FluidOptions&,
                    mva::MvaSolution<T>&) {
    throw UnsupportedError(
        "SolverLN: a fluid layer integrates its drift with LSODA, whose coefficients assume "
        "double precision; rerun with --arith double or use layer_solver 'mva'");
}

/**
 * The response-time CDF of every (station, class) of one layer, which is what
 * `SolverFluid(ensemble{e}).getCdfRespT` returns in the reference.
 *
 * The layer is integrated once to its fixed point and the converged state is
 * the marking the passage starts from, so the law reported is the STATIONARY
 * response time and not the one seen from an arbitrary initial condition. That
 * is `options.init_sol = odeStateVec` in `@@SolverFLD/getCdfRespT`.
 *
 * A pair the class never visits comes back as the degenerate curve at zero,
 * which `fluid_passage_time` already returns; the moment extraction above reads
 * it as a mean of zero and the caller drops the term, as MATLAB's
 * `m1 > CoarseTol` test does.
 */
inline std::vector<std::vector<fluid::FluidPassage>> ln_fluid_cdf_respt(
    const qn::NetworkStruct<double>& L, const fluid::FluidOptions& fo) {
    fluid::FluidOptions o = fo;
    // THE SWITCH HERE, DELIBERATELY, unlike ln_fluid_solve above. By the same
    // argument this should be the runner, but the reference for `moment3` on
    // this path is the JAR (MATLAB does not finish the method on the model the
    // golden was taken from), and the JAR's entry response times agree with the
    // first-order switch to 5e-3 and not with the resolved method, which moves
    // E:e0 from 2.2748 to 2.4964. Switching it blind would replace a measured
    // agreement with an unmeasured one; it needs a JAR-side reading first.
    const fluid::FluidSolution s = fluid::solver_fluid(L, o);
    std::vector<std::vector<fluid::FluidPassage>> out(
        L.nstations, std::vector<fluid::FluidPassage>(L.nclasses));
    for (std::size_t i = 1; i <= L.nstations; ++i) {
        if (L.stations[i - 1].nodetype == qn::NodeType::Source) continue;
        for (std::size_t r = 1; r <= L.nclasses; ++r) {
            if (L.disabled[i - 1][r - 1]) continue;
            out[i - 1][r - 1] = fluid::fluid_passage_time(L, s.xvec, i, r, o.tol, 201, s.closure);
        }
    }
    return out;
}
template <class T>
std::vector<std::vector<fluid::FluidPassage>> ln_fluid_cdf_respt(const qn::NetworkStruct<T>&,
                                                                 const fluid::FluidOptions&) {
    throw UnsupportedError(
        "SolverLN: the response-time CDF a layer contributes to the 'moment3' update is a fluid "
        "passage time, integrated by LSODA in double precision; rerun with --arith double");
}

/**
 * One layer's transient trajectory, the `getTranAvg` of a layer.
 *
 * `out_grid` REPLACES the uniform grid, and it is not a convenience: SolverENV
 * sums a layered stage's exit average against the environment's holding-time
 * CDF, so the trajectory has to be evaluated where that CDF puts its mass
 * rather than where the horizon happens to spread points. It also makes every
 * layer of one stage report on the SAME grid, which is what lets the layers be
 * assembled into one block-diagonal transient with a single time base.
 */
inline std::vector<fluid::FluidTranPoint> ln_fluid_transient(
    const qn::NetworkStruct<double>& L, const fluid::FluidOptions& fo, double t_end,
    std::size_t points, const std::vector<double>& out_grid = std::vector<double>()) {
    return fluid::solver_fluid_transient(L, fo, t_end, points, out_grid);
}
template <class T>
std::vector<fluid::FluidTranPoint> ln_fluid_transient(const qn::NetworkStruct<T>&,
                                                      const fluid::FluidOptions&, double,
                                                      std::size_t,
                                                      const std::vector<double>& =
                                                          std::vector<double>()) {
    throw UnsupportedError(
        "SolverLN: the layered transient integrates each layer's drift with LSODA, which is "
        "double precision by construction; rerun with --arith double");
}

}  // namespace detail

/** Options of SolverLN. Defaults are SolverOptions('LN'). */
struct LnOptions {
    int iter_max = 200;
    double iter_tol = 5e-3;
    double tol = 1e-4;
    bool interlocking = true;
    std::string relax = "fixed";
    double relax_factor = 0.5;
    /** Options handed to each layer solver; SolverMVA defaults. */
    mva::MvaOptions layer;
    /**
     * Which solver runs each layer: `mva`, `nc`, `fluid` or `ssa`.
     *
     * The reference names this by passing a solver FACTORY,
     * `LN(model, @(m) Fluid(m, opt))`, so the choice is per-ensemble and not
     * per-layer; this field is the same choice under a name, because a C++
     * template cannot take a MATLAB function handle.
     *
     * THEY DO NOT CONVERGE TO THE SAME PLACE. Each engine feeds different
     * service demands back into the next outer iteration, so the ensembles
     * reach different fixed points rather than the same one by different
     * routes. `ssa` in particular is a NOISY layer solver: the deterministic
     * convergence test cannot terminate against its standard error, so it is
     * driven with `LnStochController` (lqn_analyzers.h) instead.
     */
    std::string layer_solver = "mva";
    /** Options handed to each layer when `layer_solver` is `fluid`. */
    fluid::FluidOptions layer_fluid;
    /** Options handed to each layer when `layer_solver` is `nc`. */
    nc::NcSolverOptions layer_nc;
    /** Options handed to each layer when `layer_solver` is `ssa`. */
    ssa::SsaOptions layer_ssa;
    /**
     * `options.method`, which selects WHAT is reported and not merely how:
     *
     *   `default`     the mean-based update (updateMetricsDefault)
     *   `moment3`     the APH moment-based update (updateMetricsMomentBased),
     *                 which additionally produces a per-entry response-time
     *                 distribution, `get_cdf_respt()`
     *   `mwba.upper`  Majumdar-Woodside robust box bounds INSTEAD of a fixed
     *   `mwba.lower`  point: no layer is ever solved, and every metric the
     *                 bound does not define is reported as undefined
     */
    std::string method = "default";
    /**
     * `options.config.ln_transient`: how the per-layer transients are coupled.
     * `decoupled` freezes the inter-layer demands at the converged fixed point;
     * `coupled` reconciles them by waveform relaxation. Iteration 0 of the
     * coupled relaxation IS the decoupled answer.
     */
    std::string ln_transient = "coupled";
    /** `options.config.ln_transient_iter_max` and `..._tol` of the relaxation. */
    long ln_transient_iter_max = 20;
    double ln_transient_tol = 1e-2;
    /**
     * `options.config.ln_transient_channels`: which inter-layer coupling is
     * injected, `both`, `thinkt` (client delay only) or `callservt`
     * (synchronous-call service only). Used to isolate each channel's share of
     * the coupled transient.
     */
    std::string ln_transient_channels = "both";
    /** `options.timespan(2)`: the transient horizon; infinite means none is set. */
    double timespan_end = std::numeric_limits<double>::infinity();
    /** Output points per layer trajectory, and the relaxation's own grid size. */
    std::size_t tran_points = 101;
    /**
     * An EXPLICIT output grid for the layered transient, replacing the uniform
     * `tran_points` one when it is not empty.
     *
     * SolverENV is what asks for it: a layered stage's exit metric is a
     * Riemann-Stieltjes sum against the environment's holding-time CDF, and a
     * uniform grid over the horizon resolves the horizon rather than the
     * sojourn. Interpolating afterwards cannot recover resolution the
     * trajectory never had, so the points are asked for instead --
     * `SolverEnv::stage_grid` builds them and `set_tran_grid` installs them.
     * The grid must be increasing and start at zero; LSODA takes it as given.
     */
    std::vector<double> tran_grid;
};

/** Per-layer results of one iteration, the [QN,UN,RN,TN,AN,WN] of getAvg. */
template <class T>
struct LayerResult {
    Matrix<T> QN, UN, RN, TN, WN;
};

/** The LQN-level answer, indexed by element 1..nidx. */
template <class T>
struct LnSolution {
    std::vector<T> QN, UN, RN, TN, AN, WN;
    std::vector<bool> defined_Q, defined_U, defined_R, defined_T, defined_A, defined_W;
    int iterations = 0;
    bool converged = false;
    /**
     * True when the numbers are a BOUND (`method` = `mwba.upper` / `mwba.lower`)
     * rather than the fixed point. A bound defines throughput and processor
     * utilization and nothing else, so the other four measures come back
     * undefined; reporting zeros there would be a claim.
     */
    bool is_bound = false;
};

/** A CDF sampled on a grid, the [F, t] pair MATLAB's evalCDF returns. */
struct LnCdf {
    std::vector<double> t, cdf;
    bool empty() const { return t.empty(); }
};

/**
 * One layer's block of the layered transient.
 *
 * The reference assembles the layers into one BLOCK-DIAGONAL cell array whose
 * off-diagonal blocks are empty, so the blocks themselves carry the whole
 * answer; keeping them apart also keeps each layer's own time grid, which the
 * block-diagonal form has nowhere to put.
 */
struct LnTranLayer {
    std::vector<double> t;  ///< the output grid, shared by every series below
    /// [station][class][point]
    std::vector<std::vector<std::vector<double>>> QN, UN, TN;
};

/**
 * `LayeredNetwork.layerBlocks`: where each layer's block sits in the aggregate.
 *
 * An LQN has no stations and classes of its own -- SolverLN builds them, one
 * network per layer -- so the flat (station, class) view a caller like
 * SolverENV needs is the BLOCK-DIAGONAL UNION of the layer networks, in layer
 * order. `roff[e]`/`coff[e]` are layer e's 0-based row/column offset in that
 * union and `msz[e]`/`ksz[e]` its size; `M`/`K` are the totals.
 *
 * The off-diagonal blocks pair a station of one layer with a class of another
 * and stand for nothing at all. The reference leaves them as EMPTY cells; this
 * port fills them with zeros, which is what the reference's consumers make of
 * an empty cell (`tranTimeBase_` returns nothing and the metric stays 0).
 */
struct LnLayerBlocks {
    std::vector<std::size_t> roff, coff, msz, ksz;
    std::size_t M = 0, K = 0;
};

/** The layered transient: one block per layer, plus how it was produced. */
struct LnTranSolution {
    std::vector<LnTranLayer> layers;
    std::string mode;      ///< "coupled" or "decoupled"
    long iterations = 0;   ///< waveform-relaxation sweeps; 0 when decoupled
    double gap = 0.0;      ///< final sup-norm trajectory change, coupled only
};

/** getSensitivityTable of the ensemble: the layer tables under a Layer column. */
template <class T>
struct LnSensTable {
    struct Row {
        std::string layer, station, jobclass;
        T dTput, dRespT, dQLen, dUtil;
    };
    std::vector<Row> rows;
    /** Per layer, the branch that layer took; empty for a layer with no solver. */
    std::vector<std::string> layer_methods;
    /** The summary label: the common branch, or "mixed" when they differ. */
    std::string method;
    /** Per layer, the analytic Jacobian where that layer produced one. */
    std::vector<sens::SensTable<T>> layer_tables;
};

/**
 * Overtaking probability at a server entry, defined in lqn_analyzers.h.
 *
 * DECLARED, not included: lqn_analyzers.h needs LayerResult and SolverLN
 * complete (LnStochController holds a vector of the first and two refusals take
 * the second), so it must be parsed AFTER this class. The definition arrives
 * through the include at the foot of this file, which is why the declaration
 * has to stand here -- an unqualified call from a member function would
 * otherwise find nothing, LqnStruct's associated namespace being line::lqn.
 */
template <class T>
T lqn_overtake_prob_markov(const LqnStruct<T>& lqn, const std::vector<T>& servt,
                           const std::vector<T>& callresidt, const std::vector<T>& tput,
                           std::size_t eidx, const T& xj);

/** `options.config.stochiter_*` of SolverOptions.m, with its defaults. */
struct LnStochConfig {
    long burnin = 5;      ///< Picard iterations before the step decay starts
    double a0 = 1.0;      ///< Robbins-Monro step immediately after burn-in
    double alpha = 0.6;   ///< step decay exponent, in (0.5, 1]
    long conseq = 3;      ///< consecutive sub-tolerance iterations required to stop
    double iter_tol = 5e-3;
    /** Relaxation in force during burn-in, i.e. whatever init left in place. */
    double relax_burnin = 1.0;
};

/**
 * The Robbins-Monro / Polyak-Ruppert controller, defined in lqn_analyzers.h.
 *
 * DECLARED HERE, not included: same mutual dependence as above. `iterate` holds
 * one through a shared_ptr rather than by value so that this class stays
 * complete without it -- a by-value member would need the definition at the
 * point the class template is instantiated, and the definition arrives at the
 * foot of this file.
 *
 * ITS CONFIG CANNOT BE DEFERRED THE SAME WAY. `iterate` names LnStochConfig by
 * value, and that name does not depend on T, so it is looked up and required
 * COMPLETE when the template is parsed rather than when it is instantiated --
 * the shared_ptr trick only defers the class template beside it.
 */
template <class T>
class LnStochController;

template <class T>
class SolverLN {
public:
    /**
     * Port of `SolverLN.listValidMethods`.
     *
     * Each name states the LAYERING and the ENCODING; `ln_requested_method`
     * normalises the alias spellings ("ph", "cs", "srvncs", "flatcs",
     * "squashed", "squashed.ph") onto these, and they are left out here to keep
     * the list unambiguous, exactly as the reference does.
     */
    static std::vector<std::string> list_valid_methods() {
        return {"srvn", "srvn.ph", "srvn.cs", "flat", "flat.cs", "flat.ph", "moment3", "default"};
    }

    SolverLN(const LqnStruct<T>& lqn_in, const LnOptions& options) : lqn(lqn_in), opt(options) {
        // A LAYER SOLVER NEVER NARRATES, the same invariant the MATLAB, JAR and
        // python ports stamp on each constructed layer solver. The fixed point
        // runs every layer once per iteration, so a layer left at the caller's
        // verbosity would print its own banner nlayers*iter_max times and bury
        // the layered narration the caller actually asked for. `ssa` is the one
        // layer engine here carrying a verbosity knob of its own; the level is
        // forced rather than trusted so a library caller cannot set it.
        // SolverLN's own reporting is unaffected.
        opt.layer_ssa.verbose = false;
        construct();
    }

    /** Port of getEnsembleAvg: run the iteration and aggregate onto LQN elements. */
    LnSolution<T> get_ensemble_avg() {
        // getAvgTable.m:611-624 answers the two bound requests WITHOUT running
        // the fixed point at all: a box bound is a statement about the model,
        // not about an iterate, so solving first and discarding the answer would
        // only cost time and invite the two to be confused.
        if (opt.method == "mwba.upper" || opt.method == "mwba.lower") return box_bounds();
        iterate();
        // the layers of "srvn.ph" carry one class per caller task, so the
        // per-element results are rebuilt analytically -- see aggregate_ph
        return is_ph_encoding() ? aggregate_ph() : aggregate();
    }

    /**
     * Port of @@SolverLN/getCdfRespT: the per-entry response-time distribution.
     *
     * ONLY THE `moment3` METHOD PRODUCES ONE. The reference reacts to an empty
     * table by re-running getAvg under method='moment3' and restoring the
     * caller's method afterwards, and so does this: the distribution is the
     * whole point of that method, and the mean-based update has no distribution
     * to report, not a coarser one.
     *
     * Indexed by entry NUMBER 1..nentries, as the reference's
     * `entrycdfrespt{eidx - (nhosts+ntasks)}` is.
     */
    std::vector<LnCdf> get_cdf_respt() {
        if (lqn.nentries == 0) return entrycdfrespt;
        if (entrycdfrespt.size() <= lqn.nentries || entrycdfrespt[1].empty()) {
            // The distribution pass reads the ROUTING encoding of the activity
            // graph, which srvn.ph layers do not carry: re-running the iteration
            // over them would reconstruct the wrong topology rather than a
            // coarser answer, and returning the empty table would report no
            // distribution at all. Refuse by name, as the reference does
            // (@SolverLN/getCdfRespT).
            if (is_ph_encoding())
                throw UnsupportedError(
                    "getCdfRespT needs the routing encoding of the activity graph, which "
                    "method='srvn.ph' does not build. Rebuild the solver with "
                    "method='srvn.cs' or method='moment3'.");
            // BOTH the option and the RESOLVED method have to move. update_metrics
            // dispatches on `lnmethod`, which build_layers resolved once, so
            // flipping `opt.method` alone leaves the mean-based update in place and
            // the table empty -- which is what this getter did between the alias
            // landing and 2026-08-11. The routing layers already built serve
            // moment3 unchanged, so only the update pass changes.
            const std::string saved = opt.method;
            const std::string saved_resolved = lnmethod;
            opt.method = "moment3";
            lnmethod = "moment3";
            iterate();
            aggregate();
            opt.method = saved;
            lnmethod = saved_resolved;
        }
        return entrycdfrespt;
    }

    /**
     * Port of @@SolverLN/getTranAvg: the block-diagonal aggregate transient.
     *
     * `opt.ln_transient` selects the coupling; both modes return the same
     * layout, and iteration 0 of the coupled relaxation IS the decoupled
     * answer.
     */
    LnTranSolution get_tran_avg() {
        if (opt.ln_transient == "decoupled") return tran_avg_decoupled();
        if (opt.ln_transient == "coupled") return tran_avg_coupled();
        throw InputError("SolverLN: unknown ln_transient mode '" + opt.ln_transient +
                         "' (use 'coupled' or 'decoupled')");
    }

    /**
     * Port of @@SolverLN/getSensitivityTable: solve the ensemble, then
     * concatenate each layer solver's own table under a leading Layer column.
     *
     * WHAT THESE DERIVATIVES MEAN, and the reference is emphatic about it: each
     * entry is a partial derivative WITHIN ITS LAYER, taken with the layer
     * parameters the fixed point produced held fixed. Perturbing a host demand
     * moves the think times, populations and demands of every other layer
     * through the fixed-point map, and that indirect term is NOT included. The
     * table attributes a bottleneck inside a layer; it does not predict the
     * effect of a parameter change on the solved layered model.
     */
    LnSensTable<T> get_sensitivity_table(const sens::SensOptions& sopt) {
        if (results.empty()) iterate();
        LnSensTable<T> out;
        out.layer_methods.assign(ensemble.size(), std::string());
        out.layer_tables.resize(ensemble.size());
        for (std::size_t e = 0; e < ensemble.size(); ++e) {
            sens::SensOptions so = sopt;
            so.simulation = opt.layer_solver == "ssa";
            const bool exact_available =
                opt.layer_solver == "mva" || opt.layer_solver == "nc";
            sens::SensTable<T> t = sens::solver_sensitivity_table<T>(
                ensemble[e], so, exact_available && !fj_tr[e].active(),
                [this, e]() { return this->solve_layer(e); });
            out.layer_methods[e] = t.method;
            for (const sens::SensRow<T>& r : t.rows) {
                typename LnSensTable<T>::Row row;
                row.layer = ensemble[e].name;
                row.station = r.station;
                row.jobclass = r.jobclass;
                row.dTput = r.dTput;
                row.dRespT = r.dRespT;
                row.dQLen = r.dQLen;
                row.dUtil = r.dUtil;
                out.rows.push_back(row);
            }
            out.layer_tables[e] = t;
        }
        // One branch label per layer, plus a summary that is `mixed` when the
        // layers did not all take the same branch.
        for (const std::string& m : out.layer_methods) {
            if (m.empty()) continue;
            if (out.method.empty()) out.method = m;
            else if (out.method != m) out.method = "mixed";
        }
        return out;
    }

    /**
     * Majumdar-Woodside robust box bounds, reported in the shape of a solution.
     *
     * `mwba.upper` reports the upper bound of both throughput and processor
     * utilization, `mwba.lower` the lower bound of both. Everything the bound
     * does not define -- queue length, response time, residence time, arrival
     * rate -- is left UNDEFINED rather than zeroed.
     */
    LnSolution<T> box_bounds() {
        const bool upper = opt.method != "mwba.lower";
        // Fully qualified: the member `lqn` shadows the namespace of the same
        // name inside this class, so an unqualified `lqn::` would not compile.
        const ::line::lqn::LqnBoxBounds<T> b = ::line::lqn::lqn_boxbounds(lqn);
        LnSolution<T> s;
        const std::size_t N = lqn.nidx;
        auto blank = [&](std::vector<T>& v, std::vector<bool>& d) {
            v.assign(N + 1, Tzero());
            d.assign(N + 1, false);
        };
        blank(s.QN, s.defined_Q);
        blank(s.UN, s.defined_U);
        blank(s.RN, s.defined_R);
        blank(s.TN, s.defined_T);
        blank(s.AN, s.defined_A);
        blank(s.WN, s.defined_W);
        for (std::size_t i = 1; i <= N; ++i) {
            if (b.defined_T[i]) {
                s.TN[i] = upper ? b.TN_up[i] : b.TN_lo[i];
                s.defined_T[i] = true;
            }
            if (b.defined_U[i]) {
                s.UN[i] = upper ? b.UN_up[i] : b.UN_lo[i];
                s.defined_U[i] = true;
            }
        }
        s.iterations = 0;
        s.converged = true;  // a bound needs no fixed point to have converged to
        s.is_bound = true;
        return s;
    }

    std::size_t nlayers() const { return ensemble.size(); }
    const std::vector<qn::Layer<T>>& layers() const { return ensemble; }

    /**
     * Port of `LayeredNetwork.layerBlocks`: the block-diagonal layout of the
     * layers in the aggregate (station x class) view.
     *
     * Available as soon as the solver is constructed -- `build_layers` runs in
     * the constructor -- which is what lets SolverENV compare the shapes of its
     * stages before solving any of them.
     */
    LnLayerBlocks layer_blocks() const {
        LnLayerBlocks b;
        const std::size_t E = ensemble.size();
        b.roff.assign(E, 0);
        b.coff.assign(E, 0);
        b.msz.assign(E, 0);
        b.ksz.assign(E, 0);
        for (std::size_t e = 0; e < E; ++e) {
            b.roff[e] = b.M;
            b.coff[e] = b.K;
            b.msz[e] = ensemble[e].nstations;
            b.ksz[e] = ensemble[e].nclasses;
            b.M += b.msz[e];
            b.K += b.ksz[e];
        }
        return b;
    }

    /**
     * Port of `LayeredNetwork.initFromMarginal`: split an aggregate (M x K)
     * mean queue-length matrix into per-layer blocks and warm-start each layer
     * from its own.
     *
     * WHY THIS IS RECORDED AND NOT APPLIED HERE, which is the trap the JAR and
     * python ports both hit: the layered fixed point RESETS its layers as it
     * converges, so a warm start installed before the solve does not survive
     * it. The blocks are therefore kept and replayed by `run_layer_transient`,
     * i.e. after the fixed point and immediately before each layer's transient
     * -- the steady solve ignores an initial state, the transient does not.
     *
     * THE BLOCK IS NOT ROUNDED. The reference rounds a replayed block onto the
     * integer lattice for every layer engine EXCEPT the fluid one, and the
     * layered transient exists only over fluid layers (`require_transient_ready`
     * refuses the rest), so the branch that would round has no reachable case
     * here. Rounding a fluid state would quantize the very quantity SolverENV's
     * fixed point is iterating on.
     *
     * An empty matrix clears the warm start, so a caller can put the layers
     * back on their default initial state without rebuilding the solver.
     */
    void init_from_marginal(const Matrix<double>& n) {
        const std::size_t E = ensemble.size();
        layer_tran_init.assign(E, std::vector<double>());
        if (n.empty()) return;
        const LnLayerBlocks b = layer_blocks();
        if (n.rows() != b.M || n.cols() != b.K)
            throw InputError(
                "SolverLN::init_from_marginal: the marginal is " + std::to_string(n.rows()) + "x" +
                std::to_string(n.cols()) + " where the block-diagonal union of the layers is " +
                std::to_string(b.M) + "x" + std::to_string(b.K) +
                "; the aggregate view is layerBlocks', not the LQN element count");
        for (std::size_t e = 0; e < E; ++e) {
            const qn::Layer<T>& L = ensemble[e];
            const fluid::FluidLayout lay = fluid::fluid_layout(L);
            std::vector<double> y(lay.nstates, 0.0);
            bool any = false;
            for (std::size_t i = 0; i < L.nstations && i < b.msz[e]; ++i)
                for (std::size_t r = 0; r < L.nclasses && r < b.ksz[e]; ++r) {
                    if (!lay.enabled[i][r]) continue;
                    const double q = std::max(0.0, n(b.roff[e] + i, b.coff[e] + r));
                    y[lay.qidx[i][r]] = q;
                    if (q > 0.0) any = true;
                }
            // An all-zero block is NOT a warm start: it is the absence of one,
            // and installing it would empty a closed layer whose population the
            // fixed point conserves. The reference reaches the same place by
            // never calling initFromMarginal before the first stage solve.
            if (any) layer_tran_init[e] = y;
        }
    }

    /** Install the explicit output grid of the layered transient; see `LnOptions::tran_grid`. */
    void set_tran_grid(const std::vector<double>& g) { opt.tran_grid = g; }

    /** Diagnostic access to the per-iteration layer results, for the regression. */
    const std::vector<std::vector<LayerResult<T>>>& iteration_results() const { return results; }
    const std::vector<T>& state_servt() const { return servt; }
    const std::vector<T>& state_residt() const { return residt; }
    const std::vector<T>& state_tput() const { return tput; }
    const std::vector<T>& state_thinkt() const { return thinkt; }
    const std::vector<T>& state_callservt() const { return callservt; }
    const std::vector<T>& state_callresidt() const { return callresidt; }
    const std::vector<Distrib<T>>& state_servtproc() const { return servtproc; }
    const std::vector<Distrib<T>>& state_thinktproc() const { return thinktproc; }
    const std::vector<Distrib<T>>& state_callservtproc() const { return callservtproc; }
    const std::vector<T>& state_util() const { return util; }
    /** The method the layers were built for: "srvn.ph", "srvn.cs" or "moment3". */
    const std::string& state_lnmethod() const { return lnmethod; }

private:
    // state
    LqnStruct<T> lqn;
    LnOptions opt;

    std::vector<qn::Layer<T>> ensemble;   ///< compacted, layer e is ensemble[e]
    std::vector<long> idxhash;            ///< (nidx+1) element -> layer index, -1 if none
    std::vector<bool> ignore;             ///< (nidx+1) element in a component with no ref task

    // update maps, rows of [layerElement, elementOrCall, node, class]
    struct UpdRow {
        std::size_t idx, aidx, node, cls;
    };
    std::vector<UpdRow> servt_map, thinkt_map, call_map, actthinkt_map;
    /**
     * Source arrival of an async call class, `aidx` holding the call index.
     *
     * Kept apart from call_map because the two update different things about
     * the same class: call_map carries the SERVICE at the server replicas,
     * this one the arrival RATE at the Source, which is the caller activity's
     * throughput and so moves with the fixed point. Entry open arrivals have
     * no row here at all -- theirs is exogenous and never reseeded.
     */
    std::vector<UpdRow> arv_call_map;
    /** [idx, tidx_caller, eidx, nodefrom, nodeto, classfrom, classto] */
    struct RouteRow {
        std::size_t idx, tidx_caller, eidx, nodefrom, nodeto, cfrom, cto;
    };
    std::vector<RouteRow> route_map;
    std::vector<std::size_t> unique_route_idx;
    std::vector<std::size_t> route_reset, svc_reset;

    Matrix<double> njobs;  ///< (NT+1 x NT+1) population of caller tidx in layer idx

    // per-element metric state
    std::vector<T> servt, residt, tput, util, thinkt;
    std::vector<T> callservt, callresidt;
    /**
     * Phase-2 state, the reference's hasPhase2 / servt_ph1 / servt_ph2 /
     * prOvertake. The split is of `servt` and exists only while `has_phase2`;
     * `prOvertake` is indexed by entry NUMBER 1..nentries, as the reference
     * indexes it, not by element index.
     */
    bool has_phase2 = false;
    std::vector<T> servt_ph1, servt_ph2, prOvertake;
    /** Tasks whose own layer was collapsed to one representative replica. */
    std::set<std::size_t> single_replica_tasks;
    std::vector<Distrib<T>> servtproc, thinkproc, thinktproc, tputproc, callservtproc;
    /**
     * `moment3` state: the response-time CDF of each activity (by element) and
     * of each call, the fitted per-entry law and the CDF read off it. All empty
     * under the default method, which never forms a distribution at all.
     */
    std::vector<fluid::FluidPassage> servtcdf, callservtcdf;
    std::vector<LnCdf> entrycdfrespt;
    std::vector<mam::AphPair<T>> entryproc;
    /** Per layer, the response-time CDFs of the converged pass; a lazy cache. */
    std::vector<std::vector<std::vector<fluid::FluidPassage>>> cdf_repo;
    std::vector<double> servt_prev, residt_prev, tput_prev, thinkt_prev;
    std::vector<double> callservt_prev, callresidt_prev;
    std::vector<T> servt_prev_v, residt_prev_v, tput_prev_v, thinkt_prev_v, callservt_prev_v;
    Matrix<T> servtmatrix;  ///< (nidx+ncalls) x (nidx+ncalls) entry reachability

    // interlock tables
    Matrix<T> il_all, il_ph1;
    std::vector<std::vector<std::size_t>> il_common_entries, il_src_all, il_src_ph2;
    /**
     * Per-layer interlock matrix of Franks (1999), Eq. (4.7), CLASS-indexed. A host layer
     * whose engine carries the correction inside its own MVA takes the matrix instead of
     * having its residence times scaled after the fact; empty for every other layer.
     */
    std::vector<std::vector<std::vector<double>>> layer_interlock;
    std::vector<double> il_num_sources;

    std::vector<std::vector<LayerResult<T>>> results;  ///< [iteration][layer]
    std::vector<Matrix<T>> layer_init_sol;
    /**
     * Per layer, the ODE state `init_from_marginal` recorded, replayed by
     * `run_layer_transient`. Empty for a layer with no warm start.
     */
    std::vector<std::vector<double>> layer_tran_init;
    /**
     * Per layer, the fork-join transform SolverMVA solves instead of the layer
     * itself (inactive for a layer with no fork), and the auxiliary arrival
     * rates its fixed point iterates on. See fj_mmt.h.
     *
     * `fj_lambda` is the `self.fjForkLambda` of the reference and is deliberately
     * NOT reset between outer iterations: the transform is rebuilt cold on every
     * outer pass but the iterate it converged to last time is a far better
     * starting point than FineTol, and the reference warm-starts it the same way
     * (options.config.fj_warmstart, true by default).
     */
    std::vector<mva::FjMmt<T>> fj_tr;
    std::vector<std::vector<T>> fj_lambda;
    double relax_omega = 1.0;
    /** Live only while a NOISY layer engine is in force; null otherwise. */
    std::shared_ptr<LnStochController<T>> stoch_ctl;
    long averagingstart = -1;
    bool hasconverged = false;
    /** method='moment3': true once the distribution pass has formed the entry laws. */
    bool moment_pass_done = false;
    std::vector<double> maxitererr;
    int iterations_done = 0;
    bool did_converge = false;

    std::size_t NT() const { return lqn.tshift + lqn.ntasks; }

    static T Tzero() { return num_traits<T>::from_int(0); }
    static T Tone() { return num_traits<T>::from_int(1); }
    static double dbl(const T& x) { return num_traits<T>::to_double(x); }

    // -----------------------------------------------------------------------
    // construct
    // -----------------------------------------------------------------------
    void construct() {
        // Forwarding rewrite, SolverLN.m:162-164: flatten every forwarding
        // chain into caller-side pseudo rendezvous calls BEFORE anything else
        // reads lqn, so build_layers/detect_phase2/reject_unsupported all see
        // plain rendezvous arcs. The raw FWD calls survive in the struct but
        // must not contribute blocking after this point (lqn_helpers.h).
        lqn_fwd_rendezvous(lqn);
        reject_unsupported();
        detect_phase2();

        const std::size_t N = lqn.nidx;
        ignore.assign(N + 1, false);
        // weakly connected components of graph + graph'; a component with no
        // reference task is unreachable and its elements are ignored
        std::vector<long> comp(N + 1, -1);
        long ncomp = 0;
        for (std::size_t v = 1; v <= N; ++v) {
            if (comp[v] >= 0) continue;
            std::vector<std::size_t> stack{v};
            comp[v] = ncomp;
            while (!stack.empty()) {
                const std::size_t u = stack.back();
                stack.pop_back();
                for (std::size_t w : lqn.graph.succ(u))
                    if (comp[w] < 0) {
                        comp[w] = ncomp;
                        stack.push_back(w);
                    }
                for (std::size_t w : lqn.graph.pred(u))
                    if (comp[w] < 0) {
                        comp[w] = ncomp;
                        stack.push_back(w);
                    }
            }
            ++ncomp;
        }
        if (ncomp > 1) {
            std::vector<bool> has_ref(ncomp, false);
            for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
                const std::size_t tidx = lqn.tshift + t;
                if (lqn.sched[tidx] == SchedStrategy::REF) has_ref[comp[tidx]] = true;
            }
            for (std::size_t v = 1; v <= N; ++v)
                if (!has_ref[comp[v]]) ignore[v] = true;
        }

        servtproc.assign(N + 1, Distrib<T>::disabled_dist());
        thinkproc.assign(N + 1, Distrib<T>::disabled_dist());
        thinktproc.assign(N + 1, Distrib<T>::disabled_dist());
        tputproc.assign(N + 1, Distrib<T>::disabled_dist());
        for (std::size_t i = 1; i <= N; ++i) {
            servtproc[i] = lqn.hostdem[i];
            thinkproc[i] = lqn.think[i];
        }
        callservtproc.assign(lqn.ncalls + 1, Distrib<T>::disabled_dist());
        for (std::size_t c = 1; c <= lqn.ncalls; ++c)
            callservtproc[c] = lqn.hostdem[lqn.callpair_dst[c]];

        njobs = Matrix<double>(NT() + 1, NT() + 1, 0.0);
        build_layers();

        // layers whose routing or service must be re-derived after an update
        std::vector<bool> rr(NT() + 1, false), sr(NT() + 1, false);
        for (const RouteRow& r : route_map) rr[r.idx] = true;
        for (const UpdRow& r : thinkt_map) sr[r.idx] = true;
        for (const UpdRow& r : call_map) sr[r.idx] = true;
        // a moved arrival rate reweights the Sink closure, which only
        // refresh_chains rebuilds, so these layers need the routing reset
        for (const UpdRow& r : arv_call_map) rr[r.idx] = true;
        for (std::size_t i = 1; i <= NT(); ++i) {
            if (rr[i] && idxhash[i] >= 0) route_reset.push_back(std::size_t(idxhash[i]));
            if (sr[i] && idxhash[i] >= 0) svc_reset.push_back(std::size_t(idxhash[i]));
        }
        std::sort(route_reset.begin(), route_reset.end());
        route_reset.erase(std::unique(route_reset.begin(), route_reset.end()), route_reset.end());
        std::sort(svc_reset.begin(), svc_reset.end());
        svc_reset.erase(std::unique(svc_reset.begin(), svc_reset.end()), svc_reset.end());
    }

    /** Reject, by name, every construct this port does not implement. */
    // Phase-2 detection, SolverLN.m:165-173. Deliberately NOT inside
    // reject_unsupported: that method is const because it only refuses, and
    // detection sets state. Marking the flag mutable would have compiled and
    // left a validator that silently mutates the solver.
    void detect_phase2() {
        has_phase2 = false;
        for (std::size_t a = 1; a <= lqn.nacts; ++a)
            if (lqn.actphase[a] > 1) has_phase2 = true;
    }

    /**
     * True when `eidx` is the destination of a raw forwarding call.
     *
     * buildLayersRecursive.m:186-192 keeps this test even though
     * lqn_fwd_rendezvous has already run, and so must this port:
     * the rewrite walks chains out of SYNC calls only, so an entry whose
     * forwarder is reached asynchronously gets no pseudo arc and would
     * otherwise be dropped from the layer along with all of its work.
     */
    bool is_fwd_target(std::size_t eidx) const {
        for (std::size_t c = 1; c <= lqn.ncalls; ++c)
            if (lqn.calltype[c] == CallType::FWD && lqn.callpair_dst[c] == eidx) return true;
        return false;
    }

    /**
     * Total exogenous rate into the entries of task `tidx`, zero unless the arrival is
     * the ONLY way in.
     *
     * A task nobody calls has no task layer, so update_think_times never sets its
     * surrogate delay and its caller class cycles against an Immediate one. Carrying the
     * arrival as an open stream ON TOP of that unthrottled chain loads the host twice:
     * lqn_open_arrival read the processor at 0.68 where lqns, lqsim and LDES all give
     * 0.32. The chain is the representation that honours the thread pool, so build_layer
     * drops the stream for these tasks and update_think_times closes the chain on this
     * rate, as the reference does for a forwarding target. With a caller or a forwarding
     * source the stream needs a class of its own and this returns 0.
     */
    double open_arrival_rate_of(std::size_t tidx) const {
        if (lqn.isref[tidx]) return 0.0;
        for (std::size_t e : lqn.entriesof[tidx])
            if (lqn.issynccaller.any_col(e) || lqn.isasynccaller.any_col(e) || is_fwd_target(e))
                return 0.0;
        double rate = 0.0;
        for (std::size_t e : lqn.entriesof[tidx]) {
            if (!lqn.has_arrival[e]) continue;
            const double m = dbl(lqn.arrival[e].mean);
            if (std::isfinite(m) && m > GlobalConstants::FineTol) rate += 1.0 / m;
        }
        return rate;
    }

    /** True when layer `e` carries a Cache node. */
    bool has_cache_node(std::size_t e) const {
        for (std::size_t n = 0; n < ensemble[e].nodes.size(); ++n)
            if (ensemble[e].nodes[n].nodetype == NodeType::Cache) return true;
        return false;
    }

    /** True when `aidx` is an activity bound to an entry nobody calls synchronously. */
    bool async_only_activity(std::size_t aidx) const {
        if (aidx <= lqn.ashift || aidx > lqn.ashift + lqn.nacts) return false;
        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            if (lqn.graph.get(eidx, aidx) == Tzero()) continue;
            return lqn.isasynccaller.any_col(eidx) && !lqn.issynccaller.any_col(eidx);
        }
        return false;
    }

    /**
     * Reject, by name, every construct this port does not implement.
     *
     * NOTHING IS LEFT. Forwarding, asynchronous calls, entry open arrivals,
     * admission constraints, cache tasks and setup tasks are each solved.
     * The method is kept rather than deleted because it is the place a new
     * refusal belongs, and because a layer-build refusal must be raised HERE,
     * before construct() has read anything out of the struct, rather than
     * halfway through building a layer.
     */
    void reject_unsupported() const {}

    // -----------------------------------------------------------------------
    // buildLayers
    // -----------------------------------------------------------------------
    void build_layers() {
        // Method resolution. A method name carries both the LAYERING and the
        // ENCODING: "srvn.ph" replaces the routing encoding of the activity graph
        // by a composed phase-type server law, "srvn" is the alias that takes it
        // where it can serve the model and "srvn.cs" otherwise. The choice is
        // made ONCE, here, and every later dispatch reads lnmethod.
        // See _kb/06-solver-catalog.md (LN section).
        const std::string requested = ln_requested_method(opt.method);
        assert_call_groups(requested == "flat.cs");
        if (requested == "flat.cs") {
            lnmethod = "flat.cs";
            build_flat_layer();
            return;
        }
        if (requested == "flat.ph") {
            // the squashed layering with the composed law: ONE submodel holding
            // every server, and a caller visiting each of them once per
            // invocation. The feature gate is the srvn.ph one plus the refusals a
            // single submodel carries -- see ph_flat_server_set.
            ph_laws_ready = false;
            lnmethod = "flat.ph";
            build_layers_ph(true);
            return;
        }
        if (requested == "srvn.ph" || requested == "srvn") {
            ph_laws_ready = false;
            if (requested == "srvn.ph" || probe_srvn_ph()) {
                lnmethod = "srvn.ph";
                build_layers_ph();
                return;
            }
        }
        lnmethod = (requested == "moment3") ? "moment3" : "srvn.cs";
        std::vector<qn::Layer<T>> raw(NT() + 1);
        std::vector<bool> present(NT() + 1, false);

        for (std::size_t hidx = 1; hidx <= lqn.nhosts; ++hidx) {
            if (ignore[hidx]) continue;
            build_layer(raw[hidx], {hidx}, lqn.tasksof[hidx], true, false);
            present[hidx] = true;
        }
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (ignore[tidx] || lqn.isref[tidx]) continue;
            bool any_caller = lqn.iscaller.any_row(tidx) || lqn.iscaller.any_col(tidx);
            if (!any_caller) continue;
            // tasks that call some entry of tidx
            std::vector<std::size_t> callers;
            for (std::size_t ct = 1; ct <= lqn.ntasks; ++ct) {
                const std::size_t c = lqn.tshift + ct;
                bool calls = false;
                for (std::size_t e : lqn.entriesof[tidx])
                    if (lqn.iscaller.get(c, e)) calls = true;
                if (calls) callers.push_back(c);
            }
            if (callers.empty()) continue;
            build_layer(raw[tidx], {tidx}, callers, false, false);
            present[tidx] = true;
        }

        idxhash.assign(lqn.nidx + 1, -1);
        long next = 0;
        for (std::size_t i = 1; i <= NT(); ++i)
            if (present[i]) {
                idxhash[i] = next++;
                ensemble.push_back(std::move(raw[i]));
            }
        layer_init_sol.assign(ensemble.size(), Matrix<T>());
        build_fork_views();
    }

    /**
     * The processors and called tasks that become stations of the flat layer.
     *
     * Squashing is refused rather than approximated where an element carries
     * state that only a submodel of its own can hold: a REPLICATED element
     * would need one station per copy inside a layer whose routing addresses it
     * once, a CACHE task needs the Cache node in the host layer its reads queue
     * at, and a SETUP task's delay-off belongs to the station that powers down.
     */
    /**
     * Reject a routed call group under any layering or layer solver that cannot
     * carry it.
     *
     * Two conditions, and both are refusals rather than degradations. The
     * squashed layering is needed because under `srvn` each target lives in a
     * submodel of its own and is replaced, in the caller's submodel, by a
     * surrogate delay -- no node ever has arcs to more than one of them, so
     * there is nothing to dispatch among. A layer solver that resolves the
     * strategy from the STATE is needed because `refresh_routing` expands
     * RROBIN and JSQ into a uniform probability split for the matrix solvers,
     * and returning that split under a round-robin label misreports a
     * deterministic policy as a coin. In this port only `ssa` resolves them.
     */
    void assert_call_groups(bool flat) const {
        if (lqn.callgroups.empty()) return;
        if (!flat)
            throw UnsupportedError(
                "Call groups routed by a routing strategy require the squashed layering; use "
                "method='flat'. Under srvn the targets never share a submodel, so the dispatch "
                "order cannot be represented.");
        if (opt.layer_solver != "ssa")
            throw UnsupportedError(
                "Routed call groups need a layer solver that resolves the strategy from the "
                "state; set layer_solver='ssa'. MVA, NC and FLD read the routing matrix, into "
                "which refresh_routing has expanded the strategy as a uniform split, and would "
                "return that split under a round-robin or JSQ label.");
    }

    std::vector<std::size_t> flat_server_set() const {
        std::vector<std::size_t> servers;
        for (std::size_t i = 1; i <= NT(); ++i) {
            if (lqn.repl[i] > 1.0)
                throw UnsupportedError(
                    "Flat layering does not support replicated processors or tasks, use the "
                    "default 'srvn' layering.");
            if (lqn.iscache[i])
                throw UnsupportedError(
                    "Flat layering does not support cache tasks, use the default 'srvn' "
                    "layering.");
            if (lqn.hassetup[i])
                throw UnsupportedError(
                    "Flat layering does not support setup tasks, use the default 'srvn' "
                    "layering.");
        }
        for (std::size_t hidx = 1; hidx <= lqn.nhosts; ++hidx)
            if (!ignore[hidx] && !lqn.tasksof[hidx].empty()) servers.push_back(hidx);
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (ignore[tidx] || lqn.isref[tidx]) continue;
            if (!lqn.iscaller.any_row(tidx) && !lqn.iscaller.any_col(tidx)) continue;
            bool has_task_caller = false;
            for (std::size_t eidx : lqn.entriesof[tidx])
                for (std::size_t c = 1; c <= lqn.ntasks; ++c)
                    if (lqn.iscaller.get(lqn.tshift + c, eidx)) has_task_caller = true;
            if (has_task_caller) servers.push_back(tidx);
        }
        if (servers.empty())
            throw InputError(
                "Flat layering found no server: the model has no processor with tasks.");
        return servers;
    }

    /**
     * Build the single squashed layer holding every processor and called task.
     *
     * Every served element resolves to layer 0, which is at once the host layer
     * and the task layer, so a consumer that walks `idxhash` finds the same
     * network for a processor and for a task and must read the station off
     * `server_idx_of` rather than off `serverIdx`.
     */
    void build_flat_layer() {
        const std::vector<std::size_t> flat_servers = flat_server_set();
        std::vector<std::size_t> flat_callers;
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (!ignore[tidx]) flat_callers.push_back(tidx);
        }
        qn::Layer<T> layer;
        build_layer(layer, flat_servers, flat_callers, false, true);
        ensemble.clear();
        ensemble.push_back(std::move(layer));
        idxhash.assign(lqn.nidx + 1, -1);
        for (std::size_t s : flat_servers) idxhash[s] = 0;
        layer_init_sol.assign(ensemble.size(), Matrix<T>());
        build_fork_views();
    }

    /**
     * Build, for every layer that has a fork, the transformed model SolverMVA
     * actually solves, and seed its auxiliary arrival rates.
     */
    void build_fork_views() {
        fj_tr.assign(ensemble.size(), mva::FjMmt<T>());
        fj_lambda.assign(ensemble.size(), {});
        for (std::size_t e = 0; e < ensemble.size(); ++e) {
            if (!ensemble[e].has_fork()) continue;
            // fj_mmt mints its own Source/Sink pair and would silently overwrite
            // sourceIdx/sinkNode, detaching the open stream already routed there.
            // NOT a parity gap: the reference does not solve this combination
            // either, it fails inside the layer solve with "Arrays have
            // incompatible sizes" (checked 2026-07-29 on a fork+async model).
            // Refusing by name is the better answer, so this stays.
            if (ensemble[e].sourceIdx != 0)
                throw UnsupportedError(
                    "SolverLN: layer '" + ensemble[e].name +
                    "' carries both an AND fork and an open stream (an async call or an entry "
                    "arrival); the fork-join transform needs a Source of its own");
            fj_tr[e] = mva::fj_mmt(ensemble[e]);
            fj_lambda[e].assign(fj_tr[e].V.classes.size() + 1,
                                num_traits<T>::from_double(GlobalConstants::FineTol));
        }
    }

    /**
     * Port of matlab/src/lang/layered/lqn_dep_layer_handle.m: lift a
     * per-operand rate handle onto the classes of a layer station.
     *
     * COLS[j] lists the layer classes (0-based) through which operand j occupies
     * the station; CHAINCOLS is the same list in the CHAIN index space, filled
     * once refresh_chains has run. Solvers evaluate the handle in either space
     * -- CTMC and the exact recursions pass a per-class vector, the AMVA and NC
     * chain recursions a per-chain one -- so the handle reads the LENGTH of what
     * it is given to decide, aggregates the operand populations in that space,
     * and answers a vector of the same length, since the caller indexes the
     * answer with the index it passed in. An index belonging to no operand keeps
     * the neutral scaling 1.
     */
    static lang::CdScaling<T> layer_dep_handle(
        const lang::CdScaling<T>& f, const std::vector<std::vector<std::size_t>>& cols,
        const std::shared_ptr<std::vector<std::vector<std::size_t>>>& chaincols, std::size_t R) {
        return [f, cols, chaincols, R](const std::vector<T>& n) -> std::vector<T> {
            const std::size_t L = n.size();
            const std::vector<std::vector<std::size_t>>& use =
                (L == R || chaincols->empty()) ? cols : *chaincols;
            std::vector<T> nop(use.size(), num_traits<T>::from_int(0));
            for (std::size_t j = 0; j < use.size(); ++j)
                for (std::size_t k = 0; k < use[j].size(); ++k)
                    if (use[j][k] < L) nop[j] = T(nop[j] + n[use[j][k]]);
            const std::vector<T> w = f(nop);
            std::vector<T> v(L, num_traits<T>::from_int(1));
            if (w.empty()) return v;
            for (std::size_t j = 0; j < use.size(); ++j) {
                const T wj = w[j < w.size() ? j : w.size() - 1];
                for (std::size_t k = 0; k < use[j].size(); ++k)
                    if (use[j][k] < L) v[use[j][k]] = wj;
            }
            return v;
        };
    }

    /** Spread a per-operand peak rate onto the layer classes. Twin of layerPeak. */
    static std::vector<T> layer_peak(const std::vector<T>& peakPerOperand,
                                     const std::vector<std::vector<std::size_t>>& cols,
                                     std::size_t R) {
        std::vector<T> peak(R, num_traits<T>::from_int(1));
        if (peakPerOperand.empty()) return peak;
        for (std::size_t j = 0; j < cols.size(); ++j) {
            const T pj = peakPerOperand[j < peakPerOperand.size() ? j : peakPerOperand.size() - 1];
            for (std::size_t k = 0; k < cols[j].size(); ++k)
                if (cols[j][k] < R) peak[cols[j][k]] = pj;
        }
        return peak;
    }

    /**
     * Port of buildLayersRecursive.
     *
     * The layer holds a client Delay carrying the callers' think times and the
     * time they spend elsewhere, and `nreplicas` copies of the server station.
     * A class is created for every task, entry, activity and call that the
     * callers' activity graphs reach, and the graph traversal lays down the
     * routing that moves a job between those classes.
     */
    void build_layer(qn::Layer<T>& m, const std::vector<std::size_t>& idxSet,
                     const std::vector<std::size_t>& callers, bool ishostlayer, bool flat) {
        const T one = Tone();
        // The layer key: the model name, the ensemble slot and the column every
        // update map is written under. Under `srvn` it is the layer's only
        // served element; under `flat` it is the first of them.
        const std::size_t idx = idxSet[0];
        m.name = flat ? lqn.hashnames[idx] + ".Flat" : lqn.hashnames[idx];
        m.flat = flat;

        const double rawrepl = lqn.repl[idx];
        std::size_t nreplicas = 1;
        if (!flat && rawrepl > 1.0 && !callers.empty()) {
            // A replicated server is pooled into ONE station instead of one per
            // replica when every caller already addresses all of its replicas:
            // on a host layer that is the callers replicating in step with it,
            // on a task layer it is a declared fan-out covering every replica.
            // With no fan-out declared the lookup is 0 < rawrepl and the layer
            // materialises the replicas, which is the reference's answer too.
            bool reduce = true;
            if (ishostlayer) {
                for (std::size_t c : callers)
                    if (lqn.repl[c] != rawrepl) reduce = false;
            } else {
                for (std::size_t c : callers)
                    if (lqn.fanout_at(c, idx) < rawrepl) reduce = false;
            }
            nreplicas = reduce ? 1 : static_cast<std::size_t>(std::llround(rawrepl));
            if (reduce && !ishostlayer) single_replica_tasks.insert(idx);
        }
        const bool reduce_fanout = (nreplicas == 1 && rawrepl > 1.0 && !callers.empty());
        const std::vector<double>& mult = lqn.maxmult;

        // A client Delay exists unless the layer is a task layer reached only
        // by asynchronous callers, which this port refuses upstream.
        m.clientIdx = m.add_station(qn::Station<T>{"Clients", NodeType::Delay, SchedStrategy::INF,
                                                   std::numeric_limits<double>::infinity(), false, 0});
        const std::size_t clientNode = m.node_of_station(m.clientIdx);
        m.serverIdx = m.clientIdx + 1;
        // One station (times its replicas) per SERVED ELEMENT of this layer.
        // Under `srvn` that is one element and `srv[idx] == server`; under
        // `flat` every processor and called task of the model is here.
        std::vector<std::vector<std::size_t>> srv(lqn.nidx + 1);
        std::vector<std::vector<std::size_t>> srvnode(lqn.nidx + 1);
        m.server_idx_of.assign(lqn.nidx + 1, 0);
        for (std::size_t sidx : idxSet) {
            const bool sishost = sidx <= lqn.nhosts;
            srv[sidx].resize(nreplicas);
            for (std::size_t r = 0; r < nreplicas; ++r) {
                qn::Station<T> st;
                st.name = r == 0 ? lqn.hashnames[sidx]
                                 : lqn.hashnames[sidx] + "." + std::to_string(r + 1);
                // inf-scheduled Queue as Delay node rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_ln.h)
                st.nodetype =
                    lqn.sched[sidx] == SchedStrategy::INF ? NodeType::Delay : NodeType::Queue;
                st.sched = lqn.sched[sidx];
                // setNumberOfServers no-op rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_ln.h)
                st.nservers = lqn.sched[sidx] == SchedStrategy::INF
                                  ? std::numeric_limits<double>::infinity()
                                  : mult[sidx];
                st.attr_ishost = flat ? sishost : ishostlayer;
                st.attr_idx = sidx;
                srv[sidx][r] = m.add_station(st);
                srvnode[sidx].push_back(m.node_of_station(srv[sidx][r]));
            }
            m.server_idx_of[sidx] = srv[sidx][0];
            if (sishost)
                m.host_stations.push_back(srv[sidx][0]);
            else
                m.task_stations.push_back(srv[sidx][0]);
        }
        // the layer's own server, the sole server under `srvn`
        const std::vector<std::size_t>& server = srv[idx];
        const std::vector<std::size_t>& serverNode = srvnode[idx];

        /** Stations of ELEM when it is served in this layer, empty otherwise. */
        auto servers_for = [&](std::size_t elem) -> const std::vector<std::size_t>& {
            static const std::vector<std::size_t> none;
            if (elem >= 1 && elem <= lqn.nidx && !srv[elem].empty()) return srv[elem];
            return none;
        };
        auto server_nodes_for = [&](std::size_t elem) -> const std::vector<std::size_t>& {
            static const std::vector<std::size_t> none;
            if (elem >= 1 && elem <= lqn.nidx && !srvnode[elem].empty()) return srvnode[elem];
            return none;
        };
        /** True when the processor of task TIDX is a server of this layer. */
        auto host_is_server = [&](std::size_t tidx_) {
            return !servers_for(lqn.parent[tidx_]).empty();
        };
        /** True when TIDX calls an entry served by this layer. */
        auto is_layer_client = [&](std::size_t tidx_) {
            for (std::size_t s : idxSet)
                for (std::size_t e : lqn.entriesof[s])
                    if (lqn.issynccaller.get(tidx_, e)) return true;
            return false;
        };
        /** Declare a class at every server station of this layer. */
        auto set_all_servers = [&](std::size_t cl, const Distrib<T>& d) {
            for (std::size_t s : idxSet)
                for (std::size_t st : srv[s]) m.set_service(st, cl, d);
        };

        // Routed call groups, resolved from target ENTRIES to the call indices
        // that reach them: a group is ONE dispatch with n destinations, so its
        // members share a call class and a dispatch class that holds the job at
        // the client while the target is picked (callGroupsByCidx in
        // buildLayersRecursive.m).
        std::vector<std::size_t> group_of_call(lqn.ncalls + 1, 0);
        std::vector<std::vector<std::size_t>> group_members(lqn.callgroups.size() + 1);
        for (std::size_t g = 0; g < lqn.callgroups.size(); ++g) {
            const LqnCallGroup& grp = lqn.callgroups[g];
            for (std::size_t tgt : grp.targets)
                for (std::size_t cidx : lqn.callsof[grp.caller])
                    if (lqn.callpair_dst[cidx] == tgt && lqn.calltype[cidx] == CallType::SYNC &&
                        group_of_call[cidx] == 0) {
                        group_of_call[cidx] = g + 1;
                        group_members[g + 1].push_back(cidx);
                        break;
                    }
        }
        // per group: the Router node, the dispatch class and the group class
        std::vector<std::size_t> grp_router(lqn.callgroups.size() + 1, 0);
        std::vector<std::size_t> grp_dispatch(lqn.callgroups.size() + 1, 0);
        std::vector<std::size_t> grp_class(lqn.callgroups.size() + 1, 0);
        // (router node, dispatch class, group) triples whose strategy is
        // installed once the routing is written
        std::vector<std::array<std::size_t, 3>> routed_group_sites;

        // Fork/Router/Join construction rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_ln.h)
        std::vector<std::size_t> acts_in_caller;
        for (std::size_t c : callers)
            for (std::size_t a : lqn.actsof[c]) acts_in_caller.push_back(a);
        bool hasfork = false, hasjoin = false;
        std::size_t maxfanout = 1;
        for (std::size_t a : acts_in_caller) {
            if (lqn.actposttype[a] == PrecedenceType::POST_AND) hasfork = true;
            if (lqn.actpretype[a] == PrecedenceType::PRE_AND) hasjoin = true;
            std::size_t nand = 0;
            for (std::size_t sx : lqn.graph.succ(a))
                if (lqn.actposttype[sx] == PrecedenceType::POST_AND) ++nand;
            if (nand > maxfanout) maxfanout = nand;
        }
        std::size_t forkNode = 0, joinNode = 0, joinStation = 0;
        std::vector<std::size_t> forkRouter;
        if (hasfork) {
            forkNode = m.add_node("Fork_PostAnd", NodeType::Fork, false);
            for (std::size_t f = 1; f <= maxfanout; ++f)
                forkRouter.push_back(
                    m.add_node("Fork_PostAnd_" + std::to_string(f), NodeType::Router, true));
        }
        if (hasjoin) {
            qn::Station<T> js;
            js.name = "Join_PreAnd";
            js.nodetype = NodeType::Join;
            js.sched = SchedStrategy::INF;
            js.nservers = std::numeric_limits<double>::infinity();
            js.attr_ishost = false;
            js.attr_idx = 0;
            joinStation = m.add_station(js);
            joinNode = m.node_of_station(joinStation);
            if (forkNode) m.fj.emplace_back(forkNode, joinNode);
        }
        // A CACHE LAYER is a HOST layer all of whose callers are cache tasks,
        // buildLayersRecursive.m:88. The Cache node lives here and not in the
        // cache task's own layer: the item lookup is what the task DOES on its
        // processor, so the hit and miss branches have to queue at the same
        // server the read arrived at.
        bool iscachelayer = !flat && ishostlayer && !callers.empty();
        for (std::size_t c : callers)
            if (!lqn.iscache[c]) iscachelayer = false;
        // A SETUP NO LONGER CHANGES HOW A LAYER IS BUILT (buildLayersRecursive.m,
        // 2026-08-11). Wiring the setup/delay-off onto the server station routed
        // the layer through the open M/G/1-with-setup QBD, which reads the idle
        // period off the Poisson rate 1/X and so powers the thread down far more
        // often than a CLOSED layer does -- 12.67% below LDES on lqn_setup -- and
        // it charged the restart to the ACTIVITY, where it is not host demand.
        // The cold start is charged to the ENTRY instead, with the probability
        // that the thread was really found down: see setup_charge().
        std::size_t cacheNode = 0;
        qn::CacheParam<T> cachepar;
        if (iscachelayer) {
            const std::size_t ct = callers[0];
            cachepar.nitems = lqn.nitems[ct];
            cachepar.itemcap = lqn.itemcap[ct];
            cachepar.replacestrat = lqn.replacestrat[ct];
            cacheNode = m.add_node(lqn.hashnames[ct], NodeType::Cache, true);
        }

        // The open streams this layer carries, decided BEFORE the Source is
        // added: m.serverIdx is the client index plus one and replica r is
        // m.serverIdx + r, so a station inserted between them would silently
        // renumber the servers. Source and Sink therefore come last.
        // An async call is declared only in the layer whose server owns the
        // called entry, which is serversFor(parent(dst)) in the reference
        // (buildLayersRecursive.m:327) and the same test route_sync_call uses.
        std::vector<std::size_t> async_here;
        for (std::size_t c : callers)
            for (std::size_t aidx : lqn.actsof[c])
                for (std::size_t cidx : lqn.callsof[aidx])
                    if (lqn.calltype[cidx] == CallType::ASYNC &&
                        !servers_for(lqn.parent[lqn.callpair_dst[cidx]]).empty())
                        async_here.push_back(cidx);
        std::vector<std::size_t> open_entries;
        for (std::size_t c : callers) {
            // An arrival that is the only way into the task is carried by the caller
            // CHAIN, not by a stream: that task has no task layer, so nothing ever sets
            // its surrogate delay, and a chain cycling against an Immediate one plus a
            // stream loads the host twice -- lqn_open_arrival's processor read 0.68
            // against 0.32 from lqns, lqsim and LDES. update_think_times closes the
            // chain on the known rate instead, as it does for a forwarding target.
            if (open_arrival_rate_of(c) > GlobalConstants::FineTol) continue;
            for (std::size_t eidx : lqn.entriesof[c])
                if (lqn.has_arrival[eidx]) open_entries.push_back(eidx);
        }

        std::size_t sourceStation = 0, sourceNode = 0, sinkNode = 0;
        if (!async_here.empty() || !open_entries.empty()) {
            qn::Station<T> src;
            src.name = "Source";
            src.nodetype = NodeType::Source;
            // EXT keeps solver_mva out of both infSET and qSET: the Source is lambda, not a queue
            src.sched = SchedStrategy::EXT;
            src.nservers = 1.0;
            src.attr_ishost = false;
            src.attr_idx = 0;
            sourceStation = m.add_station(src);
            sourceNode = m.node_of_station(sourceStation);
            sinkNode = m.add_node("Sink", NodeType::Sink, false);
            m.sourceIdx = sourceStation;
            m.sinkNode = sinkNode;
        }

        /** The stack of entry classes at the forks currently open. */
        std::vector<std::size_t> forkClassStack;

        // class index of each LQN element and each call, 0 when absent
        std::vector<std::size_t> cls(lqn.nidx + 1, 0);
        std::vector<std::size_t> callcls(lqn.ncalls + 1, 0);
        // `<call>.Aux`, present only for a call whose mean count differs from the
        // replica count; 0 elsewhere. See the routing in route_sync_call.
        std::vector<std::size_t> auxcallcls(lqn.ncalls + 1, 0);

        // A caller carries a closed class when its own PROCESSOR is served here
        // and something drives it (it is a reference task, is called, is a
        // forwarding target or takes an arrival), or when it calls an entry
        // served here. Under `srvn` the first disjunct is exactly the host
        // layer and the second exactly the task layer, so this reduces to the
        // two-branch test it replaces (buildLayersRecursive.m:217).
        auto caller_needs_class = [&](std::size_t tidx_caller) {
            if (host_is_server(tidx_caller)) {
                if (lqn.isref[tidx_caller]) return true;
                for (std::size_t e : lqn.entriesof[tidx_caller])
                    if (lqn.issynccaller.any_col(e) || lqn.isasynccaller.any_col(e) ||
                        is_fwd_target(e) || lqn.has_arrival[e])
                        return true;
            }
            return is_layer_client(tidx_caller);
        };
        auto caller_acts_visible = [&](std::size_t tidx_caller) {
            return host_is_server(tidx_caller) || is_layer_client(tidx_caller);
        };

        // ---- first pass: create the classes --------------------------------
        for (std::size_t tidx_caller : callers) {
            if (caller_needs_class(tidx_caller)) {
                double nj = njobs(tidx_caller, idx);
                if (nj == 0.0) {
                    // representative-replica scaling rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_ln.h)
                    const bool caller_single_replica =
                        reduce_fanout || single_replica_tasks.count(tidx_caller) > 0;
                    nj = caller_single_replica ? mult[tidx_caller]
                                               : mult[tidx_caller] * lqn.repl[tidx_caller];
                    if (std::isinf(nj)) {
                        double s = 0.0;
                        for (std::size_t c = 1; c <= NT(); ++c)
                            if (lqn.taskgraph.get(c, tidx_caller) != Tzero()) s += mult[c];
                        nj = s;
                        if (std::isinf(nj)) {
                            double s2 = 0.0;
                            for (std::size_t c = 1; c <= NT(); ++c)
                                if (std::isfinite(mult[c])) s2 += mult[c] * lqn.repl[c];
                            nj = std::min(s2, 1000.0);
                        }
                    }
                    njobs(tidx_caller, idx) = nj;
                }
                qn::JobClass jc;
                jc.name = lqn.hashnames[tidx_caller];
                jc.type = JobClassType::CLOSED;
                jc.population = nj;
                jc.refstat = m.clientIdx;
                jc.completes = false;
                jc.is_ref_class = true;
                jc.attr_kind = int(LqnElement::TASK);
                jc.attr_idx = tidx_caller;
                cls[tidx_caller] = m.add_class(jc);
                m.attr_tasks.emplace_back(cls[tidx_caller], tidx_caller);
                if (lqn.isref[tidx_caller]) {
                    m.set_service(m.clientIdx, cls[tidx_caller], thinkproc[tidx_caller]);
                } else {
                    // a served task's declared think time is not a per-request
                    // delay, so the seed carries none either; update_think_times
                    // replaces this from the first iteration on
                    m.set_service(m.clientIdx, cls[tidx_caller], Distrib<T>::immediate());
                    thinkt_map.push_back({idx, tidx_caller, m.clientIdx, cls[tidx_caller]});
                }

                for (std::size_t eidx : lqn.entriesof[tidx_caller]) {
                    qn::JobClass ec;
                    ec.name = lqn.hashnames[eidx];
                    ec.type = JobClassType::CLOSED;
                    ec.population = 0.0;
                    ec.refstat = m.clientIdx;
                    ec.completes = false;
                    ec.attr_kind = int(LqnElement::ENTRY);
                    ec.attr_idx = eidx;
                    cls[eidx] = m.add_class(ec);
                    m.attr_entries.emplace_back(cls[eidx], eidx);
                    m.set_service(m.clientIdx, cls[eidx], Distrib<T>::immediate());
                }
            }

            for (std::size_t aidx : lqn.actsof[tidx_caller]) {
                if (caller_acts_visible(tidx_caller)) {
                    qn::JobClass ac;
                    ac.name = lqn.hashnames[aidx];
                    ac.type = JobClassType::CLOSED;
                    ac.population = 0.0;
                    ac.refstat = m.clientIdx;
                    ac.completes = false;
                    ac.attr_kind = int(LqnElement::ACTIVITY);
                    ac.attr_idx = aidx;
                    cls[aidx] = m.add_class(ac);
                    m.attr_activities.emplace_back(cls[aidx], aidx);
                    // The host demand is served at the processor's own station
                    // when this layer holds it; everywhere else the activity is
                    // a surrogate delay at the client.
                    const std::size_t hidx = lqn.parent[lqn.parent[aidx]];
                    if (servers_for(hidx).empty())
                        m.set_service(m.clientIdx, cls[aidx], servtproc[aidx]);
                }
                for (std::size_t cidx : lqn.callsof[aidx]) {
                    if (lqn.calltype[cidx] == CallType::ASYNC) {
                        // An async call is a stream, not a visit: the caller does
                        // not block, so it carries no closed class and instead
                        // drives an open chain of its own out of the Source.
                        const std::size_t adst = lqn.parent[lqn.callpair_dst[cidx]];
                        if (servers_for(adst).empty()) continue;
                        qn::JobClass oc;
                        oc.name = lqn.callhashnames[cidx];
                        oc.type = JobClassType::OPEN;
                        oc.population = std::numeric_limits<double>::infinity();
                        oc.refstat = sourceStation;
                        oc.completes = false;
                        oc.is_ref_class = false;
                        oc.attr_kind = int(LqnElement::CALL);
                        oc.attr_idx = cidx;
                        callcls[cidx] = m.add_class(oc);
                        m.attr_calls.push_back({callcls[cidx], cidx, lqn.callpair_src[cidx],
                                                lqn.callpair_dst[cidx]});
                        // A NEGLIGIBLE seed, not the reference's Immediate. The
                        // rate is replaced from the caller's throughput on the
                        // first post(), but apply_sink_closure weights the
                        // Sink->Source arcs by it in between, and Immediate is
                        // rate 1/0: the closed classes sharing the closure then
                        // get NaN visits. fj_mmt seeds its own open classes the
                        // same way, for the same reason.
                        m.set_service(sourceStation, callcls[cidx],
                                      Distrib<T>::exp_rate(
                                          num_traits<T>::from_double(GlobalConstants::FineTol)));
                        T minRespTA = Tzero();
                        for (std::size_t ta : lqn.actsof[adst]) minRespTA += lqn.hostdem[ta].mean;
                        for (std::size_t st : servers_for(adst)) {
                            m.set_service(st, callcls[cidx], Distrib<T>::exp_mean(minRespTA));
                            call_map.push_back({idx, cidx, st, callcls[cidx]});
                        }
                        arv_call_map.push_back({idx, cidx, sourceStation, callcls[cidx]});
                        continue;
                    }
                    if (lqn.calltype[cidx] != CallType::SYNC) continue;
                    const std::size_t gid = group_of_call[cidx];
                    if (gid != 0) {
                        // ONE dispatch with n destinations. The members SHARE the
                        // dispatch class, which is both the class the strategy
                        // routes and the class that visits the targets: the hop
                        // must not switch class, because a state-dependent
                        // routing function is evaluated at zero off the class
                        // diagonal. The class switch goes on the return arc
                        // instead, into a group class the job continues in.
                        if (grp_dispatch[gid] == 0) {
                            // The strategy is a property of a NODE and routes
                            // over that node's links, not over one class's arcs,
                            // so a dedicated Router whose only links are the
                            // group's targets is the only place where the choice
                            // is exactly the group's.
                            const std::string tag =
                                lqn.hashnames[aidx] + ".Dispatch" + std::to_string(gid);
                            grp_router[gid] = m.add_node(tag + ".Router", NodeType::Router, true);
                            qn::JobClass dc;
                            dc.name = tag;
                            dc.type = JobClassType::CLOSED;
                            dc.population = 0.0;
                            dc.refstat = m.clientIdx;
                            dc.completes = false;
                            dc.attr_kind = int(LqnElement::CALL);
                            dc.attr_idx = cidx;
                            grp_dispatch[gid] = m.add_class(dc);
                            m.set_service(m.clientIdx, grp_dispatch[gid], Distrib<T>::immediate());
                            qn::JobClass gc;
                            gc.name = lqn.callhashnames[cidx] + ".Group" + std::to_string(gid);
                            gc.type = JobClassType::CLOSED;
                            gc.population = 0.0;
                            gc.refstat = m.clientIdx;
                            gc.completes = false;
                            gc.attr_kind = int(LqnElement::CALL);
                            gc.attr_idx = cidx;
                            grp_class[gid] = m.add_class(gc);
                            m.set_service(m.clientIdx, grp_class[gid], Distrib<T>::immediate());
                            routed_group_sites.push_back(
                                {grp_router[gid], grp_dispatch[gid], gid});
                        }
                        callcls[cidx] = grp_dispatch[gid];
                        m.attr_calls.push_back({callcls[cidx], cidx, lqn.callpair_src[cidx],
                                                lqn.callpair_dst[cidx]});
                        for (std::size_t st2 : servers_for(lqn.parent[lqn.callpair_dst[cidx]]))
                            m.set_service(st2, callcls[cidx], callservtproc[cidx]);
                        continue;
                    }
                    qn::JobClass cc;
                    cc.name = lqn.callhashnames[cidx];
                    cc.type = JobClassType::CLOSED;
                    cc.population = 0.0;
                    cc.refstat = m.clientIdx;
                    cc.completes = false;
                    cc.attr_kind = int(LqnElement::CALL);
                    cc.attr_idx = cidx;
                    callcls[cidx] = m.add_class(cc);
                    m.attr_calls.push_back({callcls[cidx], cidx, lqn.callpair_src[cidx],
                                            lqn.callpair_dst[cidx]});
                    // An upper bound on the server's response, replaced at the
                    // first iteration by the measured one. The station seeded is
                    // the CALLEE's under `flat` and the layer's own under
                    // `srvn`, where the latter also seeds calls that leave the
                    // layer -- a phantom rate on a class with no visits here,
                    // kept because removing it moves the fixed point (see
                    // seedCallService in buildLayersRecursive.m).
                    const std::size_t seedidx = flat ? lqn.parent[lqn.callpair_dst[cidx]] : idx;
                    T minRespT = Tzero();
                    for (std::size_t ta : lqn.actsof[seedidx]) minRespT += lqn.hostdem[ta].mean;
                    for (std::size_t st : servers_for(seedidx))
                        m.set_service(st, callcls[cidx], Distrib<T>::exp_mean(minRespT));
                    // .Aux class rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_ln.h)
                    if (lqn.callproc_mean[cidx] !=
                        num_traits<T>::from_double(double(nreplicas))) {
                        qn::JobClass xc;
                        xc.name = lqn.callhashnames[cidx] + ".Aux";
                        xc.type = JobClassType::CLOSED;
                        xc.population = 0.0;
                        xc.refstat = m.clientIdx;
                        xc.completes = false;
                        xc.attr_kind = int(LqnElement::CALL);
                        xc.attr_idx = cidx;
                        auxcallcls[cidx] = m.add_class(xc);
                        m.set_service(m.clientIdx, auxcallcls[cidx], Distrib<T>::immediate());
                    }
                }
            }
        }

        // ---- second pass: routing out of the entries -------------------------
        struct Ctx {
            std::size_t curclass;
            int jobpos;  // 1 = at the client, 2 = at a server replica, 3 = at the Cache
            // The nodes the job stands at when jobpos is atServer. Under `srvn`
            // these are always the layer's own server replicas; under `flat`
            // they are whichever element's stations the last hop landed on, so
            // the position has to be carried rather than assumed.
            std::vector<std::size_t> curnodes;
        };
        const int atClient = 1, atServer = 2, atCache = 3;
        std::vector<int> jobposkey(lqn.nidx + 1, atClient);
        std::vector<std::size_t> curclasskey(lqn.nidx + 1, 0);
        std::vector<std::vector<std::size_t>> curnodeskey(lqn.nidx + 1);

        std::function<Ctx(std::size_t, std::size_t, Ctx)> recur =
            [&](std::size_t tidx_caller, std::size_t aidx, Ctx st) -> Ctx {
            jobposkey[aidx] = st.jobpos;
            curclasskey[aidx] = st.curclass;
            curnodeskey[aidx] = st.curnodes;
            const std::vector<std::size_t> nexts = lqn.graph.succ(aidx);
            std::size_t lastEntryClass = st.curclass;
            // fork pre-state save/restore rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_ln.h)
            bool next_is_fork = false;
            for (std::size_t sx : nexts)
                if (lqn.actposttype[sx] == PrecedenceType::POST_AND) next_is_fork = true;
            const Ctx preFork = st;
            std::vector<std::size_t> andSuccs;
            for (std::size_t sx : nexts)
                if (lqn.actposttype[sx] == PrecedenceType::POST_AND) andSuccs.push_back(sx);

            for (std::size_t k = 0; k < nexts.size(); ++k) {
                const std::size_t nextaidx = nexts[k];
                if (next_is_fork) st = preFork;
                bool isLoop = lqn.graph.get(aidx, nextaidx) != lqn.dag.get(aidx, nextaidx);
                if (lqn.parent[aidx] != lqn.parent[nextaidx]) {
                    // a call to an entry of another task
                    std::size_t cidx = 0;
                    for (std::size_t c : lqn.callsof[aidx])
                        if (lqn.callpair_dst[c] == nextaidx) cidx = c;
                    if (cidx == 0) continue;
                    // FWD lays down no routing: its work is already carried by the pseudo SYNC arc
                    if (lqn.calltype[cidx] != CallType::SYNC) continue;
                    const std::size_t gid = group_of_call[cidx];
                    if (gid != 0) {
                        // The whole group is routed once, at its FIRST member;
                        // the others are the same dispatch and lay down nothing.
                        if (group_members[gid].empty() || group_members[gid][0] != cidx) continue;
                        // The job is switched into the dispatch class while still
                        // at the client, so (client, dispatch) carries exactly
                        // the n arcs the strategy chooses among. The 1/n split
                        // laid down here is the probabilistic reading a solver
                        // without state-dependent routing would see; the
                        // declared strategy replaces it below.
                        std::vector<std::size_t> tnode2, tstat2, tcall2;
                        for (std::size_t mc : group_members[gid]) {
                            const std::size_t tt = lqn.parent[lqn.callpair_dst[mc]];
                            if (servers_for(tt).empty()) continue;
                            tnode2.push_back(server_nodes_for(tt)[0]);
                            tstat2.push_back(servers_for(tt)[0]);
                            tcall2.push_back(mc);
                        }
                        if (tnode2.size() < 2) continue;  // not enough of it is here
                        const std::size_t fromNode2 =
                            st.jobpos == atClient ? clientNode : st.curnodes[0];
                        const std::size_t dc = grp_dispatch[gid], gc = grp_class[gid];
                        m.set_route(st.curclass, dc, fromNode2, grp_router[gid], Tone());
                        const T share2 =
                            T(Tone() / num_traits<T>::from_int(int(tnode2.size())));
                        for (std::size_t d = 0; d < tnode2.size(); ++d) {
                            m.set_route(dc, dc, grp_router[gid], tnode2[d], share2);
                            m.set_route(dc, gc, tnode2[d], clientNode, Tone());
                            m.set_service(tstat2[d], dc, callservtproc[tcall2[d]]);
                            call_map.push_back({idx, tcall2[d], tstat2[d], dc});
                        }
                        st.curclass = gc;
                        st.jobpos = atClient;
                        st.curnodes.clear();
                        continue;
                    }
                    const std::size_t ctgt = lqn.parent[lqn.callpair_dst[cidx]];
                    st = route_sync_call(m, idx, cidx, st, server_nodes_for(ctgt),
                                         servers_for(ctgt), callcls, auxcallcls, clientNode,
                                         atClient, atServer, flat);
                    continue;
                }
                // a successor inside the same task
                bool any_entry_succ = false;
                for (std::size_t sx : nexts)
                    if (sx > lqn.eshift && sx <= lqn.eshift + lqn.nentries) any_entry_succ = true;
                if (!any_entry_succ) {
                    st.jobpos = jobposkey[aidx];
                    st.curclass = curclasskey[aidx];
                    st.curnodes = curnodeskey[aidx];
                } else {
                    if (k > 0 && nexts[k - 1] > lqn.eshift && nexts[k - 1] <= lqn.eshift + lqn.nentries)
                        lastEntryClass = st.curclass;
                    st.jobpos = atClient;
                    st.curclass = lastEntryClass;
                    st.curnodes.clear();
                }
                const T w = lqn.graph.get(aidx, nextaidx);

                // THE CACHE READ, buildLayersRecursive.m:743-765. `aidx` is the
                // item entry and `nextaidx` the activity bound to it: the job
                // goes from the client to the Cache node, which decides the hit
                // or the miss and switches the class accordingly, so the read
                // class itself is never served anywhere and the branch classes
                // pick the work up at the server.
                if (iscachelayer && lqn.nitems[aidx] > 0 && cls[nextaidx] != 0) {
                    const std::size_t readcls = cls[nextaidx];
                    m.set_route(st.curclass, readcls, clientNode, cacheNode, w);
                    if (cachepar.pread.size() < m.classes.size())
                        cachepar.pread.resize(m.classes.size());
                    cachepar.pread[readcls - 1] = lqn.itemproc[aidx];
                    const std::vector<std::size_t> hm = lqn.graph.succ(nextaidx);
                    if (hm.size() != 2)
                        throw InputError("SolverLN: the cache read '" + lqn.names[nextaidx] +
                                         "' needs exactly one hit and one miss successor");
                    if (cachepar.hitclass.size() < m.classes.size()) {
                        cachepar.hitclass.resize(m.classes.size(), 0);
                        cachepar.missclass.resize(m.classes.size(), 0);
                    }
                    cachepar.hitclass[readcls - 1] = cls[hm[0]];
                    cachepar.missclass[readcls - 1] = cls[hm[1]];
                    st.jobpos = atCache;
                    st.curclass = readcls;
                    st.curnodes.clear();
                    st = recur(tidx_caller, nextaidx, st);
                    continue;
                }

                const bool is_and_join_tail = lqn.actpretype[aidx] == PrecedenceType::PRE_AND;
                // the branch index of this successor among the fork's outputs
                std::size_t fbranch = 0;
                if (next_is_fork)
                    for (std::size_t q = 0; q < andSuccs.size(); ++q)
                        if (andSuccs[q] == nextaidx) fbranch = q + 1;

                // Where the successor's HOST DEMAND is served: at the station of
                // its own processor when this layer holds it, at the client
                // otherwise. Under `srvn` that is the host layer's own server
                // and nothing else, so this reduces to the `ishostlayer` test
                // it replaces.
                const std::size_t hidxOf = lqn.parent[lqn.parent[nextaidx]];
                const std::vector<std::size_t>& actStations = servers_for(hidxOf);
                const std::vector<std::size_t>& actNodes = server_nodes_for(hidxOf);
                const bool actAtServer = !actStations.empty();
                const std::size_t from = st.jobpos == atClient   ? clientNode
                                         : st.jobpos == atCache ? cacheNode
                                                                : st.curnodes[0];
                // the node a job continues to after this successor is served
                for (std::size_t r = 0; r < nreplicas; ++r) {
                    const std::size_t fromNode =
                        st.jobpos == atClient   ? clientNode
                        : st.jobpos == atCache ? cacheNode
                                               : st.curnodes[std::min(r, st.curnodes.size() - 1)];
                    const std::size_t toNode = actAtServer ? actNodes[r] : clientNode;
                    if (next_is_fork && fbranch > 0) {
                        m.set_route(st.curclass, st.curclass, fromNode, forkNode, Tone());
                        if (r == 0) forkClassStack.push_back(st.curclass);
                        m.set_route(st.curclass, st.curclass, forkNode, forkRouter[fbranch - 1],
                                    Tone());
                        m.set_route(st.curclass, cls[nextaidx], forkRouter[fbranch - 1], toNode,
                                    Tone());
                    } else if (is_and_join_tail) {
                        // rejoin the class the branch was forked from, then
                        // leave the Join in the successor's class
                        if (forkClassStack.empty())
                            throw InputError("SolverLN: an AND join has no matching fork in '" +
                                             lqn.names[aidx] + "'");
                        const std::size_t forkClass = forkClassStack.back();
                        if (r + 1 == nreplicas) forkClassStack.pop_back();
                        m.set_route(st.curclass, forkClass, fromNode, joinNode, Tone());
                        m.set_route(forkClass, cls[nextaidx], joinNode, toNode, Tone());
                    } else {
                        m.set_route(st.curclass, cls[nextaidx], fromNode, toNode, w);
                    }
                    if (actAtServer)
                        m.set_service(actStations[r], cls[nextaidx], lqn.hostdem[nextaidx]);
                }
                (void)from;
                if (actAtServer) {
                    st.jobpos = atServer;
                    st.curclass = cls[nextaidx];
                    st.curnodes = actNodes;
                    servt_map.push_back({idx, nextaidx, actStations[0], cls[nextaidx]});
                } else {
                    st.jobpos = atClient;
                    st.curclass = cls[nextaidx];
                    st.curnodes.clear();
                    m.set_service(m.clientIdx, cls[nextaidx], servtproc[nextaidx]);
                    thinkt_map.push_back({idx, nextaidx, m.clientIdx, cls[nextaidx]});
                }
                if (aidx != nextaidx && !isLoop) {
                    st = recur(tidx_caller, nextaidx, st);
                    // close the branch with a reply back to the caller's class
                    if (st.jobpos == atClient) {
                        m.set_route(st.curclass, cls[tidx_caller], clientNode, clientNode, Tone());
                    } else {
                        for (std::size_t nd : st.curnodes)
                            m.set_route(st.curclass, cls[tidx_caller], nd, clientNode, Tone());
                    }
                    // .Aux completion-guard rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_ln.h)
                    if (!is_aux_class(m.classes[st.curclass - 1].name))
                        m.classes[st.curclass - 1].completes = true;
                }
            }
            return st;
        };

        for (std::size_t tidx_caller : callers) {
            if (!caller_needs_class(tidx_caller)) continue;
            const std::vector<std::size_t>& ents = lqn.entriesof[tidx_caller];
            const T share = T(Tone() / num_traits<T>::from_int(int(ents.size())));
            for (std::size_t eidx : ents) {
                m.set_route(cls[tidx_caller], cls[eidx], clientNode, clientNode, share);
                if (ents.size() > 1)
                    route_map.push_back({idx, tidx_caller, eidx, m.clientIdx, m.clientIdx,
                                         cls[tidx_caller], cls[eidx]});
                Ctx st{cls[eidx], atClient, {}};
                recur(tidx_caller, eidx, st);
            }
        }

        // ---- open streams, laid down AFTER the activity-graph walk ----------
        // buildLayersRecursive.m:473 puts the entry-arrival routing here for the
        // same reason: the walk above rewrites whole rows of P and would erase
        // an arc written before it. Both kinds route Source -> server -> Sink
        // within one class, so they never join a closed class in a chain --
        // refresh_chains would reject that, the two carrying different refstats.
        if (sourceStation != 0) {
            const T nrep = num_traits<T>::from_int(int(nreplicas));
            for (std::size_t cidx : async_here) {
                const std::size_t oc = callcls[cidx];
                if (oc == 0) continue;
                // the stream enters the CALLEE's station, which under `srvn` is
                // this layer's own server and under `flat` one among many
                const std::vector<std::size_t>& anode =
                    server_nodes_for(lqn.parent[lqn.callpair_dst[cidx]]);
                const T callmean = lqn.callproc_mean[cidx];
                if (callmean < Tone()) {
                    // fewer than one call per firing: a single Bernoulli pass
                    m.set_route(oc, oc, sourceNode, sinkNode, T(Tone() - callmean));
                    for (std::size_t r = 0; r < anode.size(); ++r) {
                        m.set_route(oc, oc, sourceNode, anode[r], T(callmean / nrep));
                        m.set_route(oc, oc, anode[r], sinkNode, Tone());
                    }
                } else {
                    // callmean visits in expectation, as a geometric self-loop
                    const T p = T(Tone() / callmean);
                    for (std::size_t r = 0; r < anode.size(); ++r) {
                        m.set_route(oc, oc, sourceNode, anode[r], T(Tone() / nrep));
                        for (std::size_t q = 0; q < anode.size(); ++q)
                            m.set_route(oc, oc, anode[r], anode[q], T((Tone() - p) / nrep));
                        m.set_route(oc, oc, anode[r], sinkNode, p);
                    }
                }
            }
            for (std::size_t eidx : open_entries) {
                qn::JobClass eo;
                eo.name = lqn.hashnames[eidx] + "_Open";
                eo.type = JobClassType::OPEN;
                eo.population = std::numeric_limits<double>::infinity();
                eo.refstat = sourceStation;
                eo.completes = false;
                eo.is_ref_class = false;
                eo.attr_kind = int(LqnElement::ENTRY);
                eo.attr_idx = eidx;
                const std::size_t ec = m.add_class(eo);
                // the arrival is exogenous and fixed, so it is never reseeded
                m.set_service(sourceStation, ec, lqn.arrival[eidx]);
                // entries are Immediate; the work is the activity bound to them
                std::size_t bound = 0;
                for (std::size_t sx : lqn.graph.succ(eidx))
                    if (bound == 0 && sx > lqn.ashift) bound = sx;
                const Distrib<T>& svc = bound != 0 ? servtproc[bound] : servtproc[eidx];
                // the arrival enters the processor of the entry's task under
                // host layering, and the task's own station under `flat`
                const std::vector<std::size_t>& ostat =
                    flat ? servers_for(lqn.parent[eidx]) : server;
                const std::vector<std::size_t>& onode =
                    flat ? server_nodes_for(lqn.parent[eidx]) : serverNode;
                for (std::size_t r = 0; r < ostat.size(); ++r) {
                    m.set_service(ostat[r], ec, svc);
                    m.set_route(ec, ec, sourceNode, onode[r], T(Tone() / nrep));
                    m.set_route(ec, ec, onode[r], sinkNode, Tone());
                }
            }
        }

        // ---- admission constraint on the server station ---------------------
        // buildLayersRecursive.m:596-625. The constraint is stated over LQN
        // elements and has to be re-expressed in the layer's own classes: a
        // task occupies its host through the classes of its ACTIVITIES, an
        // entry is occupied by the classes of the CALLS that target it. The
        // columns are summed into the class, never assigned, because several
        // classes can stand for one column.
        if (idx < lqn.lincon_A.size() && lqn.lincon_A[idx].rows() > 0) {
            const Matrix<T>& Aelem = lqn.lincon_A[idx];
            Matrix<T> Alayer(Aelem.rows(), m.classes.size(), Tzero());
            const std::vector<std::size_t>& constrained =
                ishostlayer ? lqn.tasksof[idx] : lqn.entriesof[idx];
            bool any = false;
            for (std::size_t j = 0; j < constrained.size() && j < Aelem.cols(); ++j) {
                std::vector<std::size_t> layerClasses;
                if (ishostlayer) {
                    for (std::size_t a : lqn.actsof[constrained[j]])
                        if (cls[a] != 0) layerClasses.push_back(cls[a]);
                } else {
                    for (std::size_t c = 1; c <= lqn.ncalls; ++c)
                        if (lqn.callpair_dst[c] == constrained[j] && callcls[c] != 0)
                            layerClasses.push_back(callcls[c]);
                }
                for (std::size_t k = 0; k < layerClasses.size(); ++k)
                    for (std::size_t rr = 0; rr < Aelem.rows(); ++rr) {
                        Alayer(rr, layerClasses[k] - 1) =
                            T(Alayer(rr, layerClasses[k] - 1) + Aelem(rr, j));
                        if (Aelem(rr, j) != Tzero()) any = true;
                    }
            }
            if (any) {
                // ONE region over every replica, not one each: the constraint
                // models a passive resource of the server as a whole (a
                // semaphore, a connection pool), so the replicas share tokens.
                typename qn::NetworkStruct<T>::Region rg;
                const std::size_t M = m.stations.size(), K = m.classes.size();
                rg.cap.assign(M, std::vector<double>(K + 1, -1.0));
                rg.maxmem.assign(M, -1.0);
                rg.members.assign(M, false);
                rg.rule.assign(K, lang::DropStrategy::WAITQ);
                rg.weight.assign(K, Tone());
                rg.size.assign(K, Tone());
                for (std::size_t r = 0; r < nreplicas; ++r) rg.members[server[r] - 1] = true;
                rg.lincon_A = Alayer;
                rg.lincon_b = lqn.lincon_b[idx];
                m.regions.push_back(rg);
            }
        }

        // The Cache node's parameters are only complete once the walk has seen
        // every read: pread, hitclass and missclass are all per READ CLASS, and
        // the classes are created as the activity graph is traversed.
        if (cacheNode != 0) {
            cachepar.pread.resize(m.classes.size());
            cachepar.hitclass.resize(m.classes.size(), 0);
            cachepar.missclass.resize(m.classes.size(), 0);
            m.nodeparam[cacheNode] = cachepar;
        }

        // The declared dispatch replaces the probabilistic split on the router,
        // and ONLY on the (router, dispatch class) pair whose arcs are exactly
        // the group's targets. The split stays in P underneath, which is what
        // `refresh_routing` re-expands for anything that reads a matrix; the SSA
        // engine walks the arcs itself and honours the strategy.
        for (const std::array<std::size_t, 3>& site : routed_group_sites) {
            qn::NodeDef& nd = m.nodes[site[0] - 1];
            if (nd.routing.size() < m.classes.size())
                nd.routing.resize(m.classes.size(), RoutingStrategy::PROB);
            nd.routing[site[1] - 1] = lqn.callgroups[site[2] - 1].strategy;
        }

        // ---- queue-dependent service rates on the server station ------------
        // buildLayersRecursive.m:700-762. Declared over the server's OPERANDS
        // (the tasks of a host, the entries of a task) and re-expressed in the
        // layer's own classes, on the same mapping the admission constraint
        // uses above. The handles are evaluated in two index spaces -- CTMC and
        // the exact recursions pass a per-class vector, the AMVA and NC chain
        // recursions a per-chain one -- so each carries both column lists and
        // picks by the length of what it is handed. The chain list is only known
        // after refresh_chains, hence the shared slot filled just below.
        std::vector<std::pair<std::shared_ptr<std::vector<std::vector<std::size_t>>>,
                              std::vector<std::vector<std::size_t>>>> deferred_chaincols;
        for (std::size_t sidx : idxSet) {
            const bool hasld = sidx < lqn.lldscaling.size() && !lqn.lldscaling[sidx].empty();
            const bool hascd = sidx < lqn.cdscaling.size() && bool(lqn.cdscaling[sidx]);
            const bool hasjd = sidx < lqn.jdscaling.size() && bool(lqn.jdscaling[sidx]);
            const bool haspools = sidx < lqn.pools.size() && !lqn.pools[sidx].empty();
            if (!(hasld || hascd || hasjd || haspools)) continue;
            const bool sishost = sidx <= lqn.nhosts;
            const std::vector<std::size_t>& operandIdx =
                sishost ? lqn.tasksof[sidx] : lqn.entriesof[sidx];
            std::vector<std::vector<std::size_t>> cols(operandIdx.size());
            for (std::size_t j = 0; j < operandIdx.size(); ++j) {
                if (sishost) {
                    for (std::size_t a : lqn.actsof[operandIdx[j]])
                        if (cls[a] != 0) cols[j].push_back(cls[a] - 1);
                } else {
                    for (std::size_t c = 1; c <= lqn.ncalls; ++c)
                        if (lqn.callpair_dst[c] == operandIdx[j] && callcls[c] != 0)
                            cols[j].push_back(callcls[c] - 1);
                }
            }
            const std::size_t R = m.classes.size();
            std::shared_ptr<std::vector<std::vector<std::size_t>>> chaincols =
                std::make_shared<std::vector<std::vector<std::size_t>>>();
            deferred_chaincols.push_back(std::make_pair(chaincols, cols));
            for (std::size_t r = 0; r < nreplicas; ++r) {
                // The scalings are written straight onto the Station, as the
                // region block above writes m.regions: a Layer is a
                // NetworkStruct, not the NetworkBuilder that carries the
                // set_*_dependence helpers.
                qn::Station<T>& stn = m.stations[srv[sidx][r] - 1];
                if (hasld) stn.lldscaling = lqn.lldscaling[sidx];
                if (hascd) {
                    // beta_{i,r} is product form only while an operand maps to a
                    // single class; where it aggregates several, the same
                    // scaling is emitted as a joint dependence, numerically
                    // identical but no longer exact.
                    bool one_class_each = true;
                    for (std::size_t j = 0; j < cols.size(); ++j)
                        if (cols[j].size() > 1) one_class_each = false;
                    const lang::CdScaling<T> h =
                        layer_dep_handle(lqn.cdscaling[sidx], cols, chaincols, R);
                    const std::vector<T> pk = layer_peak(lqn.cdscalingpeak[sidx], cols, R);
                    if (one_class_each) {
                        stn.cdscaling = h;
                        stn.cdscalingpeak = pk;
                    } else {
                        stn.jdscaling = h;
                        stn.jdscalingpeak = pk;
                    }
                }
                if (hasjd) {
                    stn.jdscaling = layer_dep_handle(lqn.jdscaling[sidx], cols, chaincols, R);
                    stn.jdscalingpeak = layer_peak(lqn.jdscalingpeak[sidx], cols, R);
                }
                if (haspools) {
                    // A compatibility declaration IS a rate law: the pools clear
                    // the activated-server rate of api::sn_compat_rate, order
                    // independent at every integer state. sn_compat_scaling
                    // normalises it against the rate the SAME population would
                    // get under full compatibility, so eta isolates the
                    // compatibility GRAPH and a fully-compatible pool is the
                    // neutral eta == 1; the low-occupancy loss stays with the
                    // solver's own multiserver term.
                    const lqn::ServerPools<T>& pl = lqn.pools[sidx];
                    const Matrix<T> compat = pl.compat;
                    const std::vector<double> counts = pl.counts;
                    const std::vector<T> rates = pl.rates;
                    lang::CdScaling<T> etaPool = [compat, counts,
                                                  rates](const std::vector<T>& nop) {
                        return std::vector<T>(
                            1, api::sn_compat_scaling(compat, counts, rates, nop));
                    };
                    stn.jdscaling = layer_dep_handle(etaPool, cols, chaincols, R);
                    stn.jdscalingpeak = layer_peak(
                        std::vector<T>(cols.size(), num_traits<T>::from_int(1)), cols, R);
                }
            }
        }

        m.refresh_chains();
        // The chain columns of every handle emitted above, now that the chains
        // exist. The slots are shared with the handles, so filling them here
        // reaches every replica without rebuilding a single lambda.
        for (std::size_t d = 0; d < deferred_chaincols.size(); ++d) {
            const std::vector<std::vector<std::size_t>>& cols = deferred_chaincols[d].second;
            std::vector<std::vector<std::size_t>>& out = *deferred_chaincols[d].first;
            out.assign(cols.size(), std::vector<std::size_t>());
            for (std::size_t j = 0; j < cols.size(); ++j) {
                for (std::size_t k = 0; k < cols[j].size(); ++k)
                    for (std::size_t cc = 0; cc < m.chains.size(); ++cc)
                        if (cols[j][k] < m.chains[cc].size() && m.chains[cc][cols[j][k]]) {
                            bool seen = false;
                            for (std::size_t q = 0; q < out[j].size(); ++q)
                                if (out[j][q] == cc) seen = true;
                            if (!seen) out[j].push_back(cc);
                        }
            }
        }
        // A layer holding a POST_AND fork is a fork-join model like any other,
        // and the reference builds it as a `Network`, so buildLayers's
        // `getStruct()` runs the MMT node-visit pass on it: layer P:P2 of
        // lqn_workflows reports a Join visit of 4.3333, not the 1 the routing
        // solve leaves. THE LATER `refresh_chains()` CALLS DO NOT REPEAT IT, and
        // that is the reference's behaviour rather than an oversight here:
        // SolverLN.m re-runs `refreshChains()` alone after a routing reset,
        // which recomputes the visits WITHOUT the pass.
        api::sn_fj_nodevisits_mmt(m);
    }

    /**
     * Port of routeSynchCall (buildLayersRecursive.m).
     *
     * A synchronous call is made `callmean` times per visit to the calling
     * activity, against `nreplicas` server replicas. The branch is chosen by
     * callmean ALONE, against 1 and never against nreplicas: a Bernoulli pass
     * carries at most one call, more than one needs the geometric loop through
     * the `<call>.Aux` class, and the replicas only SPLIT each probability by
     * ntgt -- they never relax that bound (buildLayersRecursive.m:1028-1051).
     * Comparing against nreplicas instead agrees with the reference only at
     * nreplicas == 1; on lqn_sockshop's replicated P2_1 it put a callmean of 1
     * on the < branch and doubled T1's throughput and P2_1's utilization.
     *
     * A second class is needed because a self-loop on the call class would also
     * re-enter its service. Which of the two carries the call time depends on
     * the direction:
     *
     *   callmean < 1           the mean count is folded into the DEMAND
     *                          (callservt = callmean * W), so the call class must
     *                          be visited exactly once or the time is discounted
     *                          twice; .Aux absorbs the remaining branch.
     *   callmean > 1           the call class is re-entered, and .Aux is the
     *                          return path that closes the loop with probability
     *                          1/callmean, giving a mean of callmean visits.
     *
     * The four cases below are (job at client / at server) x (call targets an
     * entry of THIS layer's server / of some other task).
     */
    struct CtxPair {
        std::size_t curclass;
        int jobpos;
    };

    /** The `<call>.Aux` suffix, which is how the reference identifies them too. */
    static bool is_aux_class(const std::string& name) {
        return name.size() >= 4 && name.compare(name.size() - 4, 4, ".Aux") == 0;
    }
    template <class Ctx>
    Ctx route_sync_call(qn::Layer<T>& m, std::size_t idx, std::size_t cidx, Ctx st,
                        const std::vector<std::size_t>& tnode,
                        const std::vector<std::size_t>& tstat,
                        const std::vector<std::size_t>& callcls,
                        const std::vector<std::size_t>& auxcallcls, std::size_t clientNode,
                        int atClient, int atServer, bool flat) {
        const T one = Tone();
        // The callee's stations in THIS layer, empty when it is served
        // elsewhere: under `srvn` that is the layer's own server or nothing,
        // under `flat` it is whichever of the many servers the call targets.
        const std::size_t ntgt = tnode.size();
        const bool to_this_server = ntgt > 0;
        const T nrep = num_traits<T>::from_int(int(ntgt > 0 ? ntgt : 1));
        const T share = T(one / nrep);
        const T callmean = lqn.callproc_mean[cidx];
        const std::size_t cc = callcls[cidx];
        const std::size_t ax = auxcallcls[cidx];
        const bool below = callmean < one;
        const bool above = callmean > one;
        if (st.jobpos == atClient) {
            if (to_this_server) {
                if (below) {
                    m.set_route(st.curclass, ax, clientNode, clientNode, T(one - callmean));
                    for (std::size_t r = 0; r < ntgt; ++r) {
                        m.set_route(st.curclass, cc, clientNode, tnode[r], T(callmean / nrep));
                        m.set_route(cc, cc, tnode[r], clientNode, one);
                    }
                    // keeps .Aux attached to the routing graph; carries no time
                    m.set_route(ax, cc, clientNode, clientNode, one);
                } else if (above) {
                    for (std::size_t r = 0; r < ntgt; ++r) {
                        m.set_route(st.curclass, cc, clientNode, tnode[r], share);
                        m.set_route(cc, ax, tnode[r], clientNode, one);
                        m.set_route(ax, cc, clientNode, tnode[r], T((one - one / callmean) / nrep));
                    }
                    m.set_route(ax, cc, clientNode, clientNode, T(one / callmean));
                } else {
                    for (std::size_t r = 0; r < ntgt; ++r) {
                        m.set_route(st.curclass, cc, clientNode, tnode[r], share);
                        m.set_route(cc, cc, tnode[r], clientNode, one);
                    }
                }
                for (std::size_t r = 0; r < ntgt; ++r) {
                    m.set_service(tstat[r], cc, callservtproc[cidx]);
                    call_map.push_back({idx, cidx, tstat[r], cc});
                }
                m.set_service(m.clientIdx, cc, Distrib<T>::immediate());
                st.jobpos = atClient;
                st.curnodes.clear();
                st.curclass = cc;
            } else {
                m.set_route(st.curclass, cc, clientNode, clientNode, one);
                if (below || above) {
                    m.set_route(cc, ax, clientNode, clientNode, one);
                    st.curclass = ax;
                } else {
                    st.curclass = cc;
                }
                st.jobpos = atClient;
                st.curnodes.clear();
                m.set_service(m.clientIdx, cc, callservtproc[cidx]);
                call_map.push_back({idx, cidx, m.clientIdx, cc});
            }
        } else {
            // the node the job stands at, one per replica of wherever it landed
            auto fromNode = [&](std::size_t r) {
                return st.curnodes[std::min(r, st.curnodes.size() - 1)];
            };
            if (to_this_server) {
                if (below) {
                    // The skip and the call merge back in the CALL class at the
                    // client, exactly as in the atClient branch above. Routing
                    // the skip into the call class instead and leaving in .Aux
                    // gives .Aux no inbound arc at all, so its chain has no
                    // reference class (buildLayersRecursive.m:1100-1118). This
                    // port did that until 2026-08-11 and the arm is the same
                    // under both layerings, as it is in the reference.
                    for (std::size_t r = 0; r < ntgt; ++r) {
                        m.set_route(st.curclass, ax, fromNode(r), clientNode, T(one - callmean));
                        m.set_route(st.curclass, cc, fromNode(r), tnode[r], T(callmean / nrep));
                        m.set_route(cc, cc, tnode[r], clientNode, one);
                    }
                    m.set_route(ax, cc, clientNode, clientNode, one);
                    m.set_service(m.clientIdx, cc, Distrib<T>::immediate());
                    st.jobpos = atClient;
                    st.curnodes.clear();
                    st.curclass = cc;
                } else if (above) {
                    if (flat) {
                        // the geometric repeat transits the client between
                        // visits; a self-loop would merge them into one
                        for (std::size_t r = 0; r < ntgt; ++r) {
                            m.set_route(st.curclass, cc, fromNode(r), tnode[r], one);
                            m.set_route(cc, ax, tnode[r], clientNode, one);
                            m.set_route(ax, cc, clientNode, tnode[r],
                                        T((one - one / callmean) / nrep));
                        }
                        m.set_route(ax, cc, clientNode, clientNode, T(one / callmean));
                        m.set_service(m.clientIdx, cc, Distrib<T>::immediate());
                        st.curclass = cc;
                    } else {
                        for (std::size_t r = 0; r < ntgt; ++r) {
                            m.set_route(st.curclass, cc, fromNode(r), tnode[r], one);
                            m.set_route(cc, cc, tnode[r], tnode[r], T(one - one / callmean));
                            m.set_route(cc, ax, tnode[r], clientNode, T(one / callmean));
                        }
                        st.curclass = ax;
                    }
                    st.jobpos = atClient;
                    st.curnodes.clear();
                } else {
                    for (std::size_t r = 0; r < ntgt; ++r)
                        m.set_route(st.curclass, cc, fromNode(r), tnode[r], one);
                    if (flat) {
                        // the reply returns the job to the client, which is
                        // where the successor restoration expects it
                        for (std::size_t r = 0; r < ntgt; ++r)
                            m.set_route(cc, cc, tnode[r], clientNode, one);
                        m.set_service(m.clientIdx, cc, Distrib<T>::immediate());
                        st.jobpos = atClient;
                        st.curnodes.clear();
                    } else {
                        st.jobpos = atServer;
                        st.curnodes = tnode;
                    }
                    st.curclass = cc;
                }
                for (std::size_t r = 0; r < ntgt; ++r) {
                    m.set_service(tstat[r], cc, callservtproc[cidx]);
                    call_map.push_back({idx, cidx, tstat[r], cc});
                }
            } else {
                for (std::size_t nd : st.curnodes)
                    m.set_route(st.curclass, cc, nd, clientNode, one);
                if (below || above) {
                    m.set_route(cc, ax, clientNode, clientNode, one);
                    st.curclass = ax;
                } else {
                    st.curclass = cc;
                }
                st.jobpos = atClient;
                st.curnodes.clear();
                m.set_service(m.clientIdx, cc, callservtproc[cidx]);
                call_map.push_back({idx, cidx, m.clientIdx, cc});
            }
        }
        return st;
    }

    // -----------------------------------------------------------------------
    // getEntryServiceMatrix
    // -----------------------------------------------------------------------

    /**
     * Reachability of activities and calls from each entry, as a 0/1 matrix over
     * the combined element-and-call index space. Multiplying it by
     * [residt; callresidt] sums an entry's whole service into one number.
     */
    void build_entry_service_matrix() {
        const std::size_t dim = lqn.nidx + lqn.ncalls;
        servtmatrix = Matrix<T>(dim + 1, dim + 1, Tzero());
        std::function<void(std::size_t, std::size_t)> rec = [&](std::size_t aidx, std::size_t eidx) {
            for (std::size_t nextaidx : lqn.graph.succ(aidx)) {
                const bool isLoop = lqn.graph.get(aidx, nextaidx) != lqn.dag.get(aidx, nextaidx);
                if (lqn.parent[aidx] != lqn.parent[nextaidx]) {
                    for (std::size_t cidx : lqn.callsof[aidx])
                        if (lqn.calltype[cidx] == CallType::SYNC)
                            servtmatrix(eidx, lqn.nidx + cidx) = Tone();
                } else if (nextaidx != aidx && !isLoop) {
                    servtmatrix(eidx, nextaidx) = Tone();
                    rec(nextaidx, eidx);
                }
            }
        };
        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            rec(eidx, eidx);
        }
    }

    // -----------------------------------------------------------------------
    // initInterlock
    // -----------------------------------------------------------------------

    /** Port of initInterlock: the LQNS V5 static interlock analysis. */
    void init_interlock() {
        const std::size_t NE = lqn.nentries;
        il_all = Matrix<T>(NE + 1, NE + 1, Tzero());
        il_ph1 = Matrix<T>(NE + 1, NE + 1, Tzero());

        std::function<void(std::size_t, std::size_t, T, T, std::vector<bool>&, int)> trace =
            [&](std::size_t eidx, std::size_t root_e, T pall, T pph1, std::vector<bool>& visited,
                int depth) {
                if (eidx <= lqn.eshift || eidx > lqn.eshift + NE) return;
                const std::size_t e = eidx - lqn.eshift;
                if (visited[e]) return;
                visited[e] = true;
                il_all(root_e, e) = T(il_all(root_e, e) + pall);
                il_ph1(root_e, e) = T(il_ph1(root_e, e) + pph1);
                for (std::size_t aidx : lqn.actsof[eidx]) {
                    if (aidx <= lqn.ashift || aidx > lqn.ashift + lqn.nacts) continue;
                    const std::size_t a = aidx - lqn.ashift;
                    if (depth > 0 && lqn.actphase[a] > 1) continue;
                    const bool is_ph1 = lqn.actphase[a] <= 1;
                    for (std::size_t cidx : lqn.callsof[aidx]) {
                        if (lqn.calltype[cidx] != CallType::SYNC) continue;
                        if (!(lqn.callproc_mean[cidx] > Tzero())) continue;
                        const std::size_t dst = lqn.callpair_dst[cidx];
                        if (dst <= lqn.eshift || dst > lqn.eshift + NE) continue;
                        trace(dst, root_e, T(pall * lqn.callproc_mean[cidx]),
                              is_ph1 ? T(pph1 * lqn.callproc_mean[cidx]) : Tzero(), visited,
                              depth + 1);
                    }
                }
                visited[e] = false;
            };
        for (std::size_t e = 1; e <= NE; ++e) {
            std::vector<bool> visited(NE + 1, false);
            trace(lqn.eshift + e, e, Tone(), Tone(), visited, 0);
        }

        il_common_entries.assign(NT() + 1, {});
        il_src_all.assign(NT() + 1, {});
        il_src_ph2.assign(NT() + 1, {});
        il_num_sources.assign(NT() + 1, 0.0);

        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (lqn.isref[tidx] || lqn.sched[tidx] == SchedStrategy::INF) continue;
            interlock_for_server(tidx);
        }
        for (std::size_t h = 1; h <= lqn.nhosts; ++h) {
            if (lqn.sched[h] == SchedStrategy::INF) continue;
            interlock_for_server(h);
        }
    }

    std::vector<std::size_t> server_entry_nums(std::size_t serverIdx) const {
        std::vector<std::size_t> out;
        if (serverIdx <= lqn.nhosts) {
            for (std::size_t tidx : lqn.tasksof[serverIdx])
                for (std::size_t se : lqn.entriesof[tidx]) out.push_back(se - lqn.eshift);
        } else {
            for (std::size_t se : lqn.entriesof[serverIdx]) out.push_back(se - lqn.eshift);
        }
        return out;
    }

    std::vector<std::size_t> client_tasks(std::size_t serverIdx) const {
        std::vector<std::size_t> out;
        if (serverIdx <= lqn.nhosts) return lqn.tasksof[serverIdx];
        for (std::size_t se : lqn.entriesof[serverIdx])
            for (std::size_t ci : lqn.iscaller.col(se))
                if (ci > lqn.tshift && ci <= NT()) out.push_back(ci);
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    std::vector<std::size_t> call_dst_tasks(std::size_t src_eidx, std::size_t target_e) const {
        std::vector<std::size_t> out;
        for (std::size_t aidx : lqn.actsof[src_eidx]) {
            if (aidx <= lqn.ashift || aidx > lqn.ashift + lqn.nacts) continue;
            for (std::size_t cidx : lqn.callsof[aidx]) {
                if (lqn.calltype[cidx] != CallType::SYNC) continue;
                const std::size_t dst = lqn.callpair_dst[cidx];
                const std::size_t de = dst - lqn.eshift;
                if (de >= 1 && de <= lqn.nentries && il_all(de, target_e) > Tzero())
                    out.push_back(lqn.parent[dst]);
            }
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    bool is_branch_point(std::size_t srcX, std::size_t entryA, std::size_t srcY,
                         std::size_t entryB) const {
        const std::size_t taskA = lqn.parent[entryA], taskB = lqn.parent[entryB];
        const std::size_t taskX = lqn.parent[srcX];
        if (taskX == taskA && taskX == taskB) return false;
        if (srcX == entryA || srcY == entryB) return true;
        const std::vector<std::size_t> dx = call_dst_tasks(srcX, entryA - lqn.eshift);
        const std::vector<std::size_t> dy = call_dst_tasks(srcY, entryB - lqn.eshift);
        for (std::size_t a : dx)
            for (std::size_t b : dy)
                if (a != b) return true;
        return false;
    }

    void trace_to_server(std::size_t eidx, std::size_t serverIdx, std::vector<bool>& visited,
                         std::vector<std::size_t>& itasks, bool isHead) const {
        if (eidx <= lqn.eshift || eidx > lqn.eshift + lqn.nentries) return;
        const std::size_t e = eidx - lqn.eshift;
        if (visited[e]) return;
        const std::size_t owner = lqn.parent[eidx];
        if (owner == serverIdx) return;
        if (serverIdx <= lqn.nhosts && lqn.parent[owner] == serverIdx) return;
        visited[e] = true;
        bool found = false;
        const std::vector<std::size_t> sen = server_entry_nums(serverIdx);
        for (std::size_t aidx : lqn.actsof[eidx]) {
            if (aidx <= lqn.ashift || aidx > lqn.ashift + lqn.nacts) continue;
            for (std::size_t cidx : lqn.callsof[aidx]) {
                if (lqn.calltype[cidx] != CallType::SYNC) continue;
                const std::size_t dst = lqn.callpair_dst[cidx];
                const std::size_t dtask = lqn.parent[dst];
                bool reaches = dtask == serverIdx ||
                               (serverIdx <= lqn.nhosts && lqn.parent[dtask] == serverIdx);
                if (!reaches) {
                    const std::size_t de = dst - lqn.eshift;
                    for (std::size_t sn : sen)
                        if (il_all(de, sn) > Tzero()) reaches = true;
                }
                if (reaches) {
                    trace_to_server(dst, serverIdx, visited, itasks, false);
                    found = true;
                }
            }
        }
        if (found && !isHead) {
            itasks.push_back(owner);
            std::sort(itasks.begin(), itasks.end());
            itasks.erase(std::unique(itasks.begin(), itasks.end()), itasks.end());
        }
        visited[e] = false;
    }

    void interlock_for_server(std::size_t serverIdx) {
        const std::vector<std::size_t> sen = server_entry_nums(serverIdx);
        if (sen.empty()) return;
        const std::vector<std::size_t> cts = client_tasks(serverIdx);
        if (cts.empty()) return;

        std::vector<std::pair<std::size_t, std::size_t>> pairs;  // (task, entry number)
        for (std::size_t ct : cts)
            for (std::size_t ce : lqn.entriesof[ct]) {
                const std::size_t cen = ce - lqn.eshift;
                if (cen < 1 || cen > lqn.nentries) continue;
                for (std::size_t se : sen)
                    if (il_all(cen, se) > Tzero()) {
                        pairs.emplace_back(ct, cen);
                        break;
                    }
            }
        if (pairs.size() < 2) return;

        std::vector<std::size_t> common;
        for (std::size_t i = 0; i < pairs.size(); ++i)
            for (std::size_t j = i + 1; j < pairs.size(); ++j) {
                if (pairs[i].first == pairs[j].first) continue;
                const std::size_t eA = pairs[i].second, eC = pairs[j].second;
                for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
                    const std::size_t tidx = lqn.tshift + t;
                    for (std::size_t ex : lqn.entriesof[tidx])
                        for (std::size_t ey : lqn.entriesof[tidx]) {
                            const std::size_t xn = ex - lqn.eshift, yn = ey - lqn.eshift;
                            if (xn < 1 || yn < 1 || xn > lqn.nentries || yn > lqn.nentries) continue;
                            if (il_all(xn, eA) > Tzero() && il_all(yn, eC) > Tzero() &&
                                is_branch_point(ex, eA + lqn.eshift, ey, eC + lqn.eshift))
                                common.push_back(ex);
                        }
                }
            }
        std::sort(common.begin(), common.end());
        common.erase(std::unique(common.begin(), common.end()), common.end());
        if (common.empty()) return;

        std::vector<std::size_t> interlocked;
        for (std::size_t ce : common) {
            std::vector<bool> visited(lqn.nentries + 1, false);
            std::vector<std::size_t> it;
            trace_to_server(ce, serverIdx, visited, it, true);
            for (std::size_t x : it) interlocked.push_back(x);
        }
        std::sort(interlocked.begin(), interlocked.end());
        interlocked.erase(std::unique(interlocked.begin(), interlocked.end()), interlocked.end());

        std::vector<std::size_t> src_all;
        for (std::size_t ce : common) src_all.push_back(lqn.parent[ce]);
        std::sort(src_all.begin(), src_all.end());
        src_all.erase(std::unique(src_all.begin(), src_all.end()), src_all.end());
        {
            std::vector<std::size_t> diff;
            for (std::size_t x : src_all)
                if (!std::binary_search(interlocked.begin(), interlocked.end(), x)) diff.push_back(x);
            src_all = diff;
        }

        // Left empty on purpose, and NOT because phase 2 is unsupported: the
        // reference computes ph2SrcTasks (initInterlock.m:230-262) and never
        // reads it, and its phase-2 source count is behind `if false`
        // (initInterlock.m:288-304). Filling it would change no number.
        std::vector<std::size_t> src_ph2;
        for (std::size_t it : interlocked)
            for (std::size_t ie : lqn.entriesof[it])
                for (std::size_t ci : lqn.iscaller.col(ie))
                    if (ci > lqn.tshift && ci <= NT() &&
                        !std::binary_search(interlocked.begin(), interlocked.end(), ci))
                        src_all.push_back(ci);
        std::sort(src_all.begin(), src_all.end());
        src_all.erase(std::unique(src_all.begin(), src_all.end()), src_all.end());

        double nsrc = 0.0;
        for (std::size_t st : src_all) nsrc += lqn.mult[st];

        il_common_entries[serverIdx] = common;
        il_src_all[serverIdx] = src_all;
        il_src_ph2[serverIdx] = src_ph2;
        il_num_sources[serverIdx] = nsrc;
    }

    // -----------------------------------------------------------------------
    // init and the outer iteration
    // -----------------------------------------------------------------------
    void init() {
        const std::size_t N = lqn.nidx;
        tput.assign(N + 1, Tzero());
        util.assign(N + 1, Tzero());
        servt.assign(N + 1, Tzero());
        residt.assign(N + 1, Tzero());
        thinkt.assign(N + 1, Tzero());
        callservt.assign(lqn.ncalls + 1, Tzero());
        callresidt.assign(lqn.ncalls + 1, Tzero());
        tputproc.assign(N + 1, Distrib<T>::disabled_dist());
        servt_ph1.assign(N + 1, Tzero());
        servt_ph2.assign(N + 1, Tzero());
        prOvertake.assign(lqn.nentries + 1, Tzero());
        build_entry_service_matrix();

        relax_omega = (opt.relax == "fixed" || opt.relax == "adaptive") ? opt.relax_factor : 1.0;

        servt_prev.assign(N + 1, std::numeric_limits<double>::quiet_NaN());
        residt_prev.assign(N + 1, std::numeric_limits<double>::quiet_NaN());
        tput_prev.assign(N + 1, std::numeric_limits<double>::quiet_NaN());
        thinkt_prev.assign(N + 1, std::numeric_limits<double>::quiet_NaN());
        callservt_prev.assign(lqn.ncalls + 1, std::numeric_limits<double>::quiet_NaN());
        callresidt_prev.assign(lqn.ncalls + 1, std::numeric_limits<double>::quiet_NaN());
        servt_prev_v.assign(N + 1, Tzero());
        residt_prev_v.assign(N + 1, Tzero());
        tput_prev_v.assign(N + 1, Tzero());
        thinkt_prev_v.assign(N + 1, Tzero());
        callservt_prev_v.assign(lqn.ncalls + 1, Tzero());

        unique_route_idx.clear();
        for (const RouteRow& r : route_map) unique_route_idx.push_back(r.idx);
        std::sort(unique_route_idx.begin(), unique_route_idx.end());
        unique_route_idx.erase(std::unique(unique_route_idx.begin(), unique_route_idx.end()),
                               unique_route_idx.end());

        if (opt.interlocking) init_interlock();

        maxitererr.assign(opt.iter_max + 2, 0.0);
        averagingstart = -1;
        hasconverged = false;
        moment_pass_done = false;
        results.clear();

        servtcdf.assign(N + 1, fluid::FluidPassage());
        callservtcdf.assign(lqn.ncalls + 1, fluid::FluidPassage());
        entrycdfrespt.assign(lqn.nentries + 1, LnCdf());
        entryproc.assign(lqn.nentries + 1, mam::AphPair<T>());
        cdf_repo.assign(ensemble.size(), std::vector<std::vector<fluid::FluidPassage>>());
    }

    /**
     * The Picard iteration, with the stochastic controller in place of the
     * deterministic test when the layer engine is a simulator.
     *
     * The two tests are NOT interchangeable and cannot both run: `converged`
     * folds a moving average into `results` in place, which is exactly what the
     * Polyak-Ruppert average would then be taken of a second time. Whichever
     * test is in force therefore owns the results.
     */
    void iterate() {
        init();
        const bool stoch = opt.layer_solver == "ssa";
        if (stoch) {
            LnStochConfig cfg;
            cfg.iter_tol = opt.iter_tol;
            cfg.relax_burnin = relax_omega;
            stoch_ctl = std::make_shared<LnStochController<T>>(cfg);
        }
        int it = 0;
        while (it < opt.iter_max) {
            if (!stoch && converged(it)) break;
            ++it;
            results.emplace_back(ensemble.size());
            for (std::size_t e = 0; e < ensemble.size(); ++e) analyze(it, e);
            post(it);
            if (stoch) {
                std::vector<double> jobs(ensemble.size(), 0.0);
                for (std::size_t e = 0; e < ensemble.size(); ++e)
                    jobs[e] = ensemble[e].total_jobs();
                const bool stop = stoch_ctl->update(it, results.back(), jobs, servt, residt);
                // The step this sets is the one the NEXT updateMetrics applies,
                // which is why it is read back here and not before the sweep.
                relax_omega = stoch_ctl->relax_omega();
                if (stop) {
                    did_converge = true;
                    hasconverged = true;
                    break;
                }
            }
        }
        iterations_done = it;
        // finish(): in Robbins-Monro mode report the averaged iterate, not the
        // last (noisy) one.
        if (stoch && stoch_ctl && stoch_ctl->averaging_count() > 0) {
            const std::vector<LayerResult<T>>& avg = stoch_ctl->averaged_results();
            for (std::size_t e = 0; e < results.back().size() && e < avg.size(); ++e)
                results.back()[e] = avg[e];
            servt = stoch_ctl->averaged_servt();
            residt = stoch_ctl->averaged_residt();
        }
    }

    /**
     * Port of the filterMetric helper inside @@NetworkSolver/getAvg.m.
     *
     * This is not cosmetic post-processing: it is where a raw solver matrix
     * becomes the result the rest of LINE consumes, and three of its rules
     * change numbers that SolverLN then iterates on.
     *
     *   1. a (station, class) pair the class never visits is zeroed, whether
     *      because its service is disabled or because its visit ratio is zero;
     *   2. anything below FineTol is snapped to zero;
     *   3. the caller's zeroMask is applied -- and for the queue length and the
     *      utilization that mask is `RN < 10*FineTol`, which zeroes both
     *      wherever the response time is immediate.
     *
     * Rule 3 is the one that matters most here. Every LQN entry, task and call
     * class is served by an Immediate distribution somewhere, whose response
     * time is 1e-8, so without it a layer reports the Immediate classes' share
     * of the population as real queue length and the LN think-time update reads
     * a task utilization that MATLAB reports as zero.
     */
    Matrix<T> filter_metric(const qn::Layer<T>& L, const Matrix<T>& metric,
                            const std::vector<std::vector<bool>>* zero_mask) const {
        const std::size_t M = L.nstations, K = L.nclasses;
        Matrix<T> out(M, K, Tzero());
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k)
                if (!L.disabled[i][k]) out(i, k) = metric(i, k);
        if (zero_mask)
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t k = 0; k < K; ++k)
                    if ((*zero_mask)[i][k]) out(i, k) = Tzero();
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k)
                if (dbl(out(i, k)) < GlobalConstants::FineTol) out(i, k) = Tzero();
        for (std::size_t k = 0; k < K; ++k) {
            std::size_t c = L.nchains;
            for (std::size_t cc = 0; cc < L.nchains; ++cc)
                if (L.chains[cc][k]) c = cc;
            if (c == L.nchains) continue;
            for (std::size_t i = 0; i < M; ++i)
                if (L.visits[c](L.stateful_of_station(i + 1) - 1, k) == Tzero())
                    out(i, k) = Tzero();
        }
        return out;
    }

    /**
     * Solve one layer, running the fork-join fixed point when it has a fork.
     *
     * The driver itself is shared (fj_driver.h); the layer path supplies
     * solver_mva_analyzer as the inner solve. The auxiliary arrival rates
     * fj_lambda are warm-started across outer iterations, per the reference.
     *
     * `default` becomes `amva` on a fork layer, which is what mvaDispatch.m does
     * for any model whose BASE has a fork: the transformed layer is mixed with
     * auxiliary near-zero-rate open classes, and the exact mixed MVA the default
     * ladder would otherwise pick degenerates on them.
     */
    /**
     * Run one layer through the fluid analyzer instead of MVA.
     *
     * WHAT MAKES THIS SOUND AT ALL. The layer handed to a solver is an ordinary
     * closed queueing network -- the layering has already replaced every call
     * by a class with a service demand -- so any NetworkSolver can solve it,
     * and the reference says exactly that by taking a solver factory. The outer
     * Picard iteration only ever reads [QN, UN, RN, TN] back, so the fluid
     * result maps onto the same shape MVA returns.
     *
     * WHAT CHANGES, and it is not nothing: the fluid limit is exact only as the
     * populations grow, so on a layer of a few jobs it is a genuinely different
     * approximation, not a slower route to the same fixed point. It also feeds
     * DIFFERENT service demands back into the next outer iteration, so the two
     * ensembles converge to different fixed points rather than to the same one
     * by different paths.
     *
     * REFUSED BY NAME: a non-double backend (LSODA is double-only, see
     * solver_fluid.h) and a fork layer (`fj_fixed_point` supplies MVA as its
     * inner solve, and the auxiliary near-zero-rate open classes it introduces
     * are not something the fluid drift represents).
     */
    /**
     * Solve one layer with CTMC and report it in the shape the MVA path returns.
     *
     * Only reached for a layer carrying a Region. `qn::Layer<T>` derives from
     * NetworkStruct, so the layer goes to the analyzer as it stands; the
     * generator needs the routing matrix, which the layer does not build until
     * asked (see solve_layer_fluid).
     */
    mva::MvaSolution<T> solve_layer_ctmc(std::size_t e) {
        qn::Layer<T>& L = ensemble[e];
        L.refresh_rt();
        // The capacities the WAITQ generator reads to bound each station's
        // marginal are filled by solve_layer, for every engine; the routing
        // matrix is not, because only this path and the fluid one need it.
        ctmc::CtmcOptions co;
        // _run_any, not _run: the plain generator implements DROP only, and a
        // region built from an LQN admission constraint carries the WAITQ
        // default, whose per-region token FIFO lives in the waitq generator.
        const mva::AvgResult<T> a = ctmc::solver_ctmc_run_analyzer_any(L, co);
        mva::MvaSolution<T> out;
        out.method = "ctmc";
        out.iter = 1;
        out.Q = a.QN;
        out.U = a.UN;
        out.R = a.RN;
        out.Tp = a.TN;
        out.C = a.CN;
        out.X = a.XN;
        return out;
    }

    /**
     * Solve one layer with SolverMAM's `dec.poisson`, for a setup/delay-off.
     *
     * The method matters: `dec.poisson` caps the arrival superposition at one
     * phase AND skips the fork-join / retrial / ldqbd diversions, so the station
     * reaches the setup branch of solver_mam_basic rather than being taken by a
     * shape-matching special case. It is what SolverLN.m:183 asks for.
     */
    mva::MvaSolution<T> solve_layer_mam(std::size_t e) {
        qn::Layer<T>& L = ensemble[e];
        L.refresh_rt();
        L.refresh_capacity();
        mam::MamOptions mo;
        mo.method = "dec.poisson";
        mva::MvaSolution<T> out = mam::mam_dispatch(L, mo).sol;
        out.method = "mam";
        return out;
    }

    /**
     * Solve one layer with SolverNC, the reference's `LN(model, @(m) NC(m))`.
     *
     * This is the layer engine the Java CLI's bare `ln` token selects, against
     * `ln.mva` for the MVA one. NC evaluates the normalizing constant rather
     * than the MVA recursion, so on a product-form layer it is the SAME answer
     * reached exactly instead of through the AMVA approximation -- and on a
     * layer that is not product form the two differ, which is the whole reason
     * the token exists.
     */
    mva::MvaSolution<T> solve_layer_nc(std::size_t e) {
        if (fj_tr[e].active())
            throw UnsupportedError(
                "SolverLN: layer '" + ensemble[e].name +
                "' carries a fork, whose fixed point is driven by MVA in this port; solve this "
                "model with layer_solver 'mva'");
        qn::Layer<T>& L = ensemble[e];
        L.refresh_rt();
        // `solver_nc_solve`, NOT `nc_dispatch`: the latter is the INNER network
        // solve, the one the caching-queueing decomposition itself calls as its
        // `netfun`, and it carries no cache branch at all. The three cache gates
        // live in `solver_nc_solve` (solver_nc_runner.h), which is the true
        // counterpart of `mva_dispatch` on the branch below. Routed here, a
        // cache layer never reached `solver_nc_cacheqn_analyzer` and its Cache
        // node was solved as ordinary routing, so the hit/miss split came back
        // as the probability `link()` offered -- an even 1/2, independent of
        // capacity, item count and replacement strategy, with no warning.
        // `.sol` alone is taken, exactly as the MVA branch takes it: the layer
        // reads the split off its own station throughputs, and the auxiliary
        // `refreshed_struct` both dispatchers also return is unused on either.
        mva::MvaSolution<T> out = nc::solver_nc_solve(L, opt.layer_nc).sol;
        out.method = "nc";
        return out;
    }

    /**
     * Solve one layer by simulation, the reference's `LN(model, @(m) SSA(m))`.
     *
     * THE LAYER RESULTS ARE THEN NOISY, and the deterministic convergence test
     * cannot terminate against noise: its successive-difference error is
     * bounded below by the standard error of the estimates. `iterate` therefore
     * switches to `LnStochController` whenever this engine is selected, so the
     * relaxation decays as Robbins-Monro and the reported iterate is the
     * Polyak-Ruppert average. Selecting this engine and keeping the
     * deterministic test would simply run to iter_max.
     */
    mva::MvaSolution<T> solve_layer_ssa(std::size_t e) {
        if (fj_tr[e].active())
            throw UnsupportedError(
                "SolverLN: layer '" + ensemble[e].name +
                "' carries a fork, whose fixed point is driven by MVA in this port; solve this "
                "model with layer_solver 'mva'");
        qn::Layer<T>& L = ensemble[e];
        L.refresh_rt();
        const ssa::SsaSolution s = ssa::solver_ssa(L, opt.layer_ssa);
        mva::MvaSolution<T> out;
        out.method = "ssa";
        out.iter = 1;
        out.Q = Matrix<T>(L.nstations, L.nclasses, Tzero());
        out.U = Matrix<T>(L.nstations, L.nclasses, Tzero());
        out.R = Matrix<T>(L.nstations, L.nclasses, Tzero());
        out.Tp = Matrix<T>(L.nstations, L.nclasses, Tzero());
        for (std::size_t i = 0; i < L.nstations; ++i)
            for (std::size_t r = 0; r < L.nclasses; ++r) {
                out.Q(i, r) = num_traits<T>::from_double(s.QN(i, r));
                out.U(i, r) = num_traits<T>::from_double(s.UN(i, r));
                out.R(i, r) = num_traits<T>::from_double(s.RN(i, r));
                out.Tp(i, r) = num_traits<T>::from_double(s.TN(i, r));
            }
        out.C.assign(L.nclasses, Tzero());
        out.X.assign(L.nclasses, Tzero());
        for (std::size_t r = 0; r < L.nclasses && r < s.XN.size(); ++r) {
            out.X[r] = num_traits<T>::from_double(s.XN[r]);
            out.C[r] = num_traits<T>::from_double(s.CN[r]);
        }
        return out;
    }

    mva::MvaSolution<T> solve_layer_fluid(std::size_t e) {
        if (fj_tr[e].active())
            throw UnsupportedError(
                "SolverLN: layer '" + std::to_string(e) +
                "' carries a fork, whose fixed point is driven by MVA in this port; solve this "
                "model with layer_solver 'mva'");
        qn::Layer<T>& L = ensemble[e];
        // THE LAYER HAS NO ROUTING MATRIX UNTIL IT IS ASKED FOR ONE. `buildLayers`
        // records the wiring with set_route and then calls refresh_chains, which
        // reads P directly and computes the VISITS; MVA needs nothing else, so
        // `rt` stays empty. Both fluid methods route through `sn.rt` instead, and
        // an empty one is not an error anywhere -- it is read as "no transition
        // exists", the drift decays to zero and every metric comes back 0. It is
        // rebuilt on every solve because update_routing_probabilities can change
        // the wiring between outer iterations.
        L.refresh_rt();
        mva::MvaSolution<T> out;
        out.method = "fluid";
        out.iter = 1;
        out.Q = Matrix<T>(L.nstations, L.nclasses, Tzero());
        out.U = Matrix<T>(L.nstations, L.nclasses, Tzero());
        out.R = Matrix<T>(L.nstations, L.nclasses, Tzero());
        out.Tp = Matrix<T>(L.nstations, L.nclasses, Tzero());
        out.C.assign(L.nclasses, Tzero());
        out.X.assign(L.nclasses, Tzero());
        detail::ln_fluid_solve(L, opt.layer_fluid, out);
        return out;
    }

    mva::MvaSolution<T> solve_layer(std::size_t e) {
        // BEFORE ANY ENGINE, AND FOR EVERY ONE OF THEM. `build_layer` stops at
        // refresh_chains, so a layer reaches this point with `cap` and
        // `classcap` still EMPTY, while the reference hands SolverMVA a struct
        // that refreshCapacity has already filled. Every consumer of Kendall's
        // K -- `buffer_size`, and through it `has_blocking` and the BCMP gate
        // `has_product_form` -- indexes those vectors by station, so an empty
        // one is an out-of-range read and not a permissive default. It is
        // recomputed on each solve because the layer's class populations move
        // between outer iterations, and the buffers are derived from them.
        ensemble[e].refresh_capacity();
        // A layer carrying an admission constraint goes to CTMC whatever the
        // layer solver is, as adaptiveSolverFactory does (SolverLN.m:1047-1089):
        // MVA and the fluid path both refuse a Region by name, so leaving the
        // layer on them would refuse a model the reference solves. The
        // reference then falls back to LDES and SSA; neither is reachable here
        // (there is no cpp LDES, and cpp SSA refuses regions), so CTMC is the
        // only path and its state-space bound is this port's real limit.
        if (!ensemble[e].regions.empty()) return solve_layer_ctmc(e);
        // A setup no longer forces the MAM decomposition on the layer
        // (SolverLN.m, 2026-08-11): the cold start is charged to the entry by
        // setup_charge(), not wired into the station, so the layer is an
        // ordinary one and the user's own layer solver serves it.
        // A cache layer is dispatched to the integrated caching-queueing
        // analyzer, which reads the NODE routing to find what lies downstream
        // of the cache. The layer carries none until asked, and it has to be
        // rebuilt every pass because entry selection can move it. EVERY layer
        // engine below must reach that analyzer, not just the MVA one: this
        // refresh runs for all of them, and a branch that then hands the
        // rebuilt routing to a solver with no cache gate reports the offered
        // hit/miss split instead of the converged one.
        if (has_cache_node(e)) ensemble[e].refresh_rt();
        if (opt.layer_solver == "fluid") return solve_layer_fluid(e);
        if (opt.layer_solver == "nc") return solve_layer_nc(e);
        if (opt.layer_solver == "ssa") return solve_layer_ssa(e);
        if (opt.layer_solver != "mva")
            throw UnsupportedError("SolverLN: layer_solver '" + opt.layer_solver +
                                   "' is not available; use 'mva', 'nc', 'fluid' or 'ssa'");
        // Route each layer through the full mva_dispatch ladder, as MATLAB's
        // SolverLN does by handing every layer to SolverMVA.runAnalyzer: a plain
        // closed layer still lands on branch 13 (solver_mva_analyzer), but a
        // layer carrying load- or class-dependent scaling now reaches
        // solver_mvald_analyzer, matching the reference rather than silently
        // running the flat AMVA.
        mva::MvaOptions eopt = opt.layer;
        if (e < layer_interlock.size()) eopt.interlock = layer_interlock[e];
        if (!fj_tr[e].active())
            return mva::mva_dispatch(ensemble[e], eopt, layer_init_sol[e]).sol;
        mva::MvaOptions lopt = eopt;
        if (lopt.method == "default") lopt.method = "amva";
        return mva::fj_fixed_point(
            ensemble[e], fj_tr[e], fj_lambda[e], eopt,
            [&lopt](qn::NetworkStruct<T>& V) {
                return mva::mva_dispatch(V, lopt, Matrix<T>()).sol;
            });
    }

    void analyze(int it, std::size_t e) {
        qn::Layer<T>& L = ensemble[e];
        const mva::MvaSolution<T> s = solve_layer(e);
        LayerResult<T>& r = results[it - 1][e];
        r.RN = filter_metric(L, s.R, nullptr);
        std::vector<std::vector<bool>> zmask(L.nstations, std::vector<bool>(L.nclasses, false));
        for (std::size_t i = 0; i < L.nstations; ++i)
            for (std::size_t k = 0; k < L.nclasses; ++k)
                zmask[i][k] = dbl(r.RN(i, k)) < 10.0 * GlobalConstants::FineTol;
        r.QN = filter_metric(L, s.Q, &zmask);
        r.UN = filter_metric(L, s.U, &zmask);
        r.TN = filter_metric(L, s.Tp, nullptr);
        r.WN = residence_from_response(L, r.RN);
        // warm start the next solve of this layer from the chain-aggregated
        // queue lengths, as the reference does for SolverMVA layers
        // A fork layer is exempt: the reference guards this on the result having
        // the LAYER's station count, and the transformed model reports one more
        // row (its Source), so the hint is never installed there. A fluid layer
        // is exempt too, and for the reference's own reason: `analyze` guards
        // the warm start on `strcmp(self.solvers{e}.name, 'SolverMVA')`, because
        // init_sol means a queue-length vector to AMVA and a full ODE state
        // vector to the fluid solver -- the two are not interchangeable.
        if (!fj_tr[e].active() && opt.layer_solver == "mva") {
            Matrix<T> Qch(L.nstations, L.nchains, Tzero());
            for (std::size_t c = 0; c < L.nchains; ++c)
                for (std::size_t i = 0; i < L.nstations; ++i) {
                    T s2 = Tzero();
                    for (std::size_t k : L.inchain[c]) s2 += s.Q(i, k - 1);
                    Qch(i, c) = s2;
                }
            layer_init_sol[e] = Qch;
        }
    }

    /** Port of sn_get_residt_from_respt: response time scaled by the visit ratio. */
    Matrix<T> residence_from_response(const qn::Layer<T>& L, const Matrix<T>& RN) const {
        Matrix<T> V(L.nstations, L.nclasses, Tzero());
        for (std::size_t c = 0; c < L.nchains; ++c)
            for (std::size_t i = 0; i < L.nstations; ++i) {
                const std::size_t sf = L.stateful_of_station(i + 1) - 1;
                for (std::size_t k = 0; k < L.nclasses; ++k)
                    V(i, k) = T(V(i, k) + L.visits[c](sf, k));
            }
        Matrix<T> WN(L.nstations, L.nclasses, Tzero());
        for (std::size_t i = 0; i < L.nstations; ++i)
            for (std::size_t k = 0; k < L.nclasses; ++k) {
                if (L.disabled[i][k]) continue;
                if (!(RN(i, k) > Tzero())) continue;
                if (dbl(RN(i, k)) < GlobalConstants::FineTol) {
                    WN(i, k) = RN(i, k);
                    continue;
                }
                std::size_t c = 0;
                for (std::size_t cc = 0; cc < L.nchains; ++cc)
                    if (L.chains[cc][k]) c = cc;
                const std::size_t rstat = L.classes[k].refstat;
                T den = Tzero();
                if (L.refclass[c] > 0) {
                    den = V(rstat - 1, L.refclass[c] - 1);
                } else {
                    for (std::size_t kk : L.inchain[c]) den += V(rstat - 1, kk - 1);
                }
                if (den == Tzero()) continue;
                WN(i, k) = T(RN(i, k) * V(i, k) / den);
            }
        for (std::size_t i = 0; i < L.nstations; ++i)
            for (std::size_t k = 0; k < L.nclasses; ++k)
                if (dbl(WN(i, k)) < 10.0 * GlobalConstants::FineTol) WN(i, k) = Tzero();
        return WN;
    }

    void post(int it) {
        update_metrics(it);
        if (opt.interlocking) update_populations(it);
        update_think_times(it);
        update_layers(it);
        update_routing_probabilities(it);
        for (std::size_t e : route_reset) {
            ensemble[e].refresh_chains();
            layer_init_sol[e] = Matrix<T>();
        }
        // moment3 needs no refreshProcesses here: this port keeps the service
        // DISTRIBUTION on the layer and refresh_rates re-reads it, so there is no
        // second representation to fall out of step, unlike sn.proc in MATLAB.
        for (std::size_t e : svc_reset) ensemble[e].refresh_rates();
    }

    /** Port of converged.m: moving average of the layer results plus the test. */
    bool converged(int it) {
        const std::size_t E = ensemble.size();
        const int iter_min = std::max<int>(2 * int(E), int(std::ceil(opt.iter_max / 4.0)));
        const int wnd_size = std::max(5, int(std::ceil(iter_min / 5.0)));

        if (it >= iter_min && int(results.size()) >= wnd_size) {
            const T w = T(Tone() / num_traits<T>::from_int(wnd_size));
            for (std::size_t e = 0; e < E; ++e) {
                LayerResult<T>& cur = results[results.size() - 1][e];
                auto scale = [&](Matrix<T>& A) {
                    for (std::size_t i = 0; i < A.rows(); ++i)
                        for (std::size_t j = 0; j < A.cols(); ++j) A(i, j) = T(A(i, j) * w);
                };
                Matrix<T> Q = cur.QN, U = cur.UN, R = cur.RN, Tp = cur.TN, W = cur.WN;
                scale(Q); scale(U); scale(R); scale(Tp); scale(W);
                for (int k = 1; k < wnd_size; ++k) {
                    const LayerResult<T>& old = results[results.size() - 1 - k][e];
                    auto add = [&](Matrix<T>& A, const Matrix<T>& B) {
                        for (std::size_t i = 0; i < A.rows(); ++i)
                            for (std::size_t j = 0; j < A.cols(); ++j) A(i, j) = T(A(i, j) + B(i, j) * w);
                    };
                    add(Q, old.QN); add(U, old.UN); add(R, old.RN); add(Tp, old.TN); add(W, old.WN);
                }
                cur.QN = Q; cur.UN = U; cur.RN = R; cur.TN = Tp; cur.WN = W;
            }
        }

        if (it > 1) {
            double err = 0.0;
            for (std::size_t e = 0; e < E; ++e) {
                const Matrix<T>& Q = results[results.size() - 1][e].QN;
                const Matrix<T>& Q1 = results[results.size() - 2][e].QN;
                const double Njobs = ensemble[e].total_jobs();
                if (!(Njobs > 0.0)) continue;
                double mx = 0.0;
                for (std::size_t i = 0; i < Q.rows(); ++i)
                    for (std::size_t j = 0; j < Q.cols(); ++j)
                        mx = std::max(mx, std::fabs(dbl(Q(i, j)) - dbl(Q1(i, j))));
                err += mx / Njobs;
            }
            maxitererr[it] = err;
            line::util::LineConsole::iter(
                static_cast<long>(it),
                "layer iteration %zu: max queue-length change %.3e (tolerance %.3e)",
                static_cast<std::size_t>(it), err, opt.iter_tol);
            if (it == iter_min) {
                line::util::LineConsole::step("started averaging the iterates to aid convergence");
                averagingstart = it;
            }
        }

        if (it > iter_min && maxitererr[it] < opt.iter_tol && maxitererr[it - 1] < opt.iter_tol &&
            maxitererr[it - 2] < opt.iter_tol) {
            if (!hasconverged) {
                hasconverged = true;
            } else {
                did_converge = true;
                return true;
            }
        } else {
            hasconverged = false;
        }
        return false;
    }

    // -----------------------------------------------------------------------
    // updateMetrics
    // -----------------------------------------------------------------------
    /** Port of updateMetrics.m: the method selects which update runs. */
    // -----------------------------------------------------------------------
    // Method "srvn.ph": the activity graph of an entry as a phase-type server law
    //
    // Each layer is a two-station cycle, a client Delay plus the server, with
    // one closed class per caller task. The sequencing the default method
    // encodes as routing -- a class per entry, per activity and per call, plus
    // Fork, Join, Router and ClassSwitch nodes -- is composed instead into a
    // single phase-type service law per (layer, caller), by the exact
    // series-parallel reduction of Workflow. The layer therefore carries only
    // the client/server back-and-forth, and the activity graph survives as a
    // distribution.
    //
    // Port of the MATLAB @SolverLN/buildLayersPH.m, updateLayersPH.m,
    // updateMetricsPH.m, updateThinkTimesPH.m, phComposeEntryLaws.m and
    // getEnsembleAvgPH.m, of the JAR jline.solvers.ln.SolverLNPH and of the
    // Python line_solver.solvers.solver_ln.solver_ln_ph. See
    // _kb/06-solver-catalog.md (LN section) for the layering taxonomy.
    // -----------------------------------------------------------------------

    /**
     * The method the layers were actually BUILT for: "srvn.ph", "srvn.cs",
     * "flat.cs", "flat.ph" or "moment3". Resolved once in build_layers, because
     * the alias "srvn" may fall back; every dispatch reads this and not
     * opt.method, so a reconstruction can never disagree with the layers it is
     * reading.
     */
    std::string lnmethod;
    /** True once ph_init_laws has composed the per-entry workflows. */
    bool ph_laws_ready = false;

    /**
     * Normalise a method name onto one the solver dispatches on. A method name
     * carries TWO decisions: the LAYERING, which fixes what a submodel is, and
     * the ENCODING, which fixes how an activity graph is written into it.
     *
     * "srvn.cs" encodes the activity graph as ROUTING, "srvn.ph" as a composed
     * phase-type server law, "srvn" is the alias that takes "srvn.ph" where it
     * can serve the model and "srvn.cs" otherwise, "flat.cs" squashes every
     * server into one submodel with the routing encoding, "flat.ph" squashes
     * them with the composed one ("flat" is the alias of "flat.cs" and resolves
     * unconditionally rather than probing "flat.ph", because a model is squashed
     * in order to express what only the routing encoding carries), and "moment3"
     * is the three-moment distribution pass over the routing layers. "default"
     * is the srvn alias; an unrecognised method name takes "srvn.cs".
     */
    static std::string ln_requested_method(const std::string& method) {
        std::string m;
        for (char c : method) m += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (m.empty() || m == "srvn" || m == "default" || m == "auto") return "srvn";
        if (m == "srvn.ph" || m == "ph") return "srvn.ph";
        if (m == "srvn.cs" || m == "srvncs" || m == "cs") return "srvn.cs";
        if (m == "flat.cs" || m == "flatcs" || m == "flat" || m == "squashed") return "flat.cs";
        if (m == "flat.ph" || m == "flatph" || m == "squashed.ph") return "flat.ph";
        if (m == "moment3") return "moment3";
        // An unrecognised token takes the routing encoding, which is what every
        // name other than "moment3" resolved to before the alias existed.
        return "srvn.cs";
    }

    /** True when the layers are the collapsed phase-type ones of "srvn.ph". */
    bool is_srvn_ph() const { return lnmethod == "srvn.ph"; }

    /**
     * True when the layers carry the COMPOSED phase-type server law rather than
     * the routing encoding of the activity graph, under either layering. The
     * encoding, not the layering, decides which update and reconstruction passes
     * run, so every such dispatch asks this and not for one method name.
     */
    bool is_ph_encoding() const { return lnmethod == "srvn.ph" || lnmethod == "flat.ph"; }

    /**
     * Station of layer L that stands for LQN element ELEM, falling back to the
     * layer's own server when ELEM is not a server there. Under "srvn" the
     * fallback is the answer for every element; under "flat.cs" it is the map
     * that tells the many servers of one layer apart.
     */
    static std::size_t station_idx_of(const qn::Layer<T>& L, std::size_t elem) {
        if (elem >= 1 && elem < L.server_idx_of.size() && L.server_idx_of[elem] != 0)
            return L.server_idx_of[elem];
        return L.serverIdx;
    }

    /**
     * Station of layer L that class K (0-based) is served at: the PROCESSOR of
     * an activity, the CALLED TASK of a call, the layer's own server otherwise.
     */
    std::size_t station_idx_of_class(const qn::Layer<T>& L, std::size_t k) const {
        const int kind = L.classes[k].attr_kind;
        const std::size_t a = L.classes[k].attr_idx;
        if (kind == int(LqnElement::ACTIVITY))
            return station_idx_of(L, lqn.parent[lqn.parent[a]]);
        if (kind == int(LqnElement::CALL))
            return station_idx_of(L, lqn.parent[lqn.callpair_dst[a]]);
        return L.serverIdx;
    }

    /**
     * Answer whether "srvn.ph" can serve this model, without disturbing the
     * solver.
     *
     * Both the feature gate and the series-parallel reduction can refuse, and
     * the second only finds out by composing the per-entry workflows -- work the
     * build then reuses, since those laws do not depend on the iterate.
     */
    bool probe_srvn_ph() {
        try {
            ph_assert_supported();
            ph_init_laws();
            ph_laws_ready = true;
            return true;
        } catch (const std::exception&) {
            ph_laws_ready = false;
            return false;
        }
    }

    /** One two-station layer, and the caller classes that cycle through it. */
    struct PHLayer {
        std::size_t idx = 0;
        bool ishost = false;
        std::vector<std::size_t> callers;
        /** 1-based class index of each caller task, 0 when absent. */
        std::vector<std::size_t> class_of_caller;
        std::size_t nreplicas = 1;
        /** 1-based station indices of the server replicas. */
        std::vector<std::size_t> qstations;
        /** Mean of the law each class is currently served with, by class index. */
        std::vector<T> svcmean_by_class;
        /** (class index, entry index) or (class index, -call index). */
        std::vector<std::pair<std::size_t, long>> open_arrivals;
        /**
         * Closed population of the MODEL this server sits in. Under "flat.ph"
         * that is every caller of the single network, not only the callers of
         * this one station, so it is recorded here rather than recomputed.
         */
        double npop = 0.0;
    };

    // per-entry workflows and their composed laws
    std::vector<workflow::Workflow<T>> ph_wf, ph_wfhost;
    std::vector<std::unordered_map<std::size_t, T>> ph_execs;
    std::vector<bool> ph_has_wf;
    std::vector<workflow::PhLaw<T>> ph_hostlaw, ph_entrylaw;
    std::vector<T> ph_hostmean, ph_entrymean, ph_entryscv;
    std::vector<T> ph_share, ph_overlap, ph_setupshare, ph_xdemand;
    Matrix<T> ph_ncalls, ph_calltime;
    std::vector<T> ph_procresid, ph_actthinkt, ph_calltotal;
    std::vector<PHLayer> ph_layers;   ///< by element index, empty where absent
    std::vector<bool> ph_has_layer;

    /** Mean of an activity's own think time, 0 when it declares none. */
    T act_think_time(std::size_t aidx) const {
        if (aidx >= lqn.actthink.size() || lqn.actthink[aidx].disabled) return Tzero();
        const double m = dbl(lqn.actthink[aidx].mean);
        if (!std::isfinite(m) || m <= GlobalConstants::FineTol) return Tzero();
        return lqn.actthink[aidx].mean;
    }

    /**
     * Mean cold start one request of task TIDX pays, 0 when it declares none.
     *
     * A SetupTask powers a thread down when it goes idle and pays a setup before
     * it can serve again. The thread is released at a reply and starts a
     * delay-off countdown D of mean d; it powers off only if D expires before
     * the next request arrives, and a request arriving first cancels the
     * countdown and pays nothing. With the idle interval I seen by one thread
     * and exponential D, p = P(D < I) = E[I]/(E[I]+d) and the charge is p*s.
     *
     * Admission takes an ACTIVE idle thread before it wakes a sleeping one, so
     * the pool that actually cycles is only as large as the load needs: with
     * offered load b = X*S = rho*mult threads, about max(1,b) stay hot, giving
     * E[I] = (max(1,b) - b)/X. Exact at mult = 1; above it the exact answer is
     * matrix-analytic (Gandhi, Harchol-Balter and Adan, Performance Evaluation
     * 67(11), 2010). Twin of MATLAB lqn_setup_charge.m, the JAR
     * SolverLN.setupCharge and the Python SolverLN._setup_charge.
     */
    double setup_charge(std::size_t tidx) const {
        if (tidx >= lqn.hassetup.size() || !lqn.hassetup[tidx]) return 0.0;
        const double s = tidx < lqn.setuptime.size() && !lqn.setuptime[tidx].disabled
                             ? dbl(lqn.setuptime[tidx].mean) : 0.0;
        const double d = tidx < lqn.delayofftime.size() && !lqn.delayofftime[tidx].disabled
                             ? dbl(lqn.delayofftime[tidx].mean) : 0.0;
        if (!(s > GlobalConstants::FineTol) || !(d > GlobalConstants::FineTol)) return 0.0;
        const double mult = lqn.mult[tidx];
        if (!std::isfinite(mult) || mult <= 0.0) return 0.0;
        if (tidx >= tput.size() || tidx >= util.size()) return s;
        const double X = dbl(tput[tidx]);
        if (!std::isfinite(X) || X <= GlobalConstants::FineTol) return s;
        double rho = dbl(util[tidx]);
        if (!std::isfinite(rho) || rho < 0.0) rho = 0.0;
        rho = std::min(rho, 1.0 - GlobalConstants::FineTol);
        const double b = rho * mult;                     // offered load, in threads
        const double EI = (std::max(1.0, b) - b) / X;    // idle interval of a hot thread
        return s * EI / (EI + d);
    }

    /**
     * Probability that a request for entry EIDX finds its task's thread off.
     * ONE closure for both methods: setup_charge returns p*s, so p is that over
     * s. It also answers p = 1 during construction, before the first solve has
     * sized tput or util.
     */
    double ph_setup_prob(std::size_t eidx) const {
        const std::size_t tidx = lqn.parent[eidx];
        if (tidx >= lqn.hassetup.size() || !lqn.hassetup[tidx]) return 0.0;
        const double s = tidx < lqn.setuptime.size() && !lqn.setuptime[tidx].disabled
                             ? dbl(lqn.setuptime[tidx].mean) : 0.0;
        const double d = tidx < lqn.delayofftime.size() && !lqn.delayofftime[tidx].disabled
                             ? dbl(lqn.delayofftime[tidx].mean) : 0.0;
        if (!(s > GlobalConstants::FineTol) || !(d > GlobalConstants::FineTol)) return 0.0;
        return std::min(1.0, std::max(0.0, setup_charge(tidx) / s));
    }

    /** Divisor that scales a processor utilization into [0,1]. */
    double ph_host_servers(std::size_t hidx) const {
        if (lqn.sched[hidx] == SchedStrategy::INF) return 1.0;
        const double m = lqn.maxmult[hidx];
        return (std::isfinite(m) && m > 0.0) ? m : 1.0;
    }

    /** True when any task calls entry EIDX, synchronously or not. */
    bool ph_any_caller_of(std::size_t eidx) const {
        return lqn.issynccaller.any_col(eidx) || lqn.isasynccaller.any_col(eidx);
    }

    /**
     * True when an entry arrival is the ONLY way requests reach task TIDX.
     * "srvn.ph" refuses forwarding calls outright, so sync/async callers are the
     * whole test.
     */
    bool ph_open_arrival_only(std::size_t tidx) const {
        if (lqn.isref[tidx]) return false;
        for (std::size_t eidx : lqn.entriesof[tidx])
            if (ph_any_caller_of(eidx)) return false;
        for (std::size_t eidx : lqn.entriesof[tidx])
            if (lqn.has_arrival[eidx]) return true;
        return false;
    }

    /** Asynchronous calls whose target entry belongs to TIDX. */
    std::vector<std::size_t> ph_async_calls_into(std::size_t tidx) const {
        std::vector<std::size_t> out;
        for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
            if (lqn.calltype[cidx] != CallType::ASYNC) continue;
            for (std::size_t e : lqn.entriesof[tidx])
                if (lqn.callpair_dst[cidx] == e) { out.push_back(cidx); break; }
        }
        return out;
    }

    /**
     * Features the collapsed layer cannot represent are refused by name rather
     * than silently degraded.
     */
    void ph_assert_supported(bool flat = false) const {
        // The list is a property of the ENCODING, so it is the same under either
        // layering; what the squashing adds on top is refused in ph_flat_server_set.
        const std::string mname = flat ? "flat.ph" : "srvn.ph";
        // PHASE 2 IS ASKED HERE AND NOWHERE ELSE. `has_phase2` is built during
        // layering rather than being a property of the model, so the report twin
        // below cannot ask it; every other rule is shared with it.
        if (has_phase2)
            throw UnsupportedError(
                "method='" + mname + "' does not support second-phase activities: the composed "
                "entry law has no reply point. Use method='default'.");
        // ONE RULE LIST, TWO CALLERS: this throws what `ph_method_refusal`
        // returns, in the same order, so the run and a report cannot say
        // different things about one model. The squashing rules it also carries
        // for `flat.ph` are the ones `ph_flat_server_set` raises, which runs
        // before this on that path, so the message a caller sees is unchanged.
        const std::string why = ph_method_refusal(mname);
        if (!why.empty()) throw UnsupportedError(why);
    }

    /**
     * The same rules as a SENTENCE, so a report can withdraw a method it cannot run.
     *
     * ONE RULE LIST, TWO CALLERS. `ph_assert_supported` above answers the run
     * path by throwing; nothing answered the report, and `list_valid_methods`
     * returns the same eight names for every model, so every layered model was
     * offered every encoding and `srvn.ph`/`flat.ph` then threw on contact.
     *
     * PHASE 2 IS DELIBERATELY ABSENT. `has_phase2` is built during layering
     * rather than being a property of the model, so a gate cannot ask it without
     * doing the layering it is meant to precede; the run path still refuses it.
     *
     * @param method the concrete method name
     * @return empty string when the method can encode this model, else the reason
     */
    std::string ph_method_refusal(const std::string& method) const {
        if (method != "srvn.ph" && method != "flat.ph") return std::string();
        // The squashing refusals, `flat.ph` only: each carries PER-LAYER state
        // that one submodel cannot hold (ph_flat_server_set).
        if (method == "flat.ph") {
            for (std::size_t i = 1; i <= NT(); ++i) {
                if (lqn.repl[i] > 1.0)
                    return "method='flat.ph' does not support replicated processors or tasks, "
                           "whose replicas need a submodel each. Use method='srvn.ph'.";
                if (lqn.hassetup[i])
                    return "method='flat.ph' does not support setup tasks, whose powered-down "
                           "threads are per-layer state. Use method='srvn.ph'.";
            }
        }
        for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx)
            if (lqn.calltype[cidx] == CallType::FWD)
                return "method='" + method +
                       "' does not support forwarding calls, whose target is not part of the "
                       "caller's activity graph. Use method='default'.";
        for (std::size_t i = 1; i < lqn.iscache.size(); ++i)
            if (lqn.iscache[i])
                return "method='" + method +
                       "' does not support cache tasks. Use method='default'.";
        for (std::size_t i = 1; i < lqn.hassetup.size(); ++i) {
            if (!lqn.hassetup[i]) continue;
            if (lqn.sched[i] == SchedStrategy::INF || !std::isfinite(lqn.mult[i]))
                return "method='" + method + "': task '" + lqn.names[i] +
                       "' declares a setup time on an infinite-server task, which holds no "
                       "thread to power down; give it a finite multiplicity.";
        }
        for (std::size_t i = 0; i < lqn.lincon_A.size(); ++i)
            if (lqn.lincon_A[i].rows() > 0)
                return "method='" + method +
                       "' does not support admission constraints on a layer station. "
                       "Use method='default'.";
        for (std::size_t i = 1; i < lqn.lldscaling.size(); ++i) {
            const char* fname = nullptr;
            if (!lqn.lldscaling[i].empty()) fname = "a load dependence";
            else if (lqn.cdscaling[i]) fname = "a class dependence";
            else if (lqn.jdscaling[i]) fname = "a joint dependence";
            else if (!lqn.pools[i].empty()) fname = "server pools";
            if (fname != nullptr)
                return "method='" + method +
                       "' does not support queue-dependent service rates on a layer station ('" +
                       lqn.names[i] + "' declares " + fname + "). Use method='srvn.cs'.";
        }
        return std::string();
    }

    /** Build the per-entry workflows and the iteration-invariant processor law. */
    void ph_init_laws() {
        const std::size_t N = lqn.nidx;
        ph_wf.assign(N + 1, workflow::Workflow<T>("empty"));
        ph_wfhost.assign(N + 1, workflow::Workflow<T>("empty"));
        ph_execs.assign(N + 1, {});
        ph_has_wf.assign(N + 1, false);
        ph_hostlaw.assign(N + 1, workflow::PhLaw<T>());
        ph_entrylaw.assign(N + 1, workflow::PhLaw<T>());
        ph_hostmean.assign(N + 1, Tzero());
        ph_entrymean.assign(N + 1, Tzero());
        ph_entryscv.assign(N + 1, Tone());
        ph_share.assign(N + 1, Tzero());
        ph_overlap.assign(N + 1, Tone());
        ph_setupshare.assign(N + 1, Tzero());
        ph_xdemand.assign(N + 1, Tzero());
        ph_ncalls = Matrix<T>(N + 1, N + 1, Tzero());
        ph_calltime = Matrix<T>(N + 1, N + 1, Tzero());
        ph_procresid.assign(N + 1, Tzero());
        ph_actthinkt.assign(N + 1, Tzero());
        ph_calltotal.assign(N + 1, Tzero());
        ph_layers.assign(NT() + 1, PHLayer());
        ph_has_layer.assign(NT() + 1, false);

        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            const std::size_t tidx = lqn.parent[eidx];
            if (ignore[tidx]) continue;
            api::lqn::EntryWorkflow<T> ew = api::lqn::entry_workflow(lqn, eidx, true);
            ph_wf[eidx] = std::move(ew.wf);
            ph_execs[eidx] = ew.execs;
            ph_has_wf[eidx] = true;
            api::lqn::EntryWorkflow<T> eh = api::lqn::entry_workflow(lqn, eidx, false);
            ph_wfhost[eidx] = std::move(eh.wf);
            // the processor sees the WORK of concurrent branches, not their
            // elapsed time, so the host law serialises an AND fork
            ph_hostlaw[eidx] = api::lqn::serial_law(ph_wfhost[eidx]);
            ph_hostmean[eidx] =
                api::lqn::ph_moments(ph_hostlaw[eidx].alpha, ph_hostlaw[eidx].S).first;
        }

        // until the first iteration reports throughputs, a task splits its
        // requests evenly over its entries
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            const std::size_t n = lqn.entriesof[tidx].size();
            if (n == 0) continue;
            for (std::size_t eidx : lqn.entriesof[tidx])
                ph_share[eidx] = num_traits<T>::from_double(1.0 / double(n));
        }
    }

    /**
     * Law of the total time one execution of the issuing activity spends in call
     * CIDX: the geometric compound, of mean callproc_mean, of the response law
     * of the called entry. The response law is fitted to the response time
     * reported by the callee's layer and to the SCV of the callee's own composed
     * law, so no extra solver output is needed.
     */
    Distrib<T> ph_call_burst_law(std::size_t cidx) const {
        const double m = dbl(lqn.callproc_mean[cidx]);
        const std::size_t eidx = lqn.callpair_dst[cidx];
        if (m <= GlobalConstants::FineTol) return Distrib<T>::immediate();
        const T R = T(callservt[cidx] / lqn.callproc_mean[cidx]);
        double scvd = dbl(ph_entryscv[eidx]);
        if (!std::isfinite(scvd) || scvd <= GlobalConstants::FineTol) scvd = 1.0;
        const T Rf = dbl(R) > GlobalConstants::FineTol
                         ? R : num_traits<T>::from_double(GlobalConstants::FineTol);
        const Distrib<T> base = lang::aph_fit_mean_scv(Rf, num_traits<T>::from_double(scvd));
        const workflow::PhLaw<T> body = api::lqn::ph_law_of(base);
        const workflow::PhLaw<T> loop =
            workflow::Workflow<T>::compose_loop_geometric(body, lqn.callproc_mean[cidx]);
        return Distrib<T>::phase_type(loop.alpha, loop.S,
                                      workflow::Workflow<T>::is_acyclic_generator(loop.S));
    }

    /**
     * Station law of a composed workflow. A geometric loop over a body of two or
     * more phases closes a cycle in the phase graph, and a cyclic generator is a
     * PH and not an APH: no layer solver declares PH, so such a law is reduced
     * to the APH with the SAME first two moments. AMVA and NC read exactly those
     * two, so the reduction is lossless for them and is a two-moment fit for the
     * phase-aware layer solvers.
     */
    Distrib<T> ph_station_law(const workflow::PhLaw<T>& law) const {
        if (workflow::Workflow<T>::is_acyclic_generator(law.S))
            return Distrib<T>::phase_type(law.alpha, law.S, true);
        const std::pair<T, T> mm = api::lqn::ph_moments(law.alpha, law.S);
        return lang::aph_fit_mean_scv(mm.first, mm.second);
    }

    /** The law of a single immediate phase, the empty-composition answer. */
    static workflow::PhLaw<T> ph_immediate_law() {
        workflow::PhLaw<T> out;
        out.alpha.assign(1, Tone());
        out.S = Matrix<T>(1, 1, num_traits<T>::from_double(-GlobalConstants::Immediate));
        return out;
    }

    /**
     * Recompose the entry service laws from the current fixed-point iterate.
     *
     * The composed mean is NOT the sum of the leaf means when the graph forks:
     * the branches of an AND fork overlap, and the entry finishes with the last
     * of them. The ratio of the two, the overlap factor, is what the caller-side
     * aggregates are scaled by, so that the pieces of a cycle still add up to
     * the cycle.
     */
    void ph_compose_entry_laws() {
        const std::size_t N = lqn.nidx;
        std::vector<T> entry_setup_share(N + 1, Tzero());
        ph_overlap.assign(N + 1, Tone());

        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            if (!ph_has_wf[eidx]) continue;
            workflow::Workflow<T>& w = ph_wf[eidx];
            const std::unordered_map<std::size_t, T>& ex = ph_execs[eidx];
            T entrysum = Tzero(), procsum = Tzero();
            for (std::size_t aidx : lqn.actsof[eidx]) {
                T m = T(residt[aidx] + act_think_time(aidx));
                const T xa = ex.at(aidx);
                procsum += T(xa * m);
                const double md = dbl(m);
                w.set_activity_demand_mean(
                    lqn.names[aidx],
                    md > GlobalConstants::FineTol
                        ? m : num_traits<T>::from_double(GlobalConstants::FineTol));
                for (std::size_t cidx : lqn.callsof[aidx]) {
                    if (lqn.calltype[cidx] != CallType::SYNC) continue;
                    w.set_activity_demand(lqn.callhashnames[cidx], ph_call_burst_law(cidx));
                    m = T(m + callservt[cidx]);
                }
                entrysum += T(xa * m);
            }
            workflow::PhLaw<T> law = api::lqn::ph_law_of(w.refresh_ph());
            std::pair<T, T> mm = api::lqn::ph_moments(law.alpha, law.S);
            T m1 = mm.first, scv = mm.second;
            // All activities of an entry run on ONE processor, so the branches of
            // an AND fork cannot overlap the processor residence they request:
            // the composed maximum is a lower bound on the entry service time
            // only above that total. Where it falls below, the law is rescaled in
            // time to it, which keeps its shape, its SCV and its order.
            if (dbl(procsum) > dbl(m1) + GlobalConstants::FineTol) {
                const T f = T(m1 / procsum);
                for (std::size_t i = 0; i < law.S.rows(); ++i)
                    for (std::size_t j = 0; j < law.S.cols(); ++j) law.S(i, j) = T(law.S(i, j) * f);
                m1 = procsum;
            }
            // A SetupTask powers a thread down when it goes idle, so a request may
            // find it off and pay a cold start before the entry runs at all. The
            // setup is not part of the activity graph and never enters the
            // series-parallel reduction: it is prefixed to the composed law
            // afterwards as the mixture p*(setup THEN entry) + (1-p)*entry, which
            // is again phase-type.
            const double p = ph_setup_prob(eidx);
            if (p > GlobalConstants::FineTol) {
                const std::size_t tidx = lqn.parent[eidx];
                const double sm = dbl(lqn.setuptime[tidx].mean);
                double sscv = dbl(lqn.setuptime[tidx].scv);
                if (!std::isfinite(sscv) || sscv <= GlobalConstants::FineTol) sscv = 1.0;
                if (std::isfinite(sm) && sm > GlobalConstants::FineTol) {
                    const workflow::PhLaw<T> sl = api::lqn::ph_law_of(lang::aph_fit_mean_scv(
                        num_traits<T>::from_double(sm), num_traits<T>::from_double(sscv)));
                    std::vector<workflow::PhLaw<T>> mix;
                    mix.push_back(workflow::Workflow<T>::compose_serial(sl, law));
                    mix.push_back(law);
                    std::vector<T> probs;
                    probs.push_back(num_traits<T>::from_double(p));
                    probs.push_back(num_traits<T>::from_double(1.0 - p));
                    law = workflow::Workflow<T>::compose_mixture(mix, probs);
                    mm = api::lqn::ph_moments(law.alpha, law.S);
                    m1 = mm.first;
                    scv = mm.second;
                    // The share of the entry law that is cold start and not work.
                    // The surrogate-delay closure measures a thread's cycle in
                    // WORK, so it must not read a station utilization this has
                    // inflated -- see update_think_times_ph.
                    const double denom = std::max(dbl(m1), GlobalConstants::FineTol);
                    entry_setup_share[eidx] = num_traits<T>::from_double(p * sm / denom);
                }
            }
            ph_entrylaw[eidx] = law;
            ph_entrymean[eidx] = m1;
            ph_entryscv[eidx] = scv;
            if (dbl(entrysum) > GlobalConstants::FineTol) {
                const double r = std::min(1.0, dbl(m1) / dbl(entrysum));
                ph_overlap[eidx] = num_traits<T>::from_double(r);
            }
        }

        // Per task, the share-weighted fraction of its station service that is
        // cold start rather than work.
        ph_setupshare.assign(N + 1, Tzero());
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (ignore[tidx]) continue;
            for (std::size_t eidx : lqn.entriesof[tidx])
                ph_setupshare[tidx] = T(ph_setupshare[tidx] + ph_share[eidx] * entry_setup_share[eidx]);
        }

        // Expected number of calls per invocation, and the caller-side aggregates
        ph_ncalls = Matrix<T>(N + 1, N + 1, Tzero());
        ph_calltime = Matrix<T>(N + 1, N + 1, Tzero());
        ph_procresid.assign(N + 1, Tzero());
        ph_actthinkt.assign(N + 1, Tzero());
        ph_calltotal.assign(N + 1, Tzero());
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (ignore[tidx]) continue;
            for (std::size_t eidx : lqn.entriesof[tidx]) {
                if (!ph_has_wf[eidx]) continue;
                const T w = ph_share[eidx];
                if (!(dbl(w) > 0.0)) continue;
                const std::unordered_map<std::size_t, T>& ex = ph_execs[eidx];
                const T r = ph_overlap[eidx];
                for (std::size_t aidx : lqn.actsof[eidx]) {
                    const T xa = ex.at(aidx);
                    ph_procresid[tidx] = T(ph_procresid[tidx] + w * r * xa * residt[aidx]);
                    ph_actthinkt[tidx] = T(ph_actthinkt[tidx] + w * r * xa * act_think_time(aidx));
                    for (std::size_t cidx : lqn.callsof[aidx]) {
                        if (lqn.calltype[cidx] != CallType::SYNC) continue;
                        const std::size_t tgte = lqn.callpair_dst[cidx];
                        const std::size_t tgtt = lqn.parent[tgte];
                        // the COUNT of calls does not change with the overlap,
                        // only the time the caller is held by them
                        ph_ncalls(tidx, tgte) =
                            T(ph_ncalls(tidx, tgte) + w * xa * lqn.callproc_mean[cidx]);
                        ph_calltime(tidx, tgtt) =
                            T(ph_calltime(tidx, tgtt) + w * r * xa * callservt[cidx]);
                        ph_calltotal[tidx] = T(ph_calltotal[tidx] + w * r * xa * callservt[cidx]);
                    }
                }
            }
        }
    }

    /** Build the ensemble under method "srvn.ph". */
    void build_layers_ph(bool flat = false) {
        if (!ph_laws_ready) ph_assert_supported(flat);
        // The interlock correction rewrites the populations of the call classes,
        // which this method does not create: its callers reach the server in one
        // class each
        opt.interlocking = false;

        // A preceding probe has already composed the per-entry workflows; they do
        // not depend on the iterate, so they are not rebuilt here.
        if (!ph_laws_ready) ph_init_laws();

        // Seed the fixed point with the static demands, then compose the laws
        const std::size_t N = lqn.nidx;
        residt.assign(N + 1, Tzero());
        servt.assign(N + 1, Tzero());
        callservt.assign(lqn.ncalls + 1, Tzero());
        callresidt.assign(lqn.ncalls + 1, Tzero());
        for (std::size_t aidx = lqn.ashift + 1; aidx <= lqn.ashift + lqn.nacts; ++aidx)
            residt[aidx] = lqn.hostdem[aidx].disabled ? Tzero() : lqn.hostdem[aidx].mean;
        for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
            if (lqn.calltype[cidx] != CallType::SYNC && lqn.calltype[cidx] != CallType::ASYNC)
                continue;
            const std::size_t eidx = lqn.callpair_dst[cidx];
            callservt[cidx] = T(lqn.callproc_mean[cidx] * ph_hostmean[eidx]);
            callresidt[cidx] = callservt[cidx];
        }
        ph_compose_entry_laws();

        if (flat) {
            // ONE subnetwork holding every processor and every called task
            build_ph_flat_layer();
            tput.assign(N + 1, Tzero());
            util.assign(N + 1, Tzero());
            thinkt.assign(N + 1, Tzero());
            update_layers_ph(0);
            return;
        }

        std::vector<qn::Layer<T>> raw(NT() + 1);
        std::vector<bool> present(NT() + 1, false);

        for (std::size_t hidx = 1; hidx <= lqn.nhosts; ++hidx) {
            if (ignore[hidx]) continue;
            std::vector<std::size_t> callers;
            for (std::size_t tidx : lqn.tasksof[hidx]) {
                if (ignore[tidx]) continue;
                if (lqn.isref[tidx]) { callers.push_back(tidx); continue; }
                for (std::size_t eidx : lqn.entriesof[tidx])
                    if (ph_any_caller_of(eidx) || lqn.has_arrival[eidx]) {
                        callers.push_back(tidx);
                        break;
                    }
            }
            if (callers.empty()) continue;
            build_layer_ph(raw[hidx], hidx, callers, true);
            present[hidx] = true;
        }
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (ignore[tidx] || lqn.isref[tidx]) continue;
            std::vector<std::size_t> callers;
            for (std::size_t ct = 1; ct <= lqn.ntasks; ++ct) {
                const std::size_t c = lqn.tshift + ct;
                if (c == tidx || ignore[c]) continue;
                for (std::size_t e : lqn.entriesof[tidx])
                    if (lqn.issynccaller.get(c, e)) { callers.push_back(c); break; }
            }
            if (callers.empty() && ph_async_calls_into(tidx).empty()) continue;
            build_layer_ph(raw[tidx], tidx, callers, false);
            present[tidx] = true;
        }

        idxhash.assign(lqn.nidx + 1, -1);
        long next = 0;
        for (std::size_t i = 1; i <= NT(); ++i)
            if (present[i]) {
                idxhash[i] = next++;
                ensemble.push_back(std::move(raw[i]));
            }
        layer_init_sol.assign(ensemble.size(), Matrix<T>());
        build_fork_views();

        // install the initial laws, so that iteration 1 sees the seeded demands
        // rather than the placeholders the stations were created with
        tput.assign(N + 1, Tzero());
        util.assign(N + 1, Tzero());
        thinkt.assign(N + 1, Tzero());
        update_layers_ph(0);
    }

    /**
     * Build the ONE layer of method "flat.ph": a client delay plus a station for
     * every processor and every called task.
     *
     * A caller task is one closed class, and it visits each server it uses ONCE
     * per invocation, carrying there the composed law of the demand it places on
     * that server -- the same law method "srvn.ph" installs in the server's own
     * layer. What changes is that the servers now contend inside one network
     * instead of seeing each other through surrogate delays, so the client delay
     * keeps only the think times and whatever of the cycle this model does not
     * hold. That is the whole difference between the two encodings of the PH
     * composition, and it is why the reconstruction passes are shared verbatim.
     */
    void build_ph_flat_layer() {
        const std::vector<std::size_t> servers = ph_flat_server_set();
        const std::size_t nsrv = servers.size();

        qn::Layer<T> m;
        m.name = "FlatPH";
        m.clientIdx = m.add_station(qn::Station<T>{
            "Clients", NodeType::Delay, SchedStrategy::INF,
            std::numeric_limits<double>::infinity(), false, 0});
        const std::size_t clientNode = m.node_of_station(m.clientIdx);
        m.server_idx_of.assign(lqn.nidx + 1, 0);

        std::vector<std::size_t> station_of(lqn.nidx + 1, 0);
        std::vector<std::size_t> serverNode(nsrv, 0);
        std::vector<std::size_t> srvStation(nsrv, 0);
        for (std::size_t si = 0; si < nsrv; ++si) {
            const std::size_t idx = servers[si];
            const bool ishost = idx <= lqn.nhosts;
            qn::Station<T> st;
            st.name = lqn.hashnames[idx];
            st.nodetype = lqn.sched[idx] == SchedStrategy::INF ? NodeType::Delay : NodeType::Queue;
            st.sched = lqn.sched[idx];
            st.nservers = lqn.sched[idx] == SchedStrategy::INF
                              ? std::numeric_limits<double>::infinity()
                              : lqn.maxmult[idx];
            st.attr_ishost = ishost;
            st.attr_idx = idx;
            const std::size_t s = m.add_station(st);
            srvStation[si] = s;
            station_of[idx] = s;
            serverNode[si] = m.node_of_station(s);
            m.server_idx_of[idx] = s;
            if (ishost)
                m.host_stations.push_back(s);
            else
                m.task_stations.push_back(s);
        }
        // the scalar fallback of the station lookup, which no served element reaches
        m.serverIdx = station_of[servers[0]];
        m.flat = true;

        // Callers of each server, and the union of them, which becomes the class set
        std::vector<std::vector<std::size_t>> callers_of(lqn.nidx + 1);
        std::vector<std::size_t> all_callers;
        for (std::size_t si = 0; si < nsrv; ++si) {
            const std::size_t idx = servers[si];
            std::vector<std::size_t> cs;
            if (idx <= lqn.nhosts) {
                for (std::size_t tidx : lqn.tasksof[idx]) {
                    if (ignore[tidx]) continue;
                    if (lqn.isref[tidx]) { cs.push_back(tidx); continue; }
                    for (std::size_t eidx : lqn.entriesof[tidx])
                        if (ph_any_caller_of(eidx) || lqn.has_arrival[eidx]) {
                            cs.push_back(tidx);
                            break;
                        }
                }
            } else {
                for (std::size_t ct = 1; ct <= lqn.ntasks; ++ct) {
                    const std::size_t c = lqn.tshift + ct;
                    if (c == idx || ignore[c]) continue;
                    for (std::size_t e : lqn.entriesof[idx])
                        if (lqn.issynccaller.get(c, e)) { cs.push_back(c); break; }
                }
            }
            callers_of[idx] = cs;
            for (std::size_t c : cs)
                if (std::find(all_callers.begin(), all_callers.end(), c) == all_callers.end())
                    all_callers.push_back(c);
        }
        std::sort(all_callers.begin(), all_callers.end());

        // One closed class per caller task
        std::vector<std::size_t> class_of_caller(lqn.nidx + 1, 0);
        double npop = 0.0;
        for (std::size_t c : all_callers) {
            // ph_flat_server_set has refused every replicated element, so the
            // per-replica reduction the srvn builder makes is the identity here
            double nj = lqn.maxmult[c];
            if (std::isinf(nj)) {
                double sacc = 0.0;
                for (std::size_t k = 1; k <= NT(); ++k)
                    if (lqn.taskgraph.get(k, c) != Tzero()) sacc += lqn.maxmult[k];
                nj = sacc;
                if (std::isinf(nj) || nj == 0.0) {
                    double s2 = 0.0;
                    for (std::size_t k = 1; k <= NT(); ++k)
                        if (std::isfinite(lqn.maxmult[k])) s2 += lqn.maxmult[k];
                    nj = std::min(s2, 1000.0);
                }
            }
            qn::JobClass jc;
            jc.name = lqn.hashnames[c];
            jc.type = JobClassType::CLOSED;
            jc.population = nj;
            jc.refstat = m.clientIdx;
            jc.is_ref_class = true;
            jc.attr_kind = int(LqnElement::TASK);
            jc.attr_idx = c;
            const std::size_t k = m.add_class(jc);
            class_of_caller[c] = k;
            m.attr_tasks.emplace_back(k, c);
            npop += nj;
            const double zt = dbl(ref_think_time(c));
            m.set_service(m.clientIdx, k,
                          Distrib<T>::exp_mean(num_traits<T>::from_double(
                              std::max(zt, GlobalConstants::FineTol))));
            // A station this caller never reaches must say so with a DISABLED law,
            // not with a tiny placeholder. An FCFS station carries ONE service law
            // across its classes, so a placeholder is not inert there: it is mixed
            // into the multiserver correction and invents waiting where there is
            // none. Under "srvn.ph" the question never arises, since every class of
            // a layer visits that layer's single server.
            for (std::size_t si = 0; si < nsrv; ++si)
                m.set_service(srvStation[si], k, Distrib<T>::disabled_dist());
            for (std::size_t si = 0; si < nsrv; ++si) {
                const std::size_t idx = servers[si];
                const std::vector<std::size_t>& cs = callers_of[idx];
                if (std::find(cs.begin(), cs.end(), c) == cs.end()) continue;
                m.set_service(srvStation[si], k,
                              Distrib<T>::exp_mean(
                                  num_traits<T>::from_double(GlobalConstants::FineTol)));
                njobs(c, idx) = nj;
                thinkt_map.push_back({idx, c, m.clientIdx, k});
                servt_map.push_back({idx, c, station_of[idx], k});
            }
        }

        // Open classes: entry arrivals on a processor station, async calls on a task one
        std::vector<std::vector<std::pair<std::size_t, long>>> open_of(lqn.nidx + 1);
        std::size_t sourceStation = 0, sinkNode = 0;
        for (std::size_t si = 0; si < nsrv; ++si) {
            const std::size_t hidx = servers[si];
            if (hidx > lqn.nhosts) continue;
            for (std::size_t c : callers_of[hidx]) {
                // A task no other task calls has no task station, so the think-time
                // closure never gives its caller class a surrogate delay: the class
                // cycles against an Immediate one and an open stream on top of it
                // doubles the load. The chain is the representation that honours the
                // thread pool, so it is kept and closed on the arrival rate instead.
                if (ph_open_arrival_only(c)) continue;
                for (std::size_t eidx : lqn.entriesof[c]) {
                    if (!lqn.has_arrival[eidx]) continue;
                    if (sourceStation == 0) {
                        sourceStation = m.add_station(qn::Station<T>{
                            "Source", NodeType::Source, SchedStrategy::EXT,
                            std::numeric_limits<double>::infinity(), false, 0});
                        m.sourceIdx = sourceStation;
                        sinkNode = m.add_node("Sink", NodeType::Sink, false);
                        m.sinkNode = sinkNode;
                    }
                    qn::JobClass oc;
                    oc.name = lqn.hashnames[eidx] + ".Open";
                    oc.type = JobClassType::OPEN;
                    oc.population = std::numeric_limits<double>::infinity();
                    oc.refstat = sourceStation;
                    oc.attr_kind = int(LqnElement::ENTRY);
                    oc.attr_idx = eidx;
                    const std::size_t k = m.add_class(oc);
                    m.set_service(sourceStation, k, lqn.arrival[eidx]);
                    // disabled, not a placeholder, at every station this stream misses
                    for (std::size_t s2 = 0; s2 < nsrv; ++s2)
                        m.set_service(srvStation[s2], k, Distrib<T>::disabled_dist());
                    const double hm = std::max(dbl(ph_hostmean[eidx]), GlobalConstants::FineTol);
                    m.set_service(srvStation[si], k,
                                  Distrib<T>::exp_mean(num_traits<T>::from_double(hm)));
                    open_of[hidx].emplace_back(k, long(eidx));
                    m.attr_entries.emplace_back(k, eidx);
                }
            }
        }
        for (std::size_t si = 0; si < nsrv; ++si) {
            const std::size_t tidx = servers[si];
            if (tidx <= lqn.nhosts) continue;
            for (std::size_t cidx : ph_async_calls_into(tidx)) {
                if (sourceStation == 0) {
                    sourceStation = m.add_station(qn::Station<T>{
                        "Source", NodeType::Source, SchedStrategy::EXT,
                        std::numeric_limits<double>::infinity(), false, 0});
                    m.sourceIdx = sourceStation;
                    sinkNode = m.add_node("Sink", NodeType::Sink, false);
                    m.sinkNode = sinkNode;
                }
                const std::size_t eidx = lqn.callpair_dst[cidx];
                qn::JobClass oc;
                oc.name = lqn.callhashnames[cidx];
                oc.type = JobClassType::OPEN;
                oc.population = std::numeric_limits<double>::infinity();
                oc.refstat = sourceStation;
                oc.attr_kind = int(LqnElement::CALL);
                oc.attr_idx = cidx;
                const std::size_t k = m.add_class(oc);
                m.set_service(sourceStation, k, Distrib<T>::immediate());
                // disabled, not a placeholder, at every station this stream misses
                for (std::size_t s2 = 0; s2 < nsrv; ++s2)
                    m.set_service(srvStation[s2], k, Distrib<T>::disabled_dist());
                const double em = std::max(dbl(ph_entrymean[eidx]), GlobalConstants::FineTol);
                m.set_service(srvStation[si], k,
                              Distrib<T>::exp_mean(num_traits<T>::from_double(em)));
                open_of[tidx].emplace_back(k, -long(cidx));
                m.attr_calls.push_back({k, cidx, lqn.callpair_src[cidx], eidx});
                arv_call_map.push_back({tidx, cidx, sourceStation, k});
                call_map.push_back({tidx, cidx, station_of[tidx], k});
            }
        }

        // Routing: one visit per server the caller uses, in server order. The
        // number of calls is carried by the service law, not by a visit ratio, so
        // no arc ever moves.
        for (std::size_t c : all_callers) {
            const std::size_t k = class_of_caller[c];
            std::size_t prev = clientNode;
            bool visited = false;
            for (std::size_t si = 0; si < nsrv; ++si) {
                const std::vector<std::size_t>& cs = callers_of[servers[si]];
                if (std::find(cs.begin(), cs.end(), c) == cs.end()) continue;
                m.set_route(k, k, prev, serverNode[si], Tone());
                prev = serverNode[si];
                visited = true;
            }
            if (visited) m.set_route(k, k, prev, clientNode, Tone());
        }
        if (sourceStation != 0) {
            const std::size_t srcNode = m.node_of_station(sourceStation);
            for (std::size_t si = 0; si < nsrv; ++si)
                for (const auto& oa : open_of[servers[si]]) {
                    m.set_route(oa.first, oa.first, srcNode, serverNode[si], Tone());
                    m.set_route(oa.first, oa.first, serverNode[si], sinkNode, Tone());
                }
        }
        m.refresh_chains();

        idxhash.assign(lqn.nidx + 1, -1);
        for (std::size_t idx : servers) idxhash[idx] = 0;
        ensemble.clear();
        ensemble.push_back(std::move(m));
        layer_init_sol.assign(ensemble.size(), Matrix<T>());
        build_fork_views();

        const std::size_t nclasses = ensemble[0].classes.size();
        for (std::size_t idx : servers) {
            PHLayer L;
            L.idx = idx;
            L.ishost = idx <= lqn.nhosts;
            L.callers = callers_of[idx];
            L.class_of_caller = class_of_caller;
            L.nreplicas = 1;
            L.qstations.assign(1, station_of[idx]);
            L.svcmean_by_class.assign(nclasses + 1, Tzero());
            L.open_arrivals = open_of[idx];
            L.npop = npop < 1.0 ? 1.0 : npop;
            ph_layers[idx] = L;
            ph_has_layer[idx] = true;
        }
    }

    /**
     * Processors and called tasks that become stations of the flat layer.
     *
     * The set is the elements the srvn builder would have given a layer of their
     * own, so "flat.ph" and "srvn.ph" place the SAME stations and differ only in
     * how many networks hold them. The refusals are those of flat_server_set,
     * since they are properties of the squashing and not of the encoding: each of
     * these carries per-layer state that one submodel cannot hold.
     */
    std::vector<std::size_t> ph_flat_server_set() const {
        for (std::size_t i = 1; i <= NT(); ++i) {
            if (lqn.repl[i] > 1.0)
                throw UnsupportedError(
                    "method='flat.ph' does not support replicated processors or tasks, whose "
                    "replicas need a submodel each. Use method='srvn.ph'.");
            if (lqn.iscache[i])
                throw UnsupportedError(
                    "method='flat.ph' does not support cache tasks. Use method='default'.");
            if (lqn.hassetup[i])
                throw UnsupportedError(
                    "method='flat.ph' does not support setup tasks, whose powered-down threads "
                    "are per-layer state. Use method='srvn.ph'.");
        }
        std::vector<std::size_t> servers;
        for (std::size_t hidx = 1; hidx <= lqn.nhosts; ++hidx) {
            if (ignore[hidx] || lqn.tasksof[hidx].empty()) continue;
            bool any = false;
            for (std::size_t tidx : lqn.tasksof[hidx]) {
                if (ignore[tidx]) continue;
                if (lqn.isref[tidx]) { any = true; break; }
                for (std::size_t eidx : lqn.entriesof[tidx])
                    if (ph_any_caller_of(eidx) || lqn.has_arrival[eidx]) { any = true; break; }
                if (any) break;
            }
            if (any) servers.push_back(hidx);
        }
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (ignore[tidx] || lqn.isref[tidx]) continue;
            bool has_caller = false;
            for (std::size_t ct = 1; ct <= lqn.ntasks && !has_caller; ++ct) {
                const std::size_t c = lqn.tshift + ct;
                if (c == tidx || ignore[c]) continue;
                for (std::size_t e : lqn.entriesof[tidx])
                    if (lqn.issynccaller.get(c, e)) { has_caller = true; break; }
            }
            if (!has_caller && ph_async_calls_into(tidx).empty()) continue;
            servers.push_back(tidx);
        }
        if (servers.empty())
            throw InputError(
                "method='flat.ph' found no server: the model has no processor with tasks.");
        return servers;
    }

    /** Build the two-station layer of server element IDX. */
    void build_layer_ph(qn::Layer<T>& m, std::size_t idx,
                        const std::vector<std::size_t>& callers, bool ishost) {
        m.name = lqn.hashnames[idx];

        // Replicas of the server station, with the same fan-out reduction as the
        // default builder: a caller that reaches every replica sees one
        // representative
        const double rawrepl = lqn.repl[idx];
        std::size_t nreplicas = 1;
        if (rawrepl > 1.0 && !callers.empty()) {
            bool reduce = true;
            if (ishost) {
                for (std::size_t c : callers)
                    if (lqn.repl[c] != rawrepl) reduce = false;
            } else {
                for (std::size_t c : callers)
                    if (lqn.fanout_at(c, idx) < rawrepl) reduce = false;
            }
            nreplicas = reduce ? 1 : static_cast<std::size_t>(std::llround(rawrepl));
            if (reduce && !ishost) single_replica_tasks.insert(idx);
        }
        const bool reduce_fanout = (nreplicas == 1 && rawrepl > 1.0 && !callers.empty());

        m.clientIdx = m.add_station(qn::Station<T>{
            "Clients", NodeType::Delay, SchedStrategy::INF,
            std::numeric_limits<double>::infinity(), false, 0});
        const std::size_t clientNode = m.node_of_station(m.clientIdx);
        m.serverIdx = m.clientIdx + 1;
        PHLayer L;
        L.idx = idx;
        L.ishost = ishost;
        L.callers = callers;
        L.nreplicas = nreplicas;
        L.class_of_caller.assign(lqn.nidx + 1, 0);
        std::vector<std::size_t> serverNode(nreplicas);
        for (std::size_t r = 0; r < nreplicas; ++r) {
            qn::Station<T> st;
            st.name = r == 0 ? lqn.hashnames[idx] : lqn.hashnames[idx] + "." + std::to_string(r + 1);
            st.nodetype = lqn.sched[idx] == SchedStrategy::INF ? NodeType::Delay : NodeType::Queue;
            st.sched = lqn.sched[idx];
            st.nservers = lqn.sched[idx] == SchedStrategy::INF
                              ? std::numeric_limits<double>::infinity()
                              : lqn.maxmult[idx];
            st.attr_ishost = ishost;
            st.attr_idx = idx;
            const std::size_t s = m.add_station(st);
            L.qstations.push_back(s);
            serverNode[r] = m.node_of_station(s);
        }

        // --- closed class per caller task
        for (std::size_t c : callers) {
            double nj = njobs(c, idx);
            if (nj == 0.0) {
                const bool caller_single_replica =
                    reduce_fanout || single_replica_tasks.count(c) > 0;
                nj = caller_single_replica ? lqn.maxmult[c] : lqn.maxmult[c] * lqn.repl[c];
                if (std::isinf(nj)) {
                    double s = 0.0;
                    for (std::size_t k = 1; k <= NT(); ++k)
                        if (lqn.taskgraph.get(k, c) != Tzero()) s += lqn.maxmult[k];
                    nj = s;
                    if (std::isinf(nj) || nj == 0.0) {
                        double s2 = 0.0;
                        for (std::size_t k = 1; k <= NT(); ++k)
                            if (std::isfinite(lqn.maxmult[k])) s2 += lqn.maxmult[k] * lqn.repl[k];
                        nj = std::min(s2, 1000.0);
                    }
                }
                njobs(c, idx) = nj;
            }
            qn::JobClass jc;
            jc.name = lqn.hashnames[c];
            jc.type = JobClassType::CLOSED;
            jc.population = nj;
            jc.refstat = m.clientIdx;
            jc.is_ref_class = true;
            jc.attr_kind = int(LqnElement::TASK);
            jc.attr_idx = c;
            const std::size_t k = m.add_class(jc);
            L.class_of_caller[c] = k;
            m.attr_tasks.emplace_back(k, c);
            const double zt = dbl(ref_think_time(c));
            m.set_service(m.clientIdx, k,
                          Distrib<T>::exp_mean(num_traits<T>::from_double(
                              std::max(zt, GlobalConstants::FineTol))));
            for (std::size_t s : L.qstations)
                m.set_service(s, k, Distrib<T>::exp_mean(
                                        num_traits<T>::from_double(GlobalConstants::FineTol)));
            // every layer must be refreshed after a law change: post() resets the
            // layers named by the think-time map
            thinkt_map.push_back({idx, c, m.clientIdx, k});
            servt_map.push_back({idx, c, L.qstations[0], k});
        }

        // --- open classes: entry arrivals on a host layer, async calls on a task layer
        std::size_t sourceStation = 0, sinkNode = 0;
        if (ishost) {
            for (std::size_t c : callers) {
                // A task no other task calls has no task layer, so
                // update_think_times_ph never gives its caller class a surrogate
                // delay: the class cycles against an Immediate one and an open
                // stream on top of it doubles the load. The chain is the
                // representation that honours the thread pool, so it is kept and
                // closed on the arrival rate instead.
                if (ph_open_arrival_only(c)) continue;
                for (std::size_t eidx : lqn.entriesof[c]) {
                    if (!lqn.has_arrival[eidx]) continue;
                    if (sourceStation == 0) {
                        sourceStation = m.add_station(qn::Station<T>{
                            "Source", NodeType::Source, SchedStrategy::EXT,
                            std::numeric_limits<double>::infinity(), false, 0});
                        m.sourceIdx = sourceStation;
                        sinkNode = m.add_node("Sink", NodeType::Sink, false);
                        m.sinkNode = sinkNode;
                    }
                    qn::JobClass oc;
                    oc.name = lqn.hashnames[eidx] + ".Open";
                    oc.type = JobClassType::OPEN;
                    oc.population = std::numeric_limits<double>::infinity();
                    oc.refstat = sourceStation;
                    oc.attr_kind = int(LqnElement::ENTRY);
                    oc.attr_idx = eidx;
                    const std::size_t k = m.add_class(oc);
                    m.set_service(sourceStation, k, lqn.arrival[eidx]);
                    for (std::size_t s : L.qstations) {
                        const double hm = std::max(dbl(ph_hostmean[eidx]), GlobalConstants::FineTol);
                        m.set_service(s, k, Distrib<T>::exp_mean(num_traits<T>::from_double(hm)));
                    }
                    L.open_arrivals.emplace_back(k, long(eidx));
                    m.attr_entries.emplace_back(k, eidx);
                }
            }
        } else {
            for (std::size_t cidx : ph_async_calls_into(idx)) {
                if (sourceStation == 0) {
                    sourceStation = m.add_station(qn::Station<T>{
                        "Source", NodeType::Source, SchedStrategy::EXT,
                        std::numeric_limits<double>::infinity(), false, 0});
                    m.sourceIdx = sourceStation;
                    sinkNode = m.add_node("Sink", NodeType::Sink, false);
                    m.sinkNode = sinkNode;
                }
                const std::size_t eidx = lqn.callpair_dst[cidx];
                qn::JobClass oc;
                oc.name = lqn.callhashnames[cidx];
                oc.type = JobClassType::OPEN;
                oc.population = std::numeric_limits<double>::infinity();
                oc.refstat = sourceStation;
                oc.attr_kind = int(LqnElement::CALL);
                oc.attr_idx = cidx;
                const std::size_t k = m.add_class(oc);
                m.set_service(sourceStation, k, Distrib<T>::immediate());
                for (std::size_t s : L.qstations) {
                    const double em = std::max(dbl(ph_entrymean[eidx]), GlobalConstants::FineTol);
                    m.set_service(s, k, Distrib<T>::exp_mean(num_traits<T>::from_double(em)));
                }
                L.open_arrivals.emplace_back(k, -long(cidx));
                m.attr_calls.push_back({k, cidx, lqn.callpair_src[cidx], eidx});
                arv_call_map.push_back({idx, cidx, sourceStation, k});
                call_map.push_back({idx, cidx, L.qstations[0], k});
            }
        }

        // Routing: one visit to the server per client cycle. The number of calls
        // is carried by the service law, not by a visit ratio, so no arc changes
        const T share = num_traits<T>::from_double(1.0 / double(nreplicas));
        for (std::size_t c : callers) {
            const std::size_t k = L.class_of_caller[c];
            for (std::size_t r = 0; r < nreplicas; ++r) {
                m.set_route(k, k, clientNode, serverNode[r], share);
                m.set_route(k, k, serverNode[r], clientNode, Tone());
            }
        }
        for (const auto& oa : L.open_arrivals) {
            const std::size_t k = oa.first;
            const std::size_t srcNode = m.node_of_station(sourceStation);
            for (std::size_t r = 0; r < nreplicas; ++r) {
                m.set_route(k, k, srcNode, serverNode[r], share);
                m.set_route(k, k, serverNode[r], sinkNode, Tone());
            }
        }
        L.svcmean_by_class.assign(m.classes.size() + 1, Tzero());
        double np = 0.0;
        for (std::size_t c : callers) {
            const double v = njobs(c, idx);
            if (std::isfinite(v) && v > 0.0) np += v;
        }
        L.npop = np < 1.0 ? 1.0 : np;
        m.refresh_chains();
        ph_layers[idx] = L;
        ph_has_layer[idx] = true;
    }

    /** Law of the demand caller C places on the server of layer IDX per invocation. */
    workflow::PhLaw<T> ph_service_law(std::size_t idx, bool ishost, std::size_t c) const {
        if (ishost) {
            // mixture over the entries of C, weighted by their share of its requests
            std::vector<workflow::PhLaw<T>> laws;
            std::vector<T> probs;
            double tot = 0.0;
            for (std::size_t eidx : lqn.entriesof[c]) {
                if (ph_hostlaw[eidx].S.rows() == 0 || !(dbl(ph_share[eidx]) > 0.0)) continue;
                laws.push_back(ph_hostlaw[eidx]);
                probs.push_back(ph_share[eidx]);
                tot += dbl(ph_share[eidx]);
            }
            if (laws.empty()) return ph_immediate_law();
            for (T& p : probs) p = T(p / num_traits<T>::from_double(tot));
            return workflow::Workflow<T>::compose_mixture(laws, probs);
        }
        // task layer: the total demand is the sum, over the entries of the
        // server, of a geometric compound of the entry law of mean equal to the
        // number of calls
        bool started = false;
        workflow::PhLaw<T> out;
        for (std::size_t eidx : lqn.entriesof[idx]) {
            const T n = ph_ncalls(c, eidx);
            if (dbl(n) <= GlobalConstants::FineTol || ph_entrylaw[eidx].S.rows() == 0) continue;
            const workflow::PhLaw<T> lp =
                workflow::Workflow<T>::compose_loop_geometric(ph_entrylaw[eidx], n);
            out = started ? workflow::Workflow<T>::compose_serial(out, lp) : lp;
            started = true;
        }
        if (!started) return ph_immediate_law();
        return out;
    }

    /**
     * Mean time a thread of caller C spends away from the server of layer IDX per
     * invocation: idle, plus whatever of its cycle the layer does not hold.
     */
    T ph_delay_mean(std::size_t idx, std::size_t c) const {
        // ONE closure for both layerings. Under "srvn.ph" the model holds a single
        // server, so a host layer charges the whole call burst to the delay and a
        // task layer charges the caller's processor plus every other callee. Under
        // "flat.ph" the model holds every server, and only the think times are
        // left. Every term is SUMMED in rather than obtained by subtracting from a
        // total: that subtraction cancels catastrophically once a call time is
        // large, the think time falls below the ULP of the call time, and the layer
        // then sees a client delay of zero, saturates, and the fixed point runs
        // away.
        T z = T(thinkt[c] + ref_think_time(c));
        if (!std::isfinite(dbl(z)) || dbl(z) < 0.0) z = Tzero();
        z = T(z + ph_actthinkt[c]);
        const std::size_t hidx = lqn.parent[c];
        if (!ph_served_here(idx, hidx)) z = T(z + ph_procresid[c]);
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (!ph_served_here(idx, tidx)) z = T(z + ph_calltime(c, tidx));
        }
        if (!std::isfinite(dbl(z)) || dbl(z) < 0.0)
            z = num_traits<T>::from_double(GlobalConstants::FineTol);
        return z;
    }

    /**
     * True when LQN element ELEM is a station of the same model that holds server
     * IDX. Under "srvn.ph" that is ELEM == IDX, since each server has a layer of
     * its own; under "flat.ph" it is every server of the one network.
     */
    bool ph_served_here(std::size_t idx, std::size_t elem) const {
        if (elem < 1 || elem >= idxhash.size() || idx < 1 || idx >= idxhash.size())
            return false;
        return idxhash[elem] >= 0 && idxhash[elem] == idxhash[idx];
    }

    /**
     * Push the composed laws into the layers.
     *
     * A layer of this method carries no routing that depends on the iterate: the
     * number of calls a caller makes is folded into its service law rather than
     * into a visit ratio, so only two laws move per (layer, class) -- the
     * phase-type service law at the server and the mean of the surrogate delay
     * at the client.
     */
    void update_layers_ph(int) {
        for (std::size_t idx = 1; idx <= NT(); ++idx) {
            if (idxhash[idx] < 0 || !ph_has_layer[idx]) continue;
            PHLayer& L = ph_layers[idx];
            qn::Layer<T>& m = ensemble[std::size_t(idxhash[idx])];
            for (std::size_t c : L.callers) {
                const std::size_t k = L.class_of_caller[c];
                const workflow::PhLaw<T> sl = ph_service_law(idx, L.ishost, c);
                L.svcmean_by_class[k] = api::lqn::ph_moments(sl.alpha, sl.S).first;
                const Distrib<T> law = ph_station_law(sl);
                for (std::size_t s : L.qstations) m.set_service(s, k, law);
                const double zd = std::max(dbl(ph_delay_mean(idx, c)),
                                           GlobalConstants::FineTol);
                m.set_service(m.clientIdx, k,
                              Distrib<T>::exp_mean(num_traits<T>::from_double(zd)));
            }
            for (const auto& oa : L.open_arrivals) {
                const std::size_t k = oa.first;
                if (oa.second > 0) {
                    // entry arrival: the processor demand law of the entry is static
                    L.svcmean_by_class[k] = ph_hostmean[std::size_t(oa.second)];
                    continue;
                }
                const std::size_t cidx = std::size_t(-oa.second);
                const std::size_t eidx = lqn.callpair_dst[cidx];
                L.svcmean_by_class[k] = ph_entrymean[eidx];
                const Distrib<T> law = ph_station_law(ph_entrylaw[eidx]);
                for (std::size_t s : L.qstations) m.set_service(s, k, law);
                const std::size_t aidx = lqn.callpair_src[cidx];
                double rate = dbl(tput[aidx]) * dbl(lqn.callproc_mean[cidx]);
                if (!std::isfinite(rate) || rate <= GlobalConstants::FineTol)
                    rate = GlobalConstants::FineTol;
                m.set_service(m.sourceIdx, k,
                              Distrib<T>::exp_rate(num_traits<T>::from_double(rate)));
            }
            m.refresh_rt();
        }
    }

    /**
     * Residence time per visit, by Little from the queue length rather than from
     * the reported RN. A layer that saturates can come back from AMVA with an RN
     * no closed model can produce, and a reconstruction that trusts it feeds the
     * impossible value straight back into the call response times.
     */
    static T ph_residence(const T& Q, const T& X, const T& RN) {
        const double q = num_traits<T>::to_double(Q), x = num_traits<T>::to_double(X);
        if (std::isfinite(q) && q >= 0.0 && std::isfinite(x) && x > GlobalConstants::FineTol)
            return T(Q / X);
        return RN;
    }

    /**
     * Ratio of a residence time to the mean of the law it was measured against,
     * bounded above by the layer population: a job can wait behind at most every
     * other job in a closed layer.
     */
    static T ph_inflation_of(const T& R, const T& S, double npop) {
        double f = 1.0;
        const double s = num_traits<T>::to_double(S), r = num_traits<T>::to_double(R);
        if (s > GlobalConstants::FineTol && std::isfinite(r) && r > 0.0) f = r / s;
        if (!std::isfinite(f) || f < 1.0) f = 1.0;
        if (std::isfinite(npop) && npop >= 1.0 && f > npop) f = npop;
        return num_traits<T>::from_double(f);
    }

    /** Total closed population of a layer, i.e. how many jobs a job can queue behind. */
    double ph_layer_pop(const PHLayer& L, std::size_t idx) const {
        // Under "flat.ph" this is every caller of the single network and not only
        // the callers of this one station, so it is taken from the layer record.
        if (L.npop >= 1.0) return L.npop;
        double n = 0.0;
        for (std::size_t c : L.callers) {
            const double v = njobs(c, idx);
            if (std::isfinite(v) && v > 0.0) n += v;
        }
        return n < 1.0 ? 1.0 : n;
    }

    /**
     * Reconstruct the LQN metrics.
     *
     * A layer of this method reports one row per caller task, not one per entry,
     * activity and call, so the per-element quantities the rest of SolverLN reads
     * -- servt, residt, callservt, callresidt, tput -- are recovered analytically
     * from the series-parallel weights of the entry workflows.
     *
     * The split is conservative by construction. A station reports a residence
     * time R per visit against a service law of mean S, so the queueing inflation
     * R/S is attributed to every leaf of that visit in proportion to its own
     * mean: the pieces sum back to R exactly.
     */
    void update_metrics_ph(int it) {
        const std::size_t N = lqn.nidx;
        servt.assign(N + 1, Tzero());
        residt.assign(N + 1, Tzero());
        callservt.assign(lqn.ncalls + 1, Tzero());
        callresidt.assign(lqn.ncalls + 1, Tzero());

        std::vector<T> inflNum(N + 1, Tzero()), inflDen(N + 1, Tzero());
        std::vector<T> taskTput(N + 1, Tzero()), openTput(N + 1, Tzero());

        // Host layers: the queueing inflation of the processor demand
        for (std::size_t hidx = 1; hidx <= lqn.nhosts; ++hidx) {
            if (idxhash[hidx] < 0 || !ph_has_layer[hidx]) continue;
            const PHLayer& L = ph_layers[hidx];
            const LayerResult<T>& res = results.back()[std::size_t(idxhash[hidx])];
            const double npop = ph_layer_pop(L, hidx);
            for (std::size_t c : L.callers) {
                const std::size_t k = L.class_of_caller[c];
                T X = Tzero(), Q = Tzero();
                for (std::size_t s : L.qstations) {
                    X = T(X + res.TN(s - 1, k - 1));
                    Q = T(Q + res.QN(s - 1, k - 1));
                }
                const T R = ph_residence(Q, X, res.RN(L.qstations[0] - 1, k - 1));
                const T f = ph_inflation_of(R, L.svcmean_by_class[k], npop);
                if (!std::isfinite(dbl(X)) || dbl(X) < 0.0) X = Tzero();
                // TOTAL over the replicas. The processor layer of a replicated
                // element models ONE representative replica, so X is one replica's
                // rate and the element's own rate is REPL times it. The matching
                // per-replica quantity is ph_xdemand, which the think-time closure
                // divides down for the same reason.
                taskTput[c] = T(taskTput[c] + num_traits<T>::from_double(lqn.repl[c]) * X);
                for (std::size_t eidx : lqn.entriesof[c]) {
                    const T sh = dbl(ph_share[eidx]) > 0.0 ? ph_share[eidx] : Tzero();
                    const T w = T(sh * X);
                    inflNum[eidx] = T(inflNum[eidx] + w * f);
                    inflDen[eidx] = T(inflDen[eidx] + w);
                }
            }
            for (const auto& oa : L.open_arrivals) {
                if (oa.second <= 0) continue;  // an async call is served in the task layer
                const std::size_t k = oa.first;
                const std::size_t eidx = std::size_t(oa.second);
                T X = Tzero(), Q = Tzero();
                for (std::size_t s : L.qstations) {
                    X = T(X + res.TN(s - 1, k - 1));
                    Q = T(Q + res.QN(s - 1, k - 1));
                }
                if (!std::isfinite(dbl(X)) || dbl(X) <= 0.0) continue;
                const T f = ph_inflation_of(ph_residence(Q, X, res.RN(L.qstations[0] - 1, k - 1)),
                                            L.svcmean_by_class[k], npop);
                inflNum[eidx] = T(inflNum[eidx] + X * f);
                inflDen[eidx] = T(inflDen[eidx] + X);
                openTput[eidx] = T(openTput[eidx] + X);
                taskTput[lqn.parent[eidx]] = T(taskTput[lqn.parent[eidx]] + X);
            }
        }

        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            double f = 1.0;
            if (dbl(inflDen[eidx]) > GlobalConstants::FineTol)
                f = dbl(inflNum[eidx]) / dbl(inflDen[eidx]);
            if (!std::isfinite(f) || f < 1.0) f = 1.0;  // a residence cannot fall below its demand
            const T fv = num_traits<T>::from_double(f);
            for (std::size_t aidx : lqn.actsof[eidx])
                residt[aidx] = T(fv * (lqn.hostdem[aidx].disabled ? Tzero()
                                                                  : lqn.hostdem[aidx].mean));
        }

        // Task layers: the response time of every call
        std::vector<T> relw(N + 1, Tzero());
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (idxhash[tidx] < 0 || !ph_has_layer[tidx]) continue;
            const PHLayer& L = ph_layers[tidx];
            const LayerResult<T>& res = results.back()[std::size_t(idxhash[tidx])];
            const double npop = ph_layer_pop(L, tidx);
            for (std::size_t c : L.callers) {
                const std::size_t k = L.class_of_caller[c];
                T X = Tzero(), Q = Tzero();
                for (std::size_t s : L.qstations) {
                    X = T(X + res.TN(s - 1, k - 1));
                    Q = T(Q + res.QN(s - 1, k - 1));
                }
                if (!std::isfinite(dbl(X)) || dbl(X) < 0.0) X = Tzero();
                const T g = ph_inflation_of(ph_residence(Q, X, res.RN(L.qstations[0] - 1, k - 1)),
                                            L.svcmean_by_class[k], npop);
                for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
                    if (lqn.calltype[cidx] != CallType::SYNC) continue;
                    if (lqn.parent[lqn.callpair_src[cidx]] != c) continue;
                    if (lqn.parent[lqn.callpair_dst[cidx]] != tidx) continue;
                    const std::size_t eidx = lqn.callpair_dst[cidx];
                    callservt[cidx] = T(lqn.callproc_mean[cidx] * g * ph_entrymean[eidx]);
                    callresidt[cidx] = callservt[cidx];
                }
                for (std::size_t eidx : lqn.entriesof[tidx])
                    relw[eidx] = T(relw[eidx] + X * ph_ncalls(c, eidx));
            }
            for (const auto& oa : L.open_arrivals) {
                if (oa.second >= 0) continue;
                const std::size_t k = oa.first;
                const std::size_t cidx = std::size_t(-oa.second);
                const std::size_t eidx = lqn.callpair_dst[cidx];
                T X = Tzero();
                for (std::size_t s : L.qstations) X = T(X + res.TN(s - 1, k - 1));
                const T R = res.RN(L.qstations[0] - 1, k - 1);
                if (std::isfinite(dbl(R)) && dbl(R) > 0.0) {
                    callservt[cidx] = T(R * lqn.callproc_mean[cidx]);
                    callresidt[cidx] = callservt[cidx];
                }
                if (std::isfinite(dbl(X)) && dbl(X) > 0.0) relw[eidx] = T(relw[eidx] + X);
            }
        }

        // Entry shares and throughputs
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            const std::vector<std::size_t>& entries = lqn.entriesof[tidx];
            if (entries.empty()) continue;
            // How the requests SPLIT over the entries is a flow-balance question,
            // and is answered at the task layer: a caller class reaches that server
            // once per invocation of the caller, carrying its whole call burst in
            // its service law, so the station rate counts caller cycles and the
            // per-entry rate is that rate times the calls the caller makes.
            T tot = Tzero();
            for (std::size_t eidx : entries) tot = T(tot + relw[eidx] + openTput[eidx]);
            if (dbl(tot) > GlobalConstants::FineTol) {
                for (std::size_t eidx : entries)
                    ph_share[eidx] = T((relw[eidx] + openTput[eidx]) / tot);
            } else {
                for (std::size_t eidx : entries)
                    ph_share[eidx] = num_traits<T>::from_double(1.0 / double(entries.size()));
            }
            // HOW MANY requests the task completes is a different question, and the
            // flow-balance total does not answer it: that total is what the callers
            // DEMAND, not what the task's threads can deliver. A thread cycles
            // through its host demand AND then through the task think time, and only
            // the processor layer of the task carries both, so the rate is read there.
            tput[tidx] = dbl(taskTput[tidx]) > GlobalConstants::FineTol ? taskTput[tidx] : tot;
            for (std::size_t eidx : entries) tput[eidx] = T(tput[tidx] * ph_share[eidx]);
            // The DEMAND is kept apart because it, and not the rate just reported,
            // is what closes the surrogate delay: normalising the think time by a
            // rate the same think time produced makes the processor layer
            // self-referential. PER REPLICA, because the thread count it is paired
            // with there is per replica.
            const T nrep = num_traits<T>::from_double(std::max(1.0, lqn.repl[tidx]));
            ph_xdemand[tidx] = dbl(tot) > GlobalConstants::FineTol ? T(tot / nrep)
                                                                   : T(tput[tidx] / nrep);
        }

        // Recovery, under-relaxation, and the derived per-element quantities
        for (std::size_t aidx = lqn.ashift + 1; aidx <= lqn.ashift + lqn.nacts; ++aidx) {
            T v = residt[aidx];
            if (!std::isfinite(dbl(v)) && it > 1 && !std::isnan(residt_prev[aidx]))
                v = residt_prev_v[aidx];
            if (relax_omega < 1.0 && it > 1 && !std::isnan(residt_prev[aidx])) {
                const T om = num_traits<T>::from_double(relax_omega);
                const T om1 = num_traits<T>::from_double(1.0 - relax_omega);
                v = T(om * v + om1 * residt_prev_v[aidx]);
            }
            residt[aidx] = v;
            residt_prev[aidx] = dbl(v);
            residt_prev_v[aidx] = v;
        }
        for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
            T v = callservt[cidx];
            if (!std::isfinite(dbl(v)))
                v = (it > 1 && std::isfinite(callservt_prev[cidx])) ? callservt_prev_v[cidx]
                                                                    : Tzero();
            if (relax_omega < 1.0 && it > 1 && !std::isnan(callservt_prev[cidx])) {
                const T om = num_traits<T>::from_double(relax_omega);
                const T om1 = num_traits<T>::from_double(1.0 - relax_omega);
                v = T(om * v + om1 * callservt_prev_v[cidx]);
            }
            callservt[cidx] = v;
            callresidt[cidx] = v;
            callservt_prev[cidx] = dbl(v);
            callresidt_prev[cidx] = dbl(v);
            callservt_prev_v[cidx] = v;
            if (dbl(v) > 0.0) callservtproc[cidx] = Distrib<T>::exp_mean(v);
        }

        // Recompose the entry laws from the iterate just computed. The entry
        // service time is then the mean of the COMPOSED law and not the sum of the
        // parts: the branches of an AND fork overlap, so an entry that forks
        // finishes with the last of its branches and is not charged their sum.
        ph_compose_entry_laws();

        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            if (!ph_has_wf[eidx]) continue;
            const std::unordered_map<std::size_t, T>& ex = ph_execs[eidx];
            for (std::size_t aidx : lqn.actsof[eidx]) {
                T sa = T(residt[aidx] + act_think_time(aidx));
                for (std::size_t cidx : lqn.callsof[aidx])
                    if (lqn.calltype[cidx] == CallType::SYNC) sa = T(sa + callservt[cidx]);
                servt[aidx] = sa;
                servt_prev[aidx] = dbl(sa);
                servt_prev_v[aidx] = sa;
                tput[aidx] = T(tput[eidx] * ex.at(aidx));
                tput_prev[aidx] = dbl(tput[aidx]);
                tput_prev_v[aidx] = tput[aidx];
                // exp_rate admits a null rate, but a never-called activity has no arrivals: disabled says so, as the python twin does.
                tputproc[aidx] = dbl(tput[aidx]) > 0.0 ? Distrib<T>::exp_rate(tput[aidx])
                                                       : Distrib<T>::disabled_dist();
                if (dbl(sa) > 0.0) servtproc[aidx] = Distrib<T>::exp_mean(sa);
            }
            servt[eidx] = ph_entrymean[eidx];
            residt[eidx] = ph_entrymean[eidx];
            if (dbl(servt[eidx]) > 0.0) servtproc[eidx] = Distrib<T>::exp_mean(servt[eidx]);
        }
    }

    /**
     * Surrogate delay of every caller.
     *
     * Same closure as update_think_times -- a thread of the task is idle for
     * whatever of its cycle the task's own station does not hold -- but the rate
     * it is normalised by is the INVOCATION rate of the task and not the
     * throughput of its station. Under this method a caller class reaches the
     * server once per invocation of the caller, carrying its whole call burst in
     * its service law, so the station rate counts caller cycles rather than calls
     * and the two differ by the mean number of calls.
     */
    void update_think_times_ph(int it) {
        thinktproc.assign(lqn.nidx + 1, Distrib<T>::disabled_dist());
        const T floorv = num_traits<T>::from_double(GlobalConstants::Zero);
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (ignore[tidx]) continue;
            const T ztask = ref_think_time(tidx);
            if (idxhash[tidx] < 0) {
                // A task no other task calls but whose entries carry an arrival
                // still has a cycle: its threads are driven by the stream.
                // build_layers_ph drops the open class for it precisely so this
                // closure can set the rate.
                const double arvrate = open_arrival_rate_of(tidx);
                if (arvrate > GlobalConstants::FineTol) {
                    double nja = lqn.maxmult[tidx];
                    if (!std::isfinite(nja) || nja <= 0.0) {
                        nja = 0.0;
                        for (std::size_t c = 1; c <= NT(); ++c) nja = std::max(nja, njobs(tidx, c));
                    }
                    T hres = Tzero();
                    const std::size_t hidx = lqn.parent[tidx];
                    if (hidx >= 1 && hidx < idxhash.size() && idxhash[hidx] >= 0 &&
                        ph_has_layer[hidx] && ph_layers[hidx].class_of_caller[tidx] > 0) {
                        const PHLayer& HL = ph_layers[hidx];
                        const LayerResult<T>& hr = results.back()[std::size_t(idxhash[hidx])];
                        const T rr = hr.RN(HL.qstations[0] - 1, HL.class_of_caller[tidx] - 1);
                        if (!std::isnan(dbl(rr))) hres = rr;
                    }
                    T za = T(num_traits<T>::from_double(nja / arvrate) - hres - ztask);
                    if (za < floorv) za = floorv;
                    if (relax_omega < 1.0 && it > 1 && !std::isnan(thinkt_prev[tidx])) {
                        const T om = num_traits<T>::from_double(relax_omega);
                        const T om1 = num_traits<T>::from_double(1.0 - relax_omega);
                        za = T(om * za + om1 * thinkt_prev_v[tidx]);
                    }
                    tput[tidx] = num_traits<T>::from_double(arvrate);
                    thinkt[tidx] = za;
                    thinkt_prev[tidx] = dbl(za);
                    thinkt_prev_v[tidx] = za;
                    thinktproc[tidx] = Distrib<T>::exp_mean(T(za + ztask));
                    continue;
                }
                // a reference task, or one no other task calls: it has no station
                // of its own, so its only delay is the think time the user declared
                thinkt[tidx] = num_traits<T>::from_double(GlobalConstants::FineTol);
                thinktproc[tidx] = Distrib<T>::immediate();
                continue;
            }
            const PHLayer& L = ph_layers[tidx];
            const qn::Layer<T>& m = ensemble[std::size_t(idxhash[tidx])];
            const LayerResult<T>& r = results.back()[std::size_t(idxhash[tidx])];
            T U = Tzero();
            for (std::size_t k = 0; k < m.nclasses; ++k) {
                const T v = r.UN(L.qstations[0] - 1, k);
                if (!std::isnan(dbl(v))) U = T(U + v);
            }
            util[tidx] = U;
            // The closure below measures a thread's cycle in WORK: it is idle for
            // whatever of the cycle its station does not hold it working. A
            // SetupTask's station service also carries a cold start, which is time
            // the thread is unavailable but is not work, so it is taken back out of
            // U before the closure reads it. Zero for every task without a setup.
            if (dbl(ph_setupshare[tidx]) > 0.0)
                U = T(U * (Tone() - ph_setupshare[tidx]));
            // the rate the CALLERS ask of the task, not the rate its processor
            // layer reported: the latter is itself a function of this think time
            T X = ph_xdemand[tidx];
            if (!(dbl(X) > GlobalConstants::FineTol)) X = tput[tidx];
            // The thread pool of ONE replica, the convention ph_xdemand is kept in
            double nj = lqn.maxmult[tidx];
            if (!std::isfinite(nj) || nj <= 0.0) {
                nj = 0.0;
                for (std::size_t c = 1; c <= NT(); ++c) nj = std::max(nj, njobs(tidx, c));
            }
            T z;
            if (dbl(X) > GlobalConstants::FineTol) {
                if (lqn.sched[tidx] == SchedStrategy::INF) {
                    // an infinite server reports a mean number of busy threads
                    z = T((num_traits<T>::from_double(nj) - U) / X - ztask);
                } else {
                    const T om = U > Tone() ? T(U - Tone()) : T(Tone() - U);
                    z = T(num_traits<T>::from_double(nj) * om / X - ztask);
                }
            } else {
                z = thinkt[tidx];
            }
            if (z < floorv) z = floorv;
            if (it > 1 && !std::isnan(thinkt_prev[tidx]) && !std::isfinite(dbl(z)))
                z = thinkt_prev_v[tidx];
            if (relax_omega < 1.0 && it > 1 && !std::isnan(thinkt_prev[tidx])) {
                const T om = num_traits<T>::from_double(relax_omega);
                const T om1 = num_traits<T>::from_double(1.0 - relax_omega);
                z = T(om * z + om1 * thinkt_prev_v[tidx]);
            }
            thinkt[tidx] = z;
            thinkt_prev[tidx] = dbl(z);
            thinkt_prev_v[tidx] = z;
            thinktproc[tidx] = Distrib<T>::exp_mean(T(z + ztask));
        }
    }

    /**
     * LQN-level results of method "srvn.ph".
     *
     * The layers report per caller task, so every entry, activity and call figure
     * is rebuilt from the converged fixed point rather than read off a class row,
     * in the layout aggregate() returns: QN carries the entry and task
     * utilizations, UN the processor utilizations, RN the response times and WN
     * the residence times.
     */
    LnSolution<T> aggregate_ph() {
        const std::size_t N = lqn.nidx;
        LnSolution<T> out;
        out.QN.assign(N + 1, Tzero());
        out.UN.assign(N + 1, Tzero());
        out.RN.assign(N + 1, Tzero());
        out.TN.assign(N + 1, Tzero());
        out.AN.assign(N + 1, Tzero());
        out.WN.assign(N + 1, Tzero());
        out.defined_Q.assign(N + 1, false);
        out.defined_U.assign(N + 1, false);
        out.defined_R.assign(N + 1, false);
        out.defined_T.assign(N + 1, false);
        out.defined_A.assign(N + 1, false);
        out.defined_W.assign(N + 1, false);
        out.iterations = iterations_done;
        out.converged = did_converge;

        std::vector<T> PN(N + 1, Tzero()), UT(N + 1, Tzero());
        std::vector<bool> hasPN(N + 1, false), hasUT(N + 1, false);

        for (std::size_t a = 1; a <= lqn.nacts; ++a) {
            const std::size_t aidx = lqn.ashift + a;
            const std::size_t tidx = lqn.parent[aidx];
            if (ignore[tidx]) continue;
            const std::size_t hidx = lqn.parent[tidx];
            out.TN[aidx] = tput[aidx];
            out.defined_T[aidx] = true;
            out.RN[aidx] = servt[aidx];
            out.defined_R[aidx] = true;
            UT[aidx] = T(tput[aidx] * servt[aidx]);
            hasUT[aidx] = true;
            // LINE scales the utilization of a queueing station into [0,1] whatever
            // its multiplicity, and reports a mean number of busy servers at an
            // infinite server: the processor share of an activity follows the same
            // convention
            const T hd = lqn.hostdem[aidx].disabled ? Tzero() : lqn.hostdem[aidx].mean;
            PN[aidx] = T(tput[aidx] * hd / num_traits<T>::from_double(ph_host_servers(hidx)));
            hasPN[aidx] = true;
            PN[hidx] = T(PN[hidx] + PN[aidx]);
            hasPN[hidx] = true;
        }

        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            const std::size_t tidx = lqn.parent[eidx];
            if (ignore[tidx]) continue;
            out.TN[eidx] = tput[eidx];
            out.defined_T[eidx] = true;
            out.RN[eidx] = servt[eidx];
            out.defined_R[eidx] = true;
            UT[eidx] = T(tput[eidx] * servt[eidx]);
            hasUT[eidx] = true;
            for (std::size_t aidx : lqn.actsof[eidx]) {
                PN[eidx] = T(PN[eidx] + PN[aidx]);
                hasPN[eidx] = true;
            }
            // ResidT is reported per visit to the TASK, not per execution of the
            // activity: an activity of this entry runs EXECS times per invocation,
            // and the entry takes SHARE of the task's invocations. RespT stays per
            // execution.
            if (ph_has_wf[eidx]) {
                const std::unordered_map<std::size_t, T>& ex = ph_execs[eidx];
                for (std::size_t aidx : lqn.actsof[eidx]) {
                    out.WN[aidx] = T(ph_share[eidx] * ex.at(aidx) * residt[aidx]);
                    out.defined_W[aidx] = true;
                }
            }
            UT[tidx] = T(UT[tidx] + UT[eidx]);
            hasUT[tidx] = true;
        }

        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (ignore[tidx]) continue;
            out.TN[tidx] = tput[tidx];
            out.defined_T[tidx] = true;
            T w = Tzero();
            bool anyw = false;
            for (std::size_t aidx : lqn.actsof[tidx]) {
                PN[tidx] = T(PN[tidx] + PN[aidx]);
                hasPN[tidx] = true;
                if (out.defined_W[aidx]) { w = T(w + out.WN[aidx]); anyw = true; }
            }
            if (anyw) { out.WN[tidx] = w; out.defined_W[tidx] = true; }
        }

        for (std::size_t hidx = 1; hidx <= lqn.nhosts; ++hidx)
            out.defined_T[hidx] = false;  // kept undefined for consistency with LQNS

        for (std::size_t idx = 1; idx <= N; ++idx) {
            out.QN[idx] = UT[idx];
            out.defined_Q[idx] = hasUT[idx];
            out.UN[idx] = PN[idx];
            out.defined_U[idx] = hasPN[idx];
            // Idle, not undefined -- the same rule aggregate() applies, and for the
            // same reason: an unreachable element reports zero for the measures its
            // kind HAS and leaves the ones it never has undefined, so that the
            // table's NaN mask survives a disconnected component.
            if (ignore[idx]) {
                out.UN[idx] = Tzero();       // every kind reports a utilization
                out.defined_U[idx] = true;
                out.defined_A[idx] = false;  // nothing reports an arrival rate on an LQN
                const bool host = lqn.type[idx] == LqnElement::HOST;
                const bool task = lqn.type[idx] == LqnElement::TASK;
                const bool entry = lqn.type[idx] == LqnElement::ENTRY;
                out.QN[idx] = Tzero();
                out.defined_Q[idx] = !host;
                out.RN[idx] = Tzero();
                out.defined_R[idx] = !host && !task;
                out.WN[idx] = Tzero();
                out.defined_W[idx] = !host && !entry;
                out.TN[idx] = Tzero();
                out.defined_T[idx] = !host;
            }
        }
        return out;
    }

    void update_metrics(int it) {
        if (is_ph_encoding()) {
            update_metrics_ph(it);
            return;
        }
        if (lnmethod == "moment3") {
            update_metrics_moment_based(it);
            return;
        }
        update_metrics_default(it);
    }

    void update_metrics_default(int it) {
        const std::size_t N = lqn.nidx;
        servt.assign(N + 1, Tzero());
        residt.assign(N + 1, Tzero());
        const int iter_min = std::min(30, int(std::ceil(opt.iter_max / 4.0)));
        const bool averaging = averagingstart >= 0 && it >= iter_min;
        const int wnd = averaging ? int(it - averagingstart + 1) : 1;

        for (const UpdRow& row : servt_map) {
            const std::size_t e = std::size_t(idxhash[row.idx]);
            const qn::Layer<T>& L = ensemble[e];
            const std::size_t k = row.cls - 1;
            std::size_t c = 0;
            for (std::size_t cc = 0; cc < L.nchains; ++cc)
                if (L.chains[cc][k]) c = cc;
            const std::size_t refclass_c = L.refclass[c];
            const std::size_t refstat_k = L.classes[k].refstat;

            T sv = Tzero(), rs = Tzero(), tp = Tzero();
            const T wT = T(Tone() / num_traits<T>::from_int(wnd));
            for (int w = 0; w < wnd; ++w) {
                const LayerResult<T>& r = results[results.size() - 1 - w][e];
                sv += r.RN(row.node - 1, k) * wT;
                const T TN_ref = (refclass_c > 0 && refstat_k > 0) ? r.TN(refstat_k - 1, refclass_c - 1) : Tzero();
                if (dbl(TN_ref) > GlobalConstants::FineTol)
                    rs += r.QN(row.node - 1, k) / TN_ref * wT;
                else
                    rs += r.WN(row.node - 1, k) * wT;
                tp += r.TN(row.node - 1, k) * wT;
            }
            servt[row.aidx] = sv;
            residt[row.aidx] = rs;
            tput[row.aidx] = tp;

            // an activity think time is in series with the host demand
            const Distrib<T>& at = lqn.actthink[row.aidx];
            if (!at.disabled && dbl(at.mean) > GlobalConstants::FineTol) {
                servt[row.aidx] = T(servt[row.aidx] + at.mean);
                residt[row.aidx] = T(residt[row.aidx] + at.mean);
            }

            // An activity of an async-only entry takes RN, the response per
            // visit, and not the visit-weighted residence: entry selection
            // routes the task to each of its entries with a share, and there is
            // no caller-side visit ratio here to divide that share back out
            // (the sync branch below does exactly that). updateMetricsDefault.m:63-77.
            if (async_only_activity(row.aidx)) residt[row.aidx] = servt[row.aidx];

            if (relax_omega < 1.0 && it > 1) {
                const T om = num_traits<T>::from_double(relax_omega);
                const T om1 = num_traits<T>::from_double(1.0 - relax_omega);
                if (!std::isnan(servt_prev[row.aidx]))
                    servt[row.aidx] = T(om * servt[row.aidx] + om1 * servt_prev_v[row.aidx]);
                if (!std::isnan(residt_prev[row.aidx]))
                    residt[row.aidx] = T(om * residt[row.aidx] + om1 * residt_prev_v[row.aidx]);
                if (!std::isnan(tput_prev[row.aidx]))
                    tput[row.aidx] = T(om * tput[row.aidx] + om1 * tput_prev_v[row.aidx]);
            }
            servt_prev[row.aidx] = dbl(servt[row.aidx]);
            residt_prev[row.aidx] = dbl(residt[row.aidx]);
            tput_prev[row.aidx] = dbl(tput[row.aidx]);
            servt_prev_v[row.aidx] = servt[row.aidx];
            residt_prev_v[row.aidx] = residt[row.aidx];
            tput_prev_v[row.aidx] = tput[row.aidx];

            if (servt[row.aidx] > Tzero() && dbl(servt[row.aidx]) <= 1e10)
                servtproc[row.aidx] = Distrib<T>::exp_mean(servt[row.aidx]);
            tputproc[row.aidx] = Distrib<T>::exp_rate(tput[row.aidx]);
        }

        // The phase split of servt, updateMetricsDefault.m:120-151. It is
        // recomputed from scratch each iteration because servt is; the
        // overtaking probability it feeds is computed later, once the entry
        // throughputs exist.
        if (has_phase2) {
            servt_ph1.assign(N + 1, Tzero());
            servt_ph2.assign(N + 1, Tzero());
            for (std::size_t a = 1; a <= lqn.nacts; ++a) {
                const std::size_t aidx = lqn.ashift + a;
                if (lqn.actphase[a] == 1)
                    servt_ph1[aidx] = servt[aidx];
                else
                    servt_ph2[aidx] = servt[aidx];
            }
            for (std::size_t e = 1; e <= lqn.nentries; ++e) {
                const std::size_t eidx = lqn.eshift + e;
                for (std::size_t aidx : lqn.actsof[eidx]) {
                    if (aidx <= lqn.ashift) continue;
                    const std::size_t a = aidx - lqn.ashift;
                    if (a > lqn.nacts) continue;
                    if (lqn.actphase[a] == 1)
                        servt_ph1[eidx] = T(servt_ph1[eidx] + servt_ph1[aidx]);
                    else
                        servt_ph2[eidx] = T(servt_ph2[eidx] + servt_ph2[aidx]);
                }
            }
        }

        // throughput of the activities that appear only as client-side classes
        for (const UpdRow& row : thinkt_map) {
            if (!tputproc[row.aidx].disabled) continue;
            const std::size_t e = std::size_t(idxhash[row.idx]);
            T tp = Tzero();
            const T wT = T(Tone() / num_traits<T>::from_int(wnd));
            for (int w = 0; w < wnd; ++w)
                tp += results[results.size() - 1 - w][e].TN(row.node - 1, row.cls - 1) * wT;
            tput[row.aidx] = tp;
            tputproc[row.aidx] = Distrib<T>::exp_rate(tp);
        }

        // call service and residence times
        callservt.assign(lqn.ncalls + 1, Tzero());
        callresidt.assign(lqn.ncalls + 1, Tzero());
        for (const UpdRow& row : call_map) {
            if (row.node <= 1) continue;  // a client-side call class contributes none
            const std::size_t e = std::size_t(idxhash[row.idx]);
            const LayerResult<T>& r = results.back()[e];
            callservt[row.aidx] =
                T(r.RN(row.node - 1, row.cls - 1) * lqn.callproc_mean[row.aidx]);
            // Normalise per chain-reference visit, as residt does. WN divides by the
            // class's own reference rate when the layer is open (an INF client task),
            // which is per-ENTRY visit, and the entry rescaling in
            // resolve_entry_service would then charge the call once per entry.
            {
                const qn::Layer<T>& L = ensemble[e];
                const std::size_t k = row.cls - 1;
                std::size_t c = 0;
                for (std::size_t cc = 0; cc < L.nchains; ++cc)
                    if (L.chains[cc][k]) c = cc;
                const std::size_t refclass_c = L.refclass[c];
                const std::size_t refstat_k = L.classes[k].refstat;
                const T TN_ref = (refclass_c > 0 && refstat_k > 0) ? r.TN(refstat_k - 1, refclass_c - 1) : Tzero();
                callresidt[row.aidx] = dbl(TN_ref) > GlobalConstants::FineTol
                                           ? T(r.QN(row.node - 1, k) / TN_ref)
                                           : r.WN(row.node - 1, k);
            }
            const T rw = region_wait(e, row.cls);
            if (rw > Tzero()) {
                callservt[row.aidx] = T(callservt[row.aidx] + rw);
                callresidt[row.aidx] = T(callresidt[row.aidx] + rw);
            }
            if (relax_omega < 1.0 && it > 1 && !std::isnan(callservt_prev[row.aidx])) {
                const T om = num_traits<T>::from_double(relax_omega);
                const T om1 = num_traits<T>::from_double(1.0 - relax_omega);
                callservt[row.aidx] =
                    T(om * callservt[row.aidx] + om1 * callservt_prev_v[row.aidx]);
            }
            callservt_prev[row.aidx] = dbl(callservt[row.aidx]);
            callservt_prev_v[row.aidx] = callservt[row.aidx];
            callresidt_prev[row.aidx] = dbl(callresidt[row.aidx]);
        }

        resolve_entry_service();

        // The overtaking correction, updateMetricsDefault.m:313-351. It runs
        // HERE and not with the split above because it needs the entry
        // throughput resolve_entry_service has just produced.
        //
        // servt keeps both phases -- the server IS busy through phase 2, so the
        // utilization is unchanged -- while residt becomes the CALLER's view:
        // phase 1 in full, phase 2 only when the caller is actually overtaken.
        if (has_phase2) {
            for (std::size_t e = 1; e <= lqn.nentries; ++e) {
                const std::size_t eidx = lqn.eshift + e;
                const std::size_t tidx = lqn.parent[eidx];
                if (!(dbl(servt_ph2[eidx]) > GlobalConstants::FineTol)) continue;
                if (lqn.isref[tidx] || !lqn.issynccaller.any_col(eidx)) {
                    residt[eidx] = servt[eidx];
                    continue;
                }
                T entry_tput = Tzero();
                if (dbl(tput[eidx]) > GlobalConstants::FineTol)
                    entry_tput = tput[eidx];
                else if (dbl(tput[tidx]) > GlobalConstants::FineTol)
                    entry_tput = tput[tidx];
                prOvertake[e] =
                    dbl(entry_tput) > GlobalConstants::FineTol
                        ? lqn_overtake_prob_markov(lqn, servt, callresidt, tput, eidx,
                                                   servt_ph2[eidx])
                        : Tzero();
                residt[eidx] = T(servt_ph1[eidx] + prOvertake[e] * servt_ph2[eidx]);
            }
        }

        for (const UpdRow& row : call_map) {
            if (row.node <= 1) continue;
            const std::size_t eidx = lqn.callpair_dst[row.aidx];
            if (servt[eidx] > Tzero()) servtproc[eidx] = Distrib<T>::exp_mean(servt[eidx]);
        }
        for (const UpdRow& row : call_map) {
            if (row.node <= 1) continue;
            const std::size_t eidx = lqn.callpair_dst[row.aidx];
            if (it == 1) {
                callservt[row.aidx] = servt[eidx];
                callservtproc[row.aidx] = servtproc[eidx];
            } else if (callservt[row.aidx] > Tzero()) {
                callservtproc[row.aidx] = Distrib<T>::exp_mean(callservt[row.aidx]);
            }
        }

        // What a synchronous CALLER waits for at a phase-2 target is residt,
        // not servt (updateMetricsDefault.m:384-400): the loop just above set it
        // from the target's full service, which would charge the caller for
        // phase 2 it never waits through.
        if (has_phase2) {
            for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
                if (lqn.calltype[cidx] != CallType::SYNC) continue;
                const std::size_t target = lqn.callpair_dst[cidx];
                if (target <= lqn.eshift || target > lqn.eshift + lqn.nentries) continue;
                if (!(dbl(servt_ph2[target]) > GlobalConstants::FineTol)) continue;
                const T eff = residt[target];  // servt_ph1 + prOvertake * servt_ph2
                if (!(eff > Tzero())) continue;
                const T w = T(eff * lqn.callproc_mean[cidx]);
                callservt[cidx] = w;
                callresidt[cidx] = w;
                callservtproc[cidx] = Distrib<T>::exp_mean(w);
            }
        }
    }

    // -----------------------------------------------------------------------
    // updateMetricsMomentBased, the `moment3` method
    // -----------------------------------------------------------------------

    /**
     * The moments of an empirical CDF, MATLAB's `EmpiricalCDF.getMoments`.
     *
     * READ EXACTLY AS THE REFERENCE READS IT, and the reading is not the
     * obvious one: the weight of a bin is the CDF INCREMENT and its abscissa is
     * the MIDPOINT of the two t values, so this is a midpoint quadrature of
     * integral x^k dF and not a sum over grid points. Reading the grid as a pmf
     * at the right endpoint instead inflates every moment; that exact mistake
     * cost a wrong entry mean in the Python port (see _kb/06-solver-catalog.md).
     */
    static void cdf_moments(const fluid::FluidPassage& c, double& m1, double& m2, double& m3) {
        m1 = m2 = m3 = 0.0;
        for (std::size_t i = 0; i + 1 < c.t.size(); ++i) {
            const double x = 0.5 * (c.t[i + 1] + c.t[i]);
            const double w = c.cdf[i + 1] - c.cdf[i];
            m1 += x * w;
            m2 += x * x * w;
            m3 += x * x * x * w;
        }
    }

    /** The per-layer response-time CDFs, computed once and cached. */
    const std::vector<std::vector<fluid::FluidPassage>>& layer_cdf(std::size_t e) {
        if (cdf_repo[e].empty()) {
            // THE LAYER HAS NO ROUTING MATRIX UNTIL IT IS ASKED FOR ONE, and the
            // fluid drift is built from `sn.rt` alone: without this the ODE has
            // no transitions, the state decays to zero, every passage reports the
            // degenerate curve at the origin and every entry ends up with an
            // empty convolution and a service time of zero. Same reason as in
            // solve_layer_fluid, and the failure is silent in both.
            ensemble[e].refresh_rt();
            try {
                cdf_repo[e] = detail::ln_fluid_cdf_respt(ensemble[e], opt.layer_fluid);
            } catch (const std::exception&) {
                // The reference falls back to the LAYER's own solver when the
                // fluid passage fails, and for the MVA-family layers that is
                // the base-class exponential law with the layer's mean -- so
                // the fallback here is that law over the last solved averages,
                // rather than an empty repo that silently collapses every
                // entry law to its bare host demand.
                cdf_repo[e] = layer_cdf_exp_fallback(e);
            }
        }
        return cdf_repo[e];
    }

    /** The base-class exponential CDF over layer e's last solved mean response
     * times, the reference's fallback route when the fluid passage fails. */
    std::vector<std::vector<fluid::FluidPassage>> layer_cdf_exp_fallback(std::size_t e) {
        const qn::NetworkStruct<T>& L = ensemble[e];
        std::vector<std::vector<fluid::FluidPassage>> out(
            L.nstations, std::vector<fluid::FluidPassage>(L.nclasses));
        if (results.empty() || e >= results.back().size()) return out;
        const Matrix<T>& RN = results.back()[e].RN;
        const std::size_t npts = 100;
        for (std::size_t i = 0; i < L.nstations; ++i) {
            if (L.stations[i].nodetype == qn::NodeType::Source) continue;
            for (std::size_t r = 0; r < L.nclasses; ++r) {
                if (L.disabled[i][r]) continue;
                const double rn = (i < static_cast<std::size_t>(RN.rows()) &&
                                   r < static_cast<std::size_t>(RN.cols()))
                                      ? dbl(RN(i, r))
                                      : 0.0;
                if (!(std::isfinite(rn) && rn > 0.0)) continue;
                fluid::FluidPassage& cell = out[i][r];
                cell.t.reserve(npts);
                cell.cdf.reserve(npts);
                for (std::size_t j = 0; j < npts; ++j) {
                    const double q =
                        0.001 + (0.999 - 0.001) * static_cast<double>(j) / (npts - 1);
                    cell.cdf.push_back(q);
                    cell.t.push_back(-std::log(1.0 - q) * rn);
                }
            }
        }
        return out;
    }

    /**
     * task_tput / entry_tput at the host layer of `eidx`.
     *
     * This is what renormalises a residence time from "one visit to the TASK",
     * which is how the layer reports it, to "one visit to the ENTRY", which is
     * what an entry metric means. It exists as a helper because applying it
     * TWICE is a real and silent failure mode: `moment3` once summed per-visit
     * response times and then applied this ratio on top, inflating every entry
     * service time by the entries-per-task ratio.
     *
     * `state` is 0 when the entry's layers are ignored (nothing is assigned at
     * all), 1 when it has no synchronous caller (no ratio exists) and 2 when
     * `ratio` is set.
     */
    int entry_visit_ratio(std::size_t eidx, T& ratio) const {
        const std::size_t tidx = lqn.parent[eidx];
        const std::size_t hidx = lqn.parent[tidx];
        if (ignore[tidx] || ignore[hidx] || idxhash[hidx] < 0) return 0;
        if (!lqn.issynccaller.any_col(eidx)) return 1;
        const std::size_t hl = std::size_t(idxhash[hidx]);
        const qn::Layer<T>& L = ensemble[hl];
        const LayerResult<T>& r = results.back()[hl];
        T task_tput = Tzero(), entry_tput = Tzero();
        for (const auto& kv : L.attr_tasks)
            if (kv.second == tidx) task_tput += r.TN(L.clientIdx - 1, kv.first - 1);
        for (const auto& kv : L.attr_entries)
            if (kv.second == eidx) entry_tput += r.TN(L.clientIdx - 1, kv.first - 1);
        const T floor = num_traits<T>::from_double(GlobalConstants::Zero);
        ratio = T(task_tput / (entry_tput > floor ? entry_tput : floor));
        return 2;
    }

    /** (I - servtmatrix)^-1, the reference's `inv(eye - servtmatrix)`. */
    Matrix<T> entry_service_resolvent() const {
        const std::size_t dim = lqn.nidx + lqn.ncalls;
        Matrix<T> A(dim + 1, dim + 1, Tzero());
        for (std::size_t i = 0; i <= dim; ++i) {
            A(i, i) = Tone();
            for (std::size_t j = 0; j <= dim; ++j) A(i, j) = T(A(i, j) - servtmatrix(i, j));
        }
        return ::line::inverse(A);
    }

    /**
     * Port of @@SolverLN/updateMetricsMomentBased.m, the `moment3` method.
     *
     * WHAT IT DOES DIFFERENTLY from the default update. The default feeds MEANS
     * between the layers. This fits an APH to the response-time CDF of every
     * activity and every call a layer reports, convolves those fits along the
     * entry's activity sequence, and reads the entry's law off the convolution:
     * a mean AND a distribution, which is what `get_cdf_respt` returns.
     *
     * TWO PASSES, not one. While the ensemble is still moving (`!hasconverged`)
     * the update is mean-based -- forming a CDF per layer per iteration would
     * cost a fluid integration per layer per iteration and would be fitting
     * noise anyway -- and the distribution is formed ONCE, on the converged
     * ensemble. The reference splits it exactly here.
     *
     * THE NORMALISATION TRAP. Both passes build the entry service from
     * RESIDENCE times, which are normalised to one visit to the TASK, and then
     * apply the task/entry throughput ratio, which converts that to one visit
     * to the ENTRY. Summing per-visit RESPONSE times and applying the ratio as
     * well applies the normalisation twice and inflates every multi-entry task
     * by its entries-per-task ratio. That was a live defect in all three
     * reference codebases until 2026-07-31.
     */
    void update_metrics_moment_based(int it) {
        const std::size_t N = lqn.nidx;
        servt.assign(N + 1, Tzero());
        residt.assign(N + 1, Tzero());
        callservt.assign(lqn.ncalls + 1, Tzero());
        callresidt.assign(lqn.ncalls + 1, Tzero());

        // ---- what every activity's layer reports, common to both passes ----
        for (const UpdRow& row : servt_map) {
            const std::size_t e = std::size_t(idxhash[row.idx]);
            const qn::Layer<T>& L = ensemble[e];
            const LayerResult<T>& r = results.back()[e];
            const std::size_t k = row.cls - 1;
            std::size_t c = 0;
            for (std::size_t cc = 0; cc < L.nchains; ++cc)
                if (L.chains[cc][k]) c = cc;
            const std::size_t refclass_c = L.refclass[c];
            const std::size_t refstat_k = L.classes[k].refstat;
            const T TN_ref = (refclass_c > 0 && refstat_k > 0) ? r.TN(refstat_k - 1, refclass_c - 1) : Tzero();

            tput[row.aidx] = r.TN(row.node - 1, k);
            residt[row.aidx] = dbl(TN_ref) > GlobalConstants::FineTol
                                   ? T(r.QN(row.node - 1, k) / TN_ref)
                                   : r.WN(row.node - 1, k);
            if (!hasconverged) {
                servt[row.aidx] = r.RN(row.node - 1, k);
                servtproc[row.aidx] = Distrib<T>::exp_mean(servt[row.aidx]);
                const Distrib<T>& at = lqn.actthink[row.aidx];
                if (!at.disabled && dbl(at.mean) > GlobalConstants::FineTol) {
                    servt[row.aidx] = T(servt[row.aidx] + at.mean);
                    residt[row.aidx] = T(residt[row.aidx] + at.mean);
                    servtproc[row.aidx] = Distrib<T>::exp_mean(servt[row.aidx]);
                }
                // An activity of an async-only entry takes the per-visit
                // response, for the same reason as in the default update.
                if (async_only_activity(row.aidx)) residt[row.aidx] = servt[row.aidx];
            } else {
                servtcdf[row.aidx] = layer_cdf(e)[row.node - 1][k];
            }
        }

        for (const UpdRow& row : call_map) {
            if (row.node <= 1) continue;  // a client-side call class contributes none
            const std::size_t e = std::size_t(idxhash[row.idx]);
            const LayerResult<T>& r = results.back()[e];
            callresidt[row.aidx] = r.WN(row.node - 1, row.cls - 1);
            if (!hasconverged)
                callservt[row.aidx] =
                    T(r.RN(row.node - 1, row.cls - 1) * lqn.callproc_mean[row.aidx]);
            else
                callservtcdf[row.aidx] = layer_cdf(e)[row.node - 1][row.cls - 1];
        }

        if (!hasconverged)
            moment3_means_pass(it);
        else
            moment3_distribution_pass(it);
    }

    /** The mean-based pass of moment3, run while the ensemble is still moving. */
    void moment3_means_pass(int it) {
        const std::size_t dim = lqn.nidx + lqn.ncalls;
        std::vector<T> x(dim + 1, Tzero());
        for (std::size_t i = 1; i <= lqn.nidx; ++i) x[i] = residt[i];
        for (std::size_t c = 1; c <= lqn.ncalls; ++c) x[lqn.nidx + c] = callresidt[c];

        // entry_servt = (I - servtmatrix) \ [residt; callresidt]
        const Matrix<T> Rinv = entry_service_resolvent();
        std::vector<T> entry_servt(dim + 1, Tzero());
        for (std::size_t i = 1; i <= dim; ++i) {
            T s = Tzero();
            for (std::size_t j = 1; j <= dim; ++j)
                if (Rinv(i, j) != Tzero()) s += Rinv(i, j) * x[j];
            entry_servt[i] = s;
        }
        for (std::size_t i = 1; i <= lqn.eshift; ++i) entry_servt[i] = Tzero();

        // NO forwarding propagation here. `lqn_fwd_rendezvous` has already
        // reconnected every forwarding chain reachable from a synchronous call to
        // the client that issued the rendezvous (Franks 1999, Sec. 3.3.1), so the
        // forwarded service is in the caller's chain before this runs; charging it
        // again inflated the caller by exactly the forwarded entry's mean. An
        // asynchronous call into a chain is left untouched there by design -- a
        // send-no-reply does not block -- so it must not accumulate the forwarded
        // service either. See BUGS.md BUG-91.

        for (std::size_t e = 1; e <= lqn.nentries; ++e)
            servt[lqn.eshift + e] = entry_servt[lqn.eshift + e];

        // entry_residt = servtmatrix * [residt; callresidt]
        std::vector<T> entry_residt(dim + 1, Tzero());
        for (std::size_t i = lqn.eshift + 1; i <= lqn.eshift + lqn.nentries; ++i) {
            T s = Tzero();
            for (std::size_t j = 1; j <= dim; ++j)
                if (servtmatrix(i, j) != Tzero()) s += servtmatrix(i, j) * x[j];
            entry_residt[i] = s;
        }

        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            T ratio = Tone();
            const int state = entry_visit_ratio(eidx, ratio);
            if (state == 0) continue;
            if (state == 1) {
                residt[eidx] = entry_residt[eidx];
                continue;
            }
            servt[eidx] = T(entry_servt[eidx] * ratio);
            residt[eidx] = T(entry_residt[eidx] * ratio);
        }

        for (const UpdRow& row : call_map) {
            if (row.node <= 1) continue;
            const std::size_t eidx = lqn.callpair_dst[row.aidx];
            if (servt[eidx] > Tzero()) servtproc[eidx] = Distrib<T>::exp_mean(servt[eidx]);
        }
        for (const UpdRow& row : call_map) {
            if (row.node <= 1) continue;
            const std::size_t eidx = lqn.callpair_dst[row.aidx];
            if (it == 1) {
                // A response time is per visit, so the number of calls is 1 here.
                callservt[row.aidx] = servt[eidx];
                callservtproc[row.aidx] = servtproc[eidx];
            } else if (callservt[row.aidx] > Tzero()) {
                callservtproc[row.aidx] = Distrib<T>::exp_mean(callservt[row.aidx]);
            }
        }
    }

    /**
     * The distribution pass of moment3, run once on the converged ensemble.
     *
     * Every term reachable from an entry is fitted to an APH from the first
     * three moments of its response-time CDF, repeated as many times as the
     * entry-service matrix says it is visited, and the whole sequence is
     * convolved. A FRACTIONAL repetition count -- a call made 1.5 times on
     * average -- is realised as a branch between the fitted law and a point
     * mass at zero, which is the reference's `aph_simplify(..., pattern 3)`.
     */
    void moment3_distribution_pass(int it) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "SolverLN: the 'moment3' method fits an APH to a response-time CDF, which needs "
                "square roots and a matrix exponential; rerun with --arith double or real");
        } else {
            const std::size_t dim = lqn.nidx + lqn.ncalls;
            const Matrix<T> Rinv = entry_service_resolvent();

            // The point mass at zero the fractional branch mixes against.
            mam::AphPair<T> zero_law;
            zero_law.alpha.push_back(Tone());
            zero_law.S = Matrix<T>(1, 1, num_traits<T>::from_double(-GlobalConstants::Immediate));

            for (std::size_t en = 1; en <= lqn.nentries; ++en) {
                const std::size_t eidx = lqn.eshift + en;
                std::vector<mam::AphPair<T>> seq;
                for (std::size_t fitidx = 1; fitidx <= dim; ++fitidx) {
                    if (!(Rinv(eidx, fitidx) > Tzero())) continue;
                    // Entries themselves are dropped: the entry's own service is
                    // what is being assembled, and a host or task index carries
                    // no response-time law of its own.
                    if (fitidx <= lqn.eshift + lqn.nentries) continue;
                    const bool is_call = fitidx > lqn.nidx;
                    const fluid::FluidPassage& curve =
                        is_call ? callservtcdf[fitidx - lqn.nidx] : servtcdf[fitidx];
                    double m1 = 0.0, m2 = 0.0, m3 = 0.0;
                    cdf_moments(curve, m1, m2, m3);

                    // An activity think time is in series with the host demand,
                    // so its raw moments convolve with the measured ones.
                    if (!is_call && !lqn.actthink[fitidx].disabled &&
                        dbl(lqn.actthink[fitidx].mean) > GlobalConstants::FineTol) {
                        const Distrib<T>& zt = lqn.actthink[fitidx];
                        const double t1 = dbl(lang::dist_moment(zt, 1));
                        const double t2 = dbl(lang::dist_moment(zt, 2));
                        const double t3 = dbl(lang::dist_moment(zt, 3));
                        m3 = m3 + 3.0 * m2 * t1 + 3.0 * m1 * t2 + t3;
                        m2 = m2 + 2.0 * m1 * t1 + t2;
                        m1 = m1 + t1;
                    }

                    // CoarseTol and not FineTol: an Immediate activity has a
                    // near-zero mean whose APH fit has rates of order 1e8, and
                    // the matrix exponential of the convolution then does not
                    // terminate. The reference skips those terms outright.
                    if (!(m1 > GlobalConstants::CoarseTol)) continue;

                    const mam::AphFitResult<T> fit = mam::aph_fit(
                        num_traits<T>::from_double(m1), num_traits<T>::from_double(m2),
                        num_traits<T>::from_double(m3), 10u,
                        num_traits<T>::from_double(GlobalConstants::FineTol));
                    const mam::AphPair<T> law = aph_pair_of_map(fit.aph);

                    double reps = dbl(Rinv(eidx, fitidx));
                    // servtmatrix carries 1.0 for a call; the mean NUMBER of
                    // calls is what says how many times its law is convolved.
                    if (is_call) reps *= dbl(lqn.callproc_mean[fitidx - lqn.nidx]);
                    const long whole = static_cast<long>(std::floor(reps));
                    const double frac = reps - static_cast<double>(whole);
                    for (long q = 0; q < whole; ++q) seq.push_back(law);
                    if (frac > 0.0)
                        seq.push_back(mam::aph_simplify(law, zero_law,
                                                        num_traits<T>::from_double(frac),
                                                        num_traits<T>::from_double(1.0 - frac),
                                                        mam::AphPattern::Branch));

                    if (is_call) {
                        const std::size_t cidx = fitidx - lqn.nidx;
                        callservt[cidx] = num_traits<T>::from_double(m1);
                        callservtproc[cidx] = Distrib<T>::exp_mean(callservt[cidx]);
                    } else {
                        servt[fitidx] = num_traits<T>::from_double(m1);
                        servtproc[fitidx] = Distrib<T>::exp_mean(servt[fitidx]);
                    }
                }

                if (seq.empty()) {
                    servt[eidx] = Tzero();
                    continue;
                }
                const mam::AphPair<T> entry_law = mam::aph_convseq(seq);
                entryproc[en] = entry_law;
                servt[eidx] = Distrib<T>::ph_moment(entry_law.alpha, entry_law.S, 1);
                servtproc[eidx] = Distrib<T>::exp_mean(servt[eidx]);
                entrycdfrespt[en] = aph_eval_cdf(entry_law);
            }

            // NO forwarding propagation here, for the reason given at the
            // entry_servt assembly above. See BUGS.md BUG-91.

            // entry_residt = servtmatrix * [residt; callresidt], then the same
            // task/entry renormalisation the mean pass applies.
            std::vector<T> x(dim + 1, Tzero());
            for (std::size_t i = 1; i <= lqn.nidx; ++i) x[i] = residt[i];
            for (std::size_t c = 1; c <= lqn.ncalls; ++c) x[lqn.nidx + c] = callresidt[c];
            for (std::size_t en = 1; en <= lqn.nentries; ++en) {
                const std::size_t eidx = lqn.eshift + en;
                T s = Tzero();
                for (std::size_t j = 1; j <= dim; ++j)
                    if (servtmatrix(eidx, j) != Tzero()) s += servtmatrix(eidx, j) * x[j];
                T ratio = Tone();
                const int state = entry_visit_ratio(eidx, ratio);
                if (state == 0) continue;
                residt[eidx] = state == 2 ? T(s * ratio) : s;
            }

            if (it == 1)
                for (const UpdRow& row : call_map) {
                    if (row.node <= 1) continue;
                    const std::size_t eidx = lqn.callpair_dst[row.aidx];
                    callservt[row.aidx] = servt[eidx];
                    callservtproc[row.aidx] = Distrib<T>::exp_mean(servt[eidx]);
                }

            // The entry servt is now the MEAN OF AN APH CONVOLUTION, not a sum
            // of residence times, which the interlock rescale in
            // `update_populations` has to know. See BUGS.md BUG-97.
            moment_pass_done = true;
        }
    }

    /**
     * (alpha, S) of an APH handed back by aph_fit as a (D0, D1) pair.
     *
     * D1(i,j) = (-D0 e)_i alpha_j by construction, so alpha is any row of D1
     * divided by that row's exit rate; the first row with a positive exit rate
     * is taken. It is not recovered from the stationary phase distribution,
     * which is a different vector and would silently refit the law.
     */
    static mam::AphPair<T> aph_pair_of_map(const mam::Map<T>& m) {
        mam::AphPair<T> out;
        const std::size_t n = m.D0.rows();
        out.S = m.D0;
        out.alpha.assign(n, Tzero());
        for (std::size_t i = 0; i < n; ++i) {
            T rowsum = Tzero();
            for (std::size_t j = 0; j < n; ++j) rowsum += m.D1(i, j);
            if (!(dbl(rowsum) > GlobalConstants::Zero)) continue;
            for (std::size_t j = 0; j < n; ++j) out.alpha[j] = T(m.D1(i, j) / rowsum);
            return out;
        }
        // No phase can absorb: the law is degenerate, so it enters phase one and
        // stays there, which is what a zero alpha would NOT say.
        out.alpha[0] = Tone();
        return out;
    }

    /**
     * F(t) = 1 - alpha exp(S t) e on the reference's grid, `APH.evalCDF` with
     * no argument: 500 points over [0, mean + 10 sigma].
     */
    static LnCdf aph_eval_cdf(const mam::AphPair<T>& law) {
        LnCdf out;
        const double m1 = dbl(Distrib<T>::ph_moment(law.alpha, law.S, 1));
        const double m2 = dbl(Distrib<T>::ph_moment(law.alpha, law.S, 2));
        const double var = m2 - m1 * m1;
        const double sigma = var > 0.0 ? std::sqrt(var) : 0.0;
        const double tmax = m1 + 10.0 * sigma;
        const std::size_t P = 500;
        out.t.resize(P);
        out.cdf.resize(P);
        for (std::size_t k = 0; k < P; ++k) {
            const double t = tmax * static_cast<double>(k) / static_cast<double>(P - 1);
            out.t[k] = t;
            Matrix<T> St(law.S.rows(), law.S.cols(), Tzero());
            for (std::size_t i = 0; i < law.S.rows(); ++i)
                for (std::size_t j = 0; j < law.S.cols(); ++j)
                    St(i, j) = T(law.S(i, j) * num_traits<T>::from_double(t));
            const Matrix<T> E = ::line::expm(St);
            T surv = Tzero();
            for (std::size_t i = 0; i < E.rows(); ++i)
                for (std::size_t j = 0; j < E.cols(); ++j) surv += law.alpha[i] * E(i, j);
            out.cdf[k] = 1.0 - dbl(surv);
        }
        return out;
    }

    /**
     * The AND-join completion times, and how much they undercut the serial sum.
     *
     * The branches of an AND fork run concurrently, so the time to clear the
     * join is the k-th smallest of the branch completion times, k being the
     * quorum. The entry-service reachability matrix cannot express that -- it
     * charges every activity of every branch to the entry, i.e. it serialises
     * them -- so the difference is recorded per join and applied as a
     * correction to any entry that reaches it.
     *
     * Branch times are taken as exponential, so the variance is the square of
     * the mean; that is the reference's assumption, not an approximation added
     * here.
     */
    std::vector<T> join_excess() const {
        std::vector<T> excess(lqn.nidx + 1, Tzero());
        std::vector<std::size_t> joined;
        for (std::size_t tail = 1; tail <= lqn.nidx; ++tail) {
            if (lqn.actpretype[tail] != PrecedenceType::PRE_AND) continue;
            for (std::size_t sx : lqn.graph.succ(tail))
                if (sx > lqn.ashift && sx <= lqn.ashift + lqn.nacts) joined.push_back(sx);
        }
        std::sort(joined.begin(), joined.end());
        joined.erase(std::unique(joined.begin(), joined.end()), joined.end());
        if (joined.empty()) return excess;

        fj::LqnBranchView<T> view;
        view.graph = Matrix<T>(lqn.nidx, lqn.nidx, Tzero());
        for (std::size_t i = 1; i <= lqn.nidx; ++i)
            for (const auto& e : lqn.graph.row[i]) view.graph(i - 1, e.first - 1) = e.second;
        view.ashift = lqn.ashift;
        view.nacts = lqn.nacts;
        view.actposttype.assign(lqn.nidx, 0);
        for (std::size_t i = 1; i <= lqn.nidx; ++i)
            view.actposttype[i - 1] = static_cast<int>(lqn.actposttype[i]);

        for (std::size_t aidx : joined) {
            const std::vector<std::vector<std::size_t>> branches =
                fj::fj_branch_members(view, aidx);
            if (branches.empty()) continue;
            std::vector<T> means;
            for (const auto& br : branches) {
                T s = Tzero();
                for (std::size_t a : br) {
                    s += residt[a];
                    // A branch activity with an Immediate host demand does all its work
                    // in a rendezvous, so residt alone leaves this correction inert.
                    if (a < lqn.callsof.size())
                        for (std::size_t cidx : lqn.callsof[a])
                            if (cidx < lqn.calltype.size() && lqn.calltype[cidx] == CallType::SYNC
                                && cidx < callresidt.size())
                                s += callresidt[cidx];
                }
                means.push_back(s);
            }
            if (means.size() == 1) continue;  // a single branch cannot overlap
            std::size_t quorum = means.size();
            if (lqn.actquorum[aidx] >= 1 && lqn.actquorum[aidx] <= means.size())
                quorum = lqn.actquorum[aidx];
            std::vector<T> vars;
            for (const T& m2 : means) vars.push_back(T(m2 * m2));
            T serial = Tzero();
            for (const T& m2 : means) serial += m2;
            if constexpr (num_traits<T>::has_transcendental) {
                const fj::FJQuorumMomentsResult<T> q = fj::fj_quorum_moments(means, vars, quorum);
                excess[aidx] = T(q.m - serial);
            } else {
                throw UnsupportedError(
                    "SolverLN: the AND-join completion time is a k-th order statistic fitted "
                    "through a three-point distribution, which needs a square root; use the "
                    "double or real backend for a model with an AND join");
            }
        }
        return excess;
    }

    /** The block that turns activity residence times into entry service times. */
    void resolve_entry_service() {
        const std::size_t dim = lqn.nidx + lqn.ncalls;
        std::vector<T> x(dim + 1, Tzero());
        for (std::size_t i = 1; i <= lqn.nidx; ++i) x[i] = residt[i];
        for (std::size_t c = 1; c <= lqn.ncalls; ++c) x[lqn.nidx + c] = callresidt[c];
        std::vector<T> entry_servt(dim + 1, Tzero());
        for (std::size_t i = 1; i <= dim; ++i) {
            T s = Tzero();
            for (std::size_t j = 1; j <= dim; ++j)
                if (servtmatrix(i, j) != Tzero()) s += servtmatrix(i, j) * x[j];
            entry_servt[i] = s;
        }

        // replace each reachable join's serialised branch sum by its concurrent
        // completion time
        const std::vector<T> excess = join_excess();
        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            for (std::size_t aidx = 1; aidx <= lqn.nidx; ++aidx) {
                if (excess[aidx] == Tzero()) continue;
                if (servtmatrix(eidx, aidx) > Tzero())
                    entry_servt[eidx] = T(entry_servt[eidx] + excess[aidx]);
            }
            if (entry_servt[eidx] < Tzero()) entry_servt[eidx] = Tzero();
        }

        // A SetupTask's cold start is charged HERE, to the entry, and with the
        // probability that the thread was actually found powered down. It is not
        // host demand, so it belongs to no activity's residence -- reporting it
        // there put RespT(A2) at 1.29479 against 0.333178 from LDES on lqn_setup,
        // the bare demand. The probability is the one 'srvn.ph' uses through
        // ph_setup_prob, so the two encodings charge the same thing
        // (updateMetricsDefault.m, lqn_setup_charge.m).
        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            const double c = setup_charge(lqn.parent[eidx]);
            if (c != 0.0)
                entry_servt[eidx] = T(entry_servt[eidx] + num_traits<T>::from_double(c));
        }

        // ResidT is normalised so the TASK has one visit; the entries need one
        // visit each, hence the throughput ratio below
        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            const std::size_t tidx = lqn.parent[eidx];
            const std::size_t hidx = lqn.parent[tidx];
            if (ignore[tidx] || ignore[hidx]) continue;
            if (idxhash[hidx] < 0) continue;
            const bool has_sync = lqn.issynccaller.any_col(eidx);
            if (!has_sync) {
                servt[eidx] = entry_servt[eidx];
                residt[eidx] = entry_servt[eidx];
                continue;
            }
            const std::size_t hl = std::size_t(idxhash[hidx]);
            const qn::Layer<T>& L = ensemble[hl];
            const LayerResult<T>& r = results.back()[hl];
            T task_tput = Tzero(), entry_tput = Tzero();
            for (const auto& kv : L.attr_tasks)
                if (kv.second == tidx) task_tput += r.TN(L.clientIdx - 1, kv.first - 1);
            for (const auto& kv : L.attr_entries)
                if (kv.second == eidx) entry_tput += r.TN(L.clientIdx - 1, kv.first - 1);
            if (dbl(entry_tput) > GlobalConstants::Zero) {
                servt[eidx] = T(entry_servt[eidx] * task_tput / entry_tput);
                residt[eidx] = servt[eidx];
            } else {
                servt[eidx] = entry_servt[eidx];
                residt[eidx] = entry_servt[eidx];
            }
        }
    }

    // -----------------------------------------------------------------------
    // updatePopulations (interlock correction)
    // -----------------------------------------------------------------------
    /**
     * True when the layer engine applies Eq. (4.7) inside its own MVA. A layer whose MVA path
     * has no interlock term would be moved to another algorithm by the matrix alone: exact
     * multiserver MVA would become AMVA, the product-form kernels would become the
     * load-dependent forward step. That swap is worth far more than the correction it carries,
     * and on a layer sitting near a bifurcation it turns the LN iteration into a limit cycle.
     * Such a layer keeps the residt scaling instead.
     */
    bool layer_takes_interlock(std::size_t e) const {
        return opt.layer_solver == "mva" && e < ensemble.size() &&
               mva::mva_carries_interlock(ensemble[e], opt.layer);
    }

    /**
     * Class-level interlock matrix of one host layer. IL[r][s] is the share of the class-s
     * queue that a class-r arrival must not see at the host. It is kept CLASS-indexed, not
     * chain-indexed, so that a later chain refresh cannot leave it stale: the analyzer
     * aggregates it to chains against the struct it is about to solve. Two classes are
     * interlocked only if BOTH their tasks are, which is the 0/1 relation ir_mkj of Eq. (5);
     * the diagonal stays zero, since a request always sees its own class in full. The entry
     * is the Eq. (5) product Pr(IL_ms)*IR_ms*IR_mr, asymmetric in (r,s) because Pr(IL) is
     * taken from the QUEUED class s, so that the layer's ILw(r,s) = 1-IL(r,s) is the
     * lower-level adjustment rate r_lower. An empty result keeps the layer on the plain MVA
     * path.
     */
    std::vector<std::vector<double>> build_layer_interlock(std::size_t e,
                                                           const std::vector<std::size_t>& hts,
                                                           const std::vector<T>& prIL,
                                                           const std::vector<T>& PrIL) const {
        const qn::Layer<T>& L = ensemble[e];
        const std::size_t R = L.nclasses;
        std::vector<double> cls_ir(R, 0.0);   // IR
        std::vector<double> cls_pr(R, 0.0);   // Pr(IL)
        auto stamp = [&](std::size_t classIdx, std::size_t tidx) {
            if (classIdx < 1 || classIdx > R) return;
            for (std::size_t i = 0; i < hts.size(); ++i)
                if (hts[i] == tidx) {
                    cls_ir[classIdx - 1] = dbl(prIL[i]);
                    cls_pr[classIdx - 1] = dbl(PrIL[i]);
                }
        };
        for (const auto& a : L.attr_tasks) stamp(a.first, a.second);
        for (const auto& a : L.attr_entries) stamp(a.first, lqn.parent[a.second]);
        for (const auto& a : L.attr_activities) stamp(a.first, lqn.parent[a.second]);
        for (const auto& a : L.attr_calls) stamp(a[0], lqn.parent[a[2]]);

        std::vector<std::vector<double>> IL(R, std::vector<double>(R, 0.0));
        bool any = false;
        for (std::size_t r = 0; r < R; ++r) {
            if (!(cls_ir[r] > GlobalConstants::FineTol)) continue;
            for (std::size_t sIl = 0; sIl < R; ++sIl) {
                if (sIl == r || !(cls_ir[sIl] > GlobalConstants::FineTol)) continue;
                IL[r][sIl] = cls_pr[sIl] * cls_ir[sIl] * cls_ir[r];
                if (IL[r][sIl] > GlobalConstants::FineTol) any = true;
            }
        }
        return any ? IL : std::vector<std::vector<double>>();
    }

    /**
     * Interlock probability for one (client, server) pair, as {IR, Pr(IL)} of Li and Franks,
     * "An improved interlocking correction for decomposition of layered queueing networks",
     * CCECE 2015, Eqs. (3) and (4). isProcessorHost selects the m' rule of lqns
     * Interlock::ilrate_pril_flow: at a PROCESSOR the common-source population is doubled
     * above 3 customers and squared at or below it, which is what turns m = 4 into the
     * pril = 1/8 its trace reports. The two factors are multiplied into the Eq. (5) rate by
     * build_layer_interlock, so neither carries the source count on its own -- that lives in
     * m'. This replaces the superseded (n_s-1)/n_s discount of Franks (1999), Eq. (4.7).
     */
    std::pair<T, T> interlock_prob(std::size_t client_tidx, std::size_t server_idx,
                                   bool isProcessorHost) const {
        const std::vector<std::size_t>& common = il_common_entries[server_idx];
        const double nsrc = il_num_sources[server_idx];
        if (nsrc == 0.0 || common.empty()) return std::make_pair(Tzero(), Tzero());
        const std::vector<std::size_t>& allsrc = il_src_all[server_idx];
        T sum_flow = Tzero();
        T sum_pril = Tzero();
        for (std::size_t ce : common) {
            const std::size_t srcTask = lqn.parent[ce];
            const std::size_t cen = ce - lqn.eshift;
            // population of this common source, in customer copies
            double m_src = lqn.mult[srcTask];
            if (!std::isfinite(m_src) || m_src < 1.0) m_src = 1.0;
            double m_eff = m_src;
            if (isProcessorHost) m_eff = (m_src > 3.0) ? (m_src + m_src) : (m_src * m_src);
            for (std::size_t dst : lqn.entriesof[client_tidx]) {
                const std::size_t dn = dst - lqn.eshift;
                if (dn < 1 || dn > lqn.nentries) continue;
                if (!(il_all(cen, dn) > Tzero())) continue;
                T ce_tput = tput[ce];
                if (!(dbl(ce_tput) > GlobalConstants::FineTol) && !lqn.actsof[ce].empty())
                    ce_tput = tput[lqn.actsof[ce][0]];
                if (!(dbl(ce_tput) > GlobalConstants::FineTol)) ce_tput = tput[srcTask];
                if (!(dbl(ce_tput) > GlobalConstants::FineTol)) continue;
                // phase-2 entries are refused upstream, so every source is all-phase
                if (std::find(allsrc.begin(), allsrc.end(), srcTask) != allsrc.end()) {
                    const T contrib = T(ce_tput * il_all(cen, dn));
                    sum_flow += contrib;
                    sum_pril += T(contrib / num_traits<T>::from_double(m_eff));
                }
            }
        }
        T client_tput = tput[client_tidx];
        if (!(dbl(client_tput) > GlobalConstants::FineTol)) {
            for (std::size_t e : lqn.entriesof[client_tidx]) {
                T et = tput[e];
                if (!(dbl(et) > GlobalConstants::FineTol) && !lqn.actsof[e].empty())
                    et = tput[lqn.actsof[e][0]];
                client_tput = T(client_tput + et);
            }
        }
        if (!(dbl(client_tput) > GlobalConstants::FineTol))
            return std::make_pair(Tzero(), Tzero());
        const T flow = sum_flow < client_tput ? sum_flow : client_tput;
        T IR = T(flow / client_tput);
        if (IR > Tone()) IR = Tone();
        if (dbl(IR) < 0.0) IR = Tzero();
        T pr = Tzero();
        if (dbl(sum_flow) > GlobalConstants::FineTol) pr = T(sum_pril / sum_flow);
        if (pr > Tone()) pr = Tone();
        if (dbl(pr) < 0.0) pr = Tzero();
        return std::make_pair(IR, pr);
    }

    void update_populations(int) {
        const std::vector<T> residt_orig = residt;
        const std::vector<T> callresidt_orig = callresidt;
        bool adjusted = false;

        for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
            if (lqn.calltype[cidx] != CallType::SYNC) continue;
            const std::size_t dst = lqn.callpair_dst[cidx];
            const std::size_t server_tidx = lqn.parent[dst];
            std::size_t server = 0;
            if (server_tidx <= NT() && !il_common_entries[server_tidx].empty()) {
                server = server_tidx;
            } else if (server_tidx > lqn.tshift) {
                const std::size_t h = lqn.parent[server_tidx];
                if (h >= 1 && h <= NT() && !il_common_entries[h].empty()) server = h;
            }
            if (server == 0) continue;
            const std::size_t client_tidx = lqn.parent[lqn.callpair_src[cidx]];
            // This path serves a TASK, not a processor, so the m' rule of Li/lqns leaves the
            // source population alone; the product IR*Pr(IL) reproduces the scalar this
            // branch used before.
            const std::pair<T, T> ilp = interlock_prob(client_tidx, server, false);
            const T prIL = T(ilp.first * ilp.second);
            if (!(dbl(prIL) > GlobalConstants::FineTol)) continue;
            const T S = servt[dst];
            const T cm = lqn.callproc_mean[cidx];
            if (!(cm > Tzero()) || !(callservt[cidx] > Tzero())) continue;
            const T RN = T(callservt[cidx] / cm);
            const T W = RN > S ? T(RN - S) : Tzero();
            if (!(dbl(W) > GlobalConstants::FineTol)) continue;
            const T RN_adj = T(S + (Tone() - prIL) * W);
            const T scale = T(RN_adj / RN);
            callservt[cidx] = T(callservt[cidx] * scale);
            callresidt[cidx] = T(callresidt[cidx] * scale);
            if (callservt[cidx] > Tzero())
                callservtproc[cidx] = Distrib<T>::exp_mean(callservt[cidx]);
            adjusted = true;
        }

        // Every layer starts the pass without a matrix, so a host that stops being
        // interlocked does not keep the previous iteration's correction alive.
        layer_interlock.assign(ensemble.size(), {});
        for (std::size_t h = 1; h <= lqn.nhosts; ++h) {
            if (il_common_entries[h].empty()) continue;
            const std::vector<std::size_t>& hts = lqn.tasksof[h];
            std::vector<T> prIL(hts.size(), Tzero()), PrIL(hts.size(), Tzero()),
                tutil(hts.size(), Tzero());
            for (std::size_t i = 0; i < hts.size(); ++i) {
                // The host of a task layer is a PROCESSOR, which is what selects the m' rule.
                const std::pair<T, T> ilp = interlock_prob(hts[i], h, true);
                prIL[i] = ilp.first;
                PrIL[i] = ilp.second;
                for (std::size_t e : lqn.entriesof[hts[i]])
                    for (std::size_t a : lqn.actsof[e]) tutil[i] += tput[a] * lqn.hostdem[a].mean;
            }
            T Utot = Tzero(), Uil = Tzero();
            for (std::size_t i = 0; i < hts.size(); ++i) {
                Utot += tutil[i];
                if (dbl(prIL[i]) > GlobalConstants::FineTol) Uil += tutil[i];
            }
            if (!(dbl(Utot) > GlobalConstants::FineTol) || !(dbl(Uil) > GlobalConstants::FineTol))
                continue;
            const T frac = T(Uil / Utot);

            // When the layer engine carries Eq. (4.7) inside its own MVA, the interlock goes
            // to the layer as a class-level matrix and the residence times are left
            // untouched. Scaling them here as well would remove the same waiting twice, and
            // would still leave the layer's own THROUGHPUT uncorrected, which is what breaks
            // flow balance across a call: the reported task rate then comes from a cycle time
            // the correction has already shortened elsewhere.
            if (idxhash[h] >= 0 && layer_takes_interlock(std::size_t(idxhash[h]))) {
                const std::size_t e = std::size_t(idxhash[h]);
                layer_interlock[e] = build_layer_interlock(e, hts, prIL, PrIL);
                continue;
            }

            for (std::size_t i = 0; i < hts.size(); ++i) {
                if (!(dbl(prIL[i]) > GlobalConstants::FineTol)) continue;
                // Weight by the share of host utilization that is interlocked. The rate
                // is the SAME Eq. (5) product IR*Pr(IL) that pass 1 applies to a call and
                // that build_layer_interlock puts in the layer matrix -- IR alone is a
                // flow SHARE, ~1 whenever a layer has a single common source, and using
                // it here removed the whole processor queueing rather than the
                // interlocked part of it, which broke flow balance across a call.
                const T eff = T(prIL[i] * PrIL[i] * frac);
                for (std::size_t e : lqn.entriesof[hts[i]])
                    for (std::size_t a : lqn.actsof[e]) {
                        const T D = lqn.hostdem[a].mean;
                        if (D > Tzero() && dbl(residt[a] - D) > GlobalConstants::FineTol) {
                            residt[a] = T(D + (Tone() - eff) * (residt[a] - D));
                            adjusted = true;
                        }
                    }
            }
        }

        if (!adjusted) return;
        const std::size_t dim = lqn.nidx + lqn.ncalls;
        auto esum = [&](const std::vector<T>& rs, const std::vector<T>& cr, std::size_t i) {
            T s = Tzero();
            for (std::size_t j = 1; j <= dim; ++j) {
                if (servtmatrix(i, j) == Tzero()) continue;
                s += servtmatrix(i, j) * (j <= lqn.nidx ? rs[j] : cr[j - lqn.nidx]);
            }
            return s;
        };
        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            const T oldv = esum(residt_orig, callresidt_orig, eidx);
            if (!(dbl(oldv) > GlobalConstants::FineTol)) continue;
            const T ratio = T(esum(residt, callresidt, eidx) / oldv);
            // The entry servt is rescaled only when it was itself assembled from
            // these residence times, which is the mean-based path. After the
            // moment3 distribution pass it is the mean of an APH convolution of
            // the activities' own response laws, and a ratio of residence-time
            // sums is not a correction to it: applying it multiplied the entry
            // law by the entry's visit ratio and reported a service time BELOW
            // that of the single activity the entry contains. The residence
            // times keep their correction either way. See BUGS.md BUG-97.
            if (!moment_pass_done) {
                servt[eidx] = T(servt[eidx] * ratio);
                if (servt[eidx] > Tzero()) servtproc[eidx] = Distrib<T>::exp_mean(servt[eidx]);
            }
            residt[eidx] = T(residt[eidx] * ratio);
        }
    }

    // -----------------------------------------------------------------------
    // Declared think time of a task as it enters the thread cycle: the value for
    // a REFERENCE task, zero for any other.
    //
    // A think time is an attribute of the closed customer population a reference
    // task stands for, and it is what separates one request of that population
    // from the next. On a served task it has no such meaning, and charging it
    // per request throttles the task: lqn_basic's T3 has 25 threads and a
    // declared think time of 4, and reading it as a per-request delay caps it at
    // 25/(4+0.02) = 6.219 completions per second. Three independent oracles put
    // the rate at five calls per caller request instead -- lqsim 66.5, LDES
    // 66.955, lqns 75.6. See _kb/06-solver-catalog.md (LN section).
    T ref_think_time(std::size_t tidx) const {
        if (tidx >= lqn.isref.size() || !lqn.isref[tidx]) return Tzero();
        if (lqn.think[tidx].disabled) return Tzero();
        return lqn.think[tidx].mean;
    }

    // updateThinkTimes
    // -----------------------------------------------------------------------
    void update_think_times(int it) {
        // Under "srvn.ph" a caller reaches the server once per invocation, so the
        // station rate is not the task's invocation rate
        if (is_ph_encoding()) {
            update_think_times_ph(it);
            return;
        }
        thinktproc.assign(lqn.nidx + 1, Distrib<T>::disabled_dist());
        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            // only a REFERENCE task's think time is a per-request delay
            const T ztask = ref_think_time(tidx);
            if (idxhash[tidx] < 0) {
                // A task reached only by an entry arrival still has a cycle: its threads
                // are driven by the stream. build_layer drops the open class for it so
                // the chain can be closed on the known rate here.
                const double arvrate = open_arrival_rate_of(tidx);
                if (arvrate > GlobalConstants::FineTol) {
                    double nja = 0.0;
                    for (std::size_t c = 1; c <= NT(); ++c) nja = std::max(nja, njobs(tidx, c));
                    if (!(nja > 0.0)) nja = lqn.maxmult[tidx];
                    T hres = Tzero();
                    for (std::size_t eidx : lqn.entriesof[tidx])
                        for (std::size_t aidx : lqn.actsof[eidx])
                            if (!std::isnan(dbl(residt[aidx]))) hres += residt[aidx];
                    const T floora = num_traits<T>::from_double(GlobalConstants::Zero);
                    T za = T(num_traits<T>::from_double(nja / arvrate) - hres - ztask);
                    if (za < floora) za = floora;
                    if (relax_omega < 1.0 && it > 1 && !std::isnan(thinkt_prev[tidx])) {
                        const T om = num_traits<T>::from_double(relax_omega);
                        const T om1 = num_traits<T>::from_double(1.0 - relax_omega);
                        za = T(om * za + om1 * thinkt_prev_v[tidx]);
                    }
                    tput[tidx] = num_traits<T>::from_double(arvrate);
                    thinkt[tidx] = za;
                    thinkt_prev[tidx] = dbl(za);
                    thinkt_prev_v[tidx] = za;
                    thinktproc[tidx] = Distrib<T>::exp_mean(T(za + ztask));
                    continue;
                }
                thinkt[tidx] = num_traits<T>::from_double(GlobalConstants::FineTol);
                thinktproc[tidx] = Distrib<T>::immediate();
                continue;
            }
            const std::size_t e = std::size_t(idxhash[tidx]);
            const qn::Layer<T>& L = ensemble[e];
            const LayerResult<T>& r = results.back()[e];
            double nj = 0.0;
            for (std::size_t c = 1; c <= NT(); ++c) nj = std::max(nj, njobs(tidx, c));
            // the task's OWN station, which under `flat.cs` is one of many
            const std::size_t ts = station_idx_of(L, tidx);
            T tp = Tzero(), ut = Tzero();
            for (std::size_t k = 0; k < L.nclasses; ++k) {
                tp += r.TN(ts - 1, k);
                ut += r.UN(ts - 1, k);
            }
            tput[tidx] = T(num_traits<T>::from_double(lqn.repl[tidx]) * tp);
            util[tidx] = ut;
            T raw;
            if (lqn.sched[tidx] == SchedStrategy::INF) {
                // an infinite server reports utilization as a mean job count
                raw = tput[tidx] == Tzero()
                          ? Tzero()
                          : T((num_traits<T>::from_double(nj) - util[tidx]) / tput[tidx] - ztask);
            } else {
                const T om = util[tidx] > Tone() ? T(util[tidx] - Tone()) : T(Tone() - util[tidx]);
                raw = tput[tidx] == Tzero()
                          ? Tzero()
                          : T(num_traits<T>::from_double(nj) * om / tput[tidx] - ztask);
            }
            const T floorv = num_traits<T>::from_double(GlobalConstants::Zero);
            thinkt[tidx] = raw < floorv ? floorv : raw;

            // The cold start goes the OTHER way from the phase-2 tail. A caller
            // class cycles as delay plus station service, and the station serves
            // only the host demand: the charge is on the ENTRY, not on any
            // activity's demand, so the station never sees it and the delay has
            // to carry it. Without this the callee layer cycles faster than its
            // callers drive it -- 0.529412 against 0.5 on lqn_setup, with the
            // caller conserved and the callee not. Zero for a task with no setup
            // (updateThinkTimes.m).
            {
                const double c = setup_charge(tidx);
                if (c != 0.0) {
                    const T withc = T(thinkt[tidx] + num_traits<T>::from_double(c));
                    thinkt[tidx] = withc < floorv ? floorv : withc;
                }
            }

            if (relax_omega < 1.0 && it > 1 && !std::isnan(thinkt_prev[tidx])) {
                const double rawd = dbl(thinkt[tidx]);
                if (thinkt_prev[tidx] > 10.0 * rawd && rawd > GlobalConstants::FineTol) {
                    thinkt_prev[tidx] = rawd;
                    thinkt_prev_v[tidx] = thinkt[tidx];
                }
                const T om = num_traits<T>::from_double(relax_omega);
                const T om1 = num_traits<T>::from_double(1.0 - relax_omega);
                thinkt[tidx] = T(om * thinkt[tidx] + om1 * thinkt_prev_v[tidx]);
            }
            thinkt_prev[tidx] = dbl(thinkt[tidx]);
            thinkt_prev_v[tidx] = thinkt[tidx];
            thinktproc[tidx] = Distrib<T>::exp_mean(T(thinkt[tidx] + ztask));
        }
    }

    // -----------------------------------------------------------------------
    // updateLayers and updateRoutingProbabilities
    // -----------------------------------------------------------------------
    void update_layers(int it) {
        // Under "srvn.ph" the layer classes are one per caller task and their laws
        // are composed, not read off the update maps
        if (is_ph_encoding()) {
            update_layers_ph(it);
            return;
        }
        const bool elevator = (it % 2) == 1;
        const std::size_t nt = thinkt_map.size();
        for (std::size_t r = 0; r < nt; ++r) {
            const UpdRow& row = thinkt_map[elevator ? nt - 1 - r : r];
            qn::Layer<T>& L = ensemble[std::size_t(idxhash[row.idx])];
            if (row.aidx <= NT() && opt.interlocking &&
                L.classes[row.cls - 1].type == JobClassType::CLOSED)
                L.classes[row.cls - 1].population = njobs(row.aidx, row.idx);
            if (row.node == L.clientIdx) {
                if (lqn.type[row.aidx] == LqnElement::TASK) {
                    if (lqn.sched[row.aidx] != SchedStrategy::REF) {
                        if (!thinktproc[row.aidx].disabled)
                            L.set_service(row.node, row.cls, thinktproc[row.aidx]);
                    } else {
                        L.set_service(row.node, row.cls, servtproc[row.aidx]);
                    }
                } else {
                    L.set_service(row.node, row.cls, servtproc[row.aidx]);
                }
            } else {
                L.set_service(row.node, row.cls, servtproc[row.aidx]);
            }
        }
        const std::size_t nc = call_map.size();
        for (std::size_t r = 0; r < nc; ++r) {
            const UpdRow& row = call_map[elevator ? nc - 1 - r : r];
            qn::Layer<T>& L = ensemble[std::size_t(idxhash[row.idx])];
            if (row.node == L.clientIdx)
                L.set_service(row.node, row.cls, callservtproc[row.aidx]);
            else
                L.set_service(row.node, row.cls, servtproc[lqn.callpair_dst[row.aidx]]);
        }
        // Async arrival rate: the caller activity fires at its own throughput,
        // and every firing releases callmean jobs into this layer's Source.
        // updateLayers.m:63-73 replays the same map; the geometric self-loop in
        // build_layer already accounts for callmean, so the rate is the bare
        // activity throughput, not scaled by it.
        for (const UpdRow& row : arv_call_map) {
            qn::Layer<T>& L = ensemble[std::size_t(idxhash[row.idx])];
            const Distrib<T>& d = tputproc[lqn.callpair_src[row.aidx]];
            if (!d.disabled) L.set_service(row.node, row.cls, d);
        }
    }

    /**
     * Time a job of `cls`'s CHAIN spends blocked at this layer's region.
     *
     * A job held by an admission constraint is at no station at all, so its
     * wait is structurally absent from the RN and WN the layer reports back and
     * the caller and the callee end up disagreeing on throughput -- the fixed
     * point still converges, it just converges to the wrong place.
     *
     * Recovered by Little over the CHAIN, never per class: a job switches class
     * along the activity graph, so a call class carries population 0 and a
     * per-class deficit comes out negative and silently does nothing.
     * Zero for a layer with no region, which is every layer in a plain model.
     */
    T region_wait(std::size_t e, std::size_t cls) const {
        const qn::Layer<T>& L = ensemble[e];
        if (L.regions.empty() || results.empty()) return Tzero();
        const LayerResult<T>& r = results.back()[e];
        const std::size_t k = cls - 1;
        std::size_t c = L.nchains;
        for (std::size_t cc = 0; cc < L.nchains; ++cc)
            if (L.chains[cc][k]) c = cc;
        if (c == L.nchains) return Tzero();

        // the station the region is stated over, which under `flat.cs` is one
        // among many and under `srvn` is the layer's own server
        std::size_t rs = L.serverIdx;
        for (std::size_t i = 0; i < L.regions[0].members.size(); ++i)
            if (L.regions[0].members[i]) { rs = i + 1; break; }

        double pop = 0.0;
        T inside = Tzero(), tput_srv = Tzero();
        for (std::size_t j = 0; j < L.nclasses; ++j) {
            if (!L.chains[c][j]) continue;
            const double p = L.classes[j].population;
            if (!std::isfinite(p)) return Tzero();  // an open chain has no population to close on
            pop += p;
            for (std::size_t i = 0; i < L.nstations; ++i) inside += r.QN(i, j);
            tput_srv += r.TN(rs - 1, j);
        }
        const T deficit = T(num_traits<T>::from_double(pop) - inside);
        if (!(deficit > Tzero()) || !(tput_srv > Tzero())) return Tzero();
        return T(deficit / tput_srv);
    }

    /** Port of updateRoutingProbabilities: entry selection by throughput ratio. */
    void update_routing_probabilities(int) {
        for (std::size_t u = 0; u < unique_route_idx.size(); ++u) {
            // the reference always takes the reversed order here: its `mod(it,0)`
            // guard is NaN, which MATLAB reads as false
            const std::size_t idx = unique_route_idx[unique_route_idx.size() - 1 - u];
            qn::Layer<T>& L = ensemble[std::size_t(idxhash[idx])];
            bool updated = false;
            for (const RouteRow& r : route_map) {
                if (r.idx != idx) continue;
                if (idxhash[r.tidx_caller] < 0) continue;
                const std::size_t cl = std::size_t(idxhash[r.tidx_caller]);
                const qn::Layer<T>& CL = ensemble[cl];
                const LayerResult<T>& cr = results.back()[cl];
                // the CALLER's own station, and for a call class the station of
                // the entry it targets; the two coincide under `srvn`
                const std::size_t cs = station_idx_of(CL, r.tidx_caller);
                T Xtot = Tzero();
                for (std::size_t k = 0; k < CL.nclasses; ++k) Xtot += cr.TN(cs - 1, k);
                if (!(Xtot > Tzero())) continue;
                T entry_tput = Tzero();
                for (const auto& a : CL.attr_calls)
                    if (a[3] == r.eidx)
                        entry_tput += cr.TN(station_idx_of(CL, lqn.parent[a[3]]) - 1, a[0] - 1);
                L.set_route(r.cfrom, r.cto, r.nodefrom, r.nodeto, T(entry_tput / Xtot));
                updated = true;
            }
            if (updated) L.refresh_chains();
        }
    }

    // -----------------------------------------------------------------------
    // getTranAvg: getTranAvgDecoupled and getTranAvgCoupled
    // -----------------------------------------------------------------------

    /**
     * What a layered transient needs before it can mean anything.
     *
     * A FLUID ENSEMBLE IS REQUIRED, and not as an implementation shortcut: the
     * reference reaches the layer transient through `self.solvers{e}.getTranAvg`,
     * which only a transient-capable layer solver has, and the coupled mode
     * injects its time-varying demands through the fluid rate schedule
     * specifically. An MVA layer has no trajectory to report and no place to
     * receive an injection, so an ensemble built on one is refused here rather
     * than silently answered with its steady state repeated over a grid.
     */
    void require_transient_ready() const {
        if (opt.layer_solver != "fluid")
            throw UnsupportedError(
                "SolverLN: the layered transient is the transient OF EACH LAYER, which only the "
                "fluid layer solver has; rebuild the ensemble with layer_solver 'fluid'");
    }

    /** True when the caller named a horizon; otherwise each layer picks its own. */
    bool has_finite_horizon() const {
        return std::isfinite(opt.timespan_end) && opt.timespan_end > 0.0;
    }

    /** One layer's trajectory, resampled onto the shared grid. */
    struct TranTraj {
        std::vector<std::vector<std::vector<double>>> Q, U, Tp, R;  ///< [station][class][point]
    };

    /**
     * Run one layer's transient, optionally with an injected rate schedule, and
     * return both the reportable block and the trajectory the relaxation reads.
     */
    void run_layer_transient(std::size_t e, const std::vector<FluidRateSched>& sched,
                             LnTranLayer& block, TranTraj& traj) {
        qn::Layer<T>& L = ensemble[e];
        L.refresh_rt();
        fluid::FluidOptions fo = opt.layer_fluid;
        fo.rate_sched = sched;
        // THE WARM START IS REPLAYED HERE, and this is the only place it can be:
        // `converged` resets every layer the first time the fixed point settles,
        // so a state installed before the solve is gone by now. The steady solve
        // ignores an initial state and the transient does not, which is why
        // replaying it after the one and before the other loses nothing.
        if (e < layer_tran_init.size() && !layer_tran_init[e].empty())
            fo.init_sol = layer_tran_init[e];
        // No horizon named: this layer picks its own, by the analyzer's own rule
        // (thirty mean events of its slowest transition). That is what the
        // reference's decoupled path does -- each layer's SolverFluid chooses --
        // and it is why an unset timespan is a valid request rather than an error.
        const double t_end =
            has_finite_horizon() ? opt.timespan_end : fluid::fluid_default_horizon(L, fo);
        const std::vector<fluid::FluidTranPoint> pts =
            detail::ln_fluid_transient(L, fo, t_end, opt.tran_points, opt.tran_grid);
        const std::size_t M = L.nstations, K = L.nclasses, P = pts.size();
        block.t.resize(P);
        auto alloc = [&](std::vector<std::vector<std::vector<double>>>& A) {
            A.assign(M, std::vector<std::vector<double>>(K, std::vector<double>(P, 0.0)));
        };
        alloc(block.QN);
        alloc(block.UN);
        alloc(block.TN);
        alloc(traj.Q);
        alloc(traj.U);
        alloc(traj.Tp);
        alloc(traj.R);
        for (std::size_t p = 0; p < P; ++p) {
            block.t[p] = pts[p].t;
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    const double q = pts[p].QN(i, r), u = pts[p].UN(i, r), x = pts[p].TN(i, r);
                    block.QN[i][r][p] = q;
                    block.UN[i][r][p] = u;
                    block.TN[i][r][p] = x;
                    traj.Q[i][r][p] = q;
                    traj.U[i][r][p] = u;
                    traj.Tp[i][r][p] = x;
                    // Residence by Little, which is how the relaxation reads a
                    // callee's response time off a layer trajectory.
                    traj.R[i][r][p] = q / std::max(x, GlobalConstants::FineTol);
                }
        }
    }

    /**
     * Port of getTranAvgDecoupled: freeze the inter-layer demands at the
     * converged fixed point and run each layer's transient in isolation.
     */
    LnTranSolution tran_avg_decoupled() {
        require_transient_ready();
        if (results.empty()) iterate();
        LnTranSolution out;
        out.mode = "decoupled";
        out.layers.resize(ensemble.size());
        for (std::size_t e = 0; e < ensemble.size(); ++e) {
            TranTraj tj;
            run_layer_transient(e, std::vector<FluidRateSched>(), out.layers[e], tj);
        }
        return out;
    }

    /**
     * Port of getTranAvgCoupled: reconcile the per-layer transients by waveform
     * relaxation, so the layer populations and the inter-layer demands co-evolve
     * in model time.
     *
     * Each layer's fluid transient is driven by TIME-VARYING inter-layer demand
     * trajectories taken from the other layers' latest transients, and the loop
     * repeats until the trajectories stop moving in sup-norm. Iteration 0 uses
     * the frozen equilibrium demands, so it reproduces the decoupled answer
     * exactly; at convergence every layer relaxes to its own fixed point, so the
     * endpoint equals getAvg.
     *
     * TWO CHANNELS ARE COUPLED, the task think times (the client delay) and the
     * synchronous-call service demands (the caller's client station). Both are
     * dominant inter-layer couplings; the intra-layer host service stays at its
     * equilibrium value, as in the reference.
     */
    LnTranSolution tran_avg_coupled() {
        require_transient_ready();
        // Waveform relaxation reconciles the layers on ONE shared grid, so it
        // needs a horizon the layers agree on. With none named there is nothing
        // to co-evolve over, and the reference defers to the decoupled path
        // rather than inventing one; so does this.
        if (!has_finite_horizon()) return tran_avg_decoupled();
        if (results.empty()) iterate();
        const std::size_t E = ensemble.size();

        LnTranSolution out;
        out.mode = "coupled";
        out.layers.resize(E);
        std::vector<TranTraj> traj(E), prev(E);
        for (std::size_t e = 0; e < E; ++e)
            run_layer_transient(e, std::vector<FluidRateSched>(), out.layers[e], traj[e]);
        const std::vector<double> tgrid = out.layers.empty() ? std::vector<double>() : out.layers[0].t;

        for (long iter = 1; iter <= opt.ln_transient_iter_max; ++iter) {
            prev = traj;
            const std::vector<std::vector<FluidRateSched>> sched =
                build_rate_sched(recompute_demand(traj, tgrid), tgrid);
            for (std::size_t e = 0; e < E; ++e)
                run_layer_transient(e, sched[e], out.layers[e], traj[e]);
            double gap = 0.0;
            for (std::size_t e = 0; e < E; ++e)
                for (std::size_t i = 0; i < traj[e].Q.size(); ++i)
                    for (std::size_t r = 0; r < traj[e].Q[i].size(); ++r)
                        for (std::size_t p = 0; p < traj[e].Q[i][r].size(); ++p)
                            gap = std::max(gap, std::fabs(traj[e].Q[i][r][p] - prev[e].Q[i][r][p]));
            out.iterations = iter;
            out.gap = gap;
            if (gap < opt.ln_transient_tol) break;
        }
        return out;
    }

    /** The two coupled demand channels, as trajectories on the shared grid. */
    struct TranDemand {
        std::map<std::size_t, std::vector<double>> thinkt;     ///< by task element index
        std::map<std::size_t, std::vector<double>> callservt;  ///< by call index
    };

    /**
     * Recompute the inter-layer demands pointwise in t, mirroring the SCALAR
     * updateThinkTimes and updateMetricsDefault formulas.
     *
     * They are the same formulas, evaluated at each point of the grid instead of
     * at the fixed point: that is what makes the endpoint of the relaxation the
     * steady-state answer rather than something near it.
     */
    TranDemand recompute_demand(const std::vector<TranTraj>& traj,
                                const std::vector<double>& tgrid) const {
        TranDemand out;
        const std::size_t ng = tgrid.size();

        for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
            const std::size_t tidx = lqn.tshift + t;
            if (idxhash[tidx] < 0 || lqn.isref[tidx]) continue;
            const std::size_t e = std::size_t(idxhash[tidx]);
            const qn::Layer<T>& L = ensemble[e];
            const std::size_t s = station_idx_of(L, tidx) - 1;
            double nj = 0.0;
            for (std::size_t c = 1; c <= NT(); ++c) nj = std::max(nj, njobs(tidx, c));
            // same closure as update_think_times, so the same gate
            const double userthink = dbl(ref_think_time(tidx));
            std::vector<double> tk(ng, 0.0);
            for (std::size_t p = 0; p < ng; ++p) {
                double U = 0.0, X = 0.0;
                for (std::size_t r = 0; r < L.nclasses; ++r) {
                    U += traj[e].U[s][r][p];
                    X += traj[e].Tp[s][r][p];
                }
                const double Xs = std::max(X, GlobalConstants::FineTol);
                double v = lqn.sched[tidx] == SchedStrategy::INF
                               ? (nj - U) / Xs - userthink
                               : nj * std::fabs(1.0 - U) / Xs - userthink;
                if (v < GlobalConstants::Zero) v = GlobalConstants::Zero;
                tk[p] = v + userthink;  // total mean, user think included
            }
            out.thinkt[tidx] = tk;
        }

        for (std::size_t cidx = 1; cidx <= lqn.ncalls; ++cidx) {
            if (lqn.calltype[cidx] != CallType::SYNC) continue;
            const std::size_t eidx = lqn.callpair_dst[cidx];
            const std::size_t tidx = lqn.parent[eidx];
            if (tidx > NT() || idxhash[tidx] < 0) continue;
            const std::size_t e = std::size_t(idxhash[tidx]);
            const qn::Layer<T>& L = ensemble[e];
            const std::size_t s = station_idx_of(L, tidx) - 1;
            std::vector<double> Rc(ng, 0.0);
            bool any = false;
            for (std::size_t r = 0; r < L.nclasses; ++r) {
                if (L.classes[r].attr_kind != int(LqnElement::ENTRY)) continue;
                if (L.classes[r].attr_idx != eidx) continue;
                any = true;
                for (std::size_t p = 0; p < ng; ++p) Rc[p] += traj[e].R[s][r][p];
            }
            if (!any) {
                // No entry class of its own in this layer: the entry's work is
                // carried by its activities, so their residence is the answer.
                for (std::size_t r = 0; r < L.nclasses; ++r)
                    for (std::size_t p = 0; p < ng; ++p) Rc[p] += traj[e].R[s][r][p];
            }
            const double cm = dbl(lqn.callproc_mean[cidx]);
            for (std::size_t p = 0; p < ng; ++p) Rc[p] *= cm;
            out.callservt[cidx] = Rc;
        }
        return out;
    }

    /**
     * Map the demand trajectories onto per-layer rate schedules, through the
     * SAME update maps that place the scalar setService calls.
     *
     * The injected schedule MODULATES the layer's equilibrium rate by the ratio
     * of the transient demand to its steady-state value: passing rate = 1/d(t)
     * with nominal = 1/d(end) makes solver_fluid_ratemult's multiplier
     * d(end)/d(t), which is exactly 1 at the horizon end. That is what makes the
     * layer relax to its UNMODIFIED fixed point regardless of any small mismatch
     * between the fluid residence Q/T and the scalar equilibrium demand.
     */
    std::vector<std::vector<FluidRateSched>> build_rate_sched(
        const TranDemand& demand, const std::vector<double>& tgrid) const {
        std::vector<std::vector<FluidRateSched>> out(ensemble.size());
        const bool want_think =
            opt.ln_transient_channels == "both" || opt.ln_transient_channels == "thinkt";
        const bool want_call =
            opt.ln_transient_channels == "both" || opt.ln_transient_channels == "callservt";

        auto add = [&](std::size_t e, std::size_t station, std::size_t cls,
                       const std::vector<double>& d) {
            if (d.empty()) return;
            const double dend = d.back();
            if (!(dend > GlobalConstants::FineTol)) return;  // degenerate steady state
            // Bound the transient demand to a physical band around its
            // steady-state value: an early-transient throughput near zero sends
            // the reciprocal to infinity and the integrator with it.
            const double cap = 20.0;
            FluidRateSched s;
            s.station = station;
            s.cls = cls;
            s.tgrid = tgrid;
            s.rates.resize(d.size());
            for (std::size_t p = 0; p < d.size(); ++p) {
                const double v = std::min(std::max(d[p], dend / cap), dend * cap);
                s.rates[p] = 1.0 / v;
            }
            s.nominal = 1.0 / dend;
            out[e].push_back(s);
        };

        if (want_think)
            for (const UpdRow& row : thinkt_map) {
                if (idxhash[row.idx] < 0) continue;
                const std::size_t e = std::size_t(idxhash[row.idx]);
                if (row.node != ensemble[e].clientIdx) continue;
                if (lqn.type[row.aidx] != LqnElement::TASK) continue;
                if (lqn.sched[row.aidx] == SchedStrategy::REF) continue;
                const auto it = demand.thinkt.find(row.aidx);
                if (it == demand.thinkt.end()) continue;
                add(e, row.node, row.cls, it->second);
            }
        if (want_call)
            for (const UpdRow& row : call_map) {
                if (idxhash[row.idx] < 0) continue;
                const std::size_t e = std::size_t(idxhash[row.idx]);
                if (row.node != ensemble[e].clientIdx) continue;
                const auto it = demand.callservt.find(row.aidx);
                if (it == demand.callservt.end()) continue;
                add(e, row.node, row.cls, it->second);
            }
        return out;
    }

    // -----------------------------------------------------------------------
    // getEnsembleAvg
    // -----------------------------------------------------------------------
    LnSolution<T> aggregate() {
        const std::size_t N = lqn.nidx;
        LnSolution<T> s;
        auto mk = [&](std::vector<T>& v, std::vector<bool>& d) {
            v.assign(N + 1, Tzero());
            d.assign(N + 1, false);
        };
        std::vector<T> QN, UN, RN, TN, PN, SN, WN, AN;
        std::vector<bool> dQ, dU, dR, dT, dP, dS, dW, dA;
        mk(QN, dQ); mk(UN, dU); mk(RN, dR); mk(TN, dT);
        mk(PN, dP); mk(SN, dS); mk(WN, dW); mk(AN, dA);
        std::vector<bool> wn_done(N + 1, false);

        for (std::size_t e = 0; e < ensemble.size(); ++e) {
            const qn::Layer<T>& L = ensemble[e];
            const LayerResult<T>& r = results.back()[e];
            const std::size_t clientIdx = L.clientIdx;
            // The processors this layer serves: one under `srvn` (and only when
            // it is a host layer), every one of them under `flat.cs`. Each is
            // charged with the activities that actually run on it, which is
            // vacuous under `srvn` because a host layer holds no others.
            const bool has_host_server = !L.host_stations.empty();
            for (std::size_t hs : L.host_stations) {
                const std::size_t hidx = L.stations[hs - 1].attr_idx;
                dT[hidx] = true;
                TN[hidx] = Tzero();
                dP[hidx] = true;
                PN[hidx] = Tzero();
                for (std::size_t k = 0; k < L.nclasses; ++k) {
                    if (L.classes[k].completes) {
                        T t = clientIdx > 0 ? r.TN(clientIdx - 1, k) : Tzero();
                        const T ts = r.TN(hs - 1, k);
                        TN[hidx] = T(TN[hidx] + (t > ts ? t : ts));
                    }
                    if (L.classes[k].attr_kind == int(LqnElement::ACTIVITY)) {
                        // the activity does not run on this processor
                        if (station_idx_of_class(L, k) != hs) continue;
                        const std::size_t aidx = L.classes[k].attr_idx;
                        const std::size_t tidx = lqn.parent[aidx];
                        dP[aidx] = true;
                        dP[tidx] = true;
                        PN[aidx] = T(PN[aidx] + r.UN(hs - 1, k));
                        PN[tidx] = T(PN[tidx] + r.UN(hs - 1, k));
                        PN[hidx] = T(PN[hidx] + r.UN(hs - 1, k));
                    }
                }
                dT[hidx] = false;  // NaN in the reference, for consistency with LQNS
            }

            for (std::size_t k = 0; k < L.nclasses; ++k) {
                const int kind = L.classes[k].attr_kind;
                // the station this class is actually served at, which under
                // `flat.cs` is the processor of an activity or the called task
                // of a call rather than the one server the layer used to have
                const std::size_t serverIdx = station_idx_of_class(L, k);
                if (kind == int(LqnElement::TASK)) {
                    const std::size_t tidx = L.classes[k].attr_idx;
                    if (has_host_server && !dT[tidx]) {
                        dT[tidx] = true;
                        TN[tidx] = r.TN(clientIdx - 1, k);
                    }
                } else if (kind == int(LqnElement::ENTRY)) {
                    const std::size_t eidx = L.classes[k].attr_idx;
                    dS[eidx] = true;
                    // getEnsembleAvg.m:84-92: a phase-2 entry reports the
                    // CALLER's view, which is what residt now holds.
                    SN[eidx] = (has_phase2 && dbl(servt_ph2[eidx]) > GlobalConstants::FineTol)
                                   ? residt[eidx]
                                   : servt[eidx];
                    if (has_host_server && !dT[eidx]) {
                        dT[eidx] = true;
                        TN[eidx] = r.TN(clientIdx - 1, k);
                    }
                } else if (kind == int(LqnElement::CALL)) {
                    const std::size_t cidx = L.classes[k].attr_idx;
                    const std::size_t aidx = lqn.callpair_src[cidx];
                    if (lqn.calltype[cidx] == CallType::SYNC) {
                        dS[aidx] = true;
                        SN[aidx] = T(SN[aidx] + r.RN(serverIdx - 1, k) * lqn.callproc_mean[cidx]);
                    }
                    dQ[aidx] = true;
                    QN[aidx] = T(QN[aidx] + r.QN(serverIdx - 1, k));
                } else if (kind == int(LqnElement::ACTIVITY)) {
                    const std::size_t aidx = L.classes[k].attr_idx;
                    const std::size_t tidx = lqn.parent[aidx];
                    dQ[tidx] = true;
                    QN[tidx] = T(QN[tidx] + r.QN(serverIdx - 1, k));
                    dT[aidx] = true;
                    dQ[aidx] = true;
                    TN[aidx] = T(TN[aidx] + r.TN(serverIdx - 1, k));
                    dS[aidx] = true;
                    SN[aidx] = T(SN[aidx] + r.RN(serverIdx - 1, k));
                    dR[aidx] = true;
                    RN[aidx] = T(RN[aidx] + r.RN(serverIdx - 1, k));
                    dW[aidx] = true;
                    dW[tidx] = true;
                    WN[aidx] = residt[aidx];
                    if (!wn_done[aidx]) {
                        WN[tidx] = T(WN[tidx] + residt[aidx]);
                        wn_done[aidx] = true;
                    }
                    QN[aidx] = T(QN[aidx] + r.QN(serverIdx - 1, k));
                }
            }
        }

        for (std::size_t e = 1; e <= lqn.nentries; ++e) {
            const std::size_t eidx = lqn.eshift + e;
            const std::size_t tidx = lqn.parent[eidx];
            dU[tidx] = true;
            dU[eidx] = true;
            // getEnsembleAvg.m:186-198: the server is busy through BOTH phases,
            // so the utilization is the full service time and not SN, which for
            // a phase-2 entry has been cut down to the caller's view above.
            UN[eidx] = (has_phase2 && dbl(servt_ph2[eidx]) > GlobalConstants::FineTol)
                           ? T(TN[eidx] * (servt_ph1[eidx] + servt_ph2[eidx]))
                           : T(TN[eidx] * SN[eidx]);
            T ps = Tzero();
            bool any = false;
            for (std::size_t a : lqn.actsof[eidx])
                if (dP[a]) {
                    ps += PN[a];
                    any = true;
                }
            if (!lqn.actsof[eidx].empty()) {
                dP[eidx] = true;
                PN[eidx] = any ? ps : Tzero();
            }
            for (std::size_t a : lqn.actsof[tidx]) {
                dU[a] = true;
                UN[a] = T(TN[a] * SN[a]);
            }
            UN[tidx] = T(UN[tidx] + UN[eidx]);
        }

        // AN IGNORED ELEMENT IS IDLE, NOT UNDEFINED, and the two are different cells.
        // Its component holds no reference task, so nothing reaches it and every
        // measure it HAS is zero -- but the measures its kind never has stay
        // undefined, exactly as they do for a reachable element. A flat zero over
        // all six columns broke the table's NaN mask (a processor with a queue
        // length of 0, an arrival rate reported where no solver reports one), and
        // the mask is part of the answer: see _kb/06-solver-catalog.md. The
        // relabelling below reports Q from U, U from P and R from S, so QN and RN
        // are dead here and are not written.
        for (std::size_t i = 1; i <= N; ++i)
            if (ignore[i]) {
                PN[i] = Tzero();  // every kind reports a utilization
                dP[i] = true;
                dA[i] = false;    // nothing reports an arrival rate on an LQN
                const bool host = lqn.type[i] == LqnElement::HOST;
                const bool task = lqn.type[i] == LqnElement::TASK;
                const bool entry = lqn.type[i] == LqnElement::ENTRY;
                UN[i] = Tzero();
                dU[i] = !host;
                SN[i] = Tzero();
                dS[i] = !host && !task;
                WN[i] = Tzero();
                dW[i] = !host && !entry;
                TN[i] = Tzero();
                dT[i] = !host;
            }

        // the reference's final relabelling: Q <- U, U <- P, R <- S
        s.QN = UN; s.defined_Q = dU;
        s.UN = PN; s.defined_U = dP;
        s.RN = SN; s.defined_R = dS;
        s.TN = TN; s.defined_T = dT;
        s.AN = AN; s.defined_A = dA;
        s.WN = WN; s.defined_W = dW;
        s.iterations = iterations_done;
        s.converged = did_converge;
        return s;
    }
};

}  // namespace ln
}  // namespace line

// Deliberately at the FOOT of the file: lqn_analyzers.h needs LayerResult and
// SolverLN complete, and this file needs its lqn_overtake_prob_markov, so the
// two are mutually dependent. Either include order works -- whichever header a
// translation unit names first, the other is fully parsed before the templates
// above are instantiated -- and the forward declaration near the top of this
// file is what makes the call inside update_metrics resolve.
#include "line/solvers/ln/lqn_analyzers.h"

#endif  // LINE_SOLVERS_LN_SOLVER_LN_H
