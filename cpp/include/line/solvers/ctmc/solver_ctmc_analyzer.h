/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `solver_ctmc_analyzer.m` and the parts of `@@SolverCTMC/runAnalyzer.m`
 * that surround one solve: the method gate, the open-model CUTOFF, the state
 * space and synchronization construction, the stationary solve with its
 * reducible-generator handling, and the mapping onto the AvgTable columns.
 *
 * WHAT THE CUTOFF IS, AND WHY IT IS NOT A TOLERANCE. A closed class carries its
 * own population, so its lattice is finite. An open class does not, so the
 * chain is infinite and the reference TRUNCATES it at `options.cutoff` jobs per
 * class. The answer is therefore the exact stationary law of a DIFFERENT chain,
 * one the truncation defined, and it converges to the model's only as the
 * cutoff grows. That is why the reference prints a warning on every open model
 * and why the banner here reports the cutoff it used: a CTMC number for an open
 * model is not a number without it.
 *
 * REDUCIBILITY. `space_generator` enumerates every state the ENCODING admits,
 * and the dynamics need not reach all of them; the generator is then reducible
 * and pi is not unique. The reference resolves this by keeping the weakly
 * connected component of the INITIAL state, and falling back to the largest
 * component when the initial state is not in the space at all. Solving the whole
 * generator instead -- which is what `ctmc_solve` does on its own, splitting
 * per component and renormalizing -- spreads mass over states the model can
 * never occupy and moves every reported mean.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_ANALYZER_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_ANALYZER_H

#include "line/util/line_console.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_memory_gate.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_state_space_logsize.h"
#include "line/api/fes/fes_aggregate.h"
#include "line/api/fes/fes_compute_metrics.h"
#include "line/api/sn/sn_aggregate_chains.h"
#include "line/lang/qn/fj_tag.h"
#include "line/api/sn/sn_nonmarkov_toph.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/lang/qn/state.h"
#include "line/lang/qn/state_events.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/ctmc/ctmc_stationary.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_fcr.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/tr/fj_tag_transform.h"
#include "line/solvers/tr/transform_solve.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

/**
 * `SolverOptions.m:107`: the per-class state-space cutoff SolverCTMC defaults an
 * open or mixed model to. It is a SOLVER default, not a fallback computed from
 * the model, and the reference's `ceil(6000^(1/(M*K)))` is reached only from an
 * explicitly infinite request. See `resolve_cutoff`.
 */
constexpr std::size_t CTMC_DEFAULT_CUTOFF = 10;

/**
 * The SolverCTMC knobs this port honours.
 *
 * `cutoff` is per-class when `cutoff_vec` is given and uniform otherwise; a
 * negative `cutoff` means "not given" and takes the reference's solver default
 * `CTMC_DEFAULT_CUTOFF`, while an infinite one takes its automatic value,
 * `ceil(6000^(1/(M*K)))`.
 */
struct CtmcOptions {
    std::string method = "default";
    /**
     * Keep the per-synchronization EVENT FILTRATION alongside Q.
     *
     * Off by default because it costs one n x n matrix per synchronization,
     * which on a model with many routing pairs dwarfs the generator itself. The
     * response-time CDF needs it, since splitting Q on one event cannot be done
     * after the contributions have been summed.
     */
    bool keep_filtration = false;
    double cutoff = -1.0;                  ///< < 0 = not given
    std::vector<std::size_t> cutoff_vec;   ///< per-class override, or empty
    /**
     * `options.cutoff` AS A (station x class) MATRIX, or empty.
     *
     * The reference accepts a matrix wherever a model needs a different
     * truncation per station -- oqn_cs_routing writes `[1,1,0;3,3,0;0,0,3]`,
     * which bounds each queue only in the classes it actually serves. Reducing
     * that to its per-class maximum enumerates a far larger chain and answers a
     * different model, so it is carried whole and applied per station in
     * `space_capacity_c`.
     */
    std::vector<std::vector<std::size_t>> cutoff_mat;
    std::size_t state_max = 3000000;       ///< refuse a space larger than this
    /// `options.force`: downgrade the memory pre-gate's refusal to a warning.
    bool force = false;
    /// `options.memorySafetyFraction`: share of available memory a solve may target.
    double memory_safety_fraction = mc::CTMC_DEFAULT_SAFETY_FRACTION;
    /**
     * `options.config.nonmkvorder`: the phase budget `sn_nonmarkov_toph` spends
     * on a non-Markovian service law. The reference default is 20, and it is a
     * real cost here -- every phase multiplies the state space.
     */
    std::size_t nonmkv_order = 20;
    /**
     * `options.timestep`: the FIXED OUTPUT STEP of a transient analysis.
     *
     * <= 0 means adaptive, which is the reference's default and its `[]`. It
     * changes WHERE the solution is reported and not how it is computed: the
     * integrator takes the same steps either way, and the grid points are read
     * off its interpolant. A caller comparing two transients needs the second,
     * because two adaptive solves land on different time vectors.
     */
    double timestep = -1.0;
    /**
     * `options.config.transient_method`: "ode" (the default) integrates the
     * forward equation, "fau" marches fast adaptive uniformization
     * (`mc::ctmc_fau`) over the output grid.
     *
     * It is a config key rather than a value of `method` because it changes no
     * stationary answer -- it is the transient path only -- and because the
     * reference's valid-method list is enumerated by its sanity harness, which
     * then wants a recorded baseline per method.
     */
    std::string transient_method = "ode";
    /**
     * `options.config.fau_epsilon`: total probability mass the whole grid may
     * discard under "fau". It is divided by the number of steps, each step
     * removing mass and none putting any back, so the accumulated defect stays
     * below it.
     */
    double fau_epsilon = 1e-6;
    /// `options.config.fau_delta`: occupancy below which a state is dropped.
    double fau_delta = 1e-12;
    /// `options.config.fau_ngrid`: output grid size when `timestep` is unset.
    std::size_t fau_ngrid = 100;
    /**
     * `options.config.chain_aggregation`: solve the CHAIN-AGGREGATED model.
     *
     * The state space grows with the per-class populations, so collapsing every
     * chain onto a single class is the standard way to make an otherwise
     * intractable multiclass model solvable: `api::sn_aggregate_chains` builds
     * the collapsed model and `mva::sn_deaggregate_chain_results` maps the
     * chain-level metrics back through alpha. EXACT on a product-form model, an
     * approximation otherwise, since one aggregate service law replaces the
     * per-class ones. Off by default: a caller who needs the exact multiclass
     * answer must pay the state space, not discover the trade after the fact.
     */
    bool chain_aggregation = false;

    /**
     * `options.config.transform='lc'`: solve by LOAD CONCEALMENT.
     *
     * Birman-Kogan Algorithm 2 as a model transformation: each chain is solved
     * on its own against the residual capacity the others leave it, with the
     * single-chain subproblem a real struct this analyzer solves. Distinct from
     * the `pfqn_bklc` KERNEL, which sweeps a demand matrix and stays the fast
     * path; see `tr::transform_solve_lc`.
     */
    bool load_concealment = false;

    /** Sweep cap for an iterated transformation; the kernel's own default. */
    std::size_t transform_iter_max = 1000;
    /**
     * `options.config.fes_stations`: 1-BASED station indices to collapse into a
     * flow-equivalent server before the chain is enumerated. Empty by default.
     *
     * The collapsed stations' own metrics are recovered by conditioning on the
     * FES population, so the table still names every station of the model the
     * caller built; see `solver_ctmc_fes_aggregation`.
     */
    std::vector<std::size_t> fes_stations;
};

/**
 * Everything one CTMC solve produces.
 *
 * `chain` is the generator RESTRICTED to the component that was kept, together
 * with the state space and the per-state arrival and departure rates on it. The
 * rates are carried rather than recomputed because Q has already summed every
 * synchronization's contribution into one entry, and a per-class rate cannot be
 * recovered from that sum afterwards.
 */
template <class T>
struct CtmcSolution {
    CtmcResult<T> chain;
    std::vector<T> pi;               ///< stationary distribution over chain.space
    CtmcAvg<T> avg;
    std::vector<std::size_t> cutoff; ///< the per-class cutoff actually used
    std::string actualmethod = "default";
    /**
     * `fjclassmap` when the model was fork-join, empty otherwise: the ORIGINAL
     * class of each auxiliary sibling class of the AUGMENTED struct that `chain`
     * and `pi` are indexed by. `avg` is already folded back onto the original
     * classes, so this is what relates the two index spaces.
     */
    std::vector<std::size_t> fjclassmap;
    /**
     * What the chain says about the model's Cache nodes: exact, since the hit
     * and miss shares are read off the stationary law rather than approximated.
     * Empty on a model with no Cache.
     */
    solvers::CacheMetrics<T> cache;
    /** Set when the chain is a reducible mixture solved from an invented seed. */
    std::string warning;
};

/**
 * Port of `SolverCTMC.listValidMethods`.
 *
 * `exact` is an explicit ALIAS for the default state-space path: it pins the
 * intent at the call site so an example or test cannot be re-baselined by a
 * later change of what `default` selects. It must stay behaviourally identical
 * to `default` -- nothing below branches on the name -- and that equivalence is
 * the point of the alias.
 *
 * `mdd` and the `cftp` pair never build the explicit generator, so they are
 * served by `solver_ctmc_mdd_analyzer` and `solver_ctmc_cftp` and return before
 * this file's state-space machinery is reached; they appear here only because
 * this is the list the gate below is read against.
 */
inline std::vector<std::string> list_valid_methods() {
    return {"default", "exact", "gpu", "mdd", "cftp", "cftp.approx"};
}

/** True for a method whose analyzer is not the explicit-generator one. */
inline bool is_stateless_method(const std::string& method) {
    return method == "mdd" || method == "cftp" || method == "cftp.approx";
}

/**
 * Port of `runAnalyzerChecks`' method gate.
 *
 * `gpu` IS ACCEPTED AND FALLS BACK, because that is what the reference does.
 * `ctmc_solve.m:224-238` wraps the `gpuArray` solve in a try/catch and, when the
 * GPU is absent or the solve throws, warns "GPU either not available or
 * execution failed. Switching to default method." and runs the same `Qnnz' \
 * bnnz` the default takes. A host with no Parallel Computing Toolbox therefore
 * gets the EXACT answer from `SolverCTMC(model,'gpu')`, and refusing it here
 * made a model MATLAB answers unanswerable under `lang='cpp'`. The fallback is
 * reported through `CtmcSolution::warning` rather than silently, so a caller
 * still learns the named backend did not run.
 *
 * `mdd` and `cftp` ARE refused here, because reaching THIS analyzer with one of
 * those names means the caller routed a generator-free method into the
 * generator path: the answer would be the enumerated one under a reported
 * method that never ran. Their own entry points do not call this gate.
 */
inline void check_method(const std::string& method) {
    const std::vector<std::string> valid = list_valid_methods();
    if (std::find(valid.begin(), valid.end(), method) == valid.end())
        throw UnsupportedError("SolverCTMC: the '" + method +
                               "' method is unsupported by this solver");
    if (is_stateless_method(method))
        throw UnsupportedError("SolverCTMC: the '" + method +
                               "' method does not enumerate a state space and is served by its "
                               "own analyzer, not by solver_ctmc_analyzer");
}

/** The reference's fallback warning, or empty when the method needs none. */
inline std::string method_fallback_warning(const std::string& method) {
    if (method != "gpu") return std::string();
    return "ctmc_solve: GPU either not available or execution failed. Switching to default "
           "method.";
}

/**
 * Refuse the constructs this port generates a chain for but does not MODEL.
 *
 * WHY A GATE AND NOT A TODO. A construct this port generates a chain for but does
 * not MODEL comes back as a complete, plausible AvgTable computed from a chain
 * that is not the model's, and a caller has no way to tell. Where the reference
 * declares the construct supported, silence here would show two different answers
 * with no indication which is wrong.
 *
 * As of this change the list is empty of structural gaps: fork-join goes through
 * `fj_tag`, class and joint dependence through `cd_factor`, and BAS through the
 * blocked marker. What remains are the DECLARATION defects below -- a dependence
 * handle with no peak -- plus the refusals stated by name elsewhere
 * (`getSymbolicGenerator`, a matrix-exponential `getProb`, `getCdfRespT` on an
 * open model, and every region rule but DROP and WAITQ).
 */
template <class T>
void ctmc_check_support(const NetworkStruct<T>& sn) {
    // Immediate feedback is a SAMPLE-PATH property: the fed-back job keeps the
    // server it just used instead of releasing it and re-queueing. The
    // reference expresses that in `State.afterEvent` through its
    // `immfeed_selfloop` argument; nothing here does, so the generator would
    // route the completion out and back like any other visit and return a
    // confident answer for a model with no feedback in it. Unlike SolverMVA and
    // SolverNC, which WARN because a mean-value argument cannot express the
    // held server either way, an exact chain either has the arc or is wrong.
    if (sn.has_immediate_feedback())
        throw UnsupportedError(
            "SolverCTMC: immediate feedback (setImmediateFeedback) is not ported. The reference "
            "keeps the server across the fed-back service through State.afterEvent's immfeed "
            "self-loop, which this port's event generator has no arc for, so the chain would "
            "release and re-queue the job instead; use SolverMVA or SolverNC, which warn and "
            "approximate it");
    // Server breakdowns are the same shape of gap, and they became REACHABLE on
    // 2026-08-15 when `sn.breakdownparam` was added for the LDES engine. The
    // reference's SolverCTMC does model them -- `State.afterEventStation`'s
    // FAILURE and REPAIR branches, on the shared trailing marker column that
    // `refreshLocalVars` reserves -- and `ctmc_feature_set` transcribes that
    // declaration, so without this the chain would be built with the server
    // permanently up and the answer would look ordinary.
    if (sn.has_breakdown())
        throw UnsupportedError(
            "SolverCTMC: server breakdowns (setBreakdown) are not ported. The reference carries "
            "the up/down status in the station's trailing local-state column and fires FAILURE "
            "and REPAIR against it (State.afterEventStation); this port allocates that column for "
            "a cache, a polling controller or the BAS blocked marker only, so the chain would be "
            "built as though the server never failed; use SolverLDES, whose native engine "
            "simulates the outage");
    // Class dependence requires a DECLARED peak, and unlike the handle itself
    // the peak cannot be derived: `cd_factor` scales the generator correctly
    // without it, but the utilization column would then silently fall back to the
    // busy-server probability, which at a station whose beta emulates extra
    // servers is not a utilization. `getLimitedClassDependencePeak` makes it
    // mandatory in the reference, so a missing peak is a model defect.
    for (std::size_t i = 0; i < sn.stations.size(); ++i) {
        if (sn.stations[i].cdscaling && sn.stations[i].cdscalingpeak.empty())
            throw InputError(
                "SolverCTMC: station '" + sn.stations[i].name +
                "' declares class-dependent service without a peak rate. Utilization at a "
                "class-dependent station is reported as T/mu/peak, so pass the peak to "
                "setClassDependence");
        if (sn.stations[i].jdscaling && sn.stations[i].jdscalingpeak.empty())
            throw InputError(
                "SolverCTMC: station '" + sn.stations[i].name +
                "' declares joint-dependent service without a peak rate; pass the peak to "
                "setJointDependence");
    }
    // BAS / BBS / RSRD are SOLVED, not refused, and the two tiers are the
    // reference's, not a shortcut here:
    //
    //   BAS gets the blocked marker. `refresh_bas_blocking` reserves the column,
    //   the generator emits the become-blocked edge, the departure handler clears
    //   it at 1e7 and `solver_ctmc_avg_from_pi` moves the held job to its
    //   destination. The chain is exact.
    //
    //   BBS and RSRD get the DISABLED DEPARTURE only -- `arrival_is_lost` returns
    //   false, so the completion that would move the job is simply not generated
    //   until room frees. `declaresBlockedMarker` tests BAS alone, so the
    //   reference does the same and does not distinguish repetitive service from
    //   a random redraw of the destination. The server freezes with the job in
    //   it, which is the right occupancy but not the right sample path for RSRD.
    //
    // `solver_ctmc_waitq` still refuses all three on ITS path: the WAITQ
    // generator augments the state with a per-region token FIFO and would have to
    // carry the marker through that augmentation as well.
}

namespace analyzer_detail {

/**
 * The per-class cutoff of `@@SolverCTMC/runAnalyzer.m`.
 *
 * A closed class never consults it; an open one takes the caller's value.
 *
 * THE DEFAULT IS 10, NOT THE 6000-STATE BUDGET. `SolverOptions.m:107` sets
 * `options.cutoff = 10` for SolverCTMC, so `runAnalyzer`'s
 * `ceil(6000^(1/(M*K)))` fallback is reached only when the caller explicitly
 * asks for an INFINITE cutoff, which is the one value the truncation cannot
 * honour. Taking the fallback as the default instead put this port one lattice
 * step below the reference on every open model and the answers differed in the
 * third digit with nothing to show for it: on the Exp(0.5) -> three-FCFS-queue
 * tandem the budget gives 9 and QLen at Q1 came out 0.9437940 against the
 * reference's 0.9446411, while cutoff 10 reproduces it to every digit.
 */
template <class T>
std::vector<std::size_t> resolve_cutoff(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::vector<std::size_t> cut(K, 0);
    bool any_open = false;
    for (std::size_t k = 0; k < K; ++k)
        if (!std::isfinite(sn.njobs()[k])) any_open = true;
    if (!any_open) return cut;

    if (!opt.cutoff_mat.empty()) {
        // `spaceGenerator.m`: the LATTICE height of an open class is the largest
        // bound any station grants it, `max(capacityc(:,r))`; the per-station
        // entries then bound each node individually.
        if (opt.cutoff_mat.size() != M)
            throw InputError("SolverCTMC: the cutoff matrix must have one row per station");
        for (std::size_t k = 0; k < K; ++k) {
            if (std::isfinite(sn.njobs()[k])) continue;
            for (std::size_t i = 0; i < M; ++i) {
                if (opt.cutoff_mat[i].size() != K)
                    throw InputError(
                        "SolverCTMC: the cutoff matrix must have one column per class");
                cut[k] = std::max(cut[k], opt.cutoff_mat[i][k]);
            }
        }
        return cut;
    }
    if (!opt.cutoff_vec.empty()) {
        if (opt.cutoff_vec.size() != K)
            throw InputError("SolverCTMC: the per-class cutoff must have one entry per class");
        return opt.cutoff_vec;
    }
    std::size_t c;
    if (opt.cutoff > 0 && std::isfinite(opt.cutoff)) {
        c = static_cast<std::size_t>(opt.cutoff);
    } else if (opt.cutoff < 0) {
        c = CTMC_DEFAULT_CUTOFF;
    } else {
        const double e = 1.0 / static_cast<double>(M * K);
        c = static_cast<std::size_t>(std::ceil(std::pow(6000.0, e)));
        if (c < 1) c = 1;
    }
    for (std::size_t k = 0; k < K; ++k)
        if (!std::isfinite(sn.njobs()[k])) cut[k] = c;
    return cut;
}

/**
 * Port of the `maxPending` rule of `State.spaceGeneratorNodes`.
 *
 * Block B of a delayed-hit cache counts the secondary requests merged onto an
 * in-flight fetch, and an EXACT solver has to enumerate it, so it needs a
 * truncation level: the closed population bounds it where there is one, the
 * state-space cutoff where the model is open. `after_event_cache` reads the same
 * field and REFUSES a merge past the level rather than landing in a state the
 * enumeration does not hold -- which is why the two must be set together.
 *
 * @return false when no cache carries a retrieval system, leaving `out` untouched
 */
template <class T>
bool set_retrieval_truncation(const NetworkStruct<T>& sn, const std::vector<std::size_t>& cutoff,
                              NetworkStruct<T>& out) {
    bool any = false;
    for (typename std::map<std::size_t, qn::CacheParam<T>>::const_iterator ci =
             sn.nodeparam.begin();
         ci != sn.nodeparam.end(); ++ci)
        if (ci->second.retrieval_capacity > 0 && !ci->second.retrieval_classes.empty()) any = true;
    if (!any) return false;
    long lvl = 0;
    if (sn.nclosedjobs() > 0) {
        lvl = static_cast<long>(sn.nclosedjobs()) - 1;
    } else {
        std::size_t mx = 0;
        for (std::size_t k = 0; k < cutoff.size(); ++k) mx = std::max(mx, cutoff[k]);
        lvl = static_cast<long>(mx) - 1;
    }
    if (lvl < 0) lvl = 0;
    out = sn;
    for (typename std::map<std::size_t, qn::CacheParam<T>>::iterator ci = out.nodeparam.begin();
         ci != out.nodeparam.end(); ++ci)
        if (ci->second.retrieval_capacity > 0) ci->second.max_pending_retrieval = lvl;
    return true;
}

/**
 * Port of `solver_ctmc.m:47-62`: a G-network SIGNAL NEVER RESIDES at a station,
 * so its per-class capacity there is zero. A REPLY signal is exempt -- it
 * completes a synchronous call and queues like an ordinary job.
 *
 * WHY THE LATTICE CANNOT DISCOVER THIS ON ITS OWN, and why leaving it undiscovered
 * is not merely slow. `space_generator` walks every marginal admitted by
 * `classcap` and, at an FCFS-family station, enumerates the ORDERED buffer of each
 * one, so an uncapped signal class multiplies the buffer alphabet: a cutoff of c
 * over two classes gives on the order of 2^c sequences where the model has c+1
 * states. Every one of them is unreachable -- the arrival handler annihilates the
 * signal instead of seating it -- so the answer is right and only the cost is
 * wrong, which is exactly what makes it dangerous: it is invisible in the results.
 *
 * IT ALSO REPAIRS THE MEMORY GATE, which is the part that turns this from a
 * slowdown into a kill. `ctmc_state_space_logsize` already drops signal classes
 * from its buffered set (ctmc_state_space_logsize.h:176), so without this the gate
 * sizes the SMALL space and clears a generation of the large one. Measured on the
 * `test_ctmc_signal_util` G-network at the cutoff of 30 its tests ask for, the
 * unfixed enumeration was OOM-killed at 43.6 GB resident.
 *
 * @return false when the model declares no annihilating signal, leaving `out` untouched
 */
template <class T>
bool annihilate_signal_capacity(const NetworkStruct<T>& sn, NetworkStruct<T>& out) {
    std::vector<bool> annihilated(sn.nclasses, false);
    bool any = false;
    for (std::size_t k = 0; k < sn.nclasses && k < sn.issignal.size(); ++k) {
        if (!sn.issignal[k]) continue;
        if (k < sn.signaltype.size() && sn.signaltype[k] == lang::SignalType::REPLY) continue;
        annihilated[k] = true;
        any = true;
    }
    if (!any) return false;
    out = sn;
    for (std::size_t i = 0; i < out.stations.size() && i < out.classcap.size(); ++i) {
        // The Source is where a signal is BORN, so its capacity there is what the
        // arrival stream is drawn from and must not be zeroed.
        if (out.stations[i].sched == SchedStrategy::EXT) continue;
        for (std::size_t k = 0; k < annihilated.size() && k < out.classcap[i].size(); ++k)
            if (annihilated[k]) out.classcap[i][k] = 0.0;
    }
    return true;
}

/**
 * Port of the cache write-back of `solver_ctmc_analyzer.m:333-450`.
 *
 * WHAT A CACHE'S HIT AND MISS SHARES ARE, ON A CHAIN. A read leaves the cache in
 * its configured hit class or its miss class and in no other, so the two
 * departure rates of those classes at the cache node ARE the hit and the miss
 * flows, and normalizing by their sum gives the shares. Nothing here reads the
 * routing matrix, whose cache entries `refresh_routing` resolved to a uniform
 * split that is not the answer.
 *
 * THE DELAYED HIT IS A TRANSITION REWARD, NOT A STATE REWARD, and that is the
 * whole reason this is not two lines. A merged secondary request departs in the
 * HIT class, so the hit-class rate above is (true hits + delayed hits) and
 * something has to split it. A fetch of item i completes on exactly the
 * transitions that clear block A bit i, and each such transition releases the
 * block-B count of item i, so the delayed rate is the pi-weighted sum of
 * (count held) x (rate out on those transitions). The alternative identity
 * lambda_i * phi_i is only PASTA-exact and would be wrong on a closed model.
 *
 * `phi` and `d1` are ordinary state rewards of the same two blocks and give the
 * per-item DelayedHitQLen columns; `latency` is left empty, as the reference
 * leaves `actualresidt` NaN.
 */
template <class T>
solvers::CacheMetrics<T> cache_metrics(const NetworkStruct<T>& sn, const CtmcResult<T>& r,
                                       const std::vector<T>& pi) {
    solvers::CacheMetrics<T> out =
        solvers::cache_metrics_of(sn, std::vector<T>(), std::vector<T>(), std::vector<T>(),
                                  std::vector<T>(), Matrix<T>(), Matrix<T>(), std::vector<T>());
    if (out.caches.empty() || r.space.empty()) return out;
    const std::size_t K = sn.nclasses, ns = r.space.size();
    const double dnan = std::numeric_limits<double>::quiet_NaN();

    for (std::size_t c = 0; c < out.caches.size(); ++c) {
        solvers::CacheNodeMetrics<T>& m = out.caches[c];
        const std::size_t isf = sn.stateful_index(m.node);
        if (isf == 0) continue;
        const qn::CacheParam<T>& cp = sn.nodeparam.find(m.node)->second;

        std::vector<double> tn(K, 0.0);
        for (std::size_t s = 0; s < ns && s < pi.size(); ++s) {
            const double p = num_traits<T>::to_double(pi[s]);
            if (p == 0) continue;
            if (isf - 1 >= r.dep_rates[s].size()) continue;
            for (std::size_t k = 0; k < K && k < r.dep_rates[s][isf - 1].size(); ++k)
                tn[k] += p * num_traits<T>::to_double(r.dep_rates[s][isf - 1][k]);
        }

        const std::size_t n = cp.nitems;
        std::size_t tcc = 0;
        for (std::size_t l = 0; l < cp.itemcap.size(); ++l)
            if (cp.itemcap[l] > 0) tcc += static_cast<std::size_t>(cp.itemcap[l]);
        std::vector<std::size_t> rcl, rci, rco;
        if (cp.retrieval_capacity > 0) qn::cache_retrieval_class_map(cp, rcl, rci, rco);
        const std::size_t lvw = r.space[0].local[isf - 1].size();
        const bool retr = !rcl.empty() && lvw >= tcc + n + rcl.size();
        std::vector<double> drate(K, 0.0);

        // TIME-STATIONARY per-item occupancy of each list. The contents block of the
        // local vector holds the item index resident in each cache position, so
        // P(item i is held by list l) is a state reward of pi. This is the
        // TIME-WEIGHTED law, the counterpart of the EMBEDDED (per-request) one the
        // NC/MVA cache algorithms return; the two coincide only under PASTA.
        // The local row is [per-class counts | contents | block A | block B], so the
        // contents block sits at the offset the trailing retrieval blocks leave.
        const std::size_t tail = cp.retrieval_capacity > 0 ? n + rcl.size() : 0;
        if (n > 0 && tcc > 0 && lvw >= tcc + tail) {
            const std::size_t coff = lvw - (tcc + tail);
            const std::size_t h = cp.itemcap.size();
            Matrix<T> ip(n, h + 1);
            std::vector<double> acc(n * h, 0.0);
            std::vector<char> present(n, 0);
            for (std::size_t s = 0; s < ns && s < pi.size(); ++s) {
                const double p = num_traits<T>::to_double(pi[s]);
                if (p == 0) continue;
                const std::vector<T>& lv = r.space[s].local[isf - 1];
                std::size_t off = coff;
                for (std::size_t l = 0; l < h; ++l) {
                    const std::size_t cap =
                        cp.itemcap[l] > 0 ? static_cast<std::size_t>(cp.itemcap[l]) : 0;
                    std::fill(present.begin(), present.end(), 0);
                    for (std::size_t q = 0; q < cap; ++q) {
                        const long it = static_cast<long>(num_traits<T>::to_double(lv[off + q]));
                        if (it >= 1 && static_cast<std::size_t>(it) <= n)
                            present[static_cast<std::size_t>(it) - 1] = 1;
                    }
                    for (std::size_t i = 0; i < n; ++i)
                        if (present[i]) acc[i * h + l] += p;
                    off += cap;
                }
            }
            for (std::size_t i = 0; i < n; ++i) {
                double miss = 1.0;
                for (std::size_t l = 0; l < h; ++l) {
                    const double v = acc[i * h + l];
                    ip(i, l + 1) = num_traits<T>::from_double(v);
                    miss -= v;
                }
                ip(i, 0) = num_traits<T>::from_double(miss);
            }
            m.itemprob = ip;
        }

        if (retr) {
            const std::size_t aoff = lvw - (n + rcl.size()), boff = aoff + n;
            std::vector<double> phi(n, 0.0), d1(n, 0.0);
            for (std::size_t s = 0; s < ns && s < pi.size(); ++s) {
                const double p = num_traits<T>::to_double(pi[s]);
                if (p == 0) continue;
                const std::vector<T>& lv = r.space[s].local[isf - 1];
                for (std::size_t i = 0; i < n; ++i)
                    if (num_traits<T>::to_double(lv[aoff + i]) != 0) phi[i] += p;
                for (std::size_t j = 0; j < rcl.size(); ++j)
                    d1[rci[j] - 1] += p * num_traits<T>::to_double(lv[boff + j]);
            }
            m.delayedhitqlen.assign(n, num_traits<T>::from_int(0));
            m.delayedhitqlenfull.assign(n, num_traits<T>::from_int(0));
            for (std::size_t i = 0; i < n; ++i) {
                m.delayedhitqlen[i] = num_traits<T>::from_double(d1[i]);
                m.delayedhitqlenfull[i] = num_traits<T>::from_double(d1[i] + phi[i]);
            }
            for (std::size_t s = 0; s < ns && s < pi.size(); ++s) {
                const double p = num_traits<T>::to_double(pi[s]);
                if (p == 0) continue;
                const std::vector<T>& lv = r.space[s].local[isf - 1];
                for (std::size_t j = 0; j < rcl.size(); ++j) {
                    const std::size_t i = rci[j] - 1;
                    const double held = num_traits<T>::to_double(lv[boff + j]);
                    if (held <= 0 || num_traits<T>::to_double(lv[aoff + i]) == 0) continue;
                    double completes = 0.0;
                    for (std::size_t b = 0; b < ns; ++b) {
                        if (b == s) continue;
                        const double q = num_traits<T>::to_double(r.Q(s, b));
                        if (q == 0) continue;
                        if (num_traits<T>::to_double(r.space[b].local[isf - 1][aoff + i]) == 0)
                            completes += q;
                    }
                    if (rco[j] - 1 < K) drate[rco[j] - 1] += p * held * completes;
                }
            }
        }

        m.hitprob.assign(K, num_traits<T>::from_double(dnan));
        m.missprob.assign(K, num_traits<T>::from_double(dnan));
        if (retr) m.delayedprob.assign(K, num_traits<T>::from_double(dnan));
        for (std::size_t k = 0; k < K; ++k) {
            if (k >= cp.hitclass.size() || k >= cp.missclass.size()) continue;
            const std::size_t h = cp.hitclass[k], mi = cp.missclass[k];
            if (h == 0 || mi == 0 || h > K || mi > K) continue;
            const double denom = tn[h - 1] + tn[mi - 1];
            if (denom <= 0) continue;
            // Capped at the hit-class flow: the two are the same measurement of
            // the same transitions, and a floating-point excess would report a
            // negative true-hit share rather than a zero one.
            const double d = retr ? std::min(drate[k], tn[h - 1]) : 0.0;
            m.hitprob[k] = num_traits<T>::from_double((tn[h - 1] - d) / denom);
            m.missprob[k] = num_traits<T>::from_double(tn[mi - 1] / denom);
            if (retr) m.delayedprob[k] = num_traits<T>::from_double(d / denom);
        }
    }
    return out;
}

/**
 * Strip the reference's leading infinite-population column from a DECLARED
 * Source row, so a state written by MATLAB, the JAR or native python names a
 * row this port's enumerated space actually holds.
 *
 * The reference writes a Source as `[Inf | one phase block per class]` --
 * `[Inf 1 0 0 0 0 0]` for a Source generating the first of six classes -- and
 * this port's `from_marginal_core`, which the state-space walk runs, writes the
 * phase blocks alone. The two therefore differ by exactly one column, and an
 * untranslated row is one wider than every row of the space, so no lookup can
 * match it and the padding rule cannot help: it left-pads a NARROW candidate
 * and refuses a wide one. `rewardModel_mm1k` came back from its M2C row as
 * `solver_ctmc_transient_analyzer: the initial state is not contained in the
 * state space` while the JAR and native-python rows, each reading its own
 * encoding, passed.
 *
 * The phase block is KEPT rather than rebuilt: a Source with a MAP arrival has
 * several rows in its local space and which one the model declared is the
 * caller's statement, not a detail to regenerate.
 */
template <class T>
void strip_reference_source_column(const NetworkStruct<T>& sn, std::size_t ind,
                                   std::vector<T>& row) {
    const std::size_t ist = sn.nodes[ind - 1].station;
    if (ist == 0 || ist > sn.stations.size()) return;
    if (sn.stations[ist - 1].nodetype != NodeType::Source &&
        sn.stations[ist - 1].sched != SchedStrategy::EXT)
        return;
    std::size_t w = 0;
    for (std::size_t r = 0; r < sn.nclasses; ++r) w += sn.phasessz_of(ist, r + 1);
    // Width is the discriminator rather than the leading value: an exact
    // arithmetic has no infinity, so the reference's marker arrives as the
    // MaxInt clamp `to_marginal` uses and cannot be tested for.
    if (row.size() == w + 1) row.erase(row.begin());
}

/**
 * The model's default initial state: every closed class's jobs at its reference
 * station, everything else empty, which is `Network.initDefault`.
 *
 * The first row `from_marginal_node` emits for that marginal is taken. It is
 * used for two things -- seeding the reachable walk, and choosing a connected
 * component -- and for both any row with the right marginal serves, since rows
 * sharing a marginal are mutually reachable by construction.
 *
 * @return false when some node admits no state at all for that marginal
 */
template <class T>
bool default_init_state(const NetworkStruct<T>& sn, NetState<T>& init) {
    const std::size_t R = sn.nclasses;
    const std::vector<std::size_t>& sfn = sn.stateful_nodes;
    init.local.assign(sfn.size(), std::vector<T>());
    for (std::size_t f = 0; f < sfn.size(); ++f) {
        const std::size_t ind = sfn[f];
        // A DECLARED STATE WINS OVER THE DEFAULT MARKING, which is the whole
        // point of declaring one: `setState` is used precisely to move the jobs
        // off their reference stations, and rebuilding the default here would
        // answer `getProbSysAggr` for a state the caller did not name. The pair
        // (statespace, stateprior) is what the writers emit and the reader
        // stores; the FIRST row is the state when the prior is the trivial one,
        // and a genuine distribution over several rows is not a single initial
        // state at all, so only the one-row case is taken.
        const typename std::map<std::size_t, Matrix<T>>::const_iterator sp =
            sn.statespace.find(ind);
        if (sp != sn.statespace.end() && sp->second.rows() == 1 && sp->second.cols() > 0) {
            std::vector<T> row(sp->second.cols());
            for (std::size_t c = 0; c < sp->second.cols(); ++c) row[c] = sp->second(0, c);
            strip_reference_source_column(sn, ind, row);
            // A ROW OF EXACTLY ONE ENTRY PER CLASS IS A MARGINAL, NOT AN
            // ENCODED ROW, and the two differ wherever a class is DISABLED at
            // the station. The reference keeps a column for every class
            // (`sn.phases` is 1 at a class it does not serve), while this port
            // drops the block outright (`phases_of` is 0), so a Delay serving
            // one of three classes has a width-3 row there and a width-1 row
            // here, LEFT-PADDED to the encoding width. Installing the
            // reference's row verbatim then puts the jobs in another class's
            // block: on tut06_cache_lru_zipf the single client job landed where
            // nothing could serve it and the sample path deadlocked at once,
            // which is what a state nobody can leave looks like. Rebuilding it
            // through `from_marginal_node_first` is the same construction the
            // default branch below uses, so the two agree by construction.
            std::vector<std::size_t> marg(R, 0), mph(R, 1);
            const std::size_t mist = sn.nodes[ind - 1].station;
            bool rebuilt = false;
            if (row.size() == R && mist != 0) {
                for (std::size_t r = 0; r < R; ++r) {
                    const double v = num_traits<T>::to_double(row[r]);
                    marg[r] = v > 0.0 ? static_cast<std::size_t>(v + 0.5) : 0;
                    mph[r] = sn.phasessz_of(mist, r + 1);
                }
                std::vector<T> built;
                if (from_marginal_node_first(sn, ind, marg, mph, built) && !built.empty()) {
                    init.local[f] = built;
                    rebuilt = true;
                }
            }
            if (!rebuilt) init.local[f] = row;
            continue;
        }
        const std::size_t ist = sn.nodes[ind - 1].station;
        std::vector<std::size_t> nmarg(R, 0), ph(R, 1);
        if (ist != 0) {
            for (std::size_t r = 0; r < R; ++r) ph[r] = sn.phasessz_of(ist, r + 1);
            if (sn.stations[ist - 1].nodetype == NodeType::Source) {
                // The infinite reservoir: one job per class in service, and only
                // for the classes the Source generates -- the same guard
                // `space_generator` needs. Without it this returns false on any
                // model with a downstream-only class, the initial state is
                // reported as absent, and the reducible-component selection
                // silently falls back to "largest component": a plausible answer
                // for the wrong chain.
                for (std::size_t r = 0; r < R; ++r)
                    if (!sn.disabled[ist - 1][r]) nmarg[r] = 1;
            } else {
                for (std::size_t r = 0; r < R; ++r) {
                    const double nj = sn.njobs()[r];
                    if (std::isfinite(nj) && sn.classes[r].refstat == ist)
                        nmarg[r] = static_cast<std::size_t>(nj);
                }
                // A PLACE'S DECLARED MARKING IS THE INITIAL STATE, and it is not
                // derivable from a class population: an OPEN class has none, so
                // every Place of an open SPN started empty however many tokens
                // the model put there. On `spn_open_sevenplaces` that dropped
                // P1's two tokens and P5's one; T4 needs a P5 token to fire, so
                // the marking never moved past P3 and the chain reported Tput 0
                // with P3 pinned at the cutoff. `state_initial_occupancy` is the
                // same accessor `space_capacity_c` sizes the lattice with, so the
                // two agree on what the model declared.
                for (std::size_t r = 0; r < R; ++r) {
                    const std::size_t m0 = qn::state_initial_occupancy(sn, ind, r);
                    if (m0 > nmarg[r]) nmarg[r] = m0;
                }
            }
        }
        if (!from_marginal_node_first(sn, ind, nmarg, ph, init.local[f])) return false;
    }
    return true;
}

/** The index in `space` of the default initial state, or npos. */
template <class T>
std::size_t init_state_index(const NetworkStruct<T>& sn, const std::vector<NetState<T>>& space) {
    const std::size_t npos = static_cast<std::size_t>(-1);
    if (space.empty()) return npos;
    NetState<T> init;
    if (!default_init_state(sn, init)) return npos;
    for (std::size_t f = 0; f < init.local.size(); ++f) {
        // Every state of a node carries the node's widest row, left-padded with
        // zeros; a candidate built at the natural width would never match.
        const std::size_t w = space[0].local[f].size();
        if (init.local[f].size() > w) return npos;
        if (init.local[f].size() < w)
            init.local[f].insert(init.local[f].begin(), w - init.local[f].size(),
                                 num_traits<T>::from_int(0));
    }
    const std::vector<double> key = ctmc_detail::state_key(init);
    for (std::size_t s = 0; s < space.size(); ++s)
        if (ctmc_detail::state_key(space[s]) == key) return s;
    return npos;
}

/**
 * The initial DISTRIBUTION over `space`: the product of the declared per-node
 * priors, or a point mass on the default initial state where none is declared.
 *
 * Port of the initial-state loop of `@@SolverCTMC/runAnalyzer.m`, which walks
 * the cartesian product of the per-node state spaces, weights each combination
 * by the product of its per-node priors, integrates the forward equation once
 * per combination and SUMS the trajectories with those weights. Every quantity
 * the transient analyzer reports is linear in pi(t), and pi(t) is linear in
 * pi(0), so seeding the mixture and integrating ONCE gives the same answer at a
 * fraction of the cost -- and on ONE time grid, where the reference has to
 * interpolate its separate adaptive grids onto their union to add them.
 *
 * A combination the enumerated space does not contain is an error rather than a
 * dropped term: it means the declared space and the reachable one disagree, and
 * renormalizing over what is left would answer for a different prior.
 *
 * @return false when a node admits no state at all, as `default_init_state` does
 */
template <class T>
bool init_state_distribution(const NetworkStruct<T>& sn, const std::vector<NetState<T>>& space,
                             std::vector<T>& pi0) {
    const std::size_t npos = static_cast<std::size_t>(-1);
    pi0.assign(space.size(), num_traits<T>::from_int(0));
    if (space.empty()) return false;
    NetState<T> base;
    if (!default_init_state(sn, base)) return false;
    const std::vector<std::size_t>& sfn = sn.stateful_nodes;

    // Per stateful node: the rows it may start in and their probabilities. A
    // node with no declared prior contributes its default row alone.
    std::vector<std::vector<std::vector<T>>> rows(sfn.size());
    std::vector<std::vector<double>> wts(sfn.size());
    for (std::size_t f = 0; f < sfn.size(); ++f) {
        const typename std::map<std::size_t, Matrix<T>>::const_iterator ss =
            sn.statespace.find(sfn[f]);
        const typename std::map<std::size_t, std::vector<T>>::const_iterator sp =
            sn.stateprior.find(sfn[f]);
        if (ss == sn.statespace.end() || sp == sn.stateprior.end() || ss->second.rows() == 0 ||
            ss->second.rows() != sp->second.size()) {
            rows[f].push_back(base.local[f]);
            wts[f].push_back(1.0);
            continue;
        }
        for (std::size_t r = 0; r < ss->second.rows(); ++r) {
            const double w = num_traits<T>::to_double(sp->second[r]);
            if (!(w > 0.0)) continue;  // a zero-probability row is not a state to visit
            std::vector<T> row(ss->second.cols());
            for (std::size_t c = 0; c < ss->second.cols(); ++c) row[c] = ss->second(r, c);
            strip_reference_source_column(sn, sfn[f], row);
            rows[f].push_back(row);
            wts[f].push_back(w);
        }
        if (rows[f].empty()) return false;
    }

    std::vector<std::size_t> pick(sfn.size(), 0);
    double total = 0.0;
    for (;;) {
        NetState<T> cand = base;
        double w = 1.0;
        for (std::size_t f = 0; f < sfn.size(); ++f) {
            cand.local[f] = rows[f][pick[f]];
            w *= wts[f][pick[f]];
            // Every state of a node carries the node's widest row, left-padded
            // with zeros; a candidate at the natural width would never match.
            const std::size_t wd = space[0].local[f].size();
            if (cand.local[f].size() > wd) return false;
            if (cand.local[f].size() < wd)
                cand.local[f].insert(cand.local[f].begin(), wd - cand.local[f].size(),
                                     num_traits<T>::from_int(0));
        }
        std::size_t idx = npos;
        const std::vector<double> key = ctmc_detail::state_key(cand);
        for (std::size_t s = 0; s < space.size() && idx == npos; ++s)
            if (ctmc_detail::state_key(space[s]) == key) idx = s;
        if (idx == npos) return false;
        pi0[idx] = T(pi0[idx] + num_traits<T>::from_double(w));
        total += w;

        std::size_t f = 0;
        for (; f < sfn.size(); ++f) {
            if (++pick[f] < rows[f].size()) break;
            pick[f] = 0;
        }
        if (f == sfn.size()) break;
    }
    if (!(total > 0.0)) return false;
    // The declared priors are per node and need not multiply to one; the
    // reference sums the weighted trajectories without renormalizing, so a
    // prior that already sums to one is unchanged and one that does not is
    // reported as the reference reports it.
    return true;
}

/** Weakly connected components of the generator, as a per-state label. */
template <class T>
std::vector<std::size_t> weak_components(const Matrix<T>& Q, std::size_t& ncomp) {
    const std::size_t n = Q.rows();
    const std::size_t npos = static_cast<std::size_t>(-1);
    std::vector<std::size_t> comp(n, npos);
    ncomp = 0;
    for (std::size_t s = 0; s < n; ++s) {
        if (comp[s] != npos) continue;
        std::vector<std::size_t> stack{s};
        comp[s] = ncomp;
        while (!stack.empty()) {
            const std::size_t u = stack.back();
            stack.pop_back();
            for (std::size_t v = 0; v < n; ++v) {
                if (comp[v] != npos) continue;
                // The DIAGONAL is minus the row sum and is nonzero at almost
                // every state, so it must not be read as an edge to itself.
                if (v == u) continue;
                if (num_traits<T>::to_double(Q(u, v)) != 0 ||
                    num_traits<T>::to_double(Q(v, u)) != 0) {
                    comp[v] = comp[u];
                    stack.push_back(v);
                }
            }
        }
        ++ncomp;
    }
    return comp;
}

/** Restrict a CtmcResult to a subset of its states, keeping their order. */
template <class T>
CtmcResult<T> restrict_to(const CtmcResult<T>& r, const std::vector<std::size_t>& wset) {
    CtmcResult<T> out;
    out.Q = Matrix<T>(wset.size(), wset.size(), num_traits<T>::from_int(0));
    for (std::size_t a = 0; a < wset.size(); ++a)
        for (std::size_t b = 0; b < wset.size(); ++b) out.Q(a, b) = r.Q(wset[a], wset[b]);
    out.space.reserve(wset.size());
    out.arv_rates.reserve(wset.size());
    out.dep_rates.reserve(wset.size());
    for (std::size_t a = 0; a < wset.size(); ++a) {
        out.space.push_back(r.space[wset[a]]);
        out.arv_rates.push_back(r.arv_rates[wset[a]]);
        out.dep_rates.push_back(r.dep_rates[wset[a]]);
    }
    // The filtration is indexed by the SAME rows as Q, so it has to be
    // restricted with it or every CDF built from it would index the wrong
    // states. Its diagonal is NOT rebuilt: a filtration entry is a rate the
    // event contributed, not a generator, and it has no diagonal to speak of.
    out.filt.reserve(r.filt.size());
    for (std::size_t f = 0; f < r.filt.size(); ++f) {
        Matrix<T> F(wset.size(), wset.size(), num_traits<T>::from_int(0));
        for (std::size_t a = 0; a < wset.size(); ++a)
            for (std::size_t b = 0; b < wset.size(); ++b) F(a, b) = r.filt[f](wset[a], wset[b]);
        out.filt.push_back(F);
    }
    // The submatrix of a generator is not a generator: the rates that left the
    // component have to come off the diagonal, or the rows no longer sum to zero
    // and the stationary solve is of a matrix that is not a chain.
    make_infgen(out.Q);
    return out;
}

}  // namespace analyzer_detail

/**
 * Build the chain of ONE struct, solve it, reduce it.
 *
 * The state space is the FULL encoding space, as the reference's
 * `options.config.state_space_gen` default of 'default' asks for, so the
 * reducible handling below is on the normal path and not an edge case. Two
 * exceptions force the reachable walk instead: a global synchronization (an SPN)
 * and a fork firing list.
 *
 * Split out of `solver_ctmc_analyzer` so that the fork-join wrapper can run it on
 * the AUGMENTED struct without the gates firing twice.
 */
template <class T>
CtmcSolution<T> solve_struct(const NetworkStruct<T>& sn_in, const CtmcOptions& opt,
                             const std::vector<qn::FjSync<T>>& fjsync) {
    CtmcSolution<T> out;
    out.cutoff = analyzer_detail::resolve_cutoff(sn_in, opt);
    // A delayed-hit cache is enumerated under a declared merge truncation, and
    // the events are held to the same one; every other model walks `sn_in`.
    NetworkStruct<T> sn_trunc;
    const bool retr = analyzer_detail::set_retrieval_truncation(sn_in, out.cutoff, sn_trunc);
    const NetworkStruct<T>& sn_r = retr ? sn_trunc : sn_in;
    // A signal class holds no place in a station's buffer, and saying so BEFORE
    // the memory gate is what keeps the gate's estimate and the generator's
    // enumeration measuring the same space; see `annihilate_signal_capacity`.
    NetworkStruct<T> sn_sig;
    const bool sig = analyzer_detail::annihilate_signal_capacity(sn_r, sn_sig);
    const NetworkStruct<T>& sn = sig ? sn_sig : sn_r;

    // THE MEMORY PRE-GATE, `@SolverCTMC/runAnalyzer.m:150-164`. It runs BEFORE
    // any generation because `state_max` is enforced inside the generator loop,
    // so on its own it enumerates up to three million states before refusing --
    // and three million is a fixed count that knows nothing about how much
    // memory this host actually has. `mdd` and the `cftp` pair never reach here;
    // see the exemption note in ctmc_memory_gate.h, which is load-bearing.
    {
        mc::CtmcSizeOptions szopt;
        szopt.cutoff = opt.cutoff;
        const double log_nstates = mc::ctmc_state_space_logsize(sn, szopt);
        const mc::CtmcGateResult gate =
            mc::ctmc_memory_gate(log_nstates, opt.force, opt.memory_safety_fraction);
        if (!gate.ok) throw UnsupportedError(gate.message + " Stopping SolverCTMC.");
    }

    const std::vector<Sync<T>> sync = refresh_sync(sn);
    const std::vector<qn::GlobalSync<T>> gsync = refresh_global_sync(sn);

    // AN SPN NEEDS THE REACHABLE WALK, not the lattice enumeration. A firing
    // does not conserve the per-chain population, and a Transition's state is
    // per-MODE, so no population marginal produces the states in which a mode is
    // firing: `from_marginal_node` emits only the all-idle row. Enumerating the
    // lattice therefore yields a space every ENABLE lands outside of, and the
    // generator comes out empty. The reference forces
    // `options.config.state_space_gen = 'reachable'` for exactly this class of
    // model; the condition here is the presence of a global synchronization,
    // which is the same set.
    // A FORK FIRING LIST forces the same choice as an SPN, for the same reason:
    // the lattice enumeration walks the per-chain population, and a firing does
    // not conserve it -- one parent becomes B siblings in classes the lattice
    // gives population 0. The reference sets
    // `options.config.state_space_gen = 'reachable'` on every fork-join model.
    std::vector<NetState<T>> space;
    if (!gsync.empty() || !fjsync.empty()) {
        NetState<T> init;
        if (!analyzer_detail::default_init_state(sn, init))
            throw UnsupportedError(
                "SolverCTMC: the model's initial marking admits no state; check the Place "
                "populations against the class reference stations");
        // THE CUTOFF TRAVELS WITH THE WALK, or an open model never terminates:
        // the Source keeps producing and the walk runs to `state_max` instead of
        // answering. It is the same bound the lattice arm below applies.
        space = reachable_space_generator(sn, init, sync, gsync, opt.state_max, fjsync,
                                          out.cutoff, opt.cutoff_mat);
    } else {
        space = space_generator(sn, out.cutoff, opt.state_max, opt.cutoff_mat);
    }
    if (space.empty())
        throw UnsupportedError(
            "SolverCTMC: the state space is empty; no state satisfies the model's capacities");
    // A DROP region censors the chain: the states it forbids are simply never
    // occupied, so removing them and letting `make_infgen` re-close the rows IS
    // the censored chain. WAITQ and the blocking rules augment the state
    // instead and are refused by name inside.
    space = ctmc_filter_regions(sn, space);
    CtmcResult<T> r = solver_ctmc(sn, space, sync, gsync, opt.keep_filtration, fjsync);
    // ELIMINATE THE VANISHING STATES, `solver_ctmc.m:812`. Until this ran the
    // port returned the UNREDUCED chain: a Router or Fork pass-through, a firable
    // Join and an immediate SPN mode all kept their GlobalConstants::Immediate
    // row, so `-a states` and `-a gen` reported a state space the other three
    // codebases had already complemented away -- six states against four on
    // `fj_tiny_closed`, the two extras carrying 6.1e-09 of the mass each. `-a avg`
    // agreed to the digits printed for exactly that reason, which is why the gap
    // stayed invisible unless the chain itself was asked for.
    ctmc_eliminate_vanishing(r);

    const std::size_t npos = static_cast<std::size_t>(-1);
    // The initial state is looked up in the REDUCED space: the enumeration index
    // no longer addresses the same rows once the vanishing ones are gone.
    std::size_t init = analyzer_detail::init_state_index(sn, r.space);

    std::size_t ncomp = 0;
    const std::vector<std::size_t> comp = analyzer_detail::weak_components(r.Q, ncomp);
    if (ncomp > 1) {
        std::size_t pick;
        if (init != npos) {
            pick = comp[init];
        } else {
            // The initial state is not representable in this space -- an SPN
            // whose ENABLE states were folded away, for one -- so the reference
            // keeps the largest component instead of failing.
            std::vector<std::size_t> sizes(ncomp, 0);
            for (std::size_t s = 0; s < comp.size(); ++s) ++sizes[comp[s]];
            pick = static_cast<std::size_t>(
                std::max_element(sizes.begin(), sizes.end()) - sizes.begin());
        }
        std::vector<std::size_t> wset;
        for (std::size_t s = 0; s < comp.size(); ++s)
            if (comp[s] == pick) wset.push_back(s);
        // The seed moves WITH the restriction: keeping the pre-restriction row
        // would point the block decomposition at whatever state now sits there.
        std::size_t moved = npos;
        if (init != npos)
            for (std::size_t k = 0; k < wset.size(); ++k)
                if (wset[k] == init) {
                    moved = k;
                    break;
                }
        init = moved;
        r = analyzer_detail::restrict_to(r, wset);
    }

    // Every stationary solve goes through the block decomposition, as
    // `solver_ctmc_analyzer.m` does: the irreducible case is the degenerate one
    // BSCC / no transient states, so no dispatch can disagree with the
    // algorithm about whether the chain is reducible. `ctmc_solve` splits a
    // reducible generator into WEAK components and renormalizes across them,
    // which is not the answer the declared initial state selects.
    line::util::LineConsole::step("infinitesimal generator built: %zu states",
                                  static_cast<std::size_t>(r.Q.rows()));
    line::util::LineConsole::step("solving for the stationary distribution");
    const CtmcStationaryResult<T> st = ctmc_stationary(r.Q, init);
    line::util::LineConsole::step("stationary distribution obtained, computing the mean metrics");
    out.pi = st.pi;
    out.warning = st.warning;
    // `gpu` ran on the CPU here, exactly as it does on a MATLAB host with no
    // GPU. Say so, and do not lose a warning the stationary solve already left.
    {
        const std::string fb = method_fallback_warning(opt.method);
        if (!fb.empty()) out.warning = out.warning.empty() ? fb : out.warning + " " + fb;
    }
    out.avg = solver_ctmc_avg_from_pi(sn, r, out.pi);
    out.cache = analyzer_detail::cache_metrics(sn, r, out.pi);
    out.chain = r;
    out.actualmethod = opt.method;
    return out;
}


/**
 * Port of `solver_ctmc_analyzer.m` plus the fork-join wrapper of
 * `@@SolverCTMC/runAnalyzer.m`.
 *
 * A fork-join model is solved on the TAG-AUGMENTED copy and the sibling classes
 * are folded back at the end. The chain and the stationary vector returned are the
 * AUGMENTED ones, as the reference's `result.space` is: they are indexed by a
 * class set the caller did not declare, which is why `fjclassmap` comes back with
 * them rather than being discarded.
 */
template <class T>
CtmcSolution<T> solver_ctmc_analyzer(const NetworkStruct<T>& sn_in, const CtmcOptions& opt) {
    check_method(opt.method);

    // `@@SolverCTMC/runAnalyzer.m:126` converts the non-Markovian service laws
    // to a Markovian surrogate on its own copy of the struct. The conversion
    // runs before the support check, because what the check must see is the
    // struct the generator will actually be built from.
    //
    // THE PH FIT IS FORCED HERE, where the reference takes its CME default.
    // MATLAB's CTMC assembles a RATIONAL generator and can carry a matrix
    // exponential; this port assembles an ordinary one, and an ME's D0 has
    // off-diagonal entries that are not rates, so the chain it builds is not a
    // Markov chain. Measured on a closed Delay + Gamma(4, 0.5) FCFS model with
    // two jobs: the CME fit gives Util 1.075 -- IMPOSSIBLE for one server --
    // and a throughput of 0.538 at the delay against 0.400 at the queue, i.e.
    // flow that does not balance around a cycle. The PH fit gives Util 0.950,
    // balanced throughput 0.475 and QLen 1.525, against MATLAB's 0.954, 0.477
    // and 1.523.
    //
    // The residual 0.4% IS the divergence from the reference, and it is the
    // price of the honest fit: MATLAB matches two moments exactly with its ME
    // (scv 0.25) where the Bernstein PH lands at scv 0.323. Closing it needs
    // the CTMC to handle a rational generator, which is a separate change.
    NetworkStruct<T> converted;
    const NetworkStruct<T>* snp = &sn_in;
    if constexpr (num_traits<T>::has_transcendental) {
        if (api::sn_has_nonmarkov(sn_in, false)) {
            converted = sn_in;
            api::NonmarkovOptions no;
            no.order = opt.nonmkv_order;
            no.phfit = api::PhFit::Ph;
            api::sn_nonmarkov_toph(converted, no);
            snp = &converted;
        }
    }
    const NetworkStruct<T>& sn = *snp;

    ctmc_check_support(sn);
    // runAnalyzerChecks' universal feature gate, AFTER ctmc_check_support so the
    // declaration-defect messages, which the declared set also withholds, still
    // win over the gate's generic one.
    qn::feature_gate("SolverCTMC", qn::ctmc_feature_set(opt.method), sn);

    if (!tr::has_fork_join(sn))
        return solve_struct(sn, opt, std::vector<qn::FjSync<T>>());

    const qn::FjTagged<T> fjt = qn::fj_tag(sn);
    CtmcSolution<T> out = solve_struct(fjt.V, opt, fjt.fjsync);
    tr::fj_foldback(sn, out.avg, fjt.fjclassmap, fjt.korig);
    out.fjclassmap = fjt.fjclassmap;
    return out;
}

/**
 * Port of `@@SolverCTMC/runAnalyzer.m`'s result assembly: solve, then apply the
 * metric filter `@@NetworkSolver/getAvg` puts between the analyzer and the
 * caller, so the table is the same shape SolverMVA and SolverNC print.
 *
 * The response-time mask the product-form runners apply -- drop a metric whose
 * response time is below the tolerance -- is NOT applied here. A CTMC reports
 * what the chain does, and a station a class genuinely visits with a tiny
 * response time is a real measurement rather than a numerical artefact of a
 * fixed point that did not converge there.
 */
template <class T>
mva::AvgResult<T> solver_ctmc_avg_table(const NetworkStruct<T>& sn, const CtmcSolution<T>& d,
                                        const std::string& method) {
    const std::size_t M = sn.nstations, K = sn.nclasses;

    std::vector<std::vector<bool>> srcmask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        if (sn.stations[i].nodetype == NodeType::Source)
            for (std::size_t k = 0; k < K; ++k) srcmask[i][k] = true;

    mva::AvgResult<T> o;
    o.QN = mva::filter_metric(sn, d.avg.QN, mva::MetricKind::QLen, nullptr);
    o.UN = mva::filter_metric(sn, d.avg.UN, mva::MetricKind::Util, nullptr);
    o.RN = mva::filter_metric(sn, d.avg.RN, mva::MetricKind::RespT, nullptr);
    o.TN = mva::filter_metric(sn, d.avg.TN, mva::MetricKind::Tput, nullptr);
    o.WN = mva::filter_metric(sn, mva::sn_get_residt_from_respt(sn, o.RN),
                              mva::MetricKind::ResidT, nullptr);
    o.AN = mva::filter_metric(sn, mva::sn_get_arvr_from_tput(sn, o.TN), mva::MetricKind::ArvR,
                              &srcmask);
    o.CN = d.avg.CN;
    o.XN = d.avg.XN;
    o.cache = d.cache;
    o.method = method;
    o.actualmethod = d.actualmethod;
    o.iter = 1;
    return o;
}

template <class T>
mva::AvgResult<T> solver_ctmc_run_analyzer(const NetworkStruct<T>& sn, const CtmcOptions& opt);

/**
 * Solve the CHAIN-AGGREGATED model and map its metrics back to the classes.
 *
 * `api::sn_aggregate_chains` collapses every chain onto a single class, class
 * switching disappearing with it, and `mva::sn_deaggregate_chain_results` maps
 * chain-level metrics back through alpha, the per-station share of the chain's
 * visits each class carries. Both transforms existed in all four codebases with
 * no solver consumer; this is that consumer.
 *
 * WHAT IS TRADED. Exactness on a non-product-form model: one aggregate service
 * law, fitted to the alpha-weighted first two moments, replaces the per-class
 * ones. On a product-form model the chain IS the unit MVA and convolution
 * already solve in, so the answer is exact and the state space is the smaller
 * one.
 */
template <class T>
mva::AvgResult<T> solver_ctmc_chain_aggregation(const NetworkStruct<T>& sn,
                                                const CtmcOptions& opt) {
    // Driven by tr::transform_solve_chains, so the aggregate is solved through
    // an inner-solve seam rather than a hard-wired call to this analyzer.
    // Clearing the flag states that the aggregate must not be re-aggregated,
    // rather than relying on the caller's nchains < nclasses guard to decline.
    CtmcOptions sub = opt;
    sub.chain_aggregation = false;
    return tr::transform_solve_chains<T>(
        sn,
        [&sub](const NetworkStruct<T>& subsn) { return solver_ctmc_run_analyzer(subsn, sub); },
        opt.method);
}

/**
 * Solve by LOAD CONCEALMENT, the iterated transformation.
 *
 * `tr::transform_solve_lc` sweeps the chains in Gauss-Seidel order, solving each
 * concealed single-chain struct with THIS analyzer. Clearing the flag on the
 * inner options states that a subproblem must not be concealed again, rather
 * than relying on its one-chain shape to decline.
 */
template <class T>
mva::AvgResult<T> solver_ctmc_load_concealment(const NetworkStruct<T>& sn,
                                               const CtmcOptions& opt) {
    CtmcOptions sub = opt;
    sub.load_concealment = false;
    return tr::transform_solve_lc<T>(
        sn,
        [&sub](const NetworkStruct<T>& subsn) { return solver_ctmc_run_analyzer(subsn, sub); },
        opt.method, opt.transform_iter_max);
}

/**
 * Solve with a station subset replaced by a FLOW-EQUIVALENT SERVER, then
 * recover the collapsed stations' own metrics by conditioning.
 *
 * `api::fes_aggregate` has existed in all four codebases with no solver
 * consumer at all: it was exercised by examples and tests only, so nothing in
 * the solver stack depended on it. Flow-equivalent aggregation is the standard
 * route to HIERARCHICAL DECOMPOSITION -- a subnetwork is solved in isolation
 * and enters the outer chain as a single load-dependent station, which is what
 * makes an otherwise intractable state space tractable.
 *
 * The reduced model answers for the surviving stations directly. For a
 * collapsed station the answer is the Chandy-Herzog-Woo conditional sum
 * E[Q_i] = sum_n P(N_fes = n) * Q_i(n), with P read off the reduced chain's
 * stationary law and Q_i(n) from the isolated subnetwork
 * (`fes::fes_compute_metrics`). Throughput needs no conditioning: flow is fixed
 * by the routing and an exact reduction leaves the chain throughput unchanged.
 *
 * EXACT when the collapsed subnetwork is product-form, which is the condition
 * `fes_aggregate` already imposes; an approximation otherwise, and the
 * state-space saving is the reason to accept that.
 */
template <class T>
mva::AvgResult<T> solver_ctmc_fes_aggregation(const NetworkStruct<T>& sn,
                                              const CtmcOptions& opt) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const std::vector<std::size_t>& subset = opt.fes_stations;
    if (subset.size() < 2)
        throw InputError(
            "options.config.fes_stations must name at least two stations: collapsing one "
            "station into a flow-equivalent server saves nothing.");
    if (subset.size() >= M)
        throw InputError(
            "options.config.fes_stations names every station: there is no complement left "
            "to solve.");
    for (std::size_t i : subset)
        if (i < 1 || i > M)
            throw InputError("options.config.fes_stations must be 1-based station indices in 1.." +
                             std::to_string(M) + ".");

    fes::FesAggregateResult<T> agg = fes::fes_aggregate(sn, subset);
    const fes::FesDeaggInfo<T>& info = agg.deagg;
    const NetworkStruct<T>& snRed = agg.model.get_struct();

    CtmcOptions sub = opt;
    sub.fes_stations.clear();
    const CtmcSolution<T> red = solver_ctmc_analyzer(snRed, sub);
    const mva::AvgResult<T> redAvg = solver_ctmc_avg_table(snRed, red, sub.method);

    // P(N_fes = n). The aggregate state space carries K columns per STATION, at
    // (ist-1)*K + k, so the FES's block is the one at its station index.
    const Matrix<T> SSq = ctmc_state_space_aggr(snRed, red.chain.space);
    const std::size_t fesIst = snRed.nodes[info.fesNode - 1].station;
    std::size_t tableSize = 1;
    for (int c : info.cutoffs) tableSize *= static_cast<std::size_t>(c + 1);
    std::vector<double> Pn(tableSize, 0.0);
    for (std::size_t r = 0; r < SSq.rows(); ++r) {
        std::vector<int> nvec(K, 0);
        for (std::size_t k = 0; k < K; ++k)
            nvec[k] = static_cast<int>(
                std::llround(num_traits<T>::to_double(SSq(r, (fesIst - 1) * K + k))));
        Pn[fes::ljd_linearize(nvec, info.cutoffs) - 1] +=
            r < red.pi.size() ? num_traits<T>::to_double(red.pi[r]) : 0.0;
    }

    const fes::FesConditionalMetrics<T> cm = fes::fes_compute_metrics(
        info.isolatedDemands, info.isolatedServers, info.isolatedIsDelay, info.cutoffs);

    const T zero = num_traits<T>::from_int(0);
    mva::AvgResult<T> out;
    out.QN = Matrix<T>(M, K, zero);
    out.UN = Matrix<T>(M, K, zero);
    out.RN = Matrix<T>(M, K, zero);
    out.TN = Matrix<T>(M, K, zero);

    for (std::size_t a = 0; a < info.complementIndices.size(); ++a) {
        const std::size_t i = info.complementIndices[a] - 1;
        for (std::size_t k = 0; k < K; ++k) {
            out.QN(i, k) = redAvg.QN(a, k);
            out.UN(i, k) = redAvg.UN(a, k);
            out.TN(i, k) = redAvg.TN(a, k);
        }
    }

    const std::size_t Msub = info.subsetIndices.size();
    Matrix<T> Qsub(Msub, K, zero), Usub(Msub, K, zero);
    for (std::size_t idx = 0; idx < tableSize; ++idx) {
        if (!(Pn[idx] > 0.0)) continue;
        const T w = num_traits<T>::from_double(Pn[idx]);
        for (std::size_t a = 0; a < Msub && a < cm.QN[idx].rows(); ++a)
            for (std::size_t k = 0; k < K; ++k) {
                Qsub(a, k) += T(w * cm.QN[idx](a, k));
                Usub(a, k) += T(w * cm.UN[idx](a, k));
            }
    }
    for (std::size_t a = 0; a < Msub; ++a) {
        const std::size_t i = info.subsetIndices[a] - 1;
        for (std::size_t k = 0; k < K; ++k) {
            out.QN(i, k) = Qsub(a, k);
            out.UN(i, k) = Usub(a, k);
            // Flow through a station is fixed by the routing, so it is the FES's
            // throughput scaled by the ratio of ORIGINAL visit ratios.
            T ratio = zero;
            for (std::size_t c = 0; c < sn.nchains; ++c) {
                const std::size_t isf = sn.stateful_of_station(i + 1) - 1;
                const std::size_t isfFes = snRed.stateful_of_station(fesIst) - 1;
                if (isf < sn.visits[c].rows() && isfFes < snRed.visits[c].rows() &&
                    snRed.visits[c](isfFes, k) > zero)
                    ratio += T(sn.visits[c](isf, k) / snRed.visits[c](isfFes, k));
            }
            out.TN(i, k) = T(redAvg.TN(fesIst - 1, k) * ratio);
        }
    }

    out.CN.assign(K, zero);
    out.XN = redAvg.XN;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            if (out.TN(i, k) > zero) out.RN(i, k) = T(out.QN(i, k) / out.TN(i, k));
            out.CN[k] += out.RN(i, k);
        }
    out.AN = mva::sn_get_arvr_from_tput(sn, out.TN);
    out.WN = mva::sn_get_residt_from_respt(sn, out.RN);
    out.method = opt.method;
    out.actualmethod = opt.method + "/fes";
    return out;
}

/** Solve and format in one call, for a caller with no use for the chain. */
template <class T>
mva::AvgResult<T> solver_ctmc_run_analyzer(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    // Chain aggregation, opt-in and only where the transform is not the identity.
    if (opt.load_concealment) return solver_ctmc_load_concealment(sn, opt);
    if (opt.chain_aggregation && sn.nchains < sn.nclasses)
        return solver_ctmc_chain_aggregation(sn, opt);
    // Flow-equivalent server aggregation, opt-in.
    if (!opt.fes_stations.empty()) return solver_ctmc_fes_aggregation(sn, opt);
    return solver_ctmc_avg_table(sn, solver_ctmc_analyzer(sn, opt), opt.method);
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_ANALYZER_H
