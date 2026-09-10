/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SSA_SOLVER_SSA_SERIAL_H
#define LINE_SOLVERS_SSA_SOLVER_SSA_SERIAL_H

/**
 * SolverSSA, the `serial` method: a port of `solver_ssa_reachability.m`, of the
 * run loop of `solver_ssa.m`, and of `solver_ssa_analyzer_serial.m`.
 *
 * WHAT THE METHOD IS, AND HOW IT DIFFERS FROM THE NRM. The NRM rewrites the
 * network as a reaction grid over (node, class, phase) counts and never touches
 * the state ENCODING. The serial engine simulates the network in its own
 * encoding instead: at every step it applies the SAME event handlers the CTMC
 * generator applies (`after_event`, `after_global_event`) to the current network
 * state, collects every enabled synchronization with its rate, and takes one
 * Gillespie direct-method step. That is why it reaches models the NRM refuses --
 * anything the encoding can express, the handlers can move -- and why it is
 * slower: the enabled set is rebuilt from scratch at every firing, as the
 * reference's own inlined `solver_ssa_findenabled` does.
 *
 * TWO ENGINES, TWO STREAMS. A seed-fixed result from the serial engine cannot be
 * reproduced by the NRM and vice versa: they consume different draws in a
 * different order. Neither can be reproduced by MATLAB, the JAR or native
 * Python, for the reason `ssa_types.h` states at length. A cross-codebase check
 * against this engine is STATISTICAL. A simulated number without its sample
 * count and its seed is not a measurement, which is why `SsaSerialSolution`
 * carries both and why the tests compare against the exact CTMC only inside a
 * stated Monte Carlo band.
 *
 * THE JUMP CHAIN IS DRAWN IN ONE STEP, NOT TWO. MATLAB calls `State.afterEvent`
 * with `isSimulation = true`, which SAMPLES one successor row and returns its
 * probability; the direct method then picks among the sampled rows. The C++
 * `after_event` is the enumeration-mode handler and returns EVERY successor with
 * its probability, so this port flattens the (synchronization, active row,
 * passive row) triples into one weighted list and draws from it once. The two
 * are the same jump chain -- P(triple) = rate * p_active * p_route * p_passive,
 * normalized -- reached with a different number of uniforms, which is exactly
 * the stream difference above and not a modelling difference.
 *
 * ZERO-RATE ROWS ARE DROPPED, NOT FLOORED. The reference rewrites a zero or NaN
 * rate to 1e-38 "so that it is never selected", which leaves it in the enabled
 * list and in the arrival/departure rate statistics. Dropping it is the same
 * sample path to within 1e-38 and keeps `depRates` exactly the rate the CTMC
 * generator would accumulate for the same state, which is what makes the
 * throughput comparable between the two solvers.
 *
 * THE STATE ROW DOES NOT GROW HERE, SO IT STARTS WIDE. MATLAB's simulation-mode
 * handlers widen a buffer when a job arrives and no slot is free, and
 * `solver_ssa_reachability` left-pads the stored spaces to match. The C++
 * handlers derive the buffer width from the row they are handed, so a path
 * seeded at the natural width of the EMPTY marginal would silently saturate at
 * one waiting job per station. The initial state is therefore padded to the
 * WIDEST row each stateful node's local encoding admits (`serial_detail::
 * max_row_width`), which is the width `space_generator` gives every state of
 * that node, and the path then lives inside the enumerated encoding for free.
 *
 * DOUBLE ONLY, for the reason `solver_ssa_nrm.h` gives: the sample path is
 * generated from exponential clocks drawn as `-log(u)`, there is no exact value
 * to compute, and the answer's error is the Monte Carlo error rather than the
 * rounding. A non-`double` backend is refused BY NAME.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/lang/qn/fj_tag.h"
#include "line/solvers/tr/fj_tag_transform.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/lang/qn/state_events.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_fcr.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/ssa/ssa_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ssa {

/**
 * The serial engine's knobs: `SsaOptions` plus the three the serial path reads
 * and the NRM has no use for.
 *
 * `cutoff` bounds an open class exactly as SolverCTMC's does, and for the same
 * reason: it fixes how wide a station's buffer encoding is, so it decides where
 * the truncation sits. A refused arrival at the truncation boundary is a LOSS
 * (`arrival_is_lost` on an open class), which is what the CTMC truncation does
 * too, so the two solvers truncate the same model the same way.
 */
struct SsaSerialOptions : SsaOptions {
    double cutoff = -1.0;             ///< < 0 = the reference's automatic value
    std::size_t state_max = 3000000;  ///< refuse a reachable space larger than this
};

/**
 * Port of `solver_ssa_reachability.m`'s return: `[SSq, SSh, sn.space]`.
 *
 * `node_space[i]` is the reference's `space{i}`, the distinct local rows node
 * `i` was seen in; `hash` is `SSh`, one 1-BASED index into `node_space[i]` per
 * stateful node per state; `ssq` is `SSq`, the same states with their local rows
 * concatenated. The three are redundant by construction and the reference
 * returns all three because its callers index states by node (`SSh`) and read
 * them flat (`SSq`).
 *
 * THE ORDER IS THE WALK'S, NOT THE REFERENCE'S. The reference pushes and pops a
 * stack of its own; this reuses `reachable_space_generator`, whose stack order
 * differs. The SET is the same and nothing downstream indexes it positionally
 * across codebases, so the difference is not observable in a metric.
 */
template <class T>
struct SsaReachability {
    std::vector<qn::NetState<T>> space;                   ///< the reachable states
    std::vector<std::vector<std::vector<T>>> node_space;  ///< `sn.space`, per stateful node
    std::vector<std::vector<std::size_t>> hash;           ///< `SSh`, 1-based per node
    Matrix<T> ssq;                                        ///< `SSq`, states x concatenated width
};

/** One sample path, in the shape `solver_ssa.m` returns it. */
template <class T>
struct SsaSerialRun {
    /** The DISTINCT states visited, in first-visit order (the reference's `u`). */
    std::vector<qn::NetState<T>> space;
    /**
     * The region token FIFOs of each of those states, the reference's `fcrBuf`.
     *
     * Empty (one empty vector per region, or no vectors at all) on every model
     * without a WAITQ region. A parked job is in NO station, so it appears in no
     * queue length and is visible only here -- the JMT report convention, which
     * `ctmc_waitq_parked` states for the exact solver and
     * `SsaSerialSolution::parked` for this one.
     */
    std::vector<std::vector<std::vector<std::size_t>>> buf;
    /** `pi`: the fraction of simulated time spent in each of them. */
    std::vector<double> pi;
    /** `SSq`: the per-(station, class) job counts of each distinct state. */
    Matrix<T> ssq;
    /**
     * `arvRates` / `depRates`, indexed [distinct state][stateful-1][class-1].
     *
     * They are a deterministic function of the state, so one sample per state is
     * the exact value and not an estimate -- which is what lets the analyzer
     * multiply them by `pi` and get a throughput rather than a sample mean. The
     * reference says as much where it keeps `arvRatesSamples(ui(s),...)`.
     */
    std::vector<std::vector<std::vector<double>>> arv_rates, dep_rates;
    /**
     * The DERIVED rates per state, laid out like `arv_rates`: how fast the
     * transitions enabled in that state START a class-r service at a stateful
     * node, and how fast they PUSH a class-r job in service back into the
     * buffer there. Annotations on the arcs the engine already walks, so no
     * rate, probability or state depends on them.
     */
    std::vector<std::vector<std::vector<double>>> start_rates, preempt_rates;
    /**
     * The rate of the cache MERGE transitions, i.e. of the delayed hits, in the
     * same [state][stateful-1][class-1] shape and by the same argument.
     *
     * A delayed hit is invisible in `dep_rates`: the merged request is released
     * later, in the HIT class, and is indistinguishable there from a true hit.
     * The merge itself is the only cache transition that EMPTIES the node -- it
     * decrements the read class and adds nothing -- which is what identifies it.
     */
    std::vector<std::vector<std::vector<double>>> dly_rates;
    /** `tranSysState{1}`: the cumulative time at each firing. */
    std::vector<double> tran_time;
    /** `tranSync`: which synchronization fired, `sync.size() + g` for a global one. */
    std::vector<std::size_t> tran_sync;
    /**
     * The row of `space` the path OCCUPIED over `[t-dt, t]`, one per firing.
     *
     * The trace and the distinct-state table are two views of the same path and
     * the samplers need both: `sampleSys` prints the state at each event and
     * `getProb` sums the time spent in one state, so keeping only `pi` would
     * lose the order and keeping only the rows would lose the aggregation. It
     * indexes the state BEFORE the firing, exactly as `pi` weights it.
     */
    std::vector<std::size_t> tran_state;
    double simulated_time = 0.0;
    std::size_t samples = 0;      ///< firings actually performed
    std::size_t warmup = 0;       ///< leading firings excluded from `pi`
    unsigned long seed = 0;       ///< the stream this path came from
};

/** What the cache write-back of `solver_ssa_analyzer_serial.m` produces. */
struct SsaCacheRatio {
    std::size_t node = 0;                 ///< 1-based Cache node index
    std::vector<double> hitprob, missprob;  ///< per class, NaN where undefined
    /**
     * The delayed-hit share, EMPTY off a retrieval system.
     *
     * `hitprob` carries the TRUE hits alone once this is filled: the three shares
     * partition every read, hit + delayed + miss = 1, which is the convention NC
     * and LDES report and what makes `ArvR = arvr*(missprob + delayedprob)`.
     */
    std::vector<double> delayedprob;
    /**
     * `actualresidt`: NaN, and NOT a port gap. The reference warns
     * "Retrieval-system expected latency is not currently implemented; reporting
     * NaN" and reports NaN in every codebase, so reproducing the NaN IS parity.
     */
    std::vector<double> residt;
};

/** The serial analyzer's return: the metric table, the path, and the stream. */
template <class T>
struct SsaSerialSolution {
    SsaSolution avg;               ///< QN, UN, RN, TN, XN, CN; `method` = "serial"
    SsaSerialRun<T> run;
    unsigned long seed = 0;        ///< carried beside the numbers, never implied
    std::vector<SsaCacheRatio> cache;
    /**
     * `fjclassmap` of the tag augmentation, empty on a model with no Fork.
     *
     * The PATH inside `run` is the augmented one -- its classes are the sibling
     * classes `fj_tag` invented -- while `avg` has been folded back onto the
     * classes the caller declared. The map is what relates the two, and it is
     * returned rather than discarded for the reason SolverCTMC returns it: a
     * caller reading the trajectory needs to know which class a sibling came
     * from.
     */
    std::vector<std::size_t> fjclassmap;
    /**
     * Mean number of jobs parked in a region FIFO, per class of the struct that
     * ran. Zero everywhere without a WAITQ region.
     *
     * IT IS IN NO QLen, so a population check on a closed model has to add it
     * back by hand. That is the JMT convention `ctmc_waitq_parked` reports for
     * the exact solver and not an omission here.
     */
    std::vector<double> parked;
};

namespace serial_detail {

using lang::NodeType;
using lang::SchedStrategy;

/**
 * `solver_ssa.m`'s own guards, plus what this port cannot represent.
 *
 * The reference's remaining guards -- non-exponential reneging patience,
 * non-QUEUE_LENGTH balking, heterogeneous servers (`nodeparam.nservertypes`) --
 * have NO field in the C++ `NetworkStruct`, so a model declaring one cannot be
 * built and the guard would be dead code. Reneging is the sharpest case: it
 * enters the chain only through `refresh_sync`'s `impatience_classes` argument,
 * which no analyzer has a source for, so no RENEGE synchronization exists here
 * at all. The same holds for the class-switch mask violation the reference
 * raises inside its scan: `sn.csmask` is not carried, and the state-dependent
 * routing that can violate it is refused when the struct is built.
 */
template <class T>
bool serial_check(const qn::NetworkStruct<T>& sn, bool raise = true, bool skip_fork = false) {
    // A WAITQ region parks refused jobs in a FIFO the engine carries beside the
    // network state, and the combinations that FIFO has no meaning against are
    // the CTMC's own: one list, so the two solvers refuse the same models with
    // the same sentence and neither can silently accept what the other rejects.
    if (!sn.regions.empty()) {
        // The rule decides which machinery runs, so the gate is the matching
        // one: `ctmc_check_waitq_support` for the token FIFO, and the DROP-only
        // checker otherwise, which is what refuses BAS/BBS/RSRD by name.
        try {
            if (ctmc::ctmc_has_waitq_region(sn))
                ctmc::ctmc_check_waitq_support(sn);
            else
                ctmc::ctmc_check_region_rules(sn);
        } catch (const UnsupportedError&) {
            if (!raise) return false;
            throw;
        }
    }
    // A Fork fires through `fjsync`, which only the TAG-AUGMENTED struct carries
    // (`fj_tag`). A raw fork-join struct reaching the engine would find no
    // firing at all: the fork would never emit, the branches would stay empty
    // and the path would report a network that transparently swallows every
    // forked task. The analyzer augments before it runs, so this refusal is
    // reachable only by a caller who drove the engine directly.
    if (!skip_fork && sn.has_fork() && !sn.isfjaugmented) {
        if (!raise) return false;
        throw UnsupportedError(
            "SolverSSA(method='serial'): the model contains a Fork node and has not been "
            "tag-augmented. A fork fires through `sn.fjsync` / `State.afterFJEvent`, which only "
            "`fj_tag` builds; call `solver_ssa_serial_analyzer`, which augments and folds the "
            "sibling classes back, rather than the engine directly");
    }
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const typename std::map<std::size_t, qn::RetrialParam<T>>::const_iterator rit =
            sn.retrialparam.find(i + 1);
        if (rit == sn.retrialparam.end()) continue;
        bool any = false;
        std::size_t served = 0;
        for (std::size_t r = 0; r < rit->second.retrial_proc.size(); ++r) {
            if (rit->second.retrial_proc[r].disabled) continue;
            any = true;
            if (rit->second.retrial_proc[r].type != lang::ProcessType::EXP) {
                if (!raise) return false;
                throw UnsupportedError(
                    "SolverSSA(method='serial'): station '" + sn.stations[i].name +
                    "' retries with non-exponential patience. SOLVER_SSA supports only "
                    "exponential (memoryless) retrial delay in every codebase, because the orbit "
                    "carries no remaining-delay phase");
            }
            if (r < rit->second.max_attempts.size() && rit->second.max_attempts[r] > 0) {
                if (!raise) return false;
                throw UnsupportedError(
                    "SolverSSA(method='serial'): station '" + sn.stations[i].name +
                    "' declares a finite retrial max-attempts count. SOLVER_SSA supports only "
                    "unlimited retrials in every codebase, because the attempt counter is not "
                    "part of the state");
            }
        }
        if (!any) continue;
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            if (!sn.disabled[i][r] && sn.service[i][r].D0.rows() > 0) ++served;
        if (served > 1) {
            if (!raise) return false;
            throw UnsupportedError(
                "SolverSSA(method='serial'): station '" + sn.stations[i].name +
                "' is a multi-class retrial station. SOLVER_SSA supports retrial only for "
                "single-class stations in every codebase");
        }
    }
    return true;
}

/**
 * Can `solver_ssa_serial_analyzer` run this model at all?
 *
 * The dispatcher's fallback test, and it is the SAME body as the refusal above
 * so the two cannot drift. The fork guard is skipped because the analyzer
 * augments before the engine sees the struct, so a raw fork-join model IS one
 * the serial path runs -- answering otherwise would send it back to the NRM,
 * which excludes fork-join in every codebase.
 */
template <class T>
bool serial_can_run(const qn::NetworkStruct<T>& sn) {
    return serial_check(sn, false, true);
}

/**
 * The widest local row stateful node `ind` admits, which is the width
 * `space_generator` gives EVERY state of that node.
 *
 * WHY IT IS COMPUTED RATHER THAN ENUMERATED. Taking the width from
 * `space_generator` would mean building the whole cartesian product across
 * nodes, which is the cost simulation exists to avoid. The width of one node's
 * row is a property of that node alone: for an ordered (class-tag) buffer it
 * grows with the TOTAL jobs held and not with how they split across classes, and
 * for every other encoding it is constant. So one call to `from_marginal_node`
 * at a maximal admissible marginal settles it, and the split is chosen to fill
 * the classes in order precisely because a lopsided multiset has the fewest
 * permutations for `from_marginal` to enumerate.
 *
 * The descent to a smaller total is not defensive padding: `from_marginal_node`
 * returns NO rows for a marginal the station cannot hold, and the joint bound
 * (`cap` against the sum of the per-class bounds) can be met by a marginal that
 * some other constraint inside the handler still rejects.
 */
template <class T>
std::size_t max_row_width(const qn::NetworkStruct<T>& sn, std::size_t ind,
                          const std::vector<std::size_t>& cutoff) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    std::vector<std::size_t> ph(R, 1);
    if (ist != 0)
        for (std::size_t r = 0; r < R; ++r) ph[r] = sn.phasessz_of(ist, r + 1);

    // A stateful non-station (a Cache, a Join, a Transition) holds a per-class
    // count and its local variables, both of fixed width, so the empty marginal
    // already gives the final width.
    if (ist == 0 || sn.stations[ist - 1].nodetype == NodeType::Source) {
        std::vector<T> row;
        if (!qn::from_marginal_node_first(sn, ind, std::vector<std::size_t>(R, 0), ph, row))
            throw UnsupportedError("SolverSSA(method='serial'): node '" + sn.nodes[ind - 1].name +
                                   "' admits no state at all, so no sample path can start");
        return row.size();
    }

    // The per-class bound: the class population when closed, the cutoff when
    // open, never above the station's own per-class capacity.
    std::vector<std::size_t> bound(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        const double nj = sn.njobs()[r];
        double b = std::isfinite(nj) ? nj : static_cast<double>(cutoff[r]);
        const double cc = sn.classcap[ist - 1][r];
        if (cc < b) b = cc;
        if (!(b > 0)) b = 0;
        bound[r] = static_cast<std::size_t>(b);
    }
    double tcapd = sn.cap[ist - 1];
    std::size_t total = 0;
    for (std::size_t r = 0; r < R; ++r) total += bound[r];
    if (std::isfinite(tcapd) && tcapd < static_cast<double>(total))
        total = static_cast<std::size_t>(tcapd);

    // A PAS / OI station is the one encoding `from_marginal_node` cannot size:
    // its row is the ORDERED JOB LIST, one position per job, and the function
    // builds the ordinary [buffer | server] split instead -- which for a
    // single-server station returns a ONE-COLUMN row however many jobs the
    // marginal holds. A path seeded at that width would hold one job and block
    // every further arrival, so the width is taken from the job bound directly.
    if (sn.stations[ist - 1].sched == SchedStrategy::PAS ||
        sn.stations[ist - 1].sched == SchedStrategy::OI)
        return total + sn.nvars_of(ind);

    for (std::size_t t = total + 1; t-- > 0;) {
        std::vector<std::size_t> n(R, 0);
        std::size_t left = t;
        for (std::size_t r = 0; r < R && left > 0; ++r) {
            n[r] = std::min(left, bound[r]);
            left -= n[r];
        }
        if (left > 0) continue;  // this total does not fit the per-class bounds
        const std::vector<std::vector<T>> rows = qn::from_marginal_node(sn, ind, n, ph);
        if (!rows.empty()) return rows[0].size();
    }
    throw UnsupportedError("SolverSSA(method='serial'): station '" + sn.stations[ist - 1].name +
                           "' admits no state at all, so no sample path can start");
}

/**
 * The initial network state, padded to the encoding width of every node.
 *
 * The marginal is `Network.initDefault`'s -- every closed class's jobs at its
 * reference station -- taken from the CTMC analyzer so the two solvers start the
 * same model in the same state. The LEFT pad is not a convention chosen here: it
 * is what `space_generator` does to the narrow rows, and every slicer in
 * `state_events.h` measures its blocks from the RIGHT-hand end, so a right pad
 * would shift the server block and silently decode the wrong queue.
 */
template <class T>
qn::NetState<T> wide_init_state(const qn::NetworkStruct<T>& sn,
                                const std::vector<std::size_t>& cutoff) {
    qn::NetState<T> init;
    if (!ctmc::analyzer_detail::default_init_state(sn, init))
        throw UnsupportedError(
            "SolverSSA(method='serial'): the model's initial state admits no state; check the "
            "class populations against their reference stations");
    for (std::size_t f = 0; f < sn.stateful_nodes.size(); ++f) {
        const std::size_t w = max_row_width(sn, sn.stateful_nodes[f], cutoff);
        if (init.local[f].size() > w) continue;  // already at or beyond the encoding width
        init.local[f].insert(init.local[f].begin(), w - init.local[f].size(),
                             num_traits<T>::from_int(0));
    }
    return init;
}

}  // namespace serial_detail

/**
 * Port of `solver_ssa_reachability.m`: the states the DYNAMICS can occupy,
 * decomposed per stateful node.
 *
 * The walk itself is `reachable_space_generator`, which applies the same
 * handlers to the same synchronization list and additionally walks `gsync`. The
 * reference's reachability walks only `sync`, so an SPN's firing states reach it
 * through `solver_ssa`'s own global scan instead of through this function; here
 * they are in the space from the start, which is strictly more of the reachable
 * set and never less.
 *
 * WHAT THIS ADDS OVER THE WALK is the decomposition the reference returns and
 * the CTMC does not need: the per-node list of distinct local rows, and the
 * per-state index into it. That is the reference's `space{i}` / `SSh` pair, and
 * it exists so a caller can address a state by node without carrying the rows.
 */
template <class T>
SsaReachability<T> solver_ssa_reachability(const qn::NetworkStruct<T>& sn,
                                           const SsaSerialOptions& opt = SsaSerialOptions()) {
    serial_detail::serial_check(sn);
    ctmc::CtmcOptions copt;
    copt.cutoff = opt.cutoff;
    copt.state_max = opt.state_max;
    const std::vector<std::size_t> cutoff = ctmc::analyzer_detail::resolve_cutoff(sn, copt);
    const std::vector<qn::Sync<T>> sync = qn::refresh_sync(sn);
    const std::vector<qn::GlobalSync<T>> gsync = qn::refresh_global_sync(sn);
    const qn::NetState<T> init = serial_detail::wide_init_state(sn, cutoff);

    SsaReachability<T> out;
    // A WAITQ region's states are AUGMENTED with the token FIFO, which a
    // `NetState` cannot hold, so this decomposition has no row to report them
    // in and refuses rather than returning the region-free walk under the
    // region's name. The metrics do not need it: the engine carries the FIFO
    // itself and `SsaSerialRun::buf` reports it.
    if (ctmc::ctmc_has_waitq_region(sn))
        throw UnsupportedError(
            "solver_ssa_reachability: the model declares a WAITQ finite capacity region, whose "
            "states carry a per-region token FIFO outside every node; this decomposition is "
            "per-node and cannot represent it. Use `solver_ssa_serial_analyzer`, whose run carries "
            "the FIFO, or `solver_ctmc_waitq` for the exact augmented space");
    // The cutoff resolved above bounds the walk, as it bounds SolverCTMC's:
    // without it an open model's Source produces forever and the walk runs to
    // `state_max` rather than answering.
    out.space = ctmc::reachable_space_generator(sn, init, sync, gsync, opt.state_max,
                                                std::vector<qn::FjSync<T>>(), cutoff);
    // A DROP region censors the space exactly as it censors SolverCTMC's: the
    // forbidden states are never occupied, so leaving them in would report a
    // reachable set the dynamics cannot reach.
    out.space = ctmc::ctmc_filter_regions(sn, out.space);

    const std::size_t NF = sn.stateful_nodes.size();
    out.node_space.assign(NF, std::vector<std::vector<T>>());
    out.hash.assign(out.space.size(), std::vector<std::size_t>(NF, 0));
    // One map per node, keyed on the row itself: the reference's `matchrow`,
    // which is a linear scan and turns the decomposition quadratic on a space
    // large enough to be worth walking.
    std::vector<std::map<std::vector<double>, std::size_t>> seen(NF);
    for (std::size_t s = 0; s < out.space.size(); ++s)
        for (std::size_t f = 0; f < NF; ++f) {
            std::vector<double> key(out.space[s].local[f].size(), 0.0);
            for (std::size_t j = 0; j < key.size(); ++j)
                key[j] = num_traits<T>::to_double(out.space[s].local[f][j]);
            const typename std::map<std::vector<double>, std::size_t>::const_iterator it =
                seen[f].find(key);
            if (it != seen[f].end()) {
                out.hash[s][f] = it->second;
            } else {
                out.node_space[f].push_back(out.space[s].local[f]);
                seen[f][key] = out.node_space[f].size();
                out.hash[s][f] = out.node_space[f].size();
            }
        }

    std::size_t width = 0;
    for (std::size_t f = 0; f < NF; ++f)
        width += out.space.empty() ? 0 : out.space[0].local[f].size();
    out.ssq = Matrix<T>(out.space.size(), width, num_traits<T>::from_int(0));
    for (std::size_t s = 0; s < out.space.size(); ++s) {
        std::size_t col = 0;
        for (std::size_t f = 0; f < NF; ++f)
            for (std::size_t j = 0; j < out.space[s].local[f].size(); ++j)
                out.ssq(s, col++) = out.space[s].local[f][j];
    }
    return out;
}

/**
 * The serial engine: the sample path of `solver_ssa.m`'s main loop.
 *
 * It is a class for the reason `NrmEngine` is: the loop threads the current
 * state, the visited-state table, the two rate tables and the trace through
 * every step, and the reference threads the same through one long function.
 */
template <class T>
class SsaSerialEngine {
public:
    /**
     * `fjsync` is the fork firing list of the TAG-AUGMENTED struct, empty for a
     * model with no Fork. It is passed in rather than derived because the
     * augmentation produces the struct and the firing list together and `sn`
     * must be the augmented one: deriving it here would leave the two able to
     * disagree about which classes the branches carry.
     */
    SsaSerialEngine(const qn::NetworkStruct<T>& sn, const SsaSerialOptions& opt,
                    const std::vector<qn::FjSync<T>>& fjsync = std::vector<qn::FjSync<T>>())
        : sn_(sn), opt_(opt), rng_(opt.seed), fjsync_(fjsync) {
        serial_detail::serial_check(sn);
        ctmc::CtmcOptions copt;
        copt.cutoff = opt.cutoff;
        copt.state_max = opt.state_max;
        const std::vector<std::size_t> cutoff = ctmc::analyzer_detail::resolve_cutoff(sn, copt);
        sync_ = qn::refresh_sync(sn);
        gsync_ = qn::refresh_global_sync(sn);
        sdr_ = sn.has_sdr_routing();
        // WHICH REGION MACHINERY, decided exactly as SolverCTMC decides it: a
        // model whose every region class applies DROP is CENSORED, because a
        // refused job is destroyed and the chain simply never occupies the
        // forbidden states; one WAITQ class anywhere puts the whole model on the
        // token-FIFO relation, which handles its DROP classes inline. Choosing
        // differently from the CTMC here would make the simulator and the exact
        // solver answer different models under one model file.
        waitq_ = ctmc::ctmc_has_waitq_region(sn);
        if (!sn.regions.empty()) {
            if (waitq_) {
                caps_ = ctmc::waitq_detail::extract_caps(sn);
                ctmc::waitq_detail::resolve_lmax(sn, cutoff, caps_);
            }  // the DROP-only rule set is checked by `serial_check` above
        }
        init_.net = serial_detail::wide_init_state(sn, cutoff);
        init_.buf.assign(caps_.size(), std::vector<std::size_t>());
        if (!sn.regions.empty() && !region_admissible(init_.net))
            throw InputError(
                "SolverSSA(method='serial'): the model's initial state violates a finite capacity "
                "region; the region cannot hold the model's initial population, so no sample path "
                "can start");
    }

    /** Run `opt.samples` firings and return the path with its statistics. */
    SsaSerialRun<T> run();

    /** The synchronization list the trace's `tran_sync` indexes. */
    const std::vector<qn::Sync<T>>& sync() const { return sync_; }
    /** The state the path starts from, at full encoding width. */
    const qn::NetState<T>& init_state() const { return init_.net; }

private:
    /** One enabled transition: where it goes and what it contributes. */
    struct Move {
        std::size_t sync = 0;   ///< index into `sync_`, or `sync_.size() + g`
        double weight = 0.0;    ///< rate * p_active * p_route * p_passive
        ctmc::WaitqState<T> next;
    };

    const qn::NetworkStruct<T>& sn_;
    SsaSerialOptions opt_;
    SsaRng rng_;
    std::vector<qn::Sync<T>> sync_;
    std::vector<qn::GlobalSync<T>> gsync_;
    std::vector<qn::FjSync<T>> fjsync_;
    std::vector<ctmc::waitq_detail::RegionCaps<T>> caps_;
    bool waitq_ = false;
    /** The routing must be re-evaluated at every state visited (Krzesinski SDR). */
    bool sdr_ = false;
    ctmc::WaitqState<T> init_;

    /** `ctmc_region_admissible` on one network state, the DROP censoring test. */
    bool region_admissible(const qn::NetState<T>& st) const {
        if (sn_.regions.empty()) return true;
        std::vector<qn::NetState<T>> one(1, st);
        const Matrix<T> A = ctmc::ctmc_state_space_aggr(sn_, one);
        std::vector<T> nir(A.cols());
        for (std::size_t c = 0; c < A.cols(); ++c) nir[c] = A(0, c);
        return ctmc::ctmc_region_admissible(sn_, nir);
    }

    /**
     * The reference's inlined `solver_ssa_findenabled`: every synchronization
     * that can fire in `st`, with the arrival and departure rates it carries.
     */
    void enabled(const ctmc::WaitqState<T>& st, std::vector<Move>& moves,
                 std::vector<std::vector<double>>& arv, std::vector<std::vector<double>>& dep,
                 std::vector<std::vector<double>>& dly,
                 std::vector<std::vector<double>>& start,
                 std::vector<std::vector<double>>& preempt) const;

    /**
     * `phi(n)` on the CURRENT sample-path state, as an (nstations*nclasses)
     * row-major vector of scalings. The twin of `ctmc_gd_factor`'s single row.
     */
    std::vector<T> gd_factor_now(const qn::NetState<T>& ns) const {
        const std::size_t M = sn_.stations.size(), K = sn_.nclasses;
        const T zero = num_traits<T>::from_int(0);
        std::vector<T> npop(M * K, zero);
        for (std::size_t ist = 1; ist <= M; ++ist) {
            const std::size_t isf = sn_.stateful_of_station(ist);
            if (isf == 0) continue;
            if (sn_.stations[ist - 1].nodetype == lang::NodeType::Source) continue;
            const std::size_t ind = sn_.node_of_station(ist);
            std::vector<std::size_t> ph(K, 1), shift(K, 0);
            std::size_t w = 0;
            for (std::size_t k = 0; k < K; ++k) {
                ph[k] = sn_.phasessz_of(ist, k + 1);
                shift[k] = w;
                w += ph[k];
            }
            const qn::Marginal<T> m =
                qn::to_marginal(sn_, ist, ns.local[isf - 1], ph, shift, sn_.nvars_of(ind));
            for (std::size_t k = 0; k < K; ++k) npop[(ist - 1) * K + k] = m.nir[k];
        }
        const std::vector<T> v = sn_.gdscaling(npop);
        std::vector<T> out(M * K, num_traits<T>::from_int(1));
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                const T f = v.size() == 1 ? v[0] : (v.size() == M ? v[i] : v[i * K + r]);
                if (!(num_traits<T>::to_double(f) >= 0))
                    throw InputError(
                        "the global dependence handle returned a non-finite or negative scaling");
                out[i * K + r] = f;
            }
        return out;
    }
};

template <class T>
void SsaSerialEngine<T>::enabled(const ctmc::WaitqState<T>& st, std::vector<Move>& moves,
                                 std::vector<std::vector<double>>& arv,
                                 std::vector<std::vector<double>>& dep,
                                 std::vector<std::vector<double>>& dly,
                                 std::vector<std::vector<double>>& start,
                                 std::vector<std::vector<double>>& preempt) const {
    const std::size_t local = sn_.nodes.size() + 1;  // the dummy passive node
    const std::size_t R = sn_.nclasses;
    moves.clear();
    for (std::size_t f = 0; f < arv.size(); ++f)
        for (std::size_t r = 0; r < R; ++r) {
            arv[f][r] = 0.0;
            dep[f][r] = 0.0;
            dly[f][r] = 0.0;
            start[f][r] = 0.0;
            preempt[f][r] = 0.0;
        }

    // A WAITQ REGION REPLACES THE WHOLE ENUMERATION rather than filtering it.
    // The token FIFO is state the network encoding cannot hold, a refused job
    // parks instead of being lost, and every firing runs a release cascade to a
    // fixed point, so there is no per-move filter that turns the ordinary
    // relation into this one. `waitq_successors` IS that relation, shared with
    // the CTMC generator. Its own support gate has already refused an SPN and a
    // fork-join model beside a WAITQ region, so no global or fork firing can
    // reach here on this branch.
    if (waitq_) {
        std::vector<ctmc::waitq_detail::Successor<T>> succ;
        ctmc::waitq_detail::waitq_successors(sn_, sync_, caps_, st, succ);
        for (std::size_t i = 0; i < succ.size(); ++i) {
            const ctmc::waitq_detail::Successor<T>& su = succ[i];
            const double w = num_traits<T>::to_double(su.w);
            if (!(w > 0)) continue;
            Move m;
            m.sync = su.sync;
            m.weight = w;
            m.next = su.next;
            moves.push_back(m);
            if (su.dep_isf != 0) dep[su.dep_isf - 1][su.dep_cls - 1] += w;
            for (std::size_t q = 0; q < su.arv.size(); ++q)
                arv[su.arv[q].first - 1][su.arv[q].second - 1] += w;
        }
        return;
    }

    const qn::NetState<T>& base = st.net;
    // Global (Whittle) rate scaling declared through `set_global_dependence`. It
    // reads the FULL population matrix, so it is a CONSTANT within one state and
    // factors out of the per-transition rates, exactly as in `solver_ctmc`. The
    // CTMC tabulates it once per state of the enumerated space; a simulator has
    // one state at a time, so the table collapses to this single row.
    std::vector<T> gd_now;
    const bool has_gd = static_cast<bool>(sn_.gdscaling);
    if (has_gd) gd_now = gd_factor_now(base);
    // The state-dependent routing table, for the same reason and at the same
    // scope: one evaluation of eq. (10) per state, not one per synchronization.
    // The CTMC tabulates it over the enumerated space; a sample path holds one
    // state at a time, so the table collapses to this single matrix.
    Matrix<T> rt_now;
    if (sdr_) rt_now = qn::rt_state(sn_, base.local);
    for (std::size_t a = 0; a < sync_.size(); ++a) {
        const qn::Sync<T>& sy = sync_[a];
        const std::size_t isf_a = sn_.stateful_index(sy.active.node);
        if (isf_a == 0) continue;  // a stateless node schedules nothing
        const std::size_t isf_p =
            sy.passive.node == local ? 0 : sn_.stateful_index(sy.passive.node);
        if (sy.passive.node != local && isf_p == 0) continue;

        const qn::EventOutcome<T> oa = qn::after_event(sn_, sy.active.node, base.local[isf_a - 1],
                                                       sy.active.event, sy.active.cls);
        // PHASE is scaled too, or phase-type service would advance unscaled
        const bool gd_here = has_gd && sn_.nodes[sy.active.node - 1].station != 0 &&
                             (sy.active.event == lang::EventType::DEP ||
                              sy.active.event == lang::EventType::PHASE);
        const double gd_f =
            gd_here ? num_traits<T>::to_double(
                          gd_now[(sn_.nodes[sy.active.node - 1].station - 1) * R +
                                 (sy.active.cls - 1)])
                    : 1.0;
        // A cache READ whose successor holds one job FEWER is the delayed-hit
        // merge: the request joined an in-flight fetch and is held in block B, so
        // it leaves no departure to count and is released later in the hit class.
        // Every other cache READ keeps the server block flat.
        const bool cache_read = sy.active.event == lang::EventType::READ &&
                                sn_.nodes[sy.active.node - 1].nodetype == lang::NodeType::Cache;
        double srv_pre = 0.0;
        if (cache_read)
            for (std::size_t r = 0; r < R && r < base.local[isf_a - 1].size(); ++r)
                srv_pre += num_traits<T>::to_double(base.local[isf_a - 1][r]);

        double fired = 0.0;  // what this synchronization contributes from here
        for (std::size_t ia = 0; ia < oa.space.size(); ++ia) {
            const double rate = num_traits<T>::to_double(oa.rate[ia]) * gd_f;
            const double pa = num_traits<T>::to_double(oa.prob[ia]);
            if (!(rate > 0) || !(pa > 0)) continue;
            bool merged = false;
            if (cache_read) {
                double srv_post = 0.0;
                for (std::size_t r = 0; r < R && r < oa.space[ia].size(); ++r)
                    srv_post += num_traits<T>::to_double(oa.space[ia][r]);
                merged = srv_post - srv_pre == -1.0;
            }

            if (sy.passive.node == local) {
                Move m;
                m.sync = a;
                m.weight = rate * pa;
                m.next = st;
                m.next.net.local[isf_a - 1] = oa.space[ia];
                // A DROP region CENSORS the chain: a transition into a state the
                // region forbids is not taken at all, which is what deleting the
                // state from the CTMC's space and re-closing its rows amounts
                // to. The move is dropped whole, so it contributes neither a
                // departure nor an arrival, exactly as the deleted column does.
                if (!region_admissible(m.next.net)) continue;
                fired += m.weight;
                if (merged) dly[isf_a - 1][sy.active.cls - 1] += m.weight;
                // The START/PREEMPT tags of this arc, weighted like the rate it
                // carries: they annotate the transition itself. Written for
                // EVERY action, not only for departures -- a retrial or a
                // polling switchover starts service without being a DEP.
                ssa_detail::add_tag_rates(start, preempt, sn_, sy.active.node, oa, ia, m.weight);
                moves.push_back(m);
                continue;
            }
            // A self-loop synchronization reads the passive node's state AFTER
            // the active half has been applied, since they are the same node.
            const std::vector<T>& src =
                sy.passive.node == sy.active.node ? oa.space[ia] : base.local[isf_p - 1];
            const qn::EventOutcome<T> op = qn::after_event(sn_, sy.passive.node, src,
                                                           sy.passive.event, sy.passive.cls);
            // NO ROWS is a BLOCK, not a loss: the destination has no room and
            // cannot take the job, so the upstream departure is disabled and the
            // synchronization simply does not appear in the enabled list. This
            // is the reference's `prob_sync_p = 0`.
            // The routing probability, read at the state the job LEAVES from,
            // which is what `sub_sdr` reads.
            double proute =
                sy.passive.statedep
                    ? num_traits<T>::to_double(rt_now(sy.passive.rt_row, sy.passive.rt_col))
                    : num_traits<T>::to_double(sy.passive.prob);
            // ROUND-ROBIN reads the pointer the ACTIVE node carries once its own
            // departure has advanced it, exactly as the CTMC generator does; the
            // uniform expansion in `sy.passive.prob` would make the dispatcher a
            // coin. The pointer lives in the state, so the serial engine needs
            // no cursor of its own -- unlike the NRM, which walks the arcs.
            if (sy.active.event == lang::EventType::DEP &&
                sn_.rr_var_slot(sy.active.node, sy.active.cls) != 0) {
                const std::size_t w = sn_.nvars_of(sy.active.node);
                const std::vector<T>& arow = oa.space[ia];
                std::size_t dest = 0;
                if (arow.size() >= w) {
                    const std::vector<T> var(arow.end() - w, arow.end());
                    dest = sn_.rr_dest(sy.active.node, sy.active.cls, var);
                }
                proute = (dest == sy.passive.node && sy.passive.cls == sy.active.cls) ? 1.0 : 0.0;
            }
            for (std::size_t ip = 0; ip < op.space.size(); ++ip) {
                const double pp = num_traits<T>::to_double(op.prob[ip]);
                if (!(pp > 0)) continue;
                Move m;
                m.sync = a;
                m.weight = rate * pa * proute * pp;
                if (!(m.weight > 0)) continue;
                m.next = st;
                m.next.net.local[isf_a - 1] = oa.space[ia];
                m.next.net.local[isf_p - 1] = op.space[ip];
                if (!region_admissible(m.next.net)) continue;
                fired += m.weight;
                if (merged) dly[isf_a - 1][sy.active.cls - 1] += m.weight;
                // both halves are tagged: the arrival half is where most
                // service starts happen
                ssa_detail::add_tag_rates(start, preempt, sn_, sy.active.node, oa, ia, m.weight);
                ssa_detail::add_tag_rates(start, preempt, sn_, sy.passive.node, op, ip, m.weight);
                moves.push_back(m);
            }
        }
        // A DEP synchronization is one job LEAVING the active node and ENTERING
        // the passive one, so the same rate is a departure there and an arrival
        // here. The passive half of a LOCAL action is the dummy node, which is
        // nobody's arrival.
        if (sy.active.event == lang::EventType::DEP && fired > 0) {
            dep[isf_a - 1][sy.active.cls - 1] += fired;
            if (isf_p != 0) arv[isf_p - 1][sy.passive.cls - 1] += fired;
        }
    }

    for (std::size_t g = 0; g < gsync_.size(); ++g) {
        const qn::GlobalOutcome<T> go = qn::after_global_event(sn_, base, gsync_[g]);
        for (std::size_t io = 0; io < go.space.size(); ++io) {
            const double w = num_traits<T>::to_double(go.rate[io]) *
                             num_traits<T>::to_double(go.prob[io]);
            if (!(w > 0)) continue;
            Move m;
            m.sync = sync_.size() + g;
            m.weight = w;
            m.next = st;
            m.next.net = go.space[io];
            if (!region_admissible(m.next.net)) continue;
            moves.push_back(m);
            // A FIRE consumes from its PRE places and produces to its POST
            // places, which is a departure and an arrival respectively; an
            // ENABLE only reads the markings and moves nothing.
            if (gsync_[g].active.event != lang::EventType::FIRE) continue;
            for (std::size_t j = 0; j < gsync_[g].passive.size(); ++j) {
                const qn::ModeEvent<T>& pev = gsync_[g].passive[j];
                const std::size_t pisf = sn_.stateful_index(pev.node);
                if (pisf == 0 || pev.cls == 0 || pev.cls > R) continue;
                if (pev.event == lang::EventType::PRE) dep[pisf - 1][pev.cls - 1] += w;
                else if (pev.event == lang::EventType::POST) arv[pisf - 1][pev.cls - 1] += w;
            }
        }
    }

    // FORK FIRINGS, atomic across the fork and every branch head, so like an SPN
    // firing they take the whole network state and cannot be decomposed into
    // sync halves. `refresh_sync` emits no DEP for a Fork, which is why nothing
    // above has already counted them.
    //
    // ONE DEPARTURE, B ARRIVALS. The parent leaves the fork in its own class and
    // one sibling enters each branch head in the tag's auxiliary class, so the
    // rate statistics are accumulated by hand exactly as `solver_ctmc` does:
    // an ordinary synchronization has no way to express a one-to-many emission.
    for (std::size_t k = 0; k < fjsync_.size(); ++k) {
        const qn::FjSync<T>& e = fjsync_[k];
        const std::size_t isf_f = sn_.stateful_index(e.fork);
        if (isf_f == 0) continue;
        const qn::GlobalOutcome<T> fo = qn::after_fj_event(sn_, e, base);
        double fired = 0.0;
        for (std::size_t io = 0; io < fo.space.size(); ++io) {
            const double w = num_traits<T>::to_double(fo.rate[io]) *
                             num_traits<T>::to_double(fo.prob[io]);
            if (!(w > 0)) continue;
            Move m;
            m.sync = sync_.size() + gsync_.size() + k;
            m.weight = w;
            m.next = st;
            m.next.net = fo.space[io];
            if (!region_admissible(m.next.net)) continue;
            fired += w;
            moves.push_back(m);
        }
        if (!(fired > 0)) continue;
        dep[isf_f - 1][e.cls - 1] += fired;
        for (std::size_t b = 0; b < e.branchheads.size(); ++b) {
            const std::size_t isf_b = sn_.stateful_index(e.branchheads[b]);
            if (isf_b != 0) arv[isf_b - 1][e.auxclasses[b] - 1] += fired;
        }
    }
}

template <class T>
SsaSerialRun<T> SsaSerialEngine<T>::run() {
    // The exponential holding time is drawn as -log(u)/lambda, so the backend
    // must have a logarithm at all. The analyzer gates on `double` before it
    // instantiates this, so a caller reaching the assert is one that reached
    // past the gate and would otherwise get a compile error deep inside the
    // uniform draw instead of a sentence naming the reason.
    static_assert(num_traits<T>::has_transcendental,
                  "solver_ssa_serial: an SSA sample path is generated from exponential clocks "
                  "drawn as -log(u)/rate, which needs transcendental arithmetic");

    const std::size_t NF = sn_.stateful_nodes.size();
    const std::size_t R = sn_.nclasses;
    SsaSerialRun<T> out;
    out.seed = opt_.seed;
    out.warmup = static_cast<std::size_t>(
        std::floor(std::max(0.0, std::min(0.99, opt_.warmupfrac)) *
                   static_cast<double>(opt_.samples)));

    ctmc::WaitqState<T> cur = init_;
    std::vector<Move> moves;
    std::vector<std::vector<double>> arv(NF, std::vector<double>(R, 0.0));
    std::vector<std::vector<double>> dep(NF, std::vector<double>(R, 0.0));
    std::vector<std::vector<double>> dly(NF, std::vector<double>(R, 0.0));
    // Derived START/PREEMPT rates, sampled exactly like the three above: the
    // rate at which the transitions enabled in the current state start a
    // class-r service, or push a class-r job in service back into the buffer.
    std::vector<std::vector<double>> start(NF, std::vector<double>(R, 0.0));
    std::vector<std::vector<double>> preempt(NF, std::vector<double>(R, 0.0));
    std::vector<double> weights;
    std::map<std::vector<double>, std::size_t> index;  // state key -> row of `space`

    double cur_time = 0.0;
    out.tran_time.reserve(opt_.samples);
    out.tran_sync.reserve(opt_.samples);
    for (std::size_t n = 0; n < opt_.samples; ++n) {
        enabled(cur, moves, arv, dep, dly, start, preempt);
        if (moves.empty())
            throw NumericError(
                "solver_ssa_serial: the sample path entered a deadlock before collecting all "
                "samples, no synchronization is enabled");

        weights.resize(moves.size());
        double tot = 0.0;
        for (std::size_t i = 0; i < moves.size(); ++i) {
            weights[i] = moves[i].weight;
            tot += weights[i];
        }
        const std::size_t sel = rng_.draw(weights);
        // The transition is drawn BEFORE the holding time, as the reference
        // draws them: the two are independent, so the order changes only which
        // stream this engine is, and being a nameable stream is the point.
        const double dt = -std::log(rng_.uniform()) / tot;

        // The state is recorded with the time spent IN it, so the pair belongs
        // to the state before the firing, not after.
        //
        // THE KEY IS THE AUGMENTED ONE and the stored row is the network half.
        // Two states that agree on every node but differ in a region FIFO are
        // DIFFERENT states -- their enabled sets differ -- so merging them would
        // put two rate rows on one entry and report whichever was seen first.
        // `space` therefore may hold the same network row twice, once per FIFO
        // content, which every consumer here handles: each row carries its own
        // `pi` and the metrics are sums over rows, never lookups by row.
        const std::vector<double> key = ctmc::waitq_detail::waitq_key(cur);
        std::size_t si;
        const typename std::map<std::vector<double>, std::size_t>::const_iterator it =
            index.find(key);
        if (it != index.end()) {
            si = it->second;
        } else {
            si = out.space.size();
            index[key] = si;
            out.space.push_back(cur.net);
            out.buf.push_back(cur.buf);
            out.pi.push_back(0.0);
            out.arv_rates.push_back(arv);
            out.dep_rates.push_back(dep);
            out.dly_rates.push_back(dly);
            out.start_rates.push_back(start);
            out.preempt_rates.push_back(preempt);
        }
        // The warmup discard drops the transient from the TIME AVERAGE only: the
        // states themselves stay in the table, so a state visited only during
        // the transient keeps its (exact) rate row and contributes zero weight.
        if (n >= out.warmup) {
            out.pi[si] += dt;
            out.simulated_time += dt;
        }
        cur_time += dt;
        out.tran_time.push_back(cur_time);
        out.tran_sync.push_back(moves[sel].sync);
        out.tran_state.push_back(si);

        cur = moves[sel].next;
        out.samples = n + 1;
    }

    double tot_pi = 0.0;
    for (std::size_t s = 0; s < out.pi.size(); ++s) tot_pi += out.pi[s];
    if (tot_pi > 0)
        for (std::size_t s = 0; s < out.pi.size(); ++s) out.pi[s] /= tot_pi;
    out.ssq = ctmc::ctmc_state_space_aggr(sn_, out.space);
    return out;
}

namespace serial_detail {

/** `map_mean(PH{ist}{k})`, or a negative sentinel when the pair has no process. */
template <class T>
double service_mean(const qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t k) {
    const lang::Distrib<T>& d = sn.service[ist - 1][k - 1];
    if (d.disabled || d.D0.rows() == 0) return -1.0;
    mam::Map<T> m;
    m.D0 = d.D0;
    m.D1 = d.D1;
    try {
        return num_traits<T>::to_double(mam::map_mean(m));
    } catch (const Error&) {
        return -1.0;  // a zero-rate process has no mean; the reference skips it too
    }
}

}  // namespace serial_detail

/**
 * Port of `solver_ssa_analyzer_serial.m`: run the serial engine and reduce its
 * path to the metric table.
 *
 * THE UTILIZATION ESTIMATOR IS THE REFERENCE'S, discipline by discipline. An
 * INF station is utilized by every job it holds; a PS-family station takes the
 * ARRIVAL rate over rate*servers, because the offered load is what a processor
 * sharing server carries; every other discipline takes the arrival rate times
 * the mean service time over the servers. A class that can be DROPPED -- an open
 * class at a station with a finite capacity -- is measured on the CARRIED rate
 * instead, because the offered rate counts arrivals that never entered service.
 *
 * ONE DIVERGENCE FROM THE REFERENCE, stated rather than hidden:
 *
 *   THE CACHE LOOP IS INDEXED BY NODE, NOT BY STATEFUL INDEX. The reference
 *   writes `sn.nodetype(isf) == NodeType.Cache` with `isf` running over the
 *   STATEFUL nodes, so on any model whose stateful indices differ from its node
 *   indices -- one with a Source, a Router or a ClassSwitch, which is most of
 *   them -- it tests the type of the wrong node. Reproducing that would report
 *   hit ratios for a node that is not the cache.
 *
 * A PAS STATION takes the reference's `otherwise` branch, T*E[S]/c, and NOT the
 * in-service occupancy `solver_ctmc_avg_from_pi` computes for the same station.
 * The two disagree because a pass-and-swap job does not engage a single server;
 * the reference serial analyzer is what is ported here, and the disagreement is
 * named so it is not mistaken for a defect in either.
 */
template <class T>
SsaSerialSolution<T> solver_ssa_serial_on_struct(const qn::NetworkStruct<T>& sn,
                                                 const SsaSerialOptions& opt,
                                                 const std::vector<qn::FjSync<T>>& fjsync) {
    // `if constexpr`, not a run-time test: the engine reaches `map_mean` and the
    // logarithm of a uniform, so a Rational instantiation would fail to COMPILE
    // rather than refuse. The gate has to keep the body from being instantiated.
    if constexpr (!std::is_same<T, double>::value) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_ssa_serial: an SSA sample path is generated from exponential clocks, which "
            "are logarithms of uniform draws; there is no exact value to compute and a wider "
            "float carries no information the Monte Carlo error does not swamp. Rerun with "
            "--arith double");
    } else {
        using lang::SchedStrategy;
        const std::size_t M = sn.nstations, K = sn.nclasses;
        SsaSerialSolution<T> out;
        out.seed = opt.seed;

        SsaSerialEngine<T> eng(sn, opt, fjsync);
        out.run = eng.run();
        const SsaSerialRun<T>& r = out.run;

        // The parked population, per class: a token carries the class its job
        // will enter in, so the mean is the time average of the FIFO contents.
        out.parked.assign(K, 0.0);
        for (std::size_t s = 0; s < r.buf.size() && s < r.pi.size(); ++s)
            for (std::size_t f = 0; f < r.buf[s].size(); ++f)
                for (std::size_t j = 0; j < r.buf[s][f].size(); ++j) {
                    const std::size_t cls = (r.buf[s][f][j] - 1) % K + 1;
                    out.parked[cls - 1] += r.pi[s];
                }

        SsaSolution& a = out.avg;
        a.method = "serial";
        a.samples = r.samples;
        a.simulated_time = r.simulated_time;
        a.QN = Matrix<double>(M, K, 0.0);
        a.UN = Matrix<double>(M, K, 0.0);
        a.RN = Matrix<double>(M, K, 0.0);
        a.TN = Matrix<double>(M, K, 0.0);
        a.XN.assign(K, 0.0);
        a.CN.assign(K, 0.0);
        a.StartN = Matrix<double>(M, K, 0.0);
        a.PreemptN = Matrix<double>(M, K, 0.0);

        // System throughput is the DEPARTURE rate at each class's reference
        // station, which is what makes X a per-class quantity rather than a sum.
        for (std::size_t k = 1; k <= K; ++k) {
            const std::size_t refsf = sn.stateful_of_station(sn.classes[k - 1].refstat);
            if (refsf == 0) continue;
            for (std::size_t s = 0; s < r.space.size(); ++s)
                a.XN[k - 1] += r.pi[s] * r.dep_rates[s][refsf - 1][k - 1];
        }

        // The reference's `isempty(sn.lldscaling) && isempty(sn.cdscaling)` test
        // is over the WHOLE matrix, so one load-dependent station puts every
        // station on the scaling branch. That is kept: at a station with no
        // scaling the branch degenerates to T*E[S]/c, which differs from the
        // unscaled branch only in using the carried rather than the offered
        // rate, and reproducing the reference's table means reproducing that.
        bool scaled = static_cast<bool>(sn.gdscaling);
        for (std::size_t i = 0; i < M; ++i)
            if (!sn.stations[i].lldscaling.empty() || sn.stations[i].cdscaling) scaled = true;

        for (std::size_t ist = 1; ist <= M; ++ist) {
            const std::size_t isf = sn.stateful_of_station(ist);
            if (isf == 0) continue;
            const SchedStrategy sched = sn.stations[ist - 1].sched;
            const double S = sn.stations[ist - 1].nservers;
            for (std::size_t k = 1; k <= K; ++k) {
                for (std::size_t s = 0; s < r.space.size(); ++s) {
                    a.TN(ist - 1, k - 1) += r.pi[s] * r.dep_rates[s][isf - 1][k - 1];
                    a.QN(ist - 1, k - 1) +=
                        r.pi[s] * num_traits<T>::to_double(r.ssq(s, (ist - 1) * K + k - 1));
                    // same time average as TN, over the derived tag rates
                    if (s < r.start_rates.size()) {
                        a.StartN(ist - 1, k - 1) += r.pi[s] * r.start_rates[s][isf - 1][k - 1];
                        a.PreemptN(ist - 1, k - 1) += r.pi[s] * r.preempt_rates[s][isf - 1][k - 1];
                    }
                }
            }

            const bool is_ps = sched == SchedStrategy::PS || sched == SchedStrategy::DPS ||
                               sched == SchedStrategy::GPS || sched == SchedStrategy::LPS;
            // A SOURCE holds no jobs, so QLen, Util and hence RespT are zero
            // there BY DEFINITION and only its throughput is a quantity -- the
            // rule `solver_ssa_nrm.h` states at length and `solver_ctmc_avg_from_pi`
            // applies by the same `continue`. The reference serial analyzer has
            // no such branch and lets a Source fall into `otherwise`, where it
            // divides the arrival rate by a server count `solver_ssa.m` has
            // meanwhile overwritten with the station's capacity; that product is
            // not a utilization of anything.
            if (sn.stations[ist - 1].nodetype == lang::NodeType::Source ||
                sched == SchedStrategy::EXT)
                continue;
            if (sched == SchedStrategy::INF) {
                for (std::size_t k = 1; k <= K; ++k) a.UN(ist - 1, k - 1) = a.QN(ist - 1, k - 1);
                continue;
            }
            if (scaled) {
                // The EFFECTIVE server count a load-dependent station can
                // deliver: `max(c, max_n lld(n))`. A class-dependent station
                // normalizes by its DECLARED peak instead, which is the only
                // thing the utilization can be a fraction of.
                double ceff = S;
                const std::vector<T>& lld = sn.stations[ist - 1].lldscaling;
                for (std::size_t j = 0; j < lld.size(); ++j)
                    ceff = std::max(ceff, num_traits<T>::to_double(lld[j]));
                const bool is_cd = static_cast<bool>(sn.stations[ist - 1].cdscaling);
                const bool is_jd = static_cast<bool>(sn.stations[ist - 1].jdscaling);
                // A global (Whittle) dependence rescales the service rate the
                // same way, so the peak IT declares normalizes Util too.
                const bool is_gd = static_cast<bool>(sn.gdscaling);
                std::vector<T> gdpk;
                if (is_gd)
                    gdpk.assign(sn.gdscalingpeak.begin() + (ist - 1) * K,
                                sn.gdscalingpeak.begin() + ist * K);
                for (std::size_t k = 1; k <= K; ++k) {
                    const double mean = serial_detail::service_mean(sn, ist, k);
                    if (mean < 0) continue;
                    // The divisor is the PRODUCT of the declared peaks when either
                    // dependence is present, and the effective server count only
                    // otherwise: a station carrying both scales its rate by both,
                    // so normalizing by one of them alone leaves the other's factor
                    // in the reported utilization.
                    double cdiv = ceff;
                    if (is_cd || is_jd || is_gd) {
                        cdiv = 1.0;
                        const std::vector<T>* pks[3] = {&sn.stations[ist - 1].cdscalingpeak,
                                                        &sn.stations[ist - 1].jdscalingpeak,
                                                        &gdpk};
                        const char* names[3] = {"setClassDependence", "setJointDependence",
                                                "setGlobalDependence"};
                        const bool on[3] = {is_cd, is_jd, is_gd};
                        for (std::size_t h = 0; h < 3; ++h) {
                            if (!on[h]) continue;
                            const std::vector<T>& pk = *pks[h];
                            if (pk.size() < k || !(num_traits<T>::to_double(pk[k - 1]) > 0))
                                throw InputError(
                                    "SolverSSA(method='serial'): station '" +
                                    sn.stations[ist - 1].name +
                                    "' declares a dependent scaling with no declared peak rate. "
                                    "Utilization there is T*E[S]/peak, so pass the peak to " +
                                    names[h]);
                            cdiv *= num_traits<T>::to_double(pk[k - 1]);
                        }
                    }
                    a.UN(ist - 1, k - 1) = cdiv > 0 ? a.TN(ist - 1, k - 1) * mean / cdiv : 0.0;
                }
                continue;
            }

            // A class whose jobs can be lost here is measured on the carried
            // rate; everything else on the offered rate, which is exact in
            // steady state and is what the reference reports.
            for (std::size_t k = 1; k <= K; ++k) {
                const double mean = serial_detail::service_mean(sn, ist, k);
                if (mean < 0) continue;
                const bool can_drop = !std::isfinite(sn.njobs()[k - 1]) &&
                                      (std::isfinite(sn.cap[ist - 1]) ||
                                       std::isfinite(sn.classcap[ist - 1][k - 1]));
                double arv = 0.0;
                if (!can_drop)
                    for (std::size_t s = 0; s < r.space.size(); ++s)
                        arv += r.pi[s] * r.arv_rates[s][isf - 1][k - 1];
                if (is_ps) {
                    const double mu = num_traits<T>::to_double(sn.rates(ist - 1, k - 1));
                    if (!(mu > 0)) continue;
                    a.UN(ist - 1, k - 1) =
                        (can_drop ? a.TN(ist - 1, k - 1) / mu : arv / mu) / S;
                } else {
                    a.UN(ist - 1, k - 1) =
                        (can_drop ? a.TN(ist - 1, k - 1) * mean : arv * mean) / S;
                }
            }
        }

        // Little's law per station, then the per-class system response time.
        for (std::size_t k = 1; k <= K; ++k) {
            for (std::size_t ist = 1; ist <= M; ++ist)
                a.RN(ist - 1, k - 1) = a.TN(ist - 1, k - 1) > 0
                                           ? a.QN(ist - 1, k - 1) / a.TN(ist - 1, k - 1)
                                           : 0.0;
            // The reference's `CN(k) = NK(k)/XN(k)` with NK the class population:
            // infinite for an open class, which its NaN sweep does NOT clear and
            // which the NRM engine reports identically.
            if (a.XN[k - 1] > 0) a.CN[k - 1] = sn.classes[k - 1].population / a.XN[k - 1];
        }

        // The cache write-back: every read leaves as exactly one of hit or miss,
        // so the two departure streams divide the read rate between them and
        // their ratio is the realized hit probability. The reference stores it
        // into `sn.nodeparam{ind}.actualhitprob`; the struct is const here, so it
        // is returned beside the table instead.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (typename std::map<std::size_t, qn::CacheParam<T>>::const_iterator ci =
                 sn.nodeparam.begin();
             ci != sn.nodeparam.end(); ++ci) {
            const std::size_t ind = ci->first;
            if (ind == 0 || ind > sn.nodes.size()) continue;
            if (sn.nodes[ind - 1].nodetype != lang::NodeType::Cache) continue;
            const std::size_t isf = sn.stateful_index(ind);
            if (isf == 0) continue;
            SsaCacheRatio cr;
            cr.node = ind;
            cr.hitprob.assign(K, nan);
            cr.missprob.assign(K, nan);
            cr.residt.assign(K, nan);
            std::vector<double> dly(K, 0.0);
            bool any_delayed = false;
            for (std::size_t k = 1; k <= K; ++k) {
                if (ci->second.hitclass.size() < k || ci->second.missclass.size() < k) continue;
                const std::size_t h = ci->second.hitclass[k - 1];
                const std::size_t mi = ci->second.missclass[k - 1];
                if (h == 0 || mi == 0 || h > K || mi > K) continue;
                double th = 0.0, tm = 0.0, td = 0.0;
                for (std::size_t s = 0; s < r.space.size(); ++s) {
                    th += r.pi[s] * r.dep_rates[s][isf - 1][h - 1];
                    tm += r.pi[s] * r.dep_rates[s][isf - 1][mi - 1];
                    if (s < r.dly_rates.size()) td += r.pi[s] * r.dly_rates[s][isf - 1][k - 1];
                }
                if (th + tm > 0) {
                    // `th` already carries the released delayed hits, so the
                    // delayed share is CARVED OUT of it rather than added as a
                    // fourth share.
                    cr.hitprob[k - 1] = std::max(th - td, 0.0) / (th + tm);
                    cr.missprob[k - 1] = tm / (th + tm);
                    dly[k - 1] = td / (th + tm);
                    if (td > 0) any_delayed = true;
                }
            }
            if (any_delayed) cr.delayedprob = dly;
            out.cache.push_back(cr);
        }
        return out;
    }
}


/**
 * Port of `solver_ssa_analyzer_serial.m` plus the fork-join wrapper
 * `@@SolverSSA/runAnalyzer.m` puts in front of it.
 *
 * A FORK-JOIN MODEL IS SIMULATED ON THE TAG-AUGMENTED COPY, exactly as
 * SolverCTMC solves it there: the fork emits one sibling per branch in a class
 * of its own, the tag is what lets the Join recognize which siblings belong to
 * the same parent, and `fj_tag` is the only thing that builds the `fjsync`
 * firing list the engine fires. The sample path in the returned run is
 * therefore indexed by the AUGMENTED classes; only the metric table is folded
 * back, which is why `fjclassmap` travels with it.
 */
template <class T>
SsaSerialSolution<T> solver_ssa_serial_analyzer(const qn::NetworkStruct<T>& sn,
                                                const SsaSerialOptions& opt) {
    // The augmentation is skipped entirely on a non-`double` backend so the
    // refusal a caller reads is the arithmetic one, from inside the engine,
    // rather than a fork-join message about a model whose real problem is that
    // an exponential clock has no exact value.
    if constexpr (!std::is_same<T, double>::value) {
        return solver_ssa_serial_on_struct(sn, opt, std::vector<qn::FjSync<T>>());
    } else {
        if (!tr::has_fork_join(sn))
            return solver_ssa_serial_on_struct(sn, opt, std::vector<qn::FjSync<T>>());
        const qn::FjTagged<T> fjt = qn::fj_tag(sn);
        SsaSerialSolution<T> out = solver_ssa_serial_on_struct(fjt.V, opt, fjt.fjsync);
        tr::fj_foldback(sn, out.avg, fjt.fjclassmap, fjt.korig);
        out.fjclassmap = fjt.fjclassmap;
        out.parked.resize(fjt.korig);
        return out;
    }
}

/**
 * The `serial` entry of `solver_ssa_analyzer.m`.
 *
 * The reference reaches it from `default` (when the NRM eligibility gate fails),
 * from `ssa`, from `serial` and from `para`/`parallel` without the Parallel
 * Computing Toolbox. `para`/`parallel` is NOT that: it replicates the SAME
 * engine across workers and averages, so answering it with one replica would
 * report a number at a different variance from the one asked for, and it refuses
 * by name here.
 */
template <class T>
SsaSerialSolution<T> solver_ssa_serial(const qn::NetworkStruct<T>& sn,
                                       const SsaSerialOptions& opt) {
    const std::string& m = opt.method;
    if (m == "default" || m == "ssa" || m == "serial") return solver_ssa_serial_analyzer(sn, opt);
    if (m == "para" || m == "parallel")
        throw UnsupportedError(
            "SolverSSA: the '" + m +
            "' method runs the serial engine on several workers and averages the replicas "
            "(solver_ssa_analyzer_parallel.m). The engine is ported; the replication is not, and "
            "one replica has a different variance from the average of many. Use 'serial'");
    throw UnsupportedError("SolverSSA(serial): '" + m +
                           "' is not a method this entry accepts; it implements 'serial' and the "
                           "'default' and 'ssa' aliases that reach it");
}

}  // namespace ssa
}  // namespace line

#endif  // LINE_SOLVERS_SSA_SOLVER_SSA_SERIAL_H
