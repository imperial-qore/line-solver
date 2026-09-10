/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SSA_SOLVER_SSA_NRM_SPACE_H
#define LINE_SOLVERS_SSA_SOLVER_SSA_NRM_SPACE_H

/**
 * SolverSSA, the EXPLICIT STATE SPACE variant of the Next Reaction Method: a
 * port of `solver_ssa_nrm_space.m` and of the `else` branch of
 * `solver_ssa_analyzer_nrm.m` that consumes it, together with
 * `solver_ssa_findenabled.m`.
 *
 * HOW IT DIFFERS FROM THE PLAIN NRM. `solver_ssa_nrm.m` integrates the metrics
 * along the path and never stores a state, so its cost is independent of how
 * large the state space is. This variant instead TABULATES the path: it records
 * every distinct state it visits, the time spent in each, and the whole
 * propensity vector at each, and the analyzer then forms the means as `pi * A`
 * rather than as a running time integral. That is the same answer by a different
 * route -- which is exactly what makes the two testable against each other, and
 * what the reference's `options.config.state_space_gen` switch selects between
 * (`none` and `default` take the plain engine, anything else takes this one).
 *
 * WHAT THE TABLE BUYS. A propensity is a deterministic function of the state, so
 * one observation of it IS its exact value and not a sample mean. The departure
 * rates this variant reports are therefore exact per state, and all the Monte
 * Carlo error left in the answer sits in `pi`. The plain engine cannot make that
 * separation, because it never learns that two instants were the same state.
 *
 * WHAT IT COSTS, AND WHY THIS PORT REFUSES OPEN MODELS. The table has one row
 * per distinct state, so the variant is only viable when the reachable space is
 * small enough to enumerate. In the AGGREGATE reaction network an open model has
 * no such space: `rtnodes` gives a Sink no outgoing routing, so jobs accumulate
 * in its slot forever, and a Source's slot is a fictitious token that every
 * arrival decrements and nothing replenishes. Every firing therefore reaches a
 * state never seen before, the table grows one row per sample, and `pi` becomes
 * a uniform law over the trace. The reference does not guard this; here an open
 * model is REFUSED BY NAME, as is any closed model whose space exceeds the
 * declared cap, because a truncated table would report a normalized `pi` over
 * whichever states happened to fit.
 *
 * NO PHASE EXPANSION. The state is the (node, class) population vector, so a
 * non-exponential service process is not represented: `solver_ssa_nrm_space.m`
 * reads `sn.rates` alone and its propensity switch has five arms (EXT, INF, PS,
 * FCFS/LCFS, and non-station). A phase-type process is refused by name rather
 * than silently collapsed onto its mean rate, which would report the right
 * throughput and the wrong queue length.
 *
 * ONE DELIBERATE DIVERGENCE, and it is a reference defect that MATLAB has
 * already fixed on the other side. `solver_ssa_nrm_space.m` line 295 draws the
 * destination with `1+find(rand>=cdfVec{kfire},1)`. `find(...,1)` returns the
 * FIRST true index and the predicate holds on a prefix, so the expression
 * collapses to 2 for every draw above `cdfVec(1)`: with three or more
 * destinations the third onwards are never selected. `solver_ssa_nrm.m` line
 * 1464 names this misuse in a comment and uses the inverse CDF instead, and
 * `NrmEngine` is ported from the fixed form. This file uses the fixed form too,
 * because reproducing the defect would silently misroute jobs and would make the
 * two engines disagree on exactly the models they exist to cross-check.
 *
 * DOUBLE ONLY, for the reason `solver_ssa_nrm.h` gives: the clocks are `-log(u)`
 * of a uniform, there is no exact value to compute, and the answer's error is
 * the Monte Carlo error rather than the rounding.
 *
 * WHY `ssa_find_enabled` IS HERE. It is the standalone form of the scan the
 * serial engine inlines (`solver_ssa.m` carries both and switches on
 * `use_inline`), and it answers the same question this variant is built around:
 * what can fire from a given state, tabulated rather than evaluated on the fly.
 * It works at the STATE ENCODING level, not on the aggregate populations, so it
 * does not interoperate with the engine below and is not used by it; the two are
 * independent ports that share a file.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/lang/qn/state_events.h"
#include "line/solvers/ssa/ssa_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ssa {

// ---------------------------------------------------------------------------
// solver_ssa_findenabled.m
// ---------------------------------------------------------------------------

/** One transition the enabled scan found: where it goes, and at what rate. */
template <class T>
struct EnabledEvent {
    /** `enabled_sync`: index into `sync`, or `sync.size() + g` for a global one. */
    std::size_t sync = 0;
    /** `enabled_rates`: rate * p_active * p_route * p_passive. */
    double rate = 0.0;
    /** `enabled_next_states{act}`: the whole network state after it fires. */
    qn::NetState<T> next;
};

/**
 * Port of `solver_ssa_findenabled.m`: every synchronization that can fire in
 * `state`, with the arrival and departure rates each carries.
 *
 * THE ENUMERATION IS WIDER THAN THE REFERENCE'S, and deliberately so. MATLAB
 * calls `State.afterEvent` with `isSimulation = true`, which SAMPLES one
 * successor row and returns the probability it was drawn with; the C++
 * `after_event` is the enumeration-mode handler and returns every successor with
 * its probability. So one reference entry becomes one entry per (active row,
 * passive row) pair here, and the caller draws from the flattened list in one
 * step instead of two. The induced jump chain is identical; the number of
 * uniforms consumed is not, which is a statement about which random stream this
 * is and not about which model it simulates.
 *
 * A ZERO RATE IS DROPPED rather than rewritten to 1e-38 "so that it is never
 * selected". Keeping it leaves it in the arrival and departure statistics, where
 * it is not a rounding difference but a rate the CTMC generator does not have,
 * and these rates are what the analyzers integrate.
 *
 * `arv` and `dep` are indexed [stateful-1][class-1] and are OVERWRITTEN, not
 * accumulated, so a caller may reuse one pair of buffers across states.
 */
template <class T>
std::vector<EnabledEvent<T>> ssa_find_enabled(const qn::NetworkStruct<T>& sn,
                                              const std::vector<qn::Sync<T>>& sync,
                                              const std::vector<qn::GlobalSync<T>>& gsync,
                                              const qn::NetState<T>& state,
                                              std::vector<std::vector<double>>* arv = nullptr,
                                              std::vector<std::vector<double>>* dep = nullptr) {
    const std::size_t local = sn.nodes.size() + 1;  // the dummy passive node
    const std::size_t R = sn.nclasses;
    const std::size_t NF = sn.stateful_nodes.size();
    std::vector<EnabledEvent<T>> out;
    if (arv) arv->assign(NF, std::vector<double>(R, 0.0));
    if (dep) dep->assign(NF, std::vector<double>(R, 0.0));

    for (std::size_t a = 0; a < sync.size(); ++a) {
        const qn::Sync<T>& sy = sync[a];
        const std::size_t isf_a = sn.stateful_index(sy.active.node);
        if (isf_a == 0) continue;  // a stateless node schedules nothing
        const std::size_t isf_p =
            sy.passive.node == local ? 0 : sn.stateful_index(sy.passive.node);
        if (sy.passive.node != local && isf_p == 0) continue;

        const qn::EventOutcome<T> oa = qn::after_event(sn, sy.active.node, state.local[isf_a - 1],
                                                       sy.active.event, sy.active.cls);
        double fired = 0.0;
        for (std::size_t ia = 0; ia < oa.space.size(); ++ia) {
            const double rate = num_traits<T>::to_double(oa.rate[ia]);
            const double pa = num_traits<T>::to_double(oa.prob[ia]);
            if (!(rate > 0) || !(pa > 0)) continue;

            if (sy.passive.node == local) {
                EnabledEvent<T> e;
                e.sync = a;
                e.rate = rate * pa;
                e.next = state;
                e.next.local[isf_a - 1] = oa.space[ia];
                fired += e.rate;
                out.push_back(e);
                continue;
            }
            // A self-loop synchronization reads the passive node AFTER the
            // active half has been applied, since they are the same node.
            const std::vector<T>& src =
                sy.passive.node == sy.active.node ? oa.space[ia] : state.local[isf_p - 1];
            const qn::EventOutcome<T> op =
                qn::after_event(sn, sy.passive.node, src, sy.passive.event, sy.passive.cls);
            // NO ROWS is the reference's `prob_sync_p = 0`: the destination
            // cannot take the job, so the upstream departure is disabled rather
            // than fired into a state that does not exist.
            for (std::size_t ip = 0; ip < op.space.size(); ++ip) {
                const double pp = num_traits<T>::to_double(op.prob[ip]);
                if (!(pp > 0)) continue;
                EnabledEvent<T> e;
                e.sync = a;
                e.rate = rate * pa * num_traits<T>::to_double(sy.passive.prob) * pp;
                if (!(e.rate > 0)) continue;
                e.next = state;
                e.next.local[isf_a - 1] = oa.space[ia];
                e.next.local[isf_p - 1] = op.space[ip];
                fired += e.rate;
                out.push_back(e);
            }
        }
        // A DEP is one job leaving the active node and entering the passive one,
        // so the same rate is a departure there and an arrival here.
        if (sy.active.event == lang::EventType::DEP && fired > 0) {
            if (dep) (*dep)[isf_a - 1][sy.active.cls - 1] += fired;
            if (arv && isf_p != 0) (*arv)[isf_p - 1][sy.passive.cls - 1] += fired;
        }
    }

    for (std::size_t g = 0; g < gsync.size(); ++g) {
        const qn::GlobalOutcome<T> go = qn::after_global_event(sn, state, gsync[g]);
        for (std::size_t io = 0; io < go.space.size(); ++io) {
            const double w = num_traits<T>::to_double(go.rate[io]) *
                             num_traits<T>::to_double(go.prob[io]);
            if (!(w > 0)) continue;
            EnabledEvent<T> e;
            e.sync = sync.size() + g;
            e.rate = w;
            e.next = go.space[io];
            out.push_back(e);
            if (gsync[g].active.event != lang::EventType::FIRE) continue;
            for (std::size_t j = 0; j < gsync[g].passive.size(); ++j) {
                const qn::ModeEvent<T>& pev = gsync[g].passive[j];
                const std::size_t pisf = sn.stateful_index(pev.node);
                if (pisf == 0 || pev.cls == 0 || pev.cls > R) continue;
                if (pev.event == lang::EventType::PRE && dep)
                    (*dep)[pisf - 1][pev.cls - 1] += w;
                else if (pev.event == lang::EventType::POST && arv)
                    (*arv)[pisf - 1][pev.cls - 1] += w;
            }
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// solver_ssa_nrm_space.m
// ---------------------------------------------------------------------------

/**
 * The knobs the space variant reads.
 *
 * `state_max` is the size of the table this variant is willing to build. It is
 * not a truncation: the run refuses by name on reaching it, because a `pi`
 * normalized over the states that happened to fit is a distribution over the
 * wrong chain.
 */
struct SsaNrmSpaceOptions : SsaOptions {
    std::size_t state_max = 20000;
};

/**
 * One state of the aggregate chain: the (node, class) populations, and the
 * ordered buffer contents of every buffered node.
 *
 * THE BUFFERS ARE PART OF THE STATE, not bookkeeping beside it. Two FCFS states
 * with the same populations but a different waiting order have different
 * propensities, since the rate reads the jobs actually in service; the reference
 * says as much where its hash function concatenates the buffer contents onto the
 * population vector. The buffer is newest-first, so a departure takes the last
 * entry under FCFS and the first under LCFS.
 */
struct NrmSpaceState {
    std::vector<double> n;
    std::vector<std::vector<std::size_t>> buf;
};

/** What `solver_ssa_nrm_space.m` returns, plus what makes it a measurement. */
template <class T>
struct SsaNrmSpaceRun {
    /** `outspace`: the distinct states visited, in first-visit order. */
    std::vector<NrmSpaceState> space;
    /** `pi`: the fraction of simulated time spent in each of them. */
    std::vector<double> pi;
    /**
     * `depRates`, (states x nnodes*nclasses): the departure rate of each (node,
     * class) in each state, read off the cached propensity vector.
     *
     * EXACT PER STATE, not a sample mean: a propensity is a function of the
     * state, so observing it once is knowing it. The reference stores the same
     * vector in `reactCache` for the same reason.
     */
    Matrix<double> dep_rates;
    /** `t`: the cumulative time at each firing. */
    std::vector<double> tran_time;
    /** `kfires`: which reaction fired. */
    std::vector<std::size_t> tran_rx;
    double simulated_time = 0.0;
    std::size_t samples = 0;
    unsigned long seed = 0;
};

/** The metric table, the table it came from, and the stream that produced it. */
template <class T>
struct SsaNrmSpaceSolution {
    SsaSolution avg;
    SsaNrmSpaceRun<T> run;
    unsigned long seed = 0;
};

namespace space_detail {

using lang::NodeType;
using lang::SchedStrategy;

/** True for the two ordered-buffer disciplines this variant's rate law covers. */
inline bool space_sched_buffered(SchedStrategy s) {
    return s == SchedStrategy::FCFS || s == SchedStrategy::LCFS;
}

}  // namespace space_detail

/**
 * The aggregate NRM engine of `solver_ssa_nrm_space.m`.
 *
 * It is a class for the reason `NrmEngine` is, and it is a SEPARATE class rather
 * than a mode of that one because the two disagree about what a state is: this
 * one has no phase dimension, five rate laws instead of thirteen, and a table
 * keyed on (population, buffers) that the plain engine has no place to put. The
 * pieces that look alike -- the stoichiometry from `rtnodes`, the dependency
 * sets, the newest-first buffer -- are alike because both are ports of the same
 * reference construction, and are written out here rather than shared because
 * `NrmEngine` exposes none of them.
 */
template <class T>
class NrmSpaceEngine {
public:
    NrmSpaceEngine(const qn::NetworkStruct<T>& sn, const SsaNrmSpaceOptions& opt)
        : sn_(sn), opt_(opt), rng_(opt.seed) {
        check();
        build_layout();
        build_reactions();
        build_dependencies();
        build_initial_state();
    }

    /** Run `opt.samples` firings and return the tabulated path. */
    SsaNrmSpaceRun<T> run();

    /**
     * The reachable aggregate state space, closed forward from the initial state
     * over the REACTIONS alone.
     *
     * This is the space the run indexes into, and it is derived from the
     * reaction network rather than from the state handlers, so comparing it with
     * `reachable_space_generator`'s answer compares two independent accounts of
     * what the model can do.
     */
    std::vector<NrmSpaceState> enumerate_space() const;

    std::size_t nreactions() const { return nrx_; }
    std::size_t nslots() const { return NS_; }
    const NrmSpaceState& initial() const { return init_; }
    /** The whole propensity vector at a state, the reference's `reactCache` entry. */
    std::vector<double> propensities(const NrmSpaceState& s) const {
        std::vector<double> A(nrx_, 0.0);
        for (std::size_t k = 0; k < nrx_; ++k) A[k] = propensity(k, s);
        return A;
    }

private:
    using SchedStrategy = lang::SchedStrategy;

    const qn::NetworkStruct<T>& sn_;
    SsaNrmSpaceOptions opt_;
    SsaRng rng_;

    std::size_t I_ = 0, K_ = 0, M_ = 0, NS_ = 0, nrx_ = 0;
    std::vector<bool> is_station_;
    std::vector<std::size_t> to_station_;
    std::vector<SchedStrategy> sched_;
    std::vector<double> mi_;
    std::vector<std::vector<double>> rate_;

    Matrix<double> S_;
    /** Per reaction: the slot it consumes, its destination slots and their CDF. */
    std::vector<std::size_t> rx_node_, rx_class_, rx_from_, rx_det_dest_;
    std::vector<std::vector<std::size_t>> rx_to_;
    std::vector<std::vector<double>> rx_cdf_;
    std::vector<std::size_t> rx_nnzp_;
    std::vector<std::vector<std::size_t>> D_;

    NrmSpaceState init_;

    void check();
    void build_layout();
    void build_reactions();
    void build_dependencies();
    void build_initial_state();

    double class_pop(const std::vector<double>& X, std::size_t ind, std::size_t r) const {
        return X[ind * K_ + r];
    }
    double node_pop(const std::vector<double>& X, std::size_t ind) const {
        double s = 0.0;
        for (std::size_t r = 0; r < K_; ++r) s += X[ind * K_ + r];
        return s;
    }
    double propensity(std::size_t j, const NrmSpaceState& s) const;
    void update_buffers(std::size_t kfire, const std::vector<double>& n,
                        std::vector<std::vector<std::size_t>>& bufs, std::size_t dest_pos,
                        bool have_dest) const;

    static constexpr std::size_t npos = static_cast<std::size_t>(-1);
};

/**
 * What this variant cannot represent, refused before a single firing.
 *
 * The scheduling list is the propensity switch of `solver_ssa_nrm_space.m`, and
 * it is much shorter than the plain NRM's: the reference's own analyzer stub
 * narrows it further still, to INF/EXT/PS.
 */
template <class T>
void NrmSpaceEngine<T>::check() {
    for (std::size_t i = 0; i < sn_.nstations; ++i) {
        const SchedStrategy s = sn_.stations[i].sched;
        if (s == SchedStrategy::INF || s == SchedStrategy::EXT || s == SchedStrategy::PS ||
            s == SchedStrategy::FCFS || s == SchedStrategy::LCFS)
            continue;
        throw UnsupportedError(
            "solver_ssa_nrm_space: station '" + sn_.stations[i].name + "' uses '" +
            std::string(lang::sched_to_text(s)) +
            "' scheduling, which has no rate law in the explicit state-space variant. Its "
            "propensity switch covers EXT, INF, PS, FCFS and LCFS only; use method='nrm', whose "
            "engine carries the full set");
    }
    for (std::size_t i = 0; i < sn_.nstations; ++i)
        for (std::size_t r = 0; r < sn_.nclasses; ++r) {
            if (sn_.disabled[i][r]) continue;
            const lang::ProcessType pt = sn_.service[i][r].type;
            if (pt == lang::ProcessType::EXP || pt == lang::ProcessType::DISABLED ||
                pt == lang::ProcessType::IMMEDIATE)
                continue;
            throw UnsupportedError(
                "solver_ssa_nrm_space: class '" + sn_.classes[r].name +
                "' has non-exponential service at station '" + sn_.stations[i].name +
                "'. The explicit state-space variant carries no phase dimension -- its state is "
                "the (node, class) population vector and it reads sn.rates alone -- so a "
                "phase-type process cannot be represented; use method='nrm', which expands it");
        }
    for (std::size_t r = 0; r < sn_.nclasses; ++r)
        if (!std::isfinite(sn_.classes[r].population))
            throw UnsupportedError(
                "solver_ssa_nrm_space: class '" + sn_.classes[r].name +
                "' is open. In the aggregate reaction network a Sink has no outgoing routing and "
                "a Source's slot is a fictitious token that every arrival consumes, so the "
                "aggregate state space of an open model is unbounded: every firing would add a "
                "row to the table and pi would become the uniform law over the trace. Use "
                "method='nrm', which integrates the metrics and stores no states");
    for (const qn::NodeDef& nd : sn_.nodes) {
        switch (nd.nodetype) {
            case lang::NodeType::Cache:
                throw UnsupportedError(
                    "solver_ssa_nrm_space: node '" + nd.name +
                    "' is a Cache. The aggregate reaction network is built from sn.rtnodes and "
                    "sn.rates alone and carries no cache contents, so the hit/miss class switch "
                    "cannot be resolved at firing time");
            case lang::NodeType::Place:
            case lang::NodeType::Transition:
                throw UnsupportedError(
                    "solver_ssa_nrm_space: node '" + nd.name +
                    "' makes this model a stochastic Petri net. A firing is atomic across every "
                    "arc it touches and is not a (node, class) departure, which is the only "
                    "reaction shape this variant builds");
            case lang::NodeType::Fork:
            case lang::NodeType::Join:
                throw UnsupportedError(
                    "solver_ssa_nrm_space: node '" + nd.name +
                    "' makes this a fork-join model. A fork emits on several branches at once, "
                    "which no single-destination reaction can express, and the NRM does not "
                    "handle fork-join in any codebase");
            default:
                break;
        }
    }
    // `sn.isslc`, the SELF-LOOPING CLASS whose reaction consumes a job and
    // produces none, has no field in the C++ NetworkStruct, so the reference's
    // `Srow(from) = -Inf` branch is unreachable here rather than unported.
}

template <class T>
void NrmSpaceEngine<T>::build_layout() {
    I_ = sn_.nof_nodes();
    K_ = sn_.nclasses;
    M_ = sn_.nstations;
    NS_ = I_ * K_;

    is_station_.assign(I_, false);
    to_station_.assign(I_, npos);
    for (std::size_t i = 0; i < I_; ++i)
        if (sn_.nodes[i].station != 0) {
            is_station_[i] = true;
            to_station_[i] = sn_.nodes[i].station - 1;
        }
    sched_.assign(M_, SchedStrategy::FCFS);
    for (std::size_t i = 0; i < M_; ++i) sched_[i] = sn_.stations[i].sched;

    // A non-station is given the immediate rate and an unbounded server count,
    // which is what makes a Router or a ClassSwitch a pass-through in the
    // aggregate chain rather than a place a job can accumulate.
    mi_.assign(I_, lang::GlobalConstants::MaxInt);
    rate_.assign(I_, std::vector<double>(K_, 0.0));
    for (std::size_t i = 0; i < I_; ++i) {
        if (is_station_[i]) {
            const std::size_t ist = to_station_[i];
            mi_[i] = sn_.stations[ist].nservers;
            if (!std::isfinite(mi_[i])) mi_[i] = lang::GlobalConstants::MaxInt;
            for (std::size_t r = 0; r < K_; ++r)
                if (!sn_.disabled[ist][r])
                    rate_[i][r] = num_traits<T>::to_double(sn_.rates(ist, r));
        } else {
            for (std::size_t r = 0; r < K_; ++r) rate_[i][r] = lang::GlobalConstants::Immediate;
        }
    }
}

/**
 * The stoichiometry, one reaction per (node, class).
 *
 * There is no phase dimension, so a reaction IS a departure: it consumes one job
 * from its own slot and deposits the routing probability into every slot
 * `rtnodes` sends it to. The fractional entries are not a fluid relaxation --
 * they are the weights the destination draw reads, and exactly one destination
 * receives the whole job at firing time.
 */
template <class T>
void NrmSpaceEngine<T>::build_reactions() {
    nrx_ = I_ * K_;
    rx_node_.assign(nrx_, 0);
    rx_class_.assign(nrx_, 0);
    rx_from_.assign(nrx_, 0);
    rx_det_dest_.assign(nrx_, npos);
    rx_to_.assign(nrx_, std::vector<std::size_t>());
    rx_cdf_.assign(nrx_, std::vector<double>());
    rx_nnzp_.assign(nrx_, 0);
    S_ = Matrix<double>(NS_, nrx_, 0.0);

    for (std::size_t ind = 0; ind < I_; ++ind)
        for (std::size_t r = 0; r < K_; ++r) {
            const std::size_t k = ind * K_ + r;
            rx_node_[k] = ind;
            rx_class_[k] = r;
            rx_from_[k] = k;
            S_(k, k) -= 1.0;
            for (std::size_t jnd = 0; jnd < I_; ++jnd)
                for (std::size_t s = 0; s < K_; ++s) {
                    const double p =
                        num_traits<T>::to_double(sn_.rtnodes(ind * K_ + r, jnd * K_ + s));
                    if (!(p > 0.0)) continue;
                    S_(jnd * K_ + s, k) += p;
                }
        }

    // `P = S; P(P<0) = P(P<0)+1`: a class routing back into its own slot must
    // compete with the other destinations on equal footing, which subtracting
    // the consumed job first is what achieves.
    for (std::size_t k = 0; k < nrx_; ++k) {
        std::vector<double> Pcol(NS_, 0.0);
        for (std::size_t i = 0; i < NS_; ++i) {
            const double v = S_(i, k);
            Pcol[i] = v < 0.0 ? v + 1.0 : v;
            if (v > 0.0 && rx_det_dest_[k] == npos) rx_det_dest_[k] = i;
        }
        for (std::size_t i = 0; i < NS_; ++i)
            if (Pcol[i] != 0.0) ++rx_nnzp_[k];
        // A job that routes back into the slot it came from with probability one
        // leaves an ALL-ZERO stoichiometry column, so `find(S(:,k) > 0)` finds
        // nothing and the reference records no destination. At a buffered
        // station that loses the job from the buffer while leaving the
        // population alone, breaking numel(buf) == max(0, total - mi) and
        // silencing the rate law. Naming the source as the destination is the
        // same arithmetic on the population and the correct one on the buffer.
        if (rx_nnzp_[k] == 1 && rx_det_dest_[k] == npos && Pcol[rx_from_[k]] != 0.0)
            rx_det_dest_[k] = rx_from_[k];
        if (rx_nnzp_[k] > 1) {
            double acc = 0.0;
            for (std::size_t i = 0; i < NS_; ++i)
                if (Pcol[i] != 0.0) {
                    rx_to_[k].push_back(i);
                    acc += Pcol[i];
                    rx_cdf_[k].push_back(acc);
                }
        }
    }
}

/**
 * The dependency sets: which propensities a firing invalidates.
 *
 * Every rate law reads its node's WHOLE class-count vector, so a firing that
 * touches any slot of a node invalidates every reaction consuming from that
 * node, not merely the slot that moved.
 */
template <class T>
void NrmSpaceEngine<T>::build_dependencies() {
    D_.assign(nrx_, std::vector<std::size_t>());
    for (std::size_t k = 0; k < nrx_; ++k) {
        std::vector<bool> touched(I_, false);
        for (std::size_t i = 0; i < NS_; ++i)
            if (S_(i, k) != 0.0) touched[i / K_] = true;
        for (std::size_t ind = 0; ind < I_; ++ind) {
            if (!touched[ind]) continue;
            for (std::size_t r = 0; r < K_; ++r) D_[k].push_back(ind * K_ + r);
        }
        std::sort(D_[k].begin(), D_[k].end());
        D_[k].erase(std::unique(D_[k].begin(), D_[k].end()), D_[k].end());
    }
}

/**
 * The initial state.
 *
 * The reference reads `sn.state` through `State.toMarginalAggr`; this port has
 * no State package on the aggregate side, so it uses the same rule
 * `solver_ssa_nrm.h` does -- a closed class starts entirely at its reference
 * station. The chain is ergodic, so the steady-state means do not depend on the
 * choice.
 */
template <class T>
void NrmSpaceEngine<T>::build_initial_state() {
    init_.n.assign(NS_, 0.0);
    init_.buf.assign(I_, std::vector<std::size_t>());
    for (std::size_t r = 0; r < K_; ++r) {
        const double pop = sn_.classes[r].population;
        if (!(pop > 0.0)) continue;
        const std::size_t rs = sn_.classes[r].refstat;
        if (rs < 1 || rs > M_)
            throw InputError("solver_ssa_nrm_space: class '" + sn_.classes[r].name +
                             "' has no reference station");
        init_.n[(sn_.station_to_node[rs - 1] - 1) * K_ + r] = pop;
    }
    // The buffer holds exactly the waiting jobs, so its length is the excess
    // over the server count; which classes wait is immaterial to the steady
    // state, so they are taken in class order.
    for (std::size_t ind = 0; ind < I_; ++ind) {
        if (!is_station_[ind] || !space_detail::space_sched_buffered(sched_[to_station_[ind]]))
            continue;
        double waiting = std::max(0.0, node_pop(init_.n, ind) - mi_[ind]);
        for (std::size_t r = 0; r < K_ && waiting > 0.0; ++r) {
            const double take = std::min(waiting, init_.n[ind * K_ + r]);
            for (std::size_t c = 0; c < static_cast<std::size_t>(take); ++c)
                init_.buf[ind].push_back(r);
            waiting -= take;
        }
    }
}

/** The five-armed propensity switch of `solver_ssa_nrm_space.m` lines 136-161. */
template <class T>
double NrmSpaceEngine<T>::propensity(std::size_t j, const NrmSpaceState& st) const {
    const std::size_t ind = rx_node_[j], r = rx_class_[j];
    const std::vector<double>& X = st.n;
    if (!is_station_[ind]) {
        // A pass-through node moves at most one job at a time at the immediate
        // rate: the min is what keeps an empty node silent rather than firing at
        // 1e8 into a slot that holds nothing.
        return rate_[ind][r] * std::min(1.0, class_pop(X, ind, r));
    }
    const std::size_t ist = to_station_[ind];
    const double eps = lang::GlobalConstants::Zero;
    switch (sched_[ist]) {
        case SchedStrategy::EXT:
            return rate_[ind][r];
        case SchedStrategy::INF:
            return rate_[ind][r] * class_pop(X, ind, r);
        case SchedStrategy::PS: {
            if (K_ == 1) return rate_[ind][r] * std::min(mi_[ind], class_pop(X, ind, r));
            const double tot = node_pop(X, ind);
            return rate_[ind][r] * (class_pop(X, ind, r) / (eps + tot)) *
                   std::min(mi_[ind], eps + tot);
        }
        case SchedStrategy::FCFS:
        case SchedStrategy::LCFS: {
            // The jobs actually in service are the class population less the
            // class-r jobs still waiting, which is why the buffer has to be part
            // of the state and not merely of the bookkeeping.
            double waiting = 0.0;
            for (std::size_t c : st.buf[ind])
                if (c == r) waiting += 1.0;
            return rate_[ind][r] * std::max(0.0, class_pop(X, ind, r) - waiting);
        }
        default:
            throw UnsupportedError(
                "solver_ssa_nrm_space: the scheduling policy '" +
                std::string(lang::sched_to_text(sched_[ist])) + "' at station '" +
                sn_.stations[ist].name + "' has no rate law in the explicit state-space variant");
    }
}

/**
 * `updateBuffers`: the ordered buffers after reaction `kfire` has fired.
 *
 * A departure frees a server, so the discipline promotes a waiting job: FCFS
 * takes the oldest, which is the LAST entry of a newest-first buffer, and LCFS
 * the newest, which is the first. An arrival joins the buffer only when it finds
 * every server busy, and the test is made on the POST-firing population, which
 * already counts the arriving job.
 */
template <class T>
void NrmSpaceEngine<T>::update_buffers(std::size_t kfire, const std::vector<double>& n,
                                       std::vector<std::vector<std::size_t>>& bufs,
                                       std::size_t dest_pos, bool have_dest) const {
    const std::size_t ind = rx_node_[kfire];
    if (is_station_[ind] && !bufs[ind].empty()) {
        const SchedStrategy s = sched_[to_station_[ind]];
        if (s == SchedStrategy::FCFS) bufs[ind].pop_back();
        else if (s == SchedStrategy::LCFS) bufs[ind].erase(bufs[ind].begin());
    }
    if (!have_dest || dest_pos == npos) return;
    const std::size_t jnd = dest_pos / K_, s = dest_pos % K_;
    if (!is_station_[jnd] || !space_detail::space_sched_buffered(sched_[to_station_[jnd]])) return;
    if (node_pop(n, jnd) > mi_[jnd]) bufs[jnd].insert(bufs[jnd].begin(), s);
}

namespace space_detail {

/** The table key: the populations, then each node's buffer behind a separator. */
inline std::vector<double> space_key(const NrmSpaceState& s) {
    std::vector<double> key = s.n;
    for (std::size_t i = 0; i < s.buf.size(); ++i) {
        // The separator keeps two different splits of the same concatenation
        // from colliding, which a plain flatten would allow.
        key.push_back(-1.0);
        for (std::size_t j = 0; j < s.buf[i].size(); ++j)
            key.push_back(static_cast<double>(s.buf[i][j]) + 1.0);
    }
    return key;
}

}  // namespace space_detail

template <class T>
std::vector<NrmSpaceState> NrmSpaceEngine<T>::enumerate_space() const {
    std::vector<NrmSpaceState> out;
    std::map<std::vector<double>, std::size_t> seen;
    std::vector<std::size_t> stack;
    seen[space_detail::space_key(init_)] = 0;
    out.push_back(init_);
    stack.push_back(0);

    while (!stack.empty()) {
        const std::size_t si = stack.back();
        stack.pop_back();
        const NrmSpaceState st = out[si];  // by value: `out` grows inside the loop
        for (std::size_t k = 0; k < nrx_; ++k) {
            if (!(propensity(k, st) > 0.0)) continue;
            // Every destination the firing could draw is a distinct successor,
            // which is what makes this the reachable set of the CHAIN rather
            // than of one realization of it.
            std::vector<std::size_t> dests;
            if (rx_nnzp_[k] > 1) dests = rx_to_[k];
            else if (rx_det_dest_[k] != npos) dests.push_back(rx_det_dest_[k]);
            else dests.push_back(npos);
            for (std::size_t d = 0; d < dests.size(); ++d) {
                NrmSpaceState ns = st;
                ns.n[rx_from_[k]] -= 1.0;
                if (dests[d] != npos) ns.n[dests[d]] += 1.0;
                update_buffers(k, ns.n, ns.buf, dests[d], dests[d] != npos);
                const std::vector<double> key = space_detail::space_key(ns);
                if (seen.find(key) != seen.end()) continue;
                if (out.size() >= opt_.state_max)
                    throw UnsupportedError(
                        "solver_ssa_nrm_space: the reachable aggregate state space exceeds the "
                        "cap of " + std::to_string(opt_.state_max) +
                        " states. The explicit state-space variant tabulates one row per state, "
                        "so a larger space is refused rather than truncated: a pi normalized over "
                        "the states that happened to fit is the law of a different chain. Raise "
                        "state_max, or use method='nrm', which stores no states");
                seen[key] = out.size();
                out.push_back(ns);
                stack.push_back(out.size() - 1);
            }
        }
    }
    return out;
}

template <class T>
SsaNrmSpaceRun<T> NrmSpaceEngine<T>::run() {
    // The clocks are -log(u), so the backend must have a logarithm at all. The
    // analyzer gates on `double` before instantiating this; a caller reaching
    // the assert got past that gate and deserves a sentence rather than a
    // compile error inside the uniform draw.
    static_assert(num_traits<T>::has_transcendental,
                  "solver_ssa_nrm_space: the Next Reaction Method draws its clocks as -log(u), "
                  "which needs transcendental arithmetic");

    SsaNrmSpaceRun<T> out;
    out.seed = opt_.seed;
    if (nrx_ == 0) return out;

    NrmSpaceState cur = init_;
    std::vector<double> Ak(nrx_, 0.0), Pk(nrx_, 0.0), Tk(nrx_, 0.0), tau(nrx_, 0.0);
    for (std::size_t k = 0; k < nrx_; ++k) {
        Ak[k] = propensity(k, cur);
        Pk[k] = -std::log(rng_.uniform());
        tau[k] = Ak[k] > 0.0 ? (Pk[k] - Tk[k]) / Ak[k] : std::numeric_limits<double>::infinity();
    }

    // The table. `pi` accumulates holding time and the propensity vector is
    // stored the first time a state is seen, exactly as `reactCache` does.
    std::map<std::vector<double>, std::size_t> index;
    std::vector<std::vector<double>> cached;
    double total_time = 0.0;

    out.tran_time.reserve(opt_.samples);
    out.tran_rx.reserve(opt_.samples);
    for (std::size_t n = 0; n < opt_.samples; ++n) {
        std::size_t kfire = 0;
        double dt = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < nrx_; ++k)
            if (tau[k] < dt) {
                dt = tau[k];
                kfire = k;
            }
        if (std::isinf(dt))
            throw NumericError(
                "solver_ssa_nrm_space: deadlock -- every reaction has propensity zero, so the "
                "sample path cannot advance");

        // The state is recorded with the time spent IN it, so the pair belongs
        // to the state before the firing.
        const std::vector<double> key = space_detail::space_key(cur);
        const std::map<std::vector<double>, std::size_t>::const_iterator it = index.find(key);
        std::size_t si;
        if (it != index.end()) {
            si = it->second;
        } else {
            if (out.space.size() >= opt_.state_max)
                throw UnsupportedError(
                    "solver_ssa_nrm_space: the path has visited more than the cap of " +
                    std::to_string(opt_.state_max) +
                    " distinct states. The explicit state-space variant tabulates one row per "
                    "state, so it is refused rather than truncated: a pi normalized over the "
                    "states that happened to fit is the law of a different chain. Raise "
                    "state_max, or use method='nrm', which stores no states");
            si = out.space.size();
            index[key] = si;
            out.space.push_back(cur);
            out.pi.push_back(0.0);
            cached.push_back(Ak);
        }
        out.pi[si] += dt;
        total_time += dt;
        out.tran_time.push_back(total_time);
        out.tran_rx.push_back(kfire);

        // Apply the firing.
        std::size_t dest = npos;
        if (rx_nnzp_[kfire] > 1) {
            // Inverse CDF: the smallest index whose cumulative weight exceeds
            // the draw. See the divergence note at the top of this file.
            const double u = rng_.uniform();
            std::size_t sel = rx_cdf_[kfire].size() - 1;
            for (std::size_t i = 0; i < rx_cdf_[kfire].size(); ++i)
                if (rx_cdf_[kfire][i] > u) {
                    sel = i;
                    break;
                }
            dest = rx_to_[kfire][sel];
            cur.n[rx_from_[kfire]] -= 1.0;
            cur.n[dest] += 1.0;
        } else {
            for (std::size_t i = 0; i < NS_; ++i)
                if (S_(i, kfire) != 0.0) cur.n[i] += S_(i, kfire);
            dest = rx_det_dest_[kfire];
        }
        update_buffers(kfire, cur.n, cur.buf, dest, dest != npos);

        for (std::size_t k = 0; k < nrx_; ++k) Tk[k] += Ak[k] * dt;
        // D covers every node the stoichiometry touched, source and every
        // candidate destination alike, and the arrival's buffer join is already
        // applied when those are refreshed. What it does NOT cover is a
        // PROMOTION: the departing job's replacement enters service without any
        // slot moving, and at a self-routing reaction the column is all zeros
        // and D is empty outright. A firing at a buffered node therefore forces
        // a sweep, which is the same order as the clock scan above and so costs
        // nothing asymptotically.
        for (std::size_t k : D_[kfire]) Ak[k] = propensity(k, cur);
        if (is_station_[rx_node_[kfire]] &&
            space_detail::space_sched_buffered(sched_[to_station_[rx_node_[kfire]]]))
            for (std::size_t k = 0; k < nrx_; ++k) Ak[k] = propensity(k, cur);

        Pk[kfire] -= std::log(rng_.uniform());
        for (std::size_t k = 0; k < nrx_; ++k)
            tau[k] =
                Ak[k] > 0.0 ? (Pk[k] - Tk[k]) / Ak[k] : std::numeric_limits<double>::infinity();
        out.samples = n + 1;
    }

    double tot = 0.0;
    for (std::size_t s = 0; s < out.pi.size(); ++s) tot += out.pi[s];
    if (tot > 0)
        for (std::size_t s = 0; s < out.pi.size(); ++s) out.pi[s] /= tot;
    out.simulated_time = total_time;

    // The departure rate of (node, class) in a state IS the propensity of its
    // reaction there, which is why the cache is what the analyzer reads.
    out.dep_rates = Matrix<double>(out.space.size(), NS_, 0.0);
    for (std::size_t s = 0; s < out.space.size(); ++s)
        for (std::size_t k = 0; k < nrx_; ++k)
            out.dep_rates(s, rx_from_[k]) += cached[s][k];
    return out;
}

/** `solver_ssa_nrm_space.m`: run the tabulating engine. */
template <class T>
SsaNrmSpaceRun<T> solver_ssa_nrm_space(const qn::NetworkStruct<T>& sn,
                                       const SsaNrmSpaceOptions& opt) {
    NrmSpaceEngine<T> eng(sn, opt);
    return eng.run();
}

/**
 * Port of the `else` branch of `solver_ssa_analyzer_nrm.m`, the one
 * `state_space_gen` selects: the means as `pi * A` over the tabulated states.
 *
 * THE STANDALONE `solver_ssa_nrm_space_analyzer.m` IS A STUB. Its whole body is
 * the INF/EXT/PS scheduling gate and a debug line; it assigns none of its nine
 * declared outputs, so calling it in MATLAB raises "Output argument not
 * assigned". The working analyzer is the branch ported here, whose utilization
 * switch admits FCFS and LCFS as well -- which is also the set the engine has
 * rate laws for, so narrowing to the stub's three would refuse models this
 * variant can simulate. The stub's narrower gate is therefore NOT reproduced,
 * and the divergence is named rather than hidden.
 *
 * A SOURCE REPORTS ZERO QLen AND ZERO Util, not the negative number the
 * reference's `UN(ist,:) = QN(ist,:)` produces for an EXT station: the Source
 * slot is a fictitious token that arrivals consume, so its time average is
 * 1 - E[jobs in system]. `solver_ssa_nrm.h` states the rule at length and this
 * port applies it everywhere. It is moot for the models this variant accepts,
 * since an open one is refused outright, and it is kept so the rule holds
 * uniformly.
 */
template <class T>
SsaNrmSpaceSolution<T> solver_ssa_nrm_space_analyzer(const qn::NetworkStruct<T>& sn,
                                                     const SsaNrmSpaceOptions& opt) {
    // `if constexpr`, not a run-time test: the engine takes the logarithm of a
    // uniform, so a Rational instantiation would fail to COMPILE rather than
    // refuse. The gate has to keep the body from being instantiated at all.
    if constexpr (!std::is_same<T, double>::value) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_ssa_nrm_space: an SSA sample path is generated from exponential clocks, which "
            "are logarithms of uniform draws; there is no exact value to compute and a wider "
            "float carries no information the Monte Carlo error does not swamp. Rerun with "
            "--arith double");
    } else {
        using lang::SchedStrategy;
        const std::size_t M = sn.nstations, K = sn.nclasses;
        SsaNrmSpaceSolution<T> out;
        out.seed = opt.seed;

        NrmSpaceEngine<T> eng(sn, opt);
        out.run = eng.run();
        const SsaNrmSpaceRun<T>& r = out.run;

        SsaSolution& a = out.avg;
        a.method = "nrm.space";
        a.samples = r.samples;
        a.simulated_time = r.simulated_time;
        a.QN = Matrix<double>(M, K, 0.0);
        a.UN = Matrix<double>(M, K, 0.0);
        a.RN = Matrix<double>(M, K, 0.0);
        a.TN = Matrix<double>(M, K, 0.0);
        a.XN.assign(K, 0.0);
        a.CN.assign(K, 0.0);

        for (std::size_t k = 0; k < K; ++k) {
            const std::size_t refnd = sn.station_to_node[sn.classes[k].refstat - 1];
            for (std::size_t s = 0; s < r.space.size(); ++s)
                a.XN[k] += r.pi[s] * r.dep_rates(s, (refnd - 1) * K + k);
        }

        for (std::size_t ist = 0; ist < M; ++ist) {
            const std::size_t ind = sn.station_to_node[ist];
            for (std::size_t k = 0; k < K; ++k)
                for (std::size_t s = 0; s < r.space.size(); ++s) {
                    a.TN(ist, k) += r.pi[s] * r.dep_rates(s, (ind - 1) * K + k);
                    a.QN(ist, k) += r.pi[s] * r.space[s].n[(ind - 1) * K + k];
                }

            const SchedStrategy sched = sn.stations[ist].sched;
            if (sched == SchedStrategy::EXT ||
                sn.stations[ist].nodetype == lang::NodeType::Source) {
                for (std::size_t k = 0; k < K; ++k) a.QN(ist, k) = 0.0;
                continue;
            }
            if (sched == SchedStrategy::INF) {
                for (std::size_t k = 0; k < K; ++k) a.UN(ist, k) = a.QN(ist, k);
                continue;
            }
            // PS, FCFS and LCFS: the carried load T/(mu*c). A class-dependent
            // station normalizes by its DECLARED peak instead, which is the only
            // thing its utilization can be a fraction of.
            const bool is_cd = static_cast<bool>(sn.stations[ist].cdscaling);
            const bool is_jd = static_cast<bool>(sn.stations[ist].jdscaling);
            for (std::size_t k = 0; k < K; ++k) {
                if (sn.disabled[ist][k]) continue;
                const double mu = num_traits<T>::to_double(sn.rates(ist, k));
                if (!(mu > 0)) continue;
                // The divisor is the PRODUCT of whichever declared peaks are
                // present -- `solver_ssa_nrm.m:1784-1785` -- and the server count
                // when neither is.
                double sdiv = sn.stations[ist].nservers;
                if (is_cd || is_jd) {
                    sdiv = 1.0;
                    const std::vector<T>* pks[2] = {&sn.stations[ist].cdscalingpeak,
                                                    &sn.stations[ist].jdscalingpeak};
                    const char* names[2] = {"setClassDependence", "setJointDependence"};
                    const bool on[2] = {is_cd, is_jd};
                    for (std::size_t h = 0; h < 2; ++h) {
                        if (!on[h]) continue;
                        const std::vector<T>& pk = *pks[h];
                        if (pk.size() <= k || !(num_traits<T>::to_double(pk[k]) > 0))
                            throw InputError(
                                "solver_ssa_nrm_space: station '" + sn.stations[ist].name +
                                "' declares a dependent scaling with no declared peak rate. "
                                "Utilization there is T/mu/peak, so pass the peak to " + names[h]);
                        sdiv *= num_traits<T>::to_double(pk[k]);
                    }
                }
                a.UN(ist, k) = sdiv > 0 ? a.TN(ist, k) / mu / sdiv : 0.0;
            }
        }

        for (std::size_t k = 0; k < K; ++k) {
            for (std::size_t ist = 0; ist < M; ++ist)
                a.RN(ist, k) = a.TN(ist, k) > 0 ? a.QN(ist, k) / a.TN(ist, k) : 0.0;
            if (a.XN[k] > 0) a.CN[k] = sn.classes[k].population / a.XN[k];
        }
        return out;
    }
}

}  // namespace ssa
}  // namespace line

#endif  // LINE_SOLVERS_SSA_SOLVER_SSA_NRM_SPACE_H
