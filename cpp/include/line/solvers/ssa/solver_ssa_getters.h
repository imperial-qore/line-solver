/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SSA_SOLVER_SSA_GETTERS_H
#define LINE_SOLVERS_SSA_SOLVER_SSA_GETTERS_H

/**
 * The SolverSSA queries that are not the average table: `getProb`,
 * `getProbAggr`, `getProbSys`, `getProbSysAggr`, and the four samplers
 * `sample`, `sampleAggr`, `sampleSys`, `sampleSysAggr`.
 *
 * THEY ALL READ ONE SAMPLE PATH, which is why they live together. The reference
 * splits them across nine files in `@@SolverSSA/`, but every one of them runs
 * the analyzer and then reduces the SAME trajectory: a probability is the
 * fraction of simulated time spent in a state, and a sample is that trajectory
 * printed. Running the engine once and reducing it several ways is therefore
 * not an optimization -- it is what makes the four probabilities mutually
 * consistent, which they would not be if each drew its own path.
 *
 * THE SERIAL ENGINE, NEVER THE NRM, and that is the reference's own choice:
 * `@@SolverSSA/getProb.m` rewrites `default` and `nrm` to `serial` before it
 * asks for the trajectory. The NRM does not simulate the state ENCODING at all
 * -- it moves per-(node, class, phase) counts -- so it has no row to compare a
 * requested state against and no per-event trace to print. A caller who asked
 * for `nrm` and then for a probability gets the serial engine's answer, as they
 * do in MATLAB.
 *
 * TIME-WEIGHTED, NOT EVENT-COUNTED. A probability here is the fraction of
 * SIMULATED TIME the path spent in the state, not the fraction of firings that
 * landed in it. The two differ by the holding-time distribution and only the
 * first estimates the stationary law; the reference's `TSS(:,1) = [TSS(1,1);
 * diff(TSS(:,1))]` is exactly that reweighting, and `SsaSerialRun::pi` has
 * already applied it, so these functions sum `pi` and never counts.
 *
 * A STATE THAT WAS NEVER VISITED HAS PROBABILITY 0, and the reference warns
 * before returning it. That zero is a measurement of a finite path and not a
 * statement about the chain: a rare state simply did not occur in the sample.
 * `SsaProbResult::seen` carries the distinction so a caller can tell "the path
 * says this state is impossible" from "the path never got there".
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/lang/qn/state_events.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ssa/solver_ssa_serial.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ssa {

/** One probability query: the estimate, and whether the state occurred at all. */
struct SsaProbResult {
    double prob = 0.0;
    bool seen = false;
};

/** The four probabilities `-a prob` reports, over one requested state. */
struct SsaProbReport {
    SsaProbResult sys, sys_aggr;             ///< getProbSys, getProbSysAggr
    std::vector<SsaProbResult> marg, aggr;   ///< getProb, getProbAggr, per station
    std::size_t samples = 0;
    unsigned long seed = 0;
    double simulated_time = 0.0;
};

/** One trajectory, in the shape the reference's `sampleSys` returns it. */
template <class T>
struct SsaSamplePath {
    std::vector<double> t;              ///< the event times, increasing
    std::vector<std::size_t> event;     ///< which synchronization fired
    Matrix<T> state;                    ///< per event: the state OCCUPIED until then
    Matrix<T> aggr;                     ///< the same, as per-(stateful, class) counts
    std::size_t samples = 0;
    unsigned long seed = 0;
};

namespace getters_detail {

/** The per-class counts of one node's local row, `State.toMarginal`'s `nir`. */
template <class T>
std::vector<double> node_counts(const qn::NetworkStruct<T>& sn, std::size_t ind,
                                const std::vector<T>& row) {
    const std::pair<T, std::vector<T>> m = qn::to_marginal_aggr(sn, ind, row);
    std::vector<double> out(m.second.size(), 0.0);
    for (std::size_t k = 0; k < m.second.size(); ++k)
        out[k] = num_traits<T>::to_double(m.second[k]);
    return out;
}

/** Exact equality of two encoded rows, at the tolerance a count deserves. */
template <class T>
bool same_row(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t j = 0; j < a.size(); ++j)
        if (std::fabs(num_traits<T>::to_double(a[j]) - num_traits<T>::to_double(b[j])) > 1e-9)
            return false;
    return true;
}

inline bool same_counts(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t j = 0; j < a.size(); ++j)
        if (std::fabs(a[j] - b[j]) > 1e-9) return false;
    return true;
}

/**
 * The node's name for a message, or its index when there is no such node.
 *
 * A query naming a node past the end is the case that made this necessary:
 * reading `sn.nodes[ind - 1].name` to SAY SO is itself the out-of-range access,
 * so the refusal aborted instead of being thrown.
 */
template <class T>
std::string node_label(const qn::NetworkStruct<T>& sn, std::size_t ind) {
    if (ind == 0 || ind > sn.nodes.size()) return "#" + std::to_string(ind);
    return sn.nodes[ind - 1].name;
}

}  // namespace getters_detail

/**
 * `getProb(node, state)`: the time fraction the path spent with node `ind` in
 * `row`.
 *
 * The row is compared at the ENCODING width the path lives at, which is why the
 * caller is expected to pass a row taken from the same widening the engine
 * applies (`serial_detail::wide_init_state`, or a row of `run.space`). The
 * reference left-pads a short row with zeros; a row that is short here would
 * compare unequal to every state, so it is refused by name instead of silently
 * matching nothing.
 */
template <class T>
SsaProbResult ssa_prob(const qn::NetworkStruct<T>& sn, const SsaSerialRun<T>& r, std::size_t ind,
                       const std::vector<T>& row) {
    const std::size_t isf = sn.stateful_index(ind);
    if (isf == 0)
        throw InputError("ssa_prob: node '" + getters_detail::node_label(sn, ind) +
                         "' is stateless, so it has no state to report a probability of");
    SsaProbResult out;
    for (std::size_t s = 0; s < r.space.size(); ++s) {
        if (!getters_detail::same_row(r.space[s].local[isf - 1], row)) continue;
        out.prob += r.pi[s];
        out.seen = true;
    }
    return out;
}

/** `getProbAggr(node, n)`: the same, over the per-class counts alone. */
template <class T>
SsaProbResult ssa_prob_aggr(const qn::NetworkStruct<T>& sn, const SsaSerialRun<T>& r,
                            std::size_t ind, const std::vector<double>& counts) {
    const std::size_t isf = sn.stateful_index(ind);
    if (isf == 0)
        throw InputError("ssa_prob_aggr: node '" + getters_detail::node_label(sn, ind) +
                         "' is stateless, so it holds no jobs to count");
    SsaProbResult out;
    for (std::size_t s = 0; s < r.space.size(); ++s) {
        if (!getters_detail::same_counts(
                getters_detail::node_counts(sn, ind, r.space[s].local[isf - 1]), counts))
            continue;
        out.prob += r.pi[s];
        out.seen = true;
    }
    return out;
}

/** `getProbSys()`: the joint state of every stateful node at once. */
template <class T>
SsaProbResult ssa_prob_sys(const qn::NetworkStruct<T>& sn, const SsaSerialRun<T>& r,
                           const qn::NetState<T>& st) {
    SsaProbResult out;
    for (std::size_t s = 0; s < r.space.size(); ++s) {
        bool all = true;
        for (std::size_t f = 0; f < sn.stateful_nodes.size() && all; ++f)
            if (!getters_detail::same_row(r.space[s].local[f], st.local[f])) all = false;
        if (!all) continue;
        out.prob += r.pi[s];
        out.seen = true;
    }
    return out;
}

/** `getProbSysAggr()`: the joint per-class counts of every stateful node. */
template <class T>
SsaProbResult ssa_prob_sys_aggr(const qn::NetworkStruct<T>& sn, const SsaSerialRun<T>& r,
                                const qn::NetState<T>& st) {
    const std::size_t NF = sn.stateful_nodes.size();
    std::vector<std::vector<double>> want(NF);
    for (std::size_t f = 0; f < NF; ++f)
        want[f] = getters_detail::node_counts(sn, sn.stateful_nodes[f], st.local[f]);
    SsaProbResult out;
    for (std::size_t s = 0; s < r.space.size(); ++s) {
        bool all = true;
        for (std::size_t f = 0; f < NF && all; ++f)
            if (!getters_detail::same_counts(
                    getters_detail::node_counts(sn, sn.stateful_nodes[f], r.space[s].local[f]),
                    want[f]))
                all = false;
        if (!all) continue;
        out.prob += r.pi[s];
        out.seen = true;
    }
    return out;
}

/**
 * `sampleSys` and `sampleSysAggr`: the trajectory itself.
 *
 * THE STATE ON ROW i IS THE ONE OCCUPIED UNTIL `t[i]`, not the one entered at
 * it, which is the convention `run.tran_state` already records and the only one
 * under which the time weighting of `pi` and the rows printed here agree. The
 * event on the same row is the synchronization that ENDED that sojourn.
 */
template <class T>
SsaSamplePath<T> ssa_sample_sys(const qn::NetworkStruct<T>& sn, const SsaSerialRun<T>& r) {
    const std::size_t NF = sn.stateful_nodes.size();
    SsaSamplePath<T> out;
    out.t = r.tran_time;
    out.event = r.tran_sync;
    out.samples = r.samples;
    out.seed = r.seed;

    std::size_t width = 0;
    for (std::size_t f = 0; f < NF; ++f)
        width += r.space.empty() ? 0 : r.space[0].local[f].size();
    const std::size_t n = r.tran_state.size();
    out.state = Matrix<T>(n, width, num_traits<T>::from_int(0));
    out.aggr = Matrix<T>(n, r.ssq.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t s = r.tran_state[i];
        std::size_t col = 0;
        for (std::size_t f = 0; f < NF; ++f)
            for (std::size_t j = 0; j < r.space[s].local[f].size(); ++j)
                out.state(i, col++) = r.space[s].local[f][j];
        for (std::size_t c = 0; c < r.ssq.cols(); ++c) out.aggr(i, c) = r.ssq(s, c);
    }
    return out;
}

/**
 * `sample(node)` and `sampleAggr(node)`: the same trajectory, one node's block.
 *
 * The aggregate here is the node's own per-class counts and NOT a column slice
 * of `ssq`: `ssq` is indexed by STATION and a Cache, a Join or a Transition is
 * stateful without being one, so slicing would silently address the wrong node
 * on every model that has one.
 */
template <class T>
SsaSamplePath<T> ssa_sample_node(const qn::NetworkStruct<T>& sn, const SsaSerialRun<T>& r,
                                 std::size_t ind) {
    const std::size_t isf = sn.stateful_index(ind);
    if (isf == 0)
        throw InputError("ssa_sample_node: node '" + getters_detail::node_label(sn, ind) +
                         "' is stateless, so no trajectory passes through its state");
    SsaSamplePath<T> out;
    out.t = r.tran_time;
    out.event = r.tran_sync;
    out.samples = r.samples;
    out.seed = r.seed;
    const std::size_t n = r.tran_state.size();
    const std::size_t width = r.space.empty() ? 0 : r.space[0].local[isf - 1].size();
    out.state = Matrix<T>(n, width, num_traits<T>::from_int(0));
    out.aggr = Matrix<T>(n, sn.nclasses, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t s = r.tran_state[i];
        for (std::size_t j = 0; j < width; ++j) out.state(i, j) = r.space[s].local[isf - 1][j];
        const std::vector<double> c = getters_detail::node_counts(sn, ind, r.space[s].local[isf - 1]);
        for (std::size_t k = 0; k < c.size() && k < sn.nclasses; ++k)
            out.aggr(i, k) = num_traits<T>::from_double(c[k]);
    }
    return out;
}

/**
 * The whole `-a prob` report over the model's DEFAULT INITIAL STATE, which is
 * the state SolverCTMC's own `-a prob` reports on.
 *
 * SAME STATE, SAME QUESTION, DIFFERENT ESTIMATOR: the CTMC answer is the
 * stationary law of the chain and this one is a time average of a finite path,
 * so the pair is a direct measurement of the simulation error on a model small
 * enough for both. That is what makes the report worth having at the same
 * state rather than at a state of the simulator's choosing.
 */
template <class T>
SsaProbReport solver_ssa_prob(const qn::NetworkStruct<T>& sn, const SsaSerialOptions& opt) {
    const SsaSerialSolution<T> sim = solver_ssa_serial_analyzer(sn, opt);
    if (!sim.fjclassmap.empty())
        throw UnsupportedError(
            "SolverSSA: a probability query on a fork-join model would be asked in the CALLER's "
            "class space and answered in the tag-augmented one, whose sibling classes the caller "
            "never declared and cannot name a state in; SolverCTMC reports the same query on the "
            "augmented chain, so use it, or ask for -a avg here");

    ctmc::CtmcOptions copt;
    copt.cutoff = opt.cutoff;
    copt.state_max = opt.state_max;
    const qn::NetState<T> init =
        serial_detail::wide_init_state(sn, ctmc::analyzer_detail::resolve_cutoff(sn, copt));

    SsaProbReport out;
    out.samples = sim.run.samples;
    out.seed = sim.run.seed;
    out.simulated_time = sim.run.simulated_time;
    out.sys = ssa_prob_sys(sn, sim.run, init);
    out.sys_aggr = ssa_prob_sys_aggr(sn, sim.run, init);
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const std::size_t ind = sn.node_of_station(i + 1);
        const std::size_t isf = sn.stateful_index(ind);
        if (isf == 0) {
            out.marg.push_back(SsaProbResult());
            out.aggr.push_back(SsaProbResult());
            continue;
        }
        out.marg.push_back(ssa_prob(sn, sim.run, ind, init.local[isf - 1]));
        out.aggr.push_back(ssa_prob_aggr(
            sn, sim.run, ind, getters_detail::node_counts(sn, ind, init.local[isf - 1])));
    }
    return out;
}

/**
 * `getCdfRespT`: refused, and the refusal is the ANSWER rather than a gap.
 *
 * `@@SolverSSA/getCdfRespT.m` raises the same error in MATLAB, for the reason it
 * states: SSA samples state trajectories, not per-job sojourn times, so there is
 * nothing to build an empirical CDF from. The inherited implementation would
 * fabricate an exponential law with the right mean, which carries no information
 * about the tail and would be indistinguishable, to the caller, from a measured
 * distribution. Reproducing THAT would be the defect.
 */
inline void ssa_cdf_respt_refuse() {
    throw UnsupportedError(
        "SolverSSA does not record per-job response times, so it cannot return an empirical "
        "response time CDF; the inherited exponential fit would carry no information about the "
        "tail while looking measured. Use SolverJMT for a measured CDF, SolverFluid or SolverMAM "
        "for an analytical one");
}

}  // namespace ssa
}  // namespace line

#endif  // LINE_SOLVERS_SSA_SOLVER_SSA_GETTERS_H
