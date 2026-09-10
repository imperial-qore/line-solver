/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of the `@@SolverCTMC` sampling surface: `sample`, `sampleAggr`,
 * `sampleSys`, `sampleSysAggr`.
 *
 * WHAT IS BEING SAMPLED, and why it is not simulation in the SSA sense: the
 * chain has already been BUILT and solved, so a sample path here is a walk on a
 * generator that is exactly right, not a Monte Carlo estimate of one. The
 * randomness is in the path, not in the model, and a longer run buys a longer
 * trace rather than a more accurate answer.
 *
 * The reference builds the walk by treating the whole generator as a MARKED
 * Markovian arrival process -- D1 = sum_a filt[a], D0 = Q - D1, one mark per
 * synchronization -- and calling `mmap_sample`. The mark is what makes the
 * trace usable: it says WHICH event fired, which the state sequence alone does
 * not determine, because two synchronizations can carry the chain between the
 * same pair of states. This port samples that MMAP directly rather than
 * routing through a general MMAP sampler, which is the same walk with the
 * per-event decomposition already in hand.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_SAMPLE_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_SAMPLE_H

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

/** One sampled trajectory of the chain. */
template <class T>
struct CtmcSamplePath {
    std::vector<T> t;                     ///< time at which each state was ENTERED
    std::vector<std::size_t> state;       ///< 0-based index into `chain.space`
    std::vector<std::size_t> event;       ///< synchronization that fired to LEAVE it
    CtmcSolution<T> chain;
};

/**
 * Port of `@@SolverCTMC/sampleSys`: a marked walk on the whole network state.
 *
 * Requires the event filtration, so `opt.keep_filtration` is forced on -- the
 * mark cannot be recovered from Q, whose entries have already summed every
 * synchronization's contribution.
 *
 * @param nevents number of transitions to draw
 * @param seed    the stream; two runs are the same trace only if this matches
 * @param sn the refreshed network struct
 * @param opt_in CTMC options (state-space cutoff, tolerances, method)
 */
template <class T>
CtmcSamplePath<T> solver_ctmc_sample_sys(const NetworkStruct<T>& sn, const CtmcOptions& opt_in,
                                         std::size_t nevents, unsigned long seed = 23000) {
    static_assert(num_traits<T>::has_transcendental,
                  "solver_ctmc_sample_sys draws exponential holding times as -mean*log(u), which "
                  "is transcendental; use --arith double or real");
    assert_phase_type_states(sn, "sampleSys");

    CtmcOptions opt = opt_in;
    opt.keep_filtration = true;
    CtmcSamplePath<T> out;
    out.chain = solver_ctmc_analyzer(sn, opt);
    const std::size_t n = out.chain.chain.space.size();
    const std::size_t A = out.chain.chain.filt.size();

    // The walk starts at the model's default initial state, as the reference
    // does with `pi0(matchrow(stateSpace,s0)) = 1`. A missing initial state is
    // an error here rather than a fallback: a trace has to start somewhere
    // specific for its event sequence to mean anything.
    std::size_t cur = analyzer_detail::init_state_index(sn, out.chain.chain.space);
    if (cur == static_cast<std::size_t>(-1))
        throw InputError(
            "solver_ctmc_sample_sys: the initial state is not contained in the state space, so "
            "there is no state to start the trace from");

    pfqn::McRng rng(seed);
    T now = num_traits<T>::from_int(0);
    for (std::size_t step = 0; step < nevents; ++step) {
        out.t.push_back(now);
        out.state.push_back(cur);

        // The total exit rate is minus the diagonal, which is what makes the
        // holding time exponential with that mean.
        const double exit = -num_traits<T>::to_double(out.chain.chain.Q(cur, cur));
        if (!(exit > 0)) {
            // An absorbing state: the chain stays there forever, so the trace
            // ends rather than being padded with fictitious transitions.
            out.event.push_back(static_cast<std::size_t>(-1));
            break;
        }
        now = T(now + num_traits<T>::from_double(
                          -std::log(num_traits<T>::to_double(pfqn::mc_uniform<T>(rng))) / exit));

        // Draw the (event, destination) pair proportionally to the rate each
        // synchronization contributes -- the MMAP mark and the jump together,
        // because they are not independent.
        const double u = num_traits<T>::to_double(pfqn::mc_uniform<T>(rng)) * exit;
        double acc = 0;
        std::size_t pick_a = static_cast<std::size_t>(-1), pick_s = cur;
        for (std::size_t a = 0; a < A && pick_a == static_cast<std::size_t>(-1); ++a)
            for (std::size_t j = 0; j < n; ++j) {
                if (j == cur) continue;  // a self-loop leaves the chain where it is
                const double w = num_traits<T>::to_double(out.chain.chain.filt[a](cur, j));
                if (w <= 0) continue;
                acc += w;
                if (acc >= u) {
                    pick_a = a;
                    pick_s = j;
                    break;
                }
            }
        // Rounding can leave `u` just past the accumulated total; fall back to
        // the last positive entry rather than stalling the walk.
        if (pick_a == static_cast<std::size_t>(-1))
            for (std::size_t a = 0; a < A; ++a)
                for (std::size_t j = 0; j < n; ++j)
                    if (j != cur && num_traits<T>::to_double(out.chain.chain.filt[a](cur, j)) > 0) {
                        pick_a = a;
                        pick_s = j;
                    }
        out.event.push_back(pick_a);
        cur = pick_s;
    }
    return out;
}

/**
 * Port of `@@SolverCTMC/sampleSysAggr`: the same walk, reported as per-(station,
 * class) job counts rather than as detailed states.
 *
 * @return row `i` is the aggregate state at `path.t[i]`, in `(ist-1)*K + k`
 *         column order
 */
template <class T>
Matrix<T> solver_ctmc_sample_sys_aggr(const NetworkStruct<T>& sn,
                                      const CtmcSamplePath<T>& path) {
    const Matrix<T> A = ctmc_state_space_aggr(sn, path.chain.chain.space);
    Matrix<T> out(path.state.size(), A.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < path.state.size(); ++i)
        for (std::size_t c = 0; c < A.cols(); ++c) out(i, c) = A(path.state[i], c);
    return out;
}

/**
 * Port of `@@SolverCTMC/sample`: the walk restricted to ONE stateful node's
 * local block.
 *
 * @param ind 1-based node index
 * @param sn the refreshed network struct
 * @param path sample path to label
 * @return row `i` is that node's local state at `path.t[i]`
 */
template <class T>
Matrix<T> solver_ctmc_sample(const NetworkStruct<T>& sn, const CtmcSamplePath<T>& path,
                             std::size_t ind) {
    const std::size_t isf = sn.stateful_index(ind);
    if (isf == 0) throw InputError("solver_ctmc_sample: node " + std::to_string(ind) +
                                   " is not stateful, so it has no state to sample");
    const std::size_t w = path.chain.chain.space[0].local[isf - 1].size();
    Matrix<T> out(path.state.size(), w, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < path.state.size(); ++i)
        for (std::size_t c = 0; c < w; ++c)
            out(i, c) = path.chain.chain.space[path.state[i]].local[isf - 1][c];
    return out;
}

/** Port of `@@SolverCTMC/sampleAggr`: one node's per-class counts over time. */
template <class T>
Matrix<T> solver_ctmc_sample_aggr(const NetworkStruct<T>& sn, const CtmcSamplePath<T>& path,
                                  std::size_t ind) {
    const std::size_t isf = sn.stateful_index(ind);
    if (isf == 0) throw InputError("solver_ctmc_sample_aggr: node " + std::to_string(ind) +
                                   " is not stateful, so it has no state to sample");
    const std::size_t K = sn.nclasses;
    Matrix<T> out(path.state.size(), K, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < path.state.size(); ++i) {
        const std::vector<T> m = prob_detail::marginal_of(
            sn, ind, path.chain.chain.space[path.state[i]].local[isf - 1]);
        for (std::size_t k = 0; k < K && k < m.size(); ++k) out(i, k) = m[k];
    }
    return out;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_SAMPLE_H
