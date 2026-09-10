/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_MAM_TYPES_H
#define LINE_SOLVERS_MAM_MAM_TYPES_H

/**
 * The option and result types SolverMAM shares with its analyzers.
 *
 * The solution type IS `mva::MvaSolution<T>`: the MAM analyzers return the same
 * [Q,U,R,T,C,X] tuple, and `@@NetworkSolver/getAvg` applies the same
 * solver-independent metric filter afterwards, so reusing it is what keeps the
 * two runners from drifting apart.
 */

#include <cstddef>
#include <limits>
#include <string>

#include "line/solvers/mva/mva_types.h"

namespace line {
namespace mam {

/** The options SolverMAM reads. Defaults are SolverOptions('MAM'). */
struct MamOptions {
    std::string method = "default";
    double tol = 1e-4;
    /** SolverOptions('MAM') lowers this from the global 1000 to 100. */
    int iter_max = 100;
    /**
     * `options.config.space_max`: the order budget of the per-station arrival
     * superposition, handed to mmap_super_safe. `dec.poisson` is exactly this
     * set to 1, which collapses every stream to marked Poisson.
     */
    std::size_t space_max = 128;
    /**
     * `options.config.preserveDet`: keep a Det service as Det so the exact
     * MAP/D/c solver can claim it, instead of Erlang-approximating it. The
     * analyzer turns this on when the caller has not.
     */
    bool preserve_det = true;
    /**
     * `options.config.nonmkvorder`: the phase budget `sn_nonmarkov_toph` spends
     * replacing a non-Markovian service law. The reference default is 20.
     */
    std::size_t nonmkv_order = 20;
    /**
     * `options.config.phfit`: true selects the two-moment concentrated-ME fit
     * (the reference default), false the Bernstein density fit. MAM consumes a
     * matrix exponential happily, so the default stands here.
     */
    bool phfit_cme = true;
    /**
     * `options.config.num_cdf_pts`: how many points the response-time CDF is
     * evaluated at. The GLOBAL default is 200, not the 100 that
     * `solver_mam_passage_time.m`'s own fallback names -- the fallback is dead,
     * because SolverOptions always sets the field. Measured: MATLAB returns 200
     * rows for a model whose `options.config` was never touched.
     */
    std::size_t num_cdf_pts = 200;
    /**
     * `options.cutoff`: the level truncation `getProb` / `getProbMarg` use for
     * an OPEN model, where the queue length is unbounded. 0 selects the
     * reference's default of 100.
     */
    std::size_t cutoff = 0;
    /**
     * `options.config.bgstates_max`: the cap on the number of states of the
     * bgchain background chain, which enumerates the closed-class population
     * vector and so grows as nchoosek(N+Mc-1,Mc-1) per background class. 0
     * selects the default of 20000. Exceeding it is a named error, not a silent
     * degradation.
     */
    std::size_t bgstates_max = 0;
    /**
     * `options.config.bgaggr`: the number G of AGGREGATE background classes the
     * bgchain method carries, so the chain has 1 + G classes. G = 1 is the
     * classic tagged/aggregate pair and the default; G >= R-1 aggregates nothing
     * and carries every closed chain exactly. 0 selects the default of 1.
     */
    std::size_t bgaggr = 0;
    /**
     * `options.config.qbdphases_max`: the cap on the phase count of one bgchain
     * station QBD, the product (arrival order) x (environment states) x (service
     * phases). 0 selects the default of 500.
     */
    std::size_t qbdphases_max = 0;
    /**
     * `options.config.fj_accuracy`: the FJ_codes truncation C of the
     * queue-length DIFFERENCE between the two fork-join branches. Larger is
     * more accurate and costs (C + 1) m^2 ma in every matrix dimension; the
     * reference default is 100.
     */
    std::size_t fj_accuracy = 100;
    /**
     * `options.config.fj_tmode`: which route `computeT.m` takes to the T
     * matrix, 'NARE' (Riccati, the reference default) or 'Sylves' (the
     * fixed-point iteration of Section 5.1).
     */
    std::string fj_tmode = "NARE";
    /**
     * `options.config.timescale`: "auto", "discrete" or "continuous". "auto"
     * lets the distributions decide whether the model is slotted; the other two
     * force the reading and raise when the model does not admit it.
     */
    std::string timescale = "auto";
    /** `options.config.slotlength`: the slot in model time units. */
    double slotlength = 1.0;
    /**
     * `options.timespan`: the transient horizon. SolverOptions('MAM') defaults
     * it to [Inf, Inf], i.e. steady state only, and `runAnalyzer` treats a
     * finite upper bound as the request for a transient solve.
     */
    double timespan_start = 0.0;
    double timespan_end = std::numeric_limits<double>::infinity();
};

/** What the MAM dispatch returns: the metrics plus the algorithm that ran. */
template <class T>
struct MamSolution {
    mva::MvaSolution<T> sol;
    /** The concrete algorithm, as the reference's `actualmethod`. */
    std::string actualmethod;
};

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_MAM_TYPES_H
