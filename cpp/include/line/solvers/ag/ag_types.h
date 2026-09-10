/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AG_AG_TYPES_H
#define LINE_SOLVERS_AG_AG_TYPES_H

/**
 * @file ag_types.h
 * @brief Options of the agent-based (RCAT) solver.
 *
 * Beside the truncation level of an open agent, these carry the EXECUTION
 * BACKEND of the reversed-rate fixed point. The backend decides who evaluates an
 * agent, never what the agent evaluates to: agent k reads the rest of the model
 * only through the scalar reversed rates x, and the sweep is Jacobi, so the
 * agent order is immaterial and every backend walks the same iterates.
 */

#include <cstddef>
#include <string>
#include <vector>

namespace line {
namespace ag {

/** The caller's own loop, in agent order. The reference. */
inline const char* exec_serial() { return "serial"; }
/**
 * Fan the agents of a sweep out over a local thread pool.
 *
 * Named 'parallel', with 'para' accepted as an alias -- the same pair SolverSSA's
 * replica analyzer answers to (`case {'para','parallel'}` in
 * solver_ssa_analyzer.m, `m == "para" || m == "parallel"` in
 * solver_ssa_parallel.h), so one spelling convention covers both solvers. It was
 * called 'threads' until 2026-08-19; that name is no longer accepted, and
 * `exec_is_parallel` is the only place either spelling is recognised.
 */
inline const char* exec_parallel() { return "parallel"; }
/** The accepted alias of exec_parallel(), as SolverSSA spells it. */
inline const char* exec_para() { return "para"; }
/** True for either spelling of the local-thread-pool backend. */
inline bool exec_is_parallel(const std::string& mode) {
    return mode == exec_parallel() || mode == exec_para();
}
/** Partition the agents over remote ag-worker processes. */
inline const char* exec_cluster() { return "cluster"; }

struct AgOptions {
    /** 'default', 'inap', 'inapplus', 'inapinf' or the vestigial 'exact'. */
    std::string method = "default";

    /** Convergence tolerance of the reversed-rate fixed point. */
    double tol = 1e-4;

    /** SolverOptions('AG') lowers this from the global 1000 to 100. */
    int iter_max = 100;

    /**
     * Truncation level of an OPEN agent's queue-length dimension. A closed class
     * uses its own population instead, so this bounds only the open agents;
     * 'inapinf' ignores it and solves them on the infinite state space through
     * the matrix-geometric tail.
     */
    std::size_t max_states = 100;

    /**
     * `options.config.nonmkvorder`: the phase budget `sn_nonmarkov_toph` spends
     * replacing a non-Markovian service law. The reference default is 20.
     */
    std::size_t nonmkv_order = 20;

    /** One of exec_serial(), exec_parallel() (or exec_para()), exec_cluster(). */
    std::string exec = "serial";

    /**
     * Thread-pool size for 'parallel'; 0 means one thread per available
     * processor. Pinned rather than derived per sweep so a run is reproducible
     * on a machine whose load changes under it.
     */
    unsigned nworkers = 0;

    /** Worker addresses ("host:port") for 'cluster'. */
    std::vector<std::string> endpoints;

    /**
     * Seconds to wait on a worker before solving its agents locally instead. A
     * lost worker is never fatal: any agent can be solved anywhere given x, so a
     * cluster run degrades to a slower run and never to a wrong one.
     */
    double worker_timeout = 30.0;
};

}  // namespace ag
}  // namespace line

#endif  // LINE_SOLVERS_AG_AG_TYPES_H
