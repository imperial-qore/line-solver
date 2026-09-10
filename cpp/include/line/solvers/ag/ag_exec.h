/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AG_AG_EXEC_H
#define LINE_SOLVERS_AG_AG_EXEC_H

/**
 * @file ag_exec.h
 * @brief Execution backends of the reversed-rate fixed point.
 *
 * WHAT MAKES THIS SAFE IS THE DECOMPOSITION, NOT THE SCHEDULING. Agent k's
 * generator is
 *
 *     Q_k(x) = L_k + sum_{c passive at k} x_c Pb_c
 *
 * so an agent reads the rest of the model only through the scalar reversed rates
 * x, and it writes only its own slot of the sweep's output. The sweep is Jacobi
 * -- every x_a is read off the PREVIOUS sweep's stationary vectors and only then
 * do the agents re-solve -- so the agent order is immaterial, and a parallel or
 * distributed sweep produces the SAME iterates as the serial one.
 *
 * HOW FAR THAT SURVIVES FLOATING POINT DEPENDS ON WHO RUNS THE AGENT SOLVE.
 * `parallel` is BIT-IDENTICAL to `serial`: the same code in the same process,
 * differing only in an order that does not matter. `cluster` is bit-identical
 * only when the worker runs the same implementation as the coordinator -- the
 * wire is exact, since JSON round-trips a double without loss, but the
 * stationary vector comes back from the WORKER's solve, so a C++ coordinator
 * driving a Java ag-worker agrees to a few ulp rather than bit for bit. That is
 * the ordinary cross-codebase difference, not a protocol defect.
 *
 * ARITHMETIC. `parallel` carries any arithmetic the agent solve carries, because
 * it moves no numbers between representations. `cluster` is DOUBLE ONLY and
 * refuses by name otherwise: the wire is JSON, so an exact rational or a
 * 200-digit float would have to be rounded to send, and silently answering a
 * high-precision request at double precision is worse than refusing it.
 */

#include <algorithm>
#include <cstddef>
#include <string>
#include <thread>
#include <vector>

#include "line/num/number.h"
#include "line/solvers/ag/ag_types.h"
#include "line/solvers/ag/ag_worker_client.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ag {

/**
 * Evaluate every agent of one sweep.
 *
 * @param n         number of agents
 * @param gen       k -> agent k's generator at the current reversed rates
 * @param sol       (Q_k, k) -> agent k's stationary vector
 * @param Q,pi      per-agent output slots, resized by the caller
 */
template <class T, class Gen, class Sol>
void ag_sweep_serial(std::size_t n, Gen gen, Sol sol,
                     std::vector<Matrix<T>>& Q, std::vector<std::vector<T>>& pi) {
    for (std::size_t k = 0; k < n; ++k) {
        Q[k] = gen(k);
        pi[k] = sol(Q[k], k);
    }
}

/**
 * The same sweep over a thread pool. Each task owns one agent and writes only
 * its own slot, so no synchronisation beyond the join is needed and the result
 * cannot depend on the interleaving.
 */
template <class T, class Gen, class Sol>
void ag_sweep_parallel(std::size_t n, unsigned nworkers, Gen gen, Sol sol,
                      std::vector<Matrix<T>>& Q, std::vector<std::vector<T>>& pi) {
    unsigned hw = nworkers > 0 ? nworkers : std::thread::hardware_concurrency();
    if (hw == 0) hw = 1;
    const std::size_t nthreads = std::min<std::size_t>(hw, n == 0 ? 1 : n);
    if (nthreads <= 1) {
        ag_sweep_serial<T>(n, gen, sol, Q, pi);
        return;
    }

    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    // A STRIDED partition, not a contiguous one: agent cost varies with the
    // agent's state-space size, and neighbouring agents are the same station's
    // classes and so are similarly sized. Striding mixes big and small across
    // the threads instead of loading one thread with a station's whole block.
    for (std::size_t t = 0; t < nthreads; ++t) {
        pool.emplace_back([&, t]() {
            for (std::size_t k = t; k < n; k += nthreads) {
                Q[k] = gen(k);
                pi[k] = sol(Q[k], k);
            }
        });
    }
    for (std::size_t t = 0; t < pool.size(); ++t) pool[t].join();
}

/**
 * The same sweep with the agents partitioned over remote ag-worker processes.
 *
 * The generator is rebuilt here in any case -- the metrics stage reads it, and
 * assembling it is O(N^2) against the O(N^3) solve -- so only the stationary
 * vector crosses the wire back. An unreachable, slow or broken worker is not
 * fatal: its agents fall through to @p sol on this process.
 */
template <class T, class Gen, class Sol, class Payload>
void ag_sweep_cluster(std::size_t n, AgWorkerPool& workers, const std::vector<T>& x,
                      Gen gen, Sol sol, Payload payload,
                      std::vector<Matrix<T>>& Q, std::vector<std::vector<T>>& pi) {
    if constexpr (!std::is_same<T, double>::value) {
        (void)workers; (void)x; (void)payload;
        throw UnsupportedError(
            "ag: the 'cluster' execution backend is double only, because the ag-worker "
            "protocol is JSON and an exact or high-precision value would have to be rounded "
            "to send it; rerun with --arith double, or use exec 'serial' or 'parallel'");
    } else {
        workers.ensure_assigned(n, payload);

        for (std::size_t k = 0; k < n; ++k) Q[k] = gen(k);

        std::vector<bool> pending(n, true);
        std::vector<double> xd(x.size());
        for (std::size_t c = 0; c < x.size(); ++c) xd[c] = x[c];

        workers.sweep(xd, pi, pending);

        for (std::size_t k = 0; k < n; ++k) {
            if (pending[k]) pi[k] = sol(Q[k], k);
        }
    }
}

}  // namespace ag
}  // namespace line

#endif  // LINE_SOLVERS_AG_AG_EXEC_H
