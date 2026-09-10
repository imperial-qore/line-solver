/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MDD_MDD_CLOSEDQN_H
#define LINE_API_MDD_MDD_CLOSEDQN_H

/**
 * Exact solve of a single-class closed exponential queueing network whose CTMC
 * state space (reachable occupancy vectors) is stored in a Multi-valued
 * Decision Diagram instead of an explicit state list.
 *
 * Port of matlab/src/api/mdd/mdd_closedqn.m, jline.api.mdd.Mdd_closedqn and
 * python/line_solver/api/mdd/closedqn.py. The reachable set is generated with
 * `mdd_reachset` and the generator matrix is assembled using the MDD's O(K)
 * state indexing (`MDD::index`), so no explicit (|S| x width) state matrix is
 * ever materialised during assembly -- the diagram is the store. For
 * single-class exponential stations the aggregated (occupancy) chain is exact:
 * the rate from n to n-e_i+e_j is mu_i * min(n_i, c_i) * P(i,j) for n_i > 0,
 * matching SolverCTMC on the same model, which makes this the live exact
 * oracle the `mdd_mcd` aggregation is validated against.
 *
 * The reference's 'verbose' knob is NOT carried (this api layer is silent, as
 * `MddMcdOptions` documents); `stats` returns the same storage numbers. The
 * reference's 'ctmcmethod' knob is NOT carried either: the C++ `ctmc_solve`
 * has a single direct backend.
 */

#include <chrono>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_reachset.h"
#include "line/api/mdd/mdd_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mdd {

/** Result of the MDD-stored exact closed-network solve. */
template <class T>
struct MddClosedQnResult {
    /** The MDD holding the reachable occupancy set. */
    MDD mdd = MDD(std::vector<int>());
    /** CTMC generator, rows aligned to MDD::index order. */
    Matrix<T> Q;
    /** Stationary distribution over the reachable states, MDD::index order. */
    std::vector<T> pi;
    /** |S| x M occupancy states, in MDD::index order. */
    std::vector<std::vector<int>> states;
    /** Mean number of jobs per station. */
    std::vector<T> QLen;
    /** Utilisation: busy servers / servers, or mean busy jobs for a delay. */
    std::vector<T> U;
    /** Per-station throughput. */
    std::vector<T> X;
    /** MDD storage statistics of the reachable set. */
    MddStats stats;
    /** Phase timings in seconds: reachable-set build (0 when the diagram was
     * supplied), generator assembly, ctmc_solve, performance measures. */
    double time_reach = 0.0, time_gen = 0.0, time_solve = 0.0, time_metrics = 0.0;
};

/**
 * @param mu per-station exponential service rates, length M
 * @param P M x M Markovian routing matrix (row-stochastic, irreducible)
 * @param servers servers per station; infinite for a delay/IS station
 * @param N closed population
 * @param reuse an already-built reachable set (the `mdd` of a previous result
 *              on the same mu/P/servers/N) to skip regeneration, or nullptr
 */
template <class T>
MddClosedQnResult<T> mdd_closedqn(const std::vector<T>& mu, const Matrix<T>& P,
                                  const std::vector<double>& servers, int N,
                                  const MDD* reuse = nullptr) {
    typedef std::chrono::steady_clock clock;
    const std::size_t M = mu.size();
    if (P.rows() != M || P.cols() != M)
        throw InputError("mdd_closedqn: routing matrix must be M x M");
    if (servers.size() != M)
        throw InputError("mdd_closedqn: servers must have one entry per station");
    if (N < 1) throw InputError("mdd_closedqn: the closed population must be positive");
    const T zero = num_traits<T>::from_int(0);

    // events: a completion at station i (rate mu_i * min(n_i, c_i)) routes to
    // station j with probability P(i,j); self-routing i == j leaves the
    // occupancy vector unchanged and is skipped.
    std::vector<std::size_t> ii, jj;
    std::vector<T> pij;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j)
            if (i != j && P(i, j) != zero) {
                ii.push_back(i);
                jj.push_back(j);
                pij.push_back(P(i, j));
            }
    const std::size_t E = ii.size();

    // all jobs start at station 1; an irreducible routing chain makes every
    // composition of N over M stations reachable.
    const std::vector<int> domain(M, N + 1);
    std::vector<int> init(M, 0);
    init[0] = N;
    const MddNextState nextfun = [&ii, &jj, E](const std::vector<int>& s) {
        std::vector<std::vector<int>> succ;
        for (std::size_t a = 0; a < E; ++a) {
            if (s[ii[a]] > 0) {
                std::vector<int> t = s;
                --t[ii[a]];
                ++t[jj[a]];
                succ.push_back(t);
            }
        }
        return succ;
    };

    MddClosedQnResult<T> out;
    std::chrono::time_point<clock> t0 = clock::now();
    if (reuse != nullptr) {
        out.mdd = *reuse;
        out.time_reach = 0.0;
    } else {
        out.mdd = mdd_reachset(domain, init, nextfun);
        out.time_reach = std::chrono::duration<double>(clock::now() - t0).count();
    }

    t0 = clock::now();
    const long long n = out.mdd.cardinality();
    out.states = out.mdd.enumerate();

    // assemble the generator directly from the MDD state indexing
    out.Q = Matrix<T>(static_cast<std::size_t>(n), static_cast<std::size_t>(n));
    for (long long s = 0; s < n; ++s) {
        const std::vector<int>& st = out.states[static_cast<std::size_t>(s)];
        const long long row = out.mdd.index(st);
        for (std::size_t a = 0; a < E; ++a) {
            const std::size_t i = ii[a];
            if (st[i] > 0) {
                const int busy = std::isinf(servers[i])
                                     ? st[i]
                                     : std::min(st[i], static_cast<int>(servers[i]));
                const T rate = mu[i] * num_traits<T>::from_int(busy) * pij[a];
                std::vector<int> t = st;
                --t[i];
                ++t[jj[a]];
                const long long col = out.mdd.index(t);
                out.Q(static_cast<std::size_t>(row), static_cast<std::size_t>(col)) += rate;
            }
        }
    }
    out.Q = mc::ctmc_makeinfgen(out.Q);
    out.time_gen = std::chrono::duration<double>(clock::now() - t0).count();

    t0 = clock::now();
    out.pi = mc::ctmc_solve(out.Q);
    out.time_solve = std::chrono::duration<double>(clock::now() - t0).count();

    // performance metrics
    t0 = clock::now();
    out.QLen.assign(M, zero);
    out.U.assign(M, zero);
    out.X.assign(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        T qi = zero, busyMean = zero;
        for (long long s = 0; s < n; ++s) {
            const int ni = out.states[static_cast<std::size_t>(s)][i];
            const int busy = std::isinf(servers[i]) ? ni : std::min(ni, static_cast<int>(servers[i]));
            qi += out.pi[static_cast<std::size_t>(s)] * num_traits<T>::from_int(ni);
            busyMean += out.pi[static_cast<std::size_t>(s)] * num_traits<T>::from_int(busy);
        }
        out.QLen[i] = qi;
        out.X[i] = mu[i] * busyMean;  // throughput = mean completion rate
        if (std::isinf(servers[i]))
            out.U[i] = qi;  // mean number busy (IS station)
        else
            out.U[i] = busyMean / num_traits<T>::from_int(static_cast<int>(servers[i]));
    }
    out.time_metrics = std::chrono::duration<double>(clock::now() - t0).count();

    out.stats = out.mdd.stats();
    return out;
}

}  // namespace mdd
}  // namespace line

#endif  // LINE_API_MDD_MDD_CLOSEDQN_H
