/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_BND_LR_PF_H
#define LINE_API_MAPQN_MAPQN_BND_LR_PF_H

/**
 * Linear-reduction (LR) bound on the utilization of one station of a closed
 * product-form network.
 *
 * Templated port of matlab/lib/qrf/mapqn_bnd_lr_pf.m (ground truth), whose own
 * origin is the AMPL model bnd_linearreduction_pf.mod.
 *
 * WHAT MAKES IT A BOUND. The exact stationary distribution satisfies every
 * constraint assembled below -- normalization, the utilization and queue-length
 * definitions, population conservation, the global flow balance of Gordon-Newell
 * and the joint-marginal identities -- but those constraints do not pin it down.
 * The feasible set is therefore a POLYTOPE CONTAINING the exact solution, so
 * minimizing U_i over it cannot exceed the true utilization and maximizing it
 * cannot fall short. Both senses are valid bounds, which is why one routine
 * serves `lr.lower` and `lr.upper`.
 *
 * NOT TO BE CONFUSED WITH `qrf.mmi.linear`, whose name refers only to its
 * explicit Aeq/beq constraint representation: its OBJECTIVE is the nonlinear
 * MEM entropy and it needs `fmincon`. This method is an LP end to end, which is
 * exactly why it ports and the rest of the QRF family does not.
 *
 * ARITHMETIC. Assembly is additions and multiplications, and the optimum is a
 * vertex of a rational polytope, so at `Rational` the bound is EXACT and is
 * deliberately left ungated: `lp::simplex_solve` uses Bland's rule with no
 * tolerance under exact arithmetic. MATLAB reaches the same vertex through
 * `linprog`'s interior-point method, which approaches it from the interior and
 * stops a few digits short; a deviation against MATLAB is therefore expected to
 * be MATLAB's convergence gap and not this port's error.
 *
 * COST. The variable count is 2M + M^2 + 2 M^2 (N+1) and the row count is
 * O(M^2 + M N), so the LP grows quadratically in the station count and linearly
 * in the population. It is a bound, not a cheap one.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/simplex.h"
#include "line/api/mapqn/mapqn_params.h"

namespace line {
namespace mapqn {

/** Product-form parameters of the LR bound, mirroring MATLAB's `params`. */
template <class T>
struct LrPfParams {
    int M = 0;             ///< station count
    int N = 0;             ///< closed population
    std::vector<T> mu;     ///< (M) service rates
    Matrix<T> r;           ///< (M x M) routing probabilities
};

/** Return value of mapqn_bnd_lr_pf, mirroring the MATLAB `result` struct. */
template <class T>
struct LrPfResult {
    T objective = T();     ///< the bounded utilization of the objective station
    std::vector<T> U;      ///< (M) utilizations at the optimal vertex
    std::vector<T> Q;      ///< (M) queue lengths at the optimal vertex
    std::string status;    ///< simplex status name
    bool ok = false;
    std::size_t iterations = 0;
    std::size_t num_vars = 0, num_rows = 0;
};

namespace detail {

/**
 * The variable layout of the reference, in its declaration order: U, Q, C, p1,
 * p1c. The order is load bearing only in that the tests compare against MATLAB
 * through the RESULT fields, but keeping it identical makes the two assemblies
 * diffable line by line.
 */
struct LrPfIndex {
    int M = 0, N = 0;
    std::size_t u0 = 0, q0 = 0, c0 = 0, p1_0 = 0, p1c0 = 0, total = 0;

    LrPfIndex(int M_, int N_) : M(M_), N(N_) {
        const std::size_t m = static_cast<std::size_t>(M);
        const std::size_t np1 = static_cast<std::size_t>(N + 1);
        u0 = 0;
        q0 = u0 + m;
        c0 = q0 + m;
        p1_0 = c0 + m * m;
        p1c0 = p1_0 + m * m * np1;
        total = p1c0 + m * m * np1;
    }
    std::size_t U(int i) const { return u0 + static_cast<std::size_t>(i); }
    std::size_t Q(int i) const { return q0 + static_cast<std::size_t>(i); }
    std::size_t C(int j, int i) const {
        return c0 + static_cast<std::size_t>(j) * static_cast<std::size_t>(M) +
               static_cast<std::size_t>(i);
    }
    std::size_t p1(int j, int i, int n) const {
        return p1_0 +
               (static_cast<std::size_t>(j) * static_cast<std::size_t>(M) +
                static_cast<std::size_t>(i)) *
                   static_cast<std::size_t>(N + 1) +
               static_cast<std::size_t>(n);
    }
    std::size_t p1c(int j, int i, int n) const {
        return p1c0 +
               (static_cast<std::size_t>(j) * static_cast<std::size_t>(M) +
                static_cast<std::size_t>(i)) *
                   static_cast<std::size_t>(N + 1) +
               static_cast<std::size_t>(n);
    }
};

}  // namespace detail

/**
 * Port of `mapqn_bnd_lr_pf`.
 *
 * @param p               station count, population, rates and routing
 * @param objective_queue 1-BASED station whose utilization is bounded, as in the
 *                        reference
 * @param sense           Min for the `lr.lower` bound, Max for `lr.upper`
 */
template <class T>
LrPfResult<T> mapqn_bnd_lr_pf(const LrPfParams<T>& p, int objective_queue, MapqnSense sense) {
    const int M = p.M, N = p.N;
    if (M <= 0) throw InputError("mapqn_bnd_lr_pf: the station count must be positive");
    if (N <= 0) throw InputError("mapqn_bnd_lr_pf: the population must be positive");
    if (static_cast<int>(p.mu.size()) != M)
        throw InputError("mapqn_bnd_lr_pf: mu has the wrong length");
    if (static_cast<int>(p.r.rows()) != M || static_cast<int>(p.r.cols()) != M)
        throw InputError("mapqn_bnd_lr_pf: the routing matrix is not (M x M)");
    if (objective_queue < 1 || objective_queue > M)
        throw InputError("mapqn_bnd_lr_pf: the objective station is out of range");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T Nt = num_traits<T>::from_int(N);

    // q(i,j) = r(i,j) mu(i): the rate at which station i sends work to j.
    Matrix<T> q(static_cast<std::size_t>(M), static_cast<std::size_t>(M), zero);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            q(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) =
                T(p.r(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) *
                  p.mu[static_cast<std::size_t>(i)]);

    const detail::LrPfIndex ix(M, N);
    lp::LpModel<T> m(ix.total);

    for (int i = 0; i < M; ++i) {
        m.set_bounds(ix.U(i), zero, one);
        m.set_bounds(ix.Q(i), zero, Nt);
        for (int j = 0; j < M; ++j) m.set_bounds(ix.C(j, i), zero, Nt);
    }
    for (int j = 0; j < M; ++j)
        for (int i = 0; i < M; ++i)
            for (int n = 0; n <= N; ++n) {
                m.set_bounds(ix.p1(j, i, n), zero, one);
                m.set_bounds(ix.p1c(j, i, n), zero, one);
            }

    // ZER1 and ZER3 are BOUNDS in the reference, not rows: a station observed
    // from itself is never empty (the observer is there), and a station observed
    // from a different one never holds the whole population.
    for (int j = 0; j < M; ++j) m.set_upper(ix.p1(j, j, 0), zero);
    for (int j = 0; j < M; ++j)
        for (int i = 0; i < M; ++i)
            if (j != i) m.set_upper(ix.p1(j, i, N), zero);

    // CEQU: the self-conditional queue length is the queue length.
    for (int j = 0; j < M; ++j) {
        m.row_clear();
        m.row_add(ix.C(j, j), one);
        m.row_add(ix.Q(j), T(-one));
        m.emit(lp::LpSense::EQ, zero);
    }
    // ONE1: each conditional marginal and its complement sum to one.
    for (int j = 0; j < M; ++j)
        for (int i = 0; i < M; ++i) {
            m.row_clear();
            for (int n = 0; n <= N; ++n) {
                m.row_add(ix.p1(j, i, n), one);
                m.row_add(ix.p1c(j, i, n), one);
            }
            m.emit(lp::LpSense::EQ, one);
        }
    // UTIL: the utilization is the mass of the conditional marginal, whichever
    // station it is conditioned on -- which is what ties the M copies together.
    for (int i = 0; i < M; ++i)
        for (int t = 0; t < M; ++t) {
            m.row_clear();
            m.row_add(ix.U(i), one);
            for (int n = 0; n <= N; ++n) m.row_add(ix.p1(i, t, n), T(-one));
            m.emit(lp::LpSense::EQ, zero);
        }
    // QLEN and CLEN: first moments of the marginals.
    for (int i = 0; i < M; ++i) {
        m.row_clear();
        m.row_add(ix.Q(i), one);
        for (int n = 0; n <= N; ++n) m.row_add_int(ix.p1(i, i, n), -static_cast<long>(n));
        m.emit(lp::LpSense::EQ, zero);
    }
    for (int j = 0; j < M; ++j)
        for (int i = 0; i < M; ++i) {
            m.row_clear();
            m.row_add(ix.C(j, i), one);
            for (int n = 0; n <= N; ++n) m.row_add_int(ix.p1(j, i, n), -static_cast<long>(n));
            m.emit(lp::LpSense::EQ, zero);
        }
    // MPCB: the conditional queue lengths seen from j sum to N U_j.
    for (int j = 0; j < M; ++j) {
        m.row_clear();
        for (int i = 0; i < M; ++i) m.row_add(ix.C(j, i), one);
        m.row_add(ix.U(j), T(-Nt));
        m.emit(lp::LpSense::EQ, zero);
    }
    // POPC: the population is conserved.
    m.row_clear();
    for (int i = 0; i < M; ++i) m.row_add(ix.Q(i), one);
    m.emit(lp::LpSense::EQ, Nt);

    // GFFL0 and GFFL: Gordon-Newell global flow balance at each level. The
    // reference ACCUMULATES into the row, so the p1(i,i,.) coefficient collects
    // -sum_{j != i} q(i,j) rather than being overwritten; row_add mirrors that.
    for (int i = 0; i < M; ++i) {
        m.row_clear();
        for (int j = 0; j < M; ++j) {
            if (j == i) continue;
            m.row_add(ix.p1(j, i, 0), q(static_cast<std::size_t>(j), static_cast<std::size_t>(i)));
            m.row_add(ix.p1(i, i, 1),
                      T(-q(static_cast<std::size_t>(i), static_cast<std::size_t>(j))));
        }
        m.emit(lp::LpSense::EQ, zero);
    }
    for (int i = 0; i < M; ++i)
        for (int n = 1; n <= N - 1; ++n) {
            m.row_clear();
            for (int j = 0; j < M; ++j) {
                if (j == i) continue;
                m.row_add(ix.p1(j, i, n),
                          q(static_cast<std::size_t>(j), static_cast<std::size_t>(i)));
                m.row_add(ix.p1(i, i, n + 1),
                          T(-q(static_cast<std::size_t>(i), static_cast<std::size_t>(j))));
            }
            m.emit(lp::LpSense::EQ, zero);
        }
    // UJNT: the joint occupancy of (i,j) is symmetric in the two orderings.
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            m.row_clear();
            for (int n = 1; n <= N; ++n) m.row_add(ix.p1(j, i, n), one);
            for (int n = 1; n <= N; ++n) m.row_add(ix.p1(i, j, n), T(-one));
            m.emit(lp::LpSense::EQ, zero);
        }
    // QBAL: flow balance of station i against every other station.
    for (int i = 0; i < M; ++i) {
        m.row_clear();
        for (int j = 0; j < M; ++j) {
            if (j == i) continue;
            m.row_add(ix.U(i), q(static_cast<std::size_t>(i), static_cast<std::size_t>(j)));
            for (int n = 1; n <= N; ++n)
                m.row_add(ix.p1(i, j, n),
                          T(-q(static_cast<std::size_t>(j), static_cast<std::size_t>(i))));
            m.row_add(ix.p1(j, i, 0),
                      T(-q(static_cast<std::size_t>(j), static_cast<std::size_t>(i))));
        }
        m.emit(lp::LpSense::EQ, zero);
    }

    m.set_cost(ix.U(objective_queue - 1), one);
    m.set_maximize(sense == MapqnSense::Max);

    LrPfResult<T> res;
    res.num_vars = m.num_vars();
    res.num_rows = m.num_rows();
    const lp::LpSolution<T> s = lp::simplex_solve(m);
    res.status = lp::lp_status_name(s.status);
    res.iterations = s.iterations;
    res.ok = s.ok();
    if (!res.ok) return res;
    res.objective = s.objective;
    res.U.resize(static_cast<std::size_t>(M));
    res.Q.resize(static_cast<std::size_t>(M));
    for (int i = 0; i < M; ++i) {
        res.U[static_cast<std::size_t>(i)] = s.x[ix.U(i)];
        res.Q[static_cast<std::size_t>(i)] = s.x[ix.Q(i)];
    }
    return res;
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_BND_LR_PF_H
