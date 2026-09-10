/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_BND_LR_H
#define LINE_API_MAPQN_MAPQN_BND_LR_H

/**
 * General linear-reduction (LR) bound on the utilization of one queue-phase of
 * a closed MAP queueing network.
 *
 * Templated port of matlab/lib/qrf/mapqn_bnd_lr.m (ground truth), cross-checked
 * against python/line_solver/api/mapqn/bnd_lr.py. Note the MATLAB file's own
 * header records that it is a port OF the Python file, so the two are one
 * lineage rather than two independent derivations; where they disagree, neither
 * is automatically right and the AMPL model noblo_skel.mod settles it.
 *
 * WHAT DISTINGUISHES IT FROM THE QUADRATIC REDUCTION. This model keeps only the
 * singly-indexed p1(j,k,i,ni,h) and p1c, never the joint p2, so it is the
 * cheaper relaxation: the column count is linear in (N+1) sum_i K(i) rather
 * than quadratic. Fifteen of its twenty families are literally the quadratic
 * model's, and live in mapqn_p1_common.h. The five that are its own replace
 * what the p2 level would otherwise supply:
 *   MPCB   sum_i C(j,k,i) = N U(j,k), the aggregate of the QR model's THM2
 *   UJNT   AMPL SIMMETRY projected onto p1, standing in for PI23
 *   GFFL0  AMPL THM30 on p1, the empty-station level crossing
 *   GFFL   AMPL THM3 on p1, the level crossing between n_i and n_i + 1
 *   QBAL   throughput balance at each station
 *
 * THE FAILURE MODE OF THIS FAMILY IS A VACUOUS BOUND, NOT A CRASH. GFFL and
 * GFFL0 are the only families that read the rates into the p1 variables, and
 * SRVB the only one that reads them into U. Drop them and no constraint
 * distinguishes a fast station from a slow one, so U = 0 stays feasible and the
 * routine returns the [0,1] box with every subscript internally consistent.
 * That is exactly the 2026-07-20 defect in mapqn_bnd_lr.m, which built q and
 * then referenced it in NO constraint at all. Each family is therefore emitted
 * by a function named after it; do not inline them.
 *
 * ORACLE. At K(i) == 1 for every i the phase structure disappears and this
 * model must reproduce mapqn_bnd_lr_pf exactly, on the same instance and in the
 * same sense. That is the check test_mapqn_bnd_lr.cpp leads with, and it is
 * what the MATLAB reference itself uses (it holds there to 1e-7).
 *
 * ARITHMETIC. Assembly is +, -, * on the model data and lp::simplex_solve uses
 * Bland's rule with no tolerance, so at T = line::Rational the returned bound
 * is the EXACT optimum of the exact polytope. MATLAB reaches it with linprog's
 * 'interior-point' and lands a few digits short.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mapqn/mapqn_p1_common.h"
#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/simplex.h"

namespace line {
namespace mapqn {

/** Result of a general LR bound solve. */
template <class T>
struct MapqnBndLrResult {
    bool ok = false;      ///< the LP reached an optimal vertex
    std::string status;   ///< textual LP status
    T objective = T();    ///< the bound on U(objective_queue, objective_phase)
    Matrix<T> U;          ///< M x max(K) utilizations, 0 beyond K(i)
    Matrix<T> IT;         ///< M x max(K) idle times, 0 beyond K(i)
    Matrix<T> Q;          ///< M x max(K) mean queue lengths, 0 beyond K(i)
    std::vector<T> x;     ///< full solution vector, indexed by LrIndex
    std::size_t num_vars = 0;
    std::size_t num_rows = 0;
    std::size_t iterations = 0;
};

/**
 * Variable layout of the LR model: the shared p1-level blocks and nothing else.
 * Defined in mapqn_p1_common.h, constructed here with the joint block switched
 * off.
 */
using LrIndex = MapqnP1Index;

namespace detail {

/**
 * MPCB: sum over i of C(j,k,i) = N U(j,k).
 *
 * The aggregate the QR model derives from THM2 by summing over nj >= 1. Here it
 * is imposed directly, because there is no p2 to sum.
 */
template <class T>
void lr_mpcb(const MapqnParams<T>& p, const LrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int i = 0; i < p.M; ++i) m.row_add_int(x.C(j, k, i), 1);
            m.row_add_int(x.U(j, k), -p.N);
            m.emit_eq_int(0);
        }
    }
}

/**
 * UJNT: joint-probability symmetry (AMPL SIMMETRY, projected onto p1).
 *
 * sum{ni>=1} p1(j,k,i,ni,h) = sum{nj>=1} p1(i,h,j,nj,k): both sides are
 * P(n_j >= 1, phase_j = k, n_i >= 1, phase_i = h). Emitted once per unordered
 * pair, as the reference does; the reverse pair is the same row negated.
 */
template <class T>
void lr_ujnt(const MapqnParams<T>& p, const LrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int i = 0; i < p.M; ++i) {
                for (int h = 0; h < p.K[i]; ++h) {
                    if (!((j < i) || (j == i && k < h))) continue;
                    for (int ni = 1; ni <= p.N; ++ni) m.row_add_int(x.p1(j, k, i, ni, h), 1);
                    for (int nj = 1; nj <= p.N; ++nj) m.row_add_int(x.p1(i, h, j, nj, k), -1);
                    m.emit_eq_int(0);
                }
            }
        }
    }
}

/**
 * GFFL0: level-crossing balance at an empty station, per arrival phase (AMPL
 * THM30). The rate into {n_i = 0, phase_i = u} from a busy neighbour equals the
 * rate out of {n_i = 1} through a completion at i.
 */
template <class T>
void lr_gffl0(const MapqnParams<T>& p, const LrIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int u = 0; u < p.K[i]; ++u) {
            for (int j = 0; j < p.M; ++j) {
                if (j == i) continue;
                for (int k = 0; k < p.K[j]; ++k)
                    for (int h = 0; h < p.K[j]; ++h)
                        m.row_add(x.p1(j, k, i, 0, u), qr_rate(p, j, i, k, h));
                for (int k = 0; k < p.K[i]; ++k)
                    m.row_add(x.p1(i, k, i, 1, k), T(-qr_rate(p, i, j, k, u)));
            }
            m.emit_eq_int(0);
        }
    }
}

/**
 * GFFL: level-crossing balance between n_i and n_i + 1 (AMPL THM3).
 *
 * The arrival rate to station i while it holds ni jobs equals the completion
 * rate at i while it holds ni + 1. With it absent, a station idle with
 * probability one is feasible and the utilization lower bound collapses to
 * zero.
 */
template <class T>
void lr_gffl(const MapqnParams<T>& p, const LrIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int ni = 0; ni <= p.N - 1; ++ni) {
            for (int j = 0; j < p.M; ++j) {
                if (j == i) continue;
                for (int k = 0; k < p.K[j]; ++k)
                    for (int h = 0; h < p.K[j]; ++h)
                        for (int u = 0; u < p.K[i]; ++u)
                            m.row_add(x.p1(j, k, i, ni, u), qr_rate(p, j, i, k, h));
                for (int k = 0; k < p.K[i]; ++k)
                    for (int h = 0; h < p.K[i]; ++h)
                        m.row_add(x.p1(i, k, i, ni + 1, k), T(-qr_rate(p, i, j, k, h)));
            }
            m.emit_eq_int(0);
        }
    }
}

/**
 * QBAL: throughput balance at each station.
 *
 * The departure rate from i equals the arrival rate to i, with the arrival side
 * split over whether station i is busy (the p1(i,u,j,nj,k) term, nj >= 1) or
 * empty (the p1(j,k,i,0,u) term).
 */
template <class T>
void lr_qbal(const MapqnParams<T>& p, const LrIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int j = 0; j < p.M; ++j) {
            if (j == i) continue;
            for (int k = 0; k < p.K[i]; ++k)
                for (int h = 0; h < p.K[i]; ++h)
                    m.row_add(x.U(i, k), qr_rate(p, i, j, k, h));
            for (int k = 0; k < p.K[j]; ++k) {
                for (int h = 0; h < p.K[j]; ++h) {
                    const T w = qr_rate(p, j, i, k, h);
                    for (int u = 0; u < p.K[i]; ++u)
                        for (int nj = 1; nj <= p.N; ++nj)
                            m.row_add(x.p1(i, u, j, nj, k), T(-w));
                    for (int u = 0; u < p.K[i]; ++u) m.row_add(x.p1(j, k, i, 0, u), T(-w));
                }
            }
        }
        m.emit_eq_int(0);
    }
}

}  // namespace detail

/**
 * Bound U(objective_queue, objective_phase) over the linear-reduction polytope.
 *
 * @param p               network parameters; queues and phases are 0-based.
 *                        alpha is IGNORED: this model has no load dependence.
 * @param objective_queue queue index, 0..M-1
 * @param objective_phase phase index, 0..K(objective_queue)-1
 * @param sense           Max for an upper bound, Min for a lower bound
 */
template <class T>
MapqnBndLrResult<T> mapqn_bnd_lr(const MapqnParams<T>& p, int objective_queue,
                                 int objective_phase, MapqnSense sense = MapqnSense::Max) {
    p.validate();
    detail::p1_check_objective(p, objective_queue, objective_phase);

    const LrIndex x(p.M, p.N, p.K, false);
    lp::LpModel<T> m(x.num_vars());
    detail::p1_bounds(p, x, m);

    // Families, in the order the reference emits them. Named one per function
    // so the inventory is diffable against mapqn_bnd_lr.m; see the header note
    // on why an omitted family reads as a loose bound rather than an error.
    detail::p1_zer1(p, x, m);
    detail::p1_zer2(p, x, m);
    detail::p1_zer3(p, x, m);
    detail::p1_zer4(p, x, m);
    detail::p1_cequ(p, x, m);
    detail::p1_one1(p, x, m);
    detail::p1_utlb(p, x, m);
    detail::p1_utlc(p, x, m);
    detail::p1_qlen(p, x, m);
    detail::p1_clen(p, x, m);
    detail::p1_one(p, x, m);
    detail::p1_popc(p, x, m);
    detail::lr_mpcb(p, x, m);
    detail::p1_srvb(p, x, m);
    detail::lr_ujnt(p, x, m);
    detail::lr_gffl0(p, x, m);
    detail::lr_gffl(p, x, m);
    detail::lr_qbal(p, x, m);
    detail::p1_uub1(p, x, m);
    detail::p1_qub1(p, x, m);

    m.set_cost(x.U(objective_queue, objective_phase), num_traits<T>::from_int(1));
    m.set_maximize(sense == MapqnSense::Max);

    const lp::LpSolution<T> sol = lp::simplex_solve(m);

    MapqnBndLrResult<T> out;
    out.status = lp::lp_status_name(sol.status);
    out.ok = sol.ok();
    out.objective = sol.objective;
    out.x = sol.x;
    out.num_vars = m.num_vars();
    out.num_rows = m.num_rows();
    out.iterations = sol.iterations;
    if (!out.ok) return out;

    std::size_t maxK = 0;
    for (int i = 0; i < p.M; ++i)
        if (static_cast<std::size_t>(p.K[i]) > maxK) maxK = static_cast<std::size_t>(p.K[i]);
    out.U = Matrix<T>(static_cast<std::size_t>(p.M), maxK);
    out.IT = Matrix<T>(static_cast<std::size_t>(p.M), maxK);
    out.Q = Matrix<T>(static_cast<std::size_t>(p.M), maxK);
    for (int i = 0; i < p.M; ++i) {
        for (int k = 0; k < p.K[i]; ++k) {
            const std::size_t ii = static_cast<std::size_t>(i), kk = static_cast<std::size_t>(k);
            out.U(ii, kk) = sol.x[x.U(i, k)];
            out.IT(ii, kk) = sol.x[x.IT(i, k)];
            out.Q(ii, kk) = sol.x[x.Q(i, k)];
        }
    }
    return out;
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_BND_LR_H
