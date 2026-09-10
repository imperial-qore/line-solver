/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_BND_QR_H
#define LINE_API_MAPQN_MAPQN_BND_QR_H

/**
 * General quadratic-reduction (QR) bound on the utilization of one queue-phase
 * of a closed MAP queueing network.
 *
 * Templated port of matlab/lib/qrf/mapqn_bnd_qr.m (ground truth), cross-checked
 * against python/line_solver/api/mapqn/bnd_qr.py. It is NOT ported from
 * jar/src/main/java/jline/api/mapqn/Mapqn_bnd_qr.java, which hardcodes
 * GoalType.MAXIMIZE and takes no sense argument, so it cannot express the lower
 * bound at all (see _kb/03-api-layer.md).
 *
 * THIS IS NOT mapqn_bnd_qr_ld AT alpha == 1. The load-dependent model relaxes
 * onto the pairwise law p2 alone and imposes 15 families over it. This model
 * carries the singly-indexed p1 and p1c alongside p2, plus the aggregate U, IT,
 * Q and C, and ties the three levels together with the projection families
 * PI21/PI22/PI23. The extra level buys the aggregate inequalities UUB1, QUB1,
 * CUB1, CUB2 and THM4, which have no expression in a pure p2 model. The two
 * polytopes are different relaxations of the same chain and neither contains
 * the other; keep both.
 *
 * WHAT MAKES IT A BOUND. Every constraint below is satisfied by the exact
 * stationary distribution, and none of them pins it down, so the feasible set
 * is a polytope CONTAINING the exact solution. Optimizing U(i,k) over it
 * therefore brackets the true utilization from whichever side is asked for.
 * The bound is a relaxation, not an approximation: it is valid, not merely
 * close.
 *
 * THE FAILURE MODE OF THIS FAMILY IS A VACUOUS BOUND, NOT A CRASH. Omit THM30
 * and THM3 and no constraint mentions mu at all, so U = 0 and U = 1 are both
 * feasible and the routine returns the [0,1] box while every subscript in the
 * file stays internally consistent. That is exactly how the 2026-07-20 defect
 * in mapqn_bnd_qr.m survived (it carried only SRVB and returned [0,1/3] where
 * the true range is [0.1970, 0.2182]). Each family is therefore emitted by a
 * function named after it, so the inventory is diffable against the reference
 * by name. Do not inline them.
 *
 * INERT VARIABLES OMITTED. The reference registers UP(j,k,i,h), QP(j,k,i,h) and
 * I(j,k,i) and gives them upper bounds, then references them in no constraint
 * and in no objective (verified: UPidx, QPidx and Iidx appear in
 * mapqn_bnd_qr.m only at their own registration and bounding). A variable with
 * no row and no cost cannot move the optimum, so they are not allocated here.
 * This drops 2 (sum_i K(i))^2 + M sum_i K(i) columns and changes no bound.
 *
 * ARITHMETIC. Assembly is +, -, * on the model data and lp::simplex_solve uses
 * Bland's rule with no tolerance, so at T = line::Rational the returned bound
 * is the EXACT optimum of the exact polytope. The MATLAB reference reaches it
 * with linprog's 'interior-point' (its 'interior-point-legacy' declares this
 * system infeasible once the balance families are present) and lands a few
 * digits short, so a deviation against MATLAB is expected to be MATLAB's
 * convergence gap, not this port's error.
 *
 * COST. (N+1) sum_i K(i) is the pairwise half-index B; the model has
 * B^2 + 2 B sum_i K(i) + (M+3) sum_i K(i) columns, so it grows as the fourth
 * power of the population. The tableau is dense, which confines the port to
 * small and medium instances.
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

/** Result of a general QR bound solve. */
template <class T>
struct MapqnBndQrResult {
    bool ok = false;      ///< the LP reached an optimal vertex
    std::string status;   ///< textual LP status
    T objective = T();    ///< the bound on U(objective_queue, objective_phase)
    Matrix<T> U;          ///< M x max(K) utilizations, 0 beyond K(i)
    Matrix<T> IT;         ///< M x max(K) idle times, 0 beyond K(i)
    Matrix<T> Q;          ///< M x max(K) mean queue lengths, 0 beyond K(i)
    std::vector<T> x;     ///< full solution vector, indexed by QrIndex
    std::size_t num_vars = 0;
    std::size_t num_rows = 0;
    std::size_t iterations = 0;
};

/**
 * Variable layout of the general QR model: the shared p1-level blocks plus the
 * joint p2 block. Defined in mapqn_p1_common.h, which the linear reduction
 * shares; the alias keeps the name this model is documented and tested under.
 */
using QrIndex = MapqnP1Index;

namespace detail {

/** ZER5: p2(j,nj,k,j,nj,h) = 0 for h != k (AMPL ZERO1). */
template <class T>
void qr_zer5(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int nj = 0; nj <= p.N; ++nj)
            for (int k = 0; k < p.K[j]; ++k)
                for (int h = 0; h < p.K[j]; ++h)
                    if (h != k) m.fix(x.p2(j, nj, k, j, nj, h), T());
}

/** ZER6: p2(j,nj,k,j,ni,h) = 0 for ni != nj (AMPL ZERO2). */
template <class T>
void qr_zer6(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int nj = 0; nj <= p.N; ++nj)
            for (int k = 0; k < p.K[j]; ++k)
                for (int ni = 0; ni <= p.N; ++ni)
                    if (ni != nj)
                        for (int h = 0; h < p.K[j]; ++h) m.fix(x.p2(j, nj, k, j, ni, h), T());
}

/** ZER7: p2(j,nj,k,i,ni,h) = 0 for i != j and nj + ni > N (AMPL ZERO3). */
template <class T>
void qr_zer7(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j)
        for (int nj = 0; nj <= p.N; ++nj)
            for (int k = 0; k < p.K[j]; ++k)
                for (int i = 0; i < p.M; ++i)
                    if (i != j)
                        for (int ni = 0; ni <= p.N; ++ni)
                            if (nj + ni > p.N)
                                for (int h = 0; h < p.K[i]; ++h)
                                    m.fix(x.p2(j, nj, k, i, ni, h), T());
}

/** PCL2: the second moment of the population, sum ni*nj*p2 = N^2. */
template <class T>
void qr_pcl2(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i)
        for (int j = 0; j < p.M; ++j)
            for (int ni = 1; ni <= p.N; ++ni)
                for (int nj = 1; nj <= p.N; ++nj)
                    for (int h = 0; h < p.K[i]; ++h)
                        for (int k = 0; k < p.K[j]; ++k)
                            m.row_add_int(x.p2(i, ni, h, j, nj, k), ni * nj);
    m.emit_eq_int(p.N * p.N);
}

/** PI21: p1 is the projection of p2 over the busy populations nj >= 1. */
template <class T>
void qr_pi21(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int i = 0; i < p.M; ++i) {
                for (int ni = 0; ni <= p.N; ++ni) {
                    for (int h = 0; h < p.K[i]; ++h) {
                        m.row_add_int(x.p1(j, k, i, ni, h), 1);
                        for (int nj = 1; nj <= p.N; ++nj)
                            m.row_add_int(x.p2(j, nj, k, i, ni, h), -1);
                        m.emit_eq_int(0);
                    }
                }
            }
        }
    }
}

/** PI22: p1c is the nj = 0 slice of p2. */
template <class T>
void qr_pi22(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int i = 0; i < p.M; ++i) {
                for (int ni = 0; ni <= p.N; ++ni) {
                    for (int h = 0; h < p.K[i]; ++h) {
                        m.row_add_int(x.p1c(j, k, i, ni, h), 1);
                        m.row_add_int(x.p2(j, 0, k, i, ni, h), -1);
                        m.emit_eq_int(0);
                    }
                }
            }
        }
    }
}

/**
 * PI23: p2 is symmetric under swapping its two halves.
 *
 * Emitted only for lexicographically ordered pairs, as the reference does; the
 * reverse pair is the same row negated and the diagonal is the trivial 0 = 0.
 */
template <class T>
void qr_pi23(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int nj = 0; nj <= p.N; ++nj) {
            for (int k = 0; k < p.K[j]; ++k) {
                for (int i = 0; i < p.M; ++i) {
                    for (int ni = 0; ni <= p.N; ++ni) {
                        for (int h = 0; h < p.K[i]; ++h) {
                            const bool ordered = (j < i) || (j == i && nj < ni) ||
                                                 (j == i && nj == ni && k < h);
                            if (!ordered) continue;
                            const std::size_t a = x.p2(i, ni, h, j, nj, k);
                            const std::size_t b = x.p2(j, nj, k, i, ni, h);
                            if (a == b) continue;
                            m.row_add_int(a, 1);
                            m.row_add_int(b, -1);
                            m.emit_eq_int(0);
                        }
                    }
                }
            }
        }
    }
}

/**
 * MARG: the pairwise law agrees with its own marginal at every population
 * (AMPL MARGINALS). ONE1 imposes only the aggregate over nj, so this is a
 * separate family.
 */
template <class T>
void qr_marg(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int nj = 0; nj <= p.N; ++nj) {
                for (int i = 0; i < p.M; ++i) {
                    if (i == j) continue;
                    m.row_add_int(x.p2(j, nj, k, j, nj, k), 1);
                    for (int ni = 0; ni <= p.N - nj; ++ni)
                        for (int h = 0; h < p.K[i]; ++h)
                            m.row_add_int(x.p2(j, nj, k, i, ni, h), -1);
                    m.emit_eq_int(0);
                }
            }
        }
    }
}

/**
 * THM2: the queue-length theorem conditioned on (j,nj,k). Summed over nj >= 1
 * it recovers the aggregate sum_i C(j,k,i) = N U(j,k).
 */
template <class T>
void qr_thm2(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int nj = 0; nj <= p.N; ++nj) {
                for (int i = 0; i < p.M; ++i)
                    for (int ni = 1; ni <= p.N; ++ni)
                        for (int h = 0; h < p.K[i]; ++h)
                            m.row_add_int(x.p2(j, nj, k, i, ni, h), ni);
                m.row_add_int(x.p2(j, nj, k, j, nj, k), -p.N);
                m.emit_eq_int(0);
            }
        }
    }
}

/**
 * THM30: level-crossing balance at an empty station, per arrival phase (AMPL
 * THM30). The rate into {n_i = 0, phase_i = u} from a busy neighbour equals the
 * rate out of {n_i = 1} through a completion at i.
 */
template <class T>
void qr_thm30(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int u = 0; u < p.K[i]; ++u) {
            for (int j = 0; j < p.M; ++j) {
                if (j == i) continue;
                for (int nj = 1; nj <= p.N; ++nj)
                    for (int k = 0; k < p.K[j]; ++k)
                        for (int h = 0; h < p.K[j]; ++h)
                            m.row_add(x.p2(j, nj, k, i, 0, u), qr_rate(p, j, i, k, h));
                for (int nj = 0; nj <= p.N; ++nj)
                    for (int k = 0; k < p.K[i]; ++k)
                        for (int h = 0; h < p.K[j]; ++h)
                            m.row_add(x.p2(j, nj, h, i, 1, k), T(-qr_rate(p, i, j, k, u)));
            }
            m.emit_eq_int(0);
        }
    }
}

/**
 * THM3: level-crossing balance between n_i and n_i + 1 (AMPL THM3).
 *
 * This family and THM30 are the only ones that mention mu. Drop them and no
 * constraint distinguishes a fast station from a slow one, so U = 0 stays
 * feasible and the bound collapses to the [0,1] box.
 */
template <class T>
void qr_thm3(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int i = 0; i < p.M; ++i) {
        for (int ni = 0; ni <= p.N - 1; ++ni) {
            for (int j = 0; j < p.M; ++j) {
                if (j == i) continue;
                for (int nj = 1; nj <= p.N; ++nj)
                    for (int k = 0; k < p.K[j]; ++k)
                        for (int h = 0; h < p.K[j]; ++h)
                            for (int u = 0; u < p.K[i]; ++u)
                                m.row_add(x.p2(j, nj, k, i, ni, u), qr_rate(p, j, i, k, h));
                for (int nj = 0; nj <= p.N; ++nj)
                    for (int k = 0; k < p.K[i]; ++k)
                        for (int u = 0; u < p.K[j]; ++u)
                            for (int h = 0; h < p.K[i]; ++h)
                                m.row_add(x.p2(j, nj, u, i, ni + 1, k),
                                          T(-qr_rate(p, i, j, k, h)));
            }
            m.emit_eq_int(0);
        }
    }
}

/** CUB1: C(j,k,i) <= sum_h Q(i,h). A conditional length cannot exceed the mean. */
template <class T>
void qr_cub1(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int i = 0; i < p.M; ++i) {
                m.row_add_int(x.C(j, k, i), 1);
                for (int h = 0; h < p.K[i]; ++h) m.row_add_int(x.Q(i, h), -1);
                m.emit_le_int(0);
            }
        }
    }
}

/** CUB2: C(j,k,i) <= N U(j,k). The conditioning event has probability U. */
template <class T>
void qr_cub2(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int i = 0; i < p.M; ++i) {
                m.row_add_int(x.C(j, k, i), 1);
                m.row_add_int(x.U(j, k), -p.N);
                m.emit_le_int(0);
            }
        }
    }
}

/** THM4: the QMIN inequality, N P(j busy in k, i nonempty) <= sum_t C(j,k,t). */
template <class T>
void qr_thm4(const MapqnParams<T>& p, const QrIndex& x, lp::LpModel<T>& m) {
    for (int j = 0; j < p.M; ++j) {
        for (int k = 0; k < p.K[j]; ++k) {
            for (int i = 0; i < p.M; ++i) {
                for (int h = 0; h < p.K[i]; ++h)
                    for (int nj = 0; nj <= p.N; ++nj)
                        for (int ni = 1; ni <= p.N; ++ni)
                            m.row_add_int(x.p2(j, nj, k, i, ni, h), p.N);
                for (int t = 0; t < p.M; ++t)
                    for (int h = 0; h < p.K[t]; ++h)
                        for (int nj = 0; nj <= p.N; ++nj)
                            for (int nt = 0; nt <= p.N; ++nt)
                                m.row_add_int(x.p2(j, nj, k, t, nt, h), -nt);
                m.emit_le_int(0);
            }
        }
    }
}

}  // namespace detail

/**
 * Bound U(objective_queue, objective_phase) over the general QR polytope.
 *
 * @param p               network parameters; queues and phases are 0-based.
 *                        alpha is IGNORED: this model has no load dependence,
 *                        use mapqn_bnd_qr_ld for that.
 * @param objective_queue queue index, 0..M-1
 * @param objective_phase phase index, 0..K(objective_queue)-1
 * @param sense           Max for an upper bound, Min for a lower bound
 */
template <class T>
MapqnBndQrResult<T> mapqn_bnd_qr(const MapqnParams<T>& p, int objective_queue,
                                 int objective_phase, MapqnSense sense = MapqnSense::Max) {
    p.validate();
    detail::p1_check_objective(p, objective_queue, objective_phase);

    const QrIndex x(p.M, p.N, p.K, true);
    lp::LpModel<T> m(x.num_vars());
    detail::p1_bounds(p, x, m);

    // Families, in the order the reference emits them. Named one per function
    // so the inventory is diffable against mapqn_bnd_qr.m; see the header note
    // on why an omitted family reads as a loose bound rather than an error.
    detail::p1_zer1(p, x, m);
    detail::p1_zer2(p, x, m);
    detail::p1_zer3(p, x, m);
    detail::p1_zer4(p, x, m);
    detail::qr_zer5(p, x, m);
    detail::qr_zer6(p, x, m);
    detail::qr_zer7(p, x, m);
    detail::p1_cequ(p, x, m);
    detail::p1_one1(p, x, m);
    detail::p1_utlb(p, x, m);
    detail::p1_utlc(p, x, m);
    detail::p1_qlen(p, x, m);
    detail::p1_srvb(p, x, m);
    detail::p1_popc(p, x, m);
    detail::p1_one(p, x, m);
    detail::qr_pcl2(p, x, m);
    detail::qr_pi21(p, x, m);
    detail::qr_pi22(p, x, m);
    detail::qr_pi23(p, x, m);
    detail::p1_clen(p, x, m);
    detail::qr_marg(p, x, m);
    detail::qr_thm2(p, x, m);
    detail::qr_thm30(p, x, m);
    detail::qr_thm3(p, x, m);
    detail::p1_uub1(p, x, m);
    detail::p1_qub1(p, x, m);
    detail::qr_cub1(p, x, m);
    detail::qr_cub2(p, x, m);
    detail::qr_thm4(p, x, m);

    m.set_cost(x.U(objective_queue, objective_phase), num_traits<T>::from_int(1));
    m.set_maximize(sense == MapqnSense::Max);

    const lp::LpSolution<T> sol = lp::simplex_solve(m);

    MapqnBndQrResult<T> out;
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

#endif  // LINE_API_MAPQN_MAPQN_BND_QR_H
