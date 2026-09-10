/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_QRF_COMMON_H
#define LINE_API_MAPQN_MAPQN_QRF_COMMON_H

/**
 * Shared machinery of the QRF nonlinear bounds (`qrf_noblo_*`, `qrf_bas_*`).
 *
 * Port of python/line_solver/api/mapqn/qrf_noblo_common.py (the JAR twin is
 * `Mapqn_qrf_noblo_*` plus `Mapqn_nlp_solver`). `api/mapqn` has NO MATLAB
 * implementation, so Python and the JAR are the references here.
 *
 * WHAT THE FAMILY IS. The linear QRF bounds (`mapqn_qr_bounds_*`) MAXIMIZE one
 * station's utilization over a polytope of pairwise queue-phase probabilities.
 * These entry points optimize a different objective over the SAME feasible set:
 * mutual information (MMI) or negative entropy (MEM).
 *
 * NEITHER OBJECTIVE IS A CONVEX PROGRAM AS THE REFERENCE STATES IT, and a
 * caller has to know which one it asked for:
 *
 * - MEM minimizes -sum p log p over the diagonal entries. Despite the name,
 *   that MINIMIZES the entropy, and -p log p is CONCAVE, so the minimum sits at
 *   a vertex and every method reports a stationary point fixed by its start.
 *   The three codebases agree because they start from the same place: the
 *   MINIMUM-NORM feasible point, not an arbitrary phase-1 vertex. See
 *   `qrf_feasible_start_lp`, which is where that is arranged; measured on a
 *   two-station cycle at N = 2 with both one and two phases.
 * - MMI minimizes sum_{i != j} p_ij (log p_ij - log p_ii - log p_jj). The
 *   `-p_ij log p_ii` terms are NOT convex, so the problem is nonconvex and both
 *   codebases report a LOCAL optimum whose identity is fixed by the start
 *   point, i.e. by whichever vertex the phase-1 LP happens to return. Measured:
 *   on a one-phase pair the two agree exactly (0.857143 / 0.428571, itself the
 *   exact product-form answer); on a two-phase / one-phase pair they land on
 *   different vertices, 0.6931 against 0.6438 in objective, and NEITHER is
 *   wrong -- the reference's own solve does not move from its start either.
 *   This is verifiable rather than suspected: the reference's start point is
 *   FEASIBLE under this port's constraints to 8.9e-16, and from this port's
 *   start it is an ASCENT direction (a directional derivative of +0.184), which
 *   only a nonconvex objective permits.
 *
 * So MMI parity is start-point parity, and reproducing it would mean
 * reproducing HiGHS's pivoting rather than any property of the model.
 *
 * THE THREE SUBSTITUTIONS FOR SCIPY, and why each is exact rather than
 * approximate:
 *
 * 1. `scipy.linalg.qr(pivoting=True)` behind `independent_rows` is replaced by
 *    a pivoted modified Gram-Schmidt over the same columns. Column-pivoted QR
 *    selects, at each step, the remaining column of largest residual norm, and
 *    that is precisely what the sweep below does; the retained index SET is the
 *    same, and the set is all the caller uses. (The R factor itself is never
 *    read.)
 * 2. `scipy.optimize.linprog` behind `feasible_start` is replaced by
 *    `lp::simplex_solve` on the same rows with the same [0,1] box, followed by
 *    the minimum-norm refinement `qrf_noblo_start.m` performs with QUADPROG.
 *    The LP is exact in either; the refinement is what keeps the start off a
 *    vertex, and on these objectives the start decides the answer.
 * 3. `scipy.optimize.minimize(method='SLSQP')` is replaced by FRANK-WOLFE over
 *    the same polytope. See `solve_qrf_nlp` for why this is the right answer
 *    and not a compromise: the available augmented-Lagrangian path has a
 *    derivative-free inner solve, and started at a VERTEX -- which is what the
 *    phase-1 LP returns -- it does not move at all. Conditional gradient uses
 *    the LP that is already here and returns an optimality CERTIFICATE.
 *
 * THE REDUCTION IS NOT AN OPTIMIZATION, IT IS WHAT MAKES THE PROBLEM SOLVABLE.
 * The raw decision vector has M^2 (N+1)^2 Kmax^2 MR + M Kmax entries and the
 * equality block pins nearly all of them; substituting x = x0 + Z t for an
 * orthonormal basis Z of null(Aeq) leaves a few free directions. The reference
 * measured SLSQP at 14 iterations and under a second on the reduced problem
 * against a failure to move on the raw one.
 *
 * THE GUARD THAT MUST NOT BE DROPPED. `reduce_equalities` refuses an
 * INCONSISTENT system rather than dropping the offending rows -- an
 * inconsistent polytope is a modelling error, and silently discarding it
 * returns numbers for a model nobody wrote. The reference's second guard, a
 * one-shot probe asking whether a feasible descent direction exists at the
 * start (its SLSQP can silently return the phase-1 vertex), is SUBSUMED here:
 * `solve_qrf_nlp` asks that same question at EVERY iterate, as its Frank-Wolfe
 * gap, and terminates on the answer.
 *
 * ARITHMETIC: transcendental (p log p, and the orthogonalization).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/lp_highs.h"
#include "line/util/simplex.h"

namespace line {
namespace mapqn {

/**
 * log() in the working arithmetic.
 *
 * `num_traits` exposes the logarithm only as a double (`log_as_double`), which
 * is enough here: every consumer of this file is gated on
 * `has_transcendental`, so T is a floating type and the round trip loses
 * nothing it did not already lose.
 */
template <class T>
T qrf_log(const T& v) {
    return num_traits<T>::from_double(num_traits<T>::log_as_double(v));
}

/** The reference's LOGTOL: the shift that keeps log() off zero. */
template <class T>
T qrf_logtol() {
    return num_traits<T>::from_double(1e-6);
}

/** The unflattened decision vector: the pair tensor and the effective rates. */
template <class T>
struct QrfVars {
    /** p2[j][nj][k][i][ni][h][m], flattened row-major over the seven indices. */
    std::vector<T> p2;
    std::vector<T> e;  ///< e[i*Kmax + k]
    std::size_t M = 0, N = 0, Kmax = 0, MR = 0;

    std::size_t p2_at(std::size_t j, std::size_t nj, std::size_t k, std::size_t i,
                      std::size_t ni, std::size_t h, std::size_t m) const {
        return (((((j * (N + 1) + nj) * Kmax + k) * M + i) * (N + 1) + ni) * Kmax + h) * MR + m;
    }
    const T& p(std::size_t j, std::size_t nj, std::size_t k, std::size_t i, std::size_t ni,
               std::size_t h, std::size_t m) const {
        return p2[p2_at(j, nj, k, i, ni, h, m)];
    }
    T& p(std::size_t j, std::size_t nj, std::size_t k, std::size_t i, std::size_t ni,
         std::size_t h, std::size_t m) {
        return p2[p2_at(j, nj, k, i, ni, h, m)];
    }
};

/**
 * Number of decision variables the layout actually USES.
 *
 * THIS IS DELIBERATELY SMALLER THAN THE REFERENCE'S `compute_num_vars`, which
 * returns the full `M (N+1) Kmax M (N+1) Kmax MR + M Kmax` tensor. The fill
 * loops of `sub_qrfvar` run over K[j] and K[i], not over Kmax, so with
 * heterogeneous phase counts the tail of that vector is never written, never
 * read by a constraint, and never read by an objective -- and, being unread, it
 * has an all-zero column in Aeq, so every one of those coordinates lands in
 * null(Aeq) as a FLAT direction of the reduced problem.
 *
 * The reference's gradient-based SLSQP shrugs that off (a flat direction has
 * zero gradient). The augmented-Lagrangian inner solve here is Nelder-Mead,
 * which degrades sharply with dimension, and the padding is not a rounding
 * detail: on a two-phase / one-phase pair at N = 2 it is 64 flat directions
 * against 20 real ones, and the solve returns a utilization of exactly 1
 * against the reference's 0.8. Dropping the dead coordinates is EXACT, not an
 * approximation -- nothing reads them -- and it is what makes the port agree.
 *
 * The live coordinates are a PREFIX of the reference's vector, since the fill
 * order is consecutive, so the two layouts agree wherever both are defined.
 */
inline std::size_t qrf_num_vars(std::size_t M, std::size_t N, const std::vector<int>& K,
                                std::size_t MR) {
    std::size_t outer = 0, phases = 0;
    for (std::size_t i = 0; i < M; ++i) {
        outer += (N + 1) * static_cast<std::size_t>(K[i]);
        phases += static_cast<std::size_t>(K[i]);
    }
    return outer * outer * MR + phases;
}

/**
 * Flat position of every p2 entry, in the FILL ORDER of sub_qrfvar.
 *
 * The layout is not a plain strided tensor: the phase loops run over K[j] and
 * K[i], not over Kmax, so a station with fewer phases leaves GAPS. An entry
 * that carries no variable keeps -1, and every gradient scatter must skip it.
 */
inline std::vector<long> qrf_index_map(std::size_t M, std::size_t N, const std::vector<int>& K,
                                       std::size_t MR) {
    const std::size_t Kmax = static_cast<std::size_t>(*std::max_element(K.begin(), K.end()));
    QrfVars<double> shape;
    shape.M = M;
    shape.N = N;
    shape.Kmax = Kmax;
    shape.MR = MR;
    std::vector<long> idx(M * (N + 1) * Kmax * M * (N + 1) * Kmax * MR, -1);
    long ctr = 0;
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t nj = 0; nj <= N; ++nj)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[j]); ++k)
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t ni = 0; ni <= N; ++ni)
                        for (std::size_t h = 0; h < static_cast<std::size_t>(K[i]); ++h)
                            for (std::size_t m = 0; m < MR; ++m)
                                idx[shape.p2_at(j, nj, k, i, ni, h, m)] = ctr++;
    return idx;
}

/** Unflatten x into the pair tensor and the effective rates. */
template <class T>
QrfVars<T> sub_qrfvar(const std::vector<T>& x, std::size_t M, std::size_t N,
                      const std::vector<int>& K, std::size_t MR) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t Kmax = static_cast<std::size_t>(*std::max_element(K.begin(), K.end()));
    QrfVars<T> out;
    out.M = M;
    out.N = N;
    out.Kmax = Kmax;
    out.MR = MR;
    out.p2.assign(M * (N + 1) * Kmax * M * (N + 1) * Kmax * MR, zero);
    out.e.assign(M * Kmax, zero);
    std::size_t ctr = 0;
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t nj = 0; nj <= N; ++nj)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[j]); ++k)
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t ni = 0; ni <= N; ++ni)
                        for (std::size_t h = 0; h < static_cast<std::size_t>(K[i]); ++h)
                            for (std::size_t m = 0; m < MR; ++m)
                                out.p(j, nj, k, i, ni, h, m) = x[ctr++];
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k)
            out.e[i * Kmax + k] = x[ctr++];
    // x may be LONGER than the fill consumes, and that is the reference's own
    // layout rather than an error: `compute_num_vars` sizes the vector on the
    // full Kmax tensor while the fill loops run over K[j] and K[i], so a
    // station with fewer phases leaves DEAD entries that no constraint and no
    // objective ever reads. Too SHORT is a real defect and is refused.
    if (ctr > x.size())
        throw InputError("sub_qrfvar: the decision vector is shorter than the layout requires");
    return out;
}

/**
 * Mutual-information objective.
 *
 * sum over i != j of p_ij (log p_ij - log p_ii - log p_jj), each log taken on
 * the LOGTOL-shifted value so the sweep can visit the boundary.
 */
template <class T>
T mmi_objective(const std::vector<T>& x, std::size_t M, std::size_t N,
                const std::vector<int>& K, const std::vector<int>& F, std::size_t MR) {
    const T tol = qrf_logtol<T>();
    const QrfVars<T> v = sub_qrfvar(x, M, N, K, MR);
    T fobj = num_traits<T>::from_int(0);
    for (std::size_t m = 0; m < MR; ++m)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t ki = 0; ki < static_cast<std::size_t>(K[i]); ++ki)
                for (std::size_t j = 0; j < M; ++j) {
                    if (i == j) continue;
                    for (std::size_t kj = 0; kj < static_cast<std::size_t>(K[j]); ++kj)
                        for (std::size_t ni = 0; ni <= static_cast<std::size_t>(F[i]); ++ni)
                            for (std::size_t nj = 0; nj <= static_cast<std::size_t>(F[j]); ++nj) {
                                const T pij = v.p(i, ni, ki, j, nj, kj, m);
                                const T pii = v.p(i, ni, ki, i, ni, ki, m);
                                const T pjj = v.p(j, nj, kj, j, nj, kj, m);
                                fobj += pij * (qrf_log<T>(tol + pij) -
                                               qrf_log<T>(tol + pii) -
                                               qrf_log<T>(tol + pjj));
                            }
                }
    return fobj;
}

/**
 * Maximum-entropy objective, returned as the NEGATIVE entropy +sum p log p
 * over the diagonal entries, because the solver MINIMIZES and the AMPL model
 * states this objective as `maximize H`. Returning +H (as every port did until
 * 2026-08-29) selects the minimum-entropy face of the polytope instead, under
 * a method documented as maximum-entropy.
 */
template <class T>
T mem_objective(const std::vector<T>& x, std::size_t M, std::size_t N,
                const std::vector<int>& K, const std::vector<int>& F, std::size_t MR) {
    const T tol = qrf_logtol<T>();
    const QrfVars<T> v = sub_qrfvar(x, M, N, K, MR);
    T fobj = num_traits<T>::from_int(0);
    for (std::size_t m = 0; m < MR; ++m)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k)
                for (std::size_t ni = 1; ni <= static_cast<std::size_t>(F[i]); ++ni) {
                    const T pv = v.p(i, ni, k, i, ni, k, m);
                    fobj += pv * qrf_log<T>(tol + pv);
                }
    return fobj;
}

/**
 * Gradient of mmi_objective.
 *
 * For t = p_ij (log p_ij - log p_ii - log p_jj) the three partials are
 * dt/dp_ij = log p_ij - log p_ii - log p_jj + p_ij/p_ij', dt/dp_ii =
 * -p_ij/p_ii' and dt/dp_jj = -p_ij/p_jj', a primed denominator standing for the
 * shifted value. i != j throughout, so no term aliases its own partials.
 */
template <class T>
std::vector<T> mmi_gradient(const std::vector<T>& x, std::size_t M, std::size_t N,
                            const std::vector<int>& K, const std::vector<int>& F, std::size_t MR,
                            const std::vector<long>& idx) {
    const T tol = qrf_logtol<T>();
    const QrfVars<T> v = sub_qrfvar(x, M, N, K, MR);
    std::vector<T> g(x.size(), num_traits<T>::from_int(0));
    for (std::size_t m = 0; m < MR; ++m)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t ki = 0; ki < static_cast<std::size_t>(K[i]); ++ki)
                for (std::size_t j = 0; j < M; ++j) {
                    if (i == j) continue;
                    for (std::size_t kj = 0; kj < static_cast<std::size_t>(K[j]); ++kj)
                        for (std::size_t ni = 0; ni <= static_cast<std::size_t>(F[i]); ++ni) {
                            const T pii = v.p(i, ni, ki, i, ni, ki, m);
                            const long iii = idx[v.p2_at(i, ni, ki, i, ni, ki, m)];
                            for (std::size_t nj = 0; nj <= static_cast<std::size_t>(F[j]); ++nj) {
                                const T pij = v.p(i, ni, ki, j, nj, kj, m);
                                const T pjj = v.p(j, nj, kj, j, nj, kj, m);
                                const long iij = idx[v.p2_at(i, ni, ki, j, nj, kj, m)];
                                const long ijj = idx[v.p2_at(j, nj, kj, j, nj, kj, m)];
                                if (iij >= 0)
                                    g[iij] += qrf_log<T>(tol + pij) -
                                              qrf_log<T>(tol + pii) -
                                              qrf_log<T>(tol + pjj) + pij / (tol + pij);
                                if (iii >= 0) g[iii] -= pij / (tol + pii);
                                if (ijj >= 0) g[ijj] -= pij / (tol + pjj);
                            }
                        }
                }
    return g;
}

/** Gradient of mem_objective: d/dp of p log(p') is log p' + p/p'. */
template <class T>
std::vector<T> mem_gradient(const std::vector<T>& x, std::size_t M, std::size_t N,
                            const std::vector<int>& K, const std::vector<int>& F, std::size_t MR,
                            const std::vector<long>& idx) {
    const T tol = qrf_logtol<T>();
    const QrfVars<T> v = sub_qrfvar(x, M, N, K, MR);
    std::vector<T> g(x.size(), num_traits<T>::from_int(0));
    for (std::size_t m = 0; m < MR; ++m)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k)
                for (std::size_t ni = 1; ni <= static_cast<std::size_t>(F[i]); ++ni) {
                    const T pv = v.p(i, ni, k, i, ni, k, m);
                    const long ip = idx[v.p2_at(i, ni, k, i, ni, k, m)];
                    if (ip >= 0)
                        g[ip] += qrf_log<T>(tol + pv) + pv / (tol + pv);
                }
    return g;
}

/**
 * Tree-reweighted (Bethe) free entropy at the uniform spanning-tree weight,
 * the objective of `qrf.bethe`.
 *
 * With lambda = 1/M this is lambda*sum_{i!=j} I(n_i;n_j) - sum_i H(n_i), the
 * NEGATIVE of a tree-reweighted entropy with uniform edge weight
 * rho_ij = 2*lambda on the complete station graph. H_rho is a convex
 * combination of tree entropies -- hence concave on the local marginal
 * polytope -- exactly when rho lies in the spanning tree polytope of K_M,
 * whose uniform point is rho_ij = 2/M. So lambda = 1/M is the LARGEST uniform
 * weight for which minimising this is a CONVEX program: every local optimum is
 * global and the answer stops depending on the start point. The Bethe weight
 * lambda = 1/2 (total edge mass C(M,2) against the M-1 a spanning tree can
 * carry) is outside that polytope for every M > 2 and coincides with 1/M at
 * M = 2.
 *
 * TWO DIFFERENCES FROM `mmi_objective`, BOTH DELIBERATE. The population loops
 * start at n = 0, the range the AMPL source states (`ni, nj in 0..F`) and the
 * one `mmi_objective` does not use, so the idle/idle cell -- the strongest
 * correlation in a closed chain -- is inside the sum; and the entropy term is
 * `mem_objective`'s body over the same restored range, which already carries
 * the sign a minimiser needs. Neither repair touches `qrf.mmi` or `qrf.mem`,
 * whose values are pinned by tests.
 *
 * NUMERICAL NOTE. The restored n = 0 cells are structurally zero: they
 * contribute 0*log(tol) = 0 to the VALUE but log(tol) ~ -13.8 to the GRADIENT,
 * so the value is insensitive to the shift while the descent direction is not.
 */
template <class T>
T bethe_objective(const std::vector<T>& x, std::size_t M, std::size_t N,
                  const std::vector<int>& K, const std::vector<int>& F, std::size_t MR) {
    const T tol = qrf_logtol<T>();
    const T lambda = num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(M));
    const QrfVars<T> v = sub_qrfvar(x, M, N, K, MR);
    T fobj = num_traits<T>::from_int(0);
    for (std::size_t m = 0; m < MR; ++m)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t ki = 0; ki < static_cast<std::size_t>(K[i]); ++ki)
                for (std::size_t j = 0; j < M; ++j) {
                    if (i == j) continue;
                    for (std::size_t kj = 0; kj < static_cast<std::size_t>(K[j]); ++kj)
                        for (std::size_t ni = 0; ni <= static_cast<std::size_t>(F[i]); ++ni)
                            for (std::size_t nj = 0; nj <= static_cast<std::size_t>(F[j]); ++nj) {
                                const T pij = v.p(i, ni, ki, j, nj, kj, m);
                                const T pii = v.p(i, ni, ki, i, ni, ki, m);
                                const T pjj = v.p(j, nj, kj, j, nj, kj, m);
                                fobj += lambda * pij * (qrf_log<T>(tol + pij) -
                                                        qrf_log<T>(tol + pii) -
                                                        qrf_log<T>(tol + pjj));
                            }
                }
    for (std::size_t m = 0; m < MR; ++m)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k)
                for (std::size_t ni = 0; ni <= static_cast<std::size_t>(F[i]); ++ni) {
                    const T pv = v.p(i, ni, k, i, ni, k, m);
                    fobj += pv * qrf_log<T>(tol + pv);
                }
    return fobj;
}

/**
 * Gradient of bethe_objective.
 *
 * df/dp_ij = lambda*(log p_ij' + p_ij/p_ij' - log p_ii' - log p_jj') for
 * i != j, and df/dp_ii = log p_ii' + p_ii/p_ii'
 *                        - (lambda/p_ii')*sum_{j!=i,kj,nj}(p_ij + p_ji),
 * a primed denominator standing for the shifted value. The second sum is
 * accumulated by the scatter below, which visits both orderings of every pair.
 */
template <class T>
std::vector<T> bethe_gradient(const std::vector<T>& x, std::size_t M, std::size_t N,
                              const std::vector<int>& K, const std::vector<int>& F, std::size_t MR,
                              const std::vector<long>& idx) {
    const T tol = qrf_logtol<T>();
    const T lambda = num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(M));
    const QrfVars<T> v = sub_qrfvar(x, M, N, K, MR);
    std::vector<T> g(x.size(), num_traits<T>::from_int(0));
    for (std::size_t m = 0; m < MR; ++m)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t ki = 0; ki < static_cast<std::size_t>(K[i]); ++ki)
                for (std::size_t j = 0; j < M; ++j) {
                    if (i == j) continue;
                    for (std::size_t kj = 0; kj < static_cast<std::size_t>(K[j]); ++kj)
                        for (std::size_t ni = 0; ni <= static_cast<std::size_t>(F[i]); ++ni) {
                            const T pii = v.p(i, ni, ki, i, ni, ki, m);
                            const long iii = idx[v.p2_at(i, ni, ki, i, ni, ki, m)];
                            for (std::size_t nj = 0; nj <= static_cast<std::size_t>(F[j]); ++nj) {
                                const T pij = v.p(i, ni, ki, j, nj, kj, m);
                                const T pjj = v.p(j, nj, kj, j, nj, kj, m);
                                const long iij = idx[v.p2_at(i, ni, ki, j, nj, kj, m)];
                                const long ijj = idx[v.p2_at(j, nj, kj, j, nj, kj, m)];
                                if (iij >= 0)
                                    g[iij] += lambda * (qrf_log<T>(tol + pij) -
                                                        qrf_log<T>(tol + pii) -
                                                        qrf_log<T>(tol + pjj) + pij / (tol + pij));
                                if (iii >= 0) g[iii] -= lambda * pij / (tol + pii);
                                if (ijj >= 0) g[ijj] -= lambda * pij / (tol + pjj);
                            }
                        }
                }
    for (std::size_t m = 0; m < MR; ++m)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k)
                for (std::size_t ni = 0; ni <= static_cast<std::size_t>(F[i]); ++ni) {
                    const T pv = v.p(i, ni, k, i, ni, k, m);
                    const long ip = idx[v.p2_at(i, ni, k, i, ni, k, m)];
                    if (ip >= 0)
                        g[ip] += qrf_log<T>(tol + pv) + pv / (tol + pv);
                }
    return g;
}

/** The utilizations and queue lengths read off an optimal pair tensor. */
template <class T>
struct QrfMetrics {
    std::vector<T> UN, QN;
    /**
     * The ALPHA-WEIGHTED diagonal marginal mean: the mean number of jobs
     * actually in service, E[min(n,c)] at a c-server station, E[n] at a delay
     * and P(n >= 1) where alpha is 1. It is what the departure rate is
     * proportional to, since alpha(i,n) scales the completion rate, so a
     * station's throughput is BN/stime exactly at the relaxed point. Equal to
     * UN on the alpha-free arms, which is why they need no separate readout.
     */
    std::vector<T> BN;
};

/**
 * extract_results: the diagonal marginals of the optimal tensor, plus the
 * alpha-weighted mean BN. `alpha` may be empty, meaning load independent.
 */
template <class T>
QrfMetrics<T> qrf_extract_results(const QrfVars<T>& v, std::size_t M, const std::vector<int>& K,
                                  const std::vector<int>& F, std::size_t MR,
                                  const Matrix<T>* alpha = nullptr) {
    QrfMetrics<T> out;
    out.UN.assign(M, num_traits<T>::from_int(0));
    out.QN.assign(M, num_traits<T>::from_int(0));
    out.BN.assign(M, num_traits<T>::from_int(0));
    for (std::size_t ti = 0; ti < M; ++ti)
        for (std::size_t m = 0; m < MR; ++m)
            for (std::size_t ni = 1; ni <= static_cast<std::size_t>(F[ti]); ++ni) {
                T a = num_traits<T>::from_int(1);
                if (alpha != nullptr && alpha->rows() > static_cast<int>(ti) &&
                    alpha->cols() >= static_cast<int>(ni))
                    a = (*alpha)(ti, ni - 1);
                for (std::size_t ki = 0; ki < static_cast<std::size_t>(K[ti]); ++ki) {
                    const T pv = v.p(ti, ni, ki, ti, ni, ki, m);
                    out.UN[ti] += pv;
                    out.QN[ti] += num_traits<T>::from_int(static_cast<long>(ni)) * pv;
                    out.BN[ti] += a * pv;
                }
            }
    return out;
}

/**
 * Rows of a maximal linearly independent subset of A, by pivoted
 * Gram-Schmidt.
 *
 * The reference uses column-pivoted QR of A^T, which at each step retains the
 * remaining column of largest residual norm; this is the same greedy selection
 * written out, and only the index SET is used downstream. Selecting by
 * conditioning rather than by first-encountered matters: the QRF equality block
 * is heavily redundant (SYMMETRY states every pair twice, ZERO / MARGINALS /
 * UEFF overlap), carrying roughly twice as many rows as its rank.
 */
template <class T>
std::vector<std::size_t> qrf_independent_rows(const Matrix<T>& A, double tol = -1.0) {
    const std::size_t m = A.rows(), n = A.cols();
    if (m == 0) return std::vector<std::size_t>();
    std::vector<std::vector<double>> res(m, std::vector<double>(n, 0.0));
    std::vector<double> nrm(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) res[i][j] = num_traits<T>::to_double(A(i, j));
        for (std::size_t j = 0; j < n; ++j) nrm[i] += res[i][j] * res[i][j];
    }
    double first = 0.0;
    for (std::size_t i = 0; i < m; ++i) first = std::max(first, std::sqrt(nrm[i]));
    if (first == 0.0) return std::vector<std::size_t>();
    // The reference's default rank tolerance: max(shape) * eps * |R00|, and
    // |R00| is the largest column norm, which is what `first` holds.
    const double rtol = (tol > 0.0) ? tol
                                    : static_cast<double>(std::max(m, n)) * 2.220446049250313e-16 *
                                          first;

    std::vector<bool> taken(m, false);
    std::vector<std::size_t> keep;
    std::vector<std::vector<double>> basis;
    for (std::size_t step = 0; step < std::min(m, n); ++step) {
        std::size_t best = m;
        double bestn = rtol;
        for (std::size_t i = 0; i < m; ++i) {
            if (taken[i]) continue;
            const double v = std::sqrt(std::max(0.0, nrm[i]));
            if (v > bestn) {
                bestn = v;
                best = i;
            }
        }
        if (best == m) break;
        taken[best] = true;
        keep.push_back(best);
        std::vector<double> q = res[best];
        const double qn = std::sqrt(nrm[best]);
        for (std::size_t j = 0; j < n; ++j) q[j] /= qn;
        basis.push_back(q);
        for (std::size_t i = 0; i < m; ++i) {
            if (taken[i]) continue;
            double dot = 0.0;
            for (std::size_t j = 0; j < n; ++j) dot += res[i][j] * q[j];
            double newn = 0.0;
            for (std::size_t j = 0; j < n; ++j) {
                res[i][j] -= dot * q[j];
                newn += res[i][j] * res[i][j];
            }
            nrm[i] = newn;
        }
    }
    std::sort(keep.begin(), keep.end());
    return keep;
}

/** An affine residual map recovered as (A, b) with fn(x) = A x - b. */
template <class T>
struct QrfAffine {
    Matrix<T> A;
    std::vector<T> b;
};

/**
 * Recover (A, b) from an affine residual map.
 *
 * Every constraint of the QRF inventory is LINEAR in the decision vector --
 * only the objectives are nonlinear -- so the matrix form is exact, not a
 * linearization. Affinity is VERIFIED at a probe point and a violation is
 * raised rather than tolerated, because a silently non-affine callback would
 * make the recovered matrix wrong everywhere except at the probe.
 */
template <class T, class Fn>
QrfAffine<T> qrf_affine_matrices(Fn fn, std::size_t n) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> basis(n, zero);
    const std::vector<T> r0 = fn(basis);
    QrfAffine<T> out;
    out.A = Matrix<T>(r0.size(), n, zero);
    out.b.assign(r0.size(), zero);
    for (std::size_t i = 0; i < r0.size(); ++i) out.b[i] = -r0[i];
    const T one = num_traits<T>::from_int(1);
    for (std::size_t col = 0; col < n; ++col) {
        basis[col] = one;
        const std::vector<T> rc = fn(basis);
        basis[col] = zero;
        for (std::size_t i = 0; i < r0.size(); ++i) out.A(i, col) = rc[i] - r0[i];
    }
    if (!r0.empty() && n > 0) {
        std::vector<T> probe(n, zero);
        for (std::size_t j = 0; j < n; ++j)
            probe[j] = num_traits<T>::from_double(
                0.1 + 0.8 * (n == 1 ? 0.0 : static_cast<double>(j) / static_cast<double>(n - 1)));
        const std::vector<T> rp = fn(probe);
        double err = 0.0;
        for (std::size_t i = 0; i < r0.size(); ++i) {
            T v = -out.b[i];
            for (std::size_t j = 0; j < n; ++j) v += out.A(i, j) * probe[j];
            err = std::max(err, std::fabs(num_traits<T>::to_double(rp[i] - v)));
        }
        if (err > 1e-9)
            throw InputError(
                "qrf_affine_matrices: the constraint residuals are not affine in x; the matrix "
                "form recovered here would be wrong away from the probe");
    }
    return out;
}

/** The equality block with its dependent rows dropped. */
template <class T>
struct QrfReduced {
    Matrix<T> A;
    std::vector<T> b;
    std::vector<std::size_t> keep;
};

/**
 * Drop the linearly dependent equality rows, keeping the feasible set exact.
 *
 * An INCONSISTENT system is refused rather than reduced: rank([A|b]) above
 * rank(A) means the polytope is empty, which is a modelling error, and
 * discarding the offending rows would return numbers for a model nobody wrote.
 */
template <class T>
QrfReduced<T> qrf_reduce_equalities(const Matrix<T>& A, const std::vector<T>& b) {
    QrfReduced<T> out;
    out.keep = qrf_independent_rows(A);
    if (out.keep.size() == A.rows()) {
        out.A = A;
        out.b = b;
        return out;
    }
    Matrix<T> aug(A.rows(), A.cols() + 1, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i) {
        for (std::size_t j = 0; j < A.cols(); ++j) aug(i, j) = A(i, j);
        aug(i, A.cols()) = b[i];
    }
    if (qrf_independent_rows(aug).size() > out.keep.size())
        throw InputError(
            "qrf_reduce_equalities: the equality system is inconsistent (the augmented matrix has "
            "the higher rank), i.e. the polytope is empty");
    out.A = Matrix<T>(out.keep.size(), A.cols(), num_traits<T>::from_int(0));
    out.b.assign(out.keep.size(), num_traits<T>::from_int(0));
    for (std::size_t r = 0; r < out.keep.size(); ++r) {
        for (std::size_t j = 0; j < A.cols(); ++j) out.A(r, j) = A(out.keep[r], j);
        out.b[r] = b[out.keep[r]];
    }
    return out;
}

/**
 * The polytope of a QRF instance, as an LpModel over the box [0,1]^n.
 *
 * Shared by the phase 1 and by the conditional-gradient loop, so the two can
 * never disagree about which set they are working over.
 */
template <class T>
lp::LpModel<T> qrf_polytope(const Matrix<T>& Aeq, const std::vector<T>& beq,
                            const Matrix<T>& Aub, const std::vector<T>& bub, std::size_t n) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    lp::LpModel<T> model(n);
    model.set_maximize(false);  // a MINIMIZATION throughout; see solve_qrf_nlp
    for (std::size_t j = 0; j < n; ++j) {
        model.set_bounds(j, zero, one);
        model.set_cost(j, zero);
    }
    for (std::size_t i = 0; i < Aeq.rows(); ++i) {
        model.row_clear();
        for (std::size_t j = 0; j < n; ++j)
            if (Aeq(i, j) != zero) model.row_add(j, Aeq(i, j));
        model.emit_eq(beq[i]);
    }
    for (std::size_t i = 0; i < Aub.rows(); ++i) {
        model.row_clear();
        for (std::size_t j = 0; j < n; ++j)
            if (Aub(i, j) != zero) model.row_add(j, Aub(i, j));
        model.emit_le(bub[i]);
    }
    return model;
}

/**
 * Orthonormal basis of null(A), as an (n x d) matrix.
 *
 * The reference calls `scipy.linalg.null_space`, an SVD; this builds the same
 * space by modified Gram-Schmidt -- orthonormalize the rows of A, then sweep
 * the canonical directions and keep each one whose residual against the row
 * space and the basis so far is nontrivial. The basis is not the SVD's, but the
 * SPACE is, and the caller only ever uses the space: `x0 + Z t` ranges over the
 * same affine set whichever orthonormal basis Z carries. Gram-Schmidt is used
 * rather than `util/svd.h` because that header needs LAPACK, which this call
 * path must not require.
 */
template <class T>
Matrix<T> qrf_null_space(const Matrix<T>& A, std::size_t n) {
    std::vector<std::vector<double>> basis;  // orthonormal, spanning row(A)
    const double eps = 1e-10;
    for (std::size_t r = 0; r < A.rows(); ++r) {
        std::vector<double> v(n, 0.0);
        for (std::size_t j = 0; j < n; ++j) v[j] = num_traits<T>::to_double(A(r, j));
        for (std::size_t b = 0; b < basis.size(); ++b) {
            double dot = 0.0;
            for (std::size_t j = 0; j < n; ++j) dot += v[j] * basis[b][j];
            for (std::size_t j = 0; j < n; ++j) v[j] -= dot * basis[b][j];
        }
        double nv = 0.0;
        for (std::size_t j = 0; j < n; ++j) nv += v[j] * v[j];
        nv = std::sqrt(nv);
        if (nv <= eps) continue;
        for (std::size_t j = 0; j < n; ++j) v[j] /= nv;
        basis.push_back(v);
    }
    const std::size_t rank = basis.size();
    std::vector<std::vector<double>> nullb;
    for (std::size_t c = 0; c < n && rank + nullb.size() < n; ++c) {
        std::vector<double> v(n, 0.0);
        v[c] = 1.0;
        for (std::size_t b = 0; b < basis.size(); ++b) {
            double dot = 0.0;
            for (std::size_t j = 0; j < n; ++j) dot += v[j] * basis[b][j];
            for (std::size_t j = 0; j < n; ++j) v[j] -= dot * basis[b][j];
        }
        for (std::size_t b = 0; b < nullb.size(); ++b) {
            double dot = 0.0;
            for (std::size_t j = 0; j < n; ++j) dot += v[j] * nullb[b][j];
            for (std::size_t j = 0; j < n; ++j) v[j] -= dot * nullb[b][j];
        }
        double nv = 0.0;
        for (std::size_t j = 0; j < n; ++j) nv += v[j] * v[j];
        nv = std::sqrt(nv);
        if (nv <= 1e-8) continue;
        for (std::size_t j = 0; j < n; ++j) v[j] /= nv;
        nullb.push_back(v);
    }
    Matrix<T> Z(n, nullb.size(), num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < nullb.size(); ++j)
        for (std::size_t i = 0; i < n; ++i) Z(i, j) = num_traits<T>::from_double(nullb[j][i]);
    return Z;
}

/**
 * Minimize a convex objective over {Aeq x = beq, Aub x <= bub, 0 <= x <= 1},
 * starting from a feasible point, by FRANK-WOLFE.
 *
 * THIS IS A DELIBERATE SUBSTITUTION FOR THE REFERENCE'S SLSQP, and it is the
 * one place in this port where a different algorithm is the right answer rather
 * than a compromise. The reference eliminates the equalities with a null-space
 * basis and runs SLSQP on the reduced problem. That works because SLSQP is
 * gradient-based; the augmented-Lagrangian path available here has a
 * derivative-free Nelder-Mead inner solve, and started AT A VERTEX of the
 * polytope -- which is exactly what the phase-1 LP returns -- its first simplex
 * steps leave the feasible set, are penalized, and it does not move at all.
 * Measured on a two-phase / one-phase pair at N = 2: it returned its start
 * point, an objective of 0.6931 against the reference's 0.6438, and a
 * utilization of 1.0 against 0.8.
 *
 * The feasible set is a polytope, so conditional gradient applies directly and
 * gives more than convergence: at each step the linear minimization over the
 * SAME polytope yields the FRANK-WOLFE GAP <grad f(x), x - s>, which is the
 * exact stationarity measure -- zero iff no feasible descent direction exists
 * at x. On MEM, which is convex, it is additionally an upper bound on
 * f(x) - f(x*). It subsumes the reference's one-shot "did a feasible descent
 * direction exist at the start?" probe: the same question is asked at every
 * iterate, and the answer is what terminates the loop.
 *
 * No null-space reduction is needed, since the LP carries the equalities
 * itself.
 *
 * @param objective convex objective
 * @param gradient  its gradient
 * @param x0        a feasible point, from qrf_feasible_start
 */
template <class T, class Obj, class Grad>
std::vector<T> solve_qrf_nlp_lp(Obj objective, Grad gradient, const std::vector<T>& x0,
                                const lp::LpModel<T>& polytope, const std::string& name,
                                unsigned max_iter = 200, double gap_tol = 1e-10) {
    static_assert(num_traits<T>::has_transcendental, "solve_qrf_nlp needs a log()");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = x0.size();
    std::vector<T> x = x0;

    // The polytope, assembled once: only the LP costs change between rounds.
    lp::LpModel<T> base = polytope;
    // `LpModel` MAXIMIZES by default. The conditional-gradient step is a
    // MINIMIZATION of <grad f(x), d>, and taking the default silently returns
    // the ascent vertex: the Frank-Wolfe gap then comes out NEGATIVE (-0.184 on
    // the two-phase instance), the loop reads that as "already optimal" and the
    // routine hands back its phase-1 start.
    base.set_maximize(false);

    for (unsigned it = 0; it < max_iter; ++it) {
        const std::vector<T> g = gradient(x);
        lp::LpModel<T> model = base;
        for (std::size_t j = 0; j < n; ++j) model.set_cost(j, g[j]);
        // lp_solve, not simplex_solve: the BAS polytope outgrows the dense
        // tableau quickly, and this is the same dispatch mapqn_qr_bounds_bas uses.
        const lp::LpSolution<T> sol = lp::lp_solve(model);
        if (sol.status != lp::LpStatus::Optimal)
            throw InputError(name +
                             ": the linear minimization over the QRF polytope failed, so no "
                             "descent direction and no optimality certificate can be had");

        // The Frank-Wolfe gap bounds f(x) - f(x*) from above.
        double gap = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            gap += num_traits<T>::to_double(g[j] * (x[j] - sol.x[j]));
        if (gap <= gap_tol) break;

        // f is convex along the segment, so a golden-section search on [0,1]
        // is exact to the bracket it reports.
        const double invphi = 0.6180339887498949;
        double a = 0.0, b = 1.0;
        double c = b - invphi * (b - a), d = a + invphi * (b - a);
        auto at = [&](double t) {
            std::vector<T> y(n, zero);
            const T tt = num_traits<T>::from_double(t);
            for (std::size_t j = 0; j < n; ++j) y[j] = x[j] + tt * (sol.x[j] - x[j]);
            return num_traits<T>::to_double(objective(y));
        };
        const double f_cur = at(0.0);
        double fc = at(c), fd = at(d);
        for (int k = 0; k < 60 && (b - a) > 1e-12; ++k) {
            if (fc < fd) {
                b = d;
                d = c;
                fd = fc;
                c = b - invphi * (b - a);
                fc = at(c);
            } else {
                a = c;
                c = d;
                fc = fd;
                d = a + invphi * (b - a);
                fd = at(d);
            }
        }
        const double gamma = 0.5 * (a + b);
        if (gamma <= 0.0) break;  // the line search found no improvement to take

        // TERMINATE ON PROGRESS, not on the gap alone. The Frank-Wolfe gap is
        // an exact stationarity measure in exact arithmetic, but it is computed
        // as a sum of n floating-point products and settles at the NOISE FLOOR
        // of that sum rather than at zero. Measured on the ragged two-phase /
        // one-phase instance at N = 2: the gap is 1.45 at the first step, which
        // is real and drops f from 1.375 to 0.8246, and from the second step on
        // it rattles between 1e-10 and 2e-8 -- numerically zero, but above the
        // 1e-10 `gap_tol`. The loop then spent its whole 200-iteration budget
        // there, solving an LP per round and taking gamma ~ 3e-8 steps that
        // moved neither x nor f. Scaling `gap_tol` with n would trade one
        // arbitrary constant for another; a step that does not measurably
        // reduce the objective is the honest end of the useful run, whatever
        // the gap says. The returned point is unchanged -- this only stops the
        // spinning.
        const double f_new = at(gamma);
        if (!(f_new < f_cur - 1e-12 * (1.0 + std::fabs(f_cur)))) break;

        const T gt = num_traits<T>::from_double(gamma);
        for (std::size_t j = 0; j < n; ++j) x[j] = x[j] + gt * (sol.x[j] - x[j]);
    }
    return x;
}

/** The same, with the polytope given as matrices rather than as an LpModel. */
template <class T, class Obj, class Grad>
std::vector<T> solve_qrf_nlp(Obj objective, Grad gradient, const std::vector<T>& x0,
                             const Matrix<T>& Aeq, const std::vector<T>& beq,
                             const Matrix<T>& Aub, const std::vector<T>& bub,
                             const std::string& name, unsigned max_iter = 200,
                             double gap_tol = 1e-10) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t n = x0.size();
    lp::LpModel<T> base(n);
    base.set_maximize(false);
    for (std::size_t j = 0; j < n; ++j) base.set_bounds(j, zero, one);
    for (std::size_t i = 0; i < Aeq.rows(); ++i) {
        base.row_clear();
        for (std::size_t j = 0; j < n; ++j)
            if (Aeq(i, j) != zero) base.row_add(j, Aeq(i, j));
        base.emit_eq(beq[i]);
    }
    for (std::size_t i = 0; i < Aub.rows(); ++i) {
        base.row_clear();
        for (std::size_t j = 0; j < n; ++j)
            if (Aub(i, j) != zero) base.row_add(j, Aub(i, j));
        base.emit_le(bub[i]);
    }
    return solve_qrf_nlp_lp(objective, gradient, x0, base, name, max_iter, gap_tol);
}

/**
 * The minimum-norm point of the polytope, from any feasible point of it.
 *
 * `0.5||x||^2` is strongly convex, so the Frank-Wolfe gap is a true bound on
 * the optimality gap and the loop terminates on it rather than on its budget.
 */
template <class T>
std::vector<T> qrf_min_norm_point(const lp::LpModel<T>& polytope, const std::vector<T>& x0,
                                  const std::string& name) {
    const T half = num_traits<T>::from_double(0.5);
    return solve_qrf_nlp_lp(
        [&](const std::vector<T>& x) {
            T s = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < x.size(); ++j) s += x[j] * x[j];
            return T(half * s);
        },
        [](const std::vector<T>& x) { return x; }, x0, polytope, name + " phase 1", 500u, 1e-12);
}

/**
 * A feasible point of an LpModel polytope: the MINIMUM-NORM one.
 *
 * THE START IS PART OF THE ANSWER HERE, and a vertex is the wrong one. Neither
 * MEM nor MMI is a convex program as the reference states them -- the MEM
 * objective is `-sum p log p`, which is CONCAVE, and it is MINIMIZED -- so
 * every method reports a stationary point fixed by where it started. A
 * zero-cost phase-1 LP lands on an arbitrary VERTEX, and a vertex is already a
 * local minimum of a concave objective: conditional gradient reports a zero
 * gap and hands the vertex straight back. Measured on the two-phase /
 * one-phase pair at N = 2: f(x0) = f(xopt) = 1.0397 and a utilization of
 * exactly 1 against the references' 0.8.
 *
 * `qrf_noblo_start.m` minimizes `0.5||x||^2` over the same rows for exactly
 * this reason ("so the returned point is the unique minimum-norm feasible
 * point instead of an arbitrary vertex"), and native Python does the same.
 * That objective IS strongly convex, so conditional gradient solves it to
 * optimality over the polytope and needs no QP backend. From that start the
 * port returns [0.8, 0.4] / [1.4, 0.6], which is native Python to the last
 * digit printed and MATLAB to 1e-3 -- MATLAB's fmincon stops on its own
 * `MaxIter = 100` slightly short of it.
 */
template <class T>
std::vector<T> qrf_feasible_start_lp(const lp::LpModel<T>& polytope, const std::string& name) {
    const T zero = num_traits<T>::from_int(0);
    lp::LpModel<T> m = polytope;
    m.set_maximize(false);
    for (std::size_t j = 0; j < m.num_vars(); ++j) m.set_cost(j, zero);
    const lp::LpSolution<T> sol = lp::lp_solve(m);
    if (sol.status != lp::LpStatus::Optimal)
        throw InputError(name +
                         ": the polytope is infeasible; a well-posed instance always has a "
                         "feasible point, so this is a modelling error");
    return qrf_min_norm_point(m, sol.x, name);
}

/**
 * A point of the polytope, as the phase 1 of `qrf_noblo_start.m`.
 *
 * The all-zero vector violates normalization by a full unit and COR1 by N^2,
 * and from there a local method terminates AT the start point. Every
 * constraint here is linear, so an LP lands on the polytope exactly; the
 * minimum-norm refinement above is what keeps it off a vertex.
 */
template <class T>
std::vector<T> qrf_feasible_start(const Matrix<T>& Aeq, const std::vector<T>& beq,
                                  const Matrix<T>& Aub, const std::vector<T>& bub,
                                  std::size_t n) {
    const lp::LpModel<T> model = qrf_polytope(Aeq, beq, Aub, bub, n);
    const lp::LpSolution<T> sol = lp::simplex_solve(model);
    if (sol.status != lp::LpStatus::Optimal)
        throw InputError(
            "qrf_feasible_start: the QRF polytope is infeasible; a well-posed instance always has "
            "one, so an empty polytope is a modelling error rather than a numerical accident");
    return qrf_min_norm_point(model, sol.x, "qrf_feasible_start");
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_QRF_COMMON_H
