/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_QRF_BAS_NLP_H
#define LINE_API_MAPQN_MAPQN_QRF_BAS_NLP_H

/**
 * `qrf_bas_mmi`, `qrf_bas_mem` and `qrf_bas_bethe`: the nonlinear bounds on the
 * BAS-BLOCKING polytope.
 *
 * Port of python/line_solver/api/mapqn/qrf_bas_nlp.py. `api/mapqn` has no
 * MATLAB implementation, so Python and the JAR are the references.
 *
 * THE POLYTOPE IS NOT RE-DERIVED HERE. It is exactly the one
 * `mapqn_qr_bounds_bas` builds and that the LP token `qrf.bas` is validated on
 * against the AMPL model; only the OBJECTIVE differs, the LP maximizing one
 * station's utilization where these minimize mutual information, negative
 * entropy or the tree-reweighted (Bethe) free entropy over the same feasible
 * set. Re-transcribing the fifteen constraint
 * families would duplicate several hundred lines whose index conventions are
 * exactly where the reference's own twin twice went wrong -- the THM30/THM3
 * population-to-index shift, and a MARGINALS sum whose upper limit was taken on
 * the 1-based index rather than on the population. Sharing the builder means a
 * family added to the LP reaches the NLP with it.
 *
 * The decision vector is therefore the LP's own, indexed through `QrBasIndex`,
 * NOT the flat `sub_qrfvar` layout of the no-blocking family. The two are
 * different orderings of different variable sets and must not be mixed: the BAS
 * index skips a population above a station's capacity, which the no-blocking
 * layout carries and pins to zero.
 *
 * CONVEXITY, AS IN THE NO-BLOCKING TWIN: MEM is convex on the polytope and has
 * a unique optimum; MMI is NOT (its `-p_ij log p_ii` terms are not), so it
 * reports a local optimum fixed by the phase-1 LP vertex. BETHE is the convex
 * combination the weight lambda = 1/M is chosen to make convex on the LOCAL
 * MARGINAL polytope; the BAS set adds the blocking families on top of the
 * marginal ones, so the objective's convexity carries but the per-configuration
 * marginal consistency that the argument rests on is not proved under ZERO5.
 * Start-point independence here is therefore MEASURED, not assumed. See
 * `mapqn_qrf_common.h` for the measurement behind all three.
 *
 * ARITHMETIC: transcendental.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mapqn/mapqn_qr_bounds_bas.h"
#include "line/api/mapqn/mapqn_qrf_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lp_highs.h"
#include "line/util/simplex.h"

namespace line {
namespace mapqn {

namespace qrfbas {

/** The BAS feasible set, built by the same calls the LP bound makes. */
template <class T>
lp::LpModel<T> bas_polytope(const QrBasParams<T>& p, const QrBasIndex& x) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    lp::LpModel<T> m(x.num_vars());
    // Families in the reference's emission order. ZERO must precede SYMMETRY,
    // which skips pairs both of whose members it pinned.
    const std::vector<char> zeroed = detail::bas_zero(p, x, m);
    detail::bas_one(p, x, m);
    detail::bas_symmetry(p, x, m, zeroed);
    detail::bas_marginals(p, x, m);
    detail::bas_ueff(p, x, m);
    detail::bas_thm1(p, x, m);
    detail::bas_thm2(p, x, m);
    detail::bas_cor1(p, x, m);
    detail::bas_thm30(p, x, m);
    detail::bas_thm3(p, x, m);
    detail::bas_thm3f(p, x, m);
    detail::bas_thm3i(p, x, m);
    detail::bas_thm3l(p, x, m);
    detail::bas_thm4(p, x, m);
    // THE [0,1] BOX IS NOT IMPOSED, AND MUST NOT BE. The reference's phase-1
    // `linprog` passes `bounds=[(0,1)]*n`, but here that upper bound is
    // IMPLIED and imposing it explicitly breaks the solve: ONE fixes each
    // station's diagonal marginal to sum to one over nonnegative entries, so
    // every diagonal entry is at most one; MARGINALS ties each off-diagonal
    // entry to its own diagonal, so those are too; and UEFF writes `e` as a sum
    // of the same entries. Adding 246 redundant upper bounds to the dense
    // tableau instead made `lp_solve` report Optimal on a point whose largest
    // coordinate was 7.3e9 -- an out-of-bounds answer with a success status,
    // which then propagated into utilizations of 92 and queue lengths of 334.
    // The free-upper formulation is the one `mapqn_qr_bounds_bas` solves and
    // validates against MATLAB, so it is the one used here.
    (void)zero;
    (void)one;
    return m;
}

/**
 * The (ij, ii, jj) column triples of the MI objective, i != j, ni, nj >= n_from.
 *
 * n_from is 0 for both callers, the reference's own range: the idle cell
 * carries the strongest correlation in a closed chain, and BETHE's entropy term
 * spans the same range so the two blocks agree.
 */
inline void mmi_terms(const QrBasIndex& x, const std::vector<int>& F, std::vector<std::size_t>* ij,
                      std::vector<std::size_t>* ii, std::vector<std::size_t>* jj,
                      int n_from = 1) {
    for (int m = 0; m < x.MR; ++m)
        for (int i = 0; i < x.M; ++i)
            for (int ki = 0; ki < x.K[i]; ++ki)
                for (int j = 0; j < x.M; ++j) {
                    if (i == j) continue;
                    for (int kj = 0; kj < x.K[j]; ++kj)
                        for (int ni = n_from; ni <= F[i]; ++ni)
                            for (int nj = n_from; nj <= F[j]; ++nj) {
                                ij->push_back(x.p2(i, ni, ki, j, nj, kj, m));
                                ii->push_back(x.p2(i, ni, ki, i, ni, ki, m));
                                jj->push_back(x.p2(j, nj, kj, j, nj, kj, m));
                            }
                }
}

/** The diagonal columns of the MEM objective, ni >= n_from (1 for MEM, 0 for BETHE). */
inline std::vector<std::size_t> mem_terms(const QrBasIndex& x, const std::vector<int>& F,
                                          int n_from = 1) {
    std::vector<std::size_t> out;
    for (int m = 0; m < x.MR; ++m)
        for (int i = 0; i < x.M; ++i)
            for (int k = 0; k < x.K[i]; ++k)
                for (int ni = n_from; ni <= F[i]; ++ni) out.push_back(x.p2(i, ni, k, i, ni, k, m));
    return out;
}

/** Which functional `mapqn_qrf_bas` minimizes over the BAS polytope. */
enum class Objective { Mmi, Mem, Bethe };

}  // namespace qrfbas

/**
 * Solve one BAS-blocking NLP bound.
 *
 * @param p        the same parameters the LP token `qrf.bas` takes
 * @param obj      which functional to minimize over the polytope
 * @param max_iter Frank-Wolfe iteration cap. EVERY ITERATE IS FEASIBLE, since
 *                 each step is a convex combination of two points of the
 *                 polytope, so a truncated run returns a worse bound but never
 *                 an invalid one. That is what makes a small cap a legitimate
 *                 way to keep a test cheap: the BAS polytope has hundreds of
 *                 columns and each iteration costs one LP plus a line search.
 */
template <class T>
QrfMetrics<T> mapqn_qrf_bas(const QrBasParams<T>& p, qrfbas::Objective obj,
                            unsigned max_iter = 200) {
    p.validate();
    const T tol = qrf_logtol<T>();
    const QrBasIndex x(p.M, p.N, p.K, p.MR);
    const bool mem = (obj == qrfbas::Objective::Mem);
    const bool bethe = (obj == qrfbas::Objective::Bethe);
    // lambda = 1/M: half the uniform point 2/M of the spanning-tree polytope of
    // K_M, the largest uniform edge weight at which the tree-reweighted entropy
    // is concave and the program therefore convex.
    const T lam = num_traits<T>::from_int(1) / num_traits<T>::from_int(p.M);
    const std::string name =
        bethe ? "qrf_bas_bethe" : (mem ? "qrf_bas_mem" : "qrf_bas_mmi");

    const lp::LpModel<T> poly = qrfbas::bas_polytope(p, x);

    // The MI block spans 0..F, the reference's range; MEM's spans 1..F, except
    // under BETHE, which combines the two and needs one range across both.
    std::vector<std::size_t> ij, ii, jj, diag;
    if (mem || bethe) {
        diag = qrfbas::mem_terms(x, p.F, bethe ? 0 : 1);
        if (diag.empty())
            throw InputError(name +
                             ": the model has no marginal variables, so the objective is empty");
    }
    if (!mem) {
        qrfbas::mmi_terms(x, p.F, &ij, &ii, &jj, 0);
        if (ij.empty())
            throw InputError(name +
                             ": the model has no off-diagonal joint variables, so the mutual "
                             "information objective is empty");
    }

    // The MI block is scaled by lambda under BETHE and by one under MMI, and
    // is absent under MEM; the entropy block is present under MEM and BETHE.
    // Written as `-> T` explicitly: a deduced return type over an expression-
    // template Rational returns a node holding references to dead temporaries.
    const T mi_scale = bethe ? lam : num_traits<T>::from_int(1);
    auto objective = [&](const std::vector<T>& z) -> T {
        T f = num_traits<T>::from_int(0);
        for (std::size_t t = 0; t < diag.size(); ++t) {
            const T pv = z[diag[t]];
            f += T(pv * qrf_log<T>(tol + pv));
        }
        for (std::size_t t = 0; t < ij.size(); ++t) {
            const T pij = z[ij[t]], pii = z[ii[t]], pjj = z[jj[t]];
            f += T(mi_scale * pij *
                   (qrf_log<T>(tol + pij) - qrf_log<T>(tol + pii) - qrf_log<T>(tol + pjj)));
        }
        return f;
    };
    auto gradient = [&](const std::vector<T>& z) -> std::vector<T> {
        std::vector<T> g(z.size(), num_traits<T>::from_int(0));
        for (std::size_t t = 0; t < diag.size(); ++t) {
            const T pv = z[diag[t]];
            g[diag[t]] += T(qrf_log<T>(tol + pv) + pv / (tol + pv));
        }
        for (std::size_t t = 0; t < ij.size(); ++t) {
            const T pij = z[ij[t]], pii = z[ii[t]], pjj = z[jj[t]];
            g[ij[t]] += T(mi_scale * (qrf_log<T>(tol + pij) - qrf_log<T>(tol + pii) -
                                      qrf_log<T>(tol + pjj) + pij / (tol + pij)));
            g[ii[t]] -= T(mi_scale * pij / (tol + pii));
            g[jj[t]] -= T(mi_scale * pij / (tol + pjj));
        }
        return g;
    };

    const std::vector<T> x0 = qrf_feasible_start_lp(poly, name);
    const std::vector<T> xopt = solve_qrf_nlp_lp(objective, gradient, x0, poly, name, max_iter);

    // UTILIZATION, not occupancy. UN used to sum the diagonal p2 over ALL
    // blocking configurations, i.e. P(n_i >= 1) with the BLOCKED ones included.
    // A blocked BAS server holds a job it has already finished and does no work,
    // so that is occupancy: on the M=2, N=3, F=[2 3] cyclic model it reported
    // U2 = 1 where the exact utilization is 7/15, which the LP over the SAME
    // polytope already returns. e carries the right quantity -- bas_ueff pins
    // e(i,ki) to the mass with n_i >= 1 in the configurations where i is NOT
    // blocked.
    //
    // The 1/M matches THIS port's UEFF, which emits one row per (i,ki) with j
    // summed INSIDE, leaving e scaled by M; the python and JAR ports emit one
    // row per (j,i,ki) and carry no such factor. QN stays on the diagonal p2
    // over every configuration, blocked included, because a blocked job is still
    // held at the station and counts towards its population.
    QrfMetrics<T> out;
    out.UN.assign(static_cast<std::size_t>(p.M), num_traits<T>::from_int(0));
    out.QN.assign(static_cast<std::size_t>(p.M), num_traits<T>::from_int(0));
    const T inv_m = num_traits<T>::from_int(1) / num_traits<T>::from_int(p.M);
    for (int i = 0; i < p.M; ++i)
        for (int ki = 0; ki < p.K[i]; ++ki)
            out.UN[static_cast<std::size_t>(i)] += T(xopt[x.e(i, ki)] * inv_m);
    for (int i = 0; i < p.M; ++i)
        for (int m = 0; m < p.MR; ++m)
            for (int ni = 1; ni <= p.F[i]; ++ni)
                for (int ki = 0; ki < p.K[i]; ++ki)
                    out.QN[static_cast<std::size_t>(i)] +=
                        T(num_traits<T>::from_int(ni) * xopt[x.p2(i, ni, ki, i, ni, ki, m)]);
    return out;
}

/** Minimum-mutual-information bound on the BAS-blocking polytope. */
template <class T>
QrfMetrics<T> mapqn_qrf_bas_mmi(const QrBasParams<T>& p, unsigned max_iter = 200) {
    return mapqn_qrf_bas(p, qrfbas::Objective::Mmi, max_iter);
}

/** Maximum-entropy bound on the BAS-blocking polytope. */
template <class T>
QrfMetrics<T> mapqn_qrf_bas_mem(const QrBasParams<T>& p, unsigned max_iter = 200) {
    return mapqn_qrf_bas(p, qrfbas::Objective::Mem, max_iter);
}

/**
 * Tree-reweighted (Bethe) free entropy bound on the BAS-blocking polytope.
 *
 * `lambda*sum_{i!=j} I(n_i;n_j) - sum_i H(n_i)` at lambda = 1/M, the objective
 * of `qrf_noblo_bethe` evaluated over the BAS decision vector: the blocking
 * configurations and the per-station capacities enter through the ranges alone.
 */
template <class T>
QrfMetrics<T> mapqn_qrf_bas_bethe(const QrBasParams<T>& p, unsigned max_iter = 200) {
    return mapqn_qrf_bas(p, qrfbas::Objective::Bethe, max_iter);
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_QRF_BAS_NLP_H
