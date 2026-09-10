/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_DMAP_OPTIM_DIST_H
#define LINE_API_MAM_DMAP_OPTIM_DIST_H

/**
 * The DISCRETE twins of `map_optim_dist` / `map_optim_dist_acf`: fit a D-MAP's
 * D1 by minimizing a distance to a reference, with D0 held fixed.
 *
 * Port of matlab/lib/kpctoolbox/dmap/dmap_optim_dist.m and
 * dmap_optim_dist_acf.m (twins in
 * python/line_solver/api/mapdist/discrete.py). The distances themselves --
 * `dmap_dist`, `dmap_dist_acf` -- are already in `dmap.h`.
 *
 * EVERYTHING IS THE CONTINUOUS CASE WITH `(-D0)` REPLACED BY `(I - D0)`, and
 * that substitution is the whole of the discrete/continuous difference here. In
 * continuous time the embedded kernel is `(-D0)^-1 D1` and the row sums of D0
 * and D1 cancel; in discrete time the phase either moves without an arrival
 * (D0) or with one (D1) at every SLOT, so the two together form a stochastic
 * matrix, the kernel is `(I - D0)^-1 D1`, and the row sums are ONE rather than
 * zero. Both constraint blocks change accordingly:
 *
 *   - `alB (I - B0)^-1 B1 = alB`, alB stationary at arrivals;
 *   - the row sums of B1 are those of `(I - B0)`.
 *
 * There is no convex quadratic branch here. The continuous `map_optim_dist`
 * takes one at L = 1 because `quadprog` applies; the discrete reference calls
 * `fmincon` at every lag, so every result is a LOCAL optimum and none claims
 * otherwise.
 *
 * ARITHMETIC: transcendental, inherited from the distances.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/dmap.h"
#include "line/api/mam/map_optim_dist.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace dmapoptdetail {

/**
 * The discrete constraint block.
 *
 * `(I - B0)` where the continuous form has `(-B0)`, and the row-sum target is
 * the row sums of `(I - B0)` rather than of `-B0`.
 */
template <class T>
void build_constraints_d(const Matrix<T>& B0, const std::vector<T>& alB, Matrix<T>* Aeq,
                         std::vector<T>* beq) {
    const std::size_t n = B0.rows();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    Matrix<T> ImB0(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) ImB0(i, j) = (i == j ? one : zero) - B0(i, j);

    // (alB (I-B0)^-1)' solves (I-B0)' x = alB'.
    Matrix<T> tr(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) tr(i, j) = ImB0(j, i);
    const std::vector<T> row = solve(tr, alB);

    Matrix<T> rowM(1, n, zero);
    for (std::size_t j = 0; j < n; ++j) rowM(0, j) = row[j];
    Matrix<T> I(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) I(i, i) = one;
    Matrix<T> ones1(1, n, one);

    const Matrix<T> top = optdistdetail::mkron(I, rowM);
    const Matrix<T> bot = optdistdetail::mkron(ones1, I);
    *Aeq = Matrix<T>(2 * n, n * n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n * n; ++j) {
            (*Aeq)(i, j) = top(i, j);
            (*Aeq)(n + i, j) = bot(i, j);
        }
    beq->assign(2 * n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        (*beq)[i] = alB[i];
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += ImB0(i, j);
        (*beq)[n + i] = s;
    }
}

}  // namespace dmapoptdetail

/**
 * Fit B1 minimizing the lag-L joint-PMF distance to `a`, with B0 fixed.
 *
 * Always a LOCAL optimum: the discrete reference has no convex branch.
 */
template <class T>
MapOptimDist<T> dmap_optim_dist(const Dmap<T>& a, const std::vector<T>& alA, const Matrix<T>& B0,
                                const std::vector<T>& alB, unsigned L) {
    static_assert(num_traits<T>::has_transcendental, "dmap_optim_dist needs the D-MAP distances");
    const std::size_t n = B0.rows();
    if (n == 0 || B0.cols() != n) throw InputError("dmap_optim_dist: B0 must be square");
    if (alB.size() != n) throw InputError("dmap_optim_dist: alB has the wrong length");
    if (alA.size() != a.D0.rows()) throw InputError("dmap_optim_dist: alA has the wrong length");
    if (L == 0) throw InputError("dmap_optim_dist: at least one lag is required");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T lo = num_traits<T>::from_double(1e-6);
    Matrix<T> Aeq;
    std::vector<T> beq;
    dmapoptdetail::build_constraints_d(B0, alB, &Aeq, &beq);

    auto obj = [&a, &B0, &alA, &alB, L, n](const std::vector<T>& x) {
        Dmap<T> b;
        b.D0 = B0;
        b.D1 = optdistdetail::unvec(x, n);
        return dmap_dist(a, b, L, alA, alB);
    };
    // Start on the row-sum constraint: (I - B0)'s row sums, spread evenly.
    Matrix<T> seed(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += (i == j ? one : zero) - B0(i, j);
        for (std::size_t j = 0; j < n; ++j)
            seed(i, j) = s / num_traits<T>::from_int(static_cast<long>(n));
    }

    MapOptimDist<T> out;
    const std::vector<T> sol =
        optdistdetail::constrained_min(obj, Aeq, beq, optdistdetail::vec(seed), lo);
    out.B1 = optdistdetail::unvec(sol, n);
    out.d = obj(sol);
    out.global = false;
    return out;
}

/**
 * Fit B1 minimizing the AUTOCORRELATION distance, with B0 fixed.
 *
 * As in the continuous twin, the distance is RE-EVALUATED at the returned B1
 * rather than taken from the optimizer, so the number reported is the distance
 * of the D-MAP the caller was handed.
 */
template <class T>
MapOptimDist<T> dmap_optim_dist_acf(const Dmap<T>& a, const std::vector<T>& alA,
                                    const Matrix<T>& B0, const std::vector<T>& alB) {
    static_assert(num_traits<T>::has_transcendental,
                  "dmap_optim_dist_acf needs the D-MAP distances");
    const std::size_t n = B0.rows();
    if (n == 0 || B0.cols() != n) throw InputError("dmap_optim_dist_acf: B0 must be square");
    if (alB.size() != n) throw InputError("dmap_optim_dist_acf: alB has the wrong length");
    if (alA.size() != a.D0.rows())
        throw InputError("dmap_optim_dist_acf: alA has the wrong length");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T lo = num_traits<T>::from_double(1e-6);
    Matrix<T> Aeq;
    std::vector<T> beq;
    dmapoptdetail::build_constraints_d(B0, alB, &Aeq, &beq);

    auto obj = [&a, &B0, &alA, &alB, n](const std::vector<T>& x) {
        Dmap<T> b;
        b.D0 = B0;
        b.D1 = optdistdetail::unvec(x, n);
        return dmap_dist_acf(a, b, alA, alB);
    };
    Matrix<T> seed(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += (i == j ? one : zero) - B0(i, j);
        for (std::size_t j = 0; j < n; ++j)
            seed(i, j) = s / num_traits<T>::from_int(static_cast<long>(n));
    }

    MapOptimDist<T> out;
    const std::vector<T> sol =
        optdistdetail::constrained_min(obj, Aeq, beq, optdistdetail::vec(seed), lo);
    out.B1 = optdistdetail::unvec(sol, n);
    Dmap<T> b;
    b.D0 = B0;
    b.D1 = out.B1;
    out.d = dmap_dist_acf(a, b, alA, alB);
    out.global = false;
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_DMAP_OPTIM_DIST_H
