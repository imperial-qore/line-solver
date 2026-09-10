/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_OPTIM_DIST_H
#define LINE_API_MAM_MAP_OPTIM_DIST_H

/**
 * Fit the D1 of a MAP by MINIMIZING a distance to a reference MAP.
 *
 * Port of matlab/lib/kpctoolbox/map/map_optim_dist.m and map_optim_dist_acf.m
 * (twins in python/line_solver/api/mapdist/continuous.py). The distances
 * themselves -- `map_dist`, `map_dist_acf` -- are already in `map_dist.h`; what
 * is here is the OPTIMIZATION over D1 with D0 held fixed.
 *
 * Reference: G. Horvath, "Measuring the distance between MAPs and some
 * applications", ASMTA 2015, LNCS 9081, pp. 95-109.
 *
 * WHAT IS BEING FITTED, AND WHY D0 IS FIXED. Given a reference MAP A and a
 * chosen D0 for the approximation B, the free parameter is B's D1. The
 * constraints are exactly the two that make (B0, B1) a MAP with the DECLARED
 * embedded distribution alB:
 *
 *   - `alB (-B0)^-1 B1 = alB`, i.e. alB is stationary at arrivals. This is the
 *     `kron(I, alB inv(-B0))` block;
 *   - each row of B0 + B1 sums to zero, i.e. B1's row sums are `-B0`'s. This is
 *     the `kron(ones(1,N), I)` block.
 *
 * Every entry is bounded below by 1e-6 rather than by 0, which is the
 * reference's own floor: an exactly zero entry makes the fitted MAP reducible,
 * and the distance is then defined on a different state space than the one the
 * constraints were written for.
 *
 * THE LAG-1 CASE IS A CONVEX QP AND IS SOLVED AS ONE. At L = 1 the distance is
 * a quadratic form in vec(B1),
 *
 *     d(vB1) = vB1' H vB1 + vA1' H_AA vA1 - 2 vA1' H_AB vB1,
 *
 * with `H = mkron(X_BB, Z_BB)` positive semidefinite, so the reference calls
 * `quadprog` and gets a global optimum. Every other case is a general nonlinear
 * program over the same affine set and the reference calls `fmincon`, which is
 * local. Both are served here by `auglag`, and the DISTINCTION IS PRESERVED IN
 * WHAT THE RESULT PROMISES: `global` is set only on the L = 1 path.
 *
 * The six Lyapunov solves are the reference's; `lyap_solve` is the tree's
 * MATLAB-compatible `lyap(A,B,C)`, solving `A X + X B + C = 0`.
 *
 * COLUMN-MAJOR THROUGHOUT. `vec` here is MATLAB's, stacking COLUMNS, because
 * the Kronecker identity the quadratic form rests on (`vec(A X B) = kron(B', A)
 * vec(X)`) holds in that convention and in no other. Getting this wrong
 * transposes the fitted D1 without changing its row sums, so the constraints
 * still pass and only the distance is wrong.
 *
 * ARITHMETIC: transcendental, inherited from the distances.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_dist.h"
#include "line/api/mam/map_moment.h"
#include "line/util/lu.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/sylvester.h"

namespace line {
namespace mam {

/** What an optimizing distance fit returns. */
template <class T>
struct MapOptimDist {
    Matrix<T> B1;          ///< the fitted D1
    T d;                   ///< the distance achieved
    bool global = false;   ///< true only on the convex lag-1 path
};

namespace optdistdetail {

/** MATLAB's `vec`: the COLUMNS stacked. */
template <class T>
std::vector<T> vec(const Matrix<T>& A) {
    std::vector<T> v(A.rows() * A.cols(), num_traits<T>::from_int(0));
    std::size_t c = 0;
    for (std::size_t j = 0; j < A.cols(); ++j)
        for (std::size_t i = 0; i < A.rows(); ++i) v[c++] = A(i, j);
    return v;
}

/** The inverse of `vec` at a known shape. */
template <class T>
Matrix<T> unvec(const std::vector<T>& v, std::size_t n) {
    Matrix<T> A(n, n, num_traits<T>::from_int(0));
    std::size_t c = 0;
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i) A(i, j) = v[c++];
    return A;
}

/**
 * `kron(A, B)` in MATLAB's ordering.
 *
 * Spelled locally rather than taken from `util/` because the name is already
 * declared in the enclosing namespace and an unqualified call there is
 * ambiguous; the two agree, this one just cannot be confused for the other.
 */
template <class T>
Matrix<T> mkron(const Matrix<T>& A, const Matrix<T>& B) {
    const std::size_t ar = A.rows(), ac = A.cols(), br = B.rows(), bc = B.cols();
    Matrix<T> K(ar * br, ac * bc, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < ar; ++i)
        for (std::size_t j = 0; j < ac; ++j)
            for (std::size_t k = 0; k < br; ++k)
                for (std::size_t l = 0; l < bc; ++l)
                    K(i * br + k, j * bc + l) = A(i, j) * B(k, l);
    return K;
}

/** The row sums of -M, as a column. */
template <class T>
Matrix<T> neg_row_sums(const Matrix<T>& M) {
    Matrix<T> v(M.rows(), 1, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < M.rows(); ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < M.cols(); ++j) s += M(i, j);
        v(i, 0) = -s;
    }
    return v;
}

/**
 * The equality block: alB stationary at arrivals, then the row sums.
 *
 * Returned as (Aeq, beq) so both fitters state the SAME constraints; a fitter
 * that rebuilt them would be free to drift from the other.
 */
template <class T>
void build_constraints(const Matrix<T>& B0, const std::vector<T>& alB, Matrix<T>* Aeq,
                       std::vector<T>* beq) {
    const std::size_t n = B0.rows();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    Matrix<T> negB0(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negB0(i, j) = -B0(i, j);
    // alB * inv(-B0), by solving rather than inverting.
    Matrix<T> negB0T(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negB0T(i, j) = negB0(j, i);
    const std::vector<T> row = solve(negB0T, alB);  // (alB inv(-B0))' solves (-B0)' x = alB'

    Matrix<T> rowM(1, n, zero);
    for (std::size_t j = 0; j < n; ++j) rowM(0, j) = row[j];
    Matrix<T> I(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) I(i, i) = one;
    Matrix<T> ones1(1, n, one);

    const Matrix<T> top = mkron(I, rowM);      // (n x n^2)
    const Matrix<T> bot = mkron(ones1, I);     // (n x n^2)
    *Aeq = Matrix<T>(2 * n, n * n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n * n; ++j) {
            (*Aeq)(i, j) = top(i, j);
            (*Aeq)(n + i, j) = bot(i, j);
        }
    const Matrix<T> b = neg_row_sums(B0);
    beq->assign(2 * n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        (*beq)[i] = alB[i];
        (*beq)[n + i] = b(i, 0);
    }
}

/** Minimize `f` over {Aeq x = beq, x >= lo} from a feasible-ish start. */
template <class T, class F>
std::vector<T> constrained_min(F f, const Matrix<T>& Aeq, const std::vector<T>& beq,
                               const std::vector<T>& x0, const T& lo) {
    const std::size_t m = x0.size();
    auto h = [&Aeq, &beq, m](const std::vector<T>& x) {
        std::vector<T> r(Aeq.rows(), num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < Aeq.rows(); ++i) {
            T s = -beq[i];
            for (std::size_t j = 0; j < m; ++j) s += Aeq(i, j) * x[j];
            r[i] = s;
        }
        return r;
    };
    auto g = [](const std::vector<T>&) { return std::vector<T>(); };
    std::vector<Bound<T>> bounds(m);
    for (std::size_t j = 0; j < m; ++j) {
        bounds[j].has_lo = true;
        bounds[j].lo = lo;
    }
    return auglag(f, h, g, x0, bounds).x;
}

}  // namespace optdistdetail

/**
 * Fit B1 minimizing the lag-L joint-density distance to `a`, with B0 fixed.
 *
 * @param a   the reference MAP
 * @param alA its embedded (at-arrivals) distribution
 * @param B0  the approximation's D0, held fixed
 * @param alB the approximation's declared embedded distribution
 * @param L   number of lags; L = 1 takes the convex quadratic path
 */
template <class T>
MapOptimDist<T> map_optim_dist(const Map<T>& a, const std::vector<T>& alA, const Matrix<T>& B0,
                               const std::vector<T>& alB, unsigned L) {
    using namespace optdistdetail;
    static_assert(num_traits<T>::has_transcendental, "map_optim_dist needs the MAP distances");
    const std::size_t n = B0.rows();
    if (n == 0 || B0.cols() != n) throw InputError("map_optim_dist: B0 must be square");
    if (alB.size() != n) throw InputError("map_optim_dist: alB has the wrong length");
    if (alA.size() != a.D0.rows()) throw InputError("map_optim_dist: alA has the wrong length");
    if (L == 0) throw InputError("map_optim_dist: at least one lag is required");

    const T zero = num_traits<T>::from_int(0);
    const T lo = num_traits<T>::from_double(1e-6);
    Matrix<T> Aeq;
    std::vector<T> beq;
    build_constraints(B0, alB, &Aeq, &beq);

    MapOptimDist<T> out;
    out.d = zero;

    if (L == 1) {
        // The six Lyapunov solves of the reference, then the quadratic form.
        const std::size_t na = a.D0.rows();
        Matrix<T> A0t(na, na, zero), B0t(n, n, zero);
        for (std::size_t i = 0; i < na; ++i)
            for (std::size_t j = 0; j < na; ++j) A0t(i, j) = a.D0(j, i);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) B0t(i, j) = B0(j, i);

        Matrix<T> alAc(na, 1, zero), alBc(n, 1, zero);
        for (std::size_t i = 0; i < na; ++i) alAc(i, 0) = alA[i];
        for (std::size_t i = 0; i < n; ++i) alBc(i, 0) = alB[i];
        const Matrix<T> av = neg_row_sums(a.D0), bv = neg_row_sums(B0);

        auto outer = [zero](const Matrix<T>& u, const Matrix<T>& v) {
            Matrix<T> M(u.rows(), v.rows(), zero);
            for (std::size_t i = 0; i < u.rows(); ++i)
                for (std::size_t j = 0; j < v.rows(); ++j) M(i, j) = u(i, 0) * v(j, 0);
            return M;
        };

        const Matrix<T> Z_AB = lyap_solve(A0t, B0, outer(alAc, alBc));
        const Matrix<T> Z_AA = lyap_solve(A0t, a.D0, outer(alAc, alAc));
        const Matrix<T> Z_BB = lyap_solve(B0t, B0, outer(alBc, alBc));
        const Matrix<T> X_AB = lyap_solve(a.D0, B0t, outer(av, bv));
        const Matrix<T> X_AA = lyap_solve(a.D0, A0t, outer(av, av));
        const Matrix<T> X_BB = lyap_solve(B0, B0t, outer(bv, bv));

        const Matrix<T> H = mkron(X_BB, Z_BB);
        const Matrix<T> HAB = mkron(X_AB, Z_AB);
        const Matrix<T> HAA = mkron(X_AA, Z_AA);
        const std::vector<T> vA1 = vec(a.D1);

        // f(x) = x' H x - 2 (vA1' HAB) x, plus the constant vA1' HAA vA1.
        std::vector<T> lin(n * n, zero);
        for (std::size_t j = 0; j < n * n; ++j) {
            T s = zero;
            for (std::size_t i = 0; i < vA1.size(); ++i) s += vA1[i] * HAB(i, j);
            lin[j] = s;
        }
        T cst = zero;
        for (std::size_t i = 0; i < vA1.size(); ++i)
            for (std::size_t j = 0; j < vA1.size(); ++j) cst += vA1[i] * HAA(i, j) * vA1[j];

        auto quad = [&H, &lin, cst, zero](const std::vector<T>& x) {
            T v = cst;
            for (std::size_t i = 0; i < x.size(); ++i) {
                T row = zero;
                for (std::size_t j = 0; j < x.size(); ++j) row += H(i, j) * x[j];
                v += x[i] * row;
                v -= num_traits<T>::from_int(2) * lin[i] * x[i];
            }
            return v;
        };

        // The start point satisfies the row sums by construction, which keeps
        // the multiplier iteration from beginning far outside the affine set.
        std::vector<T> x0(n * n, lo);
        Matrix<T> seed(n, n, zero);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                seed(i, j) = bv(i, 0) / num_traits<T>::from_int(static_cast<long>(n));
        x0 = vec(seed);

        const std::vector<T> sol = constrained_min(quad, Aeq, beq, x0, lo);
        out.B1 = unvec(sol, n);
        out.d = quad(sol);
        out.global = true;  // the quadratic form is convex, so this is the optimum
        return out;
    }

    // Every other lag count: the distance itself, minimized over the same set.
    auto obj = [&a, &B0, &alA, &alB, L, n](const std::vector<T>& x) {
        Map<T> b;
        b.D0 = B0;
        b.D1 = unvec(x, n);
        return map_dist(a, b, L, alA, alB);
    };
    Matrix<T> seed(n, n, zero);
    const Matrix<T> bv = neg_row_sums(B0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            seed(i, j) = bv(i, 0) / num_traits<T>::from_int(static_cast<long>(n));
    const std::vector<T> sol = constrained_min(obj, Aeq, beq, vec(seed), lo);
    out.B1 = unvec(sol, n);
    out.d = obj(sol);
    out.global = false;  // a general nonlinear program: this is a local optimum
    return out;
}

/**
 * Fit B1 minimizing the AUTOCORRELATION distance to `a`, with B0 fixed.
 *
 * The reference re-evaluates the distance at the returned B1 rather than
 * trusting the optimizer's own objective value, and so does this: the two can
 * differ when the solve stops on its own tolerance, and the number a caller
 * reports should be the distance of the MAP it was handed.
 */
template <class T>
MapOptimDist<T> map_optim_dist_acf(const Map<T>& a, const std::vector<T>& alA,
                                   const Matrix<T>& B0, const std::vector<T>& alB) {
    using namespace optdistdetail;
    static_assert(num_traits<T>::has_transcendental, "map_optim_dist_acf needs the MAP distances");
    const std::size_t n = B0.rows();
    if (n == 0 || B0.cols() != n) throw InputError("map_optim_dist_acf: B0 must be square");
    if (alB.size() != n) throw InputError("map_optim_dist_acf: alB has the wrong length");
    if (alA.size() != a.D0.rows()) throw InputError("map_optim_dist_acf: alA has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    const T lo = num_traits<T>::from_double(1e-6);
    Matrix<T> Aeq;
    std::vector<T> beq;
    build_constraints(B0, alB, &Aeq, &beq);

    auto obj = [&a, &B0, &alA, &alB, n](const std::vector<T>& x) {
        Map<T> b;
        b.D0 = B0;
        b.D1 = unvec(x, n);
        return map_dist_acf(a, b, alA, alB);
    };
    Matrix<T> seed(n, n, zero);
    const Matrix<T> bv = neg_row_sums(B0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            seed(i, j) = bv(i, 0) / num_traits<T>::from_int(static_cast<long>(n));

    MapOptimDist<T> out;
    const std::vector<T> sol = constrained_min(obj, Aeq, beq, vec(seed), lo);
    out.B1 = unvec(sol, n);
    Map<T> b;
    b.D0 = B0;
    b.D1 = out.B1;
    out.d = map_dist_acf(a, b, alA, alB);
    out.global = false;
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_OPTIM_DIST_H
