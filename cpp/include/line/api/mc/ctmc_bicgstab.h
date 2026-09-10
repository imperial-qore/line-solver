/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_BICGSTAB_H
#define LINE_API_MC_CTMC_BICGSTAB_H

/**
 * Preconditioned stabilized biconjugate gradients, for the linear systems a
 * generator produces.
 *
 * Templated port of matlab/src/api/mc/ctmc_bicgstab.m and
 * jar/src/main/java/jline/api/mc/Ctmc_bicgstab.java. This is the
 * short-recurrence counterpart of ctmc_gmres: work and storage per iteration are
 * constant rather than growing with the Krylov dimension, so the method does not
 * restart and does not lose the optimality that restarting costs GMRES. Where
 * GMRES(m) stagnates because the useful subspace is wider than m, this
 * converges; where it does not, GMRES(m) is the more robust of the two, hence
 * the order in which ctmc_solve tries them.
 *
 * The equilibration, the reverse Cuthill-McKee reordering and the ILUT
 * preconditioner are those of ctmc_gmres, reused through detail::GmresPrepared
 * rather than reimplemented, so both methods factorize the same matrix in the
 * same order and a switch between them cannot move a reported metric for a
 * reason other than the iteration itself.
 *
 * The preconditioner is applied on the RIGHT, on the search directions p and s,
 * so the recurrence carries the residual of the ORIGINAL system and the
 * convergence test needs no unpreconditioning. This matches the choice made in
 * ctmc_gmres and, as there, means the reported relres is the true relative
 * residual rather than MATLAB's preconditioned one.
 *
 * FLAG follows the MATLAB bicgstab convention: 0 converged, 1 iteration limit,
 * 3 stagnation or divergence, 4 a scalar quantity became too small or too large
 * to continue. ITER counts matrix-vector products with A: two per complete
 * iteration, and one when the iteration converges at its half step, so an ODD
 * count is normal. Counting products rather than iterations is what makes it
 * comparable with the ITER of ctmc_gmres and across the four codebases.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC, for the reason given in ctmc_gmres: the
 * iteration stops on a residual tolerance and normalizes by a Euclidean norm,
 * so at Rational it would run to the iteration limit with exploding
 * denominators. The exact solve of the same system is ctmc_solve.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_gmres.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct BicgstabResult {
    std::vector<T> x;  ///< solution
    int flag;          ///< 0 converged, 1 iteration limit, 3 stagnation, 4 breakdown
    T relres;          ///< true relative residual norm(b - A x) / norm(b)
    long iter;         ///< matrix-vector products with A
};

template <class T>
struct BicgstabMultiResult {
    Matrix<T> X;  ///< solution block, empty unless flag == 0
    int flag;     ///< 0 all columns converged, otherwise the first failing flag
};

namespace detail {

constexpr double BICGSTAB_DEFAULT_TOL = 1e-12;
/** Complete iterations allowed by default; storage is O(n) whatever the count. */
constexpr long BICGSTAB_DEFAULT_MAXIT = 200;
/** Threshold below which rho or omega is treated as a Lanczos breakdown. */
constexpr double BICGSTAB_BREAKDOWN_TOL = 1e-14;

/** Right-preconditioned BiCGSTAB on an already prepared system. */
template <class T>
BicgstabResult<T> bicgstab_solve(const GmresPrepared<T>& prep, const std::vector<T>& rhsIn,
                                 const std::vector<T>& x0In, double tol, long maxit) {
    const std::size_t n = prep.n;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (tol <= 0.0) tol = BICGSTAB_DEFAULT_TOL;
    if (maxit <= 0) maxit = std::min(static_cast<long>(n), BICGSTAB_DEFAULT_MAXIT);
    maxit = std::max(1L, std::min(maxit, static_cast<long>(n)));
    const T tolT = num_traits<T>::from_double(tol);
    const T breakT = num_traits<T>::from_double(BICGSTAB_BREAKDOWN_TOL);

    // Row scaling, then RCM permutation, on both the right-hand side and the
    // initial guess, exactly as gmres_solve does.
    std::vector<T> rhs(n), x(n);
    for (std::size_t i = 0; i < n; ++i) rhs[i] = rhsIn[prep.perm[i]] / prep.rowScale[prep.perm[i]];
    for (std::size_t i = 0; i < n; ++i) x[i] = x0In[prep.perm[i]];

    T bnorm = vec_norm2(rhs);
    if (bnorm == zero) bnorm = one;

    BicgstabResult<T> out;
    out.flag = 1;
    out.iter = 0;

    std::vector<T> r(n);
    prep.csr.mult(x, r);
    for (std::size_t i = 0; i < n; ++i) r[i] = rhs[i] - r[i];
    out.relres = vec_norm2(r) / bnorm;
    if (out.relres <= tolT) {
        out.x.assign(n, zero);
        for (std::size_t i = 0; i < n; ++i) out.x[prep.perm[i]] = x[i];
        out.flag = 0;
        return out;
    }

    // The shadow residual is fixed at the initial residual, the standard choice:
    // any vector not orthogonal to r would do, and this one cannot be.
    const std::vector<T> rhat = r;
    std::vector<T> p(n, zero), v(n, zero), s(n, zero), t(n, zero), ph(n, zero), sh(n, zero);

    T rho = one, alpha = one, omega = one;
    T bestrelres = out.relres;

    for (long it = 0; it < maxit; ++it) {
        const T rhoNew = vec_dot(rhat, r);
        // rho vanishing is the biorthogonality breakdown of the underlying
        // Lanczos process, not slow convergence: restarting with a fresh shadow
        // vector would discard the iterate, so the caller is told to use another
        // method instead.
        if (num_abs(rhoNew) <= breakT * vec_norm2(rhat) * vec_norm2(r)) {
            out.flag = 4;
            break;
        }
        if (it == 0) {
            p = r;
        } else {
            if (omega == zero) {
                out.flag = 4;
                break;
            }
            const T beta = (rhoNew / rho) * (alpha / omega);
            for (std::size_t i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        rho = rhoNew;

        prep.M.apply(p, ph);
        prep.csr.mult(ph, v);
        ++out.iter;

        const T rhatv = vec_dot(rhat, v);
        if (rhatv == zero || !num_isfinite(rhatv)) {
            out.flag = 4;
            break;
        }
        alpha = rho / rhatv;

        for (std::size_t i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];

        // Half-step convergence: s is the residual of x + alpha*ph, so a
        // converged s reaches the answer without the second matvec.
        const T snorm = vec_norm2(s);
        if (snorm / bnorm <= tolT) {
            for (std::size_t i = 0; i < n; ++i) x[i] += alpha * ph[i];
            out.relres = snorm / bnorm;
            out.flag = 0;
            break;
        }

        prep.M.apply(s, sh);
        prep.csr.mult(sh, t);
        ++out.iter;

        const T tt = vec_dot(t, t);
        if (tt == zero || !num_isfinite(tt)) {
            out.flag = 4;
            break;
        }
        omega = vec_dot(t, s) / tt;

        for (std::size_t i = 0; i < n; ++i) x[i] += alpha * ph[i] + omega * sh[i];
        for (std::size_t i = 0; i < n; ++i) r[i] = s[i] - omega * t[i];

        out.relres = vec_norm2(r) / bnorm;
        if (out.relres <= tolT) {
            out.flag = 0;
            break;
        }
        // omega vanishing stalls the update of x while leaving r finite, so the
        // iteration would spin without progress.
        if (num_abs(omega) <= breakT) {
            out.flag = 4;
            break;
        }
        // BiCGSTAB residuals are non-monotone by construction, so an increase is
        // not by itself stagnation and the test is against the BEST residual
        // seen rather than the previous one. Growing two orders of magnitude
        // past that best is divergence.
        if (out.relres > num_traits<T>::from_int(100) * bestrelres) {
            out.flag = 3;
            break;
        }
        if (out.relres < bestrelres) bestrelres = out.relres;
    }

    out.x.assign(n, zero);
    for (std::size_t i = 0; i < n; ++i) out.x[prep.perm[i]] = x[i];
    for (std::size_t i = 0; i < n; ++i)
        if (!num_isfinite(out.x[i])) {
            out.flag = 4;
            out.relres = num_traits<T>::from_double(1e300);
            return out;
        }
    if (out.relres <= tolT) out.flag = 0;
    return out;
}

}  // namespace detail

/**
 * @param A     coefficient matrix, already assembled
 * @param b     right-hand side
 * @param tol   relative residual tolerance (default 1e-12)
 * @param maxit complete iterations; <= 0 selects min(n, 200)
 * @param x0    initial guess; empty selects the uniform vector ones(n)/n
 */
template <class T>
BicgstabResult<T> ctmc_bicgstab(const Matrix<T>& A, const std::vector<T>& b, double tol = 1e-12,
                                long maxit = 0, const std::vector<T>& x0 = std::vector<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_bicgstab requires transcendental arithmetic: the iteration stops on a "
                  "residual tolerance and normalizes by a Euclidean norm, so there is no exact "
                  "result to converge to; use ctmc_solve for an exact solve");
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("ctmc_bicgstab: matrix is not square");
    if (b.size() != n) throw InputError("ctmc_bicgstab: right-hand side has the wrong length");
    if (!x0.empty() && x0.size() != n) throw InputError("ctmc_bicgstab: initial guess has the wrong length");
    std::vector<T> guess = x0;
    if (guess.empty())
        guess.assign(n, num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(n)));
    const detail::GmresPrepared<T> prep(A);
    return detail::bicgstab_solve(prep, b, guess, tol, maxit);
}

/**
 * Every column of B solved against the SAME equilibration, reordering and ILUT
 * factorization, each column starting from the previous column's solution. FLAG
 * is zero only when every column converged; on any other value the result matrix
 * is empty and the caller must fall back, a partially converged block leaving
 * the fallback ambiguous.
 *
 * @param A     coefficient matrix
 * @param B     right-hand sides, one per column
 * @param tol   relative residual tolerance (default 1e-12)
 * @param maxit complete iterations per column; <= 0 selects min(n, 200)
 */
template <class T>
BicgstabMultiResult<T> ctmc_bicgstab_multi(const Matrix<T>& A, const Matrix<T>& B, double tol = 1e-12,
                                           long maxit = 0) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_bicgstab_multi requires transcendental arithmetic: see ctmc_bicgstab, the "
                  "iteration stops on a residual tolerance rather than reaching an exact value");
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("ctmc_bicgstab_multi: matrix is not square");
    if (B.rows() != n) throw InputError("ctmc_bicgstab_multi: right-hand side block has the wrong height");

    const detail::GmresPrepared<T> prep(A);
    const std::size_t nrhs = B.cols();
    Matrix<T> X(n, nrhs, num_traits<T>::from_int(0));
    std::vector<T> guess(n, num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(n)));
    std::vector<T> rhs(n);

    BicgstabMultiResult<T> out;
    for (std::size_t c = 0; c < nrhs; ++c) {
        for (std::size_t i = 0; i < n; ++i) rhs[i] = B(i, c);
        const BicgstabResult<T> r = detail::bicgstab_solve(prep, rhs, guess, tol, maxit);
        if (r.flag != 0) {
            out.X = Matrix<T>();
            out.flag = r.flag;
            return out;
        }
        for (std::size_t i = 0; i < n; ++i) X(i, c) = r.x[i];
        guess = r.x;
    }
    out.X = X;
    out.flag = 0;
    return out;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_BICGSTAB_H
