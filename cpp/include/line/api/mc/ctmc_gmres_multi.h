/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_GMRES_MULTI_H
#define LINE_API_MC_CTMC_GMRES_MULTI_H

/**
 * Restarted GMRES for a block of right-hand sides sharing one coefficient
 * matrix.
 *
 * Templated port of matlab/src/api/mc/ctmc_gmres_multi.m and the multi-column
 * overload of jar/src/main/java/jline/api/mc/Ctmc_gmres.java. Every column is
 * solved against the SAME equilibration, reordering and ILUT factorization, and
 * each column starts from the previous column's solution: this is the shape of
 * the stochastic complement, whose right-hand side is a whole block of the
 * generator, and refactorizing per column would cost more than the direct solve
 * the method replaces.
 *
 * FLAG is zero only when every column converged. On any other value the result
 * matrix is empty and the caller must fall back to the direct solve; returning
 * a partially converged block would leave the fallback ambiguous, which is the
 * reference's rule and is kept.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC, for the reason given in ctmc_gmres: the
 * iteration stops on a residual tolerance and normalizes by a Euclidean norm.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_gmres.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct GmresMultiResult {
    Matrix<T> X;  ///< solution block, empty unless flag == 0
    int flag;     ///< 0 all columns converged, otherwise the first failing flag
};

/**
 * @param A       coefficient matrix
 * @param B       right-hand sides, one per column
 * @param tol     relative residual tolerance (default 1e-12)
 * @param restart restart length; <= 0 selects min(n, 50)
 * @param maxit   outer cycles; <= 0 selects ceil(n / restart)
 */
template <class T>
GmresMultiResult<T> ctmc_gmres_multi(const Matrix<T>& A, const Matrix<T>& B, double tol = 1e-12,
                                     long restart = 0, long maxit = 0) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_gmres_multi requires transcendental arithmetic: see ctmc_gmres, the "
                  "iteration stops on a residual tolerance rather than reaching an exact value");
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("ctmc_gmres_multi: matrix is not square");
    if (B.rows() != n) throw InputError("ctmc_gmres_multi: right-hand side block has the wrong height");

    const detail::GmresPrepared<T> prep(A);
    const std::size_t nrhs = B.cols();
    Matrix<T> X(n, nrhs, num_traits<T>::from_int(0));
    std::vector<T> guess(n, num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(n)));
    std::vector<T> rhs(n);

    GmresMultiResult<T> out;
    for (std::size_t c = 0; c < nrhs; ++c) {
        for (std::size_t i = 0; i < n; ++i) rhs[i] = B(i, c);
        const GmresResult<T> r = detail::gmres_solve(prep, rhs, guess, tol, restart, maxit);
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

#endif  // LINE_API_MC_CTMC_GMRES_MULTI_H
