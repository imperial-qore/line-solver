/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_SVD_H
#define LINE_UTIL_SVD_H

/**
 * Singular value decomposition WITH the singular vectors, and the
 * Moore-Penrose pseudo-inverse built from it.
 *
 * util/eig.h already exposes the singular VALUES (dgesvd with jobu = jobvt =
 * 'N'), which is all the rank and conditioning tests need. This header adds
 * the factors themselves, for the one place in the port that needs them: the
 * pinv fallback of ldqbd_R, taken when a level's local block is singular.
 *
 * DOUBLE ONLY, for the same reason eig.h is: singular values of a rational
 * matrix are algebraic, not rational, so there is no exact instantiation to
 * offer, and LAPACK has no multiprecision path. Callers templated on T must
 * convert and state the precision loss at the call site. Configured without
 * LAPACK the entry points throw UnsupportedError naming the dependency rather
 * than substituting anything.
 *
 * The Fortran symbol dgesvd_ is declared by eig.h under the same LAPACK guard
 * and is reused from there, so there is exactly one declaration of it in the
 * tree.
 */

#include <cstddef>
#include <vector>

#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {

/** A = U diag(s) Vt, with U (m x m), s of length min(m,n) and Vt (n x n). */
struct SvdFactors {
    Matrix<double> U;
    std::vector<double> s;
    Matrix<double> Vt;
};

/** Full SVD of a real matrix, singular values in descending order. */
inline SvdFactors svd_full(const Matrix<double>& A) {
#ifndef LINE_MP_HAVE_LAPACK
    (void)A;
    throw UnsupportedError(
        "svd_full requires LAPACK: reconfigure with -DLINE_MP_USE_LAPACK=ON and liblapack "
        "available");
#else
    const std::size_t m = A.rows(), n = A.cols();
    if (m == 0 || n == 0) throw InputError("svd_full: empty matrix");

    // LAPACK is column-major, the port is row-major, so transpose on the way in.
    std::vector<double> a(m * n);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j) a[j * m + i] = A(i, j);

    const int mi = static_cast<int>(m), nj = static_cast<int>(n);
    const std::size_t k = m < n ? m : n;
    std::vector<double> s(k), u(m * m), vt(n * n);
    int info = 0, lwork = -1;
    double wopt = 0.0;
    dgesvd_("A", "A", &mi, &nj, a.data(), &mi, s.data(), u.data(), &mi, vt.data(), &nj, &wopt,
            &lwork, &info);
    if (info != 0) throw NumericError("svd_full: LAPACK workspace query failed");
    lwork = static_cast<int>(wopt);
    std::vector<double> work(static_cast<std::size_t>(lwork));
    dgesvd_("A", "A", &mi, &nj, a.data(), &mi, s.data(), u.data(), &mi, vt.data(), &nj, work.data(),
            &lwork, &info);
    if (info != 0) throw NumericError("svd_full: LAPACK dgesvd failed to converge");

    SvdFactors out;
    out.s = s;
    out.U = Matrix<double>(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) out.U(i, j) = u[j * m + i];
    out.Vt = Matrix<double>(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) out.Vt(i, j) = vt[j * n + i];
    return out;
#endif
}

/**
 * Moore-Penrose pseudo-inverse, A^+ = V diag(1/s_i) U^T over the singular
 * values above max(m,n) eps sigma_1, which is MATLAB's default pinv tolerance.
 */
inline Matrix<double> pinv(const Matrix<double>& A) {
    const SvdFactors f = svd_full(A);
    const std::size_t m = A.rows(), n = A.cols();
    const std::size_t k = f.s.size();
    const std::size_t dim = m > n ? m : n;
    const double tol = static_cast<double>(dim) * 2.220446049250313e-16 * (k ? f.s[0] : 0.0);
    Matrix<double> X(n, m, 0.0);
    for (std::size_t r = 0; r < k; ++r) {
        if (!(f.s[r] > tol)) continue;
        const double inv = 1.0 / f.s[r];
        // X += inv * V(:,r) * U(:,r)^T, with V(:,r) the r-th ROW of Vt.
        for (std::size_t i = 0; i < n; ++i) {
            const double vi = f.Vt(r, i);
            if (vi == 0.0) continue;
            for (std::size_t j = 0; j < m; ++j) X(i, j) += inv * vi * f.U(j, r);
        }
    }
    return X;
}

}  // namespace line

#endif  // LINE_UTIL_SVD_H
