/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_EIG_H
#define LINE_UTIL_EIG_H

/**
 * Eigenvalues and singular values, backed by LAPACK.
 *
 * LAPACK is BSD-3, so it is the one external numerical dependency that costs
 * the port nothing in licensing terms: it links cleanly into a BSD-licensed
 * binary and, unlike GMP or MPFR, imposes no relinking obligation. It is
 * requested through the Fortran symbols directly (`dgeev_`, `dgesvd_`) rather
 * than through LAPACKE, because the reference distribution ships the library
 * without the C headers.
 *
 * DOUBLE ONLY, and deliberately so. Eigenvalues of a rational matrix are
 * algebraic numbers, not rationals: there is no exact instantiation to offer,
 * and a high-precision one would need a multiprecision QR iteration that
 * LAPACK cannot provide. Both entry points therefore take Matrix<double>,
 * and callers templated on T must convert and document the precision loss at
 * that point rather than pretending otherwise.
 *
 * When the build is configured without LAPACK the entry points throw
 * UnsupportedError naming the missing dependency, which is the same
 * refuse-by-name policy the CLI uses for unported features: a caller learns
 * what is missing instead of receiving a silently substituted answer.
 */

#include <complex>
#include <cstddef>
#include <vector>

#include "line/util/error.h"
#include "line/util/matrix.h"

#ifdef LINE_MP_HAVE_LAPACK
extern "C" {
void dgeev_(const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda,
            double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr,
            double* work, const int* lwork, int* info);
void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n, double* a,
             const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt,
             double* work, const int* lwork, int* info);
void dgees_(const char* jobvs, const char* sort, int (*select)(const double*, const double*),
            const int* n, double* a, const int* lda, int* sdim, double* wr, double* wi, double* vs,
            const int* ldvs, double* work, const int* lwork, int* bwork, int* info);
void dtrexc_(const char* compq, const int* n, double* t, const int* ldt, double* q, const int* ldq,
             int* ifst, int* ilst, double* work, int* info);
}
#endif

namespace line {

/** Eigenvalues of a general real square matrix, in LAPACK's order. */
inline std::vector<std::complex<double>> eig_values(const Matrix<double>& A) {
#ifndef LINE_MP_HAVE_LAPACK
    (void)A;
    throw UnsupportedError(
        "eig_values requires LAPACK: reconfigure with -DLINE_MP_USE_LAPACK=ON and liblapack "
        "available");
#else
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("eig_values: matrix is not square");
    if (n == 0) return {};

    // LAPACK is column-major; the input is row-major, so transpose on the way in.
    std::vector<double> a(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) a[j * n + i] = A(i, j);

    const int ni = static_cast<int>(n);
    std::vector<double> wr(n), wi(n);
    double vdummy = 0.0;
    int info = 0, lwork = -1;
    double wopt = 0.0;
    const int one = 1;
    dgeev_("N", "N", &ni, a.data(), &ni, wr.data(), wi.data(), &vdummy, &one, &vdummy, &one, &wopt,
           &lwork, &info);
    if (info != 0) throw NumericError("eig_values: LAPACK workspace query failed");
    lwork = static_cast<int>(wopt);
    std::vector<double> work(static_cast<std::size_t>(lwork));
    dgeev_("N", "N", &ni, a.data(), &ni, wr.data(), wi.data(), &vdummy, &one, &vdummy, &one,
           work.data(), &lwork, &info);
    if (info != 0) throw NumericError("eig_values: LAPACK dgeev failed to converge");

    std::vector<std::complex<double>> out(n);
    for (std::size_t i = 0; i < n; ++i) out[i] = std::complex<double>(wr[i], wi[i]);
    return out;
#endif
}

/** Largest modulus over the spectrum, i.e. the spectral radius. */
inline double spectral_radius(const Matrix<double>& A) {
    double r = 0.0;
    for (const std::complex<double>& z : eig_values(A)) {
        const double m = std::abs(z);
        if (m > r) r = m;
    }
    return r;
}

/**
 * Second largest modulus over the spectrum. This is the quantity the NCD
 * machinery needs (ctmc_courtois's epsMAX is built from the subdominant
 * eigenvalue of each diagonal block); returns 0 when the matrix is 1 x 1.
 */
inline double subdominant_modulus(const Matrix<double>& A) {
    std::vector<std::complex<double>> e = eig_values(A);
    if (e.size() < 2) return 0.0;
    double first = 0.0, second = 0.0;
    for (const std::complex<double>& z : e) {
        const double m = std::abs(z);
        if (m > first) {
            second = first;
            first = m;
        } else if (m > second) {
            second = m;
        }
    }
    return second;
}

/** Singular values in descending order. */
inline std::vector<double> svd_values(const Matrix<double>& A) {
#ifndef LINE_MP_HAVE_LAPACK
    (void)A;
    throw UnsupportedError(
        "svd_values requires LAPACK: reconfigure with -DLINE_MP_USE_LAPACK=ON and liblapack "
        "available");
#else
    const std::size_t m = A.rows(), n = A.cols();
    if (m == 0 || n == 0) return {};
    std::vector<double> a(m * n);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j) a[j * m + i] = A(i, j);

    const int mi = static_cast<int>(m), nj = static_cast<int>(n);
    const std::size_t k = m < n ? m : n;
    std::vector<double> s(k);
    double udummy = 0.0, vtdummy = 0.0;
    int info = 0, lwork = -1;
    double wopt = 0.0;
    const int one = 1;
    dgesvd_("N", "N", &mi, &nj, a.data(), &mi, s.data(), &udummy, &one, &vtdummy, &one, &wopt,
            &lwork, &info);
    if (info != 0) throw NumericError("svd_values: LAPACK workspace query failed");
    lwork = static_cast<int>(wopt);
    std::vector<double> work(static_cast<std::size_t>(lwork));
    dgesvd_("N", "N", &mi, &nj, a.data(), &mi, s.data(), &udummy, &one, &vtdummy, &one, work.data(),
            &lwork, &info);
    if (info != 0) throw NumericError("svd_values: LAPACK dgesvd failed to converge");
    return s;  // dgesvd returns them descending
#endif
}

// ---------------------------------------------------------------------------
// Real Schur form, and reordering its diagonal blocks
// ---------------------------------------------------------------------------

/**
 * Real Schur factorization A = Z T Z^T, with Z orthogonal and T upper
 * quasi-triangular: 1 x 1 diagonal blocks for real eigenvalues and 2 x 2 blocks
 * for complex conjugate pairs.
 */
struct RealSchur {
    Matrix<double> Z;  ///< orthogonal Schur vectors
    Matrix<double> T;  ///< upper quasi-triangular factor
};

/**
 * Real Schur factorization of a general square matrix (LAPACK dgees, unsorted).
 *
 * The Schur form is the right basis for splitting a spectrum into invariant
 * subspaces, which is what a multi-regime fluid queue needs: the eigenvector
 * basis exists only for a diagonalizable matrix and is complex whenever the
 * spectrum is, while Z is orthogonal and real for every real A.
 */
inline RealSchur schur_decomposition(const Matrix<double>& A) {
#ifndef LINE_MP_HAVE_LAPACK
    (void)A;
    throw UnsupportedError(
        "schur_decomposition requires LAPACK: reconfigure with -DLINE_MP_USE_LAPACK=ON and "
        "liblapack available");
#else
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("schur_decomposition: matrix is not square");
    RealSchur out;
    out.Z = Matrix<double>(n, n, 0.0);
    out.T = Matrix<double>(n, n, 0.0);
    if (n == 0) return out;

    std::vector<double> a(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) a[j * n + i] = A(i, j);

    const int ni = static_cast<int>(n);
    std::vector<double> wr(n), wi(n), vs(n * n);
    int sdim = 0, info = 0, lwork = -1;
    double wopt = 0.0;
    dgees_("V", "N", nullptr, &ni, a.data(), &ni, &sdim, wr.data(), wi.data(), vs.data(), &ni,
           &wopt, &lwork, nullptr, &info);
    if (info != 0) throw NumericError("schur_decomposition: LAPACK workspace query failed");
    lwork = static_cast<int>(wopt);
    std::vector<double> work(static_cast<std::size_t>(lwork));
    dgees_("V", "N", nullptr, &ni, a.data(), &ni, &sdim, wr.data(), wi.data(), vs.data(), &ni,
           work.data(), &lwork, nullptr, &info);
    if (info != 0) throw NumericError("schur_decomposition: LAPACK dgees failed to converge");

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            out.T(i, j) = a[j * n + i];
            out.Z(i, j) = vs[j * n + i];
        }
    return out;
#endif
}

/**
 * Reorder the diagonal blocks of a real Schur form into DESCENDING key order,
 * stably, updating Z so that A = Z T Z^T still holds.
 *
 * This is MATLAB's ordschur with a CLUSTER-NUMBER select vector, whose
 * documented behaviour is that clusters appear in descending order of the
 * cluster number. The key is given per diagonal ENTRY; for a 2 x 2 block the
 * key of its first row is used and the second is ignored, which is the same
 * requirement MATLAB imposes (a select vector must be constant on a block, and
 * a caller that splits one is asking for a factorization that does not exist
 * over the reals).
 *
 * WHY dtrexc AND NOT dtrsen. dtrsen splits a spectrum into TWO clusters, the
 * selected one and the rest, so an ordering into three or more classes needs it
 * applied recursively to trailing submatrices, with the accumulated Q and the
 * shifting block boundaries tracked by hand at every level. dtrexc moves ONE
 * diagonal block from position ifst to position ilst and updates Q itself, so
 * an arbitrary key ordering is a stable selection sort over blocks with no
 * submatrix bookkeeping at all. It also refuses to split a 2 x 2 block rather
 * than silently producing a complex pair astride a boundary, which is the
 * failure that would otherwise be discovered downstream as a complex "real"
 * subspace.
 *
 * @param s   a real Schur factorization, typically from schur_decomposition
 * @param key one value per diagonal entry; blocks are sorted descending by it
 */
inline RealSchur schur_reorder(const RealSchur& s, const std::vector<double>& key) {
#ifndef LINE_MP_HAVE_LAPACK
    (void)s;
    (void)key;
    throw UnsupportedError(
        "schur_reorder requires LAPACK: reconfigure with -DLINE_MP_USE_LAPACK=ON and liblapack "
        "available");
#else
    const std::size_t n = s.T.rows();
    if (s.T.cols() != n || s.Z.rows() != n || s.Z.cols() != n)
        throw InputError("schur_reorder: the factors are not square or not conformable");
    if (key.size() != n) throw InputError("schur_reorder: one key per diagonal entry is required");
    RealSchur out;
    out.T = s.T;
    out.Z = s.Z;
    if (n <= 1) return out;

    // Column-major working copies for LAPACK.
    std::vector<double> t(n * n), q(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            t[j * n + i] = s.T(i, j);
            q[j * n + i] = s.Z(i, j);
        }
    // The key travels with its block, so it is permuted alongside T.
    std::vector<double> k = key;

    const int ni = static_cast<int>(n);
    std::vector<double> work(n);
    const double tiny = 0.0;  // dgees leaves an exact zero below a 1 x 1 block

    // Stable selection sort over blocks. pos is the first row of the region
    // still to be ordered.
    std::size_t pos = 0;
    while (pos < n) {
        // Enumerate the blocks of the remaining region and find the first one
        // carrying the largest key; "first" keeps the sort stable.
        std::size_t best = pos;
        double bestkey = k[pos];
        std::size_t i = pos;
        while (i < n) {
            const std::size_t sz = (i + 1 < n && t[i * n + (i + 1)] != tiny) ? 2u : 1u;
            if (k[i] > bestkey) {
                bestkey = k[i];
                best = i;
            }
            i += sz;
        }
        const std::size_t bsz =
            (best + 1 < n && t[best * n + (best + 1)] != tiny) ? 2u : 1u;
        if (best != pos) {
            int ifst = static_cast<int>(best) + 1;  // LAPACK is 1-based
            int ilst = static_cast<int>(pos) + 1;
            int info = 0;
            dtrexc_("V", &ni, t.data(), &ni, q.data(), &ni, &ifst, &ilst, work.data(), &info);
            if (info == 1)
                throw NumericError(
                    "schur_reorder: dtrexc could not separate two eigenvalues that are too close; "
                    "the requested ordering splits a 2 x 2 block");
            if (info != 0) throw NumericError("schur_reorder: LAPACK dtrexc failed");
            // Move the key with the block: erase it from its old position and
            // reinsert it at the new one, exactly as dtrexc permuted the rows.
            const std::vector<double> moved(k.begin() + static_cast<long>(best),
                                            k.begin() + static_cast<long>(best + bsz));
            k.erase(k.begin() + static_cast<long>(best),
                    k.begin() + static_cast<long>(best + bsz));
            k.insert(k.begin() + static_cast<long>(pos), moved.begin(), moved.end());
        }
        pos += bsz;
    }

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            out.T(i, j) = t[j * n + i];
            out.Z(i, j) = q[j * n + i];
        }
    return out;
#endif
}

/** Numerical rank at the standard max(m,n) eps sigma_1 threshold. */
inline std::size_t matrix_rank(const Matrix<double>& A) {
    const std::vector<double> s = svd_values(A);
    if (s.empty()) return 0;
    const std::size_t dim = A.rows() > A.cols() ? A.rows() : A.cols();
    const double tol = static_cast<double>(dim) * 2.220446049250313e-16 * s[0];
    std::size_t r = 0;
    for (double v : s)
        if (v > tol) ++r;
    return r;
}

}  // namespace line

#endif  // LINE_UTIL_EIG_H
