/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_SYLVESTER_H
#define LINE_UTIL_SYLVESTER_H

/**
 * The Sylvester equation A X + X B = C, and MATLAB's `lyap(A,B,C)`.
 *
 * Solved in the Kronecker form (I_m kron A + B^T kron I_n) vec(X) = vec(C) with
 * `vec` column-major, so the whole thing is one LU solve of order n*m in `T`.
 * That is the ONLY formulation available at every arithmetic the port offers:
 * Bartels-Stewart needs a real Schur factorization, `util/eig.h` is
 * deliberately double-only (eigenvalues of a rational matrix are algebraic, not
 * rational), so a Schur-based solver would silently pin every caller to double.
 * The cost is O((n m)^3) against Bartels-Stewart's O(n^3 + m^3); the callers in
 * this port pass phase-space matrices whose order is the product of an arrival
 * MMAP order and a total service-PH order, both small, and the whole point of
 * the header is that `Rational` and `Real<D>` get an answer at all.
 *
 * `mfq_multiregime.h` carries the identical Kronecker construction for
 * `Matrix<double>`; it predates this header and stays where it is because it is
 * used only from double-only code and moving it would churn two large files
 * that no test covers at another arithmetic.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

#ifdef LINE_MP_HAVE_LAPACK
extern "C" {
/** Quasi-triangular Sylvester solve, the core of Bartels-Stewart. */
void dtrsyl_(const char* trana, const char* tranb, const int* isgn, const int* m, const int* n,
             const double* a, const int* lda, const double* b, const int* ldb, double* c,
             const int* ldc, double* scale, int* info);
}
#endif

namespace line {

/**
 * The Kronecker operator of a FIXED (A,B) pair, factorized once.
 *
 * The MMAP[K]/PH[K]/1 queue-length recursion solves a chain of Sylvester
 * equations that share A and B and differ only in the right-hand side, one per
 * queue-length level. Refactorizing per level would make the recursion cubic in
 * the level count for no reason.
 */
template <class T>
class SylvesterFactor {
public:
    SylvesterFactor(const Matrix<T>& A, const Matrix<T>& B) : n_(A.rows()), m_(B.rows()) {
        const T zero = num_traits<T>::from_int(0);
        if (A.cols() != n_ || B.cols() != m_)
            throw InputError("SylvesterFactor: A and B must be square");
        if (n_ == 0 || m_ == 0) return;
        LU_ = Matrix<T>(n_ * m_, n_ * m_, zero);
        for (std::size_t j = 0; j < m_; ++j)
            for (std::size_t i = 0; i < n_; ++i) {
                const std::size_t row = j * n_ + i;
                for (std::size_t k = 0; k < n_; ++k) LU_(row, j * n_ + k) += A(i, k);
                for (std::size_t k = 0; k < m_; ++k) LU_(row, k * n_ + i) += B(k, j);
            }
        piv_ = lu_factor(LU_);
    }

    /** Solve A X + X B = C. */
    Matrix<T> solve_sylvester(const Matrix<T>& C) const {
        const T zero = num_traits<T>::from_int(0);
        Matrix<T> X(n_, m_, zero);
        if (n_ == 0 || m_ == 0) return X;
        if (C.rows() != n_ || C.cols() != m_)
            throw InputError("SylvesterFactor: the right-hand side is not conformable");
        std::vector<T> rhs(n_ * m_, zero);
        for (std::size_t j = 0; j < m_; ++j)
            for (std::size_t i = 0; i < n_; ++i) rhs[j * n_ + i] = C(i, j);
        lu_solve(LU_, piv_, rhs);
        for (std::size_t j = 0; j < m_; ++j)
            for (std::size_t i = 0; i < n_; ++i) X(i, j) = rhs[j * n_ + i];
        return X;
    }

    /** Solve A X + X B + C = 0, MATLAB's lyap(A,B,C). */
    Matrix<T> solve_lyap(const Matrix<T>& C) const {
        Matrix<T> negC(C.rows(), C.cols());
        for (std::size_t i = 0; i < C.rows(); ++i)
            for (std::size_t j = 0; j < C.cols(); ++j) negC(i, j) = -C(i, j);
        return solve_sylvester(negC);
    }

private:
    std::size_t n_, m_;
    Matrix<T> LU_;
    std::vector<std::size_t> piv_;
};

/**
 * Solve A X + X B = C for X.
 *
 * @param A n x n
 * @param B m x m
 * @param C n x m
 */
template <class T>
Matrix<T> sylvester_solve(const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C) {
    if (C.rows() != A.rows() || C.cols() != B.rows())
        throw InputError("sylvester_solve: the blocks are not conformable");
    return SylvesterFactor<T>(A, B).solve_sylvester(C);
}

/** MATLAB `lyap(A,B,C)` solves A X + X B + C = 0, i.e. sylvester_solve(A, B, -C). */
template <class T>
Matrix<T> lyap_solve(const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C) {
    if (C.rows() != A.rows() || C.cols() != B.rows())
        throw InputError("lyap_solve: the blocks are not conformable");
    return SylvesterFactor<T>(A, B).solve_lyap(C);
}

/**
 * A X + X B = C by Bartels-Stewart, at double.
 *
 * WHY THIS EXISTS BESIDE THE KRONECKER SOLVER ABOVE. The Kronecker form is the
 * only one available in an arbitrary field, and its O((n m)^3) cost is
 * acceptable for the phase-space matrices the MMAP queue recursion passes. The
 * FJ_codes fork-join engine passes matrices of order (C + 1) * m^2 * ma with C
 * defaulting to 100, i.e. n and m in the hundreds to low thousands: the
 * Kronecker operator of an 800 x 800 pair is 640000 square, which is 3 TB of
 * LU. Bartels-Stewart is O(n^3 + m^3) and factorizes the same problem in
 * fractions of a second, at the price of a real Schur factorization -- LAPACK,
 * hence double only, exactly as `util/eig.h` is.
 *
 * A X + X B = C is solved as Ta Y + Y Tb = Za' C Zb with A = Za Ta Za' and
 * B = Zb Tb Zb' real Schur, the quasi-triangular core done by dtrsyl, and
 * X = Za Y Zb'. dtrsyl's `scale` guards against overflow when the two spectra
 * nearly collide; it is divided out here, and a scale of zero means the
 * equation has no solution, which is reported by name rather than returned as
 * an infinity.
 */
inline Matrix<double> sylvester_schur(const Matrix<double>& A, const Matrix<double>& B,
                                      const Matrix<double>& C) {
#ifndef LINE_MP_HAVE_LAPACK
    (void)A;
    (void)B;
    (void)C;
    throw UnsupportedError(
        "sylvester_schur requires LAPACK: reconfigure with -DLINE_MP_USE_LAPACK=ON and liblapack "
        "available");
#else
    const std::size_t n = A.rows(), m = B.rows();
    if (A.cols() != n || B.cols() != m)
        throw InputError("sylvester_schur: A and B must be square");
    if (C.rows() != n || C.cols() != m)
        throw InputError("sylvester_schur: the right-hand side is not conformable");
    Matrix<double> X(n, m, 0.0);
    if (n == 0 || m == 0) return X;

    const RealSchur sa = schur_decomposition(A);
    const RealSchur sb = schur_decomposition(B);
    // Ct = Za' C Zb, column-major for LAPACK.
    const Matrix<double> Ct = matmul(matmul(sa.Z.transpose(), C), sb.Z);
    std::vector<double> ta(n * n), tb(m * m), c(n * m);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) ta[j * n + i] = sa.T(i, j);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) tb[j * m + i] = sb.T(i, j);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < m; ++j) c[j * n + i] = Ct(i, j);

    const int ni = static_cast<int>(n), mi = static_cast<int>(m), isgn = 1;
    double scale = 0.0;
    int info = 0;
    dtrsyl_("N", "N", &isgn, &ni, &mi, ta.data(), &ni, tb.data(), &mi, c.data(), &ni, &scale,
            &info);
    if (info < 0) throw InputError("sylvester_schur: LAPACK dtrsyl rejected an argument");
    if (info == 1)
        throw NumericError(
            "sylvester_schur: A and -B have common or very close eigenvalues, so the Sylvester "
            "equation is singular or nearly so and its solution is not determined");
    if (scale == 0.0)
        throw NumericError("sylvester_schur: LAPACK dtrsyl returned a zero scale factor");

    Matrix<double> Y(n, m, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < m; ++j) Y(i, j) = c[j * n + i] / scale;
    return matmul(matmul(sa.Z, Y), sb.Z.transpose());
#endif
}

/** MATLAB `lyap(A,B,C)` at double via Bartels-Stewart: A X + X B + C = 0. */
inline Matrix<double> lyap_schur(const Matrix<double>& A, const Matrix<double>& B,
                                 const Matrix<double>& C) {
    Matrix<double> negC(C.rows(), C.cols());
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < C.cols(); ++j) negC(i, j) = -C(i, j);
    return sylvester_schur(A, B, negC);
}

}  // namespace line

#endif  // LINE_UTIL_SYLVESTER_H
