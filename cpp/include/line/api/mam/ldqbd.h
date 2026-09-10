/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_LDQBD_H
#define LINE_API_MAM_LDQBD_H

/**
 * Level-dependent QBD processes with finitely many levels: the rate matrices
 * R^(n) by the backward matrix continued fraction, and the stationary
 * distribution.
 *
 * Templated port of matlab/src/api/mam/ldqbd_R.m, ldqbd_pi.m and ldqbd.m,
 * which implement Algorithms 1 and 3 of T. Phung-Duc, H. Masuyama,
 * S. Kasahara, Y. Takahashi, "A Simple Algorithm for the Rate Matrices of
 * Level-Dependent QBD Processes", QTNA 2010.
 *
 * The generator is block tridiagonal over levels 0..N,
 *
 *     Q = | Q1(0)  Q0(0)   0      ...            |
 *         | Q2(1)  Q1(1)  Q0(1)   ...            |
 *         |  0     Q2(2)  Q1(2)   ...            |
 *         | ...                   Q2(N)   Q1(N)  |
 *
 * with Q0[n] the up-block from level n, Q1[n] the local block at level n and
 * Q2[n] the down-block from level n. The dimensions may vary with the level,
 * so R^(n) is (order of level n-1) x (order of level n). The recursion runs
 * downwards from the top level,
 *
 *     R^(N) = Q0[N-1] (-Q1[N])^-1
 *     R^(n) = Q0[n-1] (-Q1[n] - R^(n+1) Q2[n+1])^-1,   n = N-1, ..., 1
 *
 * and the stationary vectors follow from the level-0 balance
 * pi_0 (Q1[0] + R^(1) Q2[1]) = 0 and pi_n = pi_{n-1} R^(n), normalized over
 * all levels.
 *
 * INDEXING. The reference passes Q0, Q1, Q2 as MATLAB cell arrays with three
 * different origins: Q0{k} is Q0^(k-1), Q1{k} is Q1^(k-1) and Q2{k} is
 * Q2^(k). The port takes three std::vector, all indexed by LEVEL:
 * q0[n] for n = 0..N-1, q1[n] for n = 0..N, q2[n] for n = 1..N, with q2[0]
 * required to be present but unused (a level-0 down-block does not exist), so
 * q2.size() == q1.size(). This removes the origin mismatch that makes the
 * MATLAB call sites hard to read; the tests cross-check against MATLAB with
 * the shift applied explicitly.
 *
 * ARITHMETIC. The recursion is FINITE -- N matrix inverses and N products --
 * with no tolerance and no iteration, so ldqbd_R and ldqbd_pi instantiate at
 * Rational and return the rate matrices and the stationary law as exact
 * fractions. That is unusual for a QBD: the level-independent case needs
 * cyclic reduction and can only ever be approximate (see qbd_r.h), whereas a
 * level-dependent chain with a finite top level is a finite linear algebra
 * problem. The one exception is the singular fallback, see below.
 *
 * SINGULAR LEVELS. The reference tests abs(det(U)) > 1e-14 and falls back to
 * Q0 * pinv(U) when it fails. An absolute determinant threshold is not
 * meaningful at exact arithmetic (a rational matrix is singular or it is not,
 * and det scales like the n-th power of the entries), so the port splits:
 *   - inexact T keeps the reference behaviour, |det(U)| <= 1e-14 routes to the
 *     LAPACK pseudo-inverse of util/svd.h, in double, with the result lifted
 *     back into T;
 *   - exact T tests for an exactly singular U and raises NumericError naming
 *     the level, because a pseudo-inverse of a rational matrix is not rational
 *     and silently returning a rounded one would make the "exact" label false.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"
#include "line/util/svd.h"

namespace line {
namespace mam {

namespace ldqbd_detail {

/** Determinant through the same LU the port uses everywhere else. */
template <class T>
T det_lu(const Matrix<T>& A) {
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("ldqbd: determinant of a non-square matrix");
    if (n == 0) return num_traits<T>::from_int(1);
    Matrix<T> LU = A;
    std::vector<std::size_t> piv;
    try {
        piv = lu_factor(LU);
    } catch (const NumericError&) {
        return num_traits<T>::from_int(0);  // an exactly zero pivot: singular
    }
    T d = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < n; ++i) d *= LU(i, i);
    for (std::size_t k = 0; k < n; ++k)
        if (piv[k] != k) d = -d;
    return d;
}

/** X U^-1, with the reference's pseudo-inverse fallback when U is singular. */
template <class T>
Matrix<T> right_divide(const Matrix<T>& X, const Matrix<T>& U, std::size_t level) {
    const T d = det_lu(U);
    bool singular;
    if (num_traits<T>::is_exact) {
        singular = (d == num_traits<T>::from_int(0));
    } else {
        singular = !(num_abs(T(d)) > num_traits<T>::from_double(1e-14));
    }
    if (!singular) return matmul(X, inverse(U));

    if (num_traits<T>::is_exact)
        throw NumericError("ldqbd_R: the local block at level " + std::to_string(level) +
                           " is exactly singular; its pseudo-inverse is not rational, so the "
                           "exact instantiation refuses rather than returning a rounded value");
    Matrix<double> Ud(U.rows(), U.cols());
    for (std::size_t i = 0; i < U.rows(); ++i)
        for (std::size_t j = 0; j < U.cols(); ++j) Ud(i, j) = num_traits<T>::to_double(U(i, j));
    const Matrix<double> P = pinv(Ud);
    Matrix<T> Pt(P.rows(), P.cols());
    for (std::size_t i = 0; i < P.rows(); ++i)
        for (std::size_t j = 0; j < P.cols(); ++j) Pt(i, j) = num_traits<T>::from_double(P(i, j));
    return matmul(X, Pt);
}

}  // namespace ldqbd_detail

/**
 * Rate matrices R^(1), ..., R^(N) of a level-dependent QBD (ldqbd_R.m).
 *
 * @param q0 up-blocks, q0[n] from level n to n+1, n = 0..N-1
 * @param q1 local blocks, q1[n] at level n, n = 0..N
 * @param q2 down-blocks, q2[n] from level n to n-1, n = 1..N; q2[0] is
 *           required to be present so that the vectors line up by level, and
 *           is never read
 * @return R indexed by level, R[n] for n = 1..N; R[0] is empty
 */
template <class T>
std::vector<Matrix<T>> ldqbd_R(const std::vector<Matrix<T>>& q0, const std::vector<Matrix<T>>& q1,
                               const std::vector<Matrix<T>>& q2) {
    if (q1.size() < 2) throw InputError("ldqbd_R: at least two levels are required");
    const std::size_t N = q1.size() - 1;
    if (q0.size() < N) throw InputError("ldqbd_R: q0 must hold the up-blocks of levels 0..N-1");
    if (q2.size() < N + 1) throw InputError("ldqbd_R: q2 must be indexed by level, 0..N");

    std::vector<Matrix<T>> R(N + 1);
    // Top level: R^(N) = Q0[N-1] (-Q1[N])^-1.
    R[N] = ldqbd_detail::right_divide(
        q0[N - 1], qbd_detail::mscale(q1[N], T(num_traits<T>::from_int(-1))), N);

    for (std::size_t n = N - 1; n >= 1; --n) {
        const Matrix<T> RQ2 = matmul(R[n + 1], q2[n + 1]);
        const Matrix<T> U =
            qbd_detail::msub(qbd_detail::mscale(q1[n], T(num_traits<T>::from_int(-1))), RQ2);
        R[n] = ldqbd_detail::right_divide(q0[n - 1], U, n);
    }
    return R;
}

/** Stationary distribution of a level-dependent QBD, per level and per phase. */
template <class T>
struct LdqbdPi {
    std::vector<std::vector<T>> pi_level;  ///< pi_level[n], one entry per phase of level n
    std::vector<T> pi;                     ///< pi[n] = sum of pi_level[n], the level marginal
};

/**
 * Stationary distribution of a level-dependent QBD given its rate matrices
 * (ldqbd_pi.m).
 *
 * The level-0 vector solves pi_0 (Q1[0] + R^(1) Q2[1]) = 0, the higher levels
 * follow from pi_n = pi_{n-1} R^(n), and the whole family is normalized to
 * unit total mass.
 *
 * DIVERGENCE FROM THE REFERENCE, deliberate and tested. MATLAB extracts the
 * level-0 vector from an eigendecomposition: it takes the eigenvector of A'
 * whose eigenvalue has smallest modulus, then applies real(), abs() and a
 * normalization. Three problems with that, all avoided here:
 *   - abs() silently turns a genuinely signed null vector into a nonnegative
 *     one, hiding a mis-specified generator instead of reporting it;
 *   - the eigenvector is only accurate to the square root of the eigenvalue
 *     separation, whereas the null vector of an exactly known matrix is a
 *     linear solve;
 *   - it forces the whole routine through a double-precision eigensolver, so
 *     no exact instantiation would be possible.
 * The port solves x A = 0 with sum(x) = 1 directly (qbd_detail::statvec), one
 * linear system, exact at Rational. On a well-posed instance the two agree to
 * roundoff; the tests check that against MATLAB and separately check the
 * balance residual ||pi_0 A||_inf, which the eigen route cannot drive to zero.
 *
 * The scalar-level-0 special case of the reference (a level 0 of order one is
 * seeded with pi_0 = 1 instead of solving anything) is reproduced, because for
 * an order-one level the balance equation is 1 x 1 and any nonzero scalar is
 * its solution up to the final normalization.
 */
template <class T>
LdqbdPi<T> ldqbd_pi(const std::vector<Matrix<T>>& R, const std::vector<Matrix<T>>& q0,
                    const std::vector<Matrix<T>>& q1, const std::vector<Matrix<T>>& q2) {
    (void)q0;
    if (q1.size() < 2) throw InputError("ldqbd_pi: at least two levels are required");
    const std::size_t N = q1.size() - 1;
    if (R.size() != N + 1) throw InputError("ldqbd_pi: R must be indexed by level, 0..N");
    const T zero = num_traits<T>::from_int(0);

    LdqbdPi<T> out;
    out.pi_level.resize(N + 1);

    if (q1[0].rows() == 1) {
        out.pi_level[0] = std::vector<T>(1, num_traits<T>::from_int(1));
    } else {
        const Matrix<T> A = qbd_detail::madd(q1[0], matmul(R[1], q2[1]));
        out.pi_level[0] = qbd_detail::statvec(A);
    }

    for (std::size_t n = 1; n <= N; ++n) out.pi_level[n] = vecmul(out.pi_level[n - 1], R[n]);

    T total = zero;
    for (std::size_t n = 0; n <= N; ++n)
        for (const T& v : out.pi_level[n]) total += v;
    if (total == zero) throw NumericError("ldqbd_pi: the stationary vector has zero total mass");

    out.pi.assign(N + 1, zero);
    for (std::size_t n = 0; n <= N; ++n) {
        for (T& v : out.pi_level[n]) v /= total;
        T s = zero;
        for (const T& v : out.pi_level[n]) s += v;
        out.pi[n] = s;
    }
    return out;
}

/** R and the stationary distribution together (ldqbd.m). */
template <class T>
struct LdqbdResult {
    std::vector<Matrix<T>> R;
    LdqbdPi<T> pi;
};

/** Solve a level-dependent QBD: rate matrices and stationary law (ldqbd.m). */
template <class T>
LdqbdResult<T> ldqbd(const std::vector<Matrix<T>>& q0, const std::vector<Matrix<T>>& q1,
                     const std::vector<Matrix<T>>& q2) {
    LdqbdResult<T> out;
    out.R = ldqbd_R(q0, q1, q2);
    out.pi = ldqbd_pi(out.R, q0, q1, q2);
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_LDQBD_H
