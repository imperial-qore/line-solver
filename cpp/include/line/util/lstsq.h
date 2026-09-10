/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_LSTSQ_H
#define LINE_UTIL_LSTSQ_H

/**
 * Least squares for a rectangular system, exact-capable.
 *
 * MATLAB reaches for `qr(A,0)` when A has full column rank and for a truncated
 * SVD when it does not. Neither has an exact counterpart: a QR factor contains
 * square roots and a singular value is algebraic, not rational. The two
 * routines here give the SAME vectors by rational means.
 *
 * - Full column rank: the normal equations A'A x = A'b. Their solution IS the
 *   least-squares solution, identically, not an approximation of it; the usual
 *   objection to them is the squared condition number, which is a
 *   floating-point concern and disappears at T = Rational.
 * - Rank deficient: the minimum-norm least-squares solution x = A^+ b, built
 *   from a full-rank factorization A = C F obtained by exact elimination,
 *
 *     x = F' (F F')^{-1} (C'C)^{-1} C' b.
 *
 *   For a matrix of rank k this equals the truncated-SVD solution at the same
 *   k, so it reproduces MATLAB's fallback rather than approximating it.
 *
 * RANK DETECTION is where exact and inexact arithmetic genuinely part company.
 * At T = Rational a pivot is either zero or it is not, so the rank is exact and
 * `tol` is ignored. In floating point a threshold is unavoidable; the default
 * is max(m,n) * 1e-14 * max|A|, in the spirit of MATLAB's
 * max(size(A)) * eps(max(sv)) but computed from the entries rather than from
 * singular values that are not available here.
 *
 * NOT PORTED, deliberately: MATLAB's callers follow a rank-deficient solve with
 * a RANDOM perturbation retry (`rng(23000,'twister')`, demands nudged by
 * 1e-10..1e-4 times their scale) and keep whichever perturbed model looks
 * better behaved. That is a floating-point remedy for a floating-point rank
 * test. It has no exact counterpart, it changes the answer, and reproducing it
 * would require MATLAB's Mersenne-Twister double stream bit for bit. The exact
 * pseudoinverse above is what the perturbation is trying to approximate.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {

/** Outcome of lstsq: the solution and whether the system was rank deficient. */
template <class T>
struct LstsqResult {
    std::vector<T> x;
    std::size_t rank;
    bool rankdef;
};

namespace detail {

/** Pivot threshold: exactly zero in an exact field, relative otherwise. */
template <class T>
T lstsq_tolerance(const Matrix<T>& A) {
    const T zero = num_traits<T>::from_int(0);
    if (num_traits<T>::is_exact) return zero;
    T mx = zero;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            const T a = num_abs(T(A(i, j)));
            if (a > mx) mx = a;
        }
    const std::size_t d = A.rows() > A.cols() ? A.rows() : A.cols();
    return T(mx * num_traits<T>::from_double(1e-14) *
             num_traits<T>::from_int(static_cast<long>(d)));
}

/** A'A and A'b, the normal-equation pair. */
template <class T>
void normal_equations(const Matrix<T>& A, const std::vector<T>& b, Matrix<T>& G,
                      std::vector<T>& c) {
    const std::size_t m = A.rows(), n = A.cols();
    const T zero = num_traits<T>::from_int(0);
    G = Matrix<T>(n, n, zero);
    c.assign(n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) {
            T s = zero;
            for (std::size_t k = 0; k < m; ++k) s += A(k, i) * A(k, j);
            G(i, j) = s;
            G(j, i) = s;
        }
        T s = zero;
        for (std::size_t k = 0; k < m; ++k) s += A(k, i) * b[k];
        c[i] = s;
    }
}

}  // namespace detail

/**
 * Reduced row echelon form of A, in place, returning the pivot columns.
 * Rows beyond the returned rank are identically zero.
 */
template <class T>
std::vector<std::size_t> rref(Matrix<T>& A, const T& tol) {
    const std::size_t m = A.rows(), n = A.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<std::size_t> piv;
    std::size_t r = 0;
    for (std::size_t c = 0; c < n && r < m; ++c) {
        std::size_t p = r;
        T best = num_abs(T(A(r, c)));
        for (std::size_t i = r + 1; i < m; ++i) {
            const T a = num_abs(T(A(i, c)));
            if (a > best) {
                best = a;
                p = i;
            }
        }
        if (!(best > tol)) continue;
        if (p != r)
            for (std::size_t j = 0; j < n; ++j) std::swap(A(r, j), A(p, j));
        const T d = A(r, c);
        for (std::size_t j = 0; j < n; ++j) A(r, j) /= d;
        A(r, c) = one;
        for (std::size_t i = 0; i < m; ++i) {
            if (i == r) continue;
            const T f = A(i, c);
            if (f == zero) continue;
            for (std::size_t j = 0; j < n; ++j) A(i, j) -= f * A(r, j);
            A(i, c) = zero;
        }
        piv.push_back(c);
        ++r;
    }
    return piv;
}

/**
 * Least-squares solution of A x = b, minimum-norm when A is rank deficient.
 *
 * @param A (m x n), any shape
 * @param b (m)
 * @param tol pivot threshold; pass 0 for the exact rank
 */
template <class T>
LstsqResult<T> lstsq(const Matrix<T>& A, const std::vector<T>& b, const T& tol) {
    const std::size_t m = A.rows(), n = A.cols();
    if (b.size() != m) throw InputError("lstsq: rhs length mismatch");
    if (n == 0) throw InputError("lstsq: empty system");
    const T zero = num_traits<T>::from_int(0);

    Matrix<T> Rf = A;
    const std::vector<std::size_t> pivcols = rref(Rf, tol);
    const std::size_t k = pivcols.size();

    LstsqResult<T> res;
    res.rank = k;
    res.rankdef = (k < n);

    if (!res.rankdef) {
        Matrix<T> G;
        std::vector<T> c;
        detail::normal_equations(A, b, G, c);
        res.x = solve(G, c);
        return res;
    }
    if (k == 0) {
        res.x.assign(n, zero);
        return res;
    }

    // Full-rank factorization A = C F: C the pivot columns of A, F the k
    // nonzero rows of the reduced row echelon form.
    Matrix<T> C(m, k, zero), F(k, n, zero);
    for (std::size_t j = 0; j < k; ++j)
        for (std::size_t i = 0; i < m; ++i) C(i, j) = A(i, pivcols[j]);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < n; ++j) F(i, j) = Rf(i, j);

    // y = (C'C)^{-1} C' b
    Matrix<T> CtC;
    std::vector<T> Ctb;
    detail::normal_equations(C, b, CtC, Ctb);
    const std::vector<T> y = solve(CtC, Ctb);

    // x = F' (F F')^{-1} y
    Matrix<T> FFt(k, k, zero);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = i; j < k; ++j) {
            T s = zero;
            for (std::size_t t = 0; t < n; ++t) s += F(i, t) * F(j, t);
            FFt(i, j) = s;
            FFt(j, i) = s;
        }
    const std::vector<T> w = solve(FFt, y);
    res.x.assign(n, zero);
    for (std::size_t t = 0; t < n; ++t) {
        T s = zero;
        for (std::size_t i = 0; i < k; ++i) s += F(i, t) * w[i];
        res.x[t] = s;
    }
    return res;
}

/** Overload picking the default pivot threshold for the arithmetic in use. */
template <class T>
LstsqResult<T> lstsq(const Matrix<T>& A, const std::vector<T>& b) {
    return lstsq(A, b, detail::lstsq_tolerance(A));
}

}  // namespace line

#endif  // LINE_UTIL_LSTSQ_H
