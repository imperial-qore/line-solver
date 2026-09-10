/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_LU_H
#define LINE_UTIL_LU_H

/**
 * LU factorization with partial pivoting, templated on the number type.
 *
 * Port of mp_pfqn's gmpla/mpq_ludcmp.c + mpq_lubksb.c, generalized from mpq_t
 * to any field type and with the 1-based Numerical-Recipes indexing removed.
 * In exact arithmetic pivoting is not needed for stability, only to avoid a
 * zero pivot, but the same largest-magnitude rule is kept so that the exact
 * and double paths take identical elimination orders and can be diffed.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {

/**
 * The ordered type the pivot search compares in.
 *
 * It is T itself for every ordered field, so the real instantiations compile to
 * the same code as before. A complex element type has no order, and its
 * specialization (in line/num/complex_number.h) reduces to the real modulus --
 * which is what MATLAB's own pivot rule uses on a complex matrix.
 */
template <class T>
struct pivot_mag {
    using type = T;
    static type of(const T& x) { return num_abs(x); }
};

/**
 * In-place LU of A (n x n). Returns the row permutation; A holds L (unit
 * diagonal, implicit) below and U on and above the diagonal.
 */
template <class T>
std::vector<std::size_t> lu_factor(Matrix<T>& A) {
    using Mag = typename pivot_mag<T>::type;
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("lu_factor: matrix is not square");
    std::vector<std::size_t> piv(n);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        Mag amax = pivot_mag<T>::of(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            Mag a = pivot_mag<T>::of(A(i, k));
            if (a > amax) {
                amax = a;
                p = i;
            }
        }
        if (amax == num_traits<Mag>::from_int(0)) throw NumericError("lu_factor: singular matrix");
        piv[k] = p;
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
        }
        const T d = A(k, k);
        for (std::size_t i = k + 1; i < n; ++i) {
            if (A(i, k) == num_traits<T>::from_int(0)) continue;
            T f = A(i, k) / d;
            A(i, k) = f;  // store the multiplier in place, as in Crout
            for (std::size_t j = k + 1; j < n; ++j) A(i, j) -= f * A(k, j);
        }
    }
    return piv;
}

/**
 * Solve LUx = Pb in place on b, using the factors from lu_factor.
 *
 * The whole permutation must be applied to b BEFORE any elimination.
 * lu_factor swaps entire rows, multiplier columns included, so LU(i,k) is the
 * multiplier of the row that ends up at position i. Interleaving the swaps
 * with the forward updates -- swap b[k], then update b[i] with LU(i,k) --
 * pairs the row identity at step k with a final-order multiplier whenever a
 * later step moves that row, and silently returns a wrong solution. It does
 * so without any warning sign: the factorization still satisfies LU = PA and
 * the pivots are all healthy. The failure needs a pivot sequence that moves an
 * already-eliminated row, which diagonally dominant matrices rarely produce,
 * so it hides until a generator with a small leading diagonal is solved.
 */
template <class T>
void lu_solve(const Matrix<T>& LU, const std::vector<std::size_t>& piv, std::vector<T>& b) {
    const std::size_t n = LU.rows();
    if (b.size() != n) throw InputError("lu_solve: rhs length mismatch");
    for (std::size_t k = 0; k < n; ++k)
        if (piv[k] != k) std::swap(b[k], b[piv[k]]);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < i; ++k) b[i] -= LU(i, k) * b[k];
    for (std::size_t i = n; i-- > 0;) {
        T s = b[i];
        for (std::size_t j = i + 1; j < n; ++j) s -= LU(i, j) * b[j];
        b[i] = s / LU(i, i);
    }
}

/**
 * Determinant of a square matrix, by the same partial-pivoting elimination.
 *
 * SINGULAR IS A VALUE HERE, NOT AN ERROR, which is why this does not go through
 * `lu_factor`: that one throws on a zero pivot because every caller of it is
 * solving a system, where a zero pivot means the question has no answer. A
 * determinant of a singular matrix is zero, a perfectly good answer, and
 * Cramer's rule needs it -- `MarkovProcess.getProbState` forms a numerator
 * matrix that IS singular whenever the state has probability zero.
 *
 * The sign comes from the parity of the row swaps, so it is exact even where
 * the product of the pivots is not.
 */
template <class T>
T lu_det(const Matrix<T>& A) {
    using Mag = typename pivot_mag<T>::type;
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("lu_det: matrix is not square");
    if (n == 0) return num_traits<T>::from_int(1);
    Matrix<T> U = A;
    const T zero = num_traits<T>::from_int(0);
    T det = num_traits<T>::from_int(1);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        Mag amax = pivot_mag<T>::of(U(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            Mag a = pivot_mag<T>::of(U(i, k));
            if (a > amax) {
                amax = a;
                p = i;
            }
        }
        if (amax == num_traits<Mag>::from_int(0)) return zero;
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(U(k, j), U(p, j));
            det = T(zero - det);
        }
        const T d = U(k, k);
        det = T(det * d);
        for (std::size_t i = k + 1; i < n; ++i) {
            if (U(i, k) == zero) continue;
            const T f = U(i, k) / d;
            for (std::size_t j = k + 1; j < n; ++j) U(i, j) -= f * U(k, j);
        }
    }
    return det;
}

/** Convenience: solve Ax = b, leaving A and b untouched. */
template <class T>
std::vector<T> solve(const Matrix<T>& A, const std::vector<T>& b) {
    Matrix<T> LU = A;
    std::vector<T> x = b;
    std::vector<std::size_t> piv = lu_factor(LU);
    lu_solve(LU, piv, x);
    return x;
}

}  // namespace line

#endif  // LINE_UTIL_LU_H
