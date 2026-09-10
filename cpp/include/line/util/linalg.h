/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_LINALG_H
#define LINE_UTIL_LINALG_H

/**
 * Dense linear algebra over the templated number type: products, identity,
 * inverse, and powers. Built on the LU in lu.h, so the exact instantiation
 * inverts a rational matrix exactly, with no residual at all.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {

/** Identity of order n. */
template <class T>
Matrix<T> eye(std::size_t n) {
    Matrix<T> I(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) I(i, i) = num_traits<T>::from_int(1);
    return I;
}

/** Matrix product A B. */
template <class T>
Matrix<T> matmul(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.cols() != B.rows()) throw InputError("matmul: inner dimensions disagree");
    Matrix<T> C(A.rows(), B.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t k = 0; k < A.cols(); ++k) {
            const T a = A(i, k);
            if (a == num_traits<T>::from_int(0)) continue;
            for (std::size_t j = 0; j < B.cols(); ++j) C(i, j) += a * B(k, j);
        }
    return C;
}

/** Row vector times matrix, v A. */
template <class T>
std::vector<T> vecmul(const std::vector<T>& v, const Matrix<T>& A) {
    if (v.size() != A.rows()) throw InputError("vecmul: dimensions disagree");
    std::vector<T> r(A.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (v[i] == num_traits<T>::from_int(0)) continue;
        for (std::size_t j = 0; j < A.cols(); ++j) r[j] += v[i] * A(i, j);
    }
    return r;
}

/** Matrix times column vector, A v. */
template <class T>
std::vector<T> mulvec(const Matrix<T>& A, const std::vector<T>& v) {
    if (v.size() != A.cols()) throw InputError("mulvec: dimensions disagree");
    std::vector<T> r(A.rows(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) r[i] += A(i, j) * v[j];
    return r;
}

/** Inverse by LU with one factorization and n back substitutions. */
template <class T>
Matrix<T> inverse(const Matrix<T>& A) {
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("inverse: matrix is not square");
    Matrix<T> LU = A;
    const std::vector<std::size_t> piv = lu_factor(LU);
    Matrix<T> X(n, n, num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < n; ++c) {
        std::vector<T> e(n, num_traits<T>::from_int(0));
        e[c] = num_traits<T>::from_int(1);
        lu_solve(LU, piv, e);
        for (std::size_t i = 0; i < n; ++i) X(i, c) = e[i];
    }
    return X;
}

/** Integer matrix power, by repeated squaring. */
template <class T>
Matrix<T> matpow(const Matrix<T>& A, unsigned k) {
    if (A.rows() != A.cols()) throw InputError("matpow: matrix is not square");
    Matrix<T> R = eye<T>(A.rows());
    Matrix<T> B = A;
    unsigned e = k;
    while (e > 0) {
        if (e & 1u) R = matmul(R, B);
        B = matmul(B, B);
        e >>= 1;
    }
    return R;
}

/** Column vector of ones, the ubiquitous e in MAP algebra. */
template <class T>
std::vector<T> ones(std::size_t n) {
    return std::vector<T>(n, num_traits<T>::from_int(1));
}

}  // namespace line

#endif  // LINE_UTIL_LINALG_H
