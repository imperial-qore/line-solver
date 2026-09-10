/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_DIST_H
#define LINE_API_MAM_MAP_DIST_H

/**
 * Analytic distances between continuous-time MAPs.
 *
 * Templated port of matlab/lib/kpctoolbox/map: map_exp_mul_int.m,
 * map_geo_mul_sum.m, map_dist.m, map_dist_acf.m, map_dist_lag1.m and
 * map_feastol.m. These are the objective functions the KPC fitters minimise, so
 * they are evaluated once per iteration and must agree with the reference to
 * full precision, not merely to the tolerance of a fitted MAP.
 *
 * The discrete-time twins live in dmap.h and rest on Stein equations
 * A X B - X + C = 0; the continuous ones rest on Sylvester equations
 * A X + X B + C = 0, MATLAB's three-argument lyap. Both are solved here through
 * the Kronecker form, which is exact in the rational backend where a
 * Schur-based solver could not be. The single exception is map_geo_mul_sum,
 * which still needs a Stein solve because its recursion is on the embedded
 * chain rather than on the generator.
 *
 * SIGN TRAP: MATLAB's lyap(A, B, C) solves A X + X B + C = 0, so C enters with
 * a PLUS. Reading it as the control-theory Lyapunov equation A X + X A' = -C
 * flips the sign of every distance and leaves the minimiser unchanged, which is
 * why a fitter can look healthy on top of the error.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/dmap.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace detail {

/** MATLAB lyap(A, B, C): solves A X + X B + C = 0 through the Kronecker form. */
template <class T>
Matrix<T> sylv_solve(const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C) {
    const std::size_t m = A.rows(), n = B.cols();
    if (A.cols() != m || B.rows() != n) throw InputError("sylv_solve: A and B must be square");
    if (C.rows() != m || C.cols() != n)
        throw InputError("sylv_solve: the right-hand side is not conformable");
    const std::size_t N = m * n;
    Matrix<T> M(N, N, num_traits<T>::from_int(0));
    // vec(A X) = kron(I_n, A) vec(X) and vec(X B) = kron(B', I_m) vec(X)
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t ia = 0; ia < m; ++ia)
            for (std::size_t ja = 0; ja < m; ++ja) M(j * m + ia, j * m + ja) += A(ia, ja);
    for (std::size_t jb = 0; jb < n; ++jb)
        for (std::size_t ib = 0; ib < n; ++ib)
            for (std::size_t k = 0; k < m; ++k) M(jb * m + k, ib * m + k) += B(ib, jb);
    std::vector<T> rhs = vec_colmajor(C);
    for (std::size_t k = 0; k < N; ++k) rhs[k] = -rhs[k];
    const std::vector<T> x = solve(M, rhs);
    Matrix<T> X(m, n, num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < m; ++i) X(i, j) = x[j * m + i];
    return X;
}

/** Transpose. */
template <class T>
Matrix<T> tr(const Matrix<T>& A) {
    Matrix<T> B(A.cols(), A.rows());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) B(j, i) = A(i, j);
    return B;
}

/** Row sums of -D0, the vector of total exit rates. */
template <class T>
Matrix<T> negrowsum(const Matrix<T>& D0) {
    Matrix<T> v(D0.rows(), 1, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < D0.rows(); ++i)
        for (std::size_t j = 0; j < D0.cols(); ++j) v(i, 0) -= D0(i, j);
    return v;
}

/** Outer product of a row-vector pair, u' * v. */
template <class T>
Matrix<T> outer(const std::vector<T>& u, const std::vector<T>& v) {
    Matrix<T> M(u.size(), v.size(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < u.size(); ++i)
        for (std::size_t j = 0; j < v.size(); ++j) M(i, j) = u[i] * v[j];
    return M;
}

/**
 * vec(A)' kron(X, Z) vec(B) without forming the Kronecker product.
 *
 * With the column-major vec the summand is A(i,j) X(j,l) Z(i,k) B(k,l), which
 * collapses to the Frobenius inner product of Z with A X B'.
 */
template <class T>
T quad_kron(const Matrix<T>& A, const Matrix<T>& X, const Matrix<T>& Z, const Matrix<T>& B) {
    const Matrix<T> AXBt = matmul(matmul(A, X), tr(B));
    if (AXBt.rows() != Z.rows() || AXBt.cols() != Z.cols())
        throw InputError("quad_kron: the two Sylvester solutions are not conformable");
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < Z.rows(); ++i)
        for (std::size_t j = 0; j < Z.cols(); ++j) s += Z(i, j) * AXBt(i, j);
    return s;
}

/** Reciprocal 1-norm condition number, MATLAB's rcond. */
template <class T>
T rcond1(const Matrix<T>& M) {
    const T zero = num_traits<T>::from_int(0);
    T nM = zero;
    for (std::size_t j = 0; j < M.cols(); ++j) {
        T c = zero;
        for (std::size_t i = 0; i < M.rows(); ++i) c += num_abs(M(i, j));
        if (c > nM) nM = c;
    }
    if (nM == zero) return zero;
    const Matrix<T> Mi = inverse(M);
    T nI = zero;
    for (std::size_t j = 0; j < Mi.cols(); ++j) {
        T c = zero;
        for (std::size_t i = 0; i < Mi.rows(); ++i) c += num_abs(Mi(i, j));
        if (c > nI) nI = c;
    }
    if (nI == zero) return zero;
    return num_traits<T>::from_int(1) / (nM * nI);
}

}  // namespace detail

/**
 * Integral of the product of the two interarrival densities up to lag L.
 *
 * The reference recursion is Z_1 = lyap(B0', A0, alB' alA) followed by
 * Z_i = lyap(B0', A0, B1' Z_{i-1} A1), read out as (-B0 1)' Z (-A0 1).
 */
template <class T>
T map_exp_mul_int(const Map<T>& a, const Map<T>& b, unsigned L, const std::vector<T>& alA,
                  const std::vector<T>& alB) {
    if (L == 0) throw InputError("map_exp_mul_int: the lag L must be positive");
    const Matrix<T> B0t = detail::tr(b.D0), B1t = detail::tr(b.D1);
    Matrix<T> Z = detail::sylv_solve(B0t, a.D0, detail::outer(alB, alA));
    for (unsigned i = 1; i < L; ++i)
        Z = detail::sylv_solve(B0t, a.D0, matmul(matmul(B1t, Z), a.D1));
    const Matrix<T> ea = detail::negrowsum(a.D0), eb = detail::negrowsum(b.D0);
    T d = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < eb.rows(); ++i)
        for (std::size_t j = 0; j < ea.rows(); ++j) d += eb(i, 0) * Z(i, j) * ea(j, 0);
    return d;
}

/** map_exp_mul_int with both stationary embedded distributions taken from the MAPs. */
template <class T>
T map_exp_mul_int(const Map<T>& a, const Map<T>& b, unsigned L) {
    return map_exp_mul_int(a, b, L, map_pie(a), map_pie(b));
}

/**
 * Geometrically weighted sum of the cross moments of the two embedded chains.
 *
 * When the Stein operator is numerically singular the reference returns
 * 1/rcond(M) as a large penalty rather than a distance, which is what steers the
 * fitter away from that corner of the parameter space; the same guard is kept
 * here, driven by the pivot growth of the Kronecker matrix.
 */
template <class T>
T map_geo_mul_sum(const Map<T>& a, const Map<T>& b, const std::vector<T>& alA,
                  const std::vector<T>& alB) {
    const std::size_t NA = a.order(), NB = b.order();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> negA0(NA, NA, zero), negB0(NB, NB, zero);
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NA; ++j) negA0(i, j) = -a.D0(i, j);
    for (std::size_t i = 0; i < NB; ++i)
        for (std::size_t j = 0; j < NB; ++j) negB0(i, j) = -b.D0(i, j);
    const Matrix<T> A0i = inverse(negA0), B0i = inverse(negB0);
    Matrix<T> PAh = matmul(A0i, a.D1), PBh = matmul(B0i, b.D1);
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NA; ++j) PAh(i, j) -= alA[j];
    for (std::size_t i = 0; i < NB; ++i)
        for (std::size_t j = 0; j < NB; ++j) PBh(i, j) -= alB[j];
    // The Stein operator M = I - kron(PBh', PAh) is the one the reference conditions on
    Matrix<T> M(NA * NB, NA * NB, zero);
    for (std::size_t jb = 0; jb < NB; ++jb)
        for (std::size_t ib = 0; ib < NB; ++ib)
            for (std::size_t ia = 0; ia < NA; ++ia)
                for (std::size_t ja = 0; ja < NA; ++ja)
                    M(jb * NA + ia, ib * NA + ja) = -PBh(ib, jb) * PAh(ia, ja);
    for (std::size_t k = 0; k < NA * NB; ++k) M(k, k) += one;
    const T rc = detail::rcond1(M);
    if (rc == zero) throw InputError("map_geo_mul_sum: the Stein operator is singular");
    if (rc < num_traits<T>::from_double(1e-10)) return one / rc;
    // C = sum(A0i, 2) * (alB * B0i), an NA x NB rank-one right-hand side
    std::vector<T> ra(NA, zero);
    for (std::size_t i = 0; i < NA; ++i)
        for (std::size_t j = 0; j < NA; ++j) ra[i] += A0i(i, j);
    const std::vector<T> cb = vecmul(alB, B0i);
    const Matrix<T> X = detail::stein_solve(PAh, PBh, detail::outer(ra, cb));
    // d = sum(alA * A0i * X * B0i)
    const std::vector<T> u = vecmul(vecmul(alA, A0i), X);
    const std::vector<T> w = vecmul(u, B0i);
    T d = zero;
    for (std::size_t j = 0; j < NB; ++j) d += w[j];
    return d;
}

/** map_geo_mul_sum with both embedded distributions taken from the MAPs. */
template <class T>
T map_geo_mul_sum(const Map<T>& a, const Map<T>& b) {
    return map_geo_mul_sum(a, b, map_pie(a), map_pie(b));
}

/** Squared L2 distance between the two joint densities up to lag L (map_dist.m). */
template <class T>
T map_dist(const Map<T>& a, const Map<T>& b, unsigned L, const std::vector<T>& alA,
           const std::vector<T>& alB) {
    const T two = num_traits<T>::from_int(2);
    return map_exp_mul_int(a, a, L + 1, alA, alA) - two * map_exp_mul_int(a, b, L + 1, alA, alB) +
           map_exp_mul_int(b, b, L + 1, alB, alB);
}

/** map_dist with both embedded distributions taken from the MAPs. */
template <class T>
T map_dist(const Map<T>& a, const Map<T>& b, unsigned L) {
    return map_dist(a, b, L, map_pie(a), map_pie(b));
}

/** Squared L2 distance between the two autocorrelation functions (map_dist_acf.m). */
template <class T>
T map_dist_acf(const Map<T>& a, const Map<T>& b, const std::vector<T>& alA,
               const std::vector<T>& alB) {
    const T two = num_traits<T>::from_int(2), four = num_traits<T>::from_int(4);
    const T mA = map_moment(a, 1u), m2A = map_moment(a, 2u);
    const T mB = map_moment(b, 1u), m2B = map_moment(b, 2u);
    const T varA = m2A - mA * mA, varB = m2B - mB * mB;
    if (varA == num_traits<T>::from_int(0) || varB == num_traits<T>::from_int(0))
        throw InputError("map_dist_acf: a deterministic MAP has no autocorrelation");
    return (map_geo_mul_sum(a, a, alA, alA) - m2A * m2A / four) / (varA * varA) -
           two * (map_geo_mul_sum(a, b, alA, alB) - m2A * m2B / four) / (varA * varB) +
           (map_geo_mul_sum(b, b, alB, alB) - m2B * m2B / four) / (varB * varB);
}

/** map_dist_acf with both embedded distributions taken from the MAPs. */
template <class T>
T map_dist_acf(const Map<T>& a, const Map<T>& b) {
    return map_dist_acf(a, b, map_pie(a), map_pie(b));
}

/**
 * Squared L2 distance between the two lag-one joint densities (map_dist_lag1.m).
 *
 * The quadratic form is vec(D1)' kron(X, Z) vec(D1) with X and Z the two
 * Sylvester solutions; it is accumulated directly here rather than through the
 * Kronecker product, which would be an order-four matrix in the MAP order.
 */
template <class T>
T map_dist_lag1(const Map<T>& a, const Map<T>& b, const std::vector<T>& alA,
                const std::vector<T>& alB) {
    const std::size_t NA = a.order(), NB = b.order();
    const Matrix<T> A0t = detail::tr(a.D0), B0t = detail::tr(b.D0);
    const Matrix<T> ea = detail::negrowsum(a.D0), eb = detail::negrowsum(b.D0);
    std::vector<T> av(NA), bv(NB);
    for (std::size_t i = 0; i < NA; ++i) av[i] = ea(i, 0);
    for (std::size_t i = 0; i < NB; ++i) bv[i] = eb(i, 0);
    const Matrix<T> Z_AB = detail::sylv_solve(A0t, b.D0, detail::outer(alA, alB));
    const Matrix<T> Z_AA = detail::sylv_solve(A0t, a.D0, detail::outer(alA, alA));
    const Matrix<T> Z_BB = detail::sylv_solve(B0t, b.D0, detail::outer(alB, alB));
    const Matrix<T> X_AB = detail::sylv_solve(a.D0, B0t, detail::outer(av, bv));
    const Matrix<T> X_AA = detail::sylv_solve(a.D0, A0t, detail::outer(av, av));
    const Matrix<T> X_BB = detail::sylv_solve(b.D0, B0t, detail::outer(bv, bv));
    const T two = num_traits<T>::from_int(2);
    return detail::quad_kron(b.D1, X_BB, Z_BB, b.D1) + detail::quad_kron(a.D1, X_AA, Z_AA, a.D1) -
           two * detail::quad_kron(a.D1, X_AB, Z_AB, b.D1);
}

/** map_dist_lag1 with both embedded distributions taken from the MAPs. */
template <class T>
T map_dist_lag1(const Map<T>& a, const Map<T>& b) {
    return map_dist_lag1(a, b, map_pie(a), map_pie(b));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_DIST_H
