/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_ME_ME_TYPES_H
#define LINE_API_ME_ME_TYPES_H

/**
 * Shared declarations for the maximum-entropy (Kouvatsos) queueing network
 * algorithms.
 *
 * Templated port of matlab/src/api/me/, cross-checked against
 * jar/src/main/java/jline/api/nc/Me_oqn.java, Me_cqn.java and Me_mqn.java.
 *
 * Reference: D.D. Kouvatsos, "Entropy Maximisation and Queueing Network
 * Models", Annals of Operations Research 48:63-126, 1994.
 *
 * CONVENTIONS
 *  - Stations are indexed 0..M-1, classes 0..R-1.
 *  - Per-station, per-class data (arrival rates, service rates, scvs) are
 *    (M x R) matrices; routing is a vector of R (M x M) matrices with
 *    P[r](j,i) the probability that a class-r job moves from j to i.
 *  - The server count is a vector of longs with 0 standing for an
 *    INFINITE-SERVER station. MATLAB and the JAR use Inf for this; a
 *    templated port cannot rely on T having an infinity (Rational does not),
 *    and the algorithms only ever test isinf(c(i)), never arithmetic on it.
 *  - `insens` marks stations with an insensitive discipline (PS, LCFS-PR),
 *    which are solved by the product-form mean queue length instead of the
 *    GE-type FCFS formula.
 *
 * ARITHMETIC
 *  Every function in this domain is a damped fixed-point iteration stopped by
 *  a relative tolerance, and me_cqn additionally evaluates its Lagrangian
 *  coefficient functions through log-gamma and exp. They therefore all carry
 *    static_assert(num_traits<T>::has_transcendental)
 *  and are instantiated for double and Real50 only. Raising the precision is
 *  a legitimate use here: the ME coefficients of (3.8) are products of up to
 *  sum(N) factors, and the convolution that normalizes them cancels heavily
 *  at high population, which is precisely where the double solution starts to
 *  lose its population constraint sum_i L(i,r) = N(r).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace me {

/** Iteration control, mirroring the MATLAB options struct. */
struct MeOptions {
    double tol = 1e-6;
    long maxiter = 1000;
};

namespace detail {

/** exp(v), resolved by ADL. */
template <class T>
inline T num_exp(const T& v) {
    using std::exp;
    return exp(v);
}

/** log(v), resolved by ADL. */
template <class T>
inline T num_log(const T& v) {
    using std::log;
    return log(v);
}

/** sqrt(v), resolved by ADL. */
template <class T>
inline T num_sqrt(const T& v) {
    using std::sqrt;
    return sqrt(v);
}

/** log(k!), the gammaln(k+1) of the MATLAB source. */
template <class T>
inline T log_factorial(long k) {
    T s = num_traits<T>::from_int(0);
    for (long j = 2; j <= k; ++j) s += num_log(num_traits<T>::from_int(j));
    return s;
}

/** true when station i has an infinite number of servers. */
inline bool is_is(const std::vector<long>& c, std::size_t i) { return c[i] <= 0; }

/** Solve A x = b by Gaussian elimination with partial pivoting. */
template <class T>
std::vector<T> linear_solve(Matrix<T> A, std::vector<T> b) {
    const std::size_t n = A.rows();
    if (A.cols() != n || b.size() != n)
        throw InputError("me: linear solve with inconsistent dimensions");
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        T amax = num_abs(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            const T a = num_abs(A(i, k));
            if (a > amax) {
                amax = a;
                p = i;
            }
        }
        if (amax == zero) throw NumericError("me: singular flow balance equations");
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
            std::swap(b[k], b[p]);
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            const T f = A(i, k) / A(k, k);
            if (f == zero) continue;
            for (std::size_t j = k; j < n; ++j) A(i, j) -= f * A(k, j);
            b[i] -= f * b[k];
        }
    }
    std::vector<T> x(n, zero);
    for (std::size_t ii = n; ii-- > 0;) {
        T s = b[ii];
        for (std::size_t j = ii + 1; j < n; ++j) s -= A(ii, j) * x[j];
        x[ii] = s / A(ii, ii);
    }
    return x;
}

/** Validates the shared (M x R) / routing / server arguments. */
template <class T>
inline void check_dims(std::size_t M, std::size_t R, const Matrix<T>& mu, const Matrix<T>& Cs,
                       const std::vector<Matrix<T>>& P, const std::vector<long>& c,
                       const std::vector<char>& insens, const char* who) {
    if (M == 0 || R == 0) throw InputError(std::string(who) + ": M and R must be positive");
    if (mu.rows() != M || mu.cols() != R || Cs.rows() != M || Cs.cols() != R)
        throw InputError(std::string(who) + ": mu and Cs must be M x R");
    if (P.size() != R) throw InputError(std::string(who) + ": one routing matrix per class");
    for (std::size_t r = 0; r < R; ++r)
        if (P[r].rows() != M || P[r].cols() != M)
            throw InputError(std::string(who) + ": each routing matrix must be M x M");
    if (c.size() != M) throw InputError(std::string(who) + ": one server count per station");
    if (insens.size() != M) throw InputError(std::string(who) + ": one insens flag per station");
    // insens is consulted only at single-server stations, exactly as in the
    // references, so a flag set at a multiserver station is simply inert.
}

}  // namespace detail

/** Mean-value results shared by the open, closed and mixed algorithms. */
template <class T>
struct MeResult {
    Matrix<T> L;       ///< mean queue lengths (M x R)
    Matrix<T> W;       ///< mean response times (M x R)
    Matrix<T> Ca;      ///< arrival scvs (M x R)
    Matrix<T> Cd;      ///< departure scvs (M x R)
    Matrix<T> lambda;  ///< per-station throughputs (M x R)
    Matrix<T> rho;     ///< utilizations (M x R)
    std::vector<T> X;  ///< class throughputs (R)
    long iter = 0;     ///< fixed-point iterations performed
    bool converged = false;
};

}  // namespace me
}  // namespace line

#endif  // LINE_API_ME_ME_TYPES_H
