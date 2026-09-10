/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_LE_H
#define LINE_API_PFQN_PFQN_LE_H

/**
 * Logistic expansion (LE) asymptotic approximation of the normalizing constant
 * of a closed product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_le.m (Casale, "Accelerating
 * performance inference over closed systems by asymptotic methods",
 * SIGMETRICS 2017), including its four local functions pfqn_le_fpi,
 * pfqn_le_fpiZ, pfqn_le_hessian and pfqn_le_hessianZ, which are exported here
 * because pfqn_ls needs the same mode and Hessian.
 *
 * The integral representation of G is mapped to the simplex by a logistic
 * transformation and evaluated by Laplace's method at the mode u* of the
 * transformed integrand, giving
 *
 *   log G = multinomialln([N, M-1]) + factln(M-1) + (M-1) log sqrt(2 pi)
 *           - log sqrt(det A) + sum_i log u*_i + sum_r N_r log(u*' L(:,r))
 *
 * with A the Hessian at the mode, and the analogous Z > 0 form in which the
 * mode carries an extra scale variable v*. This is Cas17 eq. (34) as published;
 * pfqn_ble (pfqn_ble.h) adds the eps->0 bias correction derived there.
 *
 * ARITHMETIC. Laplace's method is an asymptotic approximation and the formula
 * itself is a sum of logarithms, so the routine is gated on
 * num_traits<T>::has_transcendental: it has no meaning in exact arithmetic,
 * and instantiating it there would silently produce a value that is not the
 * normalizing constant.
 *
 * FIXED POINT. The mode is found by the same 1-norm fixed-point iteration as
 * MATLAB, stopped at 1e-10; the tolerance is a double constant converted into
 * T, so a Real<D> instantiation iterates to the same point, not further. That
 * is deliberate: matching MATLAB is the contract, and the Laplace error
 * dominates the fixed-point residual by many orders of magnitude anyway.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Mode of the logistic-transformed integrand, Z = 0 case (pfqn_le_fpi). */
template <class T>
std::vector<T> pfqn_le_fpi(const Matrix<T>& L, const std::vector<T>& N) {
    static_assert(num_traits<T>::has_transcendental, "pfqn_le_fpi requires transcendental arithmetic");
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_le_fpi: L and N disagree on the class count");
    T Ntot = num_traits<T>::from_int(0);
    for (const T& v : N) Ntot += v;
    const T eta = T(Ntot + num_traits<T>::from_int(static_cast<long>(M)));
    std::vector<T> u(M, T(num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(M))));
    std::vector<T> u1(M);
    const T tol = num_traits<T>::from_double(1e-10);
    for (int it = 0; it < 100000; ++it) {
        u1 = u;
        std::vector<T> uL(R, num_traits<T>::from_int(0));
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i) uL[r] += u1[i] * L(i, r);
        for (std::size_t i = 0; i < M; ++i) {
            T ui = T(num_traits<T>::from_int(1) / eta);
            for (std::size_t r = 0; r < R; ++r) {
                if (uL[r] == num_traits<T>::from_int(0)) continue;
                ui += T(N[r] / eta) * L(i, r) * u1[i] / uL[r];
            }
            u[i] = ui;
        }
        T d = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < M; ++i) d += num_abs(T(u[i] - u1[i]));
        if (d <= tol) break;
    }
    return u;
}

/** Mode of the logistic-transformed integrand, Z > 0 case (pfqn_le_fpiZ). */
template <class T>
void pfqn_le_fpiZ(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                  std::vector<T>& u, T& v) {
    static_assert(num_traits<T>::has_transcendental, "pfqn_le_fpiZ requires transcendental arithmetic");
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R || Z.size() != R)
        throw InputError("pfqn_le_fpiZ: L, N and Z disagree on the class count");
    T Ntot = num_traits<T>::from_int(0);
    for (const T& x : N) Ntot += x;
    const T eta = T(Ntot + num_traits<T>::from_int(static_cast<long>(M)));
    u.assign(M, T(num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(M))));
    // Note: eq. (35) in the SIGMETRICS 2017 paper has a spurious +1 in the v
    // equation; the correct stationary point is v = eta - sum_r xi_r*Z_r.
    v = eta;
    const T tol = num_traits<T>::from_double(1e-10);
    std::vector<T> u1(M);
    for (int it = 0; it < 100000; ++it) {
        u1 = u;
        std::vector<T> uL(R, num_traits<T>::from_int(0));
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i) uL[r] += u1[i] * L(i, r);
        for (std::size_t i = 0; i < M; ++i) {
            T ui = T(num_traits<T>::from_int(1) / eta);
            for (std::size_t r = 0; r < R; ++r) {
                const T den = T(Z[r] + v * uL[r]);
                if (den == num_traits<T>::from_int(0)) continue;
                ui += T(N[r] / eta) * T(Z[r] + v * L(i, r)) * u1[i] / den;
            }
            u[i] = ui;
        }
        T vnew = eta;
        for (std::size_t r = 0; r < R; ++r) {
            const T den = T(Z[r] + v * uL[r]);
            if (den == num_traits<T>::from_int(0)) continue;
            vnew -= T(N[r] / den) * Z[r];
        }
        T d = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < M; ++i) d += num_abs(T(u[i] - u1[i]));
        const T dv = num_abs(T(vnew - v));
        v = vnew;
        if (T(d + dv) <= tol) break;
    }
}

/** Hessian of the Z = 0 logistic integrand at the mode ((M-1) x (M-1)). */
template <class T>
Matrix<T> pfqn_le_hessian(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& u0) {
    const std::size_t M = L.rows(), R = L.cols();
    if (M < 2) throw InputError("pfqn_le_hessian: at least two stations are required");
    T Ntot = num_traits<T>::from_int(0);
    for (const T& x : N) Ntot += x;
    const T eta = T(Ntot + num_traits<T>::from_int(static_cast<long>(M)));
    std::vector<T> uL(R, num_traits<T>::from_int(0));
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i) uL[r] += u0[i] * L(i, r);

    Matrix<T> H(M - 1, M - 1, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i + 1 < M; ++i) {
        for (std::size_t j = 0; j + 1 < M; ++j) {
            if (i != j) {
                T h = T(-eta * u0[i] * u0[j]);
                for (std::size_t r = 0; r < R; ++r)
                    h += N[r] * L(i, r) * L(j, r) * T(u0[i] * u0[j]) / T(uL[r] * uL[r]);
                H(i, j) = h;
            } else {
                T rest = num_traits<T>::from_int(0);
                for (std::size_t k = 0; k < M; ++k)
                    if (k != i) rest += u0[k];
                T h = T(eta * u0[i] * rest);
                for (std::size_t r = 0; r < R; ++r) {
                    T restL = num_traits<T>::from_int(0);
                    for (std::size_t k = 0; k < M; ++k)
                        if (k != i) restL += u0[k] * L(k, r);
                    h -= N[r] * L(i, r) * u0[i] * restL / T(uL[r] * uL[r]);
                }
                H(i, i) = h;
            }
        }
    }
    return H;
}

/** Hessian of the Z > 0 logistic integrand at the mode (M x M). */
template <class T>
Matrix<T> pfqn_le_hessianZ(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                           const std::vector<T>& u, const T& v) {
    const std::size_t K = L.rows(), R = L.cols();
    T Ntot = num_traits<T>::from_int(0);
    for (const T& x : N) Ntot += x;
    const T eta = T(Ntot + num_traits<T>::from_int(static_cast<long>(K)));
    std::vector<T> uL(R, num_traits<T>::from_int(0));
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < K; ++i) uL[r] += u[i] * L(i, r);
    std::vector<T> csi(R);
    // csi2N is csi[r]*csi[r]/N[r] rewritten as N[r]/c[r]^2. Identical where both are
    // defined, but 0 rather than 0/0 for an empty class, which oner() makes routine
    // in the mean-value pipeline of pfqn_nc.
    std::vector<T> csi2N(R);
    for (std::size_t r = 0; r < R; ++r) {
        const T c = T(Z[r] + v * uL[r]);
        csi[r] = T(N[r] / c);
        csi2N[r] = T(N[r] / T(c * c));
    }
    Matrix<T> Lhat(K, R);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t r = 0; r < R; ++r) Lhat(k, r) = T(Z[r] + v * L(k, r));

    Matrix<T> A(K, K, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) {
            if (i == j) continue;
            T a = T(-eta * u[i] * u[j]);
            for (std::size_t r = 0; r < R; ++r)
                a += csi2N[r] * Lhat(i, r) * Lhat(j, r) * T(u[i] * u[j]);
            A(i, j) = a;
        }
    for (std::size_t i = 0; i < K; ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < K; ++j)
            if (j != i) s += A(i, j);
        A(i, i) = T(-s);
    }
    // MATLAB truncates A to (K-1)x(K-1) and then writes row/column K, so the
    // assembled matrix keeps order K with its last row and column rebuilt.
    Matrix<T> B(K, K, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i + 1 < K; ++i)
        for (std::size_t j = 0; j + 1 < K; ++j) B(i, j) = A(i, j);
    T akk = num_traits<T>::from_int(1);
    for (std::size_t r = 0; r < R; ++r)
        akk -= csi2N[r] * Z[r] * uL[r];
    B(K - 1, K - 1) = T(v * akk);
    for (std::size_t i = 0; i + 1 < K; ++i) {
        T a = num_traits<T>::from_int(0);
        for (std::size_t r = 0; r < R; ++r)
            a += v * u[i] * T(csi2N[r] * Lhat(i, r) * uL[r] - csi[r] * L(i, r));
        B(i, K - 1) = a;
        B(K - 1, i) = a;
    }
    return B;
}

/** Return value of pfqn_le, mirroring [Gn, lGn]. */
template <class T>
struct LeResult {
    T G;
    T lG;
};

/**
 * Logistic expansion estimate of the normalizing constant.
 *
 * @param L (M x R) demands, @param N (R) population, @param Z (R) think times
 *          (pass an empty vector or all zeros for the Z = 0 branch)
 */
template <class T>
LeResult<T> pfqn_le(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_le requires transcendental arithmetic (Laplace approximation of an integral)");
    using std::exp;
    using std::log;
    using std::sqrt;
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    LeResult<T> res;

    T Ntot = zero, Lsum = zero, Zsum = zero;
    for (const T& x : N) Ntot += x;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Lsum += L(i, r);
    for (const T& x : Z) Zsum += x;

    // Degenerate branch: no queueing stations, the delay carries everything.
    if (M == 0 || N.empty() || Ntot == zero ||
        num_traits<T>::to_double(Lsum) < 1e-4) {
        T lG = zero;
        for (std::size_t r = 0; r < R && r < N.size(); ++r) {
            lG -= detail::num_factln<T>(N[r]);
            if (!Z.empty() && Z[r] > zero) lG += N[r] * log(Z[r]);
        }
        res.lG = lG;
        res.G = exp(lG);
        return res;
    }

    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    if (Z.empty() || num_traits<T>::to_double(Zsum) < lang::GlobalConstants::Zero) {
        const std::vector<T> umax = pfqn_le_fpi(L, N);
        const Matrix<T> A = pfqn_le_hessian(L, N, umax);
        T S = zero;
        for (std::size_t r = 0; r < R; ++r) {
            T uL = zero;
            for (std::size_t i = 0; i < M; ++i) uL += umax[i] * L(i, r);
            S += N[r] * log(uL);
        }
        // multinomialln([N, M-1]) = factln(sum N + M-1) - sum factln(N) - factln(M-1)
        T mln = detail::num_factln<T>(T(Ntot + num_traits<T>::from_int(static_cast<long>(M) - 1)));
        for (std::size_t r = 0; r < R; ++r) mln -= detail::num_factln<T>(N[r]);
        mln -= detail::num_factln<T>(num_traits<T>::from_int(static_cast<long>(M) - 1));
        T lG = T(mln + detail::num_factln<T>(num_traits<T>::from_int(static_cast<long>(M) - 1)));
        lG += num_traits<T>::from_int(static_cast<long>(M) - 1) * log(sqrt(twopi));
        lG -= log(sqrt(detail::pfqn_det(A)));
        for (std::size_t i = 0; i < M; ++i) lG += log(umax[i]);
        lG += S;
        res.lG = lG;
        res.G = exp(lG);
        return res;
    }

    std::vector<T> umax;
    T vmax = zero;
    pfqn_le_fpiZ(L, N, Z, umax, vmax);
    const Matrix<T> A = pfqn_le_hessianZ(L, N, Z, umax, vmax);
    T S = zero;
    for (std::size_t r = 0; r < R; ++r) {
        T uL = zero;
        for (std::size_t i = 0; i < M; ++i) uL += umax[i] * L(i, r);
        S += N[r] * log(T(Z[r] + vmax * uL));
    }
    T lG = zero;
    for (std::size_t r = 0; r < R; ++r) lG -= detail::num_factln<T>(N[r]);
    lG -= vmax;
    lG += num_traits<T>::from_int(static_cast<long>(M)) * log(vmax);
    lG += num_traits<T>::from_int(static_cast<long>(M)) * log(sqrt(twopi));
    lG -= log(sqrt(detail::pfqn_det(A)));
    for (std::size_t i = 0; i < M; ++i) lG += log(umax[i]);
    lG += S;
    res.lG = lG;
    res.G = exp(lG);
    return res;
}

template <class T>
LeResult<T> pfqn_le(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_le(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_LE_H
