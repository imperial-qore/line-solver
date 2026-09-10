/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LS_H
#define LINE_API_PFQN_LS_H

/**
 * Logistic-sampling estimate of the normalizing constant of a closed
 * product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ls.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/nc/Pfqn_ls.java. Reference: G. Casale,
 * "Accelerating performance inference over closed systems by asymptotic
 * methods", ACM SIGMETRICS 2017.
 *
 * The integral representation used by pfqn_le is estimated by importance
 * sampling rather than by a Laplace expansion: the proposal is the Gaussian
 * N(x0, A^{-1}) whose mode x0 and precision A are exactly the fixed point and
 * Hessian that pfqn_le computes (pfqn_le_fpi / pfqn_le_hessian for Z = 0,
 * pfqn_le_fpiZ / pfqn_le_hessianZ otherwise), and the estimate is the
 * self-normalized average of the log-weights lr = log h(x) - log q(x). Where
 * pfqn_le keeps only the quadratic term of the expansion, pfqn_ls corrects it
 * by Monte Carlo, so it is the same asymptotic device made consistent.
 *
 * Everything stays in the log domain: exponentiating either the integrand or
 * the normal density overflows once lG reaches a few hundred nats, and the
 * exp(-gammaln) prefactor then underflows, which is the NaN the reference's
 * comments record.
 *
 * Deviations from the reference, both deliberate:
 *
 *  - MATLAB routes an all-zero Z to the Z > 0 branch because it tests only
 *    isempty(Z). Here an all-zero Z goes to the Z = 0 branch, which is the same
 *    integral evaluated without running the fixed point pfqn_le_fpiZ on a
 *    degenerate v. This matches the convention already used by pfqn_le in this
 *    tree, so the two stay consistent with each other.
 *  - The proposal draw is z ~ N(0,I) mapped by the Cholesky factor of A^{-1},
 *    which is what MATLAB's mvnrnd does; the factor is computed here rather
 *    than delegated, so a non-positive-definite covariance raises a numeric
 *    error instead of a library-specific one.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. The estimate is a random variable; the
 * proposal is Gaussian, the mode comes from a tolerance-stopped fixed-point
 * iteration, and every weight is a log.
 *
 * RNG contract: see pfqn_mc_common.h. Comparable to MATLAB only in
 * distribution, never stream for stream; reproducible within this port only
 * when the generator is passed in the same state.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_le.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_ls, mirroring [Gn, lGn]. */
template <class T>
struct LsResult {
    T G;
    T lG;
};

namespace detail {

/**
 * Upper-triangular Cholesky factor C of a symmetric positive-definite A, so
 * that A = C^T C. Returns false, leaving C unspecified, when a pivot is
 * non-positive; that is MATLAB's chol flag.
 */
template <class T>
bool ls_chol_upper(const Matrix<T>& A, Matrix<T>& C) {
    using std::sqrt;
    const std::size_t n = A.rows();
    C = Matrix<T>(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) {
            T s = A(i, j);
            for (std::size_t k = 0; k < i; ++k) s -= C(k, i) * C(k, j);
            if (i == j) {
                if (!(s > num_traits<T>::from_int(0))) return false;
                C(i, i) = sqrt(s);
            } else {
                C(i, j) = s / C(i, i);
            }
        }
    }
    return true;
}

}  // namespace detail

/**
 * @param L0  (M x R) demands; rows whose total demand is below 1e-4 are
 *            dropped, as in the reference
 * @param N   (R) populations
 * @param Z   (R) think times; empty or all zero selects the Z = 0 branch
 * @param I   number of importance samples
 * @param rng explicit generator, advanced by the call
 */
template <class T>
LsResult<T> pfqn_ls(const Matrix<T>& L0, const std::vector<T>& N, const std::vector<T>& Z,
                    std::size_t I, McRng& rng) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_ls requires transcendental arithmetic: it importance-samples a Gaussian "
                  "proposal centred on a tolerance-stopped fixed point and averages log-weights");
    using std::exp;
    using std::log;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t R = N.size();
    if (!L0.empty() && L0.cols() != R)
        throw InputError("pfqn_ls: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_ls: Z has the wrong length");
    if (I == 0) throw InputError("pfqn_ls: at least one sample is required");

    // ---- drop the stations that carry no demand ----------------------------
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < L0.rows(); ++i) {
        T s = zero;
        for (std::size_t r = 0; r < R; ++r) s += L0(i, r);
        if (num_traits<T>::to_double(s) > 1e-4) keep.push_back(i);
    }
    Matrix<T> L(keep.size(), R, zero);
    for (std::size_t i = 0; i < keep.size(); ++i)
        for (std::size_t r = 0; r < R; ++r) L(i, r) = L0(keep[i], r);
    const std::size_t M = L.rows();

    T Ntot = zero, Lsum = zero, Zsum = zero;
    for (const T& x : N) Ntot += x;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Lsum += L(i, r);
    for (const T& x : Z) Zsum += x;

    LsResult<T> res;
    // ---- degenerate model: the delay carries everything --------------------
    if (M == 0 || R == 0 || Ntot == zero || num_traits<T>::to_double(Lsum) < 1e-4) {
        T lG = zero;
        for (std::size_t r = 0; r < R; ++r) {
            lG -= detail::num_factln<T>(N[r]);
            if (!Z.empty() && Z[r] > zero) lG += N[r] * log(Z[r]);
        }
        res.lG = lG;
        res.G = exp(lG);
        return res;
    }
    if (M < 2)
        throw InputError(
            "pfqn_ls: at least two loaded stations are required; the logistic transform maps the "
            "simplex to R^{M-1}, which is empty for a single station");

    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    const bool zeroZ = Z.empty() || Zsum == zero;

    // ---- proposal: mode x0 and precision A ---------------------------------
    std::size_t d = 0;
    std::vector<T> x0;
    Matrix<T> A;
    std::vector<T> umax;
    T vmax = zero;
    if (zeroZ) {
        umax = pfqn_le_fpi(L, N);
        A = pfqn_le_hessian(L, N, umax);
        d = M - 1;
        x0.resize(d);
        for (std::size_t i = 0; i < d; ++i) x0[i] = log(T(umax[i] / umax[M - 1]));
    } else {
        pfqn_le_fpiZ(L, N, Z, umax, vmax);
        A = pfqn_le_hessianZ(L, N, Z, umax, vmax);
        d = M;
        x0.resize(d);
        for (std::size_t i = 0; i + 1 < M; ++i) x0[i] = log(T(umax[i] / umax[M - 1]));
        x0[M - 1] = log(vmax);
    }
    // Symmetrize away the small asymmetries the closed-form Hessian carries.
    const T half = num_traits<T>::from_rational(1, 2);
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t j = i + 1; j < d; ++j) {
            const T m = T(A(i, j) + A(j, i)) * half;
            A(i, j) = m;
            A(j, i) = m;
        }
    const Matrix<T> iA = inverse(A);

    // log|A| from its Cholesky factor, falling back to the determinant.
    Matrix<T> Ca;
    T logdetA = zero;
    if (detail::ls_chol_upper(A, Ca)) {
        for (std::size_t i = 0; i < d; ++i) logdetA += log(Ca(i, i));
        logdetA *= num_traits<T>::from_int(2);
    } else {
        const T det = detail::pfqn_det(A);
        logdetA = log(T(det < zero ? T(-det) : det));
    }

    Matrix<T> Ci;
    if (!detail::ls_chol_upper(iA, Ci))
        throw NumericError(
            "pfqn_ls: the proposal covariance (inverse Hessian at the mode) is not positive "
            "definite, so no Gaussian proposal exists there");

    // ---- draw, evaluate the integrand and the proposal density -------------
    std::vector<double> lr(I);
    std::vector<T> xs(d), z(d), diff(d);
    const T eN = num_traits<T>::from_double(1e-10) * Ntot;
    const T eta = T(Ntot + num_traits<T>::from_int(static_cast<long>(M)) * T(one + eN));

    for (std::size_t s = 0; s < I; ++s) {
        for (std::size_t i = 0; i < d; ++i) z[i] = num_traits<T>::from_double(mc_normal01(rng));
        // x = x0 + z C, with C upper triangular and C^T C = A^{-1}.
        for (std::size_t j = 0; j < d; ++j) {
            T acc = x0[j];
            for (std::size_t i = 0; i <= j; ++i) acc += z[i] * Ci(i, j);
            xs[j] = acc;
        }

        // ---- log of the integrand -----------------------------------------
        T lT = zero;
        if (zeroZ) {
            // simplex_logfun: v = [exp(x), 1].
            T vsum = one;
            std::vector<T> vv(M, one);
            for (std::size_t i = 0; i < d; ++i) {
                vv[i] = exp(xs[i]);
                vsum += vv[i];
            }
            for (std::size_t r = 0; r < R; ++r) {
                T vl = zero;
                for (std::size_t i = 0; i < M; ++i) vl += vv[i] * L(i, r);
                lT += N[r] * log(vl);
            }
            for (std::size_t i = 0; i < d; ++i) lT += xs[i];
            lT -= T(Ntot + num_traits<T>::from_int(static_cast<long>(M))) * log(vsum);
        } else {
            const T v = exp(xs[M - 1]);
            T esum = zero;
            std::vector<T> e(M - 1, zero);
            for (std::size_t i = 0; i + 1 < M; ++i) {
                e[i] = exp(xs[i]);
                esum += e[i];
            }
            lT = -v;
            lT += num_traits<T>::from_int(static_cast<long>(M)) * T(one + eN) * xs[M - 1];
            for (std::size_t r = 0; r < R; ++r) {
                T inner = T(L(M - 1, r) * v + Z[r]);
                for (std::size_t i = 0; i + 1 < M; ++i) inner += e[i] * T(L(i, r) * v + Z[r]);
                lT += N[r] * log(inner);
            }
            for (std::size_t i = 0; i + 1 < M; ++i) lT += xs[i];
            lT -= eta * log(T(one + esum));
        }

        // ---- log of the proposal density, from the precision matrix --------
        for (std::size_t i = 0; i < d; ++i) diff[i] = T(xs[i] - x0[i]);
        T q = zero;
        for (std::size_t i = 0; i < d; ++i) {
            T row = zero;
            for (std::size_t j = 0; j < d; ++j) row += diff[j] * A(j, i);
            q += row * diff[i];
        }
        const T ldpdf =
            half * logdetA - num_traits<T>::from_rational(static_cast<long>(d), 2) * log(twopi) -
            half * q;
        lr[s] = num_traits<T>::to_double(T(lT - ldpdf));
    }

    // ---- self-normalized average, factored by the largest weight ----------
    const double lmean = mc_logmeanexp(lr);
    T lG = num_traits<T>::from_double(lmean);
    if (zeroZ) {
        // multinomialln([N, M-1]) + factln(M-1) = factln(sum N + M - 1) - sum factln(N)
        lG += detail::num_factln<T>(T(Ntot + num_traits<T>::from_int(static_cast<long>(M) - 1)));
        for (std::size_t r = 0; r < R; ++r) lG -= detail::num_factln<T>(N[r]);
    } else {
        for (std::size_t r = 0; r < R; ++r) lG -= detail::num_factln<T>(N[r]);
    }
    res.lG = lG;
    res.G = exp(lG);
    return res;
}

/** Reference default of 1e5 samples. */
template <class T>
LsResult<T> pfqn_ls(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                    McRng& rng) {
    return pfqn_ls(L, N, Z, static_cast<std::size_t>(100000), rng);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LS_H
