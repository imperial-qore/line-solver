/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_TRAFFIC_IDC_H
#define LINE_API_NPFQN_TRAFFIC_IDC_H

/**
 * Traffic variability equations of the Robust Queueing Network Analyzer
 * (W. Whitt and W. You 2018, "A Robust Queueing Network Analyzer Based on
 * Indices of Dispersion").
 *
 * Templated port of matlab/src/api/npfqn/npfqn_traffic_idc.m. The JAR carries
 * the same algorithm in jar/src/main/java/jline/api/npfqn/Npfqn_traffic_idc.java;
 * both were cross-read and agree on every equation.
 *
 * The model is a single-class open network of K single-server FCFS queues with
 * Markovian routing P. Two pieces are computed:
 *
 *  1. npfqn_traffic_idc assembles and solves the LIMITING variability
 *     equations (eq. 42/44) for the asymptotic arrival variability parameters
 *     c2_{a,i} = I_{a,i}(Inf), together with the flow rates (eq. 20-21), the
 *     fundamental matrix Xi = (I - P')^{-1}, the splitting correction
 *     alpha_{i,j} (eq. 34) and the superposition correction beta_i (eq. 38-39).
 *     Every operation there is an addition, a multiplication, a division or a
 *     linear solve in the field of the inputs, so this stays exact at
 *     T = Rational with no reformulation.
 *
 *  2. npfqn_traffic_idc_at solves the TIME-DEPENDENT IDC equations (eq. 40/43)
 *     at a single time t and returns I_{a,i}(t). It weights the equations by
 *     the canonical RBM correlation weight w*, which is an erfc and an exp, so
 *     that half requires transcendental arithmetic.
 *
 * Deviations from MATLAB, both mechanical:
 *   - the external-arrival and service IDC handles a0IdcFun and sIdcFun are
 *     arguments of npfqn_traffic_idc_at rather than fields of the context, so
 *     the context stays a plain value type with no std::function in it;
 *   - MATLAB signals "no weight available at this station" by setting the
 *     weight argument to Inf, which w* maps to 1. An exact field has no
 *     infinity, so the port tests the same condition (h_i > 0 and c2x_i > 0)
 *     and uses the weight 1 directly. The two are identical in double.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/npfqn/npfqn_rqna_weight.h"
#include "line/api/npfqn/npfqn_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

/** Toggle for the two correction terms, mirroring MATLAB's `corrections`. */
struct TrafficIdcCorrections {
    bool alpha = true;  ///< splitting correction alpha_{i,j}, eq. (34)
    bool beta = true;   ///< superposition correction beta_i, eqs. (38)-(39)
};

/** Context returned by npfqn_traffic_idc, mirroring the MATLAB `ctx` struct. */
template <class T>
struct TrafficIdcContext {
    std::size_t K = 0;         ///< number of queues
    std::vector<T> lambda;     ///< (K) total arrival rate at each queue
    std::vector<T> lambda0;    ///< (K) external arrival rate at each queue
    std::vector<T> mu;         ///< (K) service rate
    std::vector<T> rho;        ///< (K) utilization lambda_i / mu_i
    std::vector<T> cs2;        ///< (K) service SCV
    std::vector<T> c2a0;       ///< (K) external arrival asymptotic IDC
    std::vector<T> c2a;        ///< (K) total arrival asymptotic IDC
    std::vector<T> c2d;        ///< (K) departure asymptotic IDC
    std::vector<T> c2x;        ///< (K) c2a + cs2
    Matrix<T> P;               ///< (K x K) routing matrix
    Matrix<T> Xi;              ///< (K x K) fundamental matrix (I - P')^{-1}
    Matrix<T> lam_ji;          ///< (K x K) lam_ji(j,i) = lambda_j p_{j,i}
    Matrix<T> c2aij;           ///< (K x K) flow (i,j) asymptotic IDC
    Matrix<T> c2alpha;         ///< (K x K) splitting correction
    std::vector<Matrix<T>> zetaAll;  ///< (K) each (K x K), zetaAll[i](j,k) = zeta_{j,i;k,i}
};

namespace detail {

/** Inverse of a square matrix by LU, one solve per column of the identity. */
template <class T>
Matrix<T> matrix_inverse(const Matrix<T>& A) {
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("npfqn_traffic_idc: matrix is not square");
    Matrix<T> LU = A;
    const std::vector<std::size_t> piv = lu_factor(LU);
    Matrix<T> inv(n, n, num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < n; ++j) {
        std::vector<T> e(n, num_traits<T>::from_int(0));
        e[j] = num_traits<T>::from_int(1);
        lu_solve(LU, piv, e);
        for (std::size_t i = 0; i < n; ++i) inv(i, j) = e[i];
    }
    return inv;
}

}  // namespace detail

/**
 * @param lambda0     (K) external arrival rate into each queue
 * @param P           (K x K) routing matrix among the queues, P(i,j) = p_{i,j}
 * @param c2a0        (K) asymptotic IDC (SCV) of each external arrival process
 * @param mu          (K) service rate at each queue
 * @param cs2         (K) service SCV c2_{s,i}
 * @param corrections which of the two correction terms to include
 */
template <class T>
TrafficIdcContext<T> npfqn_traffic_idc(const std::vector<T>& lambda0, const Matrix<T>& P,
                                       const std::vector<T>& c2a0, const std::vector<T>& mu,
                                       const std::vector<T>& cs2,
                                       const TrafficIdcCorrections& corrections) {
    const std::size_t K = mu.size();
    if (P.rows() != K || P.cols() != K)
        throw InputError("npfqn_traffic_idc: routing matrix and service rates disagree on K");
    if (lambda0.size() != K || c2a0.size() != K || cs2.size() != K)
        throw InputError("npfqn_traffic_idc: input vectors disagree on K");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    TrafficIdcContext<T> ctx;
    ctx.K = K;
    ctx.lambda0 = lambda0;
    ctx.mu = mu;
    ctx.cs2 = cs2;
    ctx.c2a0 = c2a0;
    ctx.P = P;

    // ----- traffic rate equations (eq. 20-21) -----
    Matrix<T> ImPt(K, K, zero);
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) ImPt(i, j) = (i == j ? one : zero) - P(j, i);
    ctx.Xi = detail::matrix_inverse(ImPt);

    ctx.lambda.assign(K, zero);
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) ctx.lambda[i] += ctx.Xi(i, j) * lambda0[j];

    ctx.rho.assign(K, zero);
    for (std::size_t i = 0; i < K; ++i) ctx.rho[i] = ctx.lambda[i] / mu[i];

    ctx.lam_ji = Matrix<T>(K, K, zero);
    for (std::size_t j = 0; j < K; ++j)
        for (std::size_t i = 0; i < K; ++i) ctx.lam_ji(j, i) = ctx.lambda[j] * P(j, i);

    // ----- correction terms (asymptotic, w*(Inf) = 1) -----
    // alpha: c2alpha_{i,j} = 2 Xi_{i,j} p_{i,j} (1 - p_{i,j})
    ctx.c2alpha = Matrix<T>(K, K, zero);
    if (corrections.alpha) {
        for (std::size_t i = 0; i < K; ++i)
            for (std::size_t j = 0; j < K; ++j)
                ctx.c2alpha(i, j) = two * ctx.Xi(i, j) * P(i, j) * (one - P(i, j));
    }

    // beta: splitting covariance Sigma_l at each station, then the Brownian
    // variance-rate matrix Amat, then zeta_{j,i;k,i} of eq. (39)
    std::vector<Matrix<T>> Sigma(K, Matrix<T>(K, K, zero));
    for (std::size_t l = 0; l < K; ++l) {
        Matrix<T> Sl(K, K, zero);
        for (std::size_t a = 0; a < K; ++a)
            for (std::size_t b = 0; b < K; ++b) Sl(a, b) = -P(l, a) * P(l, b) * ctx.lambda[l];
        for (std::size_t a = 0; a < K; ++a) Sl(a, a) = P(l, a) * (one - P(l, a)) * ctx.lambda[l];
        Sigma[l] = Sl;
    }
    Matrix<T> Amat(K, K, zero);
    for (std::size_t i = 0; i < K; ++i) Amat(i, i) = c2a0[i] * lambda0[i];
    for (std::size_t l = 0; l < K; ++l)
        for (std::size_t a = 0; a < K; ++a)
            for (std::size_t b = 0; b < K; ++b) Amat(a, b) += Sigma[l](a, b);

    ctx.zetaAll.assign(K, Matrix<T>(K, K, zero));
    std::vector<T> c2beta(K, zero);
    for (std::size_t i = 0; i < K; ++i) {
        // nu(l,:) = p_{l,i} Xi(l,:)
        Matrix<T> nu(K, K, zero);
        for (std::size_t l = 0; l < K; ++l)
            for (std::size_t b = 0; b < K; ++b) nu(l, b) = P(l, i) * ctx.Xi(l, b);
        // Z = nu Amat nu'
        Matrix<T> nuA(K, K, zero);
        for (std::size_t a = 0; a < K; ++a)
            for (std::size_t b = 0; b < K; ++b) {
                T s = zero;
                for (std::size_t c = 0; c < K; ++c) s += nu(a, c) * Amat(c, b);
                nuA(a, b) = s;
            }
        Matrix<T> Z(K, K, zero);
        for (std::size_t a = 0; a < K; ++a)
            for (std::size_t b = 0; b < K; ++b) {
                T s = zero;
                for (std::size_t c = 0; c < K; ++c) s += nuA(a, c) * nu(b, c);
                Z(a, b) = s;
            }
        // cross terms nu_k Sigma_j e_i + nu_j Sigma_k e_i
        for (std::size_t j = 0; j < K; ++j)
            for (std::size_t k = 0; k < K; ++k) {
                T s = zero;
                for (std::size_t c = 0; c < K; ++c)
                    s += nu(k, c) * Sigma[j](c, i) + nu(j, c) * Sigma[k](c, i);
                Z(j, k) += s;
            }
        ctx.zetaAll[i] = Z;
        T s = zero;
        for (std::size_t j = 0; j < K; ++j)
            for (std::size_t k = j + 1; k < K; ++k) s += Z(j, k);
        if (ctx.lambda[i] > zero) c2beta[i] = (two / ctx.lambda[i]) * s;
    }
    if (!corrections.beta) {
        c2beta.assign(K, zero);
        ctx.zetaAll.assign(K, Matrix<T>(K, K, zero));
    }

    // ----- limiting variability equations (eq. 44): (E - Minf) c = binf -----
    // variable ordering: [c2a(1..K), c2aij(i,j) row-major, c2d(1..K)]
    const std::size_t Na = K, Naij = K * K, N = Na + Naij + Na;
    Matrix<T> E(N, N, zero);
    std::vector<T> binf(N, zero);
    for (std::size_t i = 0; i < N; ++i) E(i, i) = one;
    for (std::size_t i = 0; i < K; ++i) {
        if (ctx.lambda[i] > zero) {
            for (std::size_t j = 0; j < K; ++j)
                E(i, Na + j * K + i) -= ctx.lam_ji(j, i) / ctx.lambda[i];
            binf[i] = (lambda0[i] / ctx.lambda[i]) * c2a0[i] + c2beta[i];
        }
        for (std::size_t j = 0; j < K; ++j) {
            E(Na + i * K + j, Na + Naij + i) -= P(i, j);
            binf[Na + i * K + j] = (one - P(i, j)) + ctx.c2alpha(i, j);
        }
        E(Na + Naij + i, i) -= one;
    }
    const std::vector<T> csol = solve(E, binf);

    ctx.c2a.assign(csol.begin(), csol.begin() + static_cast<long>(K));
    ctx.c2aij = Matrix<T>(K, K, zero);
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) ctx.c2aij(i, j) = csol[Na + i * K + j];
    ctx.c2d.assign(csol.begin() + static_cast<long>(Na + Naij), csol.end());
    ctx.c2x.assign(K, zero);
    for (std::size_t i = 0; i < K; ++i) ctx.c2x[i] = ctx.c2a[i] + cs2[i];
    return ctx;
}

/** Default corrections (both on), matching MATLAB's `nargin < 8` branch. */
template <class T>
TrafficIdcContext<T> npfqn_traffic_idc(const std::vector<T>& lambda0, const Matrix<T>& P,
                                       const std::vector<T>& c2a0, const std::vector<T>& mu,
                                       const std::vector<T>& cs2) {
    return npfqn_traffic_idc(lambda0, P, c2a0, mu, cs2, TrafficIdcCorrections());
}

/**
 * Time-dependent IDC equations (eq. 43) at a single time t, the port of the
 * MATLAB nested function local_idc_at reached through ctx.IaFun.
 *
 * @param ctx      context from npfqn_traffic_idc
 * @param t        time argument
 * @param a0IdcFun external-arrival IDC, a0IdcFun(t) -> (K) vector
 * @param sIdcFun  service IDC, sIdcFun(rho .* t) -> (K) vector
 * @return (K) total arrival IDCs I_{a,i}(t)
 */
template <class T>
std::vector<T> npfqn_traffic_idc_at(const TrafficIdcContext<T>& ctx, const T& t,
                                    const std::function<std::vector<T>(const T&)>& a0IdcFun,
                                    const std::function<std::vector<T>(const std::vector<T>&)>& sIdcFun) {
    static_assert(num_traits<T>::has_transcendental,
                  "npfqn_traffic_idc_at requires transcendental arithmetic");
    const std::size_t K = ctx.K;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    // tuning function h(rho) = rho^2
    std::vector<T> h(K, zero);
    for (std::size_t i = 0; i < K; ++i) h[i] = ctx.rho[i] * ctx.rho[i];

    // departure weights w_i(t) = w*((1-rho_i)^2 lambda_i t / (h_i c2x_i))
    std::vector<T> w(K, one);
    for (std::size_t i = 0; i < K; ++i) {
        if (h[i] > zero && ctx.c2x[i] > zero) {
            const T num = (one - ctx.rho[i]) * (one - ctx.rho[i]) * ctx.lambda[i] * t;
            w[i] = npfqn_rqna_weight(T(num / (h[i] * ctx.c2x[i])));
        }  // else MATLAB passes Inf, whose weight is 1
    }

    const std::vector<T> Ia0 = a0IdcFun(t);
    std::vector<T> rhot(K, zero);
    for (std::size_t i = 0; i < K; ++i) rhot[i] = ctx.rho[i] * t;
    const std::vector<T> Is = sIdcFun(rhot);  // service IDC at scaled time rho*t
    if (Ia0.size() != K || Is.size() != K)
        throw InputError("npfqn_traffic_idc_at: an IDC handle returned the wrong length");

    // alpha_{i,j}(t) = c2alpha_{i,j} w_i(t)
    Matrix<T> alpha_t(K, K, zero);
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) alpha_t(i, j) = ctx.c2alpha(i, j) * w[i];

    // beta_i(t) = (1/lambda_i) sum_{j != k} zeta_{j,i;k,i} w*(arg_j)
    std::vector<T> beta_t(K, zero);
    for (std::size_t i = 0; i < K; ++i) {
        const Matrix<T>& Z = ctx.zetaAll[i];
        std::vector<T> wj(K, zero);
        for (std::size_t j = 0; j < K; ++j) {
            if (h[j] > zero && ctx.c2x[j] > zero && ctx.P(j, i) > zero) {
                const T aj = (one - ctx.rho[j]) * (one - ctx.rho[j]) * ctx.P(j, i) * ctx.lambda[j] * t;
                wj[j] = npfqn_rqna_weight(T(aj / (h[j] * ctx.c2x[j])));
            }
        }
        T s = zero;
        for (std::size_t j = 0; j < K; ++j)
            for (std::size_t k = 0; k < K; ++k)
                if (j != k) s += Z(j, k) * wj[j];
        if (ctx.lambda[i] > zero) beta_t[i] = s / ctx.lambda[i];
    }

    // assemble (E - M(t)) I = b(t)
    const std::size_t Na = K, Naij = K * K, N = Na + Naij + Na;
    Matrix<T> E(N, N, zero);
    std::vector<T> b(N, zero);
    for (std::size_t i = 0; i < N; ++i) E(i, i) = one;
    for (std::size_t i = 0; i < K; ++i) {
        if (ctx.lambda[i] > zero) {
            for (std::size_t j = 0; j < K; ++j)
                E(i, Na + j * K + i) -= ctx.lam_ji(j, i) / ctx.lambda[i];
            b[i] = (ctx.lambda0[i] / ctx.lambda[i]) * Ia0[i] + beta_t[i];
        }
        for (std::size_t j = 0; j < K; ++j) {
            E(Na + i * K + j, Na + Naij + i) -= ctx.P(i, j);
            b[Na + i * K + j] = (one - ctx.P(i, j)) + alpha_t(i, j);
        }
        E(Na + Naij + i, i) -= w[i];
        b[Na + Naij + i] = (one - w[i]) * Is[i];
    }
    const std::vector<T> sol = solve(E, b);
    return std::vector<T>(sol.begin(), sol.begin() + static_cast<long>(K));
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_TRAFFIC_IDC_H
