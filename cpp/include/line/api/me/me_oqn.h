/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_ME_ME_OQN_H
#define LINE_API_ME_ME_OQN_H

/**
 * Maximum-entropy algorithm for open multiclass queueing networks.
 *
 * Templated port of matlab/src/api/me/me_oqn.m, cross-checked against
 * jar/src/main/java/jline/api/nc/Me_oqn.java. Implements Kouvatsos (1994)
 * Section 3.2 with the GE/GE/c building block of eq. (3.9) and the GE/GE/inf
 * block, i.e. a GE-type decomposition fixed point:
 *
 *   1. feedback correction: a self-loop of probability p_ii gives a geometric
 *      number of passes, so mu <- mu (1-p_ii), Cs <- p_ii + (1-p_ii) Cs and
 *      the residual routing is renormalized;
 *   2. job flow balance lambda = lambda0 + P' lambda on the ORIGINAL routing,
 *      so lambda counts revisits, while lambda_eff = lambda (1 - p_ii) is the
 *      rate seen by the corrected queue;
 *   3. mean queue lengths from the multiclass GE/GE/1/FCFS formula of
 *      Section 3.1.1, the product-form formula at insensitive stations, or
 *      eq. (3.9) on the class-aggregated stream at multiserver stations;
 *   4. departure scvs from eq. (3.6) with the marginal utilization of (3.3);
 *   5. arrival scvs by the GE merging formula (3.7) applied to the thinned
 *      departure streams; iterate 3-5 to convergence.
 *
 * DIVERGENCE from the references on unstable input: MATLAB raises a warning
 * and returns L = Inf at any finite-server station with total utilization
 * >= 1, and continues to iterate the scvs on top of that. This port raises
 * NumericError instead: the templated backends have no infinity (and no
 * warning channel), and an Inf queue length propagated into me_mqn's
 * inflation step would silently poison the open-class results.
 *
 * ARITHMETIC: a tolerance-stopped fixed point.
 *   static_assert(num_traits<T>::has_transcendental)
 */

#include <cstddef>
#include <vector>

#include "line/api/me/me_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace me {

namespace detail {

/**
 * Mean queue length of a stable GE/GE/c/FCFS queue, the exact ME solution of
 * Kouvatsos (1994) eq. (3.9).
 */
template <class T>
T ge_gec_mql(const T& lambda, const T& Ca, const T& mu, const T& Cs, long c) {
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T alpha2 = two / (Cs + one);
    const T alpha1 = one - alpha2;
    const T beta2 = two / (Ca + one);
    const T beta1 = one - beta2;
    const T lambda2 = beta2 * lambda;
    const T mu2 = alpha2 * mu;

    std::vector<T> g(static_cast<std::size_t>(c));
    for (long j = 1; j <= c - 1; ++j) {
        const T jt = num_traits<T>::from_int(j);
        const T den = jt * mu2 * (one - alpha1 * beta1);
        if (den == num_traits<T>::from_int(0))
            throw NumericError("me_oqn: degenerate GE/GE/c building block");
        g[static_cast<std::size_t>(j - 1)] =
            (lambda2 + num_traits<T>::from_int(j - 1) * mu2 * beta1) * alpha2 / den;
    }
    const T ct = num_traits<T>::from_int(c);
    const T denc = lambda2 * alpha1 + ct * mu2;
    if (denc == num_traits<T>::from_int(0))
        throw NumericError("me_oqn: degenerate GE/GE/c building block");
    g[static_cast<std::size_t>(c - 1)] =
        (lambda2 + num_traits<T>::from_int(c - 1) * mu2 * beta1) * alpha2 / denc;
    const T x = (lambda2 + ct * mu2 * beta1) / denc;
    if (x >= one) throw NumericError("me_oqn: the multiserver station is unstable");

    std::vector<T> Gn(static_cast<std::size_t>(c));
    T acc = num_traits<T>::from_int(1);
    for (long j = 0; j < c; ++j) {
        acc *= g[static_cast<std::size_t>(j)];
        Gn[static_cast<std::size_t>(j)] = acc;
    }
    T Z = one;
    for (long j = 0; j < c - 1; ++j) Z += Gn[static_cast<std::size_t>(j)];
    Z += Gn[static_cast<std::size_t>(c - 1)] / (one - x);
    T S1 = num_traits<T>::from_int(0);
    for (long n = 1; n <= c - 1; ++n)
        S1 += num_traits<T>::from_int(n) * Gn[static_cast<std::size_t>(n - 1)];
    const T S2 = Gn[static_cast<std::size_t>(c - 1)] *
                 (ct / (one - x) + x / ((one - x) * (one - x)));
    return (S1 + S2) / Z;
}

}  // namespace detail

/**
 * @param M       number of stations
 * @param R       number of classes
 * @param lambda0 external arrival rates (M x R)
 * @param Ca0     external arrival scvs (M x R)
 * @param mu      service rates (M x R)
 * @param Cs      service scvs (M x R)
 * @param P       routing, R matrices (M x M), P[r](j,i)
 * @param c       servers per station, 0 for an infinite-server station
 * @param insens  insensitive discipline flags per station
 * @param opt     tolerance and iteration budget
 */
template <class T>
MeResult<T> me_oqn(std::size_t M, std::size_t R, const Matrix<T>& lambda0, const Matrix<T>& Ca0,
                   const Matrix<T>& mu, const Matrix<T>& Cs, const std::vector<Matrix<T>>& P,
                   const std::vector<long>& c, const std::vector<char>& insens,
                   const MeOptions& opt = MeOptions()) {
    static_assert(num_traits<T>::has_transcendental, "me_oqn requires transcendental arithmetic");
    detail::check_dims(M, R, mu, Cs, P, c, insens, "me_oqn");
    if (lambda0.rows() != M || lambda0.cols() != R || Ca0.rows() != M || Ca0.cols() != R)
        throw InputError("me_oqn: lambda0 and Ca0 must be M x R");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T tol = num_traits<T>::from_double(opt.tol);

    // Step 1: feedback correction
    std::vector<Matrix<T>> Peff = P;
    Matrix<T> mueff = mu, Cseff = Cs;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < R; ++r) {
            const T pii = P[r](i, i);
            if (!(pii > zero)) continue;
            if (pii >= one) throw InputError("me_oqn: a self-loop probability of one");
            mueff(i, r) = mu(i, r) * (one - pii);
            Cseff(i, r) = pii + (one - pii) * Cs(i, r);
            for (std::size_t j = 0; j < M; ++j) Peff[r](i, j) = P[r](i, j) / (one - pii);
            Peff[r](i, i) = zero;
        }
    }

    // Step 3: job flow balance on the original routing
    MeResult<T> out;
    out.lambda = Matrix<T>(M, R, zero);
    Matrix<T> lameff(M, R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        Matrix<T> A(M, M, zero);
        std::vector<T> b(M, zero);
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < M; ++j) A(i, j) = (i == j ? one : zero) - P[r](j, i);
            b[i] = lambda0(i, r);
        }
        const std::vector<T> x = detail::linear_solve(A, b);
        for (std::size_t i = 0; i < M; ++i) {
            out.lambda(i, r) = x[i];
            lameff(i, r) = x[i] * (one - P[r](i, i));
        }
    }

    // Utilizations
    out.rho = Matrix<T>(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            if (!(mu(i, r) > zero)) continue;
            if (detail::is_is(c, i))
                out.rho(i, r) = lameff(i, r) / mueff(i, r);
            else
                out.rho(i, r) = out.lambda(i, r) / (num_traits<T>::from_int(c[i]) * mu(i, r));
        }
    for (std::size_t i = 0; i < M; ++i) {
        if (detail::is_is(c, i)) continue;
        T s = zero;
        for (std::size_t r = 0; r < R; ++r) s += out.rho(i, r);
        if (s >= one) throw NumericError("me_oqn: the network is unstable, utilization >= 1");
    }

    out.Ca = Matrix<T>(M, R, one);
    out.Cd = Matrix<T>(M, R, one);
    out.L = Matrix<T>(M, R, zero);

    T delta = zero;
    for (long it = 1; it <= opt.maxiter; ++it) {
        out.iter = it;
        const Matrix<T> Caold = out.Ca;

        // Step 4: mean queue lengths
        for (std::size_t i = 0; i < M; ++i) {
            T rho_i = zero;
            for (std::size_t r = 0; r < R; ++r) rho_i += out.rho(i, r);
            if (detail::is_is(c, i)) {
                for (std::size_t r = 0; r < R; ++r)
                    if (lameff(i, r) > zero && mueff(i, r) > zero)
                        out.L(i, r) = lameff(i, r) / mueff(i, r);
            } else if (c[i] == 1) {
                if (insens[i]) {
                    for (std::size_t r = 0; r < R; ++r)
                        if (lameff(i, r) > zero && mueff(i, r) > zero)
                            out.L(i, r) = out.rho(i, r) / (one - rho_i);
                } else {
                    T resid = zero;
                    for (std::size_t u = 0; u < R; ++u)
                        if (lameff(i, u) > zero && mueff(i, u) > zero)
                            resid += lameff(i, u) * (Cseff(i, u) + out.Ca(i, u)) /
                                     (mueff(i, u) * mueff(i, u));
                    for (std::size_t r = 0; r < R; ++r)
                        if (lameff(i, r) > zero && mueff(i, r) > zero)
                            out.L(i, r) = out.rho(i, r) * (out.Ca(i, r) + one) / two +
                                          lameff(i, r) * resid / (two * (one - rho_i));
                }
            } else {
                T lam_a = zero;
                for (std::size_t u = 0; u < R; ++u)
                    if (lameff(i, u) > zero && mueff(i, u) > zero) lam_a += lameff(i, u);
                if (!(lam_a > zero)) continue;
                T inv_a = zero, ES = zero, ES2 = zero;
                for (std::size_t u = 0; u < R; ++u) {
                    if (!(lameff(i, u) > zero && mueff(i, u) > zero)) continue;
                    const T wu = lameff(i, u) / lam_a;
                    inv_a += wu / (out.Ca(i, u) + one);
                    ES += wu / mueff(i, u);
                    ES2 += wu * (Cseff(i, u) + one) / (mueff(i, u) * mueff(i, u));
                }
                const T Ca_a = -one + one / inv_a;
                const T Cs_a = ES2 / (ES * ES) - one;
                const T L_a = detail::ge_gec_mql(lam_a, Ca_a, T(one / ES), Cs_a, c[i]);
                const T Lq_a = L_a - lam_a * ES;
                for (std::size_t r = 0; r < R; ++r)
                    if (lameff(i, r) > zero && mueff(i, r) > zero)
                        out.L(i, r) = num_traits<T>::from_int(c[i]) * out.rho(i, r) +
                                      (lameff(i, r) / lam_a) * Lq_a;
            }
        }

        // Step 5a: departure scvs
        for (std::size_t j = 0; j < M; ++j) {
            T rho_j = zero;
            for (std::size_t r = 0; r < R; ++r) rho_j += out.rho(j, r);
            for (std::size_t r = 0; r < R; ++r) {
                if (!(lameff(j, r) > zero)) continue;
                if (detail::is_is(c, j)) {
                    out.Cd(j, r) = out.Ca(j, r);
                } else if (c[j] == 1) {
                    const T den = out.L(j, r) + rho_j - out.rho(j, r);
                    const T rhohat = den == zero ? zero : out.rho(j, r) * out.L(j, r) / den;
                    out.Cd(j, r) = two * out.L(j, r) * (one - rhohat) +
                                   out.Ca(j, r) * (one - two * rhohat);
                } else {
                    out.Cd(j, r) = rho_j * (one - rho_j) + (one - rho_j) * out.Ca(j, r) +
                                   rho_j * rho_j * Cseff(j, r);
                }
            }
        }

        // Step 5b: arrival scvs by GE merging with thinning
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t r = 0; r < R; ++r) {
                if (!(lameff(i, r) > zero)) continue;
                T sum_inv = zero;
                for (std::size_t j = 0; j < M; ++j) {
                    const T pji = Peff[r](j, i);
                    if (!(pji > zero) || !(lameff(j, r) > zero)) continue;
                    const T Cdji = one + pji * (out.Cd(j, r) - one);
                    sum_inv += (lameff(j, r) * pji / lameff(i, r)) / (Cdji + one);
                }
                if (lambda0(i, r) > zero)
                    sum_inv += (lambda0(i, r) / lameff(i, r)) / (Ca0(i, r) + one);
                if (sum_inv > zero) out.Ca(i, r) = -one + one / sum_inv;
            }
        }

        delta = zero;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                const T d = num_abs(T(out.Ca(i, r) - Caold(i, r)));
                if (d > delta) delta = d;
            }
        if (delta < tol) {
            out.converged = true;
            break;
        }
    }

    // Step 6: response times by Little's law on the visit-inclusive rates
    out.W = Matrix<T>(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (out.lambda(i, r) > zero) out.W(i, r) = out.L(i, r) / out.lambda(i, r);

    out.X.assign(R, zero);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i) out.X[r] += lambda0(i, r);
    return out;
}

}  // namespace me
}  // namespace line

#endif  // LINE_API_ME_ME_OQN_H
