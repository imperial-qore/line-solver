/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_INFRADIUS_H_H
#define LINE_API_PFQN_INFRADIUS_H_H

/**
 * Integrands of the Norlund-Rice inversion of the normalizing constant, in
 * their two changes of variable.
 *
 * Templated port of matlab/src/api/pfqn/infradius_h.m (logistic substitution,
 * used by pfqn_nrl) and matlab/src/api/pfqn/infradius_hnorm.m (normal-CDF
 * substitution, used by pfqn_nrp). Both map the Cauchy contour integral for G
 * onto R^R and evaluate
 *
 *   h(x) = Re[ gld( sum_r L(:,r) e^{2 pi i (t_r - tbar)}, Ntot, alpha ) ] * jac(x)
 *
 * with t the substituted variable, tbar = sum_r beta_r t_r, beta = N/Ntot, and
 * the Jacobian prod_r e^{x_r}/(1+e^{x_r})^2 for the logistic map or
 * prod_r phi(x_r) for the normal one. gld is the single-class load-dependent
 * normalizing constant evaluated at COMPLEX demands.
 *
 * COMPLEX gld. MATLAB gets complex arithmetic from the language: pfqn_gld
 * detects a non-real demand vector and takes its linear-scale recursion, since
 * the logarithms it otherwise uses are not defined there. There is no
 * num_traits for a complex type here, so the recursion is carried out
 * explicitly on the Cx<T> pairs of pfqn_asympt_common.h; it is the same
 * recursion g(m,n,tm) = g(m-1,n,1) + L_m g(m,n-1,tm+1)/mu(m,tm), and it
 * remains exact in the field, only complex.
 *
 * PRECISION CEILING of the normal substitution. normcdf and normpdf need the
 * error function, which is not available uniformly across the arithmetic
 * backends, so infradius_hnorm evaluates the standard normal CDF by
 * converting to double, calling std::erfc, and converting back. That caps
 * infradius_hnorm -- and therefore pfqn_nrp -- near 1e-15 relative whatever T
 * is. infradius_h has no such ceiling: the logistic map needs only exp.
 *
 * ARITHMETIC. Both integrands are gated on num_traits<T>::has_transcendental:
 * exp, cos and sin are unavoidable, and so is erfc for the normal one.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/**
 * Single-class load-dependent normalizing constant at complex demands, the
 * linear-scale branch of pfqn_gldsingle.
 *
 * @param L  (M) complex demands, @param N total population
 * @param mu (M x >= N) load-dependent rates, real
 */
template <class T>
Cx<T> gldsingle_complex(const std::vector<Cx<T>>& L, int N, const Matrix<T>& mu) {
    const std::size_t M = L.size();
    if (M == 0) throw InputError("gldsingle_complex: empty demand vector");
    if (N < 0) throw InputError("gldsingle_complex: negative population");
    if (mu.rows() != M) throw InputError("gldsingle_complex: mu has the wrong station count");
    if (N > 0 && mu.cols() < static_cast<std::size_t>(N))
        throw InputError("gldsingle_complex: mu has fewer rate columns than the population");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t Nn = static_cast<std::size_t>(N);
    // g[m][n][tm], m = 0..M, n = 0..N, tm = 0..N+1 as in the reference.
    std::vector<std::vector<std::vector<Cx<T>>>> g(
        M + 1, std::vector<std::vector<Cx<T>>>(Nn + 1, std::vector<Cx<T>>(Nn + 2, Cx<T>())));
    for (std::size_t m = 1; m <= M; ++m) {
        for (std::size_t tm = 0; tm <= Nn; ++tm) g[m][0][tm] = Cx<T>(one);
        for (std::size_t n = 1; n <= Nn; ++n)
            for (std::size_t tm = 0; tm + n <= Nn; ++tm) {
                const T rate = mu(m - 1, tm);
                if (rate == zero)
                    throw NumericError("gldsingle_complex: a load-dependent rate is zero");
                const Cx<T> term = cx_scale(cx_mul(L[m - 1], g[m][n - 1][tm + 1]), T(one / rate));
                g[m][n][tm] = cx_add(g[m - 1][n][0], term);
            }
    }
    return g[M][Nn][0];
}

}  // namespace detail

/**
 * Logistic-substitution integrand (matlab/src/api/pfqn/infradius_h.m).
 *
 * @param x     point in R^R
 * @param L     (M x R) demands
 * @param N     (R) population
 * @param alpha (M x Ntot) load-dependent rates
 */
template <class T>
T infradius_h(const std::vector<T>& x, const Matrix<T>& L, const std::vector<T>& N,
              const Matrix<T>& alpha) {
    static_assert(num_traits<T>::has_transcendental,
                  "infradius_h requires transcendental arithmetic (logistic map and e^{2 pi i t})");
    using std::exp;
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R || x.size() != R)
        throw InputError("infradius_h: x, L and N disagree on the class count");
    const T one = num_traits<T>::from_int(1);
    T Nt = num_traits<T>::from_int(0);
    for (const T& v : N) Nt += v;

    std::vector<T> t(R);
    T tb = num_traits<T>::from_int(0), jac = one;
    for (std::size_t r = 0; r < R; ++r) {
        const T e = exp(x[r]);
        t[r] = T(e / T(one + e));
        tb += T(N[r] / Nt) * t[r];
        jac *= T(e / T(T(one + e) * T(one + e)));
    }

    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    std::vector<detail::Cx<T>> Lc(M);
    for (std::size_t m = 0; m < M; ++m) {
        detail::Cx<T> acc;
        for (std::size_t r = 0; r < R; ++r)
            acc = detail::cx_add(acc, detail::cx_scale(detail::cx_expi(T(twopi * T(t[r] - tb))),
                                                       L(m, r)));
        Lc[m] = acc;
    }
    const detail::Cx<T> G =
        detail::gldsingle_complex(Lc, static_cast<int>(num_traits<T>::to_double(Nt)), alpha);
    return T(G.re * jac);
}

/**
 * Normal-CDF substitution integrand
 * (matlab/src/api/pfqn/infradius_hnorm.m). See the precision note above.
 */
template <class T>
T infradius_hnorm(const std::vector<T>& x, const Matrix<T>& L, const std::vector<T>& N,
                  const Matrix<T>& alpha) {
    static_assert(num_traits<T>::has_transcendental,
                  "infradius_hnorm requires transcendental arithmetic (normal CDF and density)");
    using std::exp;
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R || x.size() != R)
        throw InputError("infradius_hnorm: x, L and N disagree on the class count");
    T Nt = num_traits<T>::from_int(0);
    for (const T& v : N) Nt += v;

    const double invsqrt2pi = 0.39894228040143267793994605993438;
    std::vector<T> t(R);
    T tb = num_traits<T>::from_int(0), jac = num_traits<T>::from_int(1);
    for (std::size_t r = 0; r < R; ++r) {
        const double xd = num_traits<T>::to_double(x[r]);
        t[r] = num_traits<T>::from_double(0.5 * std::erfc(-xd / std::sqrt(2.0)));
        tb += T(N[r] / Nt) * t[r];
        jac *= num_traits<T>::from_double(invsqrt2pi * std::exp(-0.5 * xd * xd));
    }

    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    std::vector<detail::Cx<T>> Lc(M);
    for (std::size_t m = 0; m < M; ++m) {
        detail::Cx<T> acc;
        for (std::size_t r = 0; r < R; ++r)
            acc = detail::cx_add(acc, detail::cx_scale(detail::cx_expi(T(twopi * T(t[r] - tb))),
                                                       L(m, r)));
        Lc[m] = acc;
    }
    const detail::Cx<T> G =
        detail::gldsingle_complex(Lc, static_cast<int>(num_traits<T>::to_double(Nt)), alpha);
    return T(G.re * jac);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_INFRADIUS_H_H
