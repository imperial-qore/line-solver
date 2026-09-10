/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1_FB_H
#define LINE_API_QSYS_QSYS_MG1_FB_H

/**
 * M/G/1 under FB (feedback), also called LAS (least attained service).
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1_fb.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mg1_fb.java.
 *
 * The job with the least attained service holds the server, so priority
 * depends on age, not on the original or remaining size. For a job of size x
 * (Wierman and Harchol-Balter, SIGMETRICS 2003, Sec. 3.3)
 *
 *   rho_x   = lambda int_0^x Fbar(t) dt
 *   num(x)  = lambda int_0^x t Fbar(t) dt
 *   E[T(x)] = num(x)/(1-rho_x)^2 + x/(1-rho_x)
 *
 * with Fbar the tail of the mixture job-size law. MATLAB takes two paths:
 *
 *   all cs_i = 1 (within 1e-6): the mixture is a mixture of exponentials, the
 *     two truncated integrals are closed forms, and the class mean is
 *     E[T_k] = int_0^{20/mu_k} E[T(x)] mu_k e^{-mu_k x} dx by quadrature;
 *   otherwise: the class is collapsed onto its mean size x = 1/mu_k and the
 *     non-exponential truncated integrals are replaced by the bounded
 *     surrogates min(x, 1/mu_i) and min(x^2/2, 1/mu_i^2).
 *
 * Both paths are reproduced verbatim, including the truncation of the outer
 * integral at 20 mean service times, which is not an implementation detail:
 * it biases the class mean low by the tail beyond 20 e-foldings, of relative
 * order 1e-8, and a port that integrated to infinity would not reproduce the
 * reference.
 *
 * ARITHMETIC. exp and the adaptive quadrature make this transcendental.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_mg1_setf.h"
#include "line/api/qsys/qsys_quadrature.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

namespace detail {

/** E[T(x)] under FB for a mixture-of-exponentials job-size law. */
template <class T>
T mg1_fb_response(const T& x, const std::vector<T>& mu, const std::vector<T>& p,
                  const T& lambda_total) {
    const T one = num_traits<T>::from_int(1);
    T rho_x = num_traits<T>::from_int(0), numer = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < mu.size(); ++i) {
        const T e = num_exp(T(-mu[i] * x));
        rho_x += p[i] * lambda_total * ((one - e) / mu[i]);
        numer += p[i] * lambda_total * ((one - e * (one + mu[i] * x)) / (mu[i] * mu[i]));
    }
    if (rho_x >= one)
        throw NumericError("qsys_mg1_fb: truncated load reaches one, the mean is infinite");
    return numer / ((one - rho_x) * (one - rho_x)) + x / (one - rho_x);
}

}  // namespace detail

/**
 * @param lambda per-class arrival rates
 * @param mu     per-class service rates
 * @param cs     per-class coefficients of variation of the service time
 */
template <class T>
Mg1DisciplineResult<T> qsys_mg1_fb(const std::vector<T>& lambda, const std::vector<T>& mu,
                                   const std::vector<T>& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mg1_fb requires transcendental arithmetic");
    detail::mg1_discipline_check("qsys_mg1_fb", lambda, mu, cs);
    const std::size_t K = lambda.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T cs_tol = T(num_traits<T>::from_double(1e-6));

    bool all_exp = true;
    for (std::size_t i = 0; i < K; ++i)
        if (!(num_abs(T(cs[i] - one)) < cs_tol)) all_exp = false;

    Mg1DisciplineResult<T> r;
    r.W.assign(K, zero);

    if (all_exp) {
        T lambda_total = zero;
        for (const T& v : lambda) lambda_total += v;
        std::vector<T> p(K);
        for (std::size_t i = 0; i < K; ++i) p[i] = lambda[i] / lambda_total;
        const T reltol = T(num_traits<T>::from_double(1e-8));
        const T abstol = T(num_traits<T>::from_double(1e-10));
        for (std::size_t k = 0; k < K; ++k) {
            const T mu_k = mu[k];
            const T x_max = num_traits<T>::from_int(20) / mu_k;
            r.W[k] = detail::num_integral<T>(
                [&](const T& x) {
                    return detail::mg1_fb_response(x, mu, p, lambda_total) * mu_k *
                           detail::num_exp(T(-mu_k * x));
                },
                zero, x_max, reltol, abstol);
        }
    } else {
        for (std::size_t k = 0; k < K; ++k) {
            const T x = one / mu[k];
            T rho_x = zero, numer = zero;
            for (std::size_t i = 0; i < K; ++i) {
                T int_Fbar, int_tFbar;
                if (num_abs(T(cs[i] - one)) < cs_tol) {
                    const T e = detail::num_exp(T(-mu[i] * x));
                    int_Fbar = (one - e) / mu[i];
                    int_tFbar = (one - e * (one + mu[i] * x)) / (mu[i] * mu[i]);
                } else {
                    int_Fbar = detail::num_min(x, T(one / mu[i]));
                    int_tFbar = detail::num_min(T(x * x / two), T(one / (mu[i] * mu[i])));
                }
                rho_x += lambda[i] * int_Fbar;
                numer += lambda[i] * int_tFbar;
            }
            if (rho_x >= one)
                throw NumericError("qsys_mg1_fb: truncated load reaches one, the mean is infinite");
            r.W[k] = numer / ((one - rho_x) * (one - rho_x)) + x / (one - rho_x);
        }
    }
    r.rhohat = detail::mg1_discipline_rhohat(lambda, r.W);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1_FB_H
