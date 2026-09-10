/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1_PSJF_H
#define LINE_API_QSYS_QSYS_MG1_PSJF_H

/**
 * M/G/1 under PSJF (preemptive shortest job first).
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1_psjf.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mg1_psjf.java.
 *
 * Priority follows the original size, so a job of size x is delayed only by
 * work of size at most x (Wierman and Harchol-Balter, SIGMETRICS 2003,
 * Sec. 3.2):
 *
 *   rho(x)  = lambda int_0^x t f(t) dt
 *   m2(x)   = lambda int_0^x t^2 f(t) dt
 *   E[T(x)] = x/(1-rho(x)) + m2(x)/(2 (1-rho(x))^2)
 *
 * As in qsys_mg1_fb, MATLAB has an exponential path -- closed-form truncated
 * moments of the mixture, class mean by quadrature over [0, 20/mu_k] -- and a
 * general path that sorts the classes by increasing mean size and evaluates
 * the same expression at x = 1/mu_k with the truncated moments replaced by the
 * per-class second moments of the classes at least as fast. Both are
 * reproduced.
 *
 * ARITHMETIC. exp and the adaptive quadrature make this transcendental.
 */

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

#include "line/api/qsys/qsys_mg1_setf.h"
#include "line/api/qsys/qsys_quadrature.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

namespace detail {

/** E[T(x)] under PSJF for a mixture-of-exponentials job-size law. */
template <class T>
T mg1_psjf_response(const T& x, const std::vector<T>& mu, const std::vector<T>& p,
                    const T& lambda_total) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    T m1_x = num_traits<T>::from_int(0), m2_x = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < mu.size(); ++i) {
        const T mi = mu[i];
        const T e = num_exp(T(-mi * x));
        m1_x += p[i] * (one / mi - (one / mi + x) * e);
        m2_x += p[i] * (two / (mi * mi) - (two / (mi * mi) + two * x / mi + x * x) * e);
    }
    const T rho_x = lambda_total * m1_x;
    if (rho_x >= one)
        throw NumericError("qsys_mg1_psjf: truncated load reaches one, the mean is infinite");
    return x / (one - rho_x) + lambda_total * m2_x / (two * (one - rho_x) * (one - rho_x));
}

}  // namespace detail

/**
 * @param lambda per-class arrival rates
 * @param mu     per-class service rates
 * @param cs     per-class coefficients of variation of the service time
 */
template <class T>
Mg1DisciplineResult<T> qsys_mg1_psjf(const std::vector<T>& lambda, const std::vector<T>& mu,
                                     const std::vector<T>& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mg1_psjf requires transcendental arithmetic");
    detail::mg1_discipline_check("qsys_mg1_psjf", lambda, mu, cs);
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
                    return detail::mg1_psjf_response(x, mu, p, lambda_total) * mu_k *
                           detail::num_exp(T(-mu_k * x));
                },
                zero, x_max, reltol, abstol);
        }
    } else {
        // Classes sorted by increasing mean service time; MATLAB's sort is
        // stable, so ties keep the original class order.
        std::vector<std::size_t> idx(K);
        std::iota(idx.begin(), idx.end(), std::size_t(0));
        std::stable_sort(idx.begin(), idx.end(),
                         [&](std::size_t i, std::size_t j) { return one / mu[i] < one / mu[j]; });
        T rho_cum = zero;
        for (std::size_t pos = 0; pos < K; ++pos) {
            const std::size_t k = idx[pos];
            const T x = one / mu[k];
            rho_cum += lambda[k] / mu[k];
            T m2_x = zero;
            for (std::size_t q = 0; q <= pos; ++q) {
                const std::size_t i = idx[q];
                m2_x += lambda[i] * (one + cs[i] * cs[i]) / (mu[i] * mu[i]);
            }
            if (rho_cum >= one)
                throw NumericError(
                    "qsys_mg1_psjf: truncated load reaches one, the mean is infinite");
            r.W[k] = m2_x / (two * (one - rho_cum) * (one - rho_cum)) + x / (one - rho_cum);
        }
    }
    r.rhohat = detail::mg1_discipline_rhohat(lambda, r.W);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1_PSJF_H
