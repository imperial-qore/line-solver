/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1_LRPT_H
#define LINE_API_QSYS_QSYS_MG1_LRPT_H

/**
 * M/G/1 under LRPT (longest remaining processing time).
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1_lrpt.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mg1_lrpt.java.
 *
 * LRPT gives the server to the job with the most work left, so the slowdown is
 * the same for every job in the system (Wierman and Harchol-Balter,
 * SIGMETRICS 2003, Sec. 3.2):
 *
 *   E[T(x)] = x/(1-rho) + lambda E[X^2] / (2 (1-rho)^2)
 *
 * with rho the total load and E[X^2] the second moment of the mixture size.
 * Only the first term depends on x, so under the exponential path the class
 * mean is
 *
 *   E[T_k] = int_0^{20/mu_k} E[T(x)] mu_k e^{-mu_k x} dx
 *
 * and this integral is elementary. MATLAB evaluates it by adaptive quadrature;
 * the port uses the closed form of the same truncated integral,
 *
 *   int_0^X x mu e^{-mu x} dx = (1 - e^{-muX}(1 + muX))/mu,
 *   int_0^X mu e^{-mu x} dx   = 1 - e^{-muX},        muX = 20,
 *
 * which is the value MATLAB's quadrature converges to, so the two agree to
 * MATLAB's RelTol of 1e-8 and in practice to round-off. The truncation at 20
 * mean service times is kept: dropping it would change the answer in the
 * eighth digit and no longer reproduce the reference.
 *
 * The general path (some cs != 1) is MATLAB's class-based preemptive-priority
 * surrogate, with the classes sorted by decreasing mean size. Note that it
 * ignores cs entirely -- the residual term is sum_{i<=k} lambda_i/mu_i^2,
 * i.e. the second moment of an exponential -- so for a non-exponential input
 * the reference answer does not depend on the variability it was given. That
 * is a defect of the reference, reproduced here rather than silently fixed.
 *
 * ARITHMETIC. exp appears in the exponential path, so the function is gated.
 */

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

#include "line/api/qsys/qsys_mg1_setf.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/**
 * @param lambda per-class arrival rates
 * @param mu     per-class service rates
 * @param cs     per-class coefficients of variation of the service time
 */
template <class T>
Mg1DisciplineResult<T> qsys_mg1_lrpt(const std::vector<T>& lambda, const std::vector<T>& mu,
                                     const std::vector<T>& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mg1_lrpt requires transcendental arithmetic");
    const T rho_total = detail::mg1_discipline_check("qsys_mg1_lrpt", lambda, mu, cs);
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
        T E_X2 = zero;
        for (std::size_t i = 0; i < K; ++i)
            E_X2 += (lambda[i] / lambda_total) * two / (mu[i] * mu[i]);
        const T flat = lambda_total * E_X2 / (two * (one - rho_total) * (one - rho_total));
        // muX = 20 for every class, so the two truncation factors are shared.
        const T e20 = detail::num_exp(T(num_traits<T>::from_int(-20)));
        const T mass = one - e20;
        const T first = one - e20 * num_traits<T>::from_int(21);
        for (std::size_t k = 0; k < K; ++k)
            r.W[k] = (first / mu[k]) / (one - rho_total) + flat * mass;
    } else {
        // Classes sorted by decreasing mean service time; MATLAB's sort is
        // stable, so ties keep the original class order.
        std::vector<std::size_t> idx(K);
        std::iota(idx.begin(), idx.end(), std::size_t(0));
        std::stable_sort(idx.begin(), idx.end(),
                         [&](std::size_t i, std::size_t j) { return one / mu[j] < one / mu[i]; });
        T rho_prev = zero, E_R_k = zero;
        for (std::size_t pos = 0; pos < K; ++pos) {
            const std::size_t k = idx[pos];
            const T rho_curr = rho_prev + lambda[k] / mu[k];
            E_R_k += lambda[k] / (mu[k] * mu[k]);
            r.W[k] = E_R_k / ((one - rho_prev) * (one - rho_curr)) + one / mu[k];
            rho_prev = rho_curr;
        }
    }
    r.rhohat = detail::mg1_discipline_rhohat(lambda, r.W);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1_LRPT_H
