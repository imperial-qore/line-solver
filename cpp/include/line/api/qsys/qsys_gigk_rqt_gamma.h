/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_GIGK_RQT_GAMMA_H
#define LINE_API_QSYS_QSYS_GIGK_RQT_GAMMA_H

/**
 * Service variability parameter of the Robust Queueing Theory (RQT) framework.
 *
 * Templated port of matlab/src/api/qsys/qsys_gigk_rqt_gamma.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_gigk_rqt_gamma.java.
 *
 * The adaptation of Section 7.1 turns the first two moments into an uncertainty
 * set:
 *
 *   Gamma_s = (2 (theta0 + theta1 sigma_s^2/k + theta2 Gamma_a^2 rho^2 k))^((a-1)/a)
 *             - Gamma_a k^((a-1)/a),
 *
 * with (theta0,theta1,theta2) regressed so that the worst-case system time of
 * Theorem 3 approximates the MEAN system time of the corresponding stochastic
 * queue. The arrival side needs no adaptation: Gamma_a = sigma_a for an external
 * renewal stream. Since the last term cancels Gamma_a at alpha = 2, the
 * adaptation acts on the sum Gamma_a + Gamma_s/k^(1/alpha) that Theorem 3 reads.
 *
 * THE FACTOR 2 IS NOT IN THE PRINTED FORMULA and is restored here. Section 7.1
 * states that the functional form is motivated by Kingman's bound, which the
 * alpha=2 bound of Theorem 3 reproduces when (Gamma_a+Gamma_s)^2 =
 * 2(sigma_a^2+sigma_s^2); the published thetas are all near unity, i.e.
 * corrections to that bound rather than a substitute for its factor 2. Dropping
 * the factor puts M/M/1 about 40% BELOW its exact mean system time at rho = 0.9,
 * contradicting the errors of at most 9.5% that Tables 2-3 report; restoring it
 * gives +4.7%.
 *
 * CAUTION: the form is not dimensionally homogeneous, since theta0 is an
 * additive constant on a scale of variances, so it is only valid in the time
 * unit the regression was run in. It is evaluated here in units of the mean
 * service time, 1/mu = 1, and converted back.
 *
 * ARITHMETIC. A real exponent makes this transcendental.
 *
 * Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
 * Operations Research 63(3), 676-700, Section 7.1 and Table 1.
 */

#include <cstddef>
#include <string>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/**
 * @param rho     traffic intensity lambda/(k mu)
 * @param mu      service rate of each server, which sets the time unit
 * @param Gamma_a variability parameter of the arrival uncertainty set
 * @param sigma_s standard deviation of the service time
 * @param k       number of servers
 * @param alpha_a effective arrival tail coefficient in (1,2]
 * @param regime  adaptation regime of Table 1: "independent" (service
 *                distribution unknown), "normal" or "pareto"
 */
template <class T>
T qsys_gigk_rqt_gamma(const T& rho, const T& mu, const T& Gamma_a, const T& sigma_s, std::size_t k,
                      const T& alpha_a, const std::string& regime = "independent") {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gigk_rqt_gamma requires transcendental arithmetic");
    double t0 = 0.0, t1 = 0.0, t2 = 0.0;
    if (regime == "pareto") {
        t0 = -0.05; t1 = 1.09; t2 = 1.11;
    } else if (regime == "normal") {
        t0 = -0.02; t1 = 1.03; t2 = 1.04;
    } else if (regime == "independent" || regime == "default" || regime.empty()) {
        t0 = -0.06; t1 = 1.07; t2 = 1.07;
    } else {
        throw InputError("qsys_gigk_rqt_gamma: unknown RQT adaptation regime: " + regime);
    }
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T kk = num_traits<T>::from_int(static_cast<int>(k));

    // evaluate in units of the mean service time, then convert back
    const T ga = Gamma_a * mu;
    const T ss = sigma_s * mu;
    const T e = (alpha_a - one) / alpha_a;
    T b = two * (num_traits<T>::from_double(t0) + num_traits<T>::from_double(t1) * ss * ss / kk +
                 num_traits<T>::from_double(t2) * ga * ga * rho * rho * kk);
    if (b < zero) b = zero;
    const T gs = detail::num_pow(b, e) - ga * detail::num_pow(kk, e);
    return gs / mu;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_GIGK_RQT_GAMMA_H
