/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1_SRPT_H
#define LINE_API_QSYS_QSYS_MG1_SRPT_H

/**
 * M/G/1 under SRPT (shortest remaining processing time), by the
 * Schrage-Miller formula.
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1_srpt.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mg1_srpt.java.
 *
 * For a job of size x, with f the mixture size density and Fbar its tail,
 *
 *   rho(x)  = lambda int_0^x t f(t) dt
 *   m2(x)   = int_0^x t^2 f(t) dt
 *   E[W(x)] = lambda (m2(x) + x^2 Fbar(x)) / (2 (1-rho(x))^2)
 *   E[R(x)] = int_0^x dt/(1-rho(t))
 *   E[T(x)] = E[W(x)] + E[R(x)]
 *
 * and the class mean is int_0^inf E[T(x)] f_r(x) dx. Because E[T(x)] depends
 * on the size alone -- SRPT is size-based, not class-based -- that integral is
 * exact given the size laws.
 *
 * Each class is matched to (mean 1/mu_r, scv cs_r^2) by an exponential when
 * cs_r = 1, a balanced-means two-phase hyperexponential when cs_r > 1, and a
 * Tijms mixture of Erlang-(k-1)/Erlang-k when cs_r < 1. The integrals are
 * evaluated by cumulative trapezoid quadrature on the same fixed grid MATLAB
 * builds: 40 e-foldings of the slowest phase, and at least 20000 points, or
 * 200 per unit of phase-rate spread. That grid, not the formula, sets the
 * accuracy -- the trapezoid rule on a uniform grid of N points is O(N^-2), so
 * the reference itself carries an error of order 1e-6 relative, and the port
 * matches it point for point rather than integrating better.
 *
 * ARITHMETIC. exp, log and the trapezoid rule make this transcendental.
 */

#include <algorithm>
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

/** Job-size law matched to a mean and an SCV, as in MATLAB's local srpt_fit. */
template <class T>
struct SrptFit {
    enum Kind { EXPONENTIAL, HYPEREXP2, ERLANG_MIX } kind = EXPONENTIAL;
    T rate;   ///< EXPONENTIAL and ERLANG_MIX
    T p;      ///< branch probability for HYPEREXP2 and ERLANG_MIX
    T r1, r2; ///< HYPEREXP2 phase rates
    unsigned k = 0;  ///< ERLANG_MIX larger shape
    T rate_min, rate_max;
};

template <class T>
SrptFit<T> srpt_fit(const T& mu, const T& cs) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T half = num_traits<T>::from_rational(1, 2);
    const T c2 = cs * cs;
    SrptFit<T> f;
    if (num_abs(T(c2 - one)) < T(num_traits<T>::from_double(1e-9))) {
        f.kind = SrptFit<T>::EXPONENTIAL;
        f.rate = mu;
        f.rate_min = mu;
        f.rate_max = mu;
    } else if (c2 > one) {
        const T pr = half * (one + num_sqrt(T((c2 - one) / (c2 + one))));
        f.kind = SrptFit<T>::HYPEREXP2;
        f.p = pr;
        f.r1 = two * pr * mu;
        f.r2 = two * (one - pr) * mu;
        f.rate_min = num_min(f.r1, f.r2);
        f.rate_max = f.r1 < f.r2 ? f.r2 : f.r1;
    } else {
        const double inv = 1.0 / num_traits<T>::to_double(c2);
        const unsigned k = static_cast<unsigned>(std::ceil(inv));
        const T kt = num_traits<T>::from_int(static_cast<long>(k));
        const T pr = (one / (one + c2)) * (kt * c2 - num_sqrt(T(kt * (one + c2) - kt * kt * c2)));
        f.kind = SrptFit<T>::ERLANG_MIX;
        f.k = k;
        f.p = pr;
        f.rate = (kt - pr) * mu;  // mean_x = 1/mu, so (k-p)/mean_x = (k-p) mu
        f.rate_min = f.rate;
        f.rate_max = f.rate;
    }
    return f;
}

/** log(m!) by summation, standing in for MATLAB's gammaln(m+1). */
template <class T>
T log_factorial(unsigned m) {
    T s = num_traits<T>::from_int(0);
    for (unsigned j = 2; j <= m; ++j) {
        using std::log;
        s += log(num_traits<T>::from_int(static_cast<long>(j)));
    }
    return s;
}

/** Erlang-n density with the given rate, in log space as in MATLAB. */
template <class T>
T erlang_pdf(unsigned n, const T& rate, const T& x) {
    const T zero = num_traits<T>::from_int(0);
    if (n == 0) return zero;
    const T t = rate * x;
    const unsigned m = n - 1;
    if (t <= zero) return m == 0 ? rate : zero;
    using std::log;
    const T logp = num_traits<T>::from_int(static_cast<long>(m)) * log(t) - t - log_factorial<T>(m);
    return rate * num_exp(logp);
}

/** Erlang-n complementary CDF, the upper Poisson tail, in log space. */
template <class T>
T erlang_tail(unsigned n, const T& rate, const T& x) {
    const T zero = num_traits<T>::from_int(0);
    if (n == 0) return zero;
    const T t = rate * x;
    if (t <= zero) return num_traits<T>::from_int(1);
    using std::log;
    T y = zero;
    const T logt = log(t);
    for (unsigned j = 0; j < n; ++j)
        y += num_exp(T(num_traits<T>::from_int(static_cast<long>(j)) * logt - t -
                       log_factorial<T>(j)));
    return y;
}

template <class T>
T srpt_pdf(const SrptFit<T>& f, const T& x) {
    const T one = num_traits<T>::from_int(1);
    switch (f.kind) {
        case SrptFit<T>::EXPONENTIAL:
            return f.rate * num_exp(T(-f.rate * x));
        case SrptFit<T>::HYPEREXP2:
            return f.p * f.r1 * num_exp(T(-f.r1 * x)) +
                   (one - f.p) * f.r2 * num_exp(T(-f.r2 * x));
        default:
            return f.p * erlang_pdf(f.k - 1, f.rate, x) +
                   (one - f.p) * erlang_pdf(f.k, f.rate, x);
    }
}

template <class T>
T srpt_tail(const SrptFit<T>& f, const T& x) {
    const T one = num_traits<T>::from_int(1);
    switch (f.kind) {
        case SrptFit<T>::EXPONENTIAL:
            return num_exp(T(-f.rate * x));
        case SrptFit<T>::HYPEREXP2:
            return f.p * num_exp(T(-f.r1 * x)) + (one - f.p) * num_exp(T(-f.r2 * x));
        default:
            return f.p * erlang_tail(f.k - 1, f.rate, x) +
                   (one - f.p) * erlang_tail(f.k, f.rate, x);
    }
}

}  // namespace detail

/**
 * @param lambda per-class arrival rates
 * @param mu     per-class service rates
 * @param cs     per-class coefficients of variation of the service time
 */
template <class T>
Mg1DisciplineResult<T> qsys_mg1_srpt(const std::vector<T>& lambda, const std::vector<T>& mu,
                                     const std::vector<T>& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mg1_srpt requires transcendental arithmetic");
    detail::mg1_discipline_check("qsys_mg1_srpt", lambda, mu, cs);
    const std::size_t K = lambda.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    T lambda_total = zero;
    for (const T& v : lambda) lambda_total += v;
    std::vector<T> p(K);
    for (std::size_t i = 0; i < K; ++i) p[i] = lambda[i] / lambda_total;

    std::vector<detail::SrptFit<T>> fits;
    fits.reserve(K);
    T rate_min = zero, rate_max = zero;
    for (std::size_t r = 0; r < K; ++r) {
        fits.push_back(detail::srpt_fit(mu[r], cs[r]));
        if (r == 0 || fits[r].rate_min < rate_min) rate_min = fits[r].rate_min;
        if (r == 0 || fits[r].rate_max > rate_max) rate_max = fits[r].rate_max;
    }
    if (rate_min <= zero) throw NumericError("qsys_mg1_srpt: degenerate phase rate");

    // Grid: 40 e-foldings of the slowest phase, 200 points per rate ratio.
    const T xmax = num_traits<T>::from_int(40) / rate_min;
    const double ratio = num_traits<T>::to_double(rate_max) / num_traits<T>::to_double(rate_min);
    std::size_t N = static_cast<std::size_t>(std::ceil(200.0 * ratio));
    if (N < 20000u) N = 20000u;
    if (N > 2000000u) N = 2000000u;

    std::vector<T> x(N + 1), fmix(N + 1, zero), Fbar(N + 1, zero);
    for (std::size_t i = 0; i <= N; ++i)
        x[i] = xmax * num_traits<T>::from_int(static_cast<long>(i)) /
               num_traits<T>::from_int(static_cast<long>(N));
    for (std::size_t r = 0; r < K; ++r)
        for (std::size_t i = 0; i <= N; ++i) {
            fmix[i] += p[r] * detail::srpt_pdf(fits[r], x[i]);
            Fbar[i] += p[r] * detail::srpt_tail(fits[r], x[i]);
        }

    std::vector<T> xf(N + 1), x2f(N + 1);
    for (std::size_t i = 0; i <= N; ++i) {
        xf[i] = x[i] * fmix[i];
        x2f[i] = x[i] * x[i] * fmix[i];
    }
    std::vector<T> rho_x = detail::num_cumtrapz(x, xf);
    for (T& v : rho_x) v *= lambda_total;
    const std::vector<T> m2_x = detail::num_cumtrapz(x, x2f);

    // Guard the (1-rho(x)) factors; rho(x) -> rho < 1 as x -> inf.
    const T floor_ = T(num_traits<T>::from_double(1e-12));
    std::vector<T> denom(N + 1), invden(N + 1), ET(N + 1);
    for (std::size_t i = 0; i <= N; ++i) {
        const T d = one - rho_x[i];
        denom[i] = d > floor_ ? d : floor_;
        invden[i] = one / denom[i];
    }
    const std::vector<T> Res = detail::num_cumtrapz(x, invden);
    for (std::size_t i = 0; i <= N; ++i)
        ET[i] = lambda_total * (m2_x[i] + x[i] * x[i] * Fbar[i]) / (two * denom[i] * denom[i]) +
                Res[i];

    Mg1DisciplineResult<T> r;
    r.W.assign(K, zero);
    std::vector<T> integrand(N + 1);
    for (std::size_t c = 0; c < K; ++c) {
        for (std::size_t i = 0; i <= N; ++i) integrand[i] = ET[i] * detail::srpt_pdf(fits[c], x[i]);
        r.W[c] = detail::num_trapz(x, integrand);
    }
    r.rhohat = detail::mg1_discipline_rhohat(lambda, r.W);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1_SRPT_H
