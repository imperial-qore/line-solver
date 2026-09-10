/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_GAMMA_H
#define LINE_API_MAM_MAP_GAMMA_H

/**
 * Autocorrelation decay rate of a MAP: the gamma of the geometric model
 * rho(k) = rho0 * gamma^k with rho0 = (1 - 1/scv)/2.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_gamma.m, cross-checked
 * against jar/src/main/java/jline/api/mam/Map_gamma.java.
 *
 * Order one is Poisson and has decay rate zero. Order two has a genuinely
 * geometric acf, so the rate is the exact ratio acf(2)/acf(1); a vanishing
 * acf(1) means the MAP degenerates to a phase-type renewal process and the rate
 * is again zero. Only above order two is a fit needed, and there both references
 * evaluate the acf on the ten lags 1, 1+limit/10, ..., and regress.
 *
 * WHY THE FIT IS ROBUST AND NOT PLAIN LEAST SQUARES: MATLAB uses nlinfit with
 * RobustWgtFun 'fair'. On a higher-order MAP the acf is a sum of geometric
 * terms, so the short lags sit far off any single geometric and would otherwise
 * dominate the fit and bias gamma low. The reproduction here is the JAR's: an
 * ordinary fit, the leverage of its Jacobian held fixed, then iteratively
 * reweighted fits with w = 1/(1 + |r_adj|/(1.4 sigma)) and sigma the MAD
 * estimate floored against the spread of the response. Dropping the robustness
 * is not a simplification, it changes the answer.
 *
 * ARITHMETIC: unlike most of api/mam this is a floating-point algorithm. The
 * lags enter as real exponents and the MAD needs an order statistic, so the
 * regression runs in double whatever T is; only the acf, mean and second moment
 * that feed it are computed in T. Orders one and two return an exact T.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/levmar.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace detail {

/** One weighted Levenberg-Marquardt fit of rho_k = rho0 gamma^k. */
inline double map_gamma_ls(const std::vector<double>& lag, const std::vector<double>& rho,
                           double rho0, double start, const std::vector<double>& weights) {
    const std::size_t m = lag.size();
    std::vector<double> sw(m, 1.0);
    if (!weights.empty())
        for (std::size_t i = 0; i < m; ++i) sw[i] = std::sqrt(weights[i]);
    // The optimizer minimises sqrt(w) * (model - rho), so both sides carry sqrt(w)
    const std::vector<double>* plag = &lag;
    const std::vector<double>* prho = &rho;
    auto res = [plag, prho, &sw, rho0, m](const std::vector<double>& x) {
        std::vector<double> r(m);
        for (std::size_t i = 0; i < m; ++i)
            r[i] = sw[i] * (rho0 * std::pow(x[0], (*plag)[i]) - (*prho)[i]);
        return r;
    };
    auto jac = [plag, &sw, rho0, m](const std::vector<double>& x) {
        Matrix<double> J(m, 1, 0.0);
        for (std::size_t i = 0; i < m; ++i)
            J(i, 0) = sw[i] * rho0 * (*plag)[i] * std::pow(x[0], (*plag)[i] - 1.0);
        return J;
    };
    LevmarOptions<double> opt = levmar_defaults<double>();
    opt.max_iter = 500;
    const LevmarResult<double> out = levmar_jac(res, jac, std::vector<double>(1, start), m, opt);
    return out.x[0];
}

/** Median of a copy of the sample, the MAD ingredient of the robust fit. */
inline double map_gamma_median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/** Robust nonlinear fit of the geometric acf model, nlinfit with the fair weight. */
inline double map_gamma_fit(const std::vector<double>& lag, const std::vector<double>& rho,
                            double rho0, double start) {
    const std::size_t m = lag.size();
    double gamma = map_gamma_ls(lag, rho, rho0, start, std::vector<double>());

    // Leverage of the least-squares Jacobian, held fixed across the reweighting
    // as nlinfit does; for one parameter the QR reduces to normalising it.
    double norm2 = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        const double j = rho0 * lag[i] * std::pow(gamma, lag[i] - 1.0);
        norm2 += j * j;
    }
    std::vector<double> adjust(m);
    for (std::size_t i = 0; i < m; ++i) {
        const double j = rho0 * lag[i] * std::pow(gamma, lag[i] - 1.0);
        const double h = norm2 > 0.0 ? std::min(0.9999, j * j / norm2) : 0.0;
        adjust[i] = 1.0 / std::sqrt(1.0 - h);
    }

    // A near-perfect fit drives the MAD to zero and makes every point an
    // outlier, so the scale is floored against the spread of the response
    double mean = 0.0;
    for (std::size_t i = 0; i < m; ++i) mean += rho[i];
    mean /= static_cast<double>(m);
    double var = 0.0;
    for (std::size_t i = 0; i < m; ++i) var += (rho[i] - mean) * (rho[i] - mean);
    var /= static_cast<double>(m - 1);
    double tiny = 1e-6 * std::sqrt(var);
    if (tiny == 0.0) tiny = 1.0;

    const double tune = 1.4;  // fair
    const double delta = std::sqrt(std::numeric_limits<double>::epsilon());
    std::vector<double> weights(m);
    for (unsigned iter = 0; iter < 200; ++iter) {
        const double previous = gamma;
        std::vector<double> radj(m), absr(m);
        for (std::size_t i = 0; i < m; ++i) {
            radj[i] = (rho[i] - rho0 * std::pow(gamma, lag[i])) * adjust[i];
            absr[i] = std::fabs(radj[i]);
        }
        const double sigma = map_gamma_median(absr) / 0.6745;
        const double scale = std::max(sigma, tiny) * tune;
        for (std::size_t i = 0; i < m; ++i) weights[i] = 1.0 / (1.0 + std::fabs(radj[i] / scale));
        gamma = map_gamma_ls(lag, rho, rho0, previous, weights);
        if (std::fabs(gamma - previous) <
            delta * std::max(std::fabs(gamma), std::fabs(previous)))
            break;
    }
    return gamma;
}

}  // namespace detail

/** Result of map_gamma, mirroring the MATLAB [GAMMA, RHO0] pair. */
template <class T>
struct MapGammaResult {
    T gamma;  ///< the fitted geometric decay rate of the acf
    T rho0;   ///< (1 - 1/scv)/2, the lag-0 amplitude; zero below order three
};

/**
 * @param m     the MAP
 * @param limit largest lag considered; the reference default is 1000
 */
template <class T>
MapGammaResult<T> map_gamma_full(const Map<T>& m, long limit = 1000) {
    if (limit < 1) throw InputError("map_gamma: the lag limit must be positive");
    const T zero = num_traits<T>::from_int(0);
    MapGammaResult<T> out;
    out.gamma = zero;
    out.rho0 = zero;
    const std::size_t n = m.order();
    if (n == 1) return out;
    if (n == 2) {
        const std::vector<T> a = map_acf(m, std::vector<unsigned>{1u, 2u});
        if (num_abs(a[0]) < num_traits<T>::from_double(1e-8)) return out;
        out.gamma = a[1] / a[0];
        return out;
    }
    const long step = std::max<long>(1, limit / 10);
    std::vector<unsigned> lags;
    for (long l = 1; l <= limit; l += step) lags.push_back(static_cast<unsigned>(l));
    const T m1 = map_mean(m);
    const T m2 = map_moment(m, 2u);
    const T scv = (m2 - m1 * m1) / (m1 * m1);
    if (scv == zero) throw NumericError("map_gamma: the MAP has zero scv, rho0 is undefined");
    out.rho0 = num_traits<T>::from_rational(1, 2) * (num_traits<T>::from_int(1) -
                                                     num_traits<T>::from_int(1) / scv);
    const std::vector<T> acf = map_acf(m, lags);
    std::vector<double> dl(lags.size()), dr(lags.size());
    for (std::size_t i = 0; i < lags.size(); ++i) {
        dl[i] = static_cast<double>(lags[i]);
        dr[i] = num_traits<T>::to_double(acf[i]);
    }
    out.gamma = num_traits<T>::from_double(
        detail::map_gamma_fit(dl, dr, num_traits<T>::to_double(out.rho0), 0.99));
    return out;
}

/** Autocorrelation decay rate of a MAP (map_gamma.m). */
template <class T>
T map_gamma(const Map<T>& m, long limit = 1000) {
    return map_gamma_full(m, limit).gamma;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_GAMMA_H
