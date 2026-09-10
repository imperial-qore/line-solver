/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MTGINF_H
#define LINE_API_QSYS_MTGINF_H

/**
 * Exact time-varying analysis of the Mt/G/infinity queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_mtginf.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mtginf.java.
 *
 * THE RESULT IS EXACT, not an approximation. With infinitely many servers
 * customers never interact, so the model is a Poisson random measure and the
 * number in system at time t is POISSON with mean
 *
 *   m(t) = E[ int_{t-S}^{t} lambda(u) du ] = ES E[lambda(t - Se)]
 *        = int_0^Inf lambda(t-x) P(S > x) dx
 *
 * where Se is the STATIONARY-EXCESS (equilibrium) law of the service time, with
 * density P(S>x)/ES. Because the law is Poisson the variance equals the mean.
 *
 * THE PHYSICS. Reading m(t) as ES E[lambda(t-Se)] says the time-varying load is
 * the stationary load ES lambda(t) subjected to a TIME LAG and a SPACE SHIFT: to
 * first order m(t) ~ ES lambda(t - E[Se]) with E[Se] = E[S^2]/(2 ES), so peak
 * congestion LAGS peak arrival rate, and by more than the mean service time when
 * the service law is variable. The pointwise stationary approximation is the
 * zeroth-order term of the same expansion.
 *
 * ARITHMETIC. The age integral is a Simpson quadrature against a tail cut, so
 * the answer is a quadrature approximation whatever the arithmetic; the
 * instantiation is restricted to the transcendental types.
 *
 * Reference: S. G. Eick, W. A. Massey, W. Whitt (1993). The physics of the
 * Mt/G/infinity queue. Operations Research 41(4), 731-742.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Time-varying measures of the Mt/G/infinity queue. */
template <class T>
struct QsysMtginfResult {
    std::vector<T> times;             ///< the evaluation times
    std::vector<T> meanNumber;        ///< m(t), the Poisson mean
    std::vector<T> varNumber;         ///< equal to meanNumber, the law being Poisson
    std::vector<T> arrivalRate;       ///< lambda(t)
    std::vector<T> departureRate;     ///< delta(t) = E[lambda(t-S)]
    std::vector<T> offeredLoadPSA;    ///< ES lambda(t), the pointwise stationary approximation
    T meanLag;                        ///< E[Se], set only when ES2 was supplied
    std::vector<T> lagApproximation;  ///< ES lambda(t-E[Se]), likewise
    bool hasLag = false;              ///< whether the two fields above are set
};

namespace detail {

/** Nodes and weights of the composite Simpson rule on an even panel count. */
template <class T>
void mtginf_simpson(const T& a, const T& b, std::size_t n, std::vector<T>& x, std::vector<T>& w) {
    if (n % 2 == 1) ++n;
    if (b <= a) {
        x.assign(1, a);
        w.assign(1, num_traits<T>::from_int(0));
        return;
    }
    const T h = (b - a) / num_traits<T>::from_int(static_cast<long>(n));
    x.resize(n + 1);
    w.resize(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        x[i] = a + num_traits<T>::from_int(static_cast<long>(i)) * h;
        const long c = (i == 0 || i == n) ? 1 : (i % 2 == 1 ? 4 : 2);
        w[i] = num_traits<T>::from_int(c) * h / num_traits<T>::from_int(3);
    }
}

/** Smallest doubling point at which the service ccdf is below tol. */
template <class T, class Ccdf>
T mtginf_tail_cut(Ccdf&& ccdf, double tol, double cap) {
    T x = num_traits<T>::from_int(1);
    const T tolT = num_traits<T>::from_double(tol);
    const T capT = num_traits<T>::from_double(cap);
    while (ccdf(x) > tolT) {
        x *= num_traits<T>::from_int(2);
        if (x > capT) return capT;
    }
    return x;
}

/**
 * The Poisson mean m(t), shared by the public entry point and the
 * finite-difference departure rate so that neither re-derives the other. With an
 * infinite past the age grid does not move with t, so the service ccdf is
 * evaluated once rather than once per time point.
 */
template <class T, class Lam, class Ccdf>
std::vector<T> mtginf_mean(Lam&& lambdaFun, Ccdf&& serviceCcdf, const std::vector<T>& t,
                           double startTime, const T& cut, std::size_t panels, bool unbounded) {
    std::vector<T> xs, ws, gcs;
    if (unbounded) {
        mtginf_simpson(num_traits<T>::from_int(0), cut, panels, xs, ws);
        gcs.resize(xs.size());
        for (std::size_t j = 0; j < xs.size(); ++j) gcs[j] = serviceCcdf(xs[j]);
    }
    std::vector<T> m(t.size(), num_traits<T>::from_int(0));
    std::vector<T> x, w, gc;
    for (std::size_t i = 0; i < t.size(); ++i) {
        const std::vector<T>*px;
        const std::vector<T>*pw;
        const std::vector<T>*pgc;
        if (unbounded) {
            px = &xs;
            pw = &ws;
            pgc = &gcs;
        } else {
            T hi = t[i] - num_traits<T>::from_double(startTime);
            if (hi < num_traits<T>::from_int(0)) hi = num_traits<T>::from_int(0);
            if (hi > cut) hi = cut;
            mtginf_simpson(num_traits<T>::from_int(0), hi, panels, x, w);
            gc.resize(x.size());
            for (std::size_t j = 0; j < x.size(); ++j) gc[j] = serviceCcdf(x[j]);
            px = &x;
            pw = &w;
            pgc = &gc;
        }
        T acc = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < px->size(); ++j)
            // m(t) = int lambda(t-x) P(S>x) dx: arrivals of age x still in service.
            acc += (*pw)[j] * lambdaFun(T(t[i] - (*px)[j])) * (*pgc)[j];
        m[i] = acc;
    }
    return m;
}

}  // namespace detail

/**
 * @param lambdaFun   the arrival rate; must accept arguments in the past when
 *                    startTime is infinite
 * @param serviceCcdf G^c(x) = P(S > x)
 * @param ES          the mean service time
 * @param tvals       the times at which to evaluate
 * @param startTime   time the system started empty; -Inf assumes an infinite past
 * @param ES2         the second moment of the service time; NaN skips the lag
 * @param servicePdf  the service density for the exact departure rate, or empty
 * @param tol         service-tail cut for the age integral
 * @param panels      Simpson panels for that integral
 * @param maxAge      cap on the age integrated over
 */
template <class T>
QsysMtginfResult<T> qsys_mtginf(
    const std::function<T(const T&)>& lambdaFun, const std::function<T(const T&)>& serviceCcdf,
    const T& ES, const std::vector<T>& tvals,
    double startTime = -std::numeric_limits<double>::infinity(),
    double ES2 = std::numeric_limits<double>::quiet_NaN(),
    const std::function<T(const T&)>& servicePdf = std::function<T(const T&)>(),
    double tol = 1e-12, std::size_t panels = 4000, double maxAge = 1e12) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mtginf integrates against a tail cut, so it needs inexact arithmetic");
    if (ES <= num_traits<T>::from_int(0))
        throw InputError("qsys_mtginf: the mean service time ES must be positive");
    if (!lambdaFun || !serviceCcdf)
        throw InputError("qsys_mtginf: the arrival rate and the service ccdf must be callable");

    const T cut = detail::mtginf_tail_cut<T>(serviceCcdf, tol, maxAge);
    const bool unbounded = std::isinf(startTime);

    QsysMtginfResult<T> res;
    res.times = tvals;
    res.meanNumber = detail::mtginf_mean<T>(lambdaFun, serviceCcdf, tvals, startTime, cut, panels,
                                            unbounded);
    res.varNumber = res.meanNumber;         // Poisson: the variance is the mean
    res.arrivalRate.resize(tvals.size());
    res.offeredLoadPSA.resize(tvals.size());
    for (std::size_t i = 0; i < tvals.size(); ++i) {
        res.arrivalRate[i] = lambdaFun(tvals[i]);
        res.offeredLoadPSA[i] = ES * res.arrivalRate[i];
    }

    res.departureRate.resize(tvals.size());
    if (servicePdf) {
        std::vector<T> x, w;
        for (std::size_t i = 0; i < tvals.size(); ++i) {
            T hi = cut;
            if (!unbounded) {
                hi = tvals[i] - num_traits<T>::from_double(startTime);
                if (hi < num_traits<T>::from_int(0)) hi = num_traits<T>::from_int(0);
                if (hi > cut) hi = cut;
            }
            detail::mtginf_simpson(num_traits<T>::from_int(0), hi, panels, x, w);
            T acc = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < x.size(); ++j)
                acc += w[j] * lambdaFun(T(tvals[i] - x[j])) * servicePdf(x[j]);
            res.departureRate[i] = acc;
        }
    } else {
        // Flow balance m'(t) = lambda(t) - delta(t), differentiated centrally.
        double tmax = 1.0;
        for (const T& v : tvals) tmax = std::max(tmax, std::abs(num_traits<T>::to_double(v)));
        const T h = num_traits<T>::from_double(1e-5 * tmax);
        std::vector<T> tu(tvals.size()), td(tvals.size());
        for (std::size_t i = 0; i < tvals.size(); ++i) {
            tu[i] = tvals[i] + h;
            td[i] = tvals[i] - h;
        }
        const std::vector<T> up =
            detail::mtginf_mean<T>(lambdaFun, serviceCcdf, tu, startTime, cut, panels, unbounded);
        const std::vector<T> dn =
            detail::mtginf_mean<T>(lambdaFun, serviceCcdf, td, startTime, cut, panels, unbounded);
        for (std::size_t i = 0; i < tvals.size(); ++i)
            res.departureRate[i] =
                res.arrivalRate[i] - (up[i] - dn[i]) / (num_traits<T>::from_int(2) * h);
    }

    if (!std::isnan(ES2)) {
        res.hasLag = true;
        res.meanLag = num_traits<T>::from_double(ES2) / (num_traits<T>::from_int(2) * ES);
        res.lagApproximation.resize(tvals.size());
        for (std::size_t i = 0; i < tvals.size(); ++i)
            res.lagApproximation[i] = ES * lambdaFun(T(tvals[i] - res.meanLag));
    } else {
        res.meanLag = num_traits<T>::from_int(0);
    }
    return res;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MTGINF_H
