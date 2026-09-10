/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_HYPEREXP_FIT_LONGTAIL_H
#define LINE_API_MAM_HYPEREXP_FIT_LONGTAIL_H

/**
 * Fitting a hyperexponential to a long-tail distribution.
 *
 * Templated port of matlab/src/api/mam/hyperexp_fit_longtail.m, cross-checked
 * against jar/src/main/java/jline/api/mam/HyperexpFitLongtail.java.
 *
 * WHY MOMENTS ARE THE WRONG HANDLE. A Pareto law with tail index below 2 has
 * infinite variance, so no two- or three-moment fit exists at all; and even when
 * the moments are finite, matching them says nothing about the several ORDERS OF
 * MAGNITUDE of time scale over which a long-tail law acts. This procedure
 * matches the CCDF ITSELF at points spread across those decades.
 *
 * THE RECURSION, with lambda_1 < ... < lambda_k. In the far tail only the
 * slowest component survives, so it can be fitted there alone:
 *
 *   lambda_1 = ln(F^c(c_1)/F^c(b c_1))/((b-1)c_1)                    (4.4)
 *   p_1      = F^c(c_1) exp(lambda_1 c_1)                            (4.5)
 *
 * subtract it and repeat one decade lower (4.6)-(4.11); the last component takes
 * the remaining probability, p_k = 1 - sum_{j<k} p_j, and its rate follows from
 * the ccdf at c_k (4.12)-(4.14). This is Prony's method applied to a ccdf.
 *
 * DEFAULTS. (b, decade) = (1.5, 4) rather than the paper's illustrative (2, 10):
 * the fit is exact AT the fitting arguments and free between them, and measured
 * on a Weibull(0.3) the tighter grid cuts the worst between-point error from
 * about 54% to 12%, at the cost of more components.
 *
 * ARITHMETIC. Logarithms and exponentials throughout: transcendental only.
 *
 * Reference: A. Feldmann, W. Whitt (1998). Fitting mixtures of exponentials to
 * long-tail distributions to analyze network performance models. Performance
 * Evaluation 31, 245-279, Section 4.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** Outcome of the long-tail hyperexponential fit. */
template <class T>
struct HyperexpLongtailResult {
    std::vector<T> p;        ///< mixing probabilities, summing to 1
    std::vector<T> lambda;   ///< rates, increasing
    std::vector<T> points;   ///< the fitting arguments c_i
    T mean;                  ///< mean of the fitted law
    T targetMean;            ///< mean of the original law over the covered range
    T coverageLow;           ///< c_k, the smallest constrained argument
    T coverageHigh;          ///< b c_1, the largest
    T maxRelError;           ///< worst relative error at the fitting arguments
    T maxRelErrorGrid;       ///< worst relative error on a log grid across the coverage
};

namespace detail {

/** Smallest t with F^c(t) <= prob, by doubling then bisection. */
template <class T, class Ccdf>
T hefit_quantile(Ccdf&& ccdf, const T& prob) {
    const T two = num_traits<T>::from_int(2);
    T hi = num_traits<T>::from_int(1);
    while (ccdf(hi) > prob) {
        hi *= two;
        if (hi > num_traits<T>::from_double(1e15))
            throw InputError("hyperexp_fit_longtail: the ccdf does not decay, so there is no tail "
                             "to fit");
    }
    T lo = num_traits<T>::from_int(0);
    for (int i = 0; i < 200; ++i) {
        const T mid = (lo + hi) / two;
        if (ccdf(mid) > prob) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return (lo + hi) / two;
}

}  // namespace detail

/**
 * The recursion at a fixed component count.
 *
 * @param ccdf   F^c(t) = P(X > t)
 * @param k      number of exponential components
 * @param c1     the largest fitting argument
 * @param b      the within-scale spacing, 1 < b < decade
 * @param decade the ratio between successive fitting arguments
 */
template <class T, class Ccdf>
HyperexpLongtailResult<T> hyperexp_fit_longtail_k(Ccdf&& ccdf, std::size_t k, const T& c1,
                                                  const T& b, const T& decade) {
    static_assert(num_traits<T>::has_transcendental,
                  "hyperexp_fit_longtail needs logarithms and exponentials");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (k < 1) throw InputError("hyperexp_fit_longtail: at least one component is required");
    if (b <= one) throw InputError("hyperexp_fit_longtail: the spacing b must exceed 1");
    if (decade <= b)
        throw InputError("hyperexp_fit_longtail: the decade ratio must exceed the spacing b, or "
                         "the fitting arguments would interleave");

    std::vector<T> cs(k);
    cs[0] = c1;
    for (std::size_t i = 1; i < k; ++i) cs[i] = cs[i - 1] / decade;

    std::vector<T> p(k, zero), lam(k, zero);
    for (std::size_t i = 0; i < k; ++i) {
        const T ci = cs[i];
        // Eqs. (4.6)-(4.7): what the already-fitted, slower components leave.
        T residC = ccdf(ci), residBC = ccdf(T(b * ci));
        for (std::size_t j = 0; j < i; ++j) {
            residC -= p[j] * exp(-lam[j] * ci);
            residBC -= p[j] * exp(-lam[j] * b * ci);
        }
        if (i + 1 < k) {
            if (residC <= zero || residBC <= zero || residC <= residBC)
                throw InputError("hyperexp_fit_longtail: the residual ccdf is not positive and "
                                 "decreasing at a fitting argument. The recursion needs the "
                                 "arguments well separated, c_i/c_(i+1) >> b; widen decade, lower "
                                 "k, or move c1 further into the tail");
            lam[i] = log(residC / residBC) / ((b - one) * ci);      // eq. (4.10)
            p[i] = residC * exp(lam[i] * ci);                       // eq. (4.11)
        } else {
            // Eqs. (4.12)-(4.14): the last component takes the rest of the mass.
            T rest = one;
            for (std::size_t j = 0; j < i; ++j) rest -= p[j];
            if (rest <= zero)
                throw InputError("hyperexp_fit_longtail: the fitted components already carry all "
                                 "the probability, so the last one has none left");
            if (residC <= zero)
                throw InputError("hyperexp_fit_longtail: the residual ccdf has gone non-positive "
                                 "at the last fitting argument");
            p[i] = rest;
            lam[i] = log(p[i] / residC) / ci;                       // eq. (4.14)
        }
        if (lam[i] <= zero)
            throw InputError("hyperexp_fit_longtail: a non-positive rate came out of the fit; the "
                             "ccdf is not decaying fast enough for this many components");
    }

    auto fitted = [&](const T& t) {
        T v = zero;
        for (std::size_t j = 0; j < k; ++j) v += p[j] * exp(-lam[j] * t);
        return v;
    };
    HyperexpLongtailResult<T> r;
    r.p = p;
    r.lambda = lam;
    r.points = cs;
    r.mean = zero;
    for (std::size_t j = 0; j < k; ++j) r.mean += p[j] / lam[j];
    r.coverageLow = cs[k - 1];
    r.coverageHigh = cs[0] * b;
    // The target mean over the covered range, by the trapezoid rule: a k too
    // small to reach the body shows up here and nowhere else.
    const std::size_t gn = 20000;
    const T h = r.coverageHigh / num_traits<T>::from_int(static_cast<long>(gn));
    T acc = (ccdf(zero) + ccdf(r.coverageHigh)) / num_traits<T>::from_int(2);
    for (std::size_t i = 1; i < gn; ++i)
        acc += ccdf(T(num_traits<T>::from_int(static_cast<long>(i)) * h));
    r.targetMean = acc * h;
    r.maxRelError = zero;
    for (std::size_t i = 0; i < k; ++i) {
        for (int j = 0; j < 2; ++j) {
            const T t = j == 0 ? cs[i] : T(b * cs[i]);
            const T target = ccdf(t);
            if (target > zero) {
                const T e = num_abs(T(fitted(t) - target)) / target;
                if (e > r.maxRelError) r.maxRelError = e;
            }
        }
    }
    // The fit is exact at the fitting arguments by construction; this says
    // whether it also holds BETWEEN them.
    r.maxRelErrorGrid = zero;
    const T loLog = log(r.coverageLow), hiLog = log(r.coverageHigh);
    for (int i = 0; i < 200; ++i) {
        const T t = exp(loLog + (hiLog - loLog) * num_traits<T>::from_rational(i, 199));
        const T target = ccdf(t);
        if (target > num_traits<T>::from_double(1e-300)) {
            const T e = num_abs(T(fitted(t) - target)) / target;
            if (e > r.maxRelErrorGrid) r.maxRelErrorGrid = e;
        }
    }
    return r;
}

/**
 * The fit with the component count chosen automatically: one per decade between
 * the 0.9 quantile and the 1e-6 quantile, retrying with fewer when the
 * recursion runs out of probability near the body.
 *
 * @param ccdf   F^c(t) = P(X > t)
 * @param b      the within-scale spacing
 * @param decade the ratio between successive fitting arguments
 */
template <class T, class Ccdf>
HyperexpLongtailResult<T> hyperexp_fit_longtail(Ccdf&& ccdf,
                                                const T& b = num_traits<T>::from_rational(3, 2),
                                                const T& decade = num_traits<T>::from_int(4)) {
    using std::log;
    const T top = detail::hefit_quantile<T>(ccdf, num_traits<T>::from_double(1e-6));
    const T body = detail::hefit_quantile<T>(ccdf, num_traits<T>::from_rational(9, 10));
    if (body <= num_traits<T>::from_int(0) || top <= body)
        throw InputError("hyperexp_fit_longtail: the ccdf gives no usable range of time scales");
    long k0 = static_cast<long>(std::llround(num_traits<T>::to_double(T(log(top / body) / log(decade))))) + 1;
    if (k0 < 2) k0 = 2;
    // The recursion needs each component to dominate at its own scale. Near the
    // body of a law with a lot of mass there (a Pareto, say) that fails and the
    // remaining probability runs out; back off one component at a time.
    for (long k = k0; k >= 2; --k) {
        try {
            return hyperexp_fit_longtail_k<T>(ccdf, static_cast<std::size_t>(k), top, b, decade);
        } catch (const InputError&) {
            continue;
        }
    }
    throw InputError("hyperexp_fit_longtail: no component count admits the recursion; the ccdf may "
                     "not be long-tailed enough for this scheme");
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_HYPEREXP_FIT_LONGTAIL_H
