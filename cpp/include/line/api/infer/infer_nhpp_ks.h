/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_NHPP_KS_H
#define LINE_API_INFER_NHPP_KS_H

/**
 * Kolmogorov-Smirnov tests for a non-homogeneous Poisson arrival process.
 *
 * Templated port of matlab/src/api/infer/infer_nhpp_ks.m, cross-checked against
 * jar/src/main/java/jline/api/infer/InferNhppKs.java.
 *
 * THE CONDITIONAL-UNIFORM TRANSFORMATION. Conditional on the number of arrivals
 * in [T0,T], the arrival times of an NHPP are the order statistics of iid
 * variables with cdf Lambda(t)/Lambda(T). Mapping the data through that cdf
 * turns ANY NHPP, whatever its rate, into iid uniforms, so one KS test covers
 * every rate function.
 *
 * WHY THE PLAIN TEST IS WEAK, AND WHAT FIXES IT. The CU KS test has "remarkably
 * little power" against non-exponential interarrival times: it looks at the
 * POSITIONS of the points, and those stay nearly uniform for many non-Poisson
 * processes. Lewis (1965) applies the Durbin (1961) transformation first --
 * reorder the GAPS ascending, rescale each by how many gaps remain, cumulate --
 * which turns a difference in the gap DISTRIBUTION into a difference in
 * position. Measured on 400 replications of an Erlang-4 renewal process, the CU
 * test rejects at its own size while the Lewis test rejects essentially always.
 *
 * ARITHMETIC. exp and sqrt in the p-value: transcendental only.
 *
 * Reference: S.-H. Kim, W. Whitt (2014). Are call center and hospital arrivals
 * well modeled by nonhomogeneous Poisson processes? M&SOM 16(3), 464-480;
 * J. Durbin (1961), Biometrika 48, 41-55; P. A. W. Lewis (1965), JRSS B 27.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace infer {

/** Which test to run on the conditional-uniform data. */
enum class NhppKsMethod { Cu, Lewis };

/** Outcome of the KS test. */
template <class T>
struct NhppKsResult {
    T statistic;                ///< the KS distance
    T pvalue;                   ///< asymptotic p-value
    std::size_t n = 0;          ///< number of arrivals used
    std::vector<T> uniforms;    ///< after the conditional-uniform transformation
    std::vector<T> transformed; ///< after the Durbin step, for the Lewis test
};

namespace detail {

/** Two-sided KS distance between the sample and the uniform cdf. */
template <class T>
T nhpp_ks_stat(std::vector<T> u) {
    std::sort(u.begin(), u.end());
    const std::size_t n = u.size();
    const T nT = num_traits<T>::from_int(static_cast<long>(n));
    T d = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        const T up = num_traits<T>::from_int(static_cast<long>(i + 1)) / nT - u[i];
        const T lo = u[i] - num_traits<T>::from_int(static_cast<long>(i)) / nT;
        if (up > d) d = up;
        if (lo > d) d = lo;
    }
    return d;
}

/**
 * Asymptotic Kolmogorov p-value with the small-sample correction of Stephens:
 * the effective argument is (sqrt(n)+0.12+0.11/sqrt(n))D, accurate from n = 5.
 */
template <class T>
T nhpp_ks_pvalue(const T& d, std::size_t n) {
    using std::exp;
    using std::sqrt;
    const T one = num_traits<T>::from_int(1);
    if (n == 0) return one;
    const T nT = num_traits<T>::from_int(static_cast<long>(n));
    const T x = (sqrt(nT) + num_traits<T>::from_double(0.12) +
                 num_traits<T>::from_double(0.11) / sqrt(nT)) * d;
    if (x <= num_traits<T>::from_int(0)) return one;
    T q = num_traits<T>::from_int(0);
    for (int k = 1; k <= 100; ++k) {
        const T term = exp(-num_traits<T>::from_int(2 * k * k) * x * x);
        q += (k % 2 == 1) ? term : T(-term);
    }
    T p = num_traits<T>::from_int(2) * q;
    if (p < num_traits<T>::from_int(0)) p = num_traits<T>::from_int(0);
    if (p > one) p = one;
    return p;
}

}  // namespace detail

/**
 * @param times   the arrival times, within [T0,T]
 * @param T       the right end of the observation interval
 * @param cumRate the cumulative rate Lambda(t); a constant rate when empty
 * @param method  Cu for the plain test, Lewis for the Durbin-transformed one
 * @param T0      the left end of the observation interval
 */
template <class Tv>
NhppKsResult<Tv> infer_nhpp_ks(const std::vector<Tv>& times, const Tv& T,
                               const std::function<Tv(const Tv&)>& cumRate =
                                   std::function<Tv(const Tv&)>(),
                               NhppKsMethod method = NhppKsMethod::Lewis,
                               const Tv& T0 = num_traits<Tv>::from_int(0)) {
    static_assert(num_traits<Tv>::has_transcendental, "infer_nhpp_ks needs exp for the p-value");
    const Tv zero = num_traits<Tv>::from_int(0);
    const Tv one = num_traits<Tv>::from_int(1);
    std::vector<Tv> t;
    for (const Tv& x : times)
        if (x >= T0 && x <= T) t.push_back(x);
    std::sort(t.begin(), t.end());
    const std::size_t n = t.size();
    if (n < 2) throw InputError("infer_nhpp_ks: at least two arrivals are needed to test");

    std::vector<Tv> u(n);
    if (!cumRate) {
        for (std::size_t i = 0; i < n; ++i) u[i] = (t[i] - T0) / (T - T0);
    } else {
        const Tv lo = cumRate(T0), hi = cumRate(T);
        if (hi <= lo)
            throw InputError("infer_nhpp_ks: the cumulative rate must increase over the interval");
        for (std::size_t i = 0; i < n; ++i) u[i] = (cumRate(t[i]) - lo) / (hi - lo);
    }
    for (Tv& v : u) {
        if (v < zero) v = zero;
        if (v > one) v = one;
    }

    NhppKsResult<Tv> r;
    r.n = n;
    r.uniforms = u;
    if (method == NhppKsMethod::Cu) {
        r.transformed = u;
    } else {
        // The Durbin (1961) transformation: gaps, sorted ascending, each
        // rescaled by how many gaps remain, then cumulated.
        std::vector<Tv> v = u;
        std::sort(v.begin(), v.end());
        std::vector<Tv> gaps(n + 1);
        gaps[0] = v[0];
        for (std::size_t i = 1; i < n; ++i) gaps[i] = v[i] - v[i - 1];
        gaps[n] = one - v[n - 1];
        std::sort(gaps.begin(), gaps.end());
        std::vector<Tv> c(n + 1);
        Tv prev = zero;
        for (std::size_t i = 0; i <= n; ++i) {
            c[i] = num_traits<Tv>::from_int(static_cast<long>(n + 1 - i)) * (gaps[i] - prev);
            prev = gaps[i];
        }
        r.transformed.resize(n);
        Tv acc = zero;
        for (std::size_t i = 0; i < n; ++i) {
            acc += c[i];
            Tv s = acc;
            if (s < zero) s = zero;
            if (s > one) s = one;
            r.transformed[i] = s;
        }
    }
    r.statistic = detail::nhpp_ks_stat(r.transformed);
    r.pvalue = detail::nhpp_ks_pvalue(r.statistic, n);
    return r;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_NHPP_KS_H
