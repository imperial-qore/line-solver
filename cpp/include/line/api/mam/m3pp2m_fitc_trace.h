/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_M3PP2M_FITC_TRACE_H
#define LINE_API_MAM_M3PP2M_FITC_TRACE_H

/**
 * M3PP(2, m) fitted to the counting process of a multi-class TRACE
 * (matlab/lib/m3a/m3a/m3pp/m3pp2m_fitc_trace.m), and the analogous
 * MMPP(2) entry point for a given MAP
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fitc_theoretical.m).
 *
 * The trace is reduced to counting windows with mtrace_iat2counts at two
 * resolutions, and the aggregate characteristics (rate, IDC at t1 and at tinf,
 * third central moment at t2) plus the per-class ones are read off those
 * windows; the four methods then differ only in WHICH per-class statistic they
 * feed to the corresponding fitter:
 *   - exact_delta  -> m3pp2m_fitc,                  per-class variance DIFFERENCE
 *   - approx_delta -> m3pp2m_fitc_approx,           same, least squares
 *   - approx_cov   -> m3pp22_fitc_approx_cov,       pairwise covariance, 2 classes only
 *   - approx_ag    -> m3pp2m_fitc_approx_ag,        variance plus covariance
 *
 * Two conventions of the reference are load-bearing and are kept verbatim.
 * First, bt2 is set EQUAL to bt1 rather than measured at t2: the reference
 * counts at t1 and reuses the window (mNt2 = mNt1), so the second time scale
 * enters only through t2 = t1 + mean(T) in the fitter's own algebra. Second,
 * t3 = tinf and not t1; the reference carries a comment recording that change,
 * and it is what sets how far out the per-class characteristics are matched.
 *
 * MATLAB's `var` normalizes by n - 1 and its `cov` likewise; both are used here
 * with the same denominator, since the fitted characteristics are compared
 * against the reference's.
 *
 * Gated on transcendental arithmetic, through the fitters.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/m3pp22_fitc_cov.h"
#include "line/api/mam/m3pp2m_fitc.h"
#include "line/api/mam/m3pp2m_fitc_approx.h"
#include "line/api/mam/m3pp_superpos_fitc.h"
#include "line/api/mam/map_count_mean.h"
#include "line/api/mam/map_count_moment.h"
#include "line/api/mam/map_count_var.h"
#include "line/api/trace/mtrace_iat2counts.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of m3pp2m_fitc_trace. */
template <class T>
struct M3pp2mFitcTraceResult {
    Mmap<T> mmap;        ///< the fitted M3PP(2, m)
    std::string method;  ///< the method actually run
    T a;                 ///< the aggregate rate read off the trace
    T bt1;               ///< IDC at t1
    T binf;              ///< IDC at tinf
    T m3t2;              ///< third central moment of counts at t2
};

namespace fitdetail {

/** Sample variance, MATLAB `var`, of the row sums of a count matrix. */
template <class T>
T rowsum_var(const Matrix<long>& C) {
    const std::size_t n = C.rows();
    if (n < 2) throw InputError("m3pp2m_fitc_trace: fewer than two counting windows");
    std::vector<T> s(n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        long acc = 0;
        for (std::size_t j = 0; j < C.cols(); ++j) acc += C(i, j);
        s[i] = num_traits<T>::from_int(acc);
    }
    T mu = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) mu += s[i];
    mu /= num_traits<T>::from_int(static_cast<long>(n));
    T acc = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) acc += (s[i] - mu) * (s[i] - mu);
    return acc / num_traits<T>::from_int(static_cast<long>(n - 1));
}

/** Sample raw moment of order `power` of the row sums of a count matrix. */
template <class T>
T rowsum_moment(const Matrix<long>& C, unsigned power) {
    const std::size_t n = C.rows();
    T acc = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        long s = 0;
        for (std::size_t j = 0; j < C.cols(); ++j) s += C(i, j);
        T p = num_traits<T>::from_int(1);
        for (unsigned k = 0; k < power; ++k) p *= num_traits<T>::from_int(s);
        acc += p;
    }
    return acc / num_traits<T>::from_int(static_cast<long>(n));
}

/** Sample variance of "column c" and of "everything but column c". */
template <class T>
void split_var(const Matrix<long>& C, std::size_t c, T& vc, T& vrest) {
    const std::size_t n = C.rows();
    if (n < 2) throw InputError("m3pp2m_fitc_trace: fewer than two counting windows");
    std::vector<T> x(n), y(n);
    for (std::size_t i = 0; i < n; ++i) {
        long tot = 0;
        for (std::size_t j = 0; j < C.cols(); ++j) tot += C(i, j);
        x[i] = num_traits<T>::from_int(C(i, c));
        y[i] = num_traits<T>::from_int(tot - C(i, c));
    }
    T mx = num_traits<T>::from_int(0), my = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        mx += x[i];
        my += y[i];
    }
    const T nn = num_traits<T>::from_int(static_cast<long>(n));
    mx /= nn;
    my /= nn;
    T ax = num_traits<T>::from_int(0), ay = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        ax += (x[i] - mx) * (x[i] - mx);
        ay += (y[i] - my) * (y[i] - my);
    }
    const T nm1 = num_traits<T>::from_int(static_cast<long>(n - 1));
    vc = ax / nm1;
    vrest = ay / nm1;
}

/** Sample covariance between "column c" and "everything but column c". */
template <class T>
T split_cov(const Matrix<long>& C, std::size_t c) {
    const std::size_t n = C.rows();
    if (n < 2) throw InputError("m3pp2m_fitc_trace: fewer than two counting windows");
    std::vector<T> x(n), y(n);
    for (std::size_t i = 0; i < n; ++i) {
        long tot = 0;
        for (std::size_t j = 0; j < C.cols(); ++j) tot += C(i, j);
        x[i] = num_traits<T>::from_int(C(i, c));
        y[i] = num_traits<T>::from_int(tot - C(i, c));
    }
    T mx = num_traits<T>::from_int(0), my = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        mx += x[i];
        my += y[i];
    }
    const T nn = num_traits<T>::from_int(static_cast<long>(n));
    mx /= nn;
    my /= nn;
    T acc = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) acc += (x[i] - mx) * (y[i] - my);
    return acc / num_traits<T>::from_int(static_cast<long>(n - 1));
}

}  // namespace fitdetail

/**
 * Fit a multi-class trace with an M3PP(2, m) on its counting process.
 *
 * @param Tv     inter-arrival times
 * @param A      class labels
 * @param method 'exact_delta', 'approx_delta', 'approx_cov' or 'approx_ag'
 * @param t1     finite time scale
 * @param tinf   near-infinite time scale
 */
template <class T>
M3pp2mFitcTraceResult<T> m3pp2m_fitc_trace(const std::vector<T>& Tv, const std::vector<int>& A,
                                           const std::string& method, const T& t1, const T& tinf) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp2m_fitc_trace requires transcendental arithmetic");
    if (Tv.empty()) throw InputError("m3pp2m_fitc_trace: empty trace");
    if (A.size() != Tv.size())
        throw InputError("m3pp2m_fitc_trace: labels and inter-arrival times disagree");

    const trace::MtraceCountsResult<T> N1 = trace::mtrace_iat2counts(Tv, A, t1);
    const trace::MtraceCountsResult<T> Ninf = trace::mtrace_iat2counts(Tv, A, tinf);
    const std::size_t m = N1.labels.size();
    if (method == "approx_cov" && m > 2)
        throw InputError("m3pp2m_fitc_trace: approximate covariance fitting only supports two "
                         "classes");

    T sumT = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < Tv.size(); ++i) sumT += Tv[i];
    const T mean = sumT / num_traits<T>::from_int(static_cast<long>(Tv.size()));
    const T a = num_traits<T>::from_int(1) / mean;
    const T t2 = t1 + mean;
    const T t3 = tinf;  // the reference sets the third scale to tinf, not t1

    std::vector<T> ai(m);
    for (std::size_t i = 0; i < m; ++i) {
        long cnt = 0;
        for (std::size_t k = 0; k < A.size(); ++k)
            if (A[k] == N1.labels[i]) ++cnt;
        ai[i] = a * num_traits<T>::from_int(cnt) /
                num_traits<T>::from_int(static_cast<long>(A.size()));
    }

    const T bt1 = fitdetail::rowsum_var<T>(N1.counts) / (a * t1);
    const T bt2 = bt1;  // the reference reuses the t1 window for the second scale
    const T binf = fitdetail::rowsum_var<T>(Ninf.counts) / (a * tinf);
    const T m3t2 = fitdetail::m3_from_raw(fitdetail::rowsum_moment<T>(N1.counts, 1),
                                          fitdetail::rowsum_moment<T>(N1.counts, 2),
                                          fitdetail::rowsum_moment<T>(N1.counts, 3));

    M3pp2mFitcTraceResult<T> out;
    out.method = method;
    out.a = a;
    out.bt1 = bt1;
    out.binf = binf;
    out.m3t2 = m3t2;

    if (method == "exact_delta" || method == "approx_delta") {
        std::vector<T> dvt3(m);
        for (std::size_t i = 0; i < m; ++i) {
            T vc = num_traits<T>::from_int(0), vr = num_traits<T>::from_int(0);
            fitdetail::split_var<T>(N1.counts, i, vc, vr);
            dvt3[i] = vc - vr;
        }
        if (method == "exact_delta")
            out.mmap = m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3).mmap;
        else
            out.mmap = m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3).mmap;
        return out;
    }
    if (method == "approx_cov") {
        const T V = fitdetail::rowsum_var<T>(N1.counts);
        T v0 = num_traits<T>::from_int(0), v1 = num_traits<T>::from_int(0);
        T dummy = num_traits<T>::from_int(0);
        fitdetail::split_var<T>(N1.counts, 0, v0, dummy);
        fitdetail::split_var<T>(N1.counts, 1, v1, dummy);
        const T s = (V - v0 - v1) / num_traits<T>::from_int(2);
        out.mmap = m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, s, t3).mmap;
        return out;
    }
    if (method == "approx_ag") {
        std::vector<T> gt3(m);
        for (std::size_t i = 0; i < m; ++i) {
            T vc = num_traits<T>::from_int(0), vr = num_traits<T>::from_int(0);
            fitdetail::split_var<T>(N1.counts, i, vc, vr);
            gt3[i] = vc + fitdetail::split_cov<T>(N1.counts, i);
        }
        out.mmap = m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3).mmap;
        return out;
    }
    throw InputError("m3pp2m_fitc_trace: unknown method '" + method + "'");
}

/** m3pp2m_fitc_trace with the reference's default time scales. */
template <class T>
M3pp2mFitcTraceResult<T> m3pp2m_fitc_trace(const std::vector<T>& Tv, const std::vector<int>& A,
                                           const std::string& method) {
    if (Tv.empty()) throw InputError("m3pp2m_fitc_trace: empty trace");
    T sumT = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < Tv.size(); ++i) sumT += Tv[i];
    const T mean = sumT / num_traits<T>::from_int(static_cast<long>(Tv.size()));
    const T t1 = num_traits<T>::from_int(10) * mean;
    const T span = (sumT - Tv[0]) / num_traits<T>::from_int(100);
    const T ten_t = num_traits<T>::from_int(10) * t1;
    return m3pp2m_fitc_trace(Tv, A, method, t1, ten_t > span ? ten_t : span);
}

/**
 * MMPP(2) fitted to the counting characteristics of a GIVEN MAP
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fitc_theoretical.m). The characteristics
 * are the rate at t1, the IDC at t1, t2 and tinf, and the third central moment
 * of counts at t2, all evaluated exactly on the input MAP.
 */
template <class T>
Mmpp2FitcResult<T> mmpp2_fitc_theoretical(const Map<T>& mp, const T& t1, const T& t2,
                                          const T& tinf) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmpp2_fitc_theoretical requires transcendental arithmetic");
    std::vector<T> ts;
    ts.push_back(t1);
    ts.push_back(t2);
    ts.push_back(tinf);
    const std::vector<T> mu = map_count_mean(mp, ts);
    const std::vector<T> vr = map_count_var(mp, ts);

    const T a = mu[0] / t1;
    const T bt1 = vr[0] / (a * t1);
    const T bt2 = vr[1] / (a * t2);
    const T binf = vr[2] / (a * tinf);

    std::vector<unsigned> orders;
    orders.push_back(1);
    orders.push_back(2);
    orders.push_back(3);
    const std::vector<T> mt2 = map_count_moment(mp, t2, orders);
    const T m3t2 = fitdetail::m3_from_raw(mt2[0], mt2[1], mt2[2]);

    return mmpp2_fitc(a, bt1, bt2, binf, m3t2, t1, t2);
}

/** mmpp2_fitc_theoretical with the reference's default time scales 1, 10, 1e8. */
template <class T>
Mmpp2FitcResult<T> mmpp2_fitc_theoretical(const Map<T>& mp) {
    return mmpp2_fitc_theoretical(mp, num_traits<T>::from_int(1), num_traits<T>::from_int(10),
                                  num_traits<T>::from_double(1e8));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_M3PP2M_FITC_TRACE_H
