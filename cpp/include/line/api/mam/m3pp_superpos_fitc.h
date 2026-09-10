/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_M3PP_SUPERPOS_FITC_H
#define LINE_API_MAM_M3PP_SUPERPOS_FITC_H

/**
 * M3PP obtained by SUPERPOSING one second-order process per class
 * (matlab/lib/m3a/m3a/m3pp/m3pp_superpos_fitc.m and its two entry points
 *  m3pp_superpos_fitc_theoretical.m and m3pp_superpos_fitc_trace.m).
 *
 * Each class gets its own MMPP(2) from mmpp2_fitc on that class's own rate,
 * IDC at t, IDC at infinity and third central moment of counts; the class is
 * then a single-class MMAP {D0, D1, D1}, and the k of them are superposed.
 * The result therefore has order 2^k, not k + 1: the phase process of a
 * superposition is the PRODUCT chain, and the reference's own comment "M3PP[m]
 * of order k+1" describes the lumped interleaving of m3pp2m_interleave, not
 * this. What the superposition buys instead is that every per-class second-order
 * characteristic is matched exactly and independently, since the components do
 * not interact.
 *
 * REFERENCE DEFECT (m3pp_superpos_fitc_theoretical.m): it calls
 * mmap_count_moment, which is defined NOWHERE in the MATLAB tree -- neither in
 * m3a nor in matlab/src/api. The entry point therefore raises "Unrecognized
 * function" before any fitting happens and cannot execute. The function is
 * defined in the JAR (jline.api.mam.Mmap_count_moment) as the counting moments
 * of the per-class MARGINAL MAP, {D0 + sum_{j != k} D1_j, D1_k}, which is the
 * standard definition and is what mmap_count_moment in map_count_moment.h
 * implements; that is what this port calls. MATLAB was NOT edited.
 * m3pp_superpos_fitc_trace.m is unaffected: it reads its moments off the trace.
 *
 * Gated on transcendental arithmetic, through mmpp2_fitc.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_count_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmap_stats.h"
#include "line/api/mam/mmpp2_fitc.h"
#include "line/api/trace/mtrace_iat2counts.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of the superposition fits. */
template <class T>
struct M3ppSuperposResult {
    Mmap<T> mmap;                  ///< the superposed process
    std::vector<Mmap<T>> parts;    ///< the per-class components, before superposition
};

namespace fitdetail {

/** Sample central third moment of counts from the first three raw moments. */
template <class T>
T m3_from_raw(const T& m1, const T& m2, const T& m3) {
    const T three = num_traits<T>::from_int(3);
    const T two = num_traits<T>::from_int(2);
    return m3 - three * m2 * m1 + two * m1 * m1 * m1;
}

/** Sample mean of the first `n` entries of a column. */
template <class T>
T col_mean(const Matrix<long>& C, std::size_t col, unsigned power) {
    T acc = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < C.rows(); ++i) {
        T v = num_traits<T>::from_int(static_cast<long>(C(i, col)));
        T p = num_traits<T>::from_int(1);
        for (unsigned k = 0; k < power; ++k) p *= v;
        acc += p;
    }
    return acc / num_traits<T>::from_int(static_cast<long>(C.rows()));
}

/** Sample variance (MATLAB `var`, i.e. the n-1 denominator) of a column. */
template <class T>
T col_var(const Matrix<long>& C, std::size_t col) {
    const std::size_t n = C.rows();
    if (n < 2) throw InputError("m3pp_superpos_fitc: fewer than two counting windows");
    const T mu = col_mean<T>(C, col, 1);
    T acc = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        const T d = num_traits<T>::from_int(static_cast<long>(C(i, col))) - mu;
        acc += d * d;
    }
    return acc / num_traits<T>::from_int(static_cast<long>(n - 1));
}

}  // namespace fitdetail

/**
 * Fit one second-order M3PP per class from its counting characteristics and
 * superpose them.
 *
 * @param av    per-class rates
 * @param btv   per-class IDC(t)
 * @param binfv per-class IDC(inf)
 * @param m3tv  per-class third central moment of counts at t
 * @param t     finite time scale
 * @param tinf  near-infinite time scale
 */
template <class T>
M3ppSuperposResult<T> m3pp_superpos_fitc(const std::vector<T>& av, const std::vector<T>& btv,
                                         const std::vector<T>& binfv, const std::vector<T>& m3tv,
                                         const T& t, const T& tinf) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp_superpos_fitc requires transcendental arithmetic");
    const std::size_t m = av.size();
    if (m == 0) throw InputError("m3pp_superpos_fitc: no classes");
    if (btv.size() != m || binfv.size() != m || m3tv.size() != m)
        throw InputError("m3pp_superpos_fitc: av, btv, binfv and m3tv must have the same length");

    M3ppSuperposResult<T> res;
    for (std::size_t i = 0; i < m; ++i) {
        const Mmpp2FitcResult<T> f = mmpp2_fitc(av[i], btv[i], btv[i], binfv[i], m3tv[i], t, tinf);
        Mmap<T> comp;
        comp.D0 = f.map.D0;
        comp.D1 = f.map.D1;
        comp.Dc.push_back(f.map.D1);
        res.parts.push_back(comp);
    }

    res.mmap = res.parts[0];
    for (std::size_t i = 1; i < m; ++i) res.mmap = mmap_super(res.mmap, res.parts[i]);
    return res;
}

/**
 * Superpose one M3PP per class to fit the counting characteristics of a given
 * MMAP.
 *
 * @param mm   the process to fit
 * @param t    finite time scale
 * @param tinf near-infinite time scale
 */
template <class T>
M3ppSuperposResult<T> m3pp_superpos_fitc_theoretical(const Mmap<T>& mm, const T& t, const T& tinf) {
    const std::size_t m = mm.classes();
    if (m == 0) throw InputError("m3pp_superpos_fitc_theoretical: the MMAP has no classes");

    const std::vector<T> av = mmap_count_mean(mm, num_traits<T>::from_int(1));
    const std::vector<T> btv = mmap_count_idc(mm, t);
    const std::vector<T> binfv = mmap_count_idc(mm, tinf);

    std::vector<unsigned> orders;
    orders.push_back(1);
    orders.push_back(2);
    orders.push_back(3);
    const Matrix<T> mtv = mmap_count_moment(mm, t, orders);

    std::vector<T> m3tv(m, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < m; ++i)
        m3tv[i] = fitdetail::m3_from_raw(mtv(0, i), mtv(1, i), mtv(2, i));

    return m3pp_superpos_fitc(av, btv, binfv, m3tv, t, tinf);
}

/**
 * Superpose one M3PP per class to fit a multi-class trace.
 *
 * @param Tv   inter-arrival times
 * @param A    class labels
 * @param t    finite time scale
 * @param tinf near-infinite time scale
 */
template <class T>
M3ppSuperposResult<T> m3pp_superpos_fitc_trace(const std::vector<T>& Tv, const std::vector<int>& A,
                                               const T& t, const T& tinf) {
    if (Tv.empty()) throw InputError("m3pp_superpos_fitc_trace: empty trace");
    if (A.size() != Tv.size())
        throw InputError("m3pp_superpos_fitc_trace: labels and inter-arrival times disagree");

    T sumT = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < Tv.size(); ++i) sumT += Tv[i];
    const T a = num_traits<T>::from_int(static_cast<long>(Tv.size())) / sumT;

    const trace::MtraceCountsResult<T> Nt = trace::mtrace_iat2counts(Tv, A, t);
    const trace::MtraceCountsResult<T> Ninf = trace::mtrace_iat2counts(Tv, A, tinf);
    const std::size_t m = Nt.labels.size();

    std::vector<T> av(m), btv(m), binfv(m), m3tv(m);
    for (std::size_t i = 0; i < m; ++i) {
        long cnt = 0;
        for (std::size_t k = 0; k < A.size(); ++k)
            if (A[k] == Nt.labels[i]) ++cnt;
        av[i] = a * num_traits<T>::from_int(cnt) /
                num_traits<T>::from_int(static_cast<long>(A.size()));
        btv[i] = fitdetail::col_var<T>(Nt.counts, i) / (av[i] * t);
        binfv[i] = fitdetail::col_var<T>(Ninf.counts, i) / (av[i] * tinf);
        m3tv[i] = fitdetail::m3_from_raw(fitdetail::col_mean<T>(Nt.counts, i, 1),
                                         fitdetail::col_mean<T>(Nt.counts, i, 2),
                                         fitdetail::col_mean<T>(Nt.counts, i, 3));
    }
    return m3pp_superpos_fitc(av, btv, binfv, m3tv, t, tinf);
}

/** m3pp_superpos_fitc_trace with the reference's default time scales. */
template <class T>
M3ppSuperposResult<T> m3pp_superpos_fitc_trace(const std::vector<T>& Tv,
                                               const std::vector<int>& A) {
    if (Tv.empty()) throw InputError("m3pp_superpos_fitc_trace: empty trace");
    T sumT = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < Tv.size(); ++i) sumT += Tv[i];
    const T mean = sumT / num_traits<T>::from_int(static_cast<long>(Tv.size()));
    const T t = num_traits<T>::from_int(10) * mean;
    const T span = (sumT - Tv[0]) / num_traits<T>::from_int(100);
    const T ten_t = num_traits<T>::from_int(10) * t;
    return m3pp_superpos_fitc_trace(Tv, A, t, ten_t > span ? ten_t : span);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_M3PP_SUPERPOS_FITC_H
