/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_M3PP2M_FITC_H
#define LINE_API_MAM_M3PP2M_FITC_H

/**
 * Marked MMPP(2) with m classes, fitted to counting-process characteristics
 * (matlab/lib/m3a/m3a/m3pp/m3pp2m_fitc.m).
 *
 * The underlying MMPP(2) comes from mmpp2_fitc; the per-class split is then
 * closed form. For each of the first m-1 classes the pair (q1i, q2i) of
 * per-phase marking probabilities follows from the class rate a_i and dvt3(i),
 * the difference between the variance of class i and that of all other classes
 * combined at resolution t3. The last class absorbs the remainder,
 * Dm = diag(1 - sum_i q1i, 1 - sum_i q2i) .* D1.
 *
 * The reference expressions contain sinh(u) exp(-u) with u = (r1 + r2) t / 2
 * nine times each; that product is (1 - exp(-(r1 + r2) t))/2 exactly, and is
 * evaluated in that form here -- algebraically identical, and free of the
 * cancellation that sinh times a decaying exponential suffers at large t.
 *
 * Not ported: m3pp2m_fitc_approx and m3pp2m_fitc_approx_ag, which wrap this
 * split in a quadratic program (quadprog) over the per-class variances, and
 * mmpp2_fitc_approx, whose MMPP(2) stage is an optimproblem/solve call rather
 * than a closed form.
 *
 * Gated on transcendental arithmetic: through mmpp2_fitc, and through the
 * exponential in SH.
 */

#include <vector>

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmpp2_fitc.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of m3pp2m_fitc. */
template <class T>
struct M3pp2mFitcResult {
    Mmap<T> mmap;
    bool degenerate;  ///< the underlying MMPP(2) degenerated to a Poisson process
};

/**
 * Fit an M3PP(2, m). ai holds the per-class arrival rates (summing to a),
 * dvt3 the per-class variance differences at resolution t3.
 */
template <class T>
M3pp2mFitcResult<T> m3pp2m_fitc(const T& a, const T& bt1, const T& bt2, const T& binf,
                                const T& m3t2, const T& t1, const T& t2, const std::vector<T>& ai,
                                const std::vector<T>& dvt3, const T& t3) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp2m_fitc requires transcendental arithmetic");
    using fitdetail::num_exp;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    const std::size_t m = ai.size();
    if (m == 0) throw InputError("m3pp2m_fitc: no classes");
    if (dvt3.size() + 1 < m) throw InputError("m3pp2m_fitc: dvt3 shorter than the class count");
    T asum = zero;
    for (std::size_t i = 0; i < m; ++i) asum += ai[i];
    if (num_abs(T(a - asum)) > num_traits<T>::from_double(1e-8))
        throw InputError("m3pp2m_fitc: inconsistent per-class arrival rates");

    const Mmpp2FitcResult<T> base = mmpp2_fitc(a, bt1, bt2, binf, m3t2, t1, t2);

    M3pp2mFitcResult<T> res;
    res.degenerate = base.degenerate;
    res.mmap.D0 = base.map.D0;
    res.mmap.D1 = base.map.D1;

    if (base.degenerate) {
        // Marked Poisson process: the per-class matrices are the class rates.
        for (std::size_t i = 0; i < m; ++i) {
            Matrix<T> Dc(1, 1, ai[i]);
            res.mmap.Dc.push_back(Dc);
        }
        return res;
    }

    const T l1 = base.map.D1(0, 0);
    const T l2 = base.map.D1(1, 1);
    const T r1 = base.map.D0(0, 1);
    const T r2 = base.map.D0(1, 0);
    const T t = t3;
    // sinh(u) exp(-u) with u = (r1 + r2) t / 2
    const T SH = (one - num_exp(T(-(r1 + r2) * t))) / two;

    std::vector<T> q1(m, zero), q2(m, zero);
    for (std::size_t i = 0; i + 1 < m; ++i) {
        const T a_1 = ai[i];
        const T dv_1 = dvt3[i];
        q1[i] = -(dv_1*pw(r1,4) + dv_1*pw(r2,4) - 2*a_1*pw(r1,4)*t - 2*a_1*pw(r2,4)*t + 4*dv_1*r1*pw(r2,3) +
                         4*dv_1*pw(r1,3)*r2 + l1*pw(r2,4)*t + l2*pw(r1,4)*t + 6*dv_1*pw(r1,2)*pw(r2,2) + 4*a_1*l1*pw(r2,3)*
                        t - 4*a_1*l2*pw(r2,3)*t - 8*a_1*r1*pw(r2,3)*t - 8*a_1*pw(r1,3)*r2*t + 3*l1*r1*pw(r2,3)*t + l1*
                        pw(r1,3)*r2*t + l2*r1*pw(r2,3)*t + 3*l2*pw(r1,3)*r2*t - 12*a_1*pw(r1,2)*pw(r2,2)*t + 3*l1*pw(r1,2)*
                        pw(r2,2)*t + 2*pw(l1,2)*r1*pw(r2,2)*t + 2*pw(l1,2)*pw(r1,2)*r2*t + 3*l2*pw(r1,2)*pw(r2,2)*t +
                         2*pw(l2,2)*r1*pw(r2,2)*t + 2*pw(l2,2)*pw(r1,2)*r2*t - 8*a_1*l1*pw(r2,2)*SH + 8*a_1*l2*pw(r2,2)*
                        SH - 4*pw(l1,2)*r1*r2*SH - 4*pw(l2,2)*r1*r2*SH + 8*a_1*l1*r1*pw(r2,2)*t + 4*a_1*l1*pw(r1,2)*
                        r2*t - 8*a_1*l2*r1*pw(r2,2)*t - 4*a_1*l2*pw(r1,2)*r2*t - 4*l1*l2*r1*pw(r2,2)*t - 4*l1*l2*pw(r1,2)*
                        r2*t - 8*a_1*l1*r1*r2*SH + 8*a_1*l2*r1*r2*SH + 8*l1*l2*r1*r2*SH)/(4*l1*r2*(r1 + r2)*(2*l1*SH -
                         2*l2*SH - l1*r1*t - l1*r2*t + l2*r1*t + l2*r2*t));
        q2[i] = (dv_1*pw(r1,4) + dv_1*pw(r2,4) - 2*a_1*pw(r1,4)*t - 2*a_1*pw(r2,4)*t + 4*dv_1*r1*pw(r2,3) +
                         4*dv_1*pw(r1,3)*r2 + l1*pw(r2,4)*t + l2*pw(r1,4)*t + 6*dv_1*pw(r1,2)*pw(r2,2) - 4*a_1*l1*pw(r1,3)*
                        t + 4*a_1*l2*pw(r1,3)*t - 8*a_1*r1*pw(r2,3)*t - 8*a_1*pw(r1,3)*r2*t + 3*l1*r1*pw(r2,3)*t + l1*
                        pw(r1,3)*r2*t + l2*r1*pw(r2,3)*t + 3*l2*pw(r1,3)*r2*t - 12*a_1*pw(r1,2)*pw(r2,2)*t + 3*l1*pw(r1,2)*
                        pw(r2,2)*t + 2*pw(l1,2)*r1*pw(r2,2)*t + 2*pw(l1,2)*pw(r1,2)*r2*t + 3*l2*pw(r1,2)*pw(r2,2)*t +
                         2*pw(l2,2)*r1*pw(r2,2)*t + 2*pw(l2,2)*pw(r1,2)*r2*t + 8*a_1*l1*pw(r1,2)*SH - 8*a_1*l2*pw(r1,2)*
                        SH - 4*pw(l1,2)*r1*r2*SH - 4*pw(l2,2)*r1*r2*SH - 4*a_1*l1*r1*pw(r2,2)*t - 8*a_1*l1*pw(r1,2)*
                        r2*t + 4*a_1*l2*r1*pw(r2,2)*t + 8*a_1*l2*pw(r1,2)*r2*t - 4*l1*l2*r1*pw(r2,2)*t - 4*l1*l2*pw(r1,2)*
                        r2*t + 8*a_1*l1*r1*r2*SH - 8*a_1*l2*r1*r2*SH + 8*l1*l2*r1*r2*SH)/(4*(r1 + r2)*(pw(l2,2)*pw(r1,2)*
                        t - 2*pw(l2,2)*r1*SH - l1*l2*pw(r1,2)*t + pw(l2,2)*r1*r2*t + 2*l1*l2*r1*SH - l1*l2*r1*r2*t));
    }
    T s1 = zero, s2 = zero;
    for (std::size_t i = 0; i + 1 < m; ++i) {
        s1 += q1[i];
        s2 += q2[i];
    }
    q1[m - 1] = one - s1;
    q2[m - 1] = one - s2;

    for (std::size_t i = 0; i < m; ++i) {
        Matrix<T> Dc(2, 2, zero);
        Dc(0, 0) = q1[i] * base.map.D1(0, 0);
        Dc(1, 1) = q2[i] * base.map.D1(1, 1);
        res.mmap.Dc.push_back(Dc);
    }
    return res;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_M3PP2M_FITC_H
