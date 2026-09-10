/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_M3PP22_FITC_COV_H
#define LINE_API_MAM_M3PP22_FITC_COV_H

/**
 * M3PP(2, 2) fitted to the count COVARIANCE between its two classes
 * (matlab/lib/m3a/m3a/m3pp/m3pp22_fitc_approx_cov_multiclass.m and
 *  matlab/lib/m3a/m3a/m3pp/m3pp22_fitc_approx_cov.m).
 *
 * TWO CLASSES ONLY, by construction, and the reference refuses more by name.
 * Given the underlying MMPP(2) the per-phase marking probabilities (q1, q2) of
 * the first class satisfy two relations: the class rate a1 is affine in them,
 * and the count covariance sigma at t3 is a QUADRATIC in q2 once q1 has been
 * eliminated,
 *
 *   sigma(q2) = w0 + w1 q2 + w2 q2^2 ,
 *
 * so the inverse has two roots. Rather than picking one and repairing the
 * result, the reference derives, for EACH root separately, the interval of
 * covariances over which that root keeps both marking probabilities inside
 * [0, 1] and the discriminant non-negative; it then clamps the requested
 * covariance into whichever interval is closer to it and takes the matching
 * root. A root whose interval is provably empty is flagged rather than clamped.
 * That bound derivation is the bulk of the reference and of this port, and it
 * is transcribed relation by relation: square-root argument >= 0, q2 >= 0,
 * q1 >= 0, q2 <= 1, q1 <= 1.
 *
 * The clamp is the only approximation: the rates are matched exactly, and the
 * covariance is matched exactly whenever the request lies inside the feasible
 * interval of either root.
 *
 * Degenerate inputs short-circuit as in the reference: a Poisson underlying
 * process splits D1 in proportion to the class rates, and a single class takes
 * all of D1.
 *
 * Gated on transcendental arithmetic: the bounds carry exp(-(r1 + r2) t3).
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmpp2_fitc_approx.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of the covariance-matching M3PP(2, 2) fits. */
template <class T>
struct M3pp22FitcCovResult {
    Mmap<T> mmap;      ///< the fitted M3PP(2, 2)
    T sigma;           ///< the covariance actually realised, after clamping
    int root;          ///< 1 or 2, which root of the quadratic was taken; 0 when degenerate
    bool clamped;      ///< the requested covariance lay outside the feasible interval
    bool degenerate;   ///< the underlying process was Poisson, or there is a single class
};

/**
 * Split a GIVEN MMPP(2) into two classes, matching the per-class rates exactly
 * and the count covariance between them at t3 as closely as feasible.
 *
 * @param mmpp the underlying MAP, of order 1 (Poisson) or 2
 * @param ai   the per-class rates; at most two
 * @param st3  the requested count covariance between the two classes at t3
 * @param t3   the third time scale
 */
template <class T>
M3pp22FitcCovResult<T> m3pp22_fitc_approx_cov_multiclass(const Map<T>& mmpp,
                                                         const std::vector<T>& ai, const T& st3,
                                                         const T& t3) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp22_fitc_approx_cov_multiclass requires transcendental arithmetic");
    using fitdetail::num_exp;
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T four = num_traits<T>::from_int(4);

    const std::size_t m = ai.size();
    if (m == 0) throw InputError("m3pp22_fitc_approx_cov_multiclass: no classes");
    if (m > 2) throw InputError("m3pp22_fitc_approx_cov_multiclass: no more than two classes "
                                "supported");

    M3pp22FitcCovResult<T> res;
    res.mmap.D0 = mmpp.D0;
    res.mmap.D1 = mmpp.D1;
    res.sigma = zero;
    res.root = 0;
    res.clamped = false;
    res.degenerate = false;

    if (mmpp.D0.rows() == 1) {
        // marked Poisson process: split D1 in proportion to the class rates
        T asum = zero;
        for (std::size_t i = 0; i < m; ++i) asum += ai[i];
        for (std::size_t i = 0; i < m; ++i)
            res.mmap.Dc.push_back(Matrix<T>(1, 1, T(ai[i] / asum * mmpp.D1(0, 0))));
        res.degenerate = true;
        return res;
    }
    if (m == 1) {
        res.mmap.Dc.push_back(mmpp.D1);
        res.degenerate = true;
        return res;
    }
    if (mmpp.D0.rows() != 2)
        throw InputError("m3pp22_fitc_approx_cov_multiclass: the underlying MAP must have order 2");

    const T l1 = mmpp.D1(0, 0);
    const T l2 = mmpp.D1(1, 1);
    const T r1 = mmpp.D0(0, 1);
    const T r2 = mmpp.D0(1, 0);
    const T t = t3;
    const T a1 = ai[0];

    const T E = num_exp(T(-(r1 + r2) * t));
    const T G = one - E - (r1 + r2) * t;
    const T w0 = (two * r1 * G * (a1 * a1 * (r1 + r2) - a1 * r2 * (l1 - l2))) /
                 (r2 * pw(T(r1 + r2), 3));
    const T w1 = -(two * r1 * G * (two * a1 * l2 * (r1 + r2) - l2 * r2 * (l1 - l2))) /
                 (r2 * pw(T(r1 + r2), 3));
    const T w2 = ((two * l2 * l2 * r2 * t) * (r1 + r2) + two * l2 * l2 * r1 * (one - E)) /
                     (r2 * pw(T(r1 + r2), 2)) -
                 (two * l2 * l2 * t) / r2;
    const T w3 = (r1 + r2) / (l2 * r1);
    const T w4 = (l1 * r2) / (l2 * r1);
    if (w2 == zero)
        throw NumericError("m3pp22_fitc_approx_cov_multiclass: the covariance is linear in the "
                           "marking probability, so the reference's two-root inversion degenerates");

    const T inf = num_traits<T>::from_double(1.0 / 0.0);
    T L1 = -inf, L2 = -inf, U1 = inf, U2 = inf;
    bool infeasible1 = false, infeasible2 = false;

    const T z = w0 - w1 * w1 / (four * w2);

    // square-root argument >= 0
    if (w2 > zero) {
        if (z > L1) L1 = z;
        if (z > L2) L2 = z;
    } else {
        if (z < U1) U1 = z;
        if (z < U2) U2 = z;
    }
    // q2 >= 0
    if (w1 >= zero) {
        if (w0 > L1) L1 = w0;
    } else if (w2 < zero) {
        infeasible1 = true;
    }
    if (w1 <= zero) {
        if (w0 < U2) U2 = w0;
    } else if (w2 > zero) {
        infeasible2 = true;
    }
    // q1 >= 0
    {
        const T tmp = two * a1 * w3 * w2 + w1;
        const T bnd = z + tmp * tmp / (four * w2);
        if (tmp >= zero) {
            if (bnd < U1) U1 = bnd;
        } else if (w2 > zero) {
            infeasible1 = true;
        }
        if (tmp <= zero) {
            if (bnd > L2) L2 = bnd;
        } else if (w2 < zero) {
            infeasible2 = true;
        }
    }
    // q2 <= 1
    {
        const T tmp = two * w2 + w1;
        const T bnd = z + tmp * tmp / (four * w2);
        if (tmp >= zero) {
            if (bnd < U1) U1 = bnd;
        } else if (w2 > zero) {
            infeasible1 = true;
        }
        if (tmp <= zero) {
            if (bnd > L2) L2 = bnd;
        } else if (w2 < zero) {
            infeasible2 = true;
        }
    }
    // q1 <= 1
    {
        const T tmp = two * a1 * w2 * w3 - two * w2 * w4 + w1;
        const T bnd = z + tmp * tmp / (four * w2);
        if (tmp >= zero) {
            if (bnd > L1) L1 = bnd;
        } else if (w2 < zero) {
            infeasible1 = true;
        }
        if (tmp <= zero) {
            if (bnd < U2) U2 = bnd;
        } else if (w2 > zero) {
            infeasible2 = true;
        }
    }

    if (infeasible1 && infeasible2)
        throw NumericError("m3pp22_fitc_approx_cov_multiclass: empty feasibility region");

    T sigma = zero;
    int root = 0;
    if (infeasible2) {
        sigma = st3 < U1 ? st3 : U1;
        if (sigma < L1) sigma = L1;
        root = 1;
    } else if (infeasible1) {
        sigma = st3 < U2 ? st3 : U2;
        if (sigma < L2) sigma = L2;
        root = 2;
    } else {
        T s1 = st3 < U1 ? st3 : U1;
        if (s1 < L1) s1 = L1;
        T s2 = st3 < U2 ? st3 : U2;
        if (s2 < L2) s2 = L2;
        if (num_abs(T(s1 - st3)) < num_abs(T(s2 - st3))) {
            sigma = s1;
            root = 1;
        } else {
            sigma = s2;
            root = 2;
        }
    }

    const T disc = w1 * w1 - four * w2 * (w0 - sigma);
    if (disc < zero)
        throw NumericError("m3pp22_fitc_approx_cov_multiclass: negative discriminant after "
                           "clamping the covariance");
    const T rt = num_sqrt(disc);
    T q2 = root == 1 ? T((-w1 + rt) / (two * w2)) : T((-w1 - rt) / (two * w2));
    T q1 = (a1 * (r1 + r2) - l2 * q2 * r1) / (l1 * r2);

    const T tol = num_traits<T>::from_double(1e-8);
    if (!(q1 >= tol && q1 <= one + tol && q2 >= tol && q2 <= one + tol))
        throw NumericError("m3pp22_fitc_approx_cov_multiclass: the marking probabilities left "
                           "the unit box");
    if (q1 < zero) q1 = zero;
    if (q1 > one) q1 = one;
    if (q2 < zero) q2 = zero;
    if (q2 > one) q2 = one;

    Matrix<T> Dc1(2, 2, zero), Dc2(2, 2, zero);
    Dc1(0, 0) = q1 * l1;
    Dc1(1, 1) = q2 * l2;
    Dc2(0, 0) = (one - q1) * l1;
    Dc2(1, 1) = (one - q2) * l2;
    res.mmap.Dc.push_back(Dc1);
    res.mmap.Dc.push_back(Dc2);
    res.sigma = sigma;
    res.root = root;
    res.clamped = num_abs(T(sigma - st3)) > num_traits<T>::from_double(1e-12);
    return res;
}

/**
 * Fit the underlying MMPP(2) by optimization, then apply the covariance split.
 *
 * @param a,bt1,bt2,binf,m3t2,t1,t2 the aggregate counting characteristics
 * @param ai   the rates of the two classes, which must sum to a
 * @param st3  the requested count covariance between them at t3
 * @param t3   the third time scale
 * @param opt  tuning of the underlying MMPP(2) solve
 */
template <class T>
M3pp22FitcCovResult<T> m3pp22_fitc_approx_cov(const T& a, const T& bt1, const T& bt2, const T& binf,
                                              const T& m3t2, const T& t1, const T& t2,
                                              const std::vector<T>& ai, const T& st3, const T& t3,
                                              const AugLagOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp22_fitc_approx_cov requires transcendental arithmetic");
    T asum = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < ai.size(); ++i) asum += ai[i];
    if (num_abs(T(a - asum)) > num_traits<T>::from_double(1e-8))
        throw InputError("m3pp22_fitc_approx_cov: inconsistent per-class arrival rates");

    const Mmpp2FitcApproxResult<T> base = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, opt);
    return m3pp22_fitc_approx_cov_multiclass(base.map, ai, st3, t3);
}

/** m3pp22_fitc_approx_cov with the default tuning of the MMPP(2) solve. */
template <class T>
M3pp22FitcCovResult<T> m3pp22_fitc_approx_cov(const T& a, const T& bt1, const T& bt2, const T& binf,
                                              const T& m3t2, const T& t1, const T& t2,
                                              const std::vector<T>& ai, const T& st3,
                                              const T& t3) {
    AugLagOptions<T> opt = auglag_defaults<T>();
    opt.ctol = num_traits<T>::from_double(1e-12);
    return m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, st3, t3, opt);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_M3PP22_FITC_COV_H
