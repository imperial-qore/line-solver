/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMPP2_FITC_H
#define LINE_API_MAM_MMPP2_FITC_H

/**
 * MMPP(2) matching counting-process characteristics
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fitc.m).
 *
 * Implements Heffes and Lucantoni (1986): the arrival rate mu, the index of
 * dispersion for counts at two finite time scales and in the limit
 * (bt1, bt2, binf) and the third central moment of the counts at t2 determine
 * the four MMPP(2) parameters (r1, r2, l1, l2) in closed form. The only
 * implicit step is the total switching rate x = r1 + r2, which solves
 *   (binf - 1)/(binf - bt1) = c,   x t1 = W(-c e^-c) + c
 * with W the principal branch of the Lambert function. MATLAB solves
 * w e^w = -c e^-c with fsolve started at w = 1, i.e. it targets the principal
 * branch; this port evaluates W0 directly by Halley iteration, so the result
 * is deterministic and does not depend on a solver's tolerances.
 *
 * Gated on transcendental arithmetic: exp, log and square roots appear
 * throughout, and the Lambert step is inherently tolerance-driven.
 *
 * Degenerate inputs return a Poisson process of rate mu, as in the reference:
 * a constant unit IDC (binf = bt1 = 1) has no MMPP(2) representation, and
 * neither does an IDC profile violating binf > bt1 > 1.
 */

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace fitdetail {

/**
 * Principal branch of the Lambert W function, W0(z) for z >= -1/e.
 *
 * Halley iteration on w e^w = z. The initialization is the branch-point
 * expansion near z = -1/e, where the derivative of w e^w vanishes and a naive
 * start converges slowly, and log(1 + z) elsewhere.
 */
template <class T>
T lambertw0(const T& z, unsigned max_iter) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T e = num_exp(one);
    const T ez1 = e * z + one;
    if (ez1 < zero) throw NumericError("lambertw0: argument below the branch point -1/e");
    if (z == zero) return zero;
    T w;
    if (ez1 < num_traits<T>::from_double(0.3)) {
        const T p = num_sqrt(T(two * ez1));
        w = -one + p - p * p / num_traits<T>::from_int(3) +
            num_traits<T>::from_int(11) * p * p * p / num_traits<T>::from_int(72);
    } else if (z > zero) {
        w = num_log(T(one + z));
    } else {
        w = z;
    }
    for (unsigned it = 0; it < max_iter; ++it) {
        const T ew = num_exp(w);
        const T f = w * ew - z;
        if (f == zero) break;
        const T d1 = ew * (w + one);
        // w=-1 branch point rationale: see _kb/03-api-layer.md (cpp port notes: mam)
        const T dd = two * w + two;
        T den = d1;
        if (dd != zero) den = d1 - (w + two) * f / dd;
        if (den == zero || !(den == den)) break;
        const T step = f / den;
        w -= step;
        if (num_abs(step) <= num_traits<T>::from_double(1e-40) * (one + num_abs(w))) break;
    }
    return w;
}

}  // namespace fitdetail

/** Result of mmpp2_fitc. */
template <class T>
struct Mmpp2FitcResult {
    Map<T> map;
    bool degenerate;      ///< true when a Poisson process was returned instead
    bool third_moment_ok; ///< false when the third moment had to be dropped for feasibility
};

/**
 * MMPP(2) from the arrival rate, the IDC at t1, t2 and infinity, and the
 * third central moment of the counts at t2.
 */
template <class T>
Mmpp2FitcResult<T> mmpp2_fitc(const T& mu, const T& bt1, const T& bt2, const T& binf, const T& m3t2,
                              const T& t1, const T& t2) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmpp2_fitc requires transcendental arithmetic");
    using fitdetail::lambertw0;
    using fitdetail::num_exp;
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);
    const T half = num_traits<T>::from_rational(1, 2);
    const T tiny = num_traits<T>::from_double(1e-8);

    if (mu <= zero) throw InputError("mmpp2_fitc: non-positive arrival rate");
    if (t1 <= zero || t2 <= zero) throw InputError("mmpp2_fitc: non-positive time scale");

    Mmpp2FitcResult<T> res;
    res.degenerate = false;
    res.third_moment_ok = true;

    if ((num_abs(T(binf - one)) < tiny && num_abs(T(binf - bt1)) < tiny) ||
        !(binf > bt1 && bt1 > one)) {
        res.map = map_exponential(mu);
        res.degenerate = true;
        return res;
    }

    const T c = (binf - one) / (binf - bt1);
    const T zarg = -c * num_exp(T(-c));
    const T w0 = lambertw0(zarg, 200u);
    const T x = (w0 + c) / t1;
    if (x <= zero) throw NumericError("mmpp2_fitc: non-positive switching rate");

    const T k1 = pw(mu, 3) * pw(t2, 3);
    const T k2 = three * mu * mu * (binf - one) * t2 * t2;
    const T k3 = three * mu * (binf - one) / x * t2;
    const T k4 = three * mu / (x * x) * (binf - one) * t2 * num_exp(T(-x * t2));
    const T k5 = six * mu / pw(x, 3) * (binf - one) * (one - num_exp(T(-x * t2)));
    const T g1t2 = m3t2 + three * mu * t2 * (mu * t2 - one) * bt2 +
                   mu * t2 * (mu * t2 - one) * (mu * t2 - two);
    const T hden = (k3 / x) + k4 - k5;
    if (hden == zero) throw NumericError("mmpp2_fitc: degenerate third-moment equation");
    const T h = (g1t2 - k1 - k2 - k3 * (-mu) - k4 * mu * x) / hden;

    T r1, r2, l1, l2;
    if (num_abs(h) < num_traits<T>::from_double(1e-4)) {
        r1 = x / two;
        r2 = x / two;
        const T d = half * num_sqrt(T(two * (binf - one) * mu * x));
        l2 = mu - d;
        l1 = mu + d;
    } else {
        const T y = (binf - one) * mu * pw(x, 3) / (two * h * h);
        r1 = x / two * (one + one / num_sqrt(T(num_traits<T>::from_int(4) * y + one)));
        r2 = x - r1;
        if (r1 < r2) {
            const T tmp = r1;
            r1 = r2;
            r2 = tmp;
        }
        if (r1 == r2) throw NumericError("mmpp2_fitc: degenerate switching rates");
        const T wv = h / (r1 - r2);
        const T w_min = -mu / r1 * (r1 + r2);
        const T w_max = mu / r2 * (r1 + r2);
        if (wv < w_min || wv > w_max) {
            // Third moment dropped to keep the rates non-negative.
            res.third_moment_ok = false;
            const T zz = (binf - one) * pw(x, 3) * mu;
            const T u = x * zz / (two * mu * mu * x * x + zz);
            r1 = u + (x - u) / two;
            r2 = x - r1;
            const T delta = num_sqrt(T(zz / (two * r1 * r2)));
            l2 = mu - r2 / x * delta;
            l1 = l2 + delta;
        } else {
            l2 = mu - h / (r1 - r2) * (r2 / (r1 + r2));
            l1 = h / (r1 - r2) + l2;
        }
    }

    Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 0) = -(r1 + l1);
    m.D0(0, 1) = r1;
    m.D0(1, 0) = r2;
    m.D0(1, 1) = -(r2 + l2);
    m.D1(0, 0) = l1;
    m.D1(1, 1) = l2;
    res.map = m;
    return res;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMPP2_FITC_H
