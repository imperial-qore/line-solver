/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP2_FIT_H
#define LINE_API_MAM_MAP2_FIT_H

/**
 * Explicit inverse characterization of a second-order acyclic MAP
 * (matlab/lib/kpctoolbox/map/map2_fit.m).
 *
 * Implements A. Heindl, G. Horvath, K. Gross, "Explicit inverse
 * characterization of acyclic MAPs of second order": given the first three
 * moments e1, e2, e3 and the lag-1 autocorrelation parameter g2, the two
 * canonical forms (hyperexponential for h2 > 0, hypoexponential for
 * -1/4 <= h2 < 0) are written down directly from the normalized moments
 *   h2 = (r2 - r1^2)/r1^2,  h3 = (r3 r1 - r2^2)/r1^4,  r_k = e_k/k!.
 *
 * e3 may be one of the reference's sentinel values, which pick a third moment
 * from e1, e2 and g2 instead of matching one:
 *   e3 = -1  maximize the range of feasible correlations
 *   e3 = -2  minimum feasible e3
 *   e3 = -3  maximum feasible e3
 *   e3 in (-1, 0)  interpolate between the two extremes with weight |e3|
 * The reference's e3 = -4 draws that interpolation weight from MATLAB's global
 * random stream; it is NOT ported, because it cannot be reproduced without
 * that stream. Pass e3 = -r for the deterministic equivalent with weight r.
 *
 * Gated on transcendental arithmetic: c = sqrt(b^2 + 4 h2^3) is taken in every
 * branch, and the sentinel e3 selections take sqrt(-h2).
 *
 * Divergence from the reference, deliberate: MATLAB's h2 == 0 branch builds
 * the Poisson process for the case h3 == 0 and g2 == 0 but omits the return
 * statement, so control falls through to the final else and the function
 * returns an empty MAP with ERR = 30 -- the fitted process is discarded. This
 * port returns the Poisson process, which is what the branch clearly intends,
 * and reports err = 0.
 */

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of map2_fit: the MAP plus the reference's ERR code. */
template <class T>
struct Map2FitResult {
    Map<T> map;    ///< valid only when has_map is true
    bool has_map;  ///< false when the characteristics are infeasible
    int err;       ///< 0 ok, -1 fitted but structurally infeasible, else MATLAB ERR
    T e3_used;     ///< the third moment actually matched (differs when a sentinel was passed)
};

namespace fitdetail {

/** Shared assembly of the two-parameter (a, d1, d2) canonical AMAP(2) form. */
template <class T>
Map<T> map2_fit_form(const T& r1, const T& h2, const T& h3, const T& b, const T& c, const T& g2,
                     const T& a) {
    const T two = num_traits<T>::from_int(2);
    const T one = num_traits<T>::from_int(1);
    const T den = (one - a) * (two * h2 + b - c) + two * c;
    if (den == num_traits<T>::from_int(0)) throw NumericError("map2_fit: degenerate canonical form");
    const T d1 = ((one - a) * (two * h2 * g2 + b - c) + g2 * (b + c) - (b - c)) / den;
    const T d2 = ((g2 - one) * (b - c)) / den;
    const T s = one / (two * r1 * h3);
    const T u = two * h2 + b - c;
    const T v = two * h2 + b + c;
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, num_traits<T>::from_int(0));
    m.D1 = Matrix<T>(2, 2, num_traits<T>::from_int(0));
    m.D0(0, 0) = s * (-u);
    m.D0(0, 1) = s * (u * (one - a));
    m.D0(1, 1) = s * (-v);
    m.D1(0, 0) = s * (u * d1);
    m.D1(0, 1) = s * (u * (a - d1));
    m.D1(1, 0) = s * (v * d2);
    m.D1(1, 1) = s * (v * (one - d2));
    return m;
}

/** Diagonal (reversible) hyperexponential form of map2_fit. */
template <class T>
Map<T> map2_fit_hyper_diag(const T& r1, const T& h2, const T& h3, const T& b, const T& c,
                           const T& g2) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T four = num_traits<T>::from_int(4);
    if (c == num_traits<T>::from_int(0)) throw NumericError("map2_fit: zero discriminant root");
    const T u = two * h2 + b - c;
    const T v = two * h2 + b + c;
    const T s0 = one / (two * r1 * h3);
    const T s1 = one / (four * r1 * h3);
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, num_traits<T>::from_int(0));
    m.D1 = Matrix<T>(2, 2, num_traits<T>::from_int(0));
    m.D0(0, 0) = s0 * (-u);
    m.D0(1, 1) = s0 * (-v);
    m.D1(0, 0) = s1 * (u * (one - b / c + g2 * (one + b / c)));
    m.D1(0, 1) = s1 * (u * (one + b / c) * (one - g2));
    m.D1(1, 0) = s1 * (v * (one - b / c) * (one - g2));
    m.D1(1, 1) = s1 * (v * (one + b / c + g2 * (one - b / c)));
    return m;
}

}  // namespace fitdetail

/** Fit an AMAP(2) to (e1, e2, e3, g2); see the header comment for e3 sentinels. */
template <class T>
Map2FitResult<T> map2_fit(const T& e1, const T& e2, const T& e3_in, const T& g2) {
    static_assert(num_traits<T>::has_transcendental, "map2_fit requires transcendental arithmetic");
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);
    const T twelve = num_traits<T>::from_int(12);

    Map2FitResult<T> res;
    res.has_map = false;
    res.err = 0;

    const T r1 = e1;
    const T r2 = e2 / two;
    if (r1 == zero) throw InputError("map2_fit: zero first moment");
    const T h2 = (r2 - r1 * r1) / (r1 * r1);

    T e3 = e3_in;
    const T scv = (e2 - e1 * e1) / (e1 * e1);
    const T c32 = three / two;
    if (e3 == -one) {
        if (one <= scv && scv < three) {
            if (g2 < zero) {
                const T h3 = h2 - h2 * h2;
                e3 = twelve * pw(e1, 3) * h2 + six * pw(e1, 3) * h3 +
                     six * pw(e1, 3) * (one + h2 * h2);
            } else {
                e3 = (c32 + num_traits<T>::from_double(1e-3)) * e2 * e2 / e1;
            }
        } else if (three <= scv) {
            e3 = (c32 + num_traits<T>::from_double(1e-3)) * e2 * e2 / e1;
        } else if (zero < scv && scv < one) {
            e3 = (one + num_traits<T>::from_double(1e-10)) *
                 (twelve * pw(e1, 3) * h2 +
                  six * pw(e1, 3) * (h2 * (one - h2 - two * num_sqrt(T(-h2)))) +
                  six * pw(e1, 3) * (one + h2 * h2));
        }
    } else if (e3 == -two) {
        if (one <= scv) {
            e3 = (c32 + num_traits<T>::from_double(1e-6)) * e2 * e2 / e1;
        } else if (zero < scv && scv < one) {
            const T h3 = h2 * (one - h2 - two * num_sqrt(T(-h2)));
            e3 = six * pw(e1, 3) * (h2 * h2 + h3);
        }
    } else if (e3 == -three) {
        if (one <= scv) {
            e3 = num_traits<T>::from_double(1e6);
        } else if (zero < scv && scv < one) {
            const T h3 = h2 * h2;
            e3 = six * pw(e1, 3) * (h2 * h2 + h3);
        }
    } else if (e3 > -one && e3 < zero) {
        const T r = num_abs(e3);
        if (one <= scv) {
            e3 = r * (c32 + num_traits<T>::from_double(1e-6)) * e2 * e2 / e1 +
                 (one - r) * num_traits<T>::from_double(1e6);
        } else if (zero < scv && scv < one) {
            const T h3 = r * h2 * (one - h2 - two * num_sqrt(T(-h2))) + (one - r) * (h2 * h2);
            e3 = six * pw(e1, 3) * (h2 * h2 + h3);
        }
    }
    res.e3_used = e3;

    const T r3 = e3 / six;
    const T h3 = (r3 * r1 - r2 * r2) / pw(r1, 4);
    const T b = h3 + h2 * h2 - h2;
    const T crad = b * b + num_traits<T>::from_int(4) * pw(h2, 3);
    if (crad < zero) throw NumericError("map2_fit: negative discriminant b^2 + 4 h2^3");
    const T c = num_sqrt(crad);

    if (r1 <= zero) {
        res.err = 10;  // mean out of bounds
        return res;
    }

    if (h2 == zero) {
        if (h3 == zero && g2 == zero) {
            // See the header note: MATLAB drops this result by falling through.
            res.map = map_exponential_mean(e1);
            res.has_map = true;
            res.err = 0;
        } else {
            res.err = 20;  // correlated exponential
        }
        return res;
    }

    const T qhypo_lo = h2 * (one - h2 - two * num_sqrt(T(-h2)));
    const bool hypo_region = (num_traits<T>::from_rational(-1, 4) <= h2 && h2 < zero &&
                              qhypo_lo <= h3 && h3 <= -(h2 * h2));

    if (h2 > zero && h3 > zero) {
        if (b >= zero) {
            if ((b - c) / (b + c) <= g2 && g2 < one) {
                res.map = fitdetail::map2_fit_hyper_diag(r1, h2, h3, b, c, g2);
                res.has_map = true;
            } else {
                res.err = 51;
                return res;
            }
        } else {
            if (zero <= g2 && g2 < one) {
                res.map = fitdetail::map2_fit_hyper_diag(r1, h2, h3, b, c, g2);
                res.has_map = true;
            } else if (-(h3 + h2 * h2) / h2 <= g2 && g2 < zero) {
                const T a = (h3 + h2 * h2) / h2;
                res.map = fitdetail::map2_fit_form(r1, h2, h3, b, c, g2, a);
                res.has_map = true;
            } else {
                res.err = 52;
                return res;
            }
        }
    } else if (hypo_region) {
        if (g2 >= zero) {
            const T sq = num_sqrt(T(-h3));
            if (g2 <= -((h2 + sq) * (h2 + sq)) / h2) {
                const T a = (two * h2 + b - c) * (h2 + sq) / (two * h2 * sq);
                res.map = fitdetail::map2_fit_form(r1, h2, h3, b, T(-c), g2, a);
                res.has_map = true;
            } else {
                res.err = 53;
                return res;
            }
        } else {
            if (g2 >= -(h3 + h2 * h2) / h2) {
                const T a = (h3 + h2 * h2) / h2;
                res.map = fitdetail::map2_fit_form(r1, h2, h3, b, T(-c), g2, a);
                res.has_map = true;
            } else {
                res.err = 54;
                return res;
            }
        }
    } else {
        res.err = (h2 > zero && h3 < zero) ? 40 : 30;  // h3 / h2 out of bounds
        return res;
    }

    if (res.has_map && !map_isfeasible(res.map, T(num_traits<T>::from_double(1e-10)))) res.err = -1;
    return res;
}

/** Three-argument form: map2_fit(e1, e2, g2), i.e. e3 selected automatically. */
template <class T>
Map2FitResult<T> map2_fit(const T& e1, const T& e2, const T& g2) {
    return map2_fit(e1, e2, T(-num_traits<T>::from_int(1)), g2);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP2_FIT_H
