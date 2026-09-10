/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMPP2_FITC_APPROX_H
#define LINE_API_MAM_MMPP2_FITC_APPROX_H

/**
 * MMPP(2) fitted to counting-process characteristics by optimization
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fitc_approx.m).
 *
 * Where mmpp2_fitc.h inverts the Heffes-Lucantoni relations in closed form
 * and therefore fails outright when the requested characteristics are not
 * exactly representable, this routine minimizes the mismatch. The decision
 * variables are the two phase arrival rates l1, l2 and the two switching
 * rates r1, r2 of
 *
 *   D0 = [ -(l1 + r1)   r1 ; r2   -(l2 + r2) ],   D1 = diag(l1, l2)
 *
 * and the objective is the sum of squared RELATIVE errors of the five
 * characteristics (rate a, IDC at t1 and t2, IDC at infinity, third central
 * moment of the counts at t2), each computed from the closed forms of the
 * reference rather than from a numerical counting-process evaluation. As in
 * the reference, the characteristics are evaluated on the time scale
 * stretched by factor = a/xa, so the objective is invariant to the rate, and
 * the returned MAP is finally rescaled with map_scale to have rate exactly a.
 *
 * ACCEPTANCE CONTRACT (see line/util/levmar.h). MATLAB drives this with the
 * problem-based `solve`, i.e. fmincon interior-point over an
 * automatically differentiated fcn2optimexpr expression; the JAR uses Apache
 * Commons Math BOBYQA. This port uses the augmented Lagrangian of
 * line/util/auglag.h with a Levenberg-Marquardt inner solver, since the
 * objective is literally a sum of squares. None of the three reproduce each
 * other's iterates. What is guaranteed and tested is the specification:
 *   1. the returned (D0, D1) is a valid MAP -- row sums zero, non-negative
 *      off-diagonals of D0, non-negative D1 -- verified with map_isfeasible;
 *   2. its rate is exactly a, by construction (map_scale);
 *   3. the achieved objective, the same sum of squared relative errors the
 *      reference minimizes, is returned in Mmpp2FitcApproxResult::objective so
 *      it can be compared against any other optimizer's on the same input;
 *   4. when the target characteristics are those of an actual MMPP(2), the
 *      objective reaches zero to the stated tolerance and the achieved
 *      characteristics, recomputed independently with map_count_var and
 *      map_count_moment, reproduce the targets.
 *
 * The two bound constraints of the reference (l1, l2 >= 1e-6 and r1, r2 >= 0)
 * are passed as inequality rows. Note that the characteristics are undefined
 * at r1 = r2 = 0 and at l1 r2 + l2 r1 = 0, where the expressions divide by
 * zero; no guard is needed, because a step producing a non-finite residual
 * compares false against the incumbent sum of squares and is therefore
 * rejected by the Levenberg-Marquardt loop, which then increases its damping.
 *
 * REFERENCE DEFECT (matlab/lib/kpctoolbox/mmpp/mmpp2_fitc_approx.m): when
 * t1 == t2 the local xbt2 is never assigned, yet it is used unconditionally
 * two lines later in the expression for xm3t2. Calling the reference with
 * t1 == t2 therefore raises "Unrecognized function or variable 'xbt2'". This
 * port sets xbt2 = xbt1 in that case, which is what the value means, and drops
 * the duplicated IDC residual, matching the JAR. MATLAB was NOT edited.
 *
 * Gated on transcendental arithmetic: exponentials in the IDC and third-moment
 * expressions, and tolerance-driven optimization.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of mmpp2_fitc_approx. */
template <class T>
struct Mmpp2FitcApproxResult {
    Map<T> map;      ///< the fitted MMPP(2), rescaled to rate a
    T objective;     ///< sum of squared relative errors at the optimum
    T violation;     ///< worst bound violation at the optimum
    bool converged;  ///< the constrained solve reached feasibility within its caps
};

namespace fitdetail {

/**
 * Closed-form counting characteristics of an MMPP(2) with parameters
 * (l1, l2, r1, r2), transcribed from the reference's compute_obj. The product
 * exp(-r1 T - r2 T) is written as a single exponential of -(r1 + r2) T.
 */
template <class T>
struct Mmpp2Chars {
    T xa;     ///< arrival rate
    T xbt1;   ///< IDC at the (stretched) scale t1
    T xbt2;   ///< IDC at the (stretched) scale t2
    T xbinf;  ///< limiting IDC
    T xm3t2;  ///< third central moment of the counts at the stretched t2
};

template <class T>
T mmpp2_idc_at(const T& l1, const T& l2, const T& r1, const T& r2, const T& tt) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T four = num_traits<T>::from_int(4);
    const T E = num_exp(T(-(r1 + r2) * tt));
    const T inner = r1 * (two * l1 * l1 * r2 * r2 * tt - two * l2 * l2 * r2 - two * l1 * l1 * r2 +
                          two * l2 * l2 * r2 * r2 * tt + four * l1 * l2 * r2 +
                          two * l1 * l1 * r2 * E + two * l2 * l2 * r2 * E -
                          four * l1 * l2 * r2 * r2 * tt - four * l1 * l2 * r2 * E) +
                    r1 * r1 * (two * r2 * tt * l1 * l1 - four * r2 * tt * l1 * l2 +
                               two * r2 * tt * l2 * l2);
    return inner / (tt * pw(T(r1 + r2), 3) * (l1 * r2 + l2 * r1)) + one;
}

template <class T>
Mmpp2Chars<T> mmpp2_chars(const T& l1, const T& l2, const T& r1, const T& r2, const T& a,
                          const T& t1, const T& t2, bool same_scale) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T four = num_traits<T>::from_int(4);
    const T six = num_traits<T>::from_int(6);

    Mmpp2Chars<T> c;
    c.xa = (l1 * r2 + l2 * r1) / (r1 + r2);
    const T factor = a / c.xa;

    c.xbt1 = mmpp2_idc_at(l1, l2, r1, r2, T(t1 * factor));
    c.xbt2 = same_scale ? c.xbt1 : mmpp2_idc_at(l1, l2, r1, r2, T(t2 * factor));

    c.xbinf = ((two * r2 * l1 * l1 - four * r2 * l1 * l2 + two * r2 * l2 * l2) * r1 * r1 +
               (two * l1 * l1 * r2 * r2 - four * l1 * l2 * r2 * r2 + two * l2 * l2 * r2 * r2) * r1) /
                  (pw(T(r1 + r2), 3) * (l1 * r2 + l2 * r1)) +
              one;

    const T t = t2 * factor;
    const T d = r1 + r2;
    const T p = (l1 - l2) * (r1 - r2);
    const T Etd = num_exp(T(-t * d));
    const T bm1 = c.xbinf - one;
    const T xg3t = pw(c.xa, 3) * pw(t, 3) + three * c.xa * c.xa * bm1 * t * t +
                   three * c.xa * bm1 / d * (p / d - c.xa) * t +
                   three * c.xa / (d * d) * bm1 * (p + c.xa * d) * t * Etd -
                   six * c.xa / pw(d, 3) * bm1 * p * (one - Etd);
    c.xm3t2 = xg3t - three * c.xa * t * (c.xa * t - one) * c.xbt2 -
              c.xa * t * (c.xa * t - one) * (c.xa * t - two);
    return c;
}

}  // namespace fitdetail

/**
 * Fit an MMPP(2) to counting characteristics.
 *
 * @param a     arrival rate
 * @param bt1   IDC at time scale t1
 * @param bt2   IDC at time scale t2
 * @param binf  limiting IDC
 * @param m3t2  third central moment of the counts at t2
 * @param t1,t2 the two time scales
 * @param opt   tuning of the constrained solve
 */
template <class T>
Mmpp2FitcApproxResult<T> mmpp2_fitc_approx(const T& a, const T& bt1, const T& bt2, const T& binf,
                                           const T& m3t2, const T& t1, const T& t2,
                                           const AugLagOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmpp2_fitc_approx requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (a <= zero) throw InputError("mmpp2_fitc_approx: non-positive arrival rate");
    if (t1 <= zero || t2 <= zero) throw InputError("mmpp2_fitc_approx: non-positive time scale");
    if (bt1 == zero || bt2 == zero || binf == zero || m3t2 == zero)
        throw InputError("mmpp2_fitc_approx: a target characteristic is zero, so the relative "
                         "objective of the reference is undefined");

    const bool same_scale = t1 == t2;
    const std::size_t nres = same_scale ? 4 : 5;

    auto resid = [a, bt1, bt2, binf, m3t2, t1, t2, same_scale, one](const std::vector<T>& x) {
        const fitdetail::Mmpp2Chars<T> c =
            fitdetail::mmpp2_chars(x[0], x[1], x[2], x[3], a, t1, t2, same_scale);
        std::vector<T> r;
        r.push_back(T(c.xa / a - one));
        r.push_back(T(c.xbt1 / bt1 - one));
        if (!same_scale) r.push_back(T(c.xbt2 / bt2 - one));
        r.push_back(T(c.xbinf / binf - one));
        r.push_back(T(c.xm3t2 / m3t2 - one));
        return r;
    };

    // the reference's prob.Constraints: l1, l2 >= 1e-6 and r1, r2 >= 0
    const T lmin = num_traits<T>::from_double(1e-6);
    auto g = [lmin, zero](const std::vector<T>& x) {
        std::vector<T> gv(4);
        gv[0] = lmin - x[0];
        gv[1] = lmin - x[1];
        gv[2] = zero - x[2];
        gv[3] = zero - x[3];
        return gv;
    };

    std::vector<T> x0(4);
    x0[0] = a * num_traits<T>::from_rational(3, 4);  // the reference's guess
    x0[1] = a * num_traits<T>::from_rational(3, 2);
    x0[2] = num_traits<T>::from_rational(1, 3);
    x0[3] = num_traits<T>::from_rational(2, 3);

    const AugLagResult<T> sol = auglag_ls(resid, nres, NoConstraints<T>(), g, x0, opt);

    const T l1 = sol.x[0];
    const T l2 = sol.x[1];
    const T r1 = sol.x[2];
    const T r2 = sol.x[3];

    Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 0) = -(l1 + r1);
    m.D0(0, 1) = r1;
    m.D0(1, 0) = r2;
    m.D0(1, 1) = -(l2 + r2);
    m.D1(0, 0) = l1;
    m.D1(1, 1) = l2;

    Mmpp2FitcApproxResult<T> res;
    res.map = map_scale(m, T(one / a));  // force the rate to a, as the reference does
    res.objective = sol.fval;
    res.violation = sol.violation;
    res.converged = sol.violation <= opt.ctol;
    return res;
}

/** mmpp2_fitc_approx with the default tuning. */
template <class T>
Mmpp2FitcApproxResult<T> mmpp2_fitc_approx(const T& a, const T& bt1, const T& bt2, const T& binf,
                                           const T& m3t2, const T& t1, const T& t2) {
    AugLagOptions<T> opt = auglag_defaults<T>();
    opt.ctol = num_traits<T>::from_double(1e-12);
    return mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, opt);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMPP2_FITC_APPROX_H
