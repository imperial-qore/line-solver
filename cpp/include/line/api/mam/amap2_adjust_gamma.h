/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_AMAP2_ADJUST_GAMMA_H
#define LINE_API_MAM_AMAP2_ADJUST_GAMMA_H

/**
 * Nearest AMAP(2)-feasible (M2, M3, GAMMA)
 * (matlab/lib/m3a/m3a/amap2/amap2_adjust_gamma.m).
 *
 * M1 is always feasible on its own, so only the second and third moments and
 * the autocorrelation decay rate GAMMA are adjusted. The reference offers four
 * methods, which trade fidelity between the three characteristics:
 *
 *   1  joint search over (M2, M3, GAMMA) minimizing the weighted relative
 *      deviation, subject to the theoretical AMAP(2) feasibility region;
 *   2  M2 forced feasible in closed form, then a search over (M3, GAMMA);
 *   3  (the reference default) strict priority M2 > M3 > GAMMA: apply the
 *      closed-form APH(2) adjustment of aph2_adjust to (M2, M3), then clamp
 *      GAMMA into the interval that pair admits. No optimizer at all;
 *   4  strict priority M2 > GAMMA > M3: M2 is forced feasible, then M3 is
 *      chosen inside its feasible interval to minimize a weighted objective in
 *      which GAMMA is always clamped to the interval implied by M3.
 *
 * The feasibility region and the GAMMA interval are those of Casale, Zhang
 * and Smirni's AMAP(2) characterization, transcribed from the reference's
 * compute_gamma_bounds and nonlcon_theoretical.
 *
 * ACCEPTANCE CONTRACT (see line/util/auglag.h). Methods 1 and 2 are driven in
 * MATLAB by patternsearch (a direct-search method with a Nelder-Mead search
 * step) and method 4 by PSwarm, a stochastic particle swarm. This port uses
 * the augmented Lagrangian with a simplex inner solver for 1 and 2, and a
 * deterministic multi-start bounded simplex search for 4. None of these
 * reproduce the reference's iterates, and PSwarm in particular is not
 * reproducible even against itself. What is guaranteed and tested is the
 * specification:
 *   1. the returned triple satisfies the AMAP(2) feasibility conditions, so
 *      amap2_fitall_gamma admits at least one solution for it;
 *   2. an input that is already feasible is returned unchanged by every
 *      method (objective zero);
 *   3. the achieved objective -- the reference's own weighted relative
 *      deviation -- is returned in Amap2AdjustGammaResult::objective so it can
 *      be compared against any other optimizer's on the same input.
 *
 * ON THE `constraints` ARGUMENT. The reference selects between two constraint
 * sets: constraints = 2 (default) is the theoretical characterization used
 * here, and constraints = 1 replaces it by the indicator
 * "amap2_fit_decay returns something", a zero/one function of x. That
 * indicator is discontinuous, has zero gradient wherever it is defined, and
 * carries no information about the direction of feasibility; it is usable by
 * patternsearch only because a direct-search poll needs no derivative. It is
 * NOT ported: an augmented Lagrangian on a 0/1 constraint degenerates into an
 * unguided penalty. Requesting it raises UnsupportedError rather than silently
 * substituting the theoretical set (which is what the JAR does -- it accepts
 * the argument and ignores it).
 *
 * Gated on transcendental arithmetic: square roots in the bounds and
 * tolerance-driven optimization.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/amap2_assemble.h"
#include "line/api/mam/aph2_adjust.h"
#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/error.h"
#include "line/util/neldermead.h"

namespace line {
namespace mam {

/** Result of amap2_adjust_gamma. */
template <class T>
struct Amap2AdjustGammaResult {
    T M2a;
    T M3a;
    T GAMMAa;
    T objective;  ///< the reference's weighted relative deviation at the answer
    bool feasible;///< the returned triple satisfies the AMAP(2) conditions
};

namespace fitdetail {

/** The interval of feasible GAMMA for a given (M1, M2a, M3a); compute_gamma_bounds. */
template <class T>
void amap2_gamma_bounds(const T& M1, const T& M2a, const T& M3a, const T& tol, T& lb, T& ub) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);
    const T nine = num_traits<T>::from_int(9);
    const T twelve = num_traits<T>::from_int(12);
    const T eighteen = num_traits<T>::from_int(18);
    const T twentyfour = num_traits<T>::from_int(24);
    const T twentyseven = num_traits<T>::from_int(27);
    const T half = num_traits<T>::from_rational(1, 2);
    const T zero = num_traits<T>::from_int(0);

    const T N2a = M2a / (M1 * M1);
    const T N3a = M3a / (M2a * M1);

    if (N2a < two) {
        lb = -(N2a * (N3a - six) + six) / (three * N2a - six);
        const T disc = N2a * N2a - two * N2a * N3a / three;
        const T root = disc > zero ? num_sqrt(disc) : zero;
        const T inner = half * (N2a - two) + half * root;
        ub = -(two * inner * inner) / (N2a - two);
        ub = ub * (one - tol);
    } else if (N3a < nine - twelve / N2a) {
        lb = -(N2a * (N3a - six) + six) / (three * N2a - six);
        ub = one - tol;
    } else {
        const T inner = N2a * (eighteen * N2a + N3a * (N3a - eighteen) - twentyseven) +
                        twentyfour * N3a;
        const T arg = N2a * inner;
        const T x1 = arg > zero ? num_sqrt(arg) : zero;
        const T x2 = N2a * (N3a - nine);
        lb = (x2 - x1 + twelve) / (x2 + x1 + twelve);
        ub = one - tol;
    }
}

/**
 * The reference's nonlcon_theoretical: five inequalities c <= 0 expressing
 * AMAP(2) feasibility of (M2, M3, GAMMA) at a fixed M1. Rows 2 and 3 are the
 * third-moment bounds, whose form differs on either side of N2 = 2; rows 4 and
 * 5 are the GAMMA interval.
 *
 * REFERENCE DEFECT (matlab/lib/m3a/m3a/amap2/amap2_adjust_gamma.m,
 * nonlcon_theoretical): the third-moment rows are written as
 *     if 3/2 <= n2 && n2 < 2   ... elseif n2 > 2 ...
 * with no else, and c is preallocated to zeros(5,1). For n2 BELOW 3/2, and at
 * n2 EXACTLY 2, neither branch runs and rows 2 and 3 stay at zero -- that is,
 * the third moment is left completely unconstrained. Row 1 (3/2 - n2 <= 0) is
 * the only thing pushing n2 up, and it is satisfied as soon as the deficit is
 * below the solver's constraint tolerance, at which point any n3 whatsoever is
 * declared feasible. Reproduction with this port's optimizer, before the
 * repair below: minimizing the weighted deviation from
 * (M1, M2, M3, GAMMA) = (1, 1.3, 2, 0.5) returns
 * (M2a, M3a, GAMMAa) = (1.4999999999870, 1.19265, 0.49871) with objective
 * 1.591 -- n2 is 1.3e-11 short of 3/2, so rows 2 and 3 vanish and the third
 * moment settles at n3 = 0.795 where the only APH(2) with n2 = 3/2 has
 * n3 = 2 exactly. The objective beats every genuinely feasible point.
 *
 * The repair is to evaluate the third-moment bounds at n2 CLAMPED into the
 * region, max(n2, 3/2), and to use the n2 > 2 form at n2 = 2 as well, so that
 * the rows are defined everywhere and continuous across both junctions. Row 1
 * still carries the n2 >= 3/2 requirement itself. MATLAB was NOT edited.
 */
template <class T>
std::vector<T> amap2_feasibility(const T& M1, const T& xM2, const T& xM3, const T& xGAMMA) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);
    const T nine = num_traits<T>::from_int(9);
    const T twelve = num_traits<T>::from_int(12);
    const T eighteen = num_traits<T>::from_int(18);
    const T twentyfour = num_traits<T>::from_int(24);
    const T twentyseven = num_traits<T>::from_int(27);
    const T half = num_traits<T>::from_rational(1, 2);
    const T threehalf = num_traits<T>::from_rational(3, 2);
    const T eps = num_traits<T>::from_double(1e-8);

    const T n2 = xM2 / (M1 * M1);
    const T n3 = xM3 / (M1 * xM2);

    std::vector<T> c(5, zero);
    c[0] = threehalf - n2;

    // n2 clamped into the region, so that rows 1 and 2 below are defined for
    // every argument; see the reference-defect note above
    const T n2c = n2 < threehalf ? threehalf : n2;
    if (n2c < two) {
        // the reference's p2, a2, l2, u2; all real only for n2 < 2, which is
        // why they are evaluated inside this branch alone
        const T p2 = three * (n2c - two) / (three * n2c) *
                     (-two * num_sqrt(three) / num_sqrt(T(twelve - six * n2c)) - one);
        const T inner = p2 * p2 + two * p2 * (n2c - two);
        const T a2 = (n2c - two) / (p2 * (one - n2c) + (inner > zero ? num_sqrt(inner) : zero));
        const T l2 = three * (a2 + one) / (a2 * p2 + one) -
                     six * a2 / (two + a2 * p2 * (two * a2 + two));
        const T u2 = six * (n2c - one) / n2c;
        c[1] = l2 - n3;
        c[2] = n3 - u2;
    } else {
        c[1] = threehalf * n2c - n3 + eps;
        c[2] = zero;
    }

    const T lb1 = -(n2 * (n3 - six) + six) / (three * n2 - six);
    const T disc = n2 * n2 - two * n2 * n3 / three;
    const T root1 = disc > zero ? num_sqrt(disc) : zero;
    const T innerub = half * (n2 - two) + half * root1;
    const T ub1 = -(two * innerub * innerub) / (n2 - two);
    const T tmp1 = n2 * (n3 - nine);
    const T arg2 = n2 * (n2 * (eighteen * n2 + n3 * (n3 - eighteen) - twentyseven) +
                         twentyfour * n3);
    const T tmp2 = arg2 > zero ? num_sqrt(arg2) : zero;
    const T lb2 = (tmp1 - tmp2 + twelve) / (tmp1 + tmp2 + twelve);

    if (n2 < two) {
        c[3] = lb1 - xGAMMA;
        c[4] = xGAMMA - ub1;
    } else if (n3 < nine - twelve / n2) {
        c[3] = lb1 - xGAMMA;
        c[4] = xGAMMA - one;
    } else {
        c[3] = lb2 - xGAMMA;
        c[4] = xGAMMA - one;
    }
    return c;
}

}  // namespace fitdetail

/**
 * Is (M1, M2, M3, GAMMA) AMAP(2)-feasible? The predicate behind the
 * acceptance criterion of every method here; it is the reference's own
 * nonlcon_theoretical, evaluated rather than optimized.
 */
template <class T>
bool amap2_gamma_feasible(const T& M1, const T& M2, const T& M3, const T& GAMMA, const T& slack) {
    const std::vector<T> c = fitdetail::amap2_feasibility(M1, M2, M3, GAMMA);
    for (std::size_t i = 0; i < c.size(); ++i)
        if (c[i] > slack) return false;
    return true;
}

/**
 * Nearest AMAP(2)-feasible characteristics.
 *
 * @param M1,M2,M3,GAMMA the requested characteristics
 * @param weights        the three weights on (M2, M3, GAMMA); the reference
 *                       default is (10, 1, 10). Ignored by methods 3 and 4's
 *                       first component, exactly as in the reference
 * @param method         1..4, see the header comment; the reference default 3
 * @param constraints    2 for the theoretical region (the only one ported);
 *                       1 raises UnsupportedError
 * @param tol            the strict-inequality slack, the reference's 1e-2
 */
template <class T>
Amap2AdjustGammaResult<T> amap2_adjust_gamma(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                                             const std::vector<T>& weights, int method,
                                             int constraints, const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "amap2_adjust_gamma requires transcendental arithmetic");
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T half = num_traits<T>::from_rational(1, 2);
    const T threehalf = num_traits<T>::from_rational(3, 2);

    if (M1 <= zero) throw InputError("amap2_adjust_gamma: non-positive first moment");
    if (weights.size() != 3) throw InputError("amap2_adjust_gamma: three weights are required");
    if (constraints == 1)
        throw UnsupportedError(
            "amap2_adjust_gamma: constraints = 1 uses a 0/1 feasibility indicator that only a "
            "direct-search poll can exploit; not ported");
    if (constraints != 2) throw InputError("amap2_adjust_gamma: constraints must be 1 or 2");

    // the reference's fun / fun2 divide by the targets
    if (method != 3 && (M2 == zero || M3 == zero || GAMMA == zero))
        throw InputError(
            "amap2_adjust_gamma: methods 1, 2 and 4 minimize a RELATIVE deviation, which is "
            "undefined when a target is zero; use method 3");

    Amap2AdjustGammaResult<T> r;

    if (method == 3) {
        // priorities M2 > M3 > GAMMA; entirely closed form
        const Aph2AdjustResult<T> ad = aph2_adjust(M1, M2, M3);
        r.M2a = ad.M2a;
        r.M3a = ad.M3a;
        T lb = zero, ub = zero;
        fitdetail::amap2_gamma_bounds(M1, r.M2a, r.M3a, tol, lb, ub);
        r.GAMMAa = GAMMA < lb ? lb : (GAMMA > ub ? ub : GAMMA);
    } else if (method == 1) {
        // joint search over (M2, M3, GAMMA)
        const T w0 = weights[0], w1 = weights[1], w2 = weights[2];
        auto f = [M2, M3, GAMMA, w0, w1, w2](const std::vector<T>& x) {
            const T a0 = (x[0] - M2) / M2 * w0;
            const T a1 = (x[1] - M3) / M3 * w1;
            const T a2 = (x[2] - GAMMA) / GAMMA * w2;
            using std::sqrt;
            return sqrt(T(a0 * a0 + a1 * a1 + a2 * a2));
        };
        auto g = [M1](const std::vector<T>& x) {
            return fitdetail::amap2_feasibility(M1, x[0], x[1], x[2]);
        };

        // the reference's analytically feasible starting point
        const Map<T> feas =
            GAMMA > zero
                ? map_scale(amap2_assemble(one, T(one / three), half, T(two / three), 1), M1)
                : map_scale(amap2_assemble(one, T(two / three), half, T(two / three), 2), M1);
        std::vector<unsigned> lags(2);
        lags[0] = 3;
        lags[1] = 4;
        const std::vector<T> acf = map_acf(feas, lags);
        std::vector<T> x0(3);
        x0[0] = map_moment(feas, 2);
        x0[1] = map_moment(feas, 3);
        x0[2] = acf[1] / acf[0];

        // deterministic second start rationale: see _kb/03-api-layer.md (cpp port notes: mam)
        const Aph2AdjustResult<T> ad3 = aph2_adjust(M1, M2, M3);
        T lb3 = zero, ub3 = zero;
        fitdetail::amap2_gamma_bounds(M1, ad3.M2a, ad3.M3a, tol, lb3, ub3);
        std::vector<T> x1(3);
        x1[0] = ad3.M2a;
        x1[1] = ad3.M3a;
        x1[2] = GAMMA < lb3 ? lb3 : (GAMMA > ub3 ? ub3 : GAMMA);

        std::vector<Bound<T>> bnd(3);
        bnd[0] = bound_lower(zero);
        bnd[1] = bound_lower(zero);
        bnd[2] = bound_box(T(-one), T(one - tol));

        AugLagOptions<T> opt = auglag_defaults<T>();
        opt.ctol = num_traits<T>::from_double(1e-10);
        const AugLagResult<T> s0 = auglag(f, NoConstraints<T>(), g, x0, bnd, opt);
        const AugLagResult<T> s1 = auglag(f, NoConstraints<T>(), g, x1, bnd, opt);
        const bool f0 = s0.violation <= opt.ctol;
        const bool f1 = s1.violation <= opt.ctol;
        const bool take1 = (f1 && !f0) || (f1 == f0 && s1.fval < s0.fval);
        const AugLagResult<T>& sol = take1 ? s1 : s0;
        r.M2a = sol.x[0];
        r.M3a = sol.x[1];
        r.GAMMAa = sol.x[2];
    } else if (method == 2) {
        // M2 forced feasible in closed form, then a search over (M3, GAMMA)
        r.M2a = M2 > threehalf * M1 * M1 ? M2 : T(threehalf * M1 * M1);
        // REFERENCE DEFECT (method 2 interval): see _kb/03-api-layer.md (cpp port notes: mam)
        const T n2 = r.M2a / (M1 * M1);
        T M3_LB, M3_UB, FM3;
        bool ub_finite = true;
        if (n2 >= threehalf && n2 < two) {
            const T p2 = three * (n2 - two) / (three * n2) *
                         (-two * fitdetail::num_sqrt(three) /
                              fitdetail::num_sqrt(T(num_traits<T>::from_int(12) -
                                                    num_traits<T>::from_int(6) * n2)) -
                          one);
            const T inner = p2 * p2 + two * p2 * (n2 - two);
            const T a2 = (n2 - two) /
                         (p2 * (one - n2) + (inner > zero ? fitdetail::num_sqrt(inner) : zero));
            const T l2 = three * (a2 + one) / (a2 * p2 + one) -
                         num_traits<T>::from_int(6) * a2 /
                             (two + a2 * p2 * (two * a2 + two));
            const T u2 = num_traits<T>::from_int(6) * (n2 - one) / n2;
            FM3 = (l2 + (u2 - l2) / two) * M1 * r.M2a;
            M3_LB = l2 * M1 * r.M2a;
            M3_UB = u2 * M1 * r.M2a;
        } else {
            FM3 = threehalf * r.M2a * r.M2a / M1 + tol;
            M3_LB = FM3;
            M3_UB = FM3;  // unused, the bound is one-sided
            ub_finite = false;
        }

        const T w1 = weights[1], w2 = weights[2];
        auto f = [M3, GAMMA, w1, w2](const std::vector<T>& x) {
            const T a1 = (x[0] - M3) / M3 * w1;
            const T a2 = (x[1] - GAMMA) / GAMMA * w2;
            using std::sqrt;
            return sqrt(T(a1 * a1 + a2 * a2));
        };
        const T M2a_c = r.M2a;
        auto g = [M1, M2a_c](const std::vector<T>& x) {
            return fitdetail::amap2_feasibility(M1, M2a_c, x[0], x[1]);
        };

        std::vector<T> x0(2);
        x0[0] = FM3;
        x0[1] = zero;  // a null decay rate is always feasible
        std::vector<Bound<T>> bnd(2);
        bnd[0] = ub_finite ? bound_box(M3_LB, M3_UB) : bound_lower(M3_LB);
        bnd[1] = bound_box(T(-one), T(one - tol));

        AugLagOptions<T> opt = auglag_defaults<T>();
        opt.ctol = num_traits<T>::from_double(1e-10);
        const AugLagResult<T> sol = auglag(f, NoConstraints<T>(), g, x0, bnd, opt);
        r.M3a = sol.x[0];
        r.GAMMAa = sol.x[1];
    } else if (method == 4) {
        // priorities M2 > GAMMA > M3
        const T M1sq = M1 * M1;
        const T scv = (M2 - M1sq) / M1sq;
        T scva;
        if (scv < half) {
            r.M2a = threehalf * M1sq;
            scva = (r.M2a - M1sq) / M1sq;
        } else {
            r.M2a = M2;
            scva = scv;
        }
        T M3_lb, M3_ub;
        bool ub_finite;
        if (scva <= one) {
            const T d = one - scva;
            M3_lb = three * fitdetail::pw(M1, 3) *
                    (three * scva - one + fitdetail::num_sqrt(two) * d * fitdetail::num_sqrt(d));
            M3_ub = num_traits<T>::from_int(6) * fitdetail::pw(M1, 3) * scva;
            ub_finite = true;
        } else {
            M3_lb = threehalf * fitdetail::pw(M1, 3) * (one + scva) * (one + scva);
            M3_ub = M3_lb;
            ub_finite = false;
        }

        if (ub_finite && num_abs(T(M3_lb - M3_ub)) < tol) {
            r.M3a = (M3_lb + M3_ub) / two;
            r.GAMMAa = zero;
        } else {
            const T w1 = weights[1], w2 = weights[2];
            const T M2a_c = r.M2a;
            auto obj = [M1, M2a_c, M3, GAMMA, w1, w2, tol](const std::vector<T>& x) {
                T lb = num_traits<T>::from_int(0), ub = num_traits<T>::from_int(0);
                fitdetail::amap2_gamma_bounds(M1, M2a_c, x[0], tol, lb, ub);
                const T ga = GAMMA < lb ? lb : (GAMMA > ub ? ub : GAMMA);
                const T a1 = x[0] / M3 - num_traits<T>::from_int(1);
                const T a2 = ga / GAMMA - num_traits<T>::from_int(1);
                return T(w1 * a1 * a1 + w2 * a2 * a2);
            };

            // deterministic multi-start rationale: see _kb/03-api-layer.md (cpp port notes: mam)
            const T lo = M3_lb + tol;
            std::vector<T> starts;
            const unsigned nstart = 9;
            for (unsigned k = 0; k < nstart; ++k) {
                const T frac = num_traits<T>::from_int(long(k)) /
                               num_traits<T>::from_int(long(nstart - 1));
                if (ub_finite)
                    starts.push_back(T(lo + frac * (M3_ub - lo)));
                else
                    starts.push_back(T(lo * (one + frac * num_traits<T>::from_int(9))));
            }

            std::vector<Bound<T>> bnd(1);
            bnd[0] = ub_finite ? bound_box(lo, M3_ub) : bound_lower(lo);

            NelderMeadOptions<T> nmo = nelder_mead_defaults<T>();
            nmo.xtol = num_traits<T>::from_double(1e-12);
            nmo.ftol = num_traits<T>::from_double(1e-16);

            T bestf = zero;
            T bestx = lo;
            bool have = false;
            for (std::size_t k = 0; k < starts.size(); ++k) {
                std::vector<T> x0(1, starts[k]);
                const NelderMeadResult<T> s = nelder_mead_box(obj, x0, bnd, nmo);
                if (!have || s.fval < bestf) {
                    bestf = s.fval;
                    bestx = s.x[0];
                    have = true;
                }
            }
            r.M3a = bestx;
            T lb = zero, ub = zero;
            fitdetail::amap2_gamma_bounds(M1, r.M2a, r.M3a, tol, lb, ub);
            r.GAMMAa = GAMMA < lb ? lb : (GAMMA > ub ? ub : GAMMA);
        }
    } else {
        throw InputError("amap2_adjust_gamma: method must be 1, 2, 3 or 4");
    }

    // the reference's own objective, reported for comparison
    if (M2 == zero || M3 == zero || GAMMA == zero) {
        r.objective = zero;
    } else {
        const T a0 = (r.M2a - M2) / M2 * weights[0];
        const T a1 = (r.M3a - M3) / M3 * weights[1];
        const T a2 = (r.GAMMAa - GAMMA) / GAMMA * weights[2];
        r.objective = sqrt(T(a0 * a0 + a1 * a1 + a2 * a2));
    }
    r.feasible = amap2_gamma_feasible(M1, r.M2a, r.M3a, r.GAMMAa,
                                      T(num_traits<T>::from_double(1e-8)));
    return r;
}

/** amap2_adjust_gamma with the reference defaults: weights (10,1,10), method 3, constraints 2. */
template <class T>
Amap2AdjustGammaResult<T> amap2_adjust_gamma(const T& M1, const T& M2, const T& M3,
                                             const T& GAMMA) {
    std::vector<T> w(3);
    w[0] = num_traits<T>::from_int(10);
    w[1] = num_traits<T>::from_int(1);
    w[2] = num_traits<T>::from_int(10);
    return amap2_adjust_gamma(M1, M2, M3, GAMMA, w, 3, 2, T(num_traits<T>::from_double(1e-2)));
}

/** amap2_adjust_gamma with a chosen method and the remaining reference defaults. */
template <class T>
Amap2AdjustGammaResult<T> amap2_adjust_gamma(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                                             int method) {
    std::vector<T> w(3);
    w[0] = num_traits<T>::from_int(10);
    w[1] = num_traits<T>::from_int(1);
    w[2] = num_traits<T>::from_int(10);
    return amap2_adjust_gamma(M1, M2, M3, GAMMA, w, method, 2,
                              T(num_traits<T>::from_double(1e-2)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_AMAP2_ADJUST_GAMMA_H
