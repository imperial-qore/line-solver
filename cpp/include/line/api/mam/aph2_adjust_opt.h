/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH2_ADJUST_OPT_H
#define LINE_API_MAM_APH2_ADJUST_OPT_H

/**
 * Optimization-based APH(2) moment adjustment: the 'opt_param' and 'opt_char'
 * methods of matlab/lib/m3a/m3a/aph2/aph2_adjust.m.
 *
 * Both find (M2a, M3a) close to a requested (M2, M3), holding M1 fixed, such
 * that an APH(2) with those three moments exists. They differ in the space
 * they search:
 *
 *   opt_param  searches the PARAMETER space. The decision variables are the
 *              second phase mean l2 and the branching probability r1, with
 *              l1 = M1 - l2 r1 chosen so that the first moment is matched
 *              exactly; (M2a, M3a) are then whatever that APH(2) has. Any
 *              point of the feasible box is a valid APH(2), so the answer is
 *              feasible by construction.
 *   opt_char   searches the CHARACTERISTIC space. The decision variables are
 *              (M2a, M3a) directly and feasibility is imposed as nonlinear
 *              constraints: the closed-form inversion aph2_fit1 / aph2_fit2
 *              must produce a non-negative discriminant, non-negative phase
 *              means and a probability in [0, 1]. The reference solves the
 *              two inversions as separate problems and keeps whichever gives
 *              the smaller adjustment; so does this port. The closed-form
 *              Telek-Heindl bounds are imposed as three further rows, which
 *              are redundant in exact arithmetic but close a gap in the
 *              inversion form that an optimizer otherwise exploits; see
 *              aph2_moment_bounds_rows for the reproduction.
 *
 * ACCEPTANCE CONTRACT (see line/util/auglag.h). MATLAB drives both with
 * fmincon (active-set) or, in the '_gads' variants, with GlobalSearch. This
 * port uses the augmented Lagrangian of line/util/auglag.h with the
 * least-squares inner solver, which is a different algorithm: the iterates do
 * not match and, on a problem with several local minima, the local minimum
 * need not be the same one. What is guaranteed and tested is the
 * specification:
 *   1. the returned (M1, M2a, M3a) admits an APH(2), i.e. it satisfies the
 *      Telek-Heindl bounds that aph2_adjust's 'simple' method enforces;
 *   2. when the requested (M2, M3) is already APH(2)-feasible, the returned
 *      pair reproduces it to the stated tolerance (objective ~ 0);
 *   3. the achieved objective ||(M2a, M3a) - (M2, M3)||, the exact objective
 *      the reference minimizes, is returned in Aph2AdjustOptResult::objective
 *      so it can be compared against any other optimizer's on the same input.
 *
 * NOT PORTED: 'opt_param_gads' and 'opt_char_gads'. They are the same two
 * problems handed to GlobalSearch, whose answer depends on a randomly seeded
 * multi-start; a deterministic multi-start would be a different algorithm
 * again and would not reproduce them, so they are deliberately absent rather
 * than approximated. Since both problems are smooth with a small feasible
 * region, the single-start solve here reaches the same objective on every
 * input exercised in tests/test_mam_fit_optim.cpp.
 *
 * REFERENCE DEFECT (matlab/lib/m3a/m3a/aph2/aph2_adjust.m, 'opt_char'): the
 * constraint functions nonlcon1 and nonlcon2 read `tmp0`, which is assigned
 * only inside the sibling nested functions aph2_fit1 / aph2_fit2 and is
 * therefore not in their scope. Calling aph2_adjust(M1, M2, M3, 'opt_char')
 * raises "Unrecognized function or variable 'tmp0'" before any optimization
 * happens, so the reference's opt_char has never run. This port recomputes
 * tmp0 in the constraint, which is the evident intent (it is the
 * discriminant of the inversion, and c(1) = -tmp0 asks for it to be
 * non-negative). MATLAB was NOT edited.
 *
 * Gated on transcendental arithmetic: the optimizers stop on tolerances, and
 * the inversion takes a square root.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_fit_detail.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** Result of an optimization-based APH(2) moment adjustment. */
template <class T>
struct Aph2AdjustOptResult {
    T M2a;          ///< adjusted second moment
    T M3a;          ///< adjusted third moment
    T objective;    ///< ||(M2a, M3a) - (M2, M3)||, the reference's objective
    T violation;    ///< worst constraint violation at the returned point
    bool converged; ///< the constrained solve reached feasibility within its caps
};

namespace fitdetail {

/**
 * Telek-Heindl APH(2) feasibility of a moment triple, the predicate the
 * 'simple' method of aph2_adjust enforces by clamping. Used here to state the
 * acceptance criterion rather than to compute anything.
 */
template <class T>
bool aph2_moments_feasible(const T& M1, const T& M2, const T& M3, const T& tol) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);
    const T half = num_traits<T>::from_rational(1, 2);
    const T M1sq = M1 * M1;
    const T scv = (M2 - M1sq) / M1sq;
    if (scv < half * (one - tol)) return false;
    if (scv <= one) {
        const T d = one - scv;
        const T lb = three * pw(M1, 3) * (three * scv - one + num_sqrt(two) * d * num_sqrt(d));
        const T ub = six * pw(M1, 3) * scv;
        return M3 >= lb * (one - tol) && M3 <= ub * (one + tol);
    }
    const T lb = three / two * pw(M1, 3) * (one + scv) * (one + scv);
    return M3 >= lb * (one - tol);
}

/**
 * The Telek-Heindl region as three inequality rows c <= 0 in (M2, M3) at a
 * fixed M1: the SCV floor of 1/2, and the third-moment bounds, whose upper
 * half is absent above SCV = 1.
 *
 * These rows are mathematically REDUNDANT with the inversion constraints of
 * opt_char -- they describe the same region -- but they are imposed alongside
 * them because the inversion form does not describe it as a CLOSED set, and an
 * optimizer will find the gap. Reproduction, for M1 = 1 and the target
 * (M2, M3) = (1.6, 2.0): the inversion constraints alone are minimized at
 * (M2a, M3a) = (1.99794, 1.99983) with objective 0.398, well below the 1.4735
 * of any point of the true region, by driving the branching probability to
 * p1 = -2.5e-9 while the second phase mean diverges (l2 = 645). Every
 * inversion constraint is then satisfied to 2.5e-9, which fmincon's default
 * ConstraintTolerance of 1e-6 would accept outright, yet the returned moment
 * pair admits no APH(2) at all: at p1 = 0 exactly the distribution is
 * exponential and M3 must be 6 M1^3, not 2. The bounds below close that gap.
 */
template <class T>
std::vector<T> aph2_moment_bounds_rows(const T& M1, const T& xM2, const T& xM3) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);
    const T half = num_traits<T>::from_rational(1, 2);

    const T M1sq = M1 * M1;
    const T scv = (xM2 - M1sq) / M1sq;
    std::vector<T> c(3, zero);
    c[0] = half - scv;
    if (scv <= one) {
        const T d = one - scv;
        const T lb = three * pw(M1, 3) * (three * scv - one + num_sqrt(two) * d * num_sqrt(d));
        const T ub = six * pw(M1, 3) * scv;
        c[1] = lb - xM3;
        c[2] = xM3 - ub;
    } else {
        const T lb = three / two * pw(M1, 3) * (one + scv) * (one + scv);
        c[1] = lb - xM3;
        c[2] = zero;
    }
    return c;
}

/**
 * The APH(2) inversion of aph2_adjust's nested aph2_fit1 / aph2_fit2:
 * given (M1, M2, M3), the phase means l1, l2 and the continuation probability
 * p1. `swapped` selects fit2, which exchanges the two roots.
 *
 * Returns the discriminant tmp0 as well, because the constraint set needs it.
 * The square root is taken of max(tmp0, 0): during the search the optimizer
 * does visit points with a negative discriminant, where the inversion is not
 * real; there the constraint -tmp0 <= 0 is what carries the information, and
 * clamping keeps l1, l2, p1 finite so the other three constraints stay
 * evaluable. The denominator 6 M2 - 12 M1^2 vanishes at the exponential point
 * M2 = 2 M1^2; its magnitude is clamped to degen so that the constraint
 * functions remain bounded there, which no reference does because no
 * reference evaluates them at all (see the header note).
 */
template <class T>
struct Aph2Inversion {
    T tmp0;
    T l1;
    T l2;
    T p1;
};

template <class T>
Aph2Inversion<T> aph2_invert(const T& M1, const T& xM2, const T& xM3, bool swapped,
                             const T& degen) {
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);
    const T eight = num_traits<T>::from_int(8);
    const T nine = num_traits<T>::from_int(9);
    const T twelve = num_traits<T>::from_int(12);

    Aph2Inversion<T> r;
    r.tmp0 = eight * pw(M1, 3) * xM3 / three - three * M1 * M1 * xM2 * xM2 -
             two * M1 * xM2 * xM3 + two * pw(xM2, 3) + xM3 * xM3 / nine;
    const T root = r.tmp0 > zero ? num_sqrt(r.tmp0) : zero;
    const T tmp1 = three * root;
    const T tmp2 = xM3 - three * M1 * xM2;
    T tmp3 = six * xM2 - twelve * M1 * M1;
    if (num_abs(tmp3) < degen) tmp3 = tmp3 < zero ? T(-degen) : degen;
    if (!swapped) {
        r.l1 = (tmp2 + tmp1) / tmp3;
        r.l2 = (tmp2 - tmp1) / tmp3;
    } else {
        r.l1 = (tmp2 - tmp1) / tmp3;
        r.l2 = (tmp2 + tmp1) / tmp3;
    }
    if (r.l2 == zero)
        r.p1 = zero;
    else
        r.p1 = (M1 - r.l1) / r.l2;
    return r;
}

}  // namespace fitdetail

/**
 * aph2_adjust, method 'opt_param': adjust (M2, M3) by searching the APH(2)
 * parameter space with M1 matched exactly.
 *
 * @param M1,M2,M3 the requested moments
 * @param feastol  the strict-inequality slack of the reference (MATLAB 1e-6),
 *                 used as the lower bound on l2 and in the constraint
 *                 l2 r1 <= M1 - feastol that keeps l1 positive
 * @param degentol the slack that keeps r1 away from 0 and 1 (MATLAB 1e-8);
 *                 set it to zero to allow the degenerate exponential
 */
template <class T>
Aph2AdjustOptResult<T> aph2_adjust_opt_param(const T& M1, const T& M2, const T& M3,
                                             const T& feastol, const T& degentol) {
    static_assert(num_traits<T>::has_transcendental,
                  "aph2_adjust_opt_param requires transcendental arithmetic");
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T six = num_traits<T>::from_int(6);
    const T half = num_traits<T>::from_rational(1, 2);
    if (M1 <= zero) throw InputError("aph2_adjust_opt_param: non-positive first moment");

    // moments of the APH(2) with phase means l1 = M1 - l2 r1 and l2, and
    // continuation probability r1 (aph2_adjust's fun2)
    auto moments = [M1, two, six](const std::vector<T>& x, T& xM2, T& xM3) {
        const T l2 = x[0];
        const T r1 = x[1];
        const T l1 = M1 - l2 * r1;
        xM2 = two * l1 * l1 + two * r1 * l1 * l2 + two * r1 * l2 * l2;
        xM3 = six * fitdetail::pw(l1, 3) + six * r1 * l1 * l1 * l2 + six * r1 * l1 * l2 * l2 +
              six * r1 * fitdetail::pw(l2, 3);
    };

    auto resid = [moments, M2, M3](const std::vector<T>& x) {
        T xM2 = num_traits<T>::from_int(0);
        T xM3 = num_traits<T>::from_int(0);
        moments(x, xM2, xM3);
        std::vector<T> r(2);
        r[0] = M2 - xM2;
        r[1] = M3 - xM3;
        return r;
    };

    // box bounds and the reference's nonlcon3, all as inequality rows
    auto g = [feastol, degentol, one, M1](const std::vector<T>& x) {
        std::vector<T> gv(4);
        gv[0] = feastol - x[0];              // l2 >= feastol
        gv[1] = degentol - x[1];             // r1 >= degentol
        gv[2] = x[1] - (one - degentol);     // r1 <= 1 - degentol
        gv[3] = x[0] * x[1] - M1 + feastol;  // l2 r1 <= M1 - feastol
        return gv;
    };

    std::vector<T> x0(2);
    x0[0] = M1;    // the reference's x0 = [M1, 1/2], which fits M1 exactly
    x0[1] = half;

    AugLagOptions<T> opt = auglag_defaults<T>();
    opt.ctol = num_traits<T>::from_double(1e-12);
    const AugLagResult<T> sol = auglag_ls(resid, std::size_t(2), NoConstraints<T>(), g, x0, opt);

    Aph2AdjustOptResult<T> r;
    moments(sol.x, r.M2a, r.M3a);
    const T d2 = M2 - r.M2a;
    const T d3 = M3 - r.M3a;
    r.objective = sqrt(T(d2 * d2 + d3 * d3));
    r.violation = sol.violation;
    r.converged = sol.violation <= opt.ctol;
    return r;
}

/** aph2_adjust_opt_param with the MATLAB defaults feastol = 1e-6, degentol = 1e-8. */
template <class T>
Aph2AdjustOptResult<T> aph2_adjust_opt_param(const T& M1, const T& M2, const T& M3) {
    return aph2_adjust_opt_param(M1, M2, M3, T(num_traits<T>::from_double(1e-6)),
                                 T(num_traits<T>::from_double(1e-8)));
}

/**
 * aph2_adjust, method 'opt_char': adjust (M2, M3) by searching the moment
 * space subject to APH(2) invertibility.
 *
 * Two constrained problems are solved, one per branch of the inversion, and
 * the one with the smaller adjustment is returned, exactly as the reference
 * intends.
 *
 * @param M1,M2,M3 the requested moments
 * @param postol   the strict-positivity slack applied to the phase means
 *                 (MATLAB 1e-6)
 * @param degen    magnitude below which the inversion denominator
 *                 6 M2 - 12 M1^2 is clamped; see the note on aph2_invert
 */
template <class T>
Aph2AdjustOptResult<T> aph2_adjust_opt_char(const T& M1, const T& M2, const T& M3, const T& postol,
                                            const T& degen) {
    static_assert(num_traits<T>::has_transcendental,
                  "aph2_adjust_opt_char requires transcendental arithmetic");
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    if (M1 <= zero) throw InputError("aph2_adjust_opt_char: non-positive first moment");

    auto resid = [M2, M3](const std::vector<T>& x) {
        std::vector<T> r(2);
        r[0] = x[0] - M2;
        r[1] = x[1] - M3;
        return r;
    };

    Aph2AdjustOptResult<T> best;
    bool have_best = false;

    for (int branch = 0; branch < 2; ++branch) {
        const bool swapped = branch == 1;
        auto g = [M1, swapped, postol, degen](const std::vector<T>& x) {
            const fitdetail::Aph2Inversion<T> inv =
                fitdetail::aph2_invert(M1, x[0], x[1], swapped, degen);
            std::vector<T> gv(6);
            gv[0] = -inv.tmp0;            // non-negative discriminant
            gv[1] = -inv.l1 + postol;     // positive first phase mean
            gv[2] = -inv.l2 + postol;     // positive second phase mean
            gv[3] = -inv.p1;              // non-negative branching probability
            gv[4] = -x[0];                // M2 >= 0, the reference's lb
            gv[5] = -x[1];                // M3 >= 0, the reference's lb
            const std::vector<T> th = fitdetail::aph2_moment_bounds_rows(M1, x[0], x[1]);
            for (std::size_t k = 0; k < th.size(); ++k) gv.push_back(th[k]);
            return gv;
        };

        std::vector<T> x0(2);
        x0[0] = M2;
        x0[1] = M3;
        AugLagOptions<T> opt = auglag_defaults<T>();
        opt.ctol = num_traits<T>::from_double(1e-12);
        const AugLagResult<T> sol =
            auglag_ls(resid, std::size_t(2), NoConstraints<T>(), g, x0, opt);

        Aph2AdjustOptResult<T> cand;
        cand.M2a = sol.x[0];
        cand.M3a = sol.x[1];
        const T d2 = cand.M2a - M2;
        const T d3 = cand.M3a - M3;
        cand.objective = sqrt(T(d2 * d2 + d3 * d3));
        cand.violation = sol.violation;
        cand.converged = sol.violation <= opt.ctol;

        // a branch that could not be made feasible is not a candidate at all
        if (!cand.converged && have_best) continue;
        if (!have_best || (cand.converged && !best.converged) ||
            (cand.converged == best.converged && cand.objective < best.objective)) {
            best = cand;
            have_best = true;
        }
    }
    return best;
}

/** aph2_adjust_opt_char with the MATLAB default postol = 1e-6 and degen = 1e-10. */
template <class T>
Aph2AdjustOptResult<T> aph2_adjust_opt_char(const T& M1, const T& M2, const T& M3) {
    return aph2_adjust_opt_char(M1, M2, M3, T(num_traits<T>::from_double(1e-6)),
                                T(num_traits<T>::from_double(1e-10)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH2_ADJUST_OPT_H
