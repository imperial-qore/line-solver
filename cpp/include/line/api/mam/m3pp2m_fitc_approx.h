/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_M3PP2M_FITC_APPROX_H
#define LINE_API_MAM_M3PP2M_FITC_APPROX_H

/**
 * M3PP(2, m) fitted to counting-process characteristics with an optimized
 * per-class split
 * (matlab/lib/m3a/m3a/m3pp/m3pp2m_fitc_approx.m,
 *  matlab/lib/m3a/m3a/m3pp/m3pp2m_fitc_approx_ag.m,
 *  matlab/lib/m3a/m3a/m3pp/m3pp2m_fitc_approx_ag_multiclass.m).
 *
 * The underlying MMPP(2) comes from mmpp2_fitc_approx. The per-class split is
 * then a pair of per-phase marking probabilities (q1i, q2i) per class, each an
 * AFFINE function of the class rate ai (matched exactly) and of one free
 * per-class quantity:
 *   - the non-'ag' variant frees dvi, the difference between the variance of
 *     class i and that of all other classes combined at resolution t3;
 *   - the 'ag' variant frees gi, the sum of the variance of class i and its
 *     covariance with all other classes combined at t3.
 * Since m3pp2m_fitc.h can match only m - 1 of those exactly (the last class
 * absorbs the remainder), the reference instead asks for the least-squares
 * compromise over ALL m classes, subject to the marking probabilities being
 * non-negative and summing to one in each phase. That is a quadratic program:
 *
 *   min sum_i (x_i / target_i - 1)^2
 *   s.t.  q1i(x_i) >= 0,  q2i(x_i) >= 0   for every class i
 *         sum_i q1i(x_i) = 1,  sum_i q2i(x_i) = 1
 *
 * (the reference writes the objective as 1/2 x'Hx + f'x with H = diag(2/t_i^2)
 * and f_i = -2/t_i, which is the expansion of the sum above minus the constant
 * m; its reported fit_error = fx + m is therefore exactly the sum above, and
 * that is what M3pp2mFitcApproxResult::class_objective holds).
 *
 * Because the two coefficients q1i_dvi and q2i_dvi do not depend on i, both
 * equality rows constrain only sum_i x_i; they have rank one and, as the
 * algebra of the reference guarantees, a consistent right-hand side, so the
 * program is a single equality plus 2m inequalities.
 *
 * ACCEPTANCE CONTRACT (see line/util/auglag.h). MATLAB solves the QP with
 * quadprog (interior-point-convex) or the bundled Goldfarb-Idnani QP; the JAR
 * uses OSQP. This port uses the augmented Lagrangian with a
 * Levenberg-Marquardt inner solver on the residual (x_i/target_i - 1). The
 * iterates do not match any of them. What is guaranteed and tested is the
 * specification:
 *   1. the returned MMAP is valid -- sum_i Dc_i = D1 and every Dc_i
 *      non-negative to the toolbox feasibility tolerance 1e-8, and the
 *      underlying (D0, D1) passes map_isfeasible (see mmap_isfeasible_tol);
 *   2. the per-class rates are matched exactly, by construction;
 *   3. the achieved objective (the reference's fit_error) is returned so it
 *      can be compared against any other QP solver's on the same input;
 *   4. when the requested per-class targets are those of an actual M3PP(2, m),
 *      the objective reaches zero and the achieved targets, recomputed
 *      independently with mmap_count_var, reproduce the requested ones.
 *
 * REFERENCE DEFECTS (matlab/lib/m3a/m3a/m3pp/m3pp2m_fitc_approx.m):
 *   1. the QP call reads `[x,fx] = QP(H, h, A, b, Aeq, beq, lb, ub, options)`,
 *      but `h` is never defined in that function -- the linear term is built
 *      as `f`. Every call therefore raises "Unrecognized function or variable
 *      'h'" before the QP runs, so this MATLAB entry point cannot execute.
 *      The '_ag' path is unaffected: m3pp2m_fitc_approx_ag_multiclass.m calls
 *      quadprog(H, f, ...) correctly.
 *   2. the same dead call passes lb = 1e-6 * ones(m,1). The free variables are
 *      variance DIFFERENCES, which are routinely negative (a minority class
 *      has less variance than all others combined), and the equality row fixes
 *      their sum, which is negative on ordinary inputs. Reproduction: for the
 *      MMPP(2) l1 = 2, l2 = 0.5, r1 = 0.3, r2 = 0.7 at t3 = 1 the equality
 *      forces sum_i x_i = -1.8976, so with m = 3 and lb = 1e-6 the feasible
 *      set is empty. This port therefore leaves the box unconstrained by
 *      default -- which is what the sibling _ag_multiclass does, its
 *      quadprog call passing [] for lb and ub -- and exposes optional bounds
 *      for a caller that wants them.
 * MATLAB was NOT edited.
 *
 * DIVERGENCE FROM THE JAR: M3pp2m_fitc_approx.java and
 * M3pp2m_fitc_approx_ag_multiclass.java add 2m further inequality rows
 * (q1i <= 1, q2i <= 1) that neither MATLAB file has, and then project the
 * solution post hoc with max(0, q) followed by a row renormalization. The
 * extra rows can make a program infeasible that MATLAB solves, and the
 * projection silently changes the fitted per-class rates. This port follows
 * MATLAB: 2m rows, no post-hoc projection, and the infeasibility of the result
 * reported through M3pp2mFitcApproxResult::feasible rather than repaired.
 *
 * Gated on transcendental arithmetic.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmpp2_fitc_approx.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of the optimization-based M3PP(2, m) fits. */
template <class T>
struct M3pp2mFitcApproxResult {
    Mmap<T> mmap;         ///< the fitted M3PP(2, m)
    T mmpp_objective;     ///< objective of the underlying MMPP(2) fit (zero when it was given)
    T class_objective;    ///< sum_i (x_i/target_i - 1)^2, the reference's fit_error
    T class_violation;    ///< worst constraint violation of the per-class QP
    bool degenerate;      ///< the underlying process was a Poisson process
    bool feasible;        ///< mmap_isfeasible_tol of the result at 1e-8
};

namespace fitdetail {

/** The six affine coefficients of the per-class marking probabilities. */
template <class T>
struct M3ppSplitCoeffs {
    T q1_a;  ///< coefficient of ai in q1i
    T q1_x;  ///< coefficient of the free variable in q1i
    T q1_c;  ///< constant term of q1i
    T q2_a;
    T q2_x;
    T q2_c;
};

/**
 * Coefficients of the non-'ag' split (m3pp2m_fitc_approx.m), the free variable
 * being the per-class variance difference at t3. The reference's repeated
 * product sinh(u) exp(-u) with u = (r1 + r2) t / 2 is (1 - exp(-(r1+r2) t))/2
 * exactly; it is evaluated in that form here, as in m3pp2m_fitc.h.
 */
template <class T>
M3ppSplitCoeffs<T> m3pp_split_coeffs_dv(const T& l1, const T& l2, const T& r1, const T& r2,
                                        const T& t) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T four = num_traits<T>::from_int(4);
    const T SH = (one - num_exp(T(-(r1 + r2) * t))) / two;

    const T BR = two * l1 * SH - two * l2 * SH - l1 * r1 * t - l1 * r2 * t + l2 * r1 * t +
                 l2 * r2 * t;
    const T DEN1 = l1 * r2 * (r1 + r2) * BR;
    const T BR2 = l2 * l2 * r1 * r1 * t - two * l2 * l2 * r1 * SH - l1 * l2 * r1 * r1 * t +
                  l2 * l2 * r1 * r2 * t + two * l1 * l2 * r1 * SH - l1 * l2 * r1 * r2 * t;
    const T DEN2 = (r1 + r2) * BR2;

    const T P4 = pw(r1, 4) + four * pw(r1, 3) * r2 + num_traits<T>::from_int(6) * r1 * r1 * r2 * r2 +
                 four * r1 * pw(r2, 3) + pw(r2, 4);
    const T CN = l1 * pw(r2, 4) * t + l2 * pw(r1, 4) * t + three * l1 * r1 * pw(r2, 3) * t +
                 l1 * pw(r1, 3) * r2 * t + l2 * r1 * pw(r2, 3) * t + three * l2 * pw(r1, 3) * r2 * t +
                 three * l1 * r1 * r1 * r2 * r2 * t + two * l1 * l1 * r1 * r2 * r2 * t +
                 two * l1 * l1 * r1 * r1 * r2 * t + three * l2 * r1 * r1 * r2 * r2 * t +
                 two * l2 * l2 * r1 * r2 * r2 * t + two * l2 * l2 * r1 * r1 * r2 * t -
                 four * l1 * l1 * r1 * r2 * SH - four * l2 * l2 * r1 * r2 * SH -
                 four * l1 * l2 * r1 * r2 * r2 * t - four * l1 * l2 * r1 * r1 * r2 * t +
                 num_traits<T>::from_int(8) * l1 * l2 * r1 * r2 * SH;

    M3ppSplitCoeffs<T> q;
    q.q1_a = (pw(r1, 4) * t / two + pw(r2, 4) * t / two - l1 * pw(r2, 3) * t +
              l2 * pw(r2, 3) * t + two * r1 * pw(r2, 3) * t + two * pw(r1, 3) * r2 * t +
              three * r1 * r1 * r2 * r2 * t + two * l1 * r2 * r2 * SH - two * l2 * r2 * r2 * SH -
              two * l1 * r1 * r2 * r2 * t - l1 * r1 * r1 * r2 * t + two * l2 * r1 * r2 * r2 * t +
              l2 * r1 * r1 * r2 * t + two * l1 * r1 * r2 * SH - two * l2 * r1 * r2 * SH) /
             DEN1;
    q.q1_x = -P4 / (four * DEN1);
    q.q1_c = -CN / (four * DEN1);
    q.q2_a = -(pw(r1, 4) * t / two + pw(r2, 4) * t / two + l1 * pw(r1, 3) * t -
               l2 * pw(r1, 3) * t + two * r1 * pw(r2, 3) * t + two * pw(r1, 3) * r2 * t +
               three * r1 * r1 * r2 * r2 * t - two * l1 * r1 * r1 * SH + two * l2 * r1 * r1 * SH +
               l1 * r1 * r2 * r2 * t + two * l1 * r1 * r1 * r2 * t - l2 * r1 * r2 * r2 * t -
               two * l2 * r1 * r1 * r2 * t - two * l1 * r1 * r2 * SH + two * l2 * r1 * r2 * SH) /
             DEN2;
    q.q2_x = P4 / (four * DEN2);
    q.q2_c = CN / (four * DEN2);
    return q;
}

/**
 * Coefficients of the 'ag' split (m3pp2m_fitc_approx_ag_multiclass.m), the
 * free variable being the per-class variance-plus-covariance at t3. Both
 * constant terms are zero there.
 */
template <class T>
M3ppSplitCoeffs<T> m3pp_split_coeffs_ag(const T& l1, const T& l2, const T& r1, const T& r2,
                                        const T& t) {
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T E = num_exp(T(-(r1 + r2) * t));

    const T f1num = l1 * r2 *
                    (two * l2 * r1 - two * l1 * r1 + pw(r1, 3) * t + pw(r2, 3) * t +
                     two * l1 * r1 * r1 * t - two * l2 * r1 * r1 * t + three * r1 * r2 * r2 * t +
                     three * r1 * r1 * r2 * t + two * l1 * r1 * E - two * l2 * r1 * E +
                     two * l1 * r1 * r2 * t - two * l2 * r1 * r2 * t);
    const T f1 = f1num / pw(T(r1 + r2), 4);
    const T f2num = l2 * r1 *
                    (two * l1 * r2 - two * l2 * r2 + pw(r1, 3) * t + pw(r2, 3) * t -
                     two * l1 * r2 * r2 * t + two * l2 * r2 * r2 * t + three * r1 * r2 * r2 * t +
                     three * r1 * r1 * r2 * t - two * l1 * r2 * E + two * l2 * r2 * E -
                     two * l1 * r1 * r2 * t + two * l2 * r1 * r2 * t);
    const T f2 = f2num / pw(T(r1 + r2), 4);
    const T tmp = f1 * l2 * r1 - f2 * l1 * r2;
    if (tmp == zero)
        throw NumericError("m3pp2m_fitc_approx_ag: degenerate per-class split (f1 l2 r1 = f2 l1 r2)");

    M3ppSplitCoeffs<T> q;
    q.q1_a = -(f2 * (r1 + r2)) / tmp;
    q.q1_x = (l2 * r1) / tmp;
    q.q1_c = zero;
    q.q2_a = (f1 * (r1 + r2)) / tmp;
    q.q2_x = -(l1 * r2) / tmp;
    q.q2_c = zero;
    return q;
}

/**
 * Solve the per-class quadratic program described in the header comment.
 *
 * @param q       the six affine coefficients
 * @param ai      per-class rates
 * @param target  per-class targets (dvt3 or gt3); none may be zero, since the
 *                objective is relative to them
 * @param a       total rate
 * @param bounds  optional box on the free variables; pass an empty vector for
 *                the unconstrained default (see the reference-defect note)
 * @param opt     tuning of the constrained solve
 * @param objective  out: sum_i (x_i/target_i - 1)^2 at the answer
 * @param violation  out: worst constraint violation at the answer
 */
template <class T>
std::vector<T> m3pp_split_solve(const M3ppSplitCoeffs<T>& q, const std::vector<T>& ai,
                                const std::vector<T>& target, const T& a,
                                const std::vector<Bound<T>>& bounds, const AugLagOptions<T>& opt,
                                T& objective, T& violation) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t m = ai.size();
    if (target.size() != m) throw InputError("m3pp_split_solve: target length disagrees with ai");
    for (std::size_t i = 0; i < m; ++i)
        if (target[i] == zero)
            throw InputError("m3pp_split_solve: a per-class target is zero, so the relative "
                             "objective of the reference is undefined");
    if (!bounds.empty() && bounds.size() != m)
        throw InputError("m3pp_split_solve: one bound per class is required");

    const std::vector<T> tgt = target;
    auto resid = [tgt, m, one](const std::vector<T>& x) {
        std::vector<T> r(m);
        for (std::size_t i = 0; i < m; ++i) r[i] = x[i] / tgt[i] - one;
        return r;
    };

    const std::vector<T> rates = ai;
    const M3ppSplitCoeffs<T> qc = q;
    const T mT = num_traits<T>::from_int(long(m));
    auto heq = [qc, rates, a, m, mT, one](const std::vector<T>& x) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < m; ++i) s += x[i];
        std::vector<T> hv(2);
        hv[0] = qc.q1_a * a + qc.q1_x * s + mT * qc.q1_c - one;
        hv[1] = qc.q2_a * a + qc.q2_x * s + mT * qc.q2_c - one;
        return hv;
    };
    const std::vector<Bound<T>> bnd = bounds;
    auto gineq = [qc, rates, m, bnd](const std::vector<T>& x) {
        std::vector<T> gv;
        for (std::size_t i = 0; i < m; ++i) {
            // q1i >= 0 and q2i >= 0, the reference's A x <= b
            gv.push_back(T(-(qc.q1_a * rates[i] + qc.q1_x * x[i] + qc.q1_c)));
            gv.push_back(T(-(qc.q2_a * rates[i] + qc.q2_x * x[i] + qc.q2_c)));
        }
        for (std::size_t i = 0; i < bnd.size(); ++i) {
            if (bnd[i].has_lo) gv.push_back(T(bnd[i].lo - x[i]));
            if (bnd[i].has_hi) gv.push_back(T(x[i] - bnd[i].hi));
        }
        return gv;
    };

    // start from the requested targets themselves, the point of zero objective
    std::vector<T> x0 = target;
    const AugLagResult<T> sol = auglag_ls(resid, m, heq, gineq, x0, opt);
    objective = sol.fval;
    violation = sol.violation;
    return sol.x;
}

}  // namespace fitdetail

/**
 * Assemble the M3PP from an underlying MAP and the per-phase marking
 * probabilities. The last class is NOT special-cased: every class gets the
 * probability the QP assigned it, and the equality rows are what make the
 * columns sum to D1.
 */
template <class T>
Mmap<T> m3pp2m_assemble(const Map<T>& base, const std::vector<T>& q1, const std::vector<T>& q2) {
    const T zero = num_traits<T>::from_int(0);
    if (q1.size() != q2.size()) throw InputError("m3pp2m_assemble: q1 and q2 differ in length");
    Mmap<T> mm;
    mm.D0 = base.D0;
    mm.D1 = base.D1;
    for (std::size_t i = 0; i < q1.size(); ++i) {
        Matrix<T> Dc(2, 2, zero);
        Dc(0, 0) = q1[i] * base.D1(0, 0);
        Dc(1, 1) = q2[i] * base.D1(1, 1);
        mm.Dc.push_back(Dc);
    }
    return mm;
}

namespace fitdetail {

/** Shared body of the two multiclass splits, once the coefficients are known. */
template <class T>
M3pp2mFitcApproxResult<T> m3pp2m_split_and_assemble(const Map<T>& base,
                                                    const M3ppSplitCoeffs<T>& q,
                                                    const std::vector<T>& ai,
                                                    const std::vector<T>& target, const T& a,
                                                    const std::vector<Bound<T>>& bounds,
                                                    const AugLagOptions<T>& opt) {
    const std::size_t m = ai.size();
    M3pp2mFitcApproxResult<T> res;
    res.degenerate = false;
    res.mmpp_objective = num_traits<T>::from_int(0);

    T obj = num_traits<T>::from_int(0);
    T viol = num_traits<T>::from_int(0);
    const std::vector<T> x = m3pp_split_solve(q, ai, target, a, bounds, opt, obj, viol);

    std::vector<T> q1(m), q2(m);
    for (std::size_t i = 0; i < m; ++i) {
        q1[i] = q.q1_a * ai[i] + q.q1_x * x[i] + q.q1_c;
        q2[i] = q.q2_a * ai[i] + q.q2_x * x[i] + q.q2_c;
    }
    res.mmap = m3pp2m_assemble(base, q1, q2);
    res.class_objective = obj;
    res.class_violation = viol;
    res.feasible = mmap_isfeasible_tol(res.mmap, T(num_traits<T>::from_double(1e-8)));
    return res;
}

/** The Poisson and single-class short-circuits shared by all three entry points. */
template <class T>
bool m3pp2m_trivial_split(const Map<T>& base, const std::vector<T>& ai, const T& a,
                          M3pp2mFitcApproxResult<T>& res) {
    const std::size_t m = ai.size();
    if (base.D0.rows() == 1) {
        // marked Poisson process: split D1 in proportion to the class rates
        res.degenerate = true;
        res.mmap.D0 = base.D0;
        res.mmap.D1 = base.D1;
        for (std::size_t i = 0; i < m; ++i) {
            Matrix<T> Dc(1, 1, T(ai[i] / a * base.D1(0, 0)));
            res.mmap.Dc.push_back(Dc);
        }
        res.class_objective = num_traits<T>::from_int(0);
        res.class_violation = num_traits<T>::from_int(0);
        res.feasible = mmap_isfeasible_tol(res.mmap, T(num_traits<T>::from_double(1e-8)));
        return true;
    }
    if (m == 1) {
        res.mmap.D0 = base.D0;
        res.mmap.D1 = base.D1;
        res.mmap.Dc.push_back(base.D1);
        res.class_objective = num_traits<T>::from_int(0);
        res.class_violation = num_traits<T>::from_int(0);
        res.feasible = mmap_isfeasible_tol(res.mmap, T(num_traits<T>::from_double(1e-8)));
        return true;
    }
    return false;
}

}  // namespace fitdetail

/**
 * m3pp2m_fitc_approx: fit the underlying MMPP(2) by optimization, then split
 * the classes on their variance DIFFERENCES at t3.
 *
 * @param a,bt1,bt2,binf,m3t2,t1,t2 the aggregate counting characteristics
 * @param ai    per-class rates, which must sum to a
 * @param dvt3  per-class variance differences at t3
 * @param t3    the third time scale
 * @param bounds optional box on the free variables; empty for none, which is
 *               the default (see the reference-defect note in the header)
 * @param opt tuning of the constrained solve
 */
template <class T>
M3pp2mFitcApproxResult<T> m3pp2m_fitc_approx(const T& a, const T& bt1, const T& bt2, const T& binf,
                                             const T& m3t2, const T& t1, const T& t2,
                                             const std::vector<T>& ai, const std::vector<T>& dvt3,
                                             const T& t3, const std::vector<Bound<T>>& bounds,
                                             const AugLagOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp2m_fitc_approx requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t m = ai.size();
    if (m == 0) throw InputError("m3pp2m_fitc_approx: no classes");
    T asum = zero;
    for (std::size_t i = 0; i < m; ++i) asum += ai[i];
    if (num_abs(T(a - asum)) > num_traits<T>::from_double(1e-8))
        throw InputError("m3pp2m_fitc_approx: inconsistent per-class arrival rates");

    const Mmpp2FitcApproxResult<T> base = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, opt);

    M3pp2mFitcApproxResult<T> res;
    res.degenerate = false;
    res.mmpp_objective = base.objective;
    if (fitdetail::m3pp2m_trivial_split(base.map, ai, a, res)) return res;

    const fitdetail::M3ppSplitCoeffs<T> q = fitdetail::m3pp_split_coeffs_dv(
        base.map.D1(0, 0), base.map.D1(1, 1), base.map.D0(0, 1), base.map.D0(1, 0), t3);
    M3pp2mFitcApproxResult<T> out =
        fitdetail::m3pp2m_split_and_assemble(base.map, q, ai, dvt3, a, bounds, opt);
    out.mmpp_objective = base.objective;
    return out;
}

/** m3pp2m_fitc_approx with no box on the free variables and the default tuning. */
template <class T>
M3pp2mFitcApproxResult<T> m3pp2m_fitc_approx(const T& a, const T& bt1, const T& bt2, const T& binf,
                                             const T& m3t2, const T& t1, const T& t2,
                                             const std::vector<T>& ai, const std::vector<T>& dvt3,
                                             const T& t3) {
    AugLagOptions<T> opt = auglag_defaults<T>();
    opt.ctol = num_traits<T>::from_double(1e-12);
    return m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3,
                              std::vector<Bound<T>>(), opt);
}

/**
 * m3pp2m_fitc_approx_ag_multiclass: split a GIVEN MMPP(2) into m classes on
 * their variance-plus-covariance at t3.
 *
 * @param mmpp  the underlying MAP, of order 1 (Poisson) or 2
 * @param ai    per-class rates, which must sum to the rate of mmpp
 * @param gt3   per-class variance-plus-covariance targets at t3
 * @param t3    the third time scale
 * @param bounds optional box on the free variables; empty for none
 * @param opt tuning of the constrained solve
 */
template <class T>
M3pp2mFitcApproxResult<T> m3pp2m_fitc_approx_ag_multiclass(const Map<T>& mmpp,
                                                           const std::vector<T>& ai,
                                                           const std::vector<T>& gt3, const T& t3,
                                                           const std::vector<Bound<T>>& bounds,
                                                           const AugLagOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp2m_fitc_approx_ag_multiclass requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t m = ai.size();
    if (m == 0) throw InputError("m3pp2m_fitc_approx_ag_multiclass: no classes");
    const T a = map_lambda(mmpp);
    T asum = zero;
    for (std::size_t i = 0; i < m; ++i) asum += ai[i];
    if (num_abs(T(a - asum)) > num_traits<T>::from_double(1e-8))
        throw InputError("m3pp2m_fitc_approx_ag_multiclass: inconsistent per-class arrival rates");

    M3pp2mFitcApproxResult<T> res;
    res.degenerate = false;
    res.mmpp_objective = zero;
    if (fitdetail::m3pp2m_trivial_split(mmpp, ai, a, res)) return res;
    if (mmpp.D0.rows() != 2)
        throw InputError("m3pp2m_fitc_approx_ag_multiclass: the underlying MAP must have order 2");

    const fitdetail::M3ppSplitCoeffs<T> q = fitdetail::m3pp_split_coeffs_ag(
        mmpp.D1(0, 0), mmpp.D1(1, 1), mmpp.D0(0, 1), mmpp.D0(1, 0), t3);
    return fitdetail::m3pp2m_split_and_assemble(mmpp, q, ai, gt3, a, bounds, opt);
}

/** m3pp2m_fitc_approx_ag_multiclass with no box and the default tuning. */
template <class T>
M3pp2mFitcApproxResult<T> m3pp2m_fitc_approx_ag_multiclass(const Map<T>& mmpp,
                                                           const std::vector<T>& ai,
                                                           const std::vector<T>& gt3,
                                                           const T& t3) {
    AugLagOptions<T> opt = auglag_defaults<T>();
    opt.ctol = num_traits<T>::from_double(1e-12);
    return m3pp2m_fitc_approx_ag_multiclass(mmpp, ai, gt3, t3, std::vector<Bound<T>>(), opt);
}

/**
 * m3pp2m_fitc_approx_ag: fit the underlying MMPP(2) by optimization, then
 * apply the 'ag' per-class split.
 *
 * @param a,bt1,bt2,binf,m3t2,t1,t2 the aggregate counting characteristics
 * @param ai   per-class rates, which must sum to a
 * @param gt3  per-class variance-plus-covariance targets at t3
 * @param t3   the third time scale
 * @param bounds optional box on the free variables; empty for none
 * @param opt tuning of the constrained solve
 */
template <class T>
M3pp2mFitcApproxResult<T> m3pp2m_fitc_approx_ag(const T& a, const T& bt1, const T& bt2,
                                                const T& binf, const T& m3t2, const T& t1,
                                                const T& t2, const std::vector<T>& ai,
                                                const std::vector<T>& gt3, const T& t3,
                                                const std::vector<Bound<T>>& bounds,
                                                const AugLagOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp2m_fitc_approx_ag requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    T asum = zero;
    for (std::size_t i = 0; i < ai.size(); ++i) asum += ai[i];
    if (num_abs(T(a - asum)) > num_traits<T>::from_double(1e-8))
        throw InputError("m3pp2m_fitc_approx_ag: inconsistent per-class arrival rates");

    const Mmpp2FitcApproxResult<T> base = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, opt);
    M3pp2mFitcApproxResult<T> res =
        m3pp2m_fitc_approx_ag_multiclass(base.map, ai, gt3, t3, bounds, opt);
    res.mmpp_objective = base.objective;
    return res;
}

/** m3pp2m_fitc_approx_ag with no box and the default tuning. */
template <class T>
M3pp2mFitcApproxResult<T> m3pp2m_fitc_approx_ag(const T& a, const T& bt1, const T& bt2,
                                                const T& binf, const T& m3t2, const T& t1,
                                                const T& t2, const std::vector<T>& ai,
                                                const std::vector<T>& gt3, const T& t3) {
    AugLagOptions<T> opt = auglag_defaults<T>();
    opt.ctol = num_traits<T>::from_double(1e-12);
    return m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3,
                                 std::vector<Bound<T>>(), opt);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_M3PP2M_FITC_APPROX_H
