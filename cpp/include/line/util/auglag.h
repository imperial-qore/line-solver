/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_AUGLAG_H
#define LINE_UTIL_AUGLAG_H

/**
 * Augmented Lagrangian method for equality- and inequality-constrained
 * minimization, with line/util/neldermead.h or line/util/levmar.h as the inner
 * unconstrained solver.
 *
 * Solves
 *     min f(x)   s.t.   h_i(x) = 0,   g_j(x) <= 0,   lo <= x <= hi
 * by minimizing, for a sequence of penalty parameters rho and multiplier
 * estimates (lambda, mu),
 *
 *   L_A(x; lambda, mu, rho) = f(x)
 *                           + sum_i [ lambda_i h_i + (rho/2) h_i^2 ]
 *                           + (1/(2 rho)) sum_j [ max(0, mu_j + rho g_j)^2 - mu_j^2 ]
 *
 * which is the Hestenes-Powell-Rockafellar form: the inequality term is the
 * exact penalty of Rockafellar (1973), differentiable once, and is inactive
 * for a constraint that is strictly satisfied with a zero multiplier. After
 * each inner solve
 *     lambda_i <- lambda_i + rho h_i,      mu_j <- max(0, mu_j + rho g_j)
 * and rho is multiplied by rho_factor whenever the constraint violation did
 * not shrink by at least the factor `shrink`. This is a first-order multiplier
 * method: it converges to a KKT point without driving rho to infinity, which
 * is what keeps the inner problems well conditioned.
 *
 * WHY NOT A PENALTY-ONLY LOOP: with lambda held at zero, the minimizer of the
 * penalized problem is offset from the true solution by O(1/rho) and the only
 * way to tighten it is a large rho, whose Hessian is ill conditioned by
 * exactly that factor. The multiplier update removes the offset, so a moderate
 * rho suffices.
 *
 * ACCEPTANCE CONTRACT. This replaces MATLAB's fmincon (active-set /
 * interior-point), quadprog, patternsearch and PSwarm in the m3a and
 * kpctoolbox fitting routines. It is a different algorithm: it does not
 * reproduce their iterates, their multipliers, or, on a nonconvex problem
 * with several local minima, necessarily their local minimum. Every caller in
 * line/api/mam therefore states its acceptance in terms of the specification
 * (the fitted process reproduces the target characteristics to a stated
 * tolerance, and is a valid MAP) and returns the objective value it achieved,
 * so it can be compared against the reference's.
 *
 * Deterministic: no random multi-start, no global state, no exit(), no output.
 * Non-convergence is reported through AugLagResult, never thrown.
 *
 * Gated on transcendental arithmetic through the inner solvers.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/levmar.h"
#include "line/util/neldermead.h"

namespace line {

/** Tuning of the outer multiplier iteration. */
template <class T>
struct AugLagOptions {
    T rho0;                          ///< initial penalty parameter
    T rho_factor;                    ///< growth factor applied to rho when needed
    T rho_max;                       ///< cap on rho
    T ctol;                          ///< constraint violation accepted as feasible
    T shrink;                        ///< required violation reduction to leave rho alone
    unsigned max_outer;              ///< cap on outer iterations
    NelderMeadOptions<T> inner_nm;   ///< tuning of the simplex inner solver
    LevmarOptions<T> inner_lm;       ///< tuning of the least-squares inner solver
};

/** Defaults: rho0 = 10, growth 10, feasibility 1e-10, 50 outer iterations. */
template <class T>
AugLagOptions<T> auglag_defaults() {
    AugLagOptions<T> o;
    o.rho0 = num_traits<T>::from_int(10);
    o.rho_factor = num_traits<T>::from_int(10);
    o.rho_max = num_traits<T>::from_double(1e12);
    o.ctol = num_traits<T>::from_double(1e-10);
    o.shrink = num_traits<T>::from_double(0.25);
    o.max_outer = 50;
    o.inner_nm = nelder_mead_defaults<T>();
    o.inner_lm = levmar_defaults<T>();
    return o;
}

/** Outcome of a constrained solve. */
template <class T>
struct AugLagResult {
    std::vector<T> x;         ///< best point found
    T fval;                   ///< the ORIGINAL objective f(x), not the augmented one
    T violation;              ///< max(|h_i|, max(0, g_j)) at x
    std::vector<T> lambda;    ///< final equality multipliers
    std::vector<T> mu;        ///< final inequality multipliers, all >= 0
    unsigned outer_iterations;
    bool converged;           ///< feasible to ctol and the last inner solve converged
};

/** A constraint map that returns no constraints; the default for h or g. */
template <class T>
struct NoConstraints {
    std::vector<T> operator()(const std::vector<T>&) const { return std::vector<T>(); }
};

namespace agldetail {

/** max(|h_i|, max(0, g_j)). */
template <class T>
T violation_of(const std::vector<T>& h, const std::vector<T>& g) {
    const T zero = num_traits<T>::from_int(0);
    T v = zero;
    for (std::size_t i = 0; i < h.size(); ++i) {
        const T a = num_abs(h[i]);
        if (a > v) v = a;
    }
    for (std::size_t j = 0; j < g.size(); ++j)
        if (g[j] > v) v = g[j];
    return v;
}

}  // namespace agldetail

/**
 * Augmented Lagrangian with a scalar objective and a simplex inner solver.
 *
 * @param f      objective, x -> T
 * @param h      equality constraints, x -> vector (empty for none)
 * @param g      inequality constraints g(x) <= 0, x -> vector (empty for none)
 * @param x0     starting point
 * @param bounds one Bound per variable; pass all-free bounds for none
 * @param opt    tuning
 */
template <class T, class F, class H, class G>
AugLagResult<T> auglag(F f, H h, G g, const std::vector<T>& x0,
                       const std::vector<Bound<T>>& bounds, const AugLagOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "auglag requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    const std::size_t n = x0.size();
    if (n == 0) throw InputError("auglag: no variables");
    if (bounds.size() != n) throw InputError("auglag: one bound per variable is required");

    const std::size_t ne = h(x0).size();
    const std::size_t ni = g(x0).size();

    AugLagResult<T> res;
    res.x = x0;
    res.lambda.assign(ne, zero);
    res.mu.assign(ni, zero);
    res.outer_iterations = 0;
    res.converged = false;

    T rho = opt.rho0;
    T prev_viol = num_traits<T>::from_double(-1.0);
    bool inner_ok = false;

    for (unsigned outer = 0; outer < opt.max_outer; ++outer) {
        res.outer_iterations = outer + 1;
        const std::vector<T> lam = res.lambda;
        const std::vector<T> mu = res.mu;
        const T rho_c = rho;

        auto L = [f, h, g, lam, mu, rho_c, zero, two](const std::vector<T>& x) {
            T val = f(x);
            const std::vector<T> hv = h(x);
            for (std::size_t i = 0; i < hv.size(); ++i)
                val += lam[i] * hv[i] + rho_c / two * hv[i] * hv[i];
            const std::vector<T> gv = g(x);
            for (std::size_t j = 0; j < gv.size(); ++j) {
                const T s = mu[j] + rho_c * gv[j];
                if (s > zero) val += (s * s - mu[j] * mu[j]) / (two * rho_c);
                else val -= mu[j] * mu[j] / (two * rho_c);
            }
            return val;
        };

        const NelderMeadResult<T> in = nelder_mead_box(L, res.x, bounds, opt.inner_nm);
        res.x = in.x;
        inner_ok = in.converged;

        const std::vector<T> hv = h(res.x);
        const std::vector<T> gv = g(res.x);
        const T viol = agldetail::violation_of(hv, gv);
        res.violation = viol;

        for (std::size_t i = 0; i < ne; ++i) res.lambda[i] = res.lambda[i] + rho * hv[i];
        for (std::size_t j = 0; j < ni; ++j) {
            const T s = res.mu[j] + rho * gv[j];
            res.mu[j] = s > zero ? s : zero;
        }

        if (viol <= opt.ctol) {
            res.converged = inner_ok;
            break;
        }
        if (prev_viol >= zero && viol > opt.shrink * prev_viol && rho < opt.rho_max)
            rho *= opt.rho_factor;
        prev_viol = viol;
    }

    res.fval = f(res.x);
    return res;
}

/** auglag with the default tuning. */
template <class T, class F, class H, class G>
AugLagResult<T> auglag(F f, H h, G g, const std::vector<T>& x0,
                       const std::vector<Bound<T>>& bounds) {
    return auglag(f, h, g, x0, bounds, auglag_defaults<T>());
}

/**
 * Augmented Lagrangian with a least-squares objective and levmar as the inner
 * solver.
 *
 * The augmented Lagrangian of a sum of squares is itself a sum of squares up
 * to an additive constant, because
 *     lambda h + (rho/2) h^2       = (rho/2)(h + lambda/rho)^2 - lambda^2/(2 rho)
 *     (1/(2 rho)) max(0, mu + rho g)^2 = (rho/2) max(0, g + mu/rho)^2
 * so the inner problem is handed to levmar with the extended residual
 *     [ r(x) ; sqrt(rho/2) (h + lambda/rho) ; sqrt(rho/2) max(0, g + mu/rho) ].
 * The dropped constants do not move the minimizer. The max(.) makes the
 * extended residual only piecewise smooth, which the finite-difference
 * Jacobian tolerates because an inequality is either active or inactive over
 * a whole differencing step except on a measure-zero set of iterates.
 *
 * Bounds are NOT supported here (levmar is unconstrained): express them as
 * inequality rows of g, which is what the callers in line/api/mam do.
 *
 * @param r   residual map, x -> vector of length m; the objective is sum r_i^2
 * @param m   number of residuals
 * @param h   equality constraints
 * @param g   inequality constraints g(x) <= 0
 * @param x0  starting point
 * @param opt tuning
 */
template <class T, class R, class H, class G>
AugLagResult<T> auglag_ls(R r, std::size_t m, H h, G g, const std::vector<T>& x0,
                          const AugLagOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "auglag_ls requires transcendental arithmetic");
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    const std::size_t n = x0.size();
    if (n == 0) throw InputError("auglag_ls: no variables");

    const std::size_t ne = h(x0).size();
    const std::size_t ni = g(x0).size();

    AugLagResult<T> res;
    res.x = x0;
    res.lambda.assign(ne, zero);
    res.mu.assign(ni, zero);
    res.outer_iterations = 0;
    res.converged = false;

    T rho = opt.rho0;
    T prev_viol = num_traits<T>::from_double(-1.0);
    bool inner_ok = false;

    for (unsigned outer = 0; outer < opt.max_outer; ++outer) {
        res.outer_iterations = outer + 1;
        const std::vector<T> lam = res.lambda;
        const std::vector<T> mu = res.mu;
        const T rho_c = rho;
        const T w = sqrt(T(rho_c / two));

        auto Raug = [r, h, g, lam, mu, rho_c, w, zero](const std::vector<T>& x) {
            std::vector<T> out = r(x);
            const std::vector<T> hv = h(x);
            for (std::size_t i = 0; i < hv.size(); ++i)
                out.push_back(T(w * (hv[i] + lam[i] / rho_c)));
            const std::vector<T> gv = g(x);
            for (std::size_t j = 0; j < gv.size(); ++j) {
                const T s = gv[j] + mu[j] / rho_c;
                out.push_back(s > zero ? T(w * s) : zero);
            }
            return out;
        };

        const LevmarResult<T> in = levmar(Raug, res.x, m + ne + ni, opt.inner_lm);
        res.x = in.x;
        inner_ok = in.converged;

        const std::vector<T> hv = h(res.x);
        const std::vector<T> gv = g(res.x);
        const T viol = agldetail::violation_of(hv, gv);
        res.violation = viol;

        for (std::size_t i = 0; i < ne; ++i) res.lambda[i] = res.lambda[i] + rho * hv[i];
        for (std::size_t j = 0; j < ni; ++j) {
            const T s = res.mu[j] + rho * gv[j];
            res.mu[j] = s > zero ? s : zero;
        }

        if (viol <= opt.ctol) {
            res.converged = inner_ok;
            break;
        }
        if (prev_viol >= zero && viol > opt.shrink * prev_viol && rho < opt.rho_max)
            rho *= opt.rho_factor;
        prev_viol = viol;
    }

    res.fval = levmardetail::sum_squares(r(res.x));
    return res;
}

/** auglag_ls with the default tuning. */
template <class T, class R, class H, class G>
AugLagResult<T> auglag_ls(R r, std::size_t m, H h, G g, const std::vector<T>& x0) {
    return auglag_ls(r, m, h, g, x0, auglag_defaults<T>());
}

}  // namespace line

#endif  // LINE_UTIL_AUGLAG_H
