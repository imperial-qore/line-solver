/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_NELDERMEAD_H
#define LINE_UTIL_NELDERMEAD_H

/**
 * Derivative-free simplex minimization (Nelder and Mead, 1965), with optional
 * box bounds imposed by a change of variables.
 *
 * This is the fallback for objectives that are not a sum of squares, that are
 * only piecewise smooth (the max(.) in an augmented Lagrangian's inequality
 * term, the max/min clamping in the AMAP(2) autocorrelation bounds), or whose
 * derivative is not worth n extra evaluations. Where the objective *is* a sum
 * of squares, prefer line/util/levmar.h: it converges in far fewer
 * evaluations and reports the residual vector.
 *
 * ACCEPTANCE CONTRACT. As for levmar.h: this replaces MATLAB's fminsearch /
 * fminsearchbnd / patternsearch / PSwarm, and the iterates do not and cannot
 * agree with theirs. Correctness is judged on the returned objective value and
 * on the specification the caller states, never on iterate-by-iterate
 * agreement with a reference implementation.
 *
 * DETERMINISM. The initial simplex is constructed from x0 alone: vertex 0 is
 * x0 and vertex j+1 perturbs coordinate j by step_rel * |x0_j|, or by step_abs
 * when x0_j is zero. This is the fminsearch construction and involves no
 * random numbers, so a given (objective, x0, options) triple always produces
 * the same answer, on any platform. Ties in the vertex ordering are broken by
 * vertex index, which keeps the sort stable and the run reproducible.
 *
 * BOX BOUNDS. Bounds are enforced by an unconstrained reparameterization of
 * each variable, so the simplex itself is unconstrained and every evaluated
 * point is strictly feasible:
 *   lower and upper: x = lo + (hi - lo) (sin(u) + 1)/2
 *   lower only:      x = lo + u^2
 *   upper only:      x = hi - u^2
 * The transformation is not a bijection (it folds), which is harmless for
 * minimization but means the returned x, not the internal u, is the answer.
 * A variable whose bounds coincide is held fixed at that value and removed
 * from the search. Note the well-known cost of the technique: the objective
 * seen by the simplex is flat at an active bound (dx/du = 0 there), so
 * convergence *onto* a bound is slower than in the interior; the returned
 * point still satisfies the bound exactly.
 *
 * Gated on transcendental arithmetic: the method stops on tolerances, and the
 * two-sided bound transformation needs sin/asin.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {

/** Tuning of the simplex iteration. */
template <class T>
struct NelderMeadOptions {
    T ftol;              ///< stop when the spread of f over the simplex falls below this
    T xtol;              ///< stop when the diameter of the simplex falls below this
    T step_rel;          ///< initial perturbation, relative, for a nonzero coordinate
    T step_abs;          ///< initial perturbation, absolute, for a zero coordinate
    T alpha;             ///< reflection coefficient
    T gamma;             ///< expansion coefficient
    T rho;               ///< contraction coefficient
    T sigma;             ///< shrink coefficient
    unsigned max_iter;   ///< cap on simplex iterations
    unsigned max_eval;   ///< cap on objective evaluations
};

/** fminsearch's coefficients and initial simplex, with tighter tolerances. */
template <class T>
NelderMeadOptions<T> nelder_mead_defaults() {
    NelderMeadOptions<T> o;
    o.ftol = num_traits<T>::from_double(1e-12);
    o.xtol = num_traits<T>::from_double(1e-12);
    o.step_rel = num_traits<T>::from_double(0.05);
    o.step_abs = num_traits<T>::from_double(0.00025);
    o.alpha = num_traits<T>::from_int(1);
    o.gamma = num_traits<T>::from_int(2);
    o.rho = num_traits<T>::from_rational(1, 2);
    o.sigma = num_traits<T>::from_rational(1, 2);
    o.max_iter = 5000;
    o.max_eval = 20000;
    return o;
}

/** Outcome of a simplex minimization. */
template <class T>
struct NelderMeadResult {
    std::vector<T> x;      ///< best point found
    T fval;                ///< objective there
    unsigned iterations;   ///< simplex iterations performed
    unsigned evaluations;  ///< objective evaluations
    bool converged;        ///< both tolerances met before the caps
};

/**
 * Box constraint on one variable. Absent bounds are represented by the flags,
 * not by an infinite value, because not every supported number type has an
 * infinity.
 */
template <class T>
struct Bound {
    bool has_lo;
    bool has_hi;
    T lo;
    T hi;

    Bound() : has_lo(false), has_hi(false), lo(num_traits<T>::from_int(0)),
              hi(num_traits<T>::from_int(0)) {}
};

/** Unbounded variable. */
template <class T>
Bound<T> bound_free() {
    return Bound<T>();
}

/** lo <= x. */
template <class T>
Bound<T> bound_lower(const T& lo) {
    Bound<T> b;
    b.has_lo = true;
    b.lo = lo;
    return b;
}

/** x <= hi. */
template <class T>
Bound<T> bound_upper(const T& hi) {
    Bound<T> b;
    b.has_hi = true;
    b.hi = hi;
    return b;
}

/** lo <= x <= hi. */
template <class T>
Bound<T> bound_box(const T& lo, const T& hi) {
    Bound<T> b;
    b.has_lo = true;
    b.has_hi = true;
    b.lo = lo;
    b.hi = hi;
    return b;
}

namespace nmdetail {

template <class T>
inline T nm_sin(const T& v) {
    using std::sin;
    return sin(v);
}

template <class T>
inline T nm_asin(const T& v) {
    using std::asin;
    return asin(v);
}

template <class T>
inline T nm_sqrt(const T& v) {
    using std::sqrt;
    return sqrt(v);
}

/** u -> x for one variable. */
template <class T>
T untransform(const T& u, const Bound<T>& b) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (b.has_lo && b.has_hi) {
        if (b.hi <= b.lo) return b.lo;
        const T s = nm_sin(u);
        return b.lo + (b.hi - b.lo) * (s + one) / two;
    }
    if (b.has_lo) return b.lo + u * u;
    if (b.has_hi) return b.hi - u * u;
    return u;
}

/** x -> u for one variable, clamping x into the box first. */
template <class T>
T transform(const T& x, const Bound<T>& b) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T zero = num_traits<T>::from_int(0);
    if (b.has_lo && b.has_hi) {
        if (b.hi <= b.lo) return zero;
        T z = two * (x - b.lo) / (b.hi - b.lo) - one;
        if (z < -one) z = -one;
        if (z > one) z = one;
        return nm_asin(z);
    }
    if (b.has_lo) {
        const T d = x - b.lo;
        return d > zero ? nm_sqrt(d) : zero;
    }
    if (b.has_hi) {
        const T d = b.hi - x;
        return d > zero ? nm_sqrt(d) : zero;
    }
    return x;
}

}  // namespace nmdetail

/**
 * Unconstrained simplex minimization.
 *
 * @param f   objective, x -> T
 * @param x0  starting point, which becomes vertex 0 of the simplex
 * @param opt tuning
 */
template <class T, class F>
NelderMeadResult<T> nelder_mead(F f, const std::vector<T>& x0, const NelderMeadOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "nelder_mead requires transcendental arithmetic (it stops on tolerances)");
    const std::size_t n = x0.size();
    if (n == 0) throw InputError("nelder_mead: no variables");
    const T zero = num_traits<T>::from_int(0);

    std::vector<std::vector<T>> v(n + 1, x0);
    for (std::size_t j = 0; j < n; ++j) {
        if (x0[j] == zero)
            v[j + 1][j] = opt.step_abs;
        else
            v[j + 1][j] = x0[j] * (num_traits<T>::from_int(1) + opt.step_rel);
    }

    std::vector<T> fv(n + 1);
    unsigned evals = 0;
    for (std::size_t i = 0; i <= n; ++i) {
        fv[i] = f(v[i]);
        ++evals;
    }

    std::vector<std::size_t> ord(n + 1);
    NelderMeadResult<T> res;
    res.converged = false;
    res.iterations = 0;

    for (unsigned it = 0; it < opt.max_iter; ++it) {
        res.iterations = it + 1;
        for (std::size_t i = 0; i <= n; ++i) ord[i] = i;
        std::stable_sort(ord.begin(), ord.end(),
                         [&fv](std::size_t a, std::size_t b) { return fv[a] < fv[b]; });

        const std::size_t best = ord[0];
        const std::size_t worst = ord[n];
        const std::size_t second = ord[n - 1];

        T fspread = zero;
        T xspread = zero;
        for (std::size_t i = 0; i <= n; ++i) {
            const T df = num_abs(T(fv[i] - fv[best]));
            if (df > fspread) fspread = df;
            for (std::size_t j = 0; j < n; ++j) {
                const T dx = num_abs(T(v[i][j] - v[best][j]));
                if (dx > xspread) xspread = dx;
            }
        }
        if (fspread <= opt.ftol && xspread <= opt.xtol) {
            res.converged = true;
            break;
        }
        if (evals >= opt.max_eval) break;

        // centroid of everything but the worst vertex
        std::vector<T> c(n, zero);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) c[j] += v[ord[i]][j];
        for (std::size_t j = 0; j < n; ++j) c[j] /= num_traits<T>::from_int(long(n));

        std::vector<T> xr(n);
        for (std::size_t j = 0; j < n; ++j) xr[j] = c[j] + opt.alpha * (c[j] - v[worst][j]);
        const T fr = f(xr);
        ++evals;

        if (fr < fv[best]) {
            std::vector<T> xe(n);
            for (std::size_t j = 0; j < n; ++j) xe[j] = c[j] + opt.gamma * (xr[j] - c[j]);
            const T fe = f(xe);
            ++evals;
            if (fe < fr) {
                v[worst] = xe;
                fv[worst] = fe;
            } else {
                v[worst] = xr;
                fv[worst] = fr;
            }
            continue;
        }
        if (fr < fv[second]) {
            v[worst] = xr;
            fv[worst] = fr;
            continue;
        }

        // contraction, outside when the reflection improved on the worst
        bool shrink = false;
        if (fr < fv[worst]) {
            std::vector<T> xc(n);
            for (std::size_t j = 0; j < n; ++j) xc[j] = c[j] + opt.rho * (xr[j] - c[j]);
            const T fc = f(xc);
            ++evals;
            if (fc <= fr) {
                v[worst] = xc;
                fv[worst] = fc;
            } else {
                shrink = true;
            }
        } else {
            std::vector<T> xc(n);
            for (std::size_t j = 0; j < n; ++j) xc[j] = c[j] + opt.rho * (v[worst][j] - c[j]);
            const T fc = f(xc);
            ++evals;
            if (fc < fv[worst]) {
                v[worst] = xc;
                fv[worst] = fc;
            } else {
                shrink = true;
            }
        }

        if (shrink) {
            for (std::size_t i = 0; i <= n; ++i) {
                if (i == best) continue;
                for (std::size_t j = 0; j < n; ++j)
                    v[i][j] = v[best][j] + opt.sigma * (v[i][j] - v[best][j]);
                fv[i] = f(v[i]);
                ++evals;
            }
        }
    }

    std::size_t best = 0;
    for (std::size_t i = 1; i <= n; ++i)
        if (fv[i] < fv[best]) best = i;
    res.x = v[best];
    res.fval = fv[best];
    res.evaluations = evals;
    return res;
}

/** nelder_mead with the default tuning. */
template <class T, class F>
NelderMeadResult<T> nelder_mead(F f, const std::vector<T>& x0) {
    return nelder_mead(f, x0, nelder_mead_defaults<T>());
}

/**
 * Box-constrained simplex minimization by the transformation described in the
 * header comment. Every point at which f is evaluated satisfies the bounds.
 *
 * @param f      objective, x -> T, called only at feasible x
 * @param x0     starting point, clamped into the box if it is outside
 * @param bounds one Bound per variable
 * @param opt    tuning
 */
template <class T, class F>
NelderMeadResult<T> nelder_mead_box(F f, const std::vector<T>& x0,
                                    const std::vector<Bound<T>>& bounds,
                                    const NelderMeadOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "nelder_mead_box requires transcendental arithmetic");
    const std::size_t n = x0.size();
    if (bounds.size() != n) throw InputError("nelder_mead_box: one bound per variable is required");
    for (std::size_t j = 0; j < n; ++j)
        if (bounds[j].has_lo && bounds[j].has_hi && bounds[j].hi < bounds[j].lo)
            throw InputError("nelder_mead_box: upper bound below lower bound");

    // variables with coincident bounds are fixed and taken out of the search
    std::vector<std::size_t> free_idx;
    std::vector<T> xfix = x0;
    for (std::size_t j = 0; j < n; ++j) {
        if (bounds[j].has_lo && bounds[j].has_hi && bounds[j].hi == bounds[j].lo)
            xfix[j] = bounds[j].lo;
        else
            free_idx.push_back(j);
    }
    if (free_idx.empty()) {
        NelderMeadResult<T> r;
        r.x = xfix;
        r.fval = f(xfix);
        r.iterations = 0;
        r.evaluations = 1;
        r.converged = true;
        return r;
    }

    std::vector<T> u0(free_idx.size());
    for (std::size_t k = 0; k < free_idx.size(); ++k)
        u0[k] = nmdetail::transform(x0[free_idx[k]], bounds[free_idx[k]]);

    const std::vector<std::size_t> idx = free_idx;
    const std::vector<Bound<T>> bnd = bounds;
    std::vector<T> xbuf = xfix;
    auto expand = [idx, bnd, xbuf](const std::vector<T>& u) {
        std::vector<T> x = xbuf;
        for (std::size_t k = 0; k < idx.size(); ++k)
            x[idx[k]] = nmdetail::untransform(u[k], bnd[idx[k]]);
        return x;
    };

    const NelderMeadResult<T> ru =
        nelder_mead([f, expand](const std::vector<T>& u) { return f(expand(u)); }, u0, opt);

    NelderMeadResult<T> r;
    r.x = expand(ru.x);
    r.fval = ru.fval;
    r.iterations = ru.iterations;
    r.evaluations = ru.evaluations;
    r.converged = ru.converged;
    return r;
}

/** nelder_mead_box with the default tuning. */
template <class T, class F>
NelderMeadResult<T> nelder_mead_box(F f, const std::vector<T>& x0,
                                    const std::vector<Bound<T>>& bounds) {
    return nelder_mead_box(f, x0, bounds, nelder_mead_defaults<T>());
}

}  // namespace line

#endif  // LINE_UTIL_NELDERMEAD_H
