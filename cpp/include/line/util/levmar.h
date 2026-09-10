/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_LEVMAR_H
#define LINE_UTIL_LEVMAR_H

/**
 * Levenberg-Marquardt for nonlinear least squares.
 *
 * Minimizes S(x) = sum_i r_i(x)^2 for a caller-supplied residual map
 * r : R^n -> R^m. This is the workhorse behind the moment-matching fits in
 * line/api/mam: every one of them states its target as a vector of relative
 * errors (moment_fitted/moment_target - 1) that would be zero at an exact
 * match, which is exactly the shape LM wants.
 *
 * ACCEPTANCE CONTRACT. Substituting an optimizer is not a transcription. The
 * caller must NOT expect the iterates, the iteration count, or the last digits
 * of the answer to agree with MATLAB's fmincon / fminsearch / optimproblem
 * solve, which are different algorithms with different termination rules and,
 * in the GlobalSearch cases, a random multi-start. What a caller may rely on
 * is stated per entry point in terms of the *specification*: the achieved
 * objective value, which is returned so it can be compared against any other
 * optimizer's on the same input.
 *
 * Implementation notes:
 *   - the step solves (J^T J + mu I) dx = -J^T r with mu = lambda max_j
 *     (J^T J)_jj, the uniform damping of Madsen, Nielsen and Tingleff, which
 *     is scaled to the problem yet insensitive to a column of J that is only
 *     differencing noise (see the note in levmar_jac on why Marquardt's
 *     per-column damping cannot be used with a numeric Jacobian);
 *   - lambda is multiplied by lambda_increase on a rejected step and by
 *     lambda_decrease on an accepted one;
 *   - the Jacobian defaults to central differences with a relative step, which
 *     costs 2n residual evaluations per iteration and is second-order
 *     accurate; an analytic Jacobian is supplied through levmar_jac.
 *
 * Deterministic: no random restarts, no global state, no exit(), no output.
 * Failure to converge is reported through LevmarResult::converged, never
 * thrown, since a partially converged fit is still usable and the caller is
 * the one that knows the tolerance it needs.
 *
 * Gated on transcendental arithmetic: the method stops on tolerances, so it
 * is meaningless at exact arithmetic (which would run to the iteration cap
 * carrying ever larger rationals).
 */

#include <cstddef>
#include <memory>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {

/** Tuning of the Levenberg-Marquardt iteration. */
template <class T>
struct LevmarOptions {
    T ftol;             ///< stop when the relative decrease of S falls below this
    T xtol;             ///< stop when the relative step length falls below this
    T gtol;             ///< stop when max|J^T r| falls below this
    T lambda0;          ///< initial damping
    T lambda_increase;  ///< factor applied to lambda after a rejected step
    T lambda_decrease;  ///< factor applied to lambda after an accepted step
    T lambda_max;       ///< give up on the iteration once lambda exceeds this
    T diff_step;        ///< relative step of the central-difference Jacobian
    unsigned max_iter;  ///< cap on accepted-or-rejected outer iterations
};

/** MINPACK-like defaults, with a central-difference step of eps^(1/3). */
template <class T>
LevmarOptions<T> levmar_defaults() {
    LevmarOptions<T> o;
    o.ftol = num_traits<T>::from_double(1e-12);
    o.xtol = num_traits<T>::from_double(1e-12);
    o.gtol = num_traits<T>::from_double(1e-14);
    o.lambda0 = num_traits<T>::from_double(1e-3);
    o.lambda_increase = num_traits<T>::from_int(10);
    o.lambda_decrease = num_traits<T>::from_double(0.1);
    o.lambda_max = num_traits<T>::from_double(1e16);
    o.diff_step = num_traits<T>::from_double(1e-6);
    o.max_iter = 500;
    return o;
}

/** Outcome of a least-squares solve. */
template <class T>
struct LevmarResult {
    std::vector<T> x;         ///< best point found
    std::vector<T> residual;  ///< r(x)
    T ssq;                    ///< sum of squares at x, the objective value
    unsigned iterations;      ///< outer iterations performed
    unsigned evaluations;     ///< residual evaluations, differencing included
    bool converged;           ///< a tolerance was met before the caps
};

namespace levmardetail {

template <class T>
T sum_squares(const std::vector<T>& r) {
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < r.size(); ++i) s += r[i] * r[i];
    return s;
}

template <class T>
T norm2(const std::vector<T>& v) {
    using std::sqrt;
    return sqrt(sum_squares(v));
}

}  // namespace levmardetail

/**
 * Central-difference Jacobian of r at x.
 *
 * The step for component j is diff_step * max(|x_j|, 1), so a variable of any
 * magnitude gets a meaningful perturbation and a variable at zero still gets
 * one.
 *
 * @param f         residual map, x -> vector of length m
 * @param x         evaluation point
 * @param m         number of residuals
 * @param diff_step relative differencing step
 * @return the m x n Jacobian
 */
template <class T, class F>
Matrix<T> levmar_jacobian_fd(F f, const std::vector<T>& x, std::size_t m, const T& diff_step) {
    static_assert(num_traits<T>::has_transcendental,
                  "levmar requires transcendental arithmetic (it stops on tolerances)");
    const std::size_t n = x.size();
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    Matrix<T> J(m, n, num_traits<T>::from_int(0));
    std::vector<T> xp = x;
    for (std::size_t j = 0; j < n; ++j) {
        const T ax = num_abs(x[j]);
        const T scale = ax > one ? ax : one;
        const T h = diff_step * scale;
        xp[j] = x[j] + h;
        const std::vector<T> rp = f(xp);
        xp[j] = x[j] - h;
        const std::vector<T> rm = f(xp);
        xp[j] = x[j];
        if (rp.size() != m || rm.size() != m)
            throw InputError("levmar_jacobian_fd: residual length changed between evaluations");
        const T den = two * h;
        for (std::size_t i = 0; i < m; ++i) J(i, j) = (rp[i] - rm[i]) / den;
    }
    return J;
}

/**
 * Levenberg-Marquardt with a caller-supplied Jacobian.
 *
 * @param f    residual map, x -> vector of length m
 * @param jac  Jacobian map, x -> m x n Matrix
 * @param x0   starting point
 * @param m    number of residuals
 * @param opt  tuning
 */
template <class T, class F, class J>
LevmarResult<T> levmar_jac(F f, J jac, const std::vector<T>& x0, std::size_t m,
                           const LevmarOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "levmar requires transcendental arithmetic (it stops on tolerances)");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = x0.size();
    if (n == 0) throw InputError("levmar: no variables");
    if (m == 0) throw InputError("levmar: no residuals");

    LevmarResult<T> res;
    res.x = x0;
    res.residual = f(res.x);
    if (res.residual.size() != m) throw InputError("levmar: residual length disagrees with m");
    res.ssq = levmardetail::sum_squares(res.residual);
    res.iterations = 0;
    res.evaluations = 1;
    res.converged = false;

    T lambda = opt.lambda0;
    // ftol/xtol two-consecutive-step convergence test: see _kb/14-cpp-multiprecision.md
    unsigned tol_hits = 0;

    for (unsigned it = 0; it < opt.max_iter; ++it) {
        res.iterations = it + 1;
        const Matrix<T> Jm = jac(res.x);
        if (Jm.rows() != m || Jm.cols() != n) throw InputError("levmar: Jacobian has wrong shape");

        // gradient of S/2 and the Gauss-Newton normal matrix
        std::vector<T> g(n, zero);
        for (std::size_t j = 0; j < n; ++j) {
            T s = zero;
            for (std::size_t i = 0; i < m; ++i) s += Jm(i, j) * res.residual[i];
            g[j] = s;
        }
        T gmax = zero;
        for (std::size_t j = 0; j < n; ++j) {
            const T a = num_abs(g[j]);
            if (a > gmax) gmax = a;
        }
        if (gmax <= opt.gtol) {
            res.converged = true;
            return res;
        }

        Matrix<T> A(n, n, zero);
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t k = j; k < n; ++k) {
                T s = zero;
                for (std::size_t i = 0; i < m; ++i) s += Jm(i, j) * Jm(i, k);
                A(j, k) = s;
                A(k, j) = s;
            }

        // uniform Levenberg-Marquardt damping: see _kb/14-cpp-multiprecision.md
        T dmax = zero;
        for (std::size_t j = 0; j < n; ++j)
            if (A(j, j) > dmax) dmax = A(j, j);
        if (dmax == zero) dmax = num_traits<T>::from_int(1);

        bool accepted = false;
        while (!accepted && lambda <= opt.lambda_max) {
            Matrix<T> Aug = A;
            for (std::size_t j = 0; j < n; ++j) Aug(j, j) = A(j, j) + lambda * dmax;
            std::vector<T> rhs(n);
            for (std::size_t j = 0; j < n; ++j) rhs[j] = -g[j];

            std::vector<T> dx;
            bool solved = true;
            try {
                dx = solve(Aug, rhs);
            } catch (const NumericError&) {
                solved = false;
            }
            if (!solved) {
                lambda *= opt.lambda_increase;
                continue;
            }

            std::vector<T> xn(n);
            for (std::size_t j = 0; j < n; ++j) xn[j] = res.x[j] + dx[j];
            const std::vector<T> rn = f(xn);
            ++res.evaluations;
            if (rn.size() != m) throw InputError("levmar: residual length changed between calls");
            const T sn = levmardetail::sum_squares(rn);

            if (sn < res.ssq) {
                const T dnorm = levmardetail::norm2(dx);
                const T xnorm = levmardetail::norm2(res.x);
                const T dssq = res.ssq - sn;
                const bool ftol_met = dssq <= opt.ftol * res.ssq;
                const bool xtol_met = dnorm <= opt.xtol * (xnorm + opt.xtol);
                res.x = xn;
                res.residual = rn;
                res.ssq = sn;
                lambda *= opt.lambda_decrease;
                accepted = true;
                if (ftol_met || xtol_met) {
                    ++tol_hits;
                    if (tol_hits >= 2) {
                        res.converged = true;
                        return res;
                    }
                } else {
                    tol_hits = 0;
                }
            } else {
                lambda *= opt.lambda_increase;
            }
        }

        if (!accepted) {
            // damping saturated: no downhill step exists within the model, so
            // the point is a (local) minimum to the precision of the Jacobian.
            res.converged = true;
            return res;
        }
        if (res.ssq == zero) {
            res.converged = true;
            return res;
        }
    }
    return res;
}

/**
 * Levenberg-Marquardt with a central-difference Jacobian.
 *
 * @param f   residual map, x -> vector of length m
 * @param x0  starting point
 * @param m   number of residuals
 * @param opt tuning
 */
template <class T, class F>
LevmarResult<T> levmar(F f, const std::vector<T>& x0, std::size_t m,
                       const LevmarOptions<T>& opt) {
    const T step = opt.diff_step;
    // the differencing evaluations happen inside the Jacobian callable, so they
    // are tallied separately and folded into the reported count
    std::shared_ptr<unsigned> extra = std::make_shared<unsigned>(0u);
    LevmarResult<T> r = levmar_jac(
        f,
        [f, m, step, extra](const std::vector<T>& x) {
            *extra += 2u * static_cast<unsigned>(x.size());
            return levmar_jacobian_fd(f, x, m, step);
        },
        x0, m, opt);
    r.evaluations += *extra;
    return r;
}

/** levmar with the default tuning. */
template <class T, class F>
LevmarResult<T> levmar(F f, const std::vector<T>& x0, std::size_t m) {
    return levmar(f, x0, m, levmar_defaults<T>());
}

}  // namespace line

#endif  // LINE_UTIL_LEVMAR_H
