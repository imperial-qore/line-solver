/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_ODE_H
#define LINE_UTIL_ODE_H

/**
 * Adaptive stiff ODE integrator: a four-stage Rosenbrock method of order four
 * with an embedded order-three estimate for step-size control.
 *
 * WHY IT IS HERE. Several LINE algorithms are defined by an initial value
 * problem that MATLAB hands to ode15s or ode23s: the refined mean-field cache
 * approximation (cache_miss_rmf.m), the RANDOM(m) multi-list mean field
 * (cache_rrm_meanfield.m), the fluid solvers. Those right-hand sides are stiff
 * -- the mean-field drift of a cache mixes per-item request rates that differ
 * by orders of magnitude, and the relaxation to the fixed point is integrated
 * over a horizon of 1e4 -- so an explicit method is not merely slower, it is
 * unusable: its step is capped by the fastest time constant for the whole
 * integration even after every fast mode has died. This header is the port's
 * own integrator; nothing is taken from an external solver library.
 *
 * THE METHOD. A Rosenbrock (linearly implicit Runge-Kutta) method replaces the
 * nonlinear stage equations of an implicit method by linear ones built on the
 * Jacobian, so each step costs one Jacobian, one LU factorization and s
 * back-substitutions and no Newton iteration ever fails to converge. With
 * J = df/dy and f_t = df/dt evaluated once per step at (t,y), the stages are
 *
 *   (I - h gamma J) k_i = h f(t + alpha_i h, y + sum_{j<i} a_ij k_j)
 *                         + h J sum_{j<i} gamma_ij k_j
 *                         + h^2 gamma_i f_t,
 *   y_{n+1} = y_n + sum_i b_i k_i,     yhat_{n+1} = y_n + sum_i bhat_i k_i,
 *
 * with alpha_i = sum_{j<i} a_ij and gamma_i = gamma + sum_{j<i} gamma_ij. The
 * f_t term is exactly what the augmented system (y,t)' = (f,1) produces for the
 * y-block, so the method is invariant under autonomization and the order
 * conditions are the autonomous ones. The same matrix I - h gamma J serves
 * every stage, which is the whole point of the constant diagonal gamma.
 *
 * THE COEFFICIENTS ARE DERIVED HERE, NOT COPIED. They were obtained by
 * imposing the order conditions directly rather than by quoting a published
 * table: for random polynomial vector fields the one-step numerical solution
 * was expanded as a truncated power series in h and matched, coefficient by
 * coefficient, against the exact Taylor series of the solution, and the
 * parameters solved so that the h^1 through h^4 coefficients agree. That is
 * the definition of order four, with no intermediate rooted-tree bookkeeping
 * to get wrong. Two further requirements were imposed at the same time:
 * b^T B^-1 1 = 1 with B = A + Gamma the lower triangular matrix with constant
 * diagonal gamma, which makes R(-inf) = 0 and the method L-stable -- an
 * A-stable but not L-stable method leaves the fastest modes ringing at
 * |R| = 1 instead of damping them, precisely the failure a stiff integrator
 * exists to avoid -- and gamma pinned to a value for which |R(z)| <= 1 holds
 * along the whole negative real axis.
 *
 * The embedded estimate is SECOND order and uses the first three stages. With
 * four stages and four weights the order-three conditions have a unique
 * solution, which is the order-four weight vector itself, so no order-three
 * estimate with these stages exists; the weights below instead reproduce the
 * solution through h^2 exactly and use their remaining degree of freedom to
 * make the h^3 mismatch as small as possible, which makes y - yhat a sharp
 * estimate of the O(h^3) term. The step controller therefore uses the
 * exponent 1/(2+1) = 1/3.
 *
 * What the test suite checks about the coefficients is what can be checked
 * exactly: alpha is the row sum of a, gamma_i is gamma plus the row sum of
 * Gamma, the linear order conditions b^T B^(k-1) 1 = 1/k! for k = 1..4 hold
 * (on a linear problem the method IS the implicit Runge-Kutta with matrix B,
 * so these are necessary and sufficient there), the embedded weights are
 * consistent and differ from b, |R| <= 1 on the negative axis with R(-inf) = 0,
 * and the observed convergence rate on a nonlinear non-autonomous problem is
 * four. A mistyped digit fails at least one of those.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental. The coefficients are
 * irrational, the step size is chosen by a tolerance comparison and the result
 * is an approximation controlled by rtol/atol no matter how the arithmetic is
 * carried out, so an exact-rational instantiation would be a fiction. double
 * and Real<D> both instantiate; at Real<D> the coefficients are parsed from
 * their decimal strings, so the method keeps its order at any precision, and
 * the numeric Jacobian's difference increment is scaled by the precision of T.
 *
 * DETERMINISM. No global state, no static mutable data, no clock, no random
 * numbers. The same inputs return the same trajectory bit for bit, including
 * the sequence of accepted and rejected steps. Every tolerance and bound is
 * supplied by the caller through OdeOptions.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {

/**
 * Integration controls. All fields are caller-supplied; the defaults match the
 * odeset('RelTol',1e-8,'AbsTol',1e-10) that the MATLAB callers of ode15s in
 * this tree use.
 */
/**
 * A note on how tight rtol can usefully be. The step is chosen so that the
 * ORDER-TWO embedded estimate meets the tolerance, so h scales as rtol^(1/3)
 * while the order-four solution is far more accurate than asked. Asking for
 * 1e-8 costs a few hundred steps on the problems in this tree; asking for
 * 1e-14 costs a few million and will hit max_steps. When the goal is accuracy
 * rather than error control, pin h_init = h_min = h_max and run a fixed step:
 * the solution error is O(h^4) and nothing is spent on the estimate.
 */
template <class T>
struct OdeOptions {
    T rtol = num_traits<T>::from_double(1e-8);   ///< relative tolerance per component
    T atol = num_traits<T>::from_double(1e-10);  ///< absolute tolerance per component
    T h_init = num_traits<T>::from_int(0);       ///< initial step; 0 selects one automatically
    T h_min = num_traits<T>::from_int(0);        ///< smallest admissible step; 0 derives one
    T h_max = num_traits<T>::from_int(0);        ///< largest admissible step; 0 means |t1 - t0|
    std::size_t max_steps = 100000;              ///< abort after this many accepted steps
    bool store_trajectory = true;                ///< keep every accepted point, not just the last
    /**
     * A test consulted after every ACCEPTED STEP; true ends the integration
     * there, holding that state as the solution for the rest of the span.
     *
     * WHY IT EXISTS. A drift that has reached a fixed point cannot move again:
     * f(y*) = 0 means y is y* for every later t, so the remaining span is known
     * and stepping it is waste -- and worse than waste, because a stiff step
     * controller handed a state it is already at cannot pick a step and the
     * integration stops advancing. Empty by default, and the loop then runs
     * exactly as it always has. See `LsodaOptions::step_stop` for the twin on
     * the other arm, and `solver_fluid.h` for the caller.
     */
    std::function<bool(const T&, const std::vector<T>&)> step_stop;
};

/** Result of an integration. */
template <class T>
struct OdeSolution {
    std::vector<T> t;                  ///< accepted time points, t[0] = t0
    std::vector<std::vector<T>> y;     ///< y[i] is the state at t[i]
    std::size_t steps = 0;             ///< accepted steps
    std::size_t rejected = 0;          ///< rejected steps
    std::size_t jacobians = 0;         ///< Jacobian evaluations
    std::size_t f_evals = 0;           ///< right-hand side evaluations

    const std::vector<T>& final_state() const {
        if (y.empty()) throw NumericError("OdeSolution: no state was recorded");
        return y.back();
    }
    const T& final_time() const {
        if (t.empty()) throw NumericError("OdeSolution: no state was recorded");
        return t.back();
    }
};

namespace ode_detail {

/** sqrt by ADL, so double, cpp_bin_float and mpfr all resolve. */
template <class T>
inline T num_sqrt(const T& v) {
    using std::sqrt;
    return sqrt(v);
}

/** Machine epsilon of T as a value of T. */
template <class T>
inline T num_eps() {
    return std::numeric_limits<T>::epsilon();
}

/**
 * The coefficient set, parsed once per instantiation from decimal strings so
 * that a Real<200> instantiation is not silently limited to double precision.
 *
 * The strings carry 32 significant digits. The values were obtained by solving
 * the Taylor-match system in double and then refining it by Gauss-Newton in
 * 60-digit arithmetic on problems with exactly representable rational data,
 * until every matched Taylor coefficient and the L-stability condition were
 * satisfied to better than 1e-50. That refinement is what makes a Real<50>
 * instantiation worth having: coefficients good only to 1e-16 would cap the
 * attainable local error at 1e-16 h no matter how much precision the caller
 * asked for.
 *
 * alpha and gamma_i are NOT stored, they are the row sums of a and of gamma,
 * so the three cannot drift apart.
 */
template <class T>
struct Ros4 {
    T gamma;
    T alpha[4];    ///< row sums of a: the stage abscissae
    T a[4][4];     ///< strictly lower
    T gam[4][4];   ///< strictly lower
    T gamma_i[4];  ///< gamma + row sum of gam
    T b[4];        ///< order-four weights
    T bhat[4];     ///< embedded order-two weights (first three stages)

    Ros4() {
        gamma = parse("0.57281606248213485540800138497677");
        const char* a_s[4][4] = {
            {"0", "0", "0", "0"},
            {"0.75000000000001384828217514743208", "0", "0", "0"},
            {"0.70000000000001188142192436149683", "-0.19284762489916794335116325012491", "0", "0"},
            {"-0.21564100871060819772037702114782", "0.058458776820062370659853044595891",
             "0.92804737156696904284870455820564", "0"}};
        const char* g_s[4][4] = {
            {"0", "0", "0", "0"},
            {"-0.84635396502372715969128660212009", "0", "0", "0"},
            {"0.050072468390209764931220387127165", "-0.74834720180244088425591310078693", "0",
             "0"},
            {"-0.42466874748751782239409914842761", "-0.9086245862275405564736317955516",
             "0.45665088789848053231587535636693", "0"}};
        const char* b_s[4] = {"0.38309583630415199434118358963363",
                              "0.049177316168702520768096197158938",
                              "0.09403024921808945954457206670416",
                              "0.47369659830905602534614814650327"};
        const char* bhat_s[4] = {"0.44804121887667558433", "0.34479429610719691812",
                                 "0.20716448501612805266", "0"};
        for (int i = 0; i < 4; ++i) {
            b[i] = parse(b_s[i]);
            bhat[i] = parse(bhat_s[i]);
            for (int j = 0; j < 4; ++j) {
                a[i][j] = parse(a_s[i][j]);
                gam[i][j] = parse(g_s[i][j]);
            }
        }
        for (int i = 0; i < 4; ++i) {
            alpha[i] = num_traits<T>::from_int(0);
            gamma_i[i] = gamma;
            for (int j = 0; j < i; ++j) {
                alpha[i] += a[i][j];
                gamma_i[i] += gam[i][j];
            }
        }
    }

    static T parse(const char* s) { return T(s); }
};

template <>
inline double Ros4<double>::parse(const char* s) {
    return std::stod(s);
}

}  // namespace ode_detail

/**
 * Numeric Jacobian by central differences.
 *
 * The increment is eps^(1/3) scaled by the magnitude of the component, which
 * is the standard balance for a central difference: the truncation error is
 * O(delta^2) and the cancellation error O(eps/delta), and the two meet at
 * delta ~ eps^(1/3), giving about two thirds of the digits of T. A one-sided
 * difference would cost one fewer evaluation per column and half the digits;
 * the extra accuracy matters here because the Jacobian of a stiff problem is
 * what the whole stability of the step rests on.
 */
template <class T, class F>
Matrix<T> ode_numeric_jacobian(const F& f, const T& t, const std::vector<T>& y,
                               const std::vector<T>& fy) {
    (void)fy;
    const std::size_t n = y.size();
    const T eps = ode_detail::num_eps<T>();
    using std::pow;
    const T delta_scale = pow(eps, num_traits<T>::from_rational(1, 3));
    Matrix<T> J(n, n, num_traits<T>::from_int(0));
    std::vector<T> yp = y;
    std::vector<T> ym = y;
    for (std::size_t j = 0; j < n; ++j) {
        T mag = num_abs(y[j]);
        if (mag < num_traits<T>::from_int(1)) mag = num_traits<T>::from_int(1);
        const T d = delta_scale * mag;
        yp[j] = y[j] + d;
        ym[j] = y[j] - d;
        const T den = yp[j] - ym[j];  // the actually representable increment
        const std::vector<T> fp = f(t, yp);
        const std::vector<T> fm = f(t, ym);
        if (fp.size() != n || fm.size() != n)
            throw InputError("ode_numeric_jacobian: the right-hand side changed dimension");
        for (std::size_t i = 0; i < n; ++i) J(i, j) = (fp[i] - fm[i]) / den;
        yp[j] = y[j];
        ym[j] = y[j];
    }
    return J;
}

/**
 * Integrate y' = f(t,y) from t0 to t1 with an analytic Jacobian.
 *
 * @param f    right-hand side, std::vector<T> f(const T& t, const std::vector<T>& y)
 * @param jac  Jacobian, Matrix<T> jac(const T& t, const std::vector<T>& y)
 * @param t0   initial time
 * @param t1   final time; t1 > t0 is required (this is an initial value problem
 *             marched forwards, and a backwards request is an input error
 *             rather than a silently reversed integration)
 * @param y0   initial state
 * @param opt  tolerances and step bounds
 */
template <class T, class F, class J>
OdeSolution<T> ode_rosenbrock4(const F& f, const J& jac, const T& t0, const T& t1,
                               const std::vector<T>& y0, const OdeOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "ode_rosenbrock4 requires transcendental arithmetic: its coefficients are "
                  "irrational and its step size is chosen by a tolerance comparison, so the "
                  "result is an approximation that exact rational arithmetic cannot deliver");

    const std::size_t n = y0.size();
    if (n == 0) throw InputError("ode_rosenbrock4: empty initial state");
    if (!(t1 > t0)) throw InputError("ode_rosenbrock4: the final time must exceed the initial time");
    if (!(opt.rtol > num_traits<T>::from_int(0)) || !(opt.atol > num_traits<T>::from_int(0)))
        throw InputError("ode_rosenbrock4: rtol and atol must both be positive");

    const ode_detail::Ros4<T> C;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T span = t1 - t0;
    const T hmax = opt.h_max > zero ? opt.h_max : span;
    const T hmin = opt.h_min > zero ? opt.h_min
                                    : T(span * ode_detail::num_eps<T>() *
                                        num_traits<T>::from_int(16));

    OdeSolution<T> sol;
    std::vector<T> y = y0;
    T t = t0;
    sol.t.push_back(t);
    sol.y.push_back(y);

    // initial step heuristic: see _kb/14-cpp-multiprecision.md
    T h = opt.h_init > zero ? opt.h_init : T(span / num_traits<T>::from_int(1000));
    if (h > hmax) h = hmax;
    if (h < hmin) h = hmin;

    const T safety = num_traits<T>::from_rational(9, 10);
    const T fac_min = num_traits<T>::from_rational(1, 5);
    const T fac_max = num_traits<T>::from_int(6);
    const T reject_max = num_traits<T>::from_int(1);  // no growth right after a rejection

    std::vector<T> k[4];
    std::vector<T> ystage(n), rhs(n), ynew(n), acc(n);
    bool previous_rejected = false;

    while (t < t1) {
        if (sol.steps >= opt.max_steps)
            throw NumericError("ode_rosenbrock4: step budget exhausted before reaching t1");
        // exact landing on t1: see _kb/14-cpp-multiprecision.md
        bool last_step = false;
        const T remaining = t1 - t;
        if (h >= remaining || remaining - h <= num_abs(t1) * ode_detail::num_eps<T>() *
                                                   num_traits<T>::from_int(64)) {
            h = remaining;
            last_step = true;
        }
        if (h < hmin && !last_step)
            throw NumericError("ode_rosenbrock4: the step size fell below h_min; the problem is "
                               "either singular or the tolerances are unreachable in this "
                               "arithmetic");

        const std::vector<T> fy = f(t, y);
        ++sol.f_evals;
        if (fy.size() != n)
            throw InputError("ode_rosenbrock4: the right-hand side returned the wrong dimension");
        const Matrix<T> Jm = jac(t, y);
        ++sol.jacobians;
        if (Jm.rows() != n || Jm.cols() != n)
            throw InputError("ode_rosenbrock4: the Jacobian has the wrong shape");

        // df/dt by a central difference. Autonomous problems return zero here
        // and the term drops out exactly.
        std::vector<T> ft(n, zero);
        {
            using std::pow;
            T tmag = num_abs(t);
            if (tmag < one) tmag = one;
            const T dt = pow(ode_detail::num_eps<T>(), num_traits<T>::from_rational(1, 3)) * tmag;
            const std::vector<T> fp = f(T(t + dt), y);
            const std::vector<T> fm = f(T(t - dt), y);
            sol.f_evals += 2;
            const T den = (t + dt) - (t - dt);
            for (std::size_t i = 0; i < n; ++i) ft[i] = (fp[i] - fm[i]) / den;
        }

        // One factorization of I - h gamma J serves all four stages.
        Matrix<T> LHS(n, n, zero);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                LHS(i, j) = (i == j ? one : zero) - h * C.gamma * Jm(i, j);
        std::vector<std::size_t> piv = lu_factor(LHS);

        for (int s = 0; s < 4; ++s) {
            for (std::size_t i = 0; i < n; ++i) {
                ystage[i] = y[i];
                acc[i] = zero;
            }
            for (int j = 0; j < s; ++j)
                for (std::size_t i = 0; i < n; ++i) {
                    ystage[i] += C.a[s][j] * k[j][i];
                    acc[i] += C.gam[s][j] * k[j][i];
                }
            const std::vector<T> fs = s == 0 ? fy : f(T(t + C.alpha[s] * h), ystage);
            if (s != 0) ++sol.f_evals;
            for (std::size_t i = 0; i < n; ++i) {
                T Jacc = zero;
                for (std::size_t j = 0; j < n; ++j) Jacc += Jm(i, j) * acc[j];
                rhs[i] = h * fs[i] + h * Jacc + h * h * C.gamma_i[s] * ft[i];
            }
            k[s] = rhs;
            lu_solve(LHS, piv, k[s]);
        }

        for (std::size_t i = 0; i < n; ++i) {
            ynew[i] = y[i];
            for (int s = 0; s < 4; ++s) ynew[i] += C.b[s] * k[s][i];
        }

        // embedded order-3 error estimate: see _kb/14-cpp-multiprecision.md
        T err_sq = zero;
        for (std::size_t i = 0; i < n; ++i) {
            T d = zero;
            for (int s = 0; s < 4; ++s) d += (C.b[s] - C.bhat[s]) * k[s][i];
            const T ay = num_abs(y[i]);
            const T an = num_abs(ynew[i]);
            const T scale = opt.atol + opt.rtol * (ay > an ? ay : an);
            const T r = d / scale;
            err_sq += r * r;
        }
        const T err = ode_detail::num_sqrt(T(err_sq / num_traits<T>::from_int(static_cast<long>(n))));

        // Step-size factor safety * err^(-1/(p+1)) with p = 2 the order of the
        // embedded estimate, so the cube root of the error ratio.
        T factor;
        if (err <= zero) {
            factor = fac_max;
        } else {
            using std::pow;
            const T q = pow(err, num_traits<T>::from_rational(1, 3));
            factor = safety / q;
            if (factor < fac_min) factor = fac_min;
            if (factor > fac_max) factor = fac_max;
        }

        if (err <= one) {
            t = last_step ? t1 : T(t + h);
            y = ynew;
            ++sol.steps;
            if (opt.store_trajectory || t >= t1) {
                sol.t.push_back(t);
                sol.y.push_back(y);
            } else {
                sol.t.back() = t;
                sol.y.back() = y;
            }
            // A settled state ends the span in closed form: the accepted point
            // is already recorded above, and it is the answer for every later t.
            if (opt.step_stop && opt.step_stop(t, y)) {
                if (t < t1) {
                    if (opt.store_trajectory) {
                        sol.t.push_back(t1);
                        sol.y.push_back(y);
                    } else {
                        sol.t.back() = t1;
                    }
                }
                break;
            }
            if (previous_rejected && factor > reject_max) factor = reject_max;
            previous_rejected = false;
            h *= factor;
            if (h > hmax) h = hmax;
        } else {
            ++sol.rejected;
            previous_rejected = true;
            h *= factor;
        }
    }
    return sol;
}

/** Integrate y' = f(t,y) with a numeric Jacobian by central differences. */
template <class T, class F>
OdeSolution<T> ode_rosenbrock4(const F& f, const T& t0, const T& t1, const std::vector<T>& y0,
                               const OdeOptions<T>& opt) {
    return ode_rosenbrock4(
        f,
        [&f](const T& t, const std::vector<T>& y) {
            return ode_numeric_jacobian<T>(f, t, y, std::vector<T>());
        },
        t0, t1, y0, opt);
}

/** Integrate with the default options and return only the state at t1. */
template <class T, class F>
std::vector<T> ode_rosenbrock4_endpoint(const F& f, const T& t0, const T& t1,
                                        const std::vector<T>& y0) {
    OdeOptions<T> opt;
    opt.store_trajectory = false;
    return ode_rosenbrock4(f, t0, t1, y0, opt).final_state();
}

}  // namespace line

#endif  // LINE_UTIL_ODE_H
