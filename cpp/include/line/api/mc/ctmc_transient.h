/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_TRANSIENT_H
#define LINE_API_MC_CTMC_TRANSIENT_H

/**
 * Transient distribution of a CTMC over a time interval, by integrating the
 * forward equations d pi/dt = pi Q.
 *
 * Templated port of matlab/src/api/mc/ctmc_transient.m (and the richer
 * kpctoolbox copy). The reference integrates with ode23, MATLAB's
 * Bogacki-Shampine 3(2) pair with first-same-as-last, adaptive step and
 * default tolerances RelTol 1e-3, AbsTol 1e-6; the same pair, the same step
 * controller and the same defaults are reproduced here, so the returned time
 * grid is the solver's own accepted steps rather than a fixed grid.
 *
 * MATLAB-VS-JAVA DISAGREEMENT. jline.api.mc.Ctmc_transient integrates the same
 * equations with LSODA at absolute and relative tolerance 1e-6, so the two
 * references return DIFFERENT time grids and solutions agreeing only to the
 * looser of the two tolerances. This port follows MATLAB, the ground truth.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC. The step controller raises the error
 * ratio to the power 1/3, and adaptive integration is an approximation with a
 * tolerance rather than a finite exact computation: there is no exact value of
 * pi(t) for a general rational Q, exp(Qt) not being a rational function of t.
 * For a tightly controlled transient use ctmc_foxglynn or ctmc_uniformization
 * at high precision, which bound their truncation error explicitly.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct TransientResult {
    std::vector<T> t;  ///< accepted time points, the first being t0
    Matrix<T> pi;      ///< one row per time point
};

namespace detail {

/**
 * Bogacki-Shampine 3(2) integrator with MATLAB's ode23 step control, for the
 * autonomous system y' = f(y).
 *
 * @param f    right-hand side, called as f(y, dy) and writing dy
 * @param rtol relative tolerance (MATLAB RelTol, default 1e-3)
 * @param atol absolute tolerance (MATLAB AbsTol, default 1e-6)
 * @param t0 start of the integration horizon
 * @param t1 end of the integration horizon
 * @param y0 initial condition
 * @param tout out: the accepted time points
 * @param yout out: the solution at each accepted time point
 */
template <class T, class F>
void ode23(F f, const T& t0, const T& t1, const std::vector<T>& y0, double rtol, double atol,
           std::vector<T>& tout, std::vector<std::vector<T>>& yout) {
    const std::size_t n = y0.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T third = num_traits<T>::from_rational(1, 3);
    const T threshold = num_traits<T>::from_double(atol / rtol);
    const T rtolT = num_traits<T>::from_double(rtol);
    const T eps = std::numeric_limits<T>::epsilon();
    const T span = t1 - t0;
    if (!(span > zero)) throw InputError("ode23: the time interval must have positive length");
    const T hmax = span / num_traits<T>::from_int(10);

    std::vector<T> y = y0, f1(n), f2(n), f3(n), f4(n), ytmp(n), ynew(n);
    f(y, f1);

    // Initial step, selected as MATLAB's ode23 does.
    using std::pow;
    T h = hmax;
    {
        T rh = zero;
        for (std::size_t i = 0; i < n; ++i) {
            T d = num_abs(T(y[i]));
            if (d < threshold) d = threshold;
            const T v = num_abs(T(f1[i] / d));
            if (v > rh) rh = v;
        }
        rh /= num_traits<T>::from_rational(4, 5) * pow(rtolT, third);
        if (h * rh > one) h = one / rh;
    }

    T t = t0;
    tout.assign(1, t0);
    yout.assign(1, y0);
    bool done = false;
    while (!done) {
        const T hmin = num_traits<T>::from_int(16) * eps * (num_abs(T(t)) + one);
        if (h > hmax) h = hmax;
        if (h < hmin) h = hmin;
        if (num_traits<T>::from_rational(11, 10) * h >= t1 - t) {
            h = t1 - t;
            done = true;
        }

        bool nofailed = true;
        T err = zero;
        for (;;) {
            for (std::size_t i = 0; i < n; ++i) ytmp[i] = y[i] + h * f1[i] / num_traits<T>::from_int(2);
            f(ytmp, f2);
            for (std::size_t i = 0; i < n; ++i)
                ytmp[i] = y[i] + h * num_traits<T>::from_rational(3, 4) * f2[i];
            f(ytmp, f3);
            for (std::size_t i = 0; i < n; ++i)
                ynew[i] = y[i] + h *
                                     (num_traits<T>::from_int(2) * f1[i] + num_traits<T>::from_int(3) * f2[i] +
                                      num_traits<T>::from_int(4) * f3[i]) /
                                     num_traits<T>::from_int(9);
            f(ynew, f4);

            err = zero;
            for (std::size_t i = 0; i < n; ++i) {
                T d = num_abs(T(y[i]));
                const T dn = num_abs(T(ynew[i]));
                if (dn > d) d = dn;
                if (d < threshold) d = threshold;
                const T e = (num_traits<T>::from_int(-5) * f1[i] + num_traits<T>::from_int(6) * f2[i] +
                             num_traits<T>::from_int(8) * f3[i] - num_traits<T>::from_int(9) * f4[i]) /
                            num_traits<T>::from_int(72);
                const T v = num_abs(T(e / d));
                if (v > err) err = v;
            }
            err *= num_abs(T(h));

            if (!(err > rtolT)) break;
            if (h <= hmin)
                throw NumericError("ode23: step size underflow, the system is too stiff for ode23");
            T fac = num_traits<T>::from_rational(4, 5) * pow(T(rtolT / err), third);
            if (fac < num_traits<T>::from_rational(1, 2)) fac = num_traits<T>::from_rational(1, 2);
            h *= fac;
            if (h < hmin) h = hmin;
            nofailed = false;
            done = false;
        }

        t += h;
        y = ynew;
        f1 = f4;  // first-same-as-last
        tout.push_back(t);
        yout.push_back(y);
        if (done) break;
        // A step that was accepted without a prior rejection may grow, by at
        // most a factor of five; one that failed keeps its reduced size.
        if (nofailed) {
            if (err == zero) {
                h *= num_traits<T>::from_int(5);
            } else {
                const T temp = num_traits<T>::from_rational(5, 4) * pow(T(err / rtolT), third);
                if (temp > num_traits<T>::from_rational(1, 5))
                    h /= temp;
                else
                    h *= num_traits<T>::from_int(5);
            }
        }
    }
}

}  // namespace detail

/**
 * @param Q   generator
 * @param pi0 initial distribution (row vector)
 * @param t0  initial time
 * @param t1  final time
 * @param rtol relative tolerance of the integrator (MATLAB RelTol, 1e-3)
 * @param atol absolute tolerance of the integrator (MATLAB AbsTol, 1e-6)
 */
template <class T>
TransientResult<T> ctmc_transient(const Matrix<T>& Q, const std::vector<T>& pi0, const T& t0, const T& t1,
                                  double rtol = 1e-3, double atol = 1e-6) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_transient requires transcendental arithmetic: the ode23 step controller "
                  "raises the error ratio to the power 1/3, and the result is an approximation of "
                  "pi0 exp(Qt) governed by a tolerance rather than an exact quantity");
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_transient: generator is not square");
    if (pi0.size() != n) throw InputError("ctmc_transient: pi0 has the wrong length");

    std::vector<T> tv;
    std::vector<std::vector<T>> yv;
    detail::ode23<T>(
        [&Q, n](const std::vector<T>& y, std::vector<T>& dy) {
            for (std::size_t j = 0; j < n; ++j) {
                T s = num_traits<T>::from_int(0);
                for (std::size_t i = 0; i < n; ++i) s += y[i] * Q(i, j);
                dy[j] = s;
            }
        },
        t0, t1, pi0, rtol, atol, tv, yv);

    TransientResult<T> r;
    r.t = tv;
    r.pi = Matrix<T>(yv.size(), n);
    for (std::size_t k = 0; k < yv.size(); ++k)
        for (std::size_t j = 0; j < n; ++j) r.pi(k, j) = yv[k][j];
    return r;
}

/**
 * Overload starting from the uniform distribution, as MATLAB's short forms do.
 * The initial time stays explicit: an overload taking (Q, pi0, t1) would be
 * ambiguous with the four-argument form at T = double, since the tolerances are
 * doubles too.
 */
template <class T>
TransientResult<T> ctmc_transient(const Matrix<T>& Q, const T& t0, const T& t1, double rtol = 1e-3,
                                  double atol = 1e-6) {
    const std::size_t n = Q.rows();
    if (n == 0) throw InputError("ctmc_transient: empty generator");
    const std::vector<T> pi0(n, num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(n)));
    return ctmc_transient(Q, pi0, t0, t1, rtol, atol);
}

/**
 * Resample an adaptive transient onto the uniform grid `t0 : dt : t1`.
 *
 * `options.timestep` of the reference. `ctmc_transient.m` passes the grid
 * STRAIGHT INTO ode23 as its tspan, and MATLAB's integrator then reports the
 * solution at exactly those points by evaluating its own dense output -- it does
 * not change the steps it takes, only where it reports them. This does the same
 * thing explicitly: the adaptive solve is untouched, and the answer at a grid
 * point is the cubic Hermite interpolant of the bracketing pair.
 *
 * THE DERIVATIVES ARE EXACT AND NOT DIFFERENCED. `d pi/dt = pi Q` holds at every
 * stored point, so the Hermite data is the solution and its true derivative
 * rather than a secant estimate; that is what makes this the same order as
 * ode23's own dense output instead of a linear interpolation dressed up as one.
 *
 * `t1` IS ALWAYS THE LAST POINT even when the step does not divide the horizon,
 * as `ctmc_transient.m` appends it: a transient reported to 9.9 when 10 was
 * asked for is a different answer, not a rounded one.
 */
template <class T>
TransientResult<T> ctmc_transient_on_grid(const Matrix<T>& Q, const TransientResult<T>& r,
                                          const std::vector<T>& grid);

template <class T>
TransientResult<T> ctmc_transient_on_grid(const Matrix<T>& Q, const TransientResult<T>& r,
                                          const T& t0, const T& t1, const T& dt) {
    if (!(num_traits<T>::to_double(dt) > 0.0))
        throw InputError("ctmc_transient_on_grid: the timestep must be positive");
    if (r.t.empty()) return r;

    std::vector<T> grid;
    const double d0 = num_traits<T>::to_double(t0), d1 = num_traits<T>::to_double(t1),
                 dd = num_traits<T>::to_double(dt);
    for (double x = d0; x <= d1 + 1e-12 * (d1 - d0); x += dd)
        grid.push_back(num_traits<T>::from_double(x));
    if (grid.empty() || num_traits<T>::to_double(grid.back()) < d1) grid.push_back(t1);
    return ctmc_transient_on_grid(Q, r, grid);
}

/**
 * The same resampling onto an ARBITRARY grid, which a uniform step cannot
 * express.
 *
 * A caller that integrates a quantity AGAINST the trajectory needs the points
 * its integrand asks for, not the points the step controller happened to stop
 * at: the environment coupling forms a Riemann-Stieltjes sum of each stage's
 * transient against the holding-time CDF, and the grid that resolves that CDF
 * (`refine_grid`, 90% of its points under 5*E[S]) is not uniform. Resampling
 * here rather than re-integrating keeps ONE integration behind every grid.
 *
 * LINEAR, NOT HERMITE, AND DELIBERATELY THE LESS ACCURATE CHOICE. Every
 * codebase resamples here rather than re-integrating, but MATLAB's
 * `refineForCdf_`, the JAR's `CdfGrid.on` and native python all interpolate
 * LINEARLY, and this port used a cubic Hermite built from the exact derivative
 * `y' = yQ`. That is the better interpolant and it is what made
 * renv_threestages_repairmen read Queue1 QLen 0.83092 where the other three read
 * 0.83053: a 4.7e-4 disagreement that is entirely the difference between the two
 * quadratures, on a row whose gate is 3.2e-4. Parity measures whether the
 * codebases give the SAME answer, so the interpolant is aligned rather than the
 * golden rebased onto the more accurate one. Changed 2026-08-13 on the user's
 * decision; if this is ever revisited, note that the metrics read off pi(t) are
 * LINEAR functionals of it, so interpolating pi linearly and interpolating each
 * metric linearly are the same operation -- which is why matching the
 * reference here is enough to match it in every measure.
 */
template <class T>
TransientResult<T> ctmc_transient_on_grid(const Matrix<T>& Q, const TransientResult<T>& r,
                                          const std::vector<T>& grid) {
    if (r.t.empty() || grid.empty()) return r;
    const std::size_t n = Q.rows();
    (void)Q;  // the linear interpolant needs no derivative, so no Q

    TransientResult<T> out;
    out.t = grid;
    out.pi = Matrix<T>(grid.size(), n);
    std::size_t k = 0;
    for (std::size_t g = 0; g < grid.size(); ++g) {
        const double x = num_traits<T>::to_double(grid[g]);
        while (k + 2 < r.t.size() && num_traits<T>::to_double(r.t[k + 1]) < x) ++k;
        const double a = num_traits<T>::to_double(r.t[k]);
        const std::size_t k1 = (k + 1 < r.t.size()) ? k + 1 : k;
        const double b = num_traits<T>::to_double(r.t[k1]);
        if (k1 == k || b <= a) {
            for (std::size_t j = 0; j < n; ++j) out.pi(g, j) = r.pi(k, j);
            continue;
        }
        // CLAMPED at both ends, as the reference's own resamplers are: a grid
        // point outside the stored span takes the nearest endpoint rather than
        // an extrapolation, which on a probability vector could leave the
        // simplex.
        double u = (x - a) / (b - a);
        if (u < 0.0) u = 0.0;
        if (u > 1.0) u = 1.0;
        for (std::size_t j = 0; j < n; ++j)
            out.pi(g, j) = T(num_traits<T>::from_double(1.0 - u) * r.pi(k, j) +
                             num_traits<T>::from_double(u) * r.pi(k1, j));
    }
    return out;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_TRANSIENT_H
