/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GTMTST_FLUID_H
#define LINE_API_QSYS_GTMTST_FLUID_H

/**
 * The Gt/Mt/st+GI many-server fluid queue, and the network of them.
 *
 * Templated port of matlab/src/api/qsys/qsys_gtmtst_fluid.m and
 * matlab/src/api/npfqn/npfqn_gtmtst_fluid.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gtmtst_fluid.java.
 *
 * Time-varying arrival rate lambda(t), staffing s(t), exponential service at
 * rate mu(t), general patience with ccdf F^c, unlimited waiting room.
 *
 * THE MODEL ALTERNATES BETWEEN TWO REGIMES and the algorithm is the bookkeeping
 * of that alternation:
 *
 *   UNDERLOADED  the queue is empty, every arrival enters service at once, and
 *                B' = lambda(t) - mu(t)B(t)                          (18, Mt form)
 *                ends when B reaches s while lambda > Gamma          (15)
 *   OVERLOADED   B(t) = s(t), fluid enters service at exactly
 *                Gamma(t) = s'(t) + s(t)mu(t)                        (13)
 *                q(t,x) = lambda(t-x)F^c(x) for x <= w(t)            (20)
 *                w'(t) = 1 - Gamma(t)/[lambda(t-w(t))F^c(w(t))]      (21)
 *                ends when w returns to 0 with lambda <= Gamma       (14)
 *
 * WHY w AND NOT Q. The queue content is a functional of w, but not the other way
 * round: two systems with the same Q and different age profiles abandon at
 * different rates. Tracking the boundary keeps the age profile exact, which is
 * what makes a general patience law admissible at all.
 *
 * THE NETWORK IS A FIXED POINT. lambda_j = lambda_j^0 + sum_i sigma_i P_ij with
 * sigma_i = mu_i B_i (23)-(24), iterated from the external rates; the nth
 * iterate is the fluid that has made n transitions, and the map is a monotone
 * contraction, so the rates increase to the fixed point. Only the SERVICE
 * COMPLETION flow is routed: abandoning fluid leaves the network, which is what
 * makes the traffic equations linear in sigma.
 *
 * ARITHMETIC. RK4 on a grid against a tolerance, so nothing is exact; the
 * instantiation is restricted to the transcendental types.
 *
 * Reference: Y. Liu, W. Whitt (2012). The Gt/GI/st+GI many-server fluid queue.
 * Queueing Systems 71, 405-444; Y. Liu, W. Whitt (2014). Algorithms for
 * time-varying networks of many-server fluid queues. INFORMS J. on Computing
 * 26(1), 59-73.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Trajectory of the Gt/Mt/st+GI fluid queue; every vector is on the time grid. */
template <class T>
struct QsysTvFluidResult {
    std::vector<T> times;         ///< the time grid
    std::vector<int> regime;      ///< 1 overloaded, 0 underloaded
    std::vector<T> B;             ///< fluid in service
    std::vector<T> Q;             ///< fluid in queue
    std::vector<T> X;             ///< B + Q
    std::vector<T> w;             ///< boundary waiting time
    std::vector<T> v;             ///< potential waiting time
    std::vector<T> sigma;         ///< service completion rate mu B
    std::vector<T> alpha;         ///< abandonment rate
    std::vector<T> utilization;   ///< B/s
    std::vector<T> arrivalRate;   ///< lambda on the grid
    std::vector<T> staffing;      ///< s on the grid
    std::vector<T> capacityRate;  ///< Gamma = s' + s mu
};

/** Options of `qsys_gtmtst_fluid`, all with the MATLAB defaults. */
template <class T>
struct TvFluidOptions {
    T dt = num_traits<T>::from_int(0);   ///< grid step; non-positive takes T/2000
    T B0 = num_traits<T>::from_int(0);   ///< fluid in service at time 0
    T w0 = num_traits<T>::from_int(0);   ///< boundary waiting time at time 0
    std::function<T(const T&)> sPrime;   ///< s'(t); differentiated numerically when empty
    std::function<T(const T&)> pdf;      ///< patience density; differenced when empty
    std::function<T(const T&)> lambdaPast;  ///< arrival rate before time 0
};

namespace detail {

/** Linear interpolation on an increasing grid, clamped at both ends. */
template <class T>
T tv_interp(const std::vector<T>& xs, const std::vector<T>& ys, const T& x) {
    const std::size_t n = xs.size();
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    std::size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        const std::size_t mid = (lo + hi) / 2;
        if (xs[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    const T f = (x - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] + f * (ys[hi] - ys[lo]);
}

/**
 * int_0^w lambda(t-x) WEIGHT(x) dx by Simpson: the ccdf gives the queue content
 * that has not abandoned, the density gives the abandonment rate.
 */
template <class T, class LamOf, class Weight>
T tv_integrate(LamOf&& lamOf, Weight&& weight, const T& ti, const T& wi, const T& dt) {
    const T zero = num_traits<T>::from_int(0);
    if (wi <= zero) return zero;
    std::size_t m = static_cast<std::size_t>(
        std::max(8.0, std::ceil(num_traits<T>::to_double(wi) / num_traits<T>::to_double(dt)) + 1));
    if (m % 2 == 1) ++m;
    const T h = wi / num_traits<T>::from_int(static_cast<long>(m));
    T sum = lamOf(ti) * weight(zero) + lamOf(T(ti - wi)) * weight(wi);
    for (std::size_t j = 1; j < m; ++j) {
        const T xx = num_traits<T>::from_int(static_cast<long>(j)) * h;
        sum += num_traits<T>::from_int(j % 2 == 1 ? 4 : 2) * lamOf(T(ti - xx)) * weight(xx);
    }
    return h / num_traits<T>::from_int(3) * sum;
}

}  // namespace detail

/**
 * @param lambdaFun    arrival rate lambda(t)
 * @param sFun         staffing s(t), positive
 * @param muFun        service rate mu(t), positive
 * @param patienceCcdf F^c(x) = P(patience > x)
 * @param T            horizon; the model is solved on [0,T]
 * @param opts         grid step, initial condition and the optional derivatives
 */
template <class Tv>
QsysTvFluidResult<Tv> qsys_gtmtst_fluid(const std::function<Tv(const Tv&)>& lambdaFun,
                                        const std::function<Tv(const Tv&)>& sFun,
                                        const std::function<Tv(const Tv&)>& muFun,
                                        const std::function<Tv(const Tv&)>& patienceCcdf,
                                        const Tv& T,
                                        const TvFluidOptions<Tv>& opts = TvFluidOptions<Tv>()) {
    static_assert(num_traits<Tv>::has_transcendental,
                  "qsys_gtmtst_fluid integrates on a grid, so it needs inexact arithmetic");
    const Tv zero = num_traits<Tv>::from_int(0);
    const Tv one = num_traits<Tv>::from_int(1);
    const Tv two = num_traits<Tv>::from_int(2);
    if (T <= zero) throw InputError("qsys_gtmtst_fluid: the horizon T must be positive");
    Tv dt = opts.dt;
    if (dt <= zero) dt = T / num_traits<Tv>::from_int(2000);
    const std::size_t n =
        static_cast<std::size_t>(std::llround(num_traits<Tv>::to_double(T / dt))) + 1;
    if (n < 2) throw InputError("qsys_gtmtst_fluid: the grid needs at least two points");

    QsysTvFluidResult<Tv> r;
    r.times.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        r.times[i] = T * num_traits<Tv>::from_int(static_cast<long>(i)) /
                     num_traits<Tv>::from_int(static_cast<long>(n - 1));
    const Tv step = r.times[1] - r.times[0];

    std::vector<Tv> lam(n), s(n), mu(n), sp(n), gamma(n);
    for (std::size_t i = 0; i < n; ++i) {
        lam[i] = lambdaFun(r.times[i]);
        s[i] = sFun(r.times[i]);
        mu[i] = muFun(r.times[i]);
        if (s[i] <= zero) throw InputError("qsys_gtmtst_fluid: the staffing must be positive");
        if (mu[i] <= zero) throw InputError("qsys_gtmtst_fluid: the service rate must be positive");
    }
    if (opts.sPrime) {
        for (std::size_t i = 0; i < n; ++i) sp[i] = opts.sPrime(r.times[i]);
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            if (i == 0) {
                sp[i] = (s[1] - s[0]) / step;
            } else if (i + 1 == n) {
                sp[i] = (s[n - 1] - s[n - 2]) / step;
            } else {
                sp[i] = (s[i + 1] - s[i - 1]) / (two * step);
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i) gamma[i] = sp[i] + s[i] * mu[i];   // Gamma(t), eq. (13)

    const std::function<Tv(const Tv&)> past = opts.lambdaPast ? opts.lambdaPast : lambdaFun;
    auto lamOf = [&](const Tv& u) { return u < zero ? past(u) : lambdaFun(u); };
    std::function<Tv(const Tv&)> pdf = opts.pdf;
    if (!pdf) {
        pdf = [&patienceCcdf](const Tv& x) {
            const Tv h = num_traits<Tv>::from_double(1e-6);
            const Tv lo = x - h < num_traits<Tv>::from_int(0) ? num_traits<Tv>::from_int(0) : Tv(x - h);
            const Tv d = (patienceCcdf(lo) - patienceCcdf(Tv(x + h))) / (num_traits<Tv>::from_int(2) * h);
            return d < num_traits<Tv>::from_int(0) ? num_traits<Tv>::from_int(0) : d;
        };
    }

    r.regime.assign(n, 0);
    r.B.assign(n, zero);
    r.Q.assign(n, zero);
    r.w.assign(n, zero);
    r.alpha.assign(n, zero);
    r.B[0] = opts.B0;
    r.w[0] = opts.w0;
    const bool over0 =
        opts.w0 > zero || (opts.B0 >= s[0] - num_traits<Tv>::from_double(1e-12) && lam[0] > gamma[0]);
    r.regime[0] = over0 ? 1 : 0;
    if (over0) r.B[0] = s[0];
    r.Q[0] = detail::tv_integrate<Tv>(lamOf, patienceCcdf, r.times[0], r.w[0], step);
    r.alpha[0] = detail::tv_integrate<Tv>(lamOf, pdf, r.times[0], r.w[0], step);

    auto wdot = [&](const Tv& tt, const Tv& ww) {
        const Tv den = lamOf(Tv(tt - ww)) * patienceCcdf(ww);
        // No fluid of that age survives, so the boundary can only advance with
        // the clock.
        if (den <= zero) return one;
        return Tv(one - detail::tv_interp(r.times, gamma, tt) / den);
    };

    for (std::size_t i = 0; i + 1 < n; ++i) {
        Tv bNext, wNext;
        if (r.regime[i] == 0) {
            // Underloaded: B' = lambda - mu B, by RK4 on the grid step.
            auto f = [&](const Tv& tt, const Tv& bb) {
                return detail::tv_interp(r.times, lam, tt) - detail::tv_interp(r.times, mu, tt) * bb;
            };
            const Tv k1 = f(r.times[i], r.B[i]);
            const Tv k2 = f(Tv(r.times[i] + step / two), Tv(r.B[i] + step * k1 / two));
            const Tv k3 = f(Tv(r.times[i] + step / two), Tv(r.B[i] + step * k2 / two));
            const Tv k4 = f(Tv(r.times[i] + step), Tv(r.B[i] + step * k3));
            bNext = r.B[i] + step * (k1 + two * k2 + two * k3 + k4) / num_traits<Tv>::from_int(6);
            wNext = zero;
            if (bNext >= s[i + 1] && lam[i + 1] > gamma[i + 1]) {
                // The servers just filled and the input outruns the freed
                // capacity: eq. (15), the underloaded interval ends here.
                bNext = s[i + 1];
                r.regime[i + 1] = 1;
            } else {
                r.regime[i + 1] = 0;
                if (bNext > s[i + 1]) bNext = s[i + 1];
            }
        } else {
            // Overloaded: B = s and the boundary moves by eq. (21).
            const Tv k1 = wdot(r.times[i], r.w[i]);
            const Tv w2 = r.w[i] + step * k1 / two;
            const Tv k2 = wdot(Tv(r.times[i] + step / two), w2 < zero ? zero : w2);
            const Tv w3 = r.w[i] + step * k2 / two;
            const Tv k3 = wdot(Tv(r.times[i] + step / two), w3 < zero ? zero : w3);
            const Tv w4 = r.w[i] + step * k3;
            const Tv k4 = wdot(Tv(r.times[i] + step), w4 < zero ? zero : w4);
            wNext = r.w[i] + step * (k1 + two * k2 + two * k3 + k4) / num_traits<Tv>::from_int(6);
            bNext = s[i + 1];
            if (wNext <= zero && lam[i + 1] <= gamma[i + 1]) {
                // The queue has drained and the input no longer outruns the
                // freed capacity: eq. (14), the overloaded interval ends here.
                wNext = zero;
                r.regime[i + 1] = 0;
            } else {
                if (wNext < zero) wNext = zero;
                r.regime[i + 1] = 1;
            }
        }
        r.B[i + 1] = bNext;
        r.w[i + 1] = wNext;
        if (r.regime[i + 1] == 1) {
            r.Q[i + 1] = detail::tv_integrate<Tv>(lamOf, patienceCcdf, r.times[i + 1], wNext, step);
            r.alpha[i + 1] = detail::tv_integrate<Tv>(lamOf, pdf, r.times[i + 1], wNext, step);
        }
    }

    r.X.resize(n);
    r.sigma.resize(n);
    r.utilization.resize(n);
    std::vector<Tv> entry(n);
    for (std::size_t i = 0; i < n; ++i) {
        r.sigma[i] = mu[i] * r.B[i];             // service completion rate, eq. (3)
        r.utilization[i] = r.B[i] / s[i];
        r.X[i] = r.B[i] + r.Q[i];
        entry[i] = r.times[i] - r.w[i];
    }
    // The potential waiting time of an arrival at t is the u-t at which the
    // boundary reaches it, i.e. the solution of u - w(u) = t. That map is
    // non-decreasing, so one interpolation inverts it.
    r.v.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Tv u = detail::tv_interp(entry, r.times, r.times[i]);
        r.v[i] = u - r.times[i] < zero ? zero : Tv(u - r.times[i]);
    }
    r.arrivalRate = lam;
    r.staffing = s;
    r.capacityRate = gamma;
    return r;
}

/** Result of the network solve. */
template <class T>
struct NpfqnTvFluidResult {
    std::vector<T> times;                        ///< the time grid
    std::vector<QsysTvFluidResult<T>> queues;    ///< the per-queue trajectories
    std::vector<std::vector<T>> arrivalRates;    ///< converged total rates, one row per queue
    std::size_t iterations = 0;                  ///< iterations of the traffic-rate fixed point
    T residual;                                  ///< sup-norm change at the last iteration
};

/**
 * A time-varying open network of many-server fluid queues with abandonment.
 *
 * @param lambdaFuns    external arrival rate of each queue
 * @param sFuns         staffing of each queue
 * @param muFuns        service rate of each queue
 * @param patienceCcdfs patience ccdf of each queue
 * @param P             routing proportions, substochastic
 * @param T             horizon
 * @param dt            grid step; non-positive takes T/2000
 * @param B0            initial fluid in service, empty for an empty network
 * @param w0            initial boundary waiting times, empty for an empty network
 * @param tol           sup-norm tolerance on the arrival-rate iteration
 * @param maxIter       cap on the iterations
 */
template <class Tv>
NpfqnTvFluidResult<Tv> npfqn_gtmtst_fluid(
    const std::vector<std::function<Tv(const Tv&)>>& lambdaFuns,
    const std::vector<std::function<Tv(const Tv&)>>& sFuns,
    const std::vector<std::function<Tv(const Tv&)>>& muFuns,
    const std::vector<std::function<Tv(const Tv&)>>& patienceCcdfs,
    const std::vector<std::vector<Tv>>& P, const Tv& T, const Tv& dt = num_traits<Tv>::from_int(0),
    const std::vector<Tv>& B0 = std::vector<Tv>(), const std::vector<Tv>& w0 = std::vector<Tv>(),
    double tol = 1e-6, std::size_t maxIter = 100) {
    static_assert(num_traits<Tv>::has_transcendental,
                  "npfqn_gtmtst_fluid integrates on a grid, so it needs inexact arithmetic");
    const Tv zero = num_traits<Tv>::from_int(0);
    const Tv one = num_traits<Tv>::from_int(1);
    const std::size_t m = lambdaFuns.size();
    if (sFuns.size() != m || muFuns.size() != m || patienceCcdfs.size() != m)
        throw InputError("npfqn_gtmtst_fluid: every queue needs an arrival rate, a staffing, a "
                         "service rate and a patience law");
    if (T <= zero) throw InputError("npfqn_gtmtst_fluid: the horizon T must be positive");
    Tv step = dt;
    if (step <= zero) step = T / num_traits<Tv>::from_int(2000);
    const std::size_t n =
        static_cast<std::size_t>(std::llround(num_traits<Tv>::to_double(T / step))) + 1;
    if (P.size() != m) throw InputError("npfqn_gtmtst_fluid: the routing matrix must be m x m");
    for (std::size_t i = 0; i < m; ++i) {
        if (P[i].size() != m)
            throw InputError("npfqn_gtmtst_fluid: the routing matrix must be m x m");
        Tv row = zero;
        for (std::size_t j = 0; j < m; ++j) {
            if (P[i][j] < -num_traits<Tv>::from_double(1e-12))
                throw InputError("npfqn_gtmtst_fluid: the routing matrix must be non-negative");
            row += P[i][j];
        }
        if (row > one + num_traits<Tv>::from_double(1e-9))
            throw InputError("npfqn_gtmtst_fluid: the routing matrix must be substochastic");
    }

    NpfqnTvFluidResult<Tv> out;
    out.times.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        out.times[i] = T * num_traits<Tv>::from_int(static_cast<long>(i)) /
                       num_traits<Tv>::from_int(static_cast<long>(n - 1));

    std::vector<std::vector<Tv>> ext(m, std::vector<Tv>(n, zero));
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t k = 0; k < n; ++k) ext[i][k] = lambdaFuns[i](out.times[k]);
    std::vector<std::vector<Tv>> lam = ext;

    out.residual = num_traits<Tv>::from_double(1e300);
    for (std::size_t iter = 1; iter <= maxIter; ++iter) {
        out.iterations = iter;
        out.queues.clear();
        std::vector<std::vector<Tv>> sigma(m, std::vector<Tv>(n, zero));
        for (std::size_t i = 0; i < m; ++i) {
            const std::vector<Tv>& row = lam[i];
            const std::vector<Tv>& grid = out.times;
            std::function<Tv(const Tv&)> fun = [&grid, &row](const Tv& u) {
                return detail::tv_interp(grid, row, u);
            };
            TvFluidOptions<Tv> o;
            o.dt = step;
            o.B0 = B0.empty() ? zero : B0[i];
            o.w0 = w0.empty() ? zero : w0[i];
            QsysTvFluidResult<Tv> res =
                qsys_gtmtst_fluid<Tv>(fun, sFuns[i], muFuns[i], patienceCcdfs[i], T, o);
            sigma[i] = res.sigma;
            out.queues.push_back(std::move(res));
        }
        // lambda_j = lambda_j^0 + sum_i sigma_i P_ij, eqs. (23)-(24).
        Tv diff = zero;
        std::vector<std::vector<Tv>> newlam = ext;
        for (std::size_t j = 0; j < m; ++j) {
            for (std::size_t k = 0; k < n; ++k) {
                for (std::size_t i = 0; i < m; ++i) newlam[j][k] += sigma[i][k] * P[i][j];
                const Tv d = newlam[j][k] - lam[j][k];
                const Tv ad = d < zero ? Tv(-d) : d;
                if (ad > diff) diff = ad;
            }
        }
        lam = newlam;
        out.residual = diff;
        if (num_traits<Tv>::to_double(diff) < tol) break;
    }
    out.arrivalRates = lam;
    return out;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GTMTST_FLUID_H
