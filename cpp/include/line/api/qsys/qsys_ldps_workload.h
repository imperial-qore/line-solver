/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_LDPS_WORKLOAD_H
#define LINE_API_QSYS_QSYS_LDPS_WORKLOAD_H

/**
 * Stationary distribution of the unfinished work in a single-stage
 * load-dependent generalized-processor-sharing station with Poisson arrivals
 * and blocking. Port of matlab/src/api/qsys/qsys_ldps_workload.m, which is the
 * model of J.W. Cohen, "The multiple phase service network with generalized
 * processor sharing", Acta Informatica 12, 245-284 (1979), Sect. 9.
 *
 * MODEL. Poisson arrivals of rate lambda; blocking capacity N, an arrival
 * finding N requests present being lost without trace; each of x present
 * requests accrues service at rate f(x), so the stage completes work at total
 * rate x f(x); required service times iid with an absolutely continuous law B
 * of finite mean beta. The port is parametrized, as the reference is, by the
 * LINE load-dependent TOTAL rate scaling alpha(x) = x f(x), which is the
 * argument of setLoadDependence at a PS station.
 *
 * Cohen's eqs. (9.1)-(9.3):
 *   P{psi_t < psi} = sum_{h=0}^{N} p_h Psi^{h*}(psi)
 *   p_h propto rho^h/prod_{k=1}^{h} alpha(k),   rho = lambda beta
 *   Psi(psi)       = int_0^psi (1 - B(v))/beta dv
 * with Psi^{h*} the h-fold convolution and Psi^{0*} degenerate at zero, so the
 * workload has an ATOM of size p_0 at the origin. Substituting f(k) = alpha(k)/k
 * cancels the factorial that appears in Cohen's phi(h), leaving the familiar
 * load-dependent birth-death form for p. The state probabilities depend on B
 * only through beta; the workload depends on its shape only through the
 * equilibrium residual law Psi.
 *
 * WHAT THE PORT TAKES INSTEAD OF A Distribution OBJECT. The reference takes a
 * LINE Distribution and reads getMean, getSCV and evalCDF off it. The C++ port
 * has no model layer, so it takes exactly those three things: the mean, the
 * squared coefficient of variation (used only to size the default grid through
 * the mean equilibrium residual life beta(1+SCV)/2) and the CDF as a callable.
 * Nothing else of the Distribution interface is used by the reference either,
 * so this is the same function with its dependency made explicit.
 *
 * ACCURACY. The convolutions are formed on a uniform grid with the TRAPEZOIDAL
 * rule, not the rectangle rule that a bare convolution implies. The correction
 * matters because the equilibrium density does not vanish at the origin,
 * e(0) = 1/beta: without it each convolution over-counts by dt e(0) g_i, which
 * accumulates over h and drives the mixture CDF above one. Convergence is
 * second order in the step for an absolutely continuous B, the case Cohen
 * assumes, and falls to first order when B has an atom so that 1 - B is
 * discontinuous, Det being the extreme case; the reference measures 4.00x per
 * grid doubling for Exp and about 2x for Det, and the port reproduces both
 * regimes (see the tests).
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental. The birth-death
 * weights are accumulated in LOGS, exactly as the reference does, so a large
 * rho or a large N cannot overflow before normalization; that alone needs
 * log/exp. The quadrature is a tolerance-free fixed grid, so it introduces no
 * further requirement, but it does introduce a discretization error, which is
 * why the grid size is an explicit argument rather than hidden.
 *
 * NOT the weighted GPS/DPS discipline of SchedStrategy.GPS: this formula has no
 * per-class weights and does not represent them.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/qsys/qsys_quadrature.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/**
 * The three things Cohen's formula needs from the required-service-time law:
 * its mean, its squared coefficient of variation, and its CDF.
 */
template <class T>
struct WorkloadServiceLaw {
    T mean;
    T scv;
    std::function<T(const T&)> cdf;
};

/** Return value of qsys_ldps_workload, mirroring the three MATLAB outputs. */
template <class T>
struct LdpsWorkloadResult {
    std::vector<T> F;  ///< F(j) = P{psi <= t(j)}; F(0) = p(0) when t(0) = 0
    std::vector<T> t;  ///< grid at which F is reported
    std::vector<T> p;  ///< p(h) = P{x = h}, h = 0..N, the number in system
};

namespace ldps_detail {

/**
 * Convolution of two densities sampled on a uniform grid, trapezoidal rather
 * than rectangular:
 *   (f*g)(t_i) ~ dt [ sum_{j=0}^{i} f_j g_{i-j} - (f_0 g_i + f_i g_0)/2 ].
 */
template <class T>
std::vector<T> convtrap(const std::vector<T>& f, const std::vector<T>& g, const T& dt,
                        std::size_t n) {
    const T zero = num_traits<T>::from_int(0);
    const T half = num_traits<T>::from_rational(1, 2);
    std::vector<T> c(n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T s = zero;
        for (std::size_t j = 0; j <= i; ++j) s += f[j] * g[i - j];
        c[i] = dt * (s - half * (f[0] * g[i] + f[i] * g[0]));
    }
    return c;
}

}  // namespace ldps_detail

/**
 * Workload distribution of the load-dependent PS station with blocking.
 *
 * @param lambda Poisson arrival rate
 * @param B      required service time: mean, SCV and CDF
 * @param alpha  rate scaling alpha(n) = n f(n) for n = 1..N, at least N entries
 * @param N      blocking capacity
 * @param t      grid at which the CDF is wanted; empty for an automatic grid
 * @param ngrid  points of the internal uniform quadrature grid
 */
template <class T>
LdpsWorkloadResult<T> qsys_ldps_workload(const T& lambda, const WorkloadServiceLaw<T>& B,
                                         const std::vector<T>& alpha, std::size_t N,
                                         const std::vector<T>& t, std::size_t ngrid) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_ldps_workload accumulates the birth-death weights in logs");
    using std::exp;
    using std::log;
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    if (lambda <= zero)
        throw InputError("qsys_ldps_workload: lambda must be a positive arrival rate");
    if (B.mean <= zero)
        throw InputError("qsys_ldps_workload: the service law must have a finite positive mean");
    if (!B.cdf) throw InputError("qsys_ldps_workload: the service law has no CDF");
    if (N < 1) throw InputError("qsys_ldps_workload: N must be a positive blocking capacity");
    if (alpha.size() < N)
        throw InputError(
            "qsys_ldps_workload: alpha must supply the rate scaling for n = 1..N");
    for (std::size_t k = 0; k < N; ++k)
        if (alpha[k] <= zero)
            throw InputError(
                "qsys_ldps_workload: alpha(n) must be strictly positive, since every request in a "
                "busy stage is served at a positive rate");
    if (ngrid < 2) throw InputError("qsys_ldps_workload: ngrid must be at least 2");

    // Stationary number in system, eq. (9.1), accumulated in logs.
    const T beta = B.mean;
    const T rho = lambda * beta;
    std::vector<T> logw(N + 1, zero);
    for (std::size_t h = 1; h <= N; ++h) logw[h] = logw[h - 1] + log(rho) - log(alpha[h - 1]);
    T lmax = logw[0];
    for (const T& v : logw)
        if (v > lmax) lmax = v;
    std::vector<T> p(N + 1);
    T wsum = zero;
    for (std::size_t h = 0; h <= N; ++h) {
        p[h] = exp(logw[h] - lmax);
        wsum += p[h];
    }
    for (T& v : p) v /= wsum;

    // Grid. Psi has mean m1e, the mean equilibrium residual life of B, so the
    // largest h that carries non-negligible mass sizes the support.
    const T m1e = beta * (one + B.scv) / two;
    if (m1e <= zero)
        throw InputError(
            "qsys_ldps_workload: the equilibrium residual service time needs a finite second "
            "moment");
    const bool userGrid = !t.empty();
    T tmax;
    if (userGrid) {
        tmax = zero;
        for (const T& v : t) {
            if (v < zero) throw InputError("qsys_ldps_workload: t must be non-negative");
            if (v > tmax) tmax = v;
        }
        if (tmax <= zero) tmax = m1e;
    } else {
        std::size_t hmax = 0;
        const T floor_ = num_traits<T>::from_double(1e-12);
        for (std::size_t h = 0; h <= N; ++h)
            if (p[h] > floor_) hmax = h;
        if (hmax < 1) hmax = 1;
        const T hT = num_traits<T>::from_int(static_cast<long>(hmax));
        tmax = m1e * (hT + num_traits<T>::from_int(8) * sqrt(hT));
        const T floorT = num_traits<T>::from_int(8) * m1e;
        if (tmax < floorT) tmax = floorT;
    }

    std::vector<T> tg(ngrid);
    const T span = num_traits<T>::from_int(static_cast<long>(ngrid - 1));
    for (std::size_t i = 0; i < ngrid; ++i)
        tg[i] = tmax * num_traits<T>::from_int(static_cast<long>(i)) / span;
    const T dt = tg[1] - tg[0];

    // Equilibrium residual service density, eq. (9.3).
    std::vector<T> e(ngrid);
    for (std::size_t i = 0; i < ngrid; ++i) e[i] = (one - B.cdf(tg[i])) / beta;

    // Workload distribution, eq. (9.2). The h = 0 term is degenerate at zero
    // and contributes the atom p_0 across the whole non-negative grid.
    std::vector<T> Fg(ngrid, p[0]);
    std::vector<T> dens = e;
    for (std::size_t h = 1; h <= N; ++h) {
        if (h > 1) dens = ldps_detail::convtrap(dens, e, dt, ngrid);
        const std::vector<T> C = detail::num_cumtrapz(tg, dens);
        for (std::size_t i = 0; i < ngrid; ++i) Fg[i] += p[h] * C[i];
    }

    LdpsWorkloadResult<T> out;
    out.p = p;
    if (userGrid) {
        out.t = t;
        out.F.resize(t.size());
        for (std::size_t j = 0; j < t.size(); ++j) {
            // grid interpolation rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
            const T x = t[j];
            if (x >= tg[ngrid - 1]) {
                out.F[j] = Fg[ngrid - 1];
                continue;
            }
            const double pos = num_traits<T>::to_double(T(x / dt));
            std::size_t i = static_cast<std::size_t>(pos);
            if (i + 1 >= ngrid) i = ngrid - 2;
            const T w = (x - tg[i]) / dt;
            out.F[j] = Fg[i] + w * (Fg[i + 1] - Fg[i]);
        }
    } else {
        out.t = tg;
        out.F = Fg;
    }
    return out;
}

/** qsys_ldps_workload on the automatic grid with the reference default ngrid = 2001. */
template <class T>
LdpsWorkloadResult<T> qsys_ldps_workload(const T& lambda, const WorkloadServiceLaw<T>& B,
                                         const std::vector<T>& alpha, std::size_t N) {
    return qsys_ldps_workload(lambda, B, alpha, N, std::vector<T>(),
                              static_cast<std::size_t>(2001));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_LDPS_WORKLOAD_H
