/**
 * @file cache_sojourn.h
 * @brief Grid-independent sojourn-weighted average of a cache mean-field transient.
 *
 * Twin of MATLAB `cache_sojourn_ode.m`, JAR `CacheSojourn` and Python
 * `api/cache/sojourn.py`. For a drift dx/dt = f(x) and a phase-type clock with
 * row phase vector phi, dphi/dt = phi A, density g(t) = phi(t) c, it returns
 *
 *     xbar = int_{t0}^{t1} x g dt / int_{t0}^{t1} g dt,   wtot = int_{t0}^{t1} g dt,
 *
 * by augmenting the state with phi, int g x and int g, so the integrator itself
 * carries the integral. The value is accurate to the ODE tolerance and does NOT
 * depend on any output grid: a Riemann sum over the grid did, and that is what
 * made the adaptive (MATLAB, C++) and fixed-grid (JAR, Python) ENV cache mean
 * fields differ in the third digit.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

#ifndef LINE_API_CACHE_CACHE_SOJOURN_H
#define LINE_API_CACHE_CACHE_SOJOURN_H

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/matrix.h"
#include "line/util/ode.h"

namespace line {
namespace cache {

/** A phase-type clock: dphi/dt = phi A, phi(t0) = phi0, density phi c. Empty = no clock. */
template <class T>
struct CacheSojournClock {
    Matrix<T> A;
    std::vector<T> phi0;
    std::vector<T> c;
    bool empty() const { return phi0.empty(); }
};

/** A sojourn average: `xbar` (empty when `wtot` is not positive), clock mass, per-user miss rate. */
template <class T>
struct CacheSojournResult {
    std::vector<T> xbar;
    T wtot = num_traits<T>::from_int(0);
    std::vector<T> MU;  ///< per-user miss rate of xbar, filled by the policy
    bool empty() const { return xbar.empty(); }
};

/**
 * The holding-time clock of a random-environment stage whose holding time has
 * sub-generator D0 and entry vector pie, in a drift time unit that is `lam` times
 * real time: A = D0/lam, c = -D0 1/lam, phi0 = pie, so int_0^T g = F(T/lam).
 */
template <class T>
CacheSojournClock<T> cache_sojourn_clock(const Matrix<T>& D0, const std::vector<T>& pie, double lam) {
    CacheSojournClock<T> clk;
    const std::size_t nph = D0.rows();
    const T L = num_traits<T>::from_double(lam);
    clk.A = Matrix<T>(nph, nph, num_traits<T>::from_int(0));
    clk.c.assign(nph, num_traits<T>::from_int(0));
    clk.phi0 = pie;
    for (std::size_t i = 0; i < nph; ++i) {
        T rs = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < nph; ++j) {
            clk.A(i, j) = T(D0(i, j) / L);
            rs = T(rs + D0(i, j));
        }
        clk.c[i] = T(-rs / L);
    }
    return clk;
}

/** Integrates dx/dt = f(t, x) over [t0, t1] from x0 together with `clk`; see the file comment. */
template <class T, class F>
CacheSojournResult<T> cache_sojourn_ode(const F& f, const T& t0, const T& t1,
                                        const std::vector<T>& x0, const CacheSojournClock<T>& clk) {
    const std::size_t n = x0.size(), nph = clk.phi0.size();
    const std::size_t dim = 2 * n + nph + 1;
    const T zero = num_traits<T>::from_int(0);
    const auto rhs = [&](const T& t, const std::vector<T>& z) {
        std::vector<T> x(z.begin(), z.begin() + static_cast<std::ptrdiff_t>(n));
        const std::vector<T> dx = f(t, x);
        std::vector<T> dz(dim, zero);
        for (std::size_t i = 0; i < n; ++i) dz[i] = dx[i];
        T g = zero;
        for (std::size_t j = 0; j < nph; ++j) {
            T s = zero;
            for (std::size_t i = 0; i < nph; ++i) s = T(s + T(z[n + i] * clk.A(i, j)));
            dz[n + j] = s;
            g = T(g + T(z[n + j] * clk.c[j]));
        }
        for (std::size_t i = 0; i < n; ++i) dz[n + nph + i] = T(g * x[i]);
        dz[dim - 1] = g;
        return dz;
    };
    std::vector<T> z0(dim, zero);
    for (std::size_t i = 0; i < n; ++i) z0[i] = x0[i];
    for (std::size_t j = 0; j < nph; ++j) z0[n + j] = clk.phi0[j];
    OdeOptions<T> opt;
    opt.rtol = num_traits<T>::from_double(1e-8);
    opt.atol = num_traits<T>::from_double(1e-10);
    opt.store_trajectory = false;
    const std::vector<T> z1 = ode_rosenbrock4(rhs, t0, t1, z0, opt).final_state();
    CacheSojournResult<T> r;
    r.wtot = z1[dim - 1];
    if (num_traits<T>::to_double(r.wtot) > 0.0) {
        r.xbar.assign(n, zero);
        for (std::size_t i = 0; i < n; ++i) r.xbar[i] = T(z1[n + nph + i] / r.wtot);
    }
    return r;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_CACHE_SOJOURN_H
