/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_MAP_INTERP_H
#define LINE_API_FES_MAP_INTERP_H

/**
 * Shape-preserving interpolation of the flow-equivalent descriptors, and the
 * population grid they are evaluated on.
 *
 * Templated port of matlab/src/api/fes/fes_map_interp.m, fes_map_interp_edge.m and
 * fes_map_grid.m, mirrored by the JAR and native Python.
 *
 * Fritsch and Carlson slopes are used, with the noncentered three-point endpoint
 * rule of de Boor, so the interpolant never overshoots and a monotone sequence of
 * throughputs stays monotone. The algorithm is written out rather than delegated
 * to MATLAB's pchip so that the four codebases return identical values; the
 * MATLAB port is checked against the built-in to machine precision.
 *
 * ARITHMETIC: field operations only, exact at T = Rational.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace fes {

/** Leading populations kept in full by the default grid. */
const std::size_t FES_MAP_GRID_NHEAD = 10;
/** Equispaced points above them. */
const std::size_t FES_MAP_GRID_NTAIL = 10;

/**
 * Populations at which the inter-departure MAP is evaluated. Fitting one MAP per
 * population is wasteful because the processes of neighbouring populations are
 * similar; the reference evaluates the first ten populations and ten further
 * equispaced points.
 */
inline std::vector<std::size_t> fes_map_grid(std::size_t n, std::size_t nhead = FES_MAP_GRID_NHEAD,
                                             std::size_t ntail = FES_MAP_GRID_NTAIL) {
    std::vector<std::size_t> grid;
    if (n <= nhead + ntail) {
        for (std::size_t k = 1; k <= n; ++k) grid.push_back(k);
        return grid;
    }
    std::vector<bool> taken(n + 1, false);
    for (std::size_t k = 1; k <= nhead; ++k) taken[k] = true;
    for (std::size_t j = 0; j < ntail; ++j) {
        const double t = (ntail == 1) ? static_cast<double>(n)
                                      : static_cast<double>(nhead + 1) +
                                            static_cast<double>(j) *
                                                static_cast<double>(n - nhead - 1) /
                                                static_cast<double>(ntail - 1);
        taken[static_cast<std::size_t>(std::llround(t))] = true;
    }
    for (std::size_t k = 1; k <= n; ++k)
        if (taken[k]) grid.push_back(k);
    return grid;
}

namespace detail {

/** Noncentered three-point endpoint slope with the monotonicity clamps of de Boor. */
template <class T>
T interp_edge(const T& h1, const T& h2, const T& del1, const T& del2) {
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T d = T(((two * h1 + h2) * del1 - h1 * del2) / (h1 + h2));
    const int sd = (d > zero) - (d < zero);
    const int s1 = (del1 > zero) - (del1 < zero);
    const int s2 = (del2 > zero) - (del2 < zero);
    if (sd != s1) return zero;
    const T ad = (d < zero) ? T(-d) : d;
    const T a1 = (del1 < zero) ? T(-(three * del1)) : T(three * del1);
    if (s1 != s2 && ad > a1) return three * del1;
    return d;
}

}  // namespace detail

/**
 * Monotone piecewise cubic Hermite interpolation of one series.
 *
 * @param x  sample abscissae, strictly increasing
 * @param y  sample values
 * @param xq query abscissae
 */
template <class T>
std::vector<T> fes_map_interp(const std::vector<T>& x, const std::vector<T>& y,
                              const std::vector<T>& xq) {
    const std::size_t n = x.size();
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    std::vector<T> yq(xq.size(), zero);
    if (n == 1) {
        for (std::size_t q = 0; q < xq.size(); ++q) yq[q] = y[0];
        return yq;
    }

    std::vector<T> h(n - 1), delta(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = x[i + 1] - x[i];
        delta[i] = T((y[i + 1] - y[i]) / h[i]);
    }

    std::vector<T> d(n, zero);
    if (n == 2) {
        d[0] = delta[0];
        d[1] = delta[0];
    } else {
        for (std::size_t i = 1; i + 1 < n; ++i) {
            if (delta[i - 1] * delta[i] > zero) {
                const T w1 = two * h[i] + h[i - 1];
                const T w2 = h[i] + two * h[i - 1];
                d[i] = T((w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]));
            }
        }
        d[0] = detail::interp_edge(h[0], h[1], delta[0], delta[1]);
        d[n - 1] = detail::interp_edge(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
    }

    for (std::size_t q = 0; q < xq.size(); ++q) {
        const T t = xq[q];
        std::size_t i;
        if (!(t > x[0])) {
            i = 0;
        } else if (!(t < x[n - 1])) {
            i = n - 2;
        } else {
            i = 0;
            while (i + 2 < n && !(x[i + 1] > t)) ++i;
        }
        const T s = t - x[i];
        const T c2 = T((three * delta[i] - two * d[i] - d[i + 1]) / h[i]);
        const T c3 = T((d[i] - two * delta[i] + d[i + 1]) / (h[i] * h[i]));
        yq[q] = y[i] + s * (d[i] + s * (c2 + s * c3));
    }
    return yq;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_MAP_INTERP_H
