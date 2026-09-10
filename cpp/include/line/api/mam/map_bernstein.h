/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_BERNSTEIN_H
#define LINE_API_MAM_MAP_BERNSTEIN_H

/**
 * Acyclic phase-type approximation of an arbitrary density by Bernstein
 * exponentials.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_bernstein.m and
 * jar/src/main/java/jline/api/mam/Aph_bernstein.java, after
 * Horvath and Vicario, "Construction of Phase Type Distributions by Bernstein
 * Exponentials", EPEW 2023, 201-215.
 *
 * The construction evaluates the target density at the n Bernstein nodes
 * x_i = -log(i/n) and uses the values as the entry law of the Erlang cascade
 *   T = diag(-1..-n) + superdiag(1..n-1),
 * whose i-th phase has an Erlang(i,i) residual time. The result is a renewal
 * process, D1 = -T 1 alpha, so the fit is a PH and not a general MAP; the caller
 * rescales it to the target mean with map_scale, since the construction pins the
 * time unit at one.
 *
 * WHICH REFERENCE. The MATLAB and JAR versions differ and MATLAB is the ground
 * truth here: it skips a node where the density is not finite and positive,
 * renormalizes alpha (which the bare 1/(i c) weights only satisfy up to the
 * skipped nodes), and falls back to an Erlang-n of unit mean when the
 * normalization constant is not usable at all. The JAR divides by c
 * unconditionally, so a density that underflows at some node -- a Pareto below
 * its scale, a Uniform outside its support, any bounded law -- gives it a NaN
 * generator instead of a fit. That path is exactly the one `sn_nonmarkov_toph`
 * takes for Uniform and Pareto service.
 *
 * ARITHMETIC: transcendental, since the nodes are logarithms of i/n.
 */

#include <cmath>
#include <cstddef>
#include <functional>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * @param f     the density to approximate, evaluated at positive abscissae
 * @param order number of phases, the reference default being 20
 * @return the fitted renewal MAP, of unit time scale; rescale with map_scale
 */
template <class T>
Map<T> map_bernstein(const std::function<double(double)>& f, unsigned order = 20) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_bernstein evaluates the density at logarithmic nodes");
    if (order == 0) throw InputError("map_bernstein: the order must be positive");
    const std::size_t n = order;
    const T zero = num_traits<T>::from_int(0);

    std::vector<double> fv(n, 0.0);
    double c = 0.0;
    for (std::size_t i = 1; i <= n; ++i) {
        const double xi = -std::log(static_cast<double>(i) / static_cast<double>(n));
        const double fi = f(xi);
        if (std::isfinite(fi) && fi > 0.0) {
            fv[i - 1] = fi;
            c += fi / static_cast<double>(i);
        }
    }
    // No usable mass at any node: the reference returns a unit-mean Erlang-n and
    // leaves the rescaling to the caller, rather than emitting a NaN generator.
    if (!(c > 0.0) || !std::isfinite(c))
        return map_erlang(num_traits<T>::from_int(1), static_cast<unsigned>(n));

    std::vector<double> alpha(n, 0.0);
    double asum = 0.0;
    for (std::size_t i = 1; i <= n; ++i) {
        if (fv[i - 1] > 0.0) {
            alpha[i - 1] = fv[i - 1] / (static_cast<double>(i) * c);
            asum += alpha[i - 1];
        }
    }
    if (asum > 0.0) {
        for (std::size_t i = 0; i < n; ++i) alpha[i] /= asum;
    } else {
        alpha[0] = 1.0;  // degenerate: start in the first phase
    }

    Map<T> m;
    m.D0 = Matrix<T>(n, n, zero);
    m.D1 = Matrix<T>(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        m.D0(i, i) = num_traits<T>::from_int(-static_cast<int>(i + 1));
        if (i + 1 < n) m.D0(i, i + 1) = num_traits<T>::from_int(static_cast<int>(i + 1));
    }
    // D1 = -T P with P = 1 alpha: every phase restarts from the entry law.
    for (std::size_t i = 0; i < n; ++i) {
        T rowsum = zero;
        for (std::size_t j = 0; j < n; ++j) rowsum += m.D0(i, j);
        for (std::size_t j = 0; j < n; ++j)
            m.D1(i, j) = T(-rowsum * num_traits<T>::from_double(alpha[j]));
    }
    return m;
}

/**
 * The reference's (D0, D1) spelling of the same fit, kept because the JAR
 * exposes it under this name (Aph_bernstein) and callers ported from it expect
 * the pair rather than a Map.
 */
template <class T>
Map<T> aph_bernstein(const std::function<double(double)>& f, unsigned order = 20) {
    return map_bernstein<T>(f, order);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_BERNSTEIN_H
