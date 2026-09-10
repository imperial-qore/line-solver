/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_JOINT_DERIVATIVE_H
#define LINE_API_MAM_MAP_JOINT_DERIVATIVE_H

/**
 * Derivatives at the origin of a MAP's complementary CDF and of its joint
 * inter-arrival density.
 *
 * Templated port of matlab/src/api/mam/map_ccdf_derivative.m and
 * matlab/src/api/mam/map_jointpdf_derivative.m. Both are the building blocks
 * of the joint-moment analysis of networks of MAP/MAP/1 queues in
 * A. Horvath, G. Horvath, M. Telek, "A Joint Moments Based Analysis of
 * Networks of MAP/MAP/1 Queues".
 *
 * The CCDF of the inter-arrival time of a MAP is F^c(t) = pie exp(D0 t) e, so
 * its i-th derivative at 0 is
 *
 *     nu_i = pie D0^i e.
 *
 * The joint density of a run of consecutive inter-arrival times is
 * f(t_1,...,t_k) = pie exp(D0 t_1) D1 ... exp(D0 t_k) D1 e, so the mixed
 * partial derivative of orders (i_1,...,i_k) at the origin is
 *
 *     gamma = pie (D0^{i_1} D1) (D0^{i_2} D1) ... (D0^{i_k} D1) e.
 *
 * ARITHMETIC. Both are finite products of the descriptor matrices with the
 * embedded arrival vector pie, which is itself a linear solve, so the whole
 * computation is rational in the entries of D0 and D1. Neither is gated:
 * they instantiate at Rational, where the derivatives come out as exact
 * fractions. That matters here because these derivatives alternate in sign
 * and grow like i! ||D0||^i, so at double a moderately stiff MAP loses most
 * of its significant digits by i = 6, and only the exact instantiation can
 * tell a genuine near-cancellation from an accumulation of rounding.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Derivative of order i at 0 of the MAP's complementary CDF, nu = pie D0^i e.
 * Order 0 returns pie e = 1.
 */
template <class T>
T map_ccdf_derivative(const Map<T>& m, unsigned i) {
    const std::size_t n = m.order();
    if (n == 0) throw InputError("map_ccdf_derivative: empty MAP");
    const std::vector<T> pie = map_pie(m);
    const std::vector<T> v = vecmul(pie, matpow(m.D0, i));
    T s = num_traits<T>::from_int(0);
    for (const T& x : v) s += x;
    return s;
}

/**
 * Mixed partial derivative at the origin of the joint density of consecutive
 * inter-arrival times, gamma = pie prod_j (D0^{i_j} D1) e.
 *
 * An empty index set returns pie e = 1, matching the MATLAB loop over an
 * empty vector.
 */
template <class T>
T map_jointpdf_derivative(const Map<T>& m, const std::vector<unsigned>& iset) {
    const std::size_t n = m.order();
    if (n == 0) throw InputError("map_jointpdf_derivative: empty MAP");
    std::vector<T> g = map_pie(m);
    for (std::size_t k = 0; k < iset.size(); ++k) {
        g = vecmul(g, matpow(m.D0, iset[k]));
        g = vecmul(g, m.D1);
    }
    T s = num_traits<T>::from_int(0);
    for (const T& x : g) s += x;
    return s;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_JOINT_DERIVATIVE_H
