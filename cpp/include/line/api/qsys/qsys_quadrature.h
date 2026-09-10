/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QUADRATURE_H
#define LINE_API_QSYS_QUADRATURE_H

/**
 * Adaptive quadrature for the qsys functions whose MATLAB originals call
 * integral(), and the cumulative trapezoid rule for the one that calls
 * cumtrapz().
 *
 * Three MATLAB files in matlab/src/api/qsys reduce a per-class mean response
 * time to a definite integral of a non-elementary integrand -- qsys_mg1_fb,
 * qsys_mg1_psjf and qsys_mg1k_loss -- and one, qsys_mg1_srpt, integrates on a
 * fixed uniform grid with cumtrapz/trapz. This header carries the two rules
 * they need so that each ported function stays a 1:1 image of its MATLAB file.
 * It is the numerical counterpart of qsys_types.h, which carries the shared
 * return type; nothing here corresponds to a MATLAB file of its own.
 *
 * The adaptive rule is the Gauss-Kronrod 7/15 pair on recursively bisected
 * subintervals, with the local error estimated as |K - G| and the tolerance
 * split evenly between halves. That is the same rule MATLAB's integral() uses
 * (MATLAB applies it to a transformed interval and controls the error
 * globally), so the two agree to the requested relative tolerance on the
 * smooth integrands here; the ported functions therefore claim agreement at
 * the tolerance MATLAB was asked for, not beyond it.
 *
 * ARITHMETIC. Both rules are inherently inexact -- the Kronrod nodes are
 * irrational and the trapezoid rule has a discretization error -- so both are
 * gated on num_traits<T>::has_transcendental. The node and weight constants
 * are the QUADPACK values, but they enter through num_traits<T>::from_double
 * and are therefore carried at double precision: a Real<D> instantiation gains
 * exact accumulation and no cancellation in the sums, but the quadrature error
 * floor stays near 1e-16 relative because the nodes themselves do. Any qsys
 * function that routes through this header inherits that floor, and its tests
 * assert at the tolerance MATLAB's integral() was asked for, never below it.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {
namespace detail {

/** Kronrod 15-point abscissae on [-1,1], positive half plus the origin. */
template <class T>
const std::vector<T>& gk15_nodes() {
    static const std::vector<T> x = {
        num_traits<T>::from_double(0.991455371120812639206854697526329),
        num_traits<T>::from_double(0.949107912342758524526189684047851),
        num_traits<T>::from_double(0.864864423359769072789712788640926),
        num_traits<T>::from_double(0.741531185599394439863864773280788),
        num_traits<T>::from_double(0.586087235467691130294144838258730),
        num_traits<T>::from_double(0.405845151377397166906606412076961),
        num_traits<T>::from_double(0.207784955007898467600689403773245),
        num_traits<T>::from_double(0.0)};
    return x;
}

/** Kronrod weights matched to gk15_nodes. */
template <class T>
const std::vector<T>& gk15_weights() {
    static const std::vector<T> w = {
        num_traits<T>::from_double(0.022935322010529224963732008058970),
        num_traits<T>::from_double(0.063092092629978553290700663189204),
        num_traits<T>::from_double(0.104790010322250183839876322541518),
        num_traits<T>::from_double(0.140653259715525918745189590510238),
        num_traits<T>::from_double(0.169004726639267902826583426598550),
        num_traits<T>::from_double(0.190350578064785409913256402421014),
        num_traits<T>::from_double(0.204432940075298892414161999234649),
        num_traits<T>::from_double(0.209482141084727828012999174891714)};
    return w;
}

/** Gauss 7-point weights, applied at the odd-indexed Kronrod nodes. */
template <class T>
const std::vector<T>& g7_weights() {
    static const std::vector<T> w = {
        num_traits<T>::from_double(0.129484966168869693270611432679082),
        num_traits<T>::from_double(0.279705391489276667901467771423780),
        num_traits<T>::from_double(0.381830050505118944950369775488975),
        num_traits<T>::from_double(0.417959183673469387755102040816327)};
    return w;
}

/**
 * One Gauss-Kronrod 7/15 panel on [a,b]. Returns the Kronrod estimate and
 * writes |K - G| into err.
 */
template <class T, class F>
T gk15_panel(F&& f, const T& a, const T& b, T& err) {
    const T two = num_traits<T>::from_int(2);
    const T c = (a + b) / two;
    const T h = (b - a) / two;
    const std::vector<T>& x = gk15_nodes<T>();
    const std::vector<T>& wk = gk15_weights<T>();
    const std::vector<T>& wg = g7_weights<T>();

    T K = num_traits<T>::from_int(0);
    T G = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < 7; ++i) {
        const T d = h * x[i];
        const T fsum = f(T(c - d)) + f(T(c + d));
        K += wk[i] * fsum;
        if (i % 2 == 1) G += wg[i / 2] * fsum;
    }
    const T f0 = f(c);
    K += wk[7] * f0;
    G += wg[3] * f0;
    K *= h;
    G *= h;
    err = num_abs(T(K - G));
    return K;
}

/**
 * One adaptively bisected panel. min_depth forces the first few bisections
 * unconditionally: on a long interval where the integrand is supported on a
 * small part of it -- qsys_mg1k_loss integrates over [0, 1e4/lambda] a density
 * that has died out by t = O(1/mu) -- a single Kronrod panel can see fifteen
 * near-zero samples, estimate a near-zero error and accept a wrong answer.
 * MATLAB's integral() avoids this by subdividing the transformed interval
 * before it adapts; forcing a few levels here is the same guard.
 */
template <class T, class F>
T gk15_adapt(F&& f, const T& a, const T& b, const T& reltol, const T& abstol, unsigned depth,
             unsigned min_depth) {
    T err = num_traits<T>::from_int(0);
    const T q = gk15_panel(f, a, b, err);
    const T target = reltol * num_abs(q);
    if (min_depth == 0 && (depth == 0 || err <= (target > abstol ? target : abstol))) return q;
    if (depth == 0) return q;
    const T m = (a + b) / num_traits<T>::from_int(2);
    const T half = abstol / num_traits<T>::from_int(2);
    const unsigned md = min_depth == 0 ? 0u : min_depth - 1;
    return gk15_adapt(f, a, m, reltol, half, depth - 1, md) +
           gk15_adapt(f, m, b, reltol, half, depth - 1, md);
}

/**
 * Adaptive Gauss-Kronrod integration of f over [a,b].
 *
 * The interval is first split into init_panels equal pieces, each of which is
 * then bisected adaptively with at least min_depth forced levels.
 *
 * @param reltol relative tolerance, MATLAB's 'RelTol'
 * @param abstol absolute tolerance, MATLAB's 'AbsTol'
 * @param depth  bisection budget per panel; exhausting it returns the best
 *               estimate reached, exactly as MATLAB warns and returns rather
 *               than failing
 * @param f the integrand
 * @param a lower limit
 * @param b upper limit
 * @param init_panels number of panels of the first pass
 * @param min_depth refinement levels performed before the tolerance test
 */
template <class T, class F>
T num_integral(F&& f, const T& a, const T& b, const T& reltol, const T& abstol,
               unsigned depth = 50u, unsigned init_panels = 8u, unsigned min_depth = 3u) {
    static_assert(num_traits<T>::has_transcendental,
                  "num_integral requires transcendental arithmetic");
    if (b == a) return num_traits<T>::from_int(0);
    if (init_panels == 0) init_panels = 1;
    const T np = num_traits<T>::from_int(static_cast<long>(init_panels));
    const T panel_abstol = abstol / np;
    T total = num_traits<T>::from_int(0);
    for (unsigned i = 0; i < init_panels; ++i) {
        const T lo = a + (b - a) * num_traits<T>::from_int(static_cast<long>(i)) / np;
        const T hi = a + (b - a) * num_traits<T>::from_int(static_cast<long>(i + 1)) / np;
        total += gk15_adapt(f, lo, hi, reltol, panel_abstol, depth, min_depth);
    }
    return total;
}

/**
 * Trapezoid rule on the sample pairs (x, y), MATLAB's trapz(x, y).
 */
template <class T>
T num_trapz(const std::vector<T>& x, const std::vector<T>& y) {
    static_assert(num_traits<T>::has_transcendental,
                  "num_trapz requires transcendental arithmetic");
    if (x.size() != y.size()) throw InputError("num_trapz: x and y have different lengths");
    const T two = num_traits<T>::from_int(2);
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 1; i < x.size(); ++i) s += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / two;
    return s;
}

/**
 * Cumulative trapezoid rule, MATLAB's cumtrapz(x, y): element i is the
 * integral of y from x[0] to x[i], so the first element is zero.
 */
template <class T>
std::vector<T> num_cumtrapz(const std::vector<T>& x, const std::vector<T>& y) {
    static_assert(num_traits<T>::has_transcendental,
                  "num_cumtrapz requires transcendental arithmetic");
    if (x.size() != y.size()) throw InputError("num_cumtrapz: x and y have different lengths");
    const T two = num_traits<T>::from_int(2);
    std::vector<T> c(x.size(), num_traits<T>::from_int(0));
    for (std::size_t i = 1; i < x.size(); ++i)
        c[i] = c[i - 1] + (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / two;
    return c;
}

}  // namespace detail
}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QUADRATURE_H
