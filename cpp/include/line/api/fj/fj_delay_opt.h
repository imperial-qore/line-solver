/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_DELAY_OPT_H
#define LINE_API_FJ_DELAY_OPT_H

/**
 * Deterministic subtask delays that minimise mean dispersion.
 *
 * Templated port of matlab/src/api/fj/fj_delay_opt.m.
 *
 * Holding back a fast branch costs little at the last completion and buys a
 * great deal at the first, so the minimiser of the dispersion of fj_dispersion
 * is generally interior and strictly positive on every branch but the slowest.
 *
 * The objective is minimised by cyclic coordinate descent with a golden section
 * line search on each coordinate: deterministic, derivative-free, and the same
 * sequence of evaluations in all four codebases. Adding a constant to every
 * delay shifts both order statistics equally, so the search returns the
 * representative with min(d) = 0.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_dispersion.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [d, Edisp, Emax] of fj_delay_opt. */
template <class T>
struct FJDelayOptResult {
    std::vector<T> d;
    T Edisp;
    T Emax;
};

/**
 * @param shape     Erlang stage counts, one per branch
 * @param rate      Erlang stage rates, one per branch
 * @param maxsweeps coordinate-descent sweep cap
 * @param dtol      relative convergence tolerance
 * @param npanels   Simpson panel count passed to fj_dispersion
 * @return          the optimal delays and the dispersion and last-completion mean there
 */
template <class T>
FJDelayOptResult<T> fj_delay_opt(const std::vector<unsigned>& shape, const std::vector<T>& rate,
                                 unsigned maxsweeps = 40,
                                 const T& dtol = num_traits<T>::from_double(1e-8),
                                 unsigned npanels = 2000) {
    const std::size_t N = shape.size();
    if (rate.size() != N)
        throw InputError("fj_delay_opt: shape and rate must have the same length");
    if (N < 1) throw InputError("fj_delay_opt: at least one branch is required");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T tol = num_traits<T>::from_double(1e-10);

    FJDelayOptResult<T> out;
    out.d.assign(N, zero);
    if (N < 2) {
        const FJDispersionResult<T> r = fj_dispersion<T>(shape, rate, out.d, tol, npanels);
        out.Edisp = r.Edisp;
        out.Emax = r.Emax;
        return out;
    }

    // Delaying past the slowest branch never helps, so that is the search ceiling
    T mmax = zero, smax = zero;
    for (std::size_t i = 0; i < N; ++i) {
        const T mi = num_traits<T>::from_int(static_cast<long>(shape[i])) / rate[i];
        const T si = detail::num_sqrt<T>(num_traits<T>::from_int(static_cast<long>(shape[i]))) / rate[i];
        if (mi > mmax) mmax = mi;
        if (si > smax) smax = si;
    }
    const T ub = mmax + num_traits<T>::from_int(8) * smax;

    const T invphi = (detail::num_sqrt<T>(num_traits<T>::from_int(5)) - one) /
                     num_traits<T>::from_int(2);
    T fcur = fj_dispersion<T>(shape, rate, out.d, tol, npanels).Edisp;
    for (unsigned sweep = 0; sweep < maxsweeps; ++sweep) {
        const T fprev = fcur;
        for (std::size_t i = 0; i < N; ++i) {
            // Golden section on coordinate i, the other delays held fixed
            T a = zero, b = ub;
            T c = b - invphi * (b - a), dd = a + invphi * (b - a);
            std::vector<T> probe = out.d;
            probe[i] = c;
            T fc = fj_dispersion<T>(shape, rate, probe, tol, npanels).Edisp;
            probe[i] = dd;
            T fd = fj_dispersion<T>(shape, rate, probe, tol, npanels).Edisp;
            for (unsigned it = 0; it < 60; ++it) {
                if (fc < fd) {
                    b = dd; dd = c; fd = fc;
                    c = b - invphi * (b - a);
                    probe[i] = c;
                    fc = fj_dispersion<T>(shape, rate, probe, tol, npanels).Edisp;
                } else {
                    a = c; c = dd; fc = fd;
                    dd = a + invphi * (b - a);
                    probe[i] = dd;
                    fd = fj_dispersion<T>(shape, rate, probe, tol, npanels).Edisp;
                }
                if ((b - a) <= dtol * (ub > one ? ub : one)) break;
            }
            out.d[i] = (fc < fd) ? c : dd;
        }
        // Normalise so that the smallest delay is zero
        T dmin = out.d[0];
        for (std::size_t i = 1; i < N; ++i)
            if (out.d[i] < dmin) dmin = out.d[i];
        for (std::size_t i = 0; i < N; ++i) out.d[i] -= dmin;
        fcur = fj_dispersion<T>(shape, rate, out.d, tol, npanels).Edisp;
        const T gap = (fprev > fcur) ? (fprev - fcur) : (fcur - fprev);
        const T mag = (fprev > one || fprev < -one) ? (fprev > zero ? fprev : -fprev) : one;
        if (gap <= dtol * mag) break;
    }

    const FJDispersionResult<T> r = fj_dispersion<T>(shape, rate, out.d, tol, npanels);
    out.Edisp = r.Edisp;
    out.Emax = r.Emax;
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_DELAY_OPT_H
