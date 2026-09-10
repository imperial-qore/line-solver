/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GIG1_APPROX_WHITT_H
#define LINE_API_QSYS_GIG1_APPROX_WHITT_H

/**
 * Whitt's approximation of the G/G/1 mean response time.
 *
 * Port of `qsys_gig1_approx_whitt` in
 * python/line_solver/api/qsys/approximations.py. PYTHON-ONLY, and its own
 * docstring says so: there is no `qsys_gig1_approx_whitt.m`.
 *
 * It is the Kingman diffusion form `Lq = rho^2 (ca^2 + cs^2) / (2 (1 - rho))`
 * multiplied by a CORRECTION FACTOR phi. The correction is what distinguishes
 * it: Kingman's form is an upper bound that is loose when the arrival stream is
 * more regular than Poisson, and phi discounts it exactly there.
 *
 * PHI IS PIECEWISE AND ONLY TWO OF ITS FOUR ARMS DO ANYTHING. With both
 * variability parameters at or below one -- the regular regime where Kingman is
 * loosest -- phi is an exponential discount in `(1 - ca^2)^2`. With a bursty
 * arrival stream and regular service it is a milder discount. In the two arms
 * where the SERVICE is bursty, `cs^2 > 1`, phi is exactly one and the
 * approximation falls back to Kingman: the correction has nothing to offer
 * there, and pretending otherwise would be an invented formula.
 *
 * An unstable queue returns an infinite response time and a unit utilization
 * rather than dividing by a non-positive `1 - rho`.
 *
 * ARITHMETIC: transcendental (phi is an exponential).
 */

#include <cmath>
#include <limits>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"

namespace line {
namespace qsys {

/**
 * @param lambda arrival rate
 * @param mu     service rate
 * @param ca     coefficient of variation of the interarrival time
 * @param cs     coefficient of variation of the service time
 */
template <class T>
QsysResult<T> qsys_gig1_approx_whitt(const T& lambda, const T& mu, const T& ca, const T& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gig1_approx_whitt needs exp() for the correction factor");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T four = num_traits<T>::from_int(4);
    const T rho = lambda / mu;
    if (!(num_traits<T>::to_double(one - rho) > 0.0)) {
        // Unstable: the queue has no stationary response time, and the
        // reference reports that rather than dividing.
        return {num_traits<T>::from_double(std::numeric_limits<double>::infinity()), one};
    }

    const T ca2 = ca * ca, cs2 = cs * cs;
    const double ca2d = num_traits<T>::to_double(ca2), cs2d = num_traits<T>::to_double(cs2);
    T phi = one;
    if (ca2d <= 1.0 && cs2d <= 1.0) {
        const T e = T(-two * (one - rho) * (one - ca2) * (one - ca2) / (three * rho * (ca2 + cs2)));
        phi = num_traits<T>::from_double(std::exp(num_traits<T>::to_double(e)));
    } else if (ca2d > 1.0 && cs2d <= 1.0) {
        const T e = T(-(one - rho) * (ca2 - one) / (ca2 + four * cs2));
        phi = num_traits<T>::from_double(std::exp(num_traits<T>::to_double(e)));
    }
    // Both remaining arms leave phi at one: with bursty service the correction
    // has nothing to add and this IS Kingman.

    const T Lq = T(phi * rho * rho * (ca2 + cs2) / (two * (one - rho)));
    const T L = T(Lq + rho);
    const T W = T(L / lambda);
    return {W, detail::rhohat_from_W(W, lambda)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GIG1_APPROX_WHITT_H
