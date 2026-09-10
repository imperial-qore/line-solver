/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_GIG1_RQ_H
#define LINE_API_QSYS_QSYS_GIG1_RQ_H

/**
 * Robust Queueing (RQ) approximation of a G/GI/1 queue characterized by its
 * arrival index of dispersion and the first two service moments.
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_rq.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gig1_rq.java.
 *
 * The mean steady-state workload is the value of a one-dimensional variational
 * problem (Whitt and You 2018, eqs. (13), (16)-(18)):
 *
 *   Zstar = sup_{x>=0} [ -(1-rho) x + sqrt( 2 rho x (I_a(x) + c_s^2)/mu ) ]
 *   W  = max(0, Zstar/rho - (c_s^2+1)/(2 mu))
 *   Q  = lambda W          the mean number waiting
 *   X  = Q + rho           the mean number in system
 *
 * The objective is unimodal in practice but not guaranteed to be, so MATLAB
 * brackets it with a 200-point log-spaced scan over [1e-6, 1e8] and refines
 * the best bracket with fminbnd at TolX 1e-10, keeping the better of the scan
 * value and the refined value. The port keeps the same scan and refines with
 * golden-section search on the same bracket at the same TolX. fminbnd is
 * golden section with parabolic acceleration, so on a bracket containing a
 * single interior maximum the two locate the same point to within TolX; and
 * because MATLAB and the port both return max(scan, refined), a refinement
 * that lands short can never fall below the scan value. Where the objective
 * really is multimodal both are equally at the mercy of the scan, and neither
 * claims a global optimum.
 *
 * ARITHMETIC. The square root and the tolerance-driven search make this
 * transcendental.
 *
 * The arrival process enters only through the callable I_a(x), so the caller
 * supplies whatever index-of-dispersion model applies; for a renewal arrival
 * stream I_a is the constant c_a^2.
 *
 * At rho <= 0 all four measures are zero, as in MATLAB. At rho >= 1 MATLAB
 * returns Inf; the port raises instead, since the exact instantiations have no
 * infinity and a silent Inf propagates into whatever consumes the result.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

template <class T>
struct Gig1RqResult {
    T Z;  ///< mean steady-state workload E[Z]
    T W;  ///< mean steady-state waiting time E[W]
    T Q;  ///< mean number waiting, lambda W
    T X;  ///< mean number in system, Q + rho
};

/**
 * @param rho   traffic intensity lambda/mu
 * @param mu    service rate
 * @param cs2   squared coefficient of variation of the service time
 * @param IaFun_ callable, IaFun_(x) -> the arrival IDC I_a(x) at x > 0
 */
template <class T, class IaFun>
Gig1RqResult<T> qsys_gig1_rq(const T& rho, const T& mu, const T& cs2, IaFun&& IaFun_) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gig1_rq requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    Gig1RqResult<T> r;
    if (rho <= zero) {
        r.Z = r.W = r.Q = r.X = zero;
        return r;
    }
    if (rho >= one)
        throw InputError("qsys_gig1_rq: rho must be strictly less than 1 for a finite workload");
    const T lambda = rho * mu;

    // f(x) = -(1-rho) x + sqrt(2 rho x (I_a(x) + c_s^2)/mu), and f(x<=0) = 0.
    auto f = [&](const T& x) -> T {
        if (x <= zero) return zero;
        const T ia = IaFun_(x);
        T inner = two * rho * x * (ia + cs2) / mu;
        if (inner < zero) inner = zero;
        return -(one - rho) * x + detail::num_sqrt(inner);
    };

    // Coarse log-spaced scan over [1e-6, 1e8], 200 points, as in MATLAB.
    const std::size_t NS = 200;
    std::vector<T> xs(NS);
    for (std::size_t i = 0; i < NS; ++i) {
        const double e = -6.0 + 14.0 * static_cast<double>(i) / static_cast<double>(NS - 1);
        xs[i] = T(num_traits<T>::from_double(std::pow(10.0, e)));
    }
    std::size_t imax = 0;
    T best = f(xs[0]);
    for (std::size_t i = 1; i < NS; ++i) {
        const T v = f(xs[i]);
        if (v > best) {
            best = v;
            imax = i;
        }
    }
    T lo = xs[imax > 0 ? imax - 1 : 0];
    T hi = xs[imax + 1 < NS ? imax + 1 : NS - 1];

    // Golden-section refinement on the bracket at TolX 1e-10.
    const T tolx = T(num_traits<T>::from_double(1e-10));
    const T invphi = T(num_traits<T>::from_double(0.6180339887498949));
    T c = hi - (hi - lo) * invphi;
    T d = lo + (hi - lo) * invphi;
    T fc = f(c), fd = f(d);
    for (unsigned it = 0; it < 500u && hi - lo > tolx; ++it) {
        if (fc > fd) {
            hi = d;
            d = c;
            fd = fc;
            c = hi - (hi - lo) * invphi;
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + (hi - lo) * invphi;
            fd = f(d);
        }
    }
    const T refined = fc > fd ? fc : fd;
    T Z = best > refined ? best : refined;
    if (Z < zero) Z = zero;

    r.Z = Z;
    T W = Z / rho - (cs2 + one) / (two * mu);
    if (W < zero) W = zero;
    r.W = W;
    r.Q = lambda * W;
    r.X = r.Q + rho;
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_GIG1_RQ_H
