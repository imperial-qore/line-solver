/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_LAP_H
#define LINE_API_PFQN_PFQN_LAP_H

/**
 * Laplace approximation of the normalizing constant of a repairman
 * (single-queue, multiclass) model.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_lap.m. The McKenna-Mitra
 * integral for a single queueing station is reduced to a one-dimensional
 * Laplace integral whose saddle point u0 solves
 *
 *   f(u) = 1 - sum_r N_r L_r / (Z_r + Ntot L_r u) = 0
 *
 * after which
 *
 *   log I = log Ntot - sum_r factln(N_r) - Ntot u0 + sum_r N_r log(Z_r + L_r u0 Ntot)
 *           + (1/2) log(2 pi) - (1/2) log(sum_r (N_r/Ntot)/(Z_r/(Ntot L_r) + u0)^2)
 *           - (1/2) log Ntot.
 *
 * ROOT FINDING. MATLAB calls fzero from the initial guess 1, and falls back to
 * a 1e-4 grid scan over (0,10] when fzero returns a non-finite root. f is
 * strictly increasing in u on u > 0 (each term N_r L_r/(Z_r + Ntot L_r u) is
 * decreasing), so the port brackets the root by doubling from 1 and then
 * bisects, which lands on the same root fzero converges to but without a
 * derivative or a Newton step that could leave the domain. The fallback scan
 * is kept for the case where no sign change exists on any bracket, exactly as
 * MATLAB's is.
 *
 * MATLAB returns NaN when the root is negative. The port throws instead:
 * a NaN normalizing constant propagates silently through a solver, whereas the
 * condition it signals (no admissible saddle point) is a modelling error.
 *
 * ARITHMETIC. Laplace's method plus logarithms, so gated on
 * num_traits<T>::has_transcendental.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Assembly of the expansion at the saddle point; defined below. */
template <class T>
T finish_lap(const std::vector<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
             const T& Ntot, const T& u0);

/**
 * @param L (R) per-class demand at the single station
 * @param N (R) per-class population
 * @param Z (R) per-class think time
 * @return  log of the approximate normalizing constant
 */
template <class T>
T pfqn_lap(const std::vector<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_lap requires transcendental arithmetic (Laplace approximation of an integral)");
    using std::log;
    const std::size_t R = L.size();
    if (N.size() != R || Z.size() != R)
        throw InputError(
            "pfqn_lap expects per-class vectors for a single queueing station (repairman models)");
    const T zero = num_traits<T>::from_int(0);
    T Ntot = zero;
    for (const T& v : N) Ntot += v;
    if (Ntot <= zero) throw InputError("pfqn_lap: empty population");

    // f(u) = 1 - sum_r N_r L_r / (Z_r + Ntot L_r u), increasing on u > 0.
    const auto f = [&](const T& u) {
        T s = num_traits<T>::from_int(1);
        for (std::size_t r = 0; r < R; ++r) {
            const T den = T(Z[r] + Ntot * L[r] * u);
            if (den == zero) continue;
            s -= T(N[r] * L[r] / den);
        }
        return s;
    };

    T lo = num_traits<T>::from_double(1e-12), hi = num_traits<T>::from_int(1);
    bool bracketed = false;
    if (f(lo) <= zero) {
        for (int k = 0; k < 200; ++k) {
            if (f(hi) >= zero) {
                bracketed = true;
                break;
            }
            lo = hi;
            hi = T(hi * num_traits<T>::from_int(2));
        }
    } else {
        // f already positive at the left edge: MATLAB's grid scan finds no
        // sign change either and u0 stays at the left edge.
        return finish_lap(L, N, Z, Ntot, lo);
    }
    if (!bracketed) throw NumericError("pfqn_lap: no saddle point on (0, 2^200]");
    for (int it = 0; it < 400; ++it) {
        const T mid = T(T(lo + hi) / num_traits<T>::from_int(2));
        if (f(mid) < zero)
            lo = mid;
        else
            hi = mid;
        if (num_traits<T>::to_double(num_abs(T(hi - lo))) <=
            1e-16 * (1.0 + num_traits<T>::to_double(hi)))
            break;
    }
    const T u0 = T(T(lo + hi) / num_traits<T>::from_int(2));
    return finish_lap(L, N, Z, Ntot, u0);
}

template <class T>
T finish_lap(const std::vector<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
             const T& Ntot, const T& u0) {
    using std::log;
    const std::size_t R = L.size();
    const T zero = num_traits<T>::from_int(0);
    if (u0 < zero) throw NumericError("pfqn_lap: negative saddle point, no admissible expansion");
    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);

    T logI = T(log(Ntot));
    for (std::size_t r = 0; r < R; ++r) logI -= detail::num_factln<T>(N[r]);
    logI -= Ntot * u0;
    for (std::size_t r = 0; r < R; ++r) logI += N[r] * log(T(Z[r] + L[r] * u0 * Ntot));
    T f2 = zero;
    for (std::size_t r = 0; r < R; ++r) {
        if (L[r] == zero) continue;
        const T d = T(T(Z[r] / T(Ntot * L[r])) + u0);
        f2 += T(T(N[r] / Ntot) / T(d * d));
    }
    if (f2 <= zero) throw NumericError("pfqn_lap: non-positive curvature at the saddle point");
    logI += T(num_traits<T>::from_rational(1, 2) * log(twopi));
    logI -= T(num_traits<T>::from_rational(1, 2) * log(f2));
    logI -= T(num_traits<T>::from_rational(1, 2) * log(Ntot));
    return logI;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_LAP_H
