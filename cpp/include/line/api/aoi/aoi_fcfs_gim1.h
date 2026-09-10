/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_FCFS_GIM1_H
#define LINE_API_AOI_FCFS_GIM1_H

/**
 * Mean Age of Information, its transform and the peak age of a GI/M/1 FCFS
 * queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_fcfs_gim1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_fcfs_gim1.java (identical).
 *
 *   sigma solves Y*(mu - mu sigma) = sigma in (0,1);  eta = mu (1 - sigma)
 *   E[D]     = 1/eta
 *   E[A]     = lambda E[Y^2]/2 + 1/mu + lambda (-Y*'(eta))/eta
 *   E[Apeak] = E[Y] + E[D]
 *   A*(s)    = (lambda/s) * ( T*(s) - Apeak*(s) )                  (Theorem 3)
 *   with T*(s) = eta/(s + eta) the exponential system time and
 *   Apeak*(s) = mu/(s + mu) * ( Y*(s) - s/(s + eta) * Y*(s + eta) )
 *
 * from Inoue et al. (2019). The stationary system time is exponential with
 * rate eta and independent of the next interarrival, which is why the
 * correlation term reduces to -Y*'(eta)/eta.
 *
 * static_assert(num_traits<T>::has_transcendental) -- a bracketed root of a
 * transcendental equation plus MATLAB's finite-difference derivative.
 *
 * Note that MATLAB uses the argument E_Y2 in the mean but never validates it
 * against Y_lst, so an inconsistent pair silently produces an inconsistent
 * answer; the port keeps that behaviour and only checks E_Y2 >= E_Y^2.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @param Y_lst LST of the interarrival time
 * @param mu    service rate, > 0
 * @param E_Y   mean interarrival time, > 0
 * @param E_Y2  second raw moment of the interarrival time, >= E_Y^2
 * @return      [meanAoI, peakAoI, A*(s)]
 */
template <class T>
AoiLstResult<T> aoi_fcfs_gim1(const Lst<T>& Y_lst, const T& mu, const T& E_Y, const T& E_Y2) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_fcfs_gim1 requires transcendental arithmetic");
    detail::require_positive(mu, "aoi_fcfs_gim1", "the service rate mu");
    detail::require_positive(E_Y, "aoi_fcfs_gim1", "the mean interarrival time E_Y");
    if (!Y_lst) throw InputError("aoi_fcfs_gim1: the interarrival LST must be callable");
    if (E_Y2 < E_Y * E_Y) throw InputError("aoi_fcfs_gim1: E_Y2 must be at least E_Y^2");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T lambda = one / E_Y;
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_fcfs_gim1");

    const T sigma = detail::gim1_sigma<T>(Y_lst, mu, "aoi_fcfs_gim1");
    const T E_D = one / (mu * (one - sigma));
    const T eta = mu * (one - sigma);
    const T dYstar = detail::lst_derivative<T>(Y_lst, eta);

    const T meanAoI = lambda * E_Y2 / two + one / mu + lambda * (-dYstar) / eta;
    const T peakAoI = E_Y + E_D;

    // LST of AoI (Inoue et al. 2019, Theorem 3), via the general age formula
    //   A*(s) = (lambda/s) * ( T*(s) - Apeak*(s) )
    // the cycle average of exp(-s*age) over a departure interval. In GI/M/1 the
    // system time is EXPONENTIAL at rate eta = mu*(1-sigma), so T*(s) =
    // eta/(s+eta), and Lindley gives W' + Y = max(Y, T) with T ~ Exp(eta)
    // independent of the next interarrival Y, so
    //   E[exp(-s*max(Y,T))] = Y*(s) - (s/(s+eta)) * Y*(s+eta),
    // and the peak adds one fresh Exp(mu) service. A*(0) = 1 follows from the
    // defining relation Y*(eta) = sigma.
    //
    // THE PREVIOUS FORM WAS NOT AN LST: (mu*sigma(s))/(s+mu-mu*sigma(s))*D*(s)
    // gives sigma/(1-sigma) at s = 0 rather than 1, and it re-solved sigma(s)
    // by bisection at every point. Checked against simulation on E2/M/1: at
    // s = 0.2 the old form gave 0.23056, the form below 0.59319, and the sample
    // path 0.59332.
    Lst<T> lstAoI = [Y_lst, mu, eta, lambda, one](const T& s) {
        if (num_traits<T>::to_double(s < T(0) ? T(-s) : s) < 1e-12)
            return one;  // A*(0) = 1 for any proper LST
        const T Ts = eta / (s + eta);
        const T peak = (mu / (s + mu)) * (Y_lst(s) - (s / (s + eta)) * Y_lst(T(s + eta)));
        return T((lambda / s) * (Ts - peak));
    };
    return {meanAoI, peakAoI, lstAoI, true};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_FCFS_GIM1_H
