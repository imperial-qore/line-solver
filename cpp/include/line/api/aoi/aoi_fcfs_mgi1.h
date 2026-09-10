/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_FCFS_MGI1_H
#define LINE_API_AOI_FCFS_MGI1_H

/**
 * Mean Age of Information, its transform and the peak age of an M/GI/1 FCFS
 * queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_fcfs_mgi1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_fcfs_mgi1.java (identical).
 *
 *   W*(s)    = (1-rho) s / (s - lambda + lambda H*(s))       (Pollaczek-Khinchine)
 *   T*(s)    = H*(s) W*(s)
 *   A*(s)    = lambda H*(s) / (s + lambda - lambda H*(s)) * W*(s)   (Theorem 2)
 *   E[A]     = E[H] + E[T] + (1 - 2 rho)/lambda - (d/ds) T*(s) at s = lambda
 *   E[Apeak] = E[T] + 1/lambda
 *
 * from Inoue, Masuyama, Takine and Tanaka (IEEE Trans. IT 65(12), 2019). The
 * derivative term is what carries the correlation between the interarrival
 * time and the system time.
 *
 * static_assert(num_traits<T>::has_transcendental) -- the derivative is taken
 * by MATLAB's central difference with step 1e-6 max(1,lambda), reproduced here
 * so the two agree digit for digit. That step is a tolerance, not an exact
 * operation: in an exact field the function would return an exact value of the
 * wrong quantity, which is worse than refusing.
 *
 * Note that the transform A*(s) itself is a rational function of s whenever
 * H*(s) is, so the returned Lst is exact for its argument; only the mean is
 * approximate.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @param lambda arrival rate, > 0
 * @param H_lst  LST of the service time
 * @param E_H    mean service time, > 0
 * @param E_H2   second raw moment of the service time, >= E_H^2
 * @return       [meanAoI, peakAoI, A*(s)]
 */
template <class T>
AoiLstResult<T> aoi_fcfs_mgi1(const T& lambda, const Lst<T>& H_lst, const T& E_H, const T& E_H2) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_fcfs_mgi1 requires transcendental arithmetic");
    detail::require_positive(lambda, "aoi_fcfs_mgi1", "the arrival rate lambda");
    detail::require_positive(E_H, "aoi_fcfs_mgi1", "the mean service time E_H");
    if (!H_lst) throw InputError("aoi_fcfs_mgi1: the service-time LST must be callable");
    if (E_H2 < E_H * E_H) throw InputError("aoi_fcfs_mgi1: E_H2 must be at least E_H^2");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T rho = lambda * E_H;
    detail::require_stable(rho, "aoi_fcfs_mgi1");

    const T E_Y = one / lambda;
    const T E_W = lambda * E_H2 / (two * (one - rho));
    const T E_T = E_W + E_H;

    const Lst<T> Tstar = [lambda, H_lst, rho, one](const T& s) {
        if (num_traits<T>::to_double(s < T(0) ? T(-s) : s) < 1e-12)
            return one;  // removable 0/0 at the origin
        const T Hs = H_lst(s);
        return T(Hs * ((one - rho) * s / (s - lambda + lambda * Hs)));
    };
    const T dTstar = detail::lst_derivative<T>(Tstar, lambda);
    const T meanAoI = E_H + E_T + (one - two * rho) / lambda - dTstar;
    const T peakAoI = E_T + E_Y;

    // LST of AoI (Inoue et al. 2019, Theorem 2), via the general age formula
    //   A*(s) = (lambda/s) * ( T*(s) - Apeak*(s) )
    // the cycle average of exp(-s*age) over a departure interval: the age starts
    // each cycle at the system time T of the packet just delivered and grows to
    // the peak T + (next interarrival) at the next delivery. For M/GI/1 FCFS
    // Lindley gives W' = max(0, T - Y), so W' + Y is max(Y, T), and with
    // Y ~ Exp(lambda) independent of T,
    //   E[exp(-s*max(Y,T))] = T*(s) - (s/(s+lambda)) * T*(s+lambda).
    //
    // THE PREVIOUS FORM WAS NOT AN LST: (lambda*H*(s))/(s+lambda-lambda*H*(s))
    // diverges as s -> 0, so A*(0) was +Inf instead of 1 and the value exceeded
    // 1 for small s. Checked against simulation on M/E2/1: at s = 0.3 the old
    // form gave 1.3062, the form below 0.56935, and the sample path 0.56948.
    Lst<T> lstAoI = [lambda, H_lst, Tstar, one](const T& s) {
        if (num_traits<T>::to_double(s < T(0) ? T(-s) : s) < 1e-12)
            return one;  // A*(0) = 1 for any proper LST
        const T Hs = H_lst(s);
        const T Ts = Tstar(s);
        const T Tsl = Tstar(T(s + lambda));
        const T peak = Hs * (Ts - (s / (s + lambda)) * Tsl);
        return T((lambda / s) * (Ts - peak));
    };
    return {meanAoI, peakAoI, lstAoI, true};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_FCFS_MGI1_H
