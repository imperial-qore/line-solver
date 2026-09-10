/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_TYPES_H
#define LINE_API_AOI_TYPES_H

/**
 * Shared return types and arithmetic helpers for the templated Age of
 * Information port.
 *
 * The MATLAB AoI family in matlab/src/api/aoi/ returns one of two triples:
 * [meanAoI, varAoI, peakAoI] for the fully parameterized queues (M/M/1,
 * M/D/1, D/M/1) and [meanAoI, lstAoI, peakAoI] for the ones taking a
 * Laplace-Stieltjes transform as input (M/GI/1, GI/M/1). Both live here,
 * mirroring jline.api.aoi.AoiResult and jline.api.aoi.AoiLstResult; each
 * ported function lives in its own header named after the MATLAB file, as
 * required by the port convention.
 *
 * An LST is represented as a std::function<T(const T&)> rather than a
 * distribution object, exactly as MATLAB represents it as a function handle
 * and the JAR as a LstFunction. Nothing in the family inverts an LST: the
 * transforms are evaluated at real arguments only, and where MATLAB needs a
 * derivative it takes a central difference with a fixed step.
 */

#include <cmath>
#include <functional>
#include <string>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/** A Laplace-Stieltjes transform evaluated at real arguments. */
template <class T>
using Lst = std::function<T(const T&)>;

/** [meanAoI, varAoI, peakAoI], mirroring jline.api.aoi.AoiResult. */
template <class T>
struct AoiResult {
    T meanAoI;
    T varAoI;
    T peakAoI;
};

/**
 * [meanAoI, lstAoI, peakAoI], mirroring jline.api.aoi.AoiLstResult.
 *
 * has_lst is false for the LCFS-D and LCFS-S disciplines, where MATLAB
 * returns [] because the transform has no tractable closed form.
 */
template <class T>
struct AoiLstResult {
    T meanAoI;
    T peakAoI;
    Lst<T> lstAoI;
    bool has_lst;
};

namespace detail {

/** exp(v), resolved by ADL so double, cpp_bin_float and mpfr all work. */
template <class T>
inline T num_exp(const T& v) {
    using std::exp;
    return exp(v);
}

/** Guard on positivity of a rate or a time. */
template <class T>
inline void require_positive(const T& v, const char* fn, const char* what) {
    if (v <= num_traits<T>::from_int(0))
        throw InputError(std::string(fn) + ": " + what + " must be positive");
}

/** Guard on stability; every function in the family requires rho < 1. */
template <class T>
inline void require_stable(const T& rho, const char* fn) {
    if (rho >= num_traits<T>::from_int(1))
        throw NumericError(std::string(fn) + ": unstable system, rho >= 1");
}

/**
 * Central difference of an LST at s with MATLAB's step, h = 1e-6 max(1,|s|).
 *
 * Every MATLAB file in the family that needs d/ds H*(s) uses exactly this
 * step; reproducing it, rather than improving on it, is what makes the port
 * agree with MATLAB to the twelfth digit. The approximation is the reason
 * these functions require transcendental arithmetic even when the transform
 * itself is rational: the answer is step-dependent, so an exact field would
 * return an exact value of the wrong quantity.
 */
template <class T>
inline T lst_derivative(const Lst<T>& f, const T& s) {
    const T one = num_traits<T>::from_int(1);
    const T mag = s < num_traits<T>::from_int(0) ? T(-s) : s;
    const T scale = mag > one ? mag : one;
    const T h = num_traits<T>::from_double(1e-6) * scale;
    return (f(T(s + h)) - f(T(s - h))) / (num_traits<T>::from_int(2) * h);
}

/**
 * Bisection on a bracketed sign change, the replacement for MATLAB's fzero.
 *
 * MATLAB brackets every sigma root in this family on [0.001, 0.999] and calls
 * fzero (Brent); the JAR replaces it with bisection (Aoi_fcfs_dm1). This port
 * follows the JAR. With max_iter = 200 the bracket collapses below the
 * representable resolution of double and of the 50-digit real type alike.
 */
template <class T, class F>
inline T bisect(const F& f, T lo, T hi, const char* fn, unsigned max_iter = 200) {
    const T zero = num_traits<T>::from_int(0);
    T flo = f(lo);
    const T fhi = f(hi);
    if ((flo > zero && fhi > zero) || (flo < zero && fhi < zero))
        throw NumericError(std::string(fn) + ": the root is not bracketed by [0.001, 0.999]");
    for (unsigned it = 0; it < max_iter; ++it) {
        const T mid = (lo + hi) / num_traits<T>::from_int(2);
        if (mid == lo || mid == hi) break;
        const T fm = f(mid);
        if (fm == zero) return mid;
        if ((fm > zero) == (flo > zero)) { lo = mid; flo = fm; }
        else { hi = mid; }
    }
    return (lo + hi) / num_traits<T>::from_int(2);
}

/**
 * The sigma of a GI/M/1 queue: the unique root in (0,1) of Y*(mu - mu s) = s,
 * the probability an arriving job finds the server busy.
 */
template <class T>
inline T gim1_sigma(const Lst<T>& Y_lst, const T& mu, const char* fn) {
    return bisect<T>([&](const T& s) { return T(Y_lst(T(mu - mu * s)) - s); },
                     num_traits<T>::from_double(0.001), num_traits<T>::from_double(0.999), fn);
}

}  // namespace detail

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_TYPES_H
