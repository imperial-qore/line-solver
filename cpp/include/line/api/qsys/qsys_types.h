/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_TYPES_H
#define LINE_API_QSYS_TYPES_H

/**
 * Shared return type and arithmetic helpers for the templated qsys port.
 *
 * The MATLAB closed-form queueing-system functions in matlab/src/api/qsys/
 * almost all return the pair [W,rhohat], mirrored in the JAR by Ret.qsys.
 * This header carries that pair plus the handful of ADL wrappers the family
 * needs; each ported function lives in its own header named after the MATLAB
 * file, as required by the port convention.
 */

#include <cmath>
#include <string>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/**
 * Return value of the qsys family, mirroring MATLAB's [W,rhohat] and the
 * JAR's Ret.qsys.
 *
 * W      mean response time (time in system, service included)
 * rhohat modified utilization chosen so that the M/M/1 relations still hold
 */
template <class T>
struct QsysResult {
    T W;
    T rhohat;
};

namespace detail {

/** exp(v), resolved by ADL so double, cpp_bin_float and mpfr all work. */
template <class T>
inline T num_exp(const T& v) {
    using std::exp;
    return exp(v);
}

/** sqrt(v), resolved by ADL. */
template <class T>
inline T num_sqrt(const T& v) {
    using std::sqrt;
    return sqrt(v);
}

/** base^exponent for a real-valued exponent, resolved by ADL. */
template <class T>
inline T num_pow(const T& base, const T& exponent) {
    using std::pow;
    return pow(base, exponent);
}

template <class T>
inline T num_min(const T& a, const T& b) {
    return a < b ? a : b;
}

/**
 * Guard against the exact pole at rho == 1.
 *
 * MATLAB divides by 1-rho and returns Inf there; in an exact field that
 * division has no value at all, so the port reports it as an input error.
 * Nothing else is checked: rho > 1 still returns MATLAB's negative W, since
 * the MATLAB functions do not reject it either.
 */
template <class E>
inline void require_no_pole(const E& one_minus_rho, const char* fn) {
    // E is deduced, not fixed to T, so a Boost expression template binds here
    // without being materialized first.
    if (one_minus_rho == 0)
        throw InputError(std::string(fn) + ": rho == 1, the mean waiting time is undefined");
}

/** rhohat = W*lambda/(1+W*lambda), the closing line of most qsys functions. */
template <class T>
inline T rhohat_from_W(const T& W, const T& lambda) {
    const T Wl = W * lambda;
    return Wl / (num_traits<T>::from_int(1) + Wl);
}

}  // namespace detail

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_TYPES_H
