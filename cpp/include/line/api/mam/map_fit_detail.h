/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_FIT_DETAIL_H
#define LINE_API_MAM_MAP_FIT_DETAIL_H

/**
 * Scalar helpers shared by the MAP/PH moment-matching headers.
 *
 * These are not ports of any MATLAB function: they are the arithmetic
 * primitives the closed forms need (square root, exponential, logarithm,
 * integer power) resolved through argument-dependent lookup so that the same
 * expression compiles for double and for Boost.Multiprecision types, plus a
 * minimal complex arithmetic type.
 *
 * The complex type exists because aph_fit's second fitting case (the
 * Bobbio-Horvath-Telek chain K9..K22) takes square and cube roots of
 * quantities that are negative for many *feasible* moment sets; the imaginary
 * parts cancel in the final f. MATLAB evaluates that chain in complex
 * arithmetic implicitly and the JAR does it explicitly with
 * org.apache.commons.math3.complex.Complex, so the port must too. Roots use
 * the principal branch, matching both references.
 */

#include <cmath>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {
namespace fitdetail {

template <class T>
inline T num_sqrt(const T& v) {
    using std::sqrt;
    return sqrt(v);
}

template <class T>
inline T num_exp(const T& v) {
    using std::exp;
    return exp(v);
}

template <class T>
inline T num_log(const T& v) {
    using std::log;
    return log(v);
}

template <class T>
inline T num_atan2(const T& y, const T& x) {
    using std::atan2;
    return atan2(y, x);
}

template <class T>
inline T num_cos(const T& v) {
    using std::cos;
    return cos(v);
}

template <class T>
inline T num_sin(const T& v) {
    using std::sin;
    return sin(v);
}

/** Shorthand for the integer powers that litter the symbolic expressions. */
template <class T>
inline T pw(const T& v, unsigned e) {
    return num_pow_int(v, e);
}

/**
 * Minimal complex number over T. Only the operations the aph_fit chain needs
 * are provided. Division is the naive formula: the magnitudes involved stay
 * far from the overflow limits of any of the supported backends.
 */
template <class T>
struct Cplx {
    T re;
    T im;

    Cplx() : re(num_traits<T>::from_int(0)), im(num_traits<T>::from_int(0)) {}
    explicit Cplx(const T& r) : re(r), im(num_traits<T>::from_int(0)) {}
    Cplx(const T& r, const T& i) : re(r), im(i) {}
};

template <class T>
inline Cplx<T> operator+(const Cplx<T>& a, const Cplx<T>& b) {
    return Cplx<T>(T(a.re + b.re), T(a.im + b.im));
}

template <class T>
inline Cplx<T> operator-(const Cplx<T>& a, const Cplx<T>& b) {
    return Cplx<T>(T(a.re - b.re), T(a.im - b.im));
}

template <class T>
inline Cplx<T> operator-(const Cplx<T>& a) {
    return Cplx<T>(T(-a.re), T(-a.im));
}

template <class T>
inline Cplx<T> operator*(const Cplx<T>& a, const Cplx<T>& b) {
    return Cplx<T>(T(a.re * b.re - a.im * b.im), T(a.re * b.im + a.im * b.re));
}

template <class T>
inline Cplx<T> operator*(const Cplx<T>& a, const T& s) {
    return Cplx<T>(T(a.re * s), T(a.im * s));
}

template <class T>
inline Cplx<T> operator/(const Cplx<T>& a, const Cplx<T>& b) {
    const T d = b.re * b.re + b.im * b.im;
    if (d == num_traits<T>::from_int(0)) throw NumericError("Cplx: division by zero");
    return Cplx<T>(T((a.re * b.re + a.im * b.im) / d), T((a.im * b.re - a.re * b.im) / d));
}

template <class T>
inline Cplx<T> operator/(const Cplx<T>& a, const T& s) {
    if (s == num_traits<T>::from_int(0)) throw NumericError("Cplx: division by zero");
    return Cplx<T>(T(a.re / s), T(a.im / s));
}

template <class T>
inline T cplx_abs(const Cplx<T>& a) {
    return num_sqrt(T(a.re * a.re + a.im * a.im));
}

/** Principal square root. */
template <class T>
inline Cplx<T> cplx_sqrt(const Cplx<T>& a) {
    const T zero = num_traits<T>::from_int(0);
    if (a.im == zero) {
        if (a.re >= zero) return Cplx<T>(num_sqrt(a.re), zero);
        return Cplx<T>(zero, num_sqrt(T(-a.re)));
    }
    const T r = cplx_abs(a);
    const T two = num_traits<T>::from_int(2);
    const T u = num_sqrt(T((r + a.re) / two));
    T v = num_sqrt(T((r - a.re) / two));
    if (a.im < zero) v = -v;
    return Cplx<T>(u, v);
}

/**
 * Principal z^x for a real exponent x: exp(x log z) with the principal
 * logarithm, which is what Apache Commons Math Complex.pow(double) computes.
 */
template <class T>
inline Cplx<T> cplx_pow_real(const Cplx<T>& a, const T& x) {
    const T zero = num_traits<T>::from_int(0);
    const T r = cplx_abs(a);
    if (r == zero) return Cplx<T>(zero, zero);
    const T theta = num_atan2(a.im, a.re);
    const T lr = num_log(r);
    const T mag = num_exp(T(x * lr));
    const T ang = x * theta;
    return Cplx<T>(T(mag * num_cos(ang)), T(mag * num_sin(ang)));
}

/** 1/z. */
template <class T>
inline Cplx<T> cplx_inv(const Cplx<T>& a) {
    return Cplx<T>(num_traits<T>::from_int(1)) / a;
}

}  // namespace fitdetail
}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_FIT_DETAIL_H
