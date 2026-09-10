/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_TYPES_H
#define LINE_API_NPFQN_TYPES_H

/**
 * Shared arithmetic helpers for the templated npfqn port.
 *
 * The non-product-form traffic approximations in matlab/src/api/npfqn/ do not
 * share a single return type the way the qsys family does, so each ported
 * function declares its own result struct in its own header. What they do
 * share is a handful of ADL wrappers around the transcendental functions and
 * the finiteness test, which are collected here so that the same incantation
 * serves double, cpp_bin_float and mpfr. This header mirrors the role of
 * line/api/qsys/qsys_types.h and adds no algorithm of its own.
 */

#include <cmath>
#include <cstddef>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace npfqn {

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

/** Complementary error function, from Boost.Math so that it is type generic. */
template <class T>
inline T num_erfc(const T& v) {
    return boost::math::erfc(v);
}

/** MATLAB isfinite: false for +-Inf and NaN. Always true in an exact field. */
template <class T>
inline bool num_isfinite(const T& v) {
    return boost::math::isfinite(v);
}

template <>
inline bool num_isfinite<Rational>(const Rational&) {
    return true;  // a rational has neither an infinity nor a NaN
}

/** MATLAB isnan. Always false in an exact field. */
template <class T>
inline bool num_isnan(const T& v) {
    return boost::math::isnan(v);
}

template <>
inline bool num_isnan<Rational>(const Rational&) {
    return false;
}

template <class T>
inline const T& num_max(const T& a, const T& b) {
    return a < b ? b : a;
}

template <class T>
inline const T& num_min(const T& a, const T& b) {
    return a < b ? a : b;
}

}  // namespace detail

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_TYPES_H
