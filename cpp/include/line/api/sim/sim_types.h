/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_TYPES_H
#define LINE_API_SIM_SIM_TYPES_H

/**
 * Shared arithmetic helpers for the templated simulation output-analysis port.
 *
 * This family is the port of matlab/src/api/sim/ (jline.api.sim,
 * line_solver.api.sim). The domain is named for the statistics of a simulation
 * output process rather than for any queueing model: nothing in it takes a
 * NetworkStruct, the input is a sequence of observations such as successive
 * waiting times exported from a simulation run. It was called `oa`, for output
 * analysis, in all four codebases until 2026-07-31, so older commits, log
 * entries and manuals name it that way. This header holds the same role
 * npfqn_types.h plays for the traffic approximations; no algorithm lives here.
 *
 * ARITHMETIC. Every function in this family carries
 * `static_assert(num_traits<T>::has_transcendental)`, so exact rational
 * arithmetic is a build error rather than a silent fallback. That is not a
 * porting shortcut: the standardized time series areas carry the factor
 * weight/(m sqrt(m)) with the normalizing weight sqrt(12), and the interval
 * half width is a t quantile times a square root, so the delivered numbers are
 * irrational in the data. There is nothing for an exact field to preserve.
 *
 * The distributional quantiles (normal and Student t) are computed in double
 * and only then converted to T, exactly as lossn_mci's normal_quantile is: the
 * significance level alpha and the quantile order p reach these functions as
 * doubles, so what they determine carries double information and no more,
 * whatever the working type of the estimator is. The estimator's own
 * statistical error dwarfs the arithmetic one by many orders of magnitude.
 */

#include <cmath>
#include <cstddef>
#include <limits>

#include <boost/math/special_functions/fpclassify.hpp>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace sim {

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

/** MATLAB isfinite: false for +-Inf and NaN. */
template <class T>
inline bool num_isfinite(const T& v) {
    return boost::math::isfinite(v);
}

/** MATLAB isnan. */
template <class T>
inline bool num_isnan(const T& v) {
    return boost::math::isnan(v);
}

/**
 * MATLAB's NaN as a value of T. This family returns it where the reference
 * does, i.e. for a variance-parameter estimator with no degrees of freedom and
 * for a refused interval under force = false; a caller must test it rather than
 * read the field as a number.
 */
template <class T>
inline T num_nan() {
    return num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
}

}  // namespace detail

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_SIM_TYPES_H
