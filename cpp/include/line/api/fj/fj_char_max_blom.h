/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_CHAR_MAX_BLOM_H
#define LINE_API_FJ_CHAR_MAX_BLOM_H

/**
 * Blom-corrected plotting position for the characteristic maximum.
 *
 * Templated port of matlab/src/api/fj/fj_char_max_blom.m.
 *
 *   m_K = F^-1( (K - alpha) / (K - alpha - beta + 1) )
 *
 * which for alpha = beta = 0 falls back on the naive K/(K+1). The survey quotes
 * alpha = 0.4886 and beta = 0.3140, which are the defaults. For the standard
 * normal the position is bracketed without any inversion, for K >= 5, by
 *
 *   sqrt(2 ln K - ln ln K - 3) < m_K < sqrt(2 ln K - ln ln K).
 */

#include <cmath>
#include <functional>
#include <limits>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [mK, lo, hi] of fj_char_max_blom; lo and hi are NaN below K = 5. */
template <class T>
struct FJCharMaxBlomResult {
    T mK;
    T lo;
    T hi;
    bool bracket_available;
};

/**
 * @param Finv  quantile function; an empty target selects the standard normal
 * @param K     number of i.i.d. copies, K >= 1
 * @param alpha Blom numerator offset
 * @param beta  Blom denominator offset
 * @return      the corrected position and, for the normal, its bracket
 */
template <class T>
FJCharMaxBlomResult<T> fj_char_max_blom(unsigned K,
                                        const std::function<T(const T&)>& Finv =
                                            std::function<T(const T&)>(),
                                        const T& alpha = num_traits<T>::from_double(0.4886),
                                        const T& beta = num_traits<T>::from_double(0.3140)) {
    detail::require_positive_K(K, "fj_char_max_blom");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    const T den = Kt - alpha - beta + one;
    if (!(den > zero))
        throw NumericError("fj_char_max_blom: the Blom offsets leave a non-positive denominator");
    const T q = (Kt - alpha) / den;
    if (!(q > zero) || !(q < one))
        throw NumericError("fj_char_max_blom: the plotting position fell outside (0,1)");

    FJCharMaxBlomResult<T> out;
    if (Finv) {
        out.mK = Finv(q);
    } else {
        // Standard normal quantile through the inverse error function
        const double qd = num_traits<T>::to_double(q);
        out.mK = num_traits<T>::from_double(detail::normal_quantile(qd));
    }

    out.lo = num_traits<T>::from_int(0);
    out.hi = num_traits<T>::from_int(0);
    out.bracket_available = false;
    if (K >= 5) {
        const double z = 2.0 * std::log(static_cast<double>(K)) -
                         std::log(std::log(static_cast<double>(K)));
        if (z > 3.0) {
            out.lo = num_traits<T>::from_double(std::sqrt(z - 3.0));
            out.hi = num_traits<T>::from_double(std::sqrt(z));
            out.bracket_available = true;
        }
    }
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_CHAR_MAX_BLOM_H
