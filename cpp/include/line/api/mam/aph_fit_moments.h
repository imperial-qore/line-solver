/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH_FIT_MOMENTS_H
#define LINE_API_MAM_APH_FIT_MOMENTS_H

/**
 * Acyclic phase-type fitters from the first two moments.
 *
 * Port of BUTools' `APHFrom2Moments` and of MATLAB `APH.fitMeanAndSCV`. Both
 * lived in solver_mna.h until a second caller appeared: the closed
 * setup/delay-off branch of solver_mam_basic fits its delay-off with
 * `APH.fitMeanAndSCV`, exactly as the reference does, and solver_mna.h already
 * includes solver_mam_basic.h, so leaving them there would have been a cycle.
 *
 * These are NOT the canonical Coxian of `coxian_phase_subgen`
 * (qbd_setupdelayoff.h). The two agree in the first two moments but not in
 * shape: the Coxian is entered at phase 1, this APH is entered at phase 1 with
 * probability p and at the LAST phase otherwise. Anything reading more than the
 * first two moments -- an LST, for instance -- sees the difference, so a caller
 * must use whichever one the reference names for that call site.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Port of BUTools' `APHFrom2Moments`.
 *
 * Absorption is possible ONLY from the last phase: the first N-1 rows of the
 * generator have a zero row sum by construction, so D1 is zero everywhere
 * except its last row. That is the shape, not an artifact of the fit.
 */
template <class T>
Map<T> aph_from_2moments(const T& e1, const T& e2) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T cv2 = T(T(e2 / T(e1 * e1)) - one);
    const double cv2d = num_traits<T>::to_double(cv2);
    if (!(cv2d > 0.0))
        throw NumericError(
            "APHFrom2Moments: the moment pair implies a non-positive squared coefficient of "
            "variation (" + std::to_string(cv2d) +
            "), for which the reference's order ceil(1/cv2) is not a positive integer");
    const T lambda = T(one / e1);
    const long N = std::max(static_cast<long>(std::ceil(1.0 / cv2d)), 2L);
    const std::size_t n = static_cast<std::size_t>(N);
    const T Nt = num_traits<T>::from_int(N);
    const T p = T(one / T(cv2 + one + T(T(cv2 - one) / num_traits<T>::from_int(N - 1))));

    Matrix<T> A(n, n, zero);
    const T d = T(lambda * p * Nt);
    for (std::size_t i = 0; i < n; ++i) A(i, i) = T(-d);
    for (std::size_t i = 0; i + 1 < n; ++i) A(i, i + 1) = d;
    A(n - 1, n - 1) = T(-T(lambda * Nt));

    std::vector<T> alpha(n, zero);
    alpha[0] = p;
    alpha[n - 1] = T(one - p);

    Map<T> m;
    m.D0 = A;
    m.D1 = Matrix<T>(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T ex = zero;
        for (std::size_t j = 0; j < n; ++j) ex += A(i, j);
        for (std::size_t j = 0; j < n; ++j) m.D1(i, j) = T(-ex * alpha[j]);
    }
    return m;
}

/** Port of `APH.fitMeanAndSCV`, the entry point the analyzers fit arrivals with. */
template <class T>
Map<T> aph_fit_mean_scv(const T& mean, const T& scv) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (num_traits<T>::to_double(mean) <= lang::GlobalConstants::FineTol) {
        // A sub-tolerance mean is answered by an exponential, and a
        // non-positive one by the Zero-mean surrogate the reference names.
        const T m = (mean > zero) ? mean : num_traits<T>::from_double(lang::GlobalConstants::Zero);
        return map_exponential(T(one / m));
    }
    if (scv == one) return map_exponential(T(one / mean));
    return aph_from_2moments(mean, T(T(one + scv) * mean * mean));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH_FIT_MOMENTS_H
