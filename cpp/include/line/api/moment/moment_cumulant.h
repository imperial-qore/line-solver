/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_CUMULANT_H
#define LINE_API_MOMENT_MOMENT_CUMULANT_H

/**
 * Cumulants from raw moments and the inverse, plus the factorial-cumulant pair.
 *
 * Templated port of matlab/src/api/moment/moment_cumulant_from_raw.m,
 * moment_raw_from_cumulant.m, moment_factcumulant_from_factorial.m and
 * moment_factorial_from_factcumulant.m. The recurrence is the one obtained by
 * differentiating log M(t) once, so it is triangular and every operation is
 * integer or rational.
 *
 * The factorial-cumulant pair is the SAME recurrence read on the factorial
 * sequence, which is why MATLAB delegates rather than duplicating it.
 *
 * Entry 0 of the cumulant vector is always 0 by convention, and entry 0 of the
 * raw vector is always 1, whatever the caller passes in.
 */

#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace moment {

/** kappa_i = m_i - sum_{k=1}^{i-1} C(i-1,k-1) kappa_k m_{i-k}. */
template <class T>
std::vector<T> moment_cumulant_from_raw(const std::vector<T>& m) {
    if (m.empty()) throw InputError("moment_cumulant_from_raw: m must be nonempty");
    const int n = static_cast<int>(m.size()) - 1;
    std::vector<T> kappa(m.size(), num_traits<T>::from_int(0));
    for (int i = 1; i <= n; ++i) {
        T acc = num_traits<T>::from_int(0);
        for (int k = 1; k <= i - 1; ++k) acc += num_nck<T>(i - 1, k - 1) * kappa[k] * m[i - k];
        kappa[i] = m[i] - acc;
    }
    return kappa;
}

/** m_i = sum_{k=1}^{i} C(i-1,k-1) kappa_k m_{i-k}, with m_0 = 1. */
template <class T>
std::vector<T> moment_raw_from_cumulant(const std::vector<T>& kappa) {
    if (kappa.empty()) throw InputError("moment_raw_from_cumulant: kappa must be nonempty");
    const int n = static_cast<int>(kappa.size()) - 1;
    std::vector<T> m(kappa.size(), num_traits<T>::from_int(0));
    m[0] = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i) {
        T acc = num_traits<T>::from_int(0);
        for (int k = 1; k <= i; ++k) acc += num_nck<T>(i - 1, k - 1) * kappa[k] * m[i - k];
        m[i] = acc;
    }
    return m;
}

/** Factorial cumulants from factorial moments. */
template <class T>
std::vector<T> moment_factcumulant_from_factorial(const std::vector<T>& f) {
    return moment_cumulant_from_raw<T>(f);
}

/** Factorial moments from factorial cumulants. */
template <class T>
std::vector<T> moment_factorial_from_factcumulant(const std::vector<T>& kappa) {
    return moment_raw_from_cumulant<T>(kappa);
}

}  // namespace moment
}  // namespace line

#endif
