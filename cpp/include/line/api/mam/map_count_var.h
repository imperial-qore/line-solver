/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_COUNT_VAR_H
#define LINE_API_MAM_MAP_COUNT_VAR_H

/**
 * Variance of the counting process of a MAP at resolution t.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_count_var.m, cross-checked
 * against jar/src/main/java/jline/api/mam/Map_count_var.java. With
 * D = D0 + D1, theta the stationary phase vector, e the vector of ones and
 * the deviation matrix tmp = (e theta - D)^-1,
 *
 *   Var[N(t)] = (lambda - 2 lambda^2 + 2 theta D1 tmp D1 e) t
 *               - 2 theta D1 tmp (I - exp(D t)) tmp D1 e,
 *
 * from He and Neuts, "Markov chains with marked transitions" (1998).
 *
 * ARITHMETIC: the linear part is rational in the entries, the transient
 * correction needs exp(D t), so the function requires transcendental
 * arithmetic. The two terms have opposite signs and nearly cancel for small t,
 * where Var -> lambda t; that cancellation is the reason the high-precision
 * backends are useful here.
 *
 * map_varcount.m is the same quantity written with (e theta - D)^-2 in the
 * middle instead of tmp on both sides. The two agree identically, because
 * D (e theta) = (D e) theta = 0 and (e theta) D = e (theta D) = 0, so e theta
 * and D commute and hence tmp commutes with exp(D t). Both spellings are
 * ported (see map_varcount.h) and the tests assert their agreement.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace detail {

/**
 * Deviation matrix (e theta - D)^-1 of the phase process, MATLAB's tmp. Shared
 * by every counting-process descriptor.
 */
template <class T>
Matrix<T> map_count_deviation(const Matrix<T>& D, const std::vector<T>& theta) {
    const std::size_t n = D.rows();
    Matrix<T> A(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = theta[j] - D(i, j);
    return inverse(A);
}

}  // namespace detail

/**
 * @param m the MAP (D0, D1)
 * @param t window lengths
 * @return Var[N(t)] for each window length, in the order of t
 */
template <class T>
std::vector<T> map_count_var(const Map<T>& m, const std::vector<T>& t) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_count_var requires transcendental arithmetic");
    const std::size_t n = m.order();
    const Matrix<T> D = map_infgen(m);
    const std::vector<T> theta = map_prob(m);
    const Matrix<T> tmp = detail::map_count_deviation(D, theta);
    const std::vector<T> e = ones<T>(n);

    const std::vector<T> thetaD1 = vecmul(theta, m.D1);
    T lam = num_traits<T>::from_int(0);
    for (const T& v : thetaD1) lam += v;                 // lambda = theta D1 e
    const std::vector<T> c = vecmul(thetaD1, tmp);       // theta D1 tmp
    const std::vector<T> d = mulvec(tmp, mulvec(m.D1, e));  // tmp D1 e

    const std::vector<T> cD1 = vecmul(c, m.D1);
    T cD1e = num_traits<T>::from_int(0);
    for (const T& v : cD1) cD1e += v;
    const T two = num_traits<T>::from_int(2);
    const T linear = lam - two * lam * lam + two * cD1e;

    std::vector<T> out;
    out.reserve(t.size());
    for (std::size_t k = 0; k < t.size(); ++k) {
        if (t[k] < num_traits<T>::from_int(0))
            throw InputError("map_count_var: negative window length");
        // c (I - exp(D t)) d
        const Matrix<T> E = expm(D, t[k]);
        const std::vector<T> cE = vecmul(c, E);
        T corr = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i) corr += (c[i] - cE[i]) * d[i];
        out.push_back(linear * t[k] - two * corr);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_COUNT_VAR_H
