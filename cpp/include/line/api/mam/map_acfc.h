/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_ACFC_H
#define LINE_API_MAM_MAP_ACFC_H

/**
 * Autocorrelation of the counting process of a MAP at a given timescale.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_acfc.m. With Q = D0 + D1, pi
 * the stationary phase vector, u the slot length and tmp = (e pi - Q)^-1,
 *
 *   rho(k) = PRE exp(Q (k-1) u) POST / Var[N(u)],
 *   PRE  = pi D1 (I - exp(Q u)),
 *   POST = (I - exp(Q u)) tmp^2 D1 e,
 *
 * i.e. the lag-k covariance of the numbers of arrivals in consecutive windows
 * of length u, normalized by their common variance (map_varcount). This is the
 * counting-process autocorrelation, not the inter-arrival autocorrelation
 * map_acf, and the two are different functions of the same MAP: a renewal
 * process has zero inter-arrival autocorrelation at every lag but a nonzero
 * count autocorrelation unless it is Poisson.
 *
 * ARITHMETIC: transcendental, three matrix exponentials per lag.
 *
 * The JAR has no counterpart of this function.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_count_var.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_varcount.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * @param m    the MAP (D0, D1)
 * @param kset lags, each >= 1
 * @param u    length of the counting window (the timescale)
 * @return rho(k) for each lag, in the order of kset
 */
template <class T>
std::vector<T> map_acfc(const Map<T>& m, const std::vector<unsigned>& kset, const T& u) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_acfc requires transcendental arithmetic");
    const std::size_t n = m.order();
    if (!(u > num_traits<T>::from_int(0))) throw InputError("map_acfc: the timescale must be positive");
    const Matrix<T> Q = map_infgen(m);
    const std::vector<T> piq = map_prob(m);
    const Matrix<T> tmp = detail::map_count_deviation(Q, piq);
    const Matrix<T> tmp2 = matmul(tmp, tmp);
    const std::vector<T> e = ones<T>(n);

    const Matrix<T> Eu = expm(Q, u);
    const std::vector<T> piD1 = vecmul(piq, m.D1);
    const std::vector<T> piD1Eu = vecmul(piD1, Eu);
    std::vector<T> pre(n);
    for (std::size_t i = 0; i < n; ++i) pre[i] = piD1[i] - piD1Eu[i];  // pi D1 (I - e^{Qu})

    const std::vector<T> tail = mulvec(tmp2, mulvec(m.D1, e));  // tmp^2 D1 e
    const std::vector<T> Eutail = mulvec(Eu, tail);
    std::vector<T> post(n);
    for (std::size_t i = 0; i < n; ++i) post[i] = tail[i] - Eutail[i];  // (I - e^{Qu}) tmp^2 D1 e

    const std::vector<T> uv(1, u);
    const T vart = map_varcount(m, uv)[0];
    if (vart == num_traits<T>::from_int(0))
        throw NumericError("map_acfc: the counts have zero variance at this timescale");

    std::vector<T> out;
    out.reserve(kset.size());
    for (std::size_t j = 0; j < kset.size(); ++j) {
        if (kset[j] < 1u) throw InputError("map_acfc: lags must be at least 1");
        const T lagtime = num_traits<T>::from_int(static_cast<long>(kset[j] - 1u)) * u;
        std::vector<T> v = pre;
        if (!(lagtime == num_traits<T>::from_int(0))) v = vecmul(pre, expm(Q, lagtime));
        T s = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i) s += v[i] * post[i];
        out.push_back(s / vart);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_ACFC_H
