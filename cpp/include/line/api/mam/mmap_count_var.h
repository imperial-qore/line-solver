/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_COUNT_VAR_H
#define LINE_API_MAM_MMAP_COUNT_VAR_H

/**
 * Per-class variance of the counting process of a marked MAP.
 *
 * Templated port of matlab/lib/m3a/m3a/mmap/mmap_count_var.m. It is the He and
 * Neuts variance of map_count_var applied class by class: for class k, with
 * D = D0 + D1 the generator of the phase process (the AGGREGATE one, not the
 * class-k one), theta its stationary vector and tmp = (e theta - D)^-1,
 *
 *   Var[N_k(t)] = (lambda_k - 2 lambda_k^2 + 2 theta Dk tmp Dk e) t
 *                 - 2 theta Dk tmp (I - exp(D t)) tmp Dk e,
 *
 * with lambda_k = theta Dk e. The phase process is shared across classes, so
 * the exponential is computed once and reused for every class.
 *
 * ARITHMETIC: transcendental, as map_count_var.
 *
 * This function was previously out of reach of the port for want of expm only;
 * the other MMAP counting descriptors (mmap_count_mean, mmap_count_lambda) are
 * rational and are already in mmap_lambda.h.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_count_var.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * @param mm the marked MAP
 * @param t  window length
 * @return Var[N_k(t)] for each class k
 */
template <class T>
std::vector<T> mmap_count_var(const Mmap<T>& mm, const T& t) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmap_count_var requires transcendental arithmetic");
    if (t < num_traits<T>::from_int(0)) throw InputError("mmap_count_var: negative window length");
    const std::size_t n = mm.order();
    const std::size_t K = mm.classes();
    if (K == 0) throw InputError("mmap_count_var: the MMAP has no classes");

    const Map<T> base = mm.map();
    const Matrix<T> D = map_infgen(base);
    const std::vector<T> theta = map_prob(base);
    const Matrix<T> tmp = detail::map_count_deviation(D, theta);
    const std::vector<T> e = ones<T>(n);
    const Matrix<T> E = expm(D, t);
    const T two = num_traits<T>::from_int(2);

    std::vector<T> out;
    out.reserve(K);
    for (std::size_t k = 0; k < K; ++k) {
        const Matrix<T>& Dk = mm.Dc[k];
        const std::vector<T> thetaDk = vecmul(theta, Dk);
        T lam = num_traits<T>::from_int(0);
        for (const T& v : thetaDk) lam += v;
        const std::vector<T> c = vecmul(thetaDk, tmp);         // theta Dk tmp
        const std::vector<T> d = mulvec(tmp, mulvec(Dk, e));   // tmp Dk e
        const std::vector<T> cDk = vecmul(c, Dk);
        T cDke = num_traits<T>::from_int(0);
        for (const T& v : cDk) cDke += v;

        const std::vector<T> cE = vecmul(c, E);
        T corr = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i) corr += (c[i] - cE[i]) * d[i];
        out.push_back((lam - two * lam * lam + two * cDke) * t - two * corr);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAP_COUNT_VAR_H
