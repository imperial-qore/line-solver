/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CDFUN_H
#define LINE_API_PFQN_CDFUN_H

/**
 * AMVA-QD class-dependence function.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_cdfun.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/ld/Pfqn_cdfun.java.
 *
 * Returns, for every station i, the reciprocal of the class-dependent scaling
 * beta_{i,r}(n_{i1}, ..., n_{iR}) evaluated at the per-class population row
 * nvec(i,:), for the requested class r.
 *
 * cdscaling[i] is a callable of that row. It may return a single value, i.e. a
 * chain-independent scaling beta_i(n) shared by every class (the historical
 * contract), or R values, i.e. the per-class scalings, of which element
 * classIdx is taken. The per-class form expresses Sauer's chain-dependent
 * service rates mu_{r,i}(n) (Sauer 1983, eq. (40)). An absent (empty) callable
 * leaves the station at 1.
 *
 * Arithmetic: EXACT-CAPABLE. The function itself performs one reciprocal per
 * station and nothing else, so it carries no transcendental gate; whether the
 * result is exact is entirely a property of the supplied callables.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** A per-station class-dependence callable: the population row -> 1 or R rates. */
template <class T>
using CdScaling = std::function<std::vector<T>(const std::vector<T>&)>;

/**
 * @param nvec      (M x R) per-station, per-class populations
 * @param cdscaling (M) callables; entries may be empty for "no scaling"
 * @param classIdx  0-based class whose scaling is selected from a vector result
 * @return (M) reciprocals of the scalings
 */
template <class T>
std::vector<T> pfqn_cdfun(const Matrix<T>& nvec, const std::vector<CdScaling<T>>& cdscaling,
                          std::size_t classIdx) {
    const std::size_t M = nvec.rows();
    const std::size_t R = nvec.cols();
    const T one = num_traits<T>::from_int(1);
    std::vector<T> r(M, one);
    if (cdscaling.empty()) return r;
    if (cdscaling.size() != M)
        throw InputError("pfqn_cdfun: scaling vector has the wrong station count");
    for (std::size_t i = 0; i < M; ++i) {
        if (!cdscaling[i]) continue;
        std::vector<T> row(R);
        for (std::size_t s = 0; s < R; ++s) row[s] = nvec(i, s);
        const std::vector<T> v = cdscaling[i](row);
        if (v.empty()) throw InputError("pfqn_cdfun: a scaling callable returned nothing");
        const T& chosen = v.size() > 1 ? v.at(classIdx) : v[0];
        if (chosen == num_traits<T>::from_int(0))
            throw NumericError("pfqn_cdfun: a class-dependent scaling is zero");
        r[i] = one / chosen;
    }
    return r;
}

/** MATLAB default: classIdx = 1, i.e. the first class. */
template <class T>
std::vector<T> pfqn_cdfun(const Matrix<T>& nvec, const std::vector<CdScaling<T>>& cdscaling) {
    return pfqn_cdfun(nvec, cdscaling, static_cast<std::size_t>(0));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CDFUN_H
