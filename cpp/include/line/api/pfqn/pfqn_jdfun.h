/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_JDFUN_H
#define LINE_API_PFQN_JDFUN_H

/**
 * AMVA joint-dependence function for non-product-form scaling.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_jdfun.m.
 *
 * Returns, for every station i, the reciprocal of the joint-dependent scaling
 * eta_i(n_{i1}, ..., n_{iR}) evaluated at the per-station population row
 * nvec(i,:), for the requested class.
 *
 * jdscaling[i] is a callable of that row. It may return a single value, i.e. a
 * scaling eta_i(n) shared by every class (broadcast, as in the flagship
 * min(n_i1, c)), or R values, i.e. the per-class scalings eta_{i,r}(n), of
 * which element classIdx is taken (Sauer chain-dependent rate mu_{r,i}(n)).
 *
 * Unlike pfqn_cdfun -- whose beta_{i,r} is a product-form scaling reading only
 * the own-class marginal n_{i,r} -- eta may read the joint vector arbitrarily
 * and is therefore NON-product-form: the AMVA result is an approximation with
 * no exactness or uniqueness guarantee. The numerical evaluation is identical
 * to pfqn_cdfun; the distinction is semantic and is carried by the separate
 * sn.jdscaling field.
 *
 * Arithmetic: EXACT-CAPABLE. One reciprocal per station and nothing else, so
 * no transcendental gate; exactness is a property of the supplied callables.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** A per-station joint-dependence callable: the population row -> 1 or R rates. */
template <class T>
using JdScaling = std::function<std::vector<T>(const std::vector<T>&)>;

/**
 * @param nvec      (M x R) per-station, per-class populations
 * @param jdscaling (M) callables; entries may be empty for "no scaling"
 * @param classIdx  0-based class whose scaling is selected from a vector result
 * @return (M) reciprocals of the scalings
 */
template <class T>
std::vector<T> pfqn_jdfun(const Matrix<T>& nvec, const std::vector<JdScaling<T>>& jdscaling,
                          std::size_t classIdx) {
    const std::size_t M = nvec.rows();
    const std::size_t R = nvec.cols();
    const T one = num_traits<T>::from_int(1);
    std::vector<T> r(M, one);
    if (jdscaling.empty()) return r;
    if (jdscaling.size() != M)
        throw InputError("pfqn_jdfun: scaling vector has the wrong station count");
    for (std::size_t i = 0; i < M; ++i) {
        if (!jdscaling[i]) continue;
        std::vector<T> row(R);
        for (std::size_t s = 0; s < R; ++s) row[s] = nvec(i, s);
        const std::vector<T> v = jdscaling[i](row);
        if (v.empty()) throw InputError("pfqn_jdfun: a scaling callable returned nothing");
        const T& chosen = v.size() > 1 ? v.at(classIdx) : v[0];
        if (chosen == num_traits<T>::from_int(0))
            throw NumericError("pfqn_jdfun: a joint-dependent scaling is zero");
        r[i] = one / chosen;
    }
    return r;
}

/** MATLAB default: classIdx = 1, i.e. the first class. */
template <class T>
std::vector<T> pfqn_jdfun(const Matrix<T>& nvec, const std::vector<JdScaling<T>>& jdscaling) {
    return pfqn_jdfun(nvec, jdscaling, static_cast<std::size_t>(0));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_JDFUN_H
