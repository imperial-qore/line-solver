/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_HZ_HET_H
#define LINE_API_FJ_XMAX_HZ_HET_H

/**
 * Harrison-Zertal approximation of the maximum of general variables.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_hz_het.m.
 *
 *   I(S) = (1/|S|) sum_{i in S} [ I(S \ i) + (m2_i/(2 m1_i)) L*_{S\i}(alpha_i) ]
 *
 * anchored at I({i}) = m1_i, with alpha_i = 1/m1_i. The transform of the
 * maximum over a sub-collection is recovered from the distribution functions,
 *
 *   L*_T(s) = s integral_0^inf exp(-s t) prod_{j in T} F_j(t) dt,
 *
 * by composite Simpson quadrature on a horizon widened until the product of the
 * distribution functions is within TOL of one. For identical branches the
 * recurrence collapses onto fj_xmax_hz, and for identical exponential branches
 * it is exact at H_K/lambda.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

namespace detail {

/** L*_T(s) for the sub-collection selected by mask, by composite Simpson. */
template <class T>
inline T hz_lst_max(const std::vector<std::function<T(const T&)> >& cdf, std::size_t mask,
                    std::size_t K, const T& s, const T& U, unsigned npanels) {
    if (mask == 0) return num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);
    const T h = (U - zero) / num_traits<T>::from_int(static_cast<long>(npanels));
    const T two = num_traits<T>::from_int(2), four = num_traits<T>::from_int(4);
    T acc = zero;
    for (unsigned i = 0; i <= npanels; ++i) {
        const T t = h * num_traits<T>::from_int(static_cast<long>(i));
        T g = num_exp<T>(-s * t);
        for (std::size_t j = 0; j < K; ++j)
            if (mask & (static_cast<std::size_t>(1) << j)) g *= cdf[j](t);
        T w;
        if (i == 0 || i == npanels) w = num_traits<T>::from_int(1);
        else w = (i % 2 == 1) ? four : two;
        acc += w * g;
    }
    return s * (h / num_traits<T>::from_int(3)) * acc;
}

}  // namespace detail

/**
 * @param m1      the K branch means, all positive
 * @param m2      the K branch second moments, m2[i] >= m1[i]^2
 * @param cdf     the K distribution functions, cdf[i](t) = P(X_i <= t)
 * @param tol     completion tolerance used to pick the quadrature horizon
 * @param npanels Simpson panel count, forced even
 * @return        the approximate expected maximum
 */
template <class T>
T fj_xmax_hz_het(const std::vector<T>& m1, const std::vector<T>& m2,
                 const std::vector<std::function<T(const T&)> >& cdf,
                 const T& tol = num_traits<T>::from_double(1e-10), unsigned npanels = 2000) {
    const std::size_t K = m1.size();
    if (m2.size() != K || cdf.size() != K)
        throw InputError("fj_xmax_hz_het: m1, m2 and cdf must have the same length");
    if (K < 1) throw InputError("fj_xmax_hz_het: at least one branch is required");
    if (K > 14)
        throw InputError(
            "fj_xmax_hz_het: the recurrence enumerates 2^K sub-collections with a quadrature each");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1),
            two = num_traits<T>::from_int(2);
    for (std::size_t i = 0; i < K; ++i) {
        if (!(m1[i] > zero)) throw InputError("fj_xmax_hz_het: all branch means must be positive");
        if (m2[i] < m1[i] * m1[i])
            throw InputError("fj_xmax_hz_het: a second moment is below the square of its mean");
    }
    if (npanels % 2 != 0) ++npanels;

    std::vector<T> alpha(K), resid(K);
    T mmax = m1[0];
    for (std::size_t i = 0; i < K; ++i) {
        alpha[i] = one / m1[i];
        resid[i] = m2[i] / (two * m1[i]);
        if (m1[i] > mmax) mmax = m1[i];
    }

    // Horizon: widen until every branch is essentially complete
    T U = num_traits<T>::from_int(8) * mmax;
    for (unsigned it = 0; it < 60; ++it) {
        T prodF = one;
        for (std::size_t j = 0; j < K; ++j) prodF *= cdf[j](U);
        if (one - prodF < tol) break;
        U = two * U;
    }

    const std::size_t nmask = static_cast<std::size_t>(1) << K;
    std::vector<T> Ival(nmask, zero);
    for (std::size_t mask = 1; mask < nmask; ++mask) {
        std::vector<std::size_t> members;
        for (std::size_t i = 0; i < K; ++i)
            if (mask & (static_cast<std::size_t>(1) << i)) members.push_back(i);
        if (members.size() == 1) {
            Ival[mask] = m1[members[0]];
            continue;
        }
        T acc = zero;
        for (std::size_t idx = 0; idx < members.size(); ++idx) {
            const std::size_t i = members[idx];
            const std::size_t rest = mask ^ (static_cast<std::size_t>(1) << i);
            acc += Ival[rest] +
                   resid[i] * detail::hz_lst_max<T>(cdf, rest, K, alpha[i], U, npanels);
        }
        Ival[mask] = acc / num_traits<T>::from_int(static_cast<long>(members.size()));
    }
    return Ival[nmask - 1];
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_HZ_HET_H
