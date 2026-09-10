/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_DISPERSION_H
#define LINE_API_FJ_DISPERSION_H

/**
 * Mean subtask dispersion of a split-merge system with Erlang branches.
 *
 * Templated port of matlab/src/api/fj/fj_dispersion.m.
 *
 *   E[X_(N)] = integral_0^inf [ 1 - prod_i F_i(x - d_i) ] dx
 *   E[X_(1)] = integral_0^inf prod_i [ 1 - F_i(x - d_i) ] dx
 *   E[D_d]   = E[X_(N)] - E[X_(1)]
 *
 * The integrand actually evaluated is 1 - prod F_i - prod (1-F_i), which is
 * non-negative and vanishes at both ends; the difference of the two products
 * printed in the survey is not the dispersion and can go negative.
 *
 * Branch i is an Erlang with shape(i) stages of rate rate(i), the split-merge
 * equivalent used in the delay-scheduling construction: a subtask with q others
 * ahead of it in its parallel queue behaves as an Erlang(q+1, mu).
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/api/fj/fj_xmax_erlang.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [Edisp, Emax, Emin] of fj_dispersion. */
template <class T>
struct FJDispersionResult {
    T Edisp;
    T Emax;
    T Emin;
};

/**
 * @param shape   Erlang stage counts, one per branch
 * @param rate    Erlang stage rates, one per branch
 * @param d       deterministic delays, one per branch, all non-negative
 * @param tol     completion tolerance used to pick the quadrature horizon
 * @param npanels Simpson panel count, forced even
 * @return        the mean dispersion with the two order-statistic means
 */
template <class T>
FJDispersionResult<T> fj_dispersion(const std::vector<unsigned>& shape,
                                    const std::vector<T>& rate, const std::vector<T>& d,
                                    const T& tol = num_traits<T>::from_double(1e-10),
                                    unsigned npanels = 4000) {
    const std::size_t N = shape.size();
    if (rate.size() != N || d.size() != N)
        throw InputError("fj_dispersion: shape, rate and d must have the same length");
    if (N < 1) throw InputError("fj_dispersion: at least one branch is required");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1),
            two = num_traits<T>::from_int(2);
    for (std::size_t i = 0; i < N; ++i) {
        if (shape[i] < 1) throw InputError("fj_dispersion: Erlang stage counts must be positive");
        if (!(rate[i] > zero)) throw InputError("fj_dispersion: Erlang stage rates must be positive");
        if (d[i] < zero) throw InputError("fj_dispersion: delays must be non-negative");
    }
    if (npanels % 2 != 0) ++npanels;

    T dmax = d[0], mmax = num_traits<T>::from_int(static_cast<long>(shape[0])) / rate[0];
    for (std::size_t i = 0; i < N; ++i) {
        if (d[i] > dmax) dmax = d[i];
        const T mi = num_traits<T>::from_int(static_cast<long>(shape[i])) / rate[i];
        if (mi > mmax) mmax = mi;
    }
    T U = dmax + num_traits<T>::from_int(8) * mmax;
    for (unsigned it = 0; it < 60; ++it) {
        T prodF = one;
        for (std::size_t i = 0; i < N; ++i)
            prodF *= detail::erlang_cdf<T>(U - d[i], shape[i], rate[i]);
        if (one - prodF < tol) break;
        U = two * U;
    }

    const T h = U / num_traits<T>::from_int(static_cast<long>(npanels));
    const T four = num_traits<T>::from_int(4);
    T accMax = zero, accMin = zero;
    for (unsigned i = 0; i <= npanels; ++i) {
        const T x = h * num_traits<T>::from_int(static_cast<long>(i));
        T Fprod = one, Sprod = one;
        for (std::size_t j = 0; j < N; ++j) {
            const T Fj = detail::erlang_cdf<T>(x - d[j], shape[j], rate[j]);
            Fprod *= Fj;
            Sprod *= (one - Fj);
        }
        T w;
        if (i == 0 || i == npanels) w = one;
        else w = (i % 2 == 1) ? four : two;
        accMax += w * (one - Fprod);
        accMin += w * Sprod;
    }
    const T scale = h / num_traits<T>::from_int(3);

    FJDispersionResult<T> out;
    out.Emax = scale * accMax;
    out.Emin = scale * accMin;
    out.Edisp = out.Emax - out.Emin;
    return out;
}

/** The undelayed system, d = 0. */
template <class T>
FJDispersionResult<T> fj_dispersion(const std::vector<unsigned>& shape,
                                    const std::vector<T>& rate) {
    return fj_dispersion<T>(shape, rate, std::vector<T>(shape.size(), num_traits<T>::from_int(0)));
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_DISPERSION_H
