/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_COXIAN_H
#define LINE_API_FJ_XMAX_COXIAN_H

/**
 * Expected maximum of K i.i.d. two-stage Coxian variables.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_coxian.m.
 *
 * With X = T1 + B T2, T1 ~ Exp(mu1), T2 ~ Exp(mu2) and B ~ Bernoulli(q), the
 * survival function is a two-term exponential mixture
 *
 *   S(t) = A exp(-mu1 t) + B exp(-mu2 t),
 *   A = (1-q) + q mu2/(mu2-mu1),   B = -q mu1/(mu2-mu1),
 *
 * so expanding 1 - (1-S)^K binomially and integrating term by term gives
 *
 *   E[Y_K] = sum_{j=1..K} (-1)^(j+1) C(K,j)
 *              sum_{i=0..j} C(j,i) A^(j-i) B^i / ((j-i) mu1 + i mu2).
 *
 * At coincident stage rates the mixture degenerates into
 * S(t) = (1 + q mu t) exp(-mu t) and the same expansion is carried out with
 * integral t^m exp(-j mu t) dt = m!/(j mu)^(m+1), which is selected
 * automatically.
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [Xmax, m1, c2] of fj_xmax_coxian. */
template <class T>
struct FJXmaxCoxianResult {
    T Xmax;
    T m1;
    T c2;
};

/**
 * @param K   number of branches, K >= 1
 * @param mu1 rate of the first stage
 * @param mu2 rate of the second stage
 * @param q   probability that the second stage is visited, in [0,1]
 * @return    the exact expected maximum with the branch mean and SCV
 */
template <class T>
FJXmaxCoxianResult<T> fj_xmax_coxian(unsigned K, const T& mu1, const T& mu2, const T& q) {
    detail::require_positive_K(K, "fj_xmax_coxian");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1),
            two = num_traits<T>::from_int(2);
    if (!(mu1 > zero) || !(mu2 > zero))
        throw InputError("fj_xmax_coxian: both stage rates must be positive");
    if (q < zero || q > one)
        throw InputError("fj_xmax_coxian: the branching probability must lie in [0,1]");
    if (K > 60) throw InputError("fj_xmax_coxian: the binomial expansion loses precision past K=60");

    FJXmaxCoxianResult<T> out;
    out.m1 = one / mu1 + q / mu2;
    const T var1 = one / (mu1 * mu1) + q * (two - q) / (mu2 * mu2);
    out.c2 = var1 / (out.m1 * out.m1);

    T acc = zero;
    const double gap = num_traits<T>::to_double(mu2 - mu1);
    const double scale = num_traits<T>::to_double(mu1 > mu2 ? mu1 : mu2);
    if ((gap > 0 ? gap : -gap) > 1e-12 * scale) {
        const T A = (one - q) + q * mu2 / (mu2 - mu1);
        const T B = -q * mu1 / (mu2 - mu1);
        for (unsigned j = 1; j <= K; ++j) {
            T inner = zero;
            for (unsigned i = 0; i <= j; ++i) {
                const T rate = num_traits<T>::from_int(static_cast<long>(j - i)) * mu1 +
                               num_traits<T>::from_int(static_cast<long>(i)) * mu2;
                T Apw = one, Bpw = one;
                for (unsigned e = 0; e < j - i; ++e) Apw *= A;
                for (unsigned e = 0; e < i; ++e) Bpw *= B;
                inner += detail::fj_binom<T>(j, i) * Apw * Bpw / rate;
            }
            const T term = detail::fj_binom<T>(K, j) * inner;
            if (j % 2 == 1) acc += term; else acc -= term;
        }
    } else {
        // Coincident stage rates: S(t) = (1 + q mu t) exp(-mu t)
        const T mu = mu1;
        for (unsigned j = 1; j <= K; ++j) {
            T inner = zero;
            const T jmu = num_traits<T>::from_int(static_cast<long>(j)) * mu;
            for (unsigned i = 0; i <= j; ++i) {
                T num = one, den = one, fact = one;
                for (unsigned e = 0; e < i; ++e) num *= q * mu;
                for (unsigned e = 2; e <= i; ++e) fact *= num_traits<T>::from_int(static_cast<long>(e));
                for (unsigned e = 0; e <= i; ++e) den *= jmu;
                inner += detail::fj_binom<T>(j, i) * num * fact / den;
            }
            const T term = detail::fj_binom<T>(K, j) * inner;
            if (j % 2 == 1) acc += term; else acc -= term;
        }
    }
    out.Xmax = acc;
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_COXIAN_H
