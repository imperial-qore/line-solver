/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_LINDLEY_MOMENT_H
#define LINE_API_QSYS_QSYS_LINDLEY_MOMENT_H

/**
 * One conditional Lindley moment for exponential primitives.
 *
 * Templated port of matlab/src/api/qsys/qsys_lindley_moment.m. No JAR
 * counterpart. Returns E[max(Wn + S - A, 0)^m] with A ~ Exp(lambda),
 * S ~ Exp(mu) and m >= 1, evaluated at every entry of Wn.
 *
 * This is the algorithm shared by qsys_mm1_lindley, which calls it once per
 * moment order, and qsys_hh1_lindley, which mixes it over the arrival and
 * service phases. The density of S - A is the asymmetric Laplace density
 * lambda mu/(lambda+mu) times e^{-mu x} on x > 0 and e^{lambda x} on x < 0,
 * which splits the expectation into
 *
 *   S = sum_{k=0}^{m} C(m,k) w^k (m-k)! / mu^(m-k+1)
 *   T = (-1)^m m! ( sum_{k=0}^{m} (-lambda w)^k/k! - e^{-lambda w} ) / lambda^(m+1)
 *
 * with the value lambda mu/(lambda+mu) (S + T). The T term is the upper
 * incomplete gamma Gamma(m+1, -lambda w), which for integer m+1 has the finite
 * form m! e^{-x} sum_k x^k/k! valid at the negative argument needed here.
 * Substituting it cancels the growing exponential, so no incomplete gamma
 * routine is needed and nothing overflows at large w.
 *
 * ARITHMETIC: exp() is the only transcendental step, and it appears once per
 * evaluation point rather than inside the sum, so the exact instantiation is
 * refused rather than silently rounded. Everything else is a finite sum of
 * rational terms and is exact whenever the arithmetic is.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/**
 * @param lambda arrival rate, positive
 * @param mu     service rate, positive
 * @param Wn     current waiting times, nonnegative
 * @param m      moment order, at least 1
 * @return       E[max(Wn + S - A, 0)^m] evaluated at every entry of Wn
 */
template <class T>
std::vector<T> qsys_lindley_moment(const T& lambda, const T& mu, const std::vector<T>& Wn,
                                   unsigned m) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_lindley_moment requires transcendental arithmetic: the incomplete "
                  "gamma term carries an exp(-lambda w) that no finite field evaluates");
    const T zero = num_traits<T>::from_int(0);
    if (lambda <= zero) throw InputError("qsys_lindley_moment: lambda must be positive");
    if (mu <= zero) throw InputError("qsys_lindley_moment: mu must be positive");
    if (m < 1) throw InputError("qsys_lindley_moment: the moment order must be at least 1");

    // C(m,k) built by the multiplicative recurrence so it stays exact in T
    std::vector<T> binom(m + 1, num_traits<T>::from_int(1));
    for (unsigned k = 1; k <= m; ++k)
        binom[k] = binom[k - 1] * num_traits<T>::from_int(static_cast<long>(m - k + 1)) /
                   num_traits<T>::from_int(static_cast<long>(k));

    std::vector<T> out;
    out.reserve(Wn.size());
    const T pref = lambda * mu / (lambda + mu);
    const T sign = (m % 2 == 0) ? num_traits<T>::from_int(1) : num_traits<T>::from_int(-1);
    const T mfact = num_factorial<T>(m);
    for (std::size_t i = 0; i < Wn.size(); ++i) {
        const T& w = Wn[i];
        if (w < zero) throw InputError("qsys_lindley_moment: Wn must be nonnegative");
        T sTerm = zero;
        for (unsigned k = 0; k <= m; ++k)
            sTerm += binom[k] * num_pow_int(w, k) * num_factorial<T>(m - k) /
                     num_pow_int(mu, m - k + 1);
        const T x = -lambda * w;
        T inner = zero;
        for (unsigned k = 0; k <= m; ++k) inner += num_pow_int(x, k) / num_factorial<T>(k);
        const T tTerm = sign * mfact * (inner - detail::num_exp(x)) / num_pow_int(lambda, m + 1);
        out.push_back(pref * (sTerm + tTerm));
    }
    return out;
}

/** Scalar overload of the same expression. */
template <class T>
T qsys_lindley_moment(const T& lambda, const T& mu, const T& Wn, unsigned m) {
    return qsys_lindley_moment(lambda, mu, std::vector<T>(1, Wn), m)[0];
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_LINDLEY_MOMENT_H
