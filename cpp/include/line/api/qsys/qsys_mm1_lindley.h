/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MM1_LINDLEY_H
#define LINE_API_QSYS_QSYS_MM1_LINDLEY_H

/**
 * Conditional waiting-time moments of the M/M/1 Lindley recursion.
 *
 * Templated port of matlab/src/api/qsys/qsys_mm1_lindley.m. No JAR
 * counterpart. One step of W_{n+1} = max(W_n + S_n - A_n, 0) with
 * A_n ~ Exp(lambda) and S_n ~ Exp(mu), returning the conditional raw moments
 * of orders 1..mmax given W_n.
 *
 * Unlike every other qsys_* function the quantities here are conditional on
 * the current state rather than stationary, so they are defined and finite at
 * any load, lambda >= mu included. There is no stability check for that
 * reason, and adding one would reject the very regime the recursion is used
 * to study.
 *
 * The mean is returned from the equivalent explicit form
 *   E[W_{n+1} | W_n] = W_n + (lambda-mu)/(lambda mu)
 *                      + mu e^{-lambda W_n} / (lambda (lambda+mu))
 * rather than from moments(:,1), which is what MATLAB does: the two agree in
 * exact arithmetic, and the explicit form is the better conditioned of the
 * two at small lambda where the moment expression differences two large
 * quantities.
 *
 * Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
 * Recursion", Proc. Winter Simulation Conference, 2021, theorem 1 and
 * corollary 2. Registered in .citations() as 'condlindley'.
 *
 * ARITHMETIC: see qsys_lindley_moment.h; exp() makes the exact instantiation
 * unavailable.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/qsys/qsys_lindley_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Mirrors the struct MATLAB returns from qsys_mm1_lindley. */
template <class T>
struct LindleyResult {
    std::vector<T> mean;     ///< conditional mean, one per Wn entry
    std::vector<T> var;      ///< conditional variance, one per Wn entry
    Matrix<T> moments;       ///< (numel(Wn) x mmax) conditional raw moments
    unsigned mmax = 2;       ///< highest moment order computed
    std::string analyzer;    ///< identifier string, as in MATLAB
};

/**
 * @param lambda arrival rate, positive
 * @param mu     service rate, positive
 * @param Wn     waiting times of customer n, finite and nonnegative
 * @param mmax   highest moment order; raised to 2 when smaller, as in MATLAB
 */
template <class T>
LindleyResult<T> qsys_mm1_lindley(const T& lambda, const T& mu, const std::vector<T>& Wn,
                                  unsigned mmax = 2) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mm1_lindley requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    if (lambda <= zero) throw InputError("qsys_mm1_lindley: lambda must be a positive real");
    if (mu <= zero) throw InputError("qsys_mm1_lindley: mu must be a positive real");
    if (mmax < 1) throw InputError("qsys_mm1_lindley: mmax must be a positive integer");
    for (std::size_t i = 0; i < Wn.size(); ++i)
        if (Wn[i] < zero) throw InputError("qsys_mm1_lindley: Wn must be nonnegative");

    if (mmax < 2) mmax = 2;
    const std::size_t nw = Wn.size();
    LindleyResult<T> r;
    r.mmax = mmax;
    r.analyzer = "qsys_mm1_lindley";
    r.moments = Matrix<T>(nw, mmax, zero);
    for (unsigned m = 1; m <= mmax; ++m) {
        const std::vector<T> col = qsys_lindley_moment(lambda, mu, Wn, m);
        for (std::size_t i = 0; i < nw; ++i) r.moments(i, m - 1) = col[i];
    }

    r.mean.reserve(nw);
    r.var.reserve(nw);
    for (std::size_t i = 0; i < nw; ++i) {
        const T m1 = Wn[i] + (lambda - mu) / (lambda * mu) +
                     mu * detail::num_exp(T(-lambda * Wn[i])) / (lambda * (lambda + mu));
        r.mean.push_back(m1);
        r.var.push_back(r.moments(i, 1) - r.moments(i, 0) * r.moments(i, 0));
    }
    return r;
}

/** Scalar overload. */
template <class T>
LindleyResult<T> qsys_mm1_lindley(const T& lambda, const T& mu, const T& Wn,
                                  unsigned mmax = 2) {
    return qsys_mm1_lindley(lambda, mu, std::vector<T>(1, Wn), mmax);
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MM1_LINDLEY_H
