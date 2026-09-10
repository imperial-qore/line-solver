/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_HH1_LINDLEY_H
#define LINE_API_QSYS_QSYS_HH1_LINDLEY_H

/**
 * Conditional waiting-time moments of the Hl/Hn/1 Lindley recursion.
 *
 * Templated port of matlab/src/api/qsys/qsys_hh1_lindley.m. No JAR
 * counterpart. Hyperexponential primitives are mixtures of exponentials, so
 * conditioning on the arrival phase i and the service phase j reduces one
 * Lindley step to the M/M/1 step at rates lambda(i) and mu(j), and
 *
 *   E[W_{n+1}^m | W_n] = sum_i sum_j pa(i) ps(j) E_ij[W_{n+1}^m | W_n].
 *
 * Phases are drawn independently for each customer, which is what makes the
 * mixture exact rather than an approximation; a Markov-modulated arrival
 * stream would not decompose this way.
 *
 * The variance is NOT the mixture of the per-phase variances, because the
 * phase is itself random. It is recovered from the first two MIXED raw
 * moments, which adds the between-phase spread of the means, and mixing the
 * variances instead would understate it.
 *
 * Unlike qsys_mm1_lindley the mean is read from moments(:,1) rather than from
 * an explicit closed form, since the mixture has none; that is also what
 * MATLAB does.
 *
 * Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
 * Recursion", Proc. Winter Simulation Conference, 2021, theorem 3.
 * Registered in .citations() as 'condlindley'.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/qsys/qsys_lindley_moment.h"
#include "line/api/qsys/qsys_mm1_lindley.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/**
 * @param lambda arrival phase rates, positive
 * @param pa     arrival phase probabilities, nonnegative and summing to 1
 * @param mu     service phase rates, positive
 * @param ps     service phase probabilities, nonnegative and summing to 1
 * @param Wn     waiting times of customer n, finite and nonnegative
 * @param mmax   highest moment order; raised to 2 when smaller, as in MATLAB
 */
template <class T>
LindleyResult<T> qsys_hh1_lindley(const std::vector<T>& lambda, const std::vector<T>& pa,
                                  const std::vector<T>& mu, const std::vector<T>& ps,
                                  const std::vector<T>& Wn, unsigned mmax = 2) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_hh1_lindley requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T tol = num_traits<T>::from_double(1e-10);
    if (lambda.empty() || lambda.size() != pa.size())
        throw InputError("qsys_hh1_lindley: lambda and pa must be nonempty and of equal length");
    if (mu.empty() || mu.size() != ps.size())
        throw InputError("qsys_hh1_lindley: mu and ps must be nonempty and of equal length");
    if (mmax < 1) throw InputError("qsys_hh1_lindley: mmax must be a positive integer");

    T sa = zero, ss = zero;
    for (std::size_t i = 0; i < lambda.size(); ++i) {
        if (lambda[i] <= zero)
            throw InputError("qsys_hh1_lindley: the arrival rates lambda must be positive real");
        if (pa[i] < zero) throw InputError("qsys_hh1_lindley: pa must be nonnegative");
        sa += pa[i];
    }
    for (std::size_t j = 0; j < mu.size(); ++j) {
        if (mu[j] <= zero)
            throw InputError("qsys_hh1_lindley: the service rates mu must be positive real");
        if (ps[j] < zero) throw InputError("qsys_hh1_lindley: ps must be nonnegative");
        ss += ps[j];
    }
    if (num_abs(T(sa - one)) > tol) throw InputError("qsys_hh1_lindley: pa must sum to 1");
    if (num_abs(T(ss - one)) > tol) throw InputError("qsys_hh1_lindley: ps must sum to 1");
    for (std::size_t i = 0; i < Wn.size(); ++i)
        if (Wn[i] < zero) throw InputError("qsys_hh1_lindley: Wn must be nonnegative");

    if (mmax < 2) mmax = 2;
    const std::size_t nw = Wn.size();
    LindleyResult<T> r;
    r.mmax = mmax;
    r.analyzer = "qsys_hh1_lindley";
    r.moments = Matrix<T>(nw, mmax, zero);
    for (std::size_t i = 0; i < lambda.size(); ++i)
        for (std::size_t j = 0; j < mu.size(); ++j) {
            const T weight = pa[i] * ps[j];
            if (weight == zero) continue;
            for (unsigned m = 1; m <= mmax; ++m) {
                const std::vector<T> col = qsys_lindley_moment(lambda[i], mu[j], Wn, m);
                for (std::size_t k = 0; k < nw; ++k) r.moments(k, m - 1) += weight * col[k];
            }
        }

    r.mean.reserve(nw);
    r.var.reserve(nw);
    for (std::size_t k = 0; k < nw; ++k) {
        r.mean.push_back(r.moments(k, 0));
        r.var.push_back(r.moments(k, 1) - r.moments(k, 0) * r.moments(k, 0));
    }
    return r;
}

/** Scalar-Wn overload. */
template <class T>
LindleyResult<T> qsys_hh1_lindley(const std::vector<T>& lambda, const std::vector<T>& pa,
                                  const std::vector<T>& mu, const std::vector<T>& ps,
                                  const T& Wn, unsigned mmax = 2) {
    return qsys_hh1_lindley(lambda, pa, mu, ps, std::vector<T>(1, Wn), mmax);
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_HH1_LINDLEY_H
