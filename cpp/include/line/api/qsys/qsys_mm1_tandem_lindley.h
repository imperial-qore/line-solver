/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MM1_TANDEM_LINDLEY_H
#define LINE_API_QSYS_QSYS_MM1_TANDEM_LINDLEY_H

/**
 * Conditional waiting time at the SECOND station of an M/M/1 -> /M/1 tandem.
 *
 * Templated port of matlab/src/api/qsys/qsys_mm1_tandem_lindley.m. No JAR
 * counterpart. Station 1 is M/M/1 at rates (lambda, mu1) and station 2 is a
 * single server at rate mu2 fed by its departures. Given that customer n
 * waited Wk at station 1 and Wk1 at station 2, this returns
 * E[W_{n+1} at station 2 | Wk, Wk1] together with the mean interdeparture
 * time of station 1 and the probability that station 1 is idle when customer
 * n+1 arrives.
 *
 * The interdeparture time of an M/M/1 queue is a MIXTURE, not an exponential
 * of a single rate: with probability q = e^{-lambda Wk} mu1/(lambda+mu1) the
 * server empties before the next arrival, and the gap is then the arrival
 * time PLUS a service, a convolution of Exp(lambda) and Exp(mu1); otherwise
 * the gap is the service Exp(mu1) alone. J below is the corresponding
 * conditional Lindley step at station 2,
 *
 *   J(c, y, mu2) = E[max(y + S2 - X_c, 0)] with X_c ~ Exp(c), S2 ~ Exp(mu2),
 *
 * so the mean is (1-q) mu1 J(mu1,.) + q * (the convolution mixture). At
 * lambda == mu1 the two-term partial fraction of the convolution degenerates
 * and the gap becomes Erlang(2, mu1); that branch is taken from a relative
 * tolerance of 1e-9, which is MATLAB's, and evaluates the Erlang form Jw
 * directly rather than differencing two nearly equal terms.
 *
 * Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
 * Recursion", Proc. Winter Simulation Conference, 2021. Registered in
 * .citations() as 'tandemlindley'.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Mirrors the struct MATLAB returns from qsys_mm1_tandem_lindley. */
template <class T>
struct Mm1TandemLindleyResult {
    std::vector<T> mean;          ///< conditional mean waiting time at station 2
    std::vector<T> interdepMean;  ///< mean interdeparture time of station 1
    std::vector<T> idleProb;      ///< probability station 1 empties first
    std::string analyzer;
};

namespace detail {

/** E[max(y + Exp(mu2)^{-1} - Exp(c)^{-1}, 0)], the station-2 Lindley step. */
template <class T>
T tandem_j(const T& c, const T& y, const T& mu2) {
    const T one = num_traits<T>::from_int(1);
    const T e = num_exp(T(-c * y));
    return (y + one / mu2) * (one - e) / c - (one - e * (one + c * y)) / (c * c) +
           e / (mu2 * (c + mu2));
}

/** The same expectation when the gap is Erlang(2,c) rather than Exp(c). */
template <class T>
T tandem_jw(const T& c, const T& y, const T& mu2) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T e = num_exp(T(-c * y));
    const T d = c + mu2;
    return (y + one / mu2) * (one - e * (one + c * y)) / (c * c) -
           (two - e * (two + two * c * y + c * c * y * y)) / (c * c * c) +
           e * (y / d + one / (d * d)) / mu2;
}

}  // namespace detail

/**
 * @param lambda arrival rate at station 1, positive
 * @param mu1    service rate at station 1, positive
 * @param mu2    service rate at station 2, positive
 * @param Wk     waiting times at station 1, finite nonnegative
 * @param Wk1    waiting times at station 2, same length as Wk
 */
template <class T>
Mm1TandemLindleyResult<T> qsys_mm1_tandem_lindley(const T& lambda, const T& mu1, const T& mu2,
                                                  const std::vector<T>& Wk,
                                                  const std::vector<T>& Wk1) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mm1_tandem_lindley requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (lambda <= zero) throw InputError("qsys_mm1_tandem_lindley: lambda must be positive");
    if (mu1 <= zero) throw InputError("qsys_mm1_tandem_lindley: mu1 must be positive");
    if (mu2 <= zero) throw InputError("qsys_mm1_tandem_lindley: mu2 must be positive");
    if (Wk.size() != Wk1.size())
        throw InputError("qsys_mm1_tandem_lindley: Wk and Wk1 must have the same size");
    for (std::size_t i = 0; i < Wk.size(); ++i) {
        if (Wk[i] < zero) throw InputError("qsys_mm1_tandem_lindley: Wk must be nonnegative");
        if (Wk1[i] < zero) throw InputError("qsys_mm1_tandem_lindley: Wk1 must be nonnegative");
    }

    // MATLAB's relative test abs(lambda-mu1) > 1e-9*max(lambda,mu1)
    const T scale = (lambda > mu1) ? lambda : mu1;
    const bool distinct = num_abs(T(lambda - mu1)) > num_traits<T>::from_double(1e-9) * scale;

    Mm1TandemLindleyResult<T> r;
    r.analyzer = "qsys_mm1_tandem_lindley";
    r.mean.reserve(Wk.size());
    r.interdepMean.reserve(Wk.size());
    r.idleProb.reserve(Wk.size());
    for (std::size_t i = 0; i < Wk.size(); ++i) {
        const T q = detail::num_exp(T(-lambda * Wk[i])) * mu1 / (lambda + mu1);
        const T base = mu1 * detail::tandem_j(mu1, Wk1[i], mu2);
        T conv;
        if (distinct)
            conv = lambda * mu1 / (lambda - mu1) *
                   (detail::tandem_j(mu1, Wk1[i], mu2) - detail::tandem_j(lambda, Wk1[i], mu2));
        else
            conv = mu1 * mu1 * detail::tandem_jw(mu1, Wk1[i], mu2);
        r.mean.push_back((one - q) * base + q * conv);
        r.interdepMean.push_back(one / mu1 + q / lambda);
        r.idleProb.push_back(q);
    }
    return r;
}

/** Scalar overload. */
template <class T>
Mm1TandemLindleyResult<T> qsys_mm1_tandem_lindley(const T& lambda, const T& mu1, const T& mu2,
                                                  const T& Wk, const T& Wk1) {
    return qsys_mm1_tandem_lindley(lambda, mu1, mu2, std::vector<T>(1, Wk),
                                   std::vector<T>(1, Wk1));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MM1_TANDEM_LINDLEY_H
