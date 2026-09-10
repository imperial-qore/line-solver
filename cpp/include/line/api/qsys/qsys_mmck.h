/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MMCK_H
#define LINE_API_QSYS_QSYS_MMCK_H

/**
 * Exact analysis of the M/M/c/K queue (truncated Erlang form).
 *
 * Templated port of matlab/src/api/qsys/qsys_mmck.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mmck.java.
 *
 *   a = lambda/mu,  rho = a/c
 *   p_n = a^n/n! p_0                    0 <= n <= c
 *   p_n = a^c/c! rho^(n-c) p_0          c <= n <= K
 *   p_0 = 1 / [ sum_{n<c} a^n/n! + (a^c/c!) sum_{n=0}^{K-c} rho^n ]
 *
 * All exponents are integers and the normalization is a finite sum, so the
 * whole computation stays in the field of the inputs and is exact for
 * T = Rational. The exact instantiation is not academic here: MATLAB builds
 * the unnormalized vector as a^n/n! and warns about overflow, which is the
 * failure mode for large a and K; in rational arithmetic the intermediate
 * magnitudes are irrelevant.
 *
 * Unlike the unbounded M/M/c, no stability condition is needed: a finite
 * capacity makes every load admissible, rho >= 1 included.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

template <class T>
struct MmckResult {
    T meanQueueLength;    ///< L, mean number in system
    T meanQueueLengthQ;   ///< Lq, mean number waiting
    T meanWaitingTime;    ///< Wq = Lq/lambda_eff
    T meanSojournTime;    ///< W  = L /lambda_eff
    T utilization;        ///< per-server utilization lambda_eff/(c mu)
    T throughput;         ///< lambda_eff = lambda (1 - p_K)
    T lossProbability;    ///< p_K
    std::vector<T> queueLengthDist;  ///< p_0 ... p_K
};

/**
 * @param lambda Poisson arrival rate
 * @param mu     exponential service rate of one server
 * @param c      number of servers, c >= 1
 * @param K      system capacity, K >= c
 */
template <class T>
MmckResult<T> qsys_mmck(const T& lambda, const T& mu, unsigned c, unsigned K) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (lambda <= zero) throw InputError("qsys_mmck: lambda must be a positive real scalar");
    if (mu <= zero) throw InputError("qsys_mmck: mu must be a positive real scalar");
    if (c < 1) throw InputError("qsys_mmck: c must be a positive integer");
    if (K < c) throw InputError("qsys_mmck: K must be an integer >= c");

    const T a = lambda / mu;
    const T rho = a / num_traits<T>::from_int(static_cast<long>(c));

    std::vector<T> p(K + 1, zero);
    for (unsigned n = 0; n < c; ++n) p[n] = num_pow_int(a, n) / num_factorial<T>(n);
    const T ac_over_cfact = num_pow_int(a, c) / num_factorial<T>(c);
    for (unsigned n = c; n <= K; ++n) p[n] = ac_over_cfact * num_pow_int(rho, n - c);

    T S = zero;
    for (const T& v : p) S += v;
    if (S <= zero) throw NumericError("qsys_mmck: stationary distribution failed to normalize");
    for (T& v : p) v /= S;

    MmckResult<T> r;
    r.queueLengthDist = p;
    T L = zero, Lq = zero;
    for (unsigned n = 0; n <= K; ++n) {
        L += num_traits<T>::from_int(static_cast<long>(n)) * p[n];
        if (n > c) Lq += num_traits<T>::from_int(static_cast<long>(n - c)) * p[n];
    }
    r.meanQueueLength = L;
    r.meanQueueLengthQ = Lq;
    r.lossProbability = p[K];
    const T lambdaEff = lambda * (one - p[K]);
    r.throughput = lambdaEff;
    r.utilization = lambdaEff / (num_traits<T>::from_int(static_cast<long>(c)) * mu);
    if (lambdaEff > zero) {
        r.meanWaitingTime = Lq / lambdaEff;
        r.meanSojournTime = L / lambdaEff;
    } else {
        r.meanWaitingTime = zero;
        r.meanSojournTime = zero;
    }
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MMCK_H
