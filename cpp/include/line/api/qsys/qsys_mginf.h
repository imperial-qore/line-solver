/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MGINF_H
#define LINE_API_QSYS_MGINF_H

/**
 * Exact solution of the M/G/infinity queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_mginf.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mginf.java. The JAR carries an extra
 * cv2 argument that it never reads, and returns a HashMap rather than a tuple;
 * the numbers are identical.
 *
 *   L = rho,  Lq = 0,  W = 1/mu,  Wq = 0,  p0 = exp(-rho)
 *   pk = exp(-rho) rho^k / k!
 *
 * The number in system is Poisson(rho), so p0 and pk carry exp(-rho). That is
 * a genuine transcendental: the function therefore requires
 * num_traits<T>::has_transcendental and cannot be instantiated at exact
 * arithmetic. L, Lq, W and Wq are field values, but they are not separable
 * from the struct, so the whole function is gated.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"

namespace line {
namespace qsys {

/** Return value of qsys_mginf, mirroring MATLAB's [L,Lq,W,Wq,p0,pk]. */
template <class T>
struct MginfResult {
    T L;         ///< mean number in system
    T Lq;        ///< mean number queueing, always 0
    T W;         ///< mean response time, always 1/mu
    T Wq;        ///< mean waiting time, always 0
    T p0;        ///< probability of an empty system
    T pk;        ///< probability of exactly k in system, valid iff has_pk
    bool has_pk; ///< true when the k overload was called
};

/**
 * @param lambda arrival rate
 * @param mu     service rate, mean service time 1/mu
 */
template <class T>
MginfResult<T> qsys_mginf(const T& lambda, const T& mu) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mginf requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    MginfResult<T> r;
    r.L = rho;
    r.Lq = zero;
    r.W = one / mu;
    r.Wq = zero;
    r.p0 = detail::num_exp(T(-rho));
    r.pk = zero;
    r.has_pk = false;
    return r;
}

/**
 * @param k state whose probability is additionally returned
 * @param lambda arrival rate
 * @param mu service rate
 */
template <class T>
MginfResult<T> qsys_mginf(const T& lambda, const T& mu, unsigned k) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mginf requires transcendental arithmetic");
    MginfResult<T> r = qsys_mginf(lambda, mu);
    const T rho = lambda / mu;
    r.pk = detail::num_exp(T(-rho)) * num_pow_int(rho, k) / num_factorial<T>(k);
    r.has_pk = true;
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MGINF_H
