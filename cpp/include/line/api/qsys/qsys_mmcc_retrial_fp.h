/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MMCC_RETRIAL_FP_H
#define LINE_API_QSYS_QSYS_MMCC_RETRIAL_FP_H

/**
 * Fixed-point approximation for the M/M/c/c retrial queue.
 *
 * Port of matlab/src/api/qsys/qsys_mmcc_retrial_fp.m. Blocked customers join an
 * orbit and retry; under the assumption that the retrial rate is small relative
 * to the service rate the superposition of fresh and retrial arrivals is
 * approximated by a Poisson stream of rate lambda + r, with r the solution of
 *
 *     r = (lambda + r) B((lambda + r)/mu, c),
 *
 * B being Erlang's loss formula. Cohen (1957); Phung-Duc, "Retrial Queueing
 * Models: A Survey on Theory and Applications" (2019), eq. (1).
 *
 * The reference evaluates B by the rational recursion
 * B_k = a B_{k-1}/(k + a B_{k-1}), which is exactly what line::lossn::erlang_b
 * declines to do (it goes through log/exp and says so). This port keeps the
 * reference's recursion, so B itself is a rational function of a and needs no
 * transcendental arithmetic. The static_assert is nonetheless present because
 * the OUTER iteration is a fixed point tested against a tolerance: it converges
 * geometrically but not in a finite number of field operations, and at Rational
 * the iterates would grow without bound in representation size while never
 * meeting an exact stopping rule.
 *
 * The iteration is a monotone increasing map of r started from 0 and bounded
 * above by lambda/(1-B) at the fixed point, so it converges from below and the
 * successive-difference test is a genuine stopping criterion rather than a
 * heuristic. The reference does NOT report non-convergence; when maxiter is
 * exhausted it silently returns the last iterate, and the port reports the
 * iteration count so the caller can tell the two apart, exactly as the third
 * MATLAB output niter does.
 */

#include <cstddef>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Return value of qsys_mmcc_retrial_fp, mirroring the three MATLAB outputs. */
template <class T>
struct MmccRetrialFpResult {
    T blockingProbability;  ///< B((lambda + r)/mu, c)
    T retrialRate;          ///< r, the extra arrival rate contributed by the orbit
    std::size_t iterations;  ///< iterations performed, = maxiter when unconverged
    bool converged;          ///< whether the successive-difference test was met
};

/**
 * Erlang's loss formula B(a, c) by the numerically stable rational recursion
 * B_0 = 1, B_k = a B_{k-1}/(k + a B_{k-1}). Exact in any field: no logs, no
 * factorials, no cancellation.
 */
template <class T>
T erlang_b_recursive(const T& a, unsigned c) {
    T b = num_traits<T>::from_int(1);
    for (unsigned i = 1; i <= c; ++i) b = a * b / (num_traits<T>::from_int(static_cast<long>(i)) + a * b);
    return b;
}

/**
 * M/M/c/c with retrials by the Cohen fixed point.
 *
 * @param lambda  fresh arrival rate
 * @param mu      service rate of one server
 * @param c       number of servers, which is also the capacity
 * @param tol     stopping tolerance on |r_{k+1} - r_k|
 * @param maxiter iteration cap
 */
template <class T>
MmccRetrialFpResult<T> qsys_mmcc_retrial_fp(const T& lambda, const T& mu, unsigned c, const T& tol,
                                            std::size_t maxiter) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mmcc_retrial_fp is a tolerance-terminated fixed point");
    const T zero = num_traits<T>::from_int(0);
    if (lambda <= zero) throw InputError("qsys_mmcc_retrial_fp: arrival rate must be positive");
    if (mu <= zero) throw InputError("qsys_mmcc_retrial_fp: service rate must be positive");
    if (c == 0) throw InputError("qsys_mmcc_retrial_fp: at least one server is required");
    if (tol <= zero) throw InputError("qsys_mmcc_retrial_fp: tolerance must be positive");
    if (maxiter == 0) throw InputError("qsys_mmcc_retrial_fp: maxiter must be positive");

    T r = zero;
    std::size_t it = 0;
    bool converged = false;
    for (; it < maxiter;) {
        ++it;
        const T a = (lambda + r) / mu;
        const T b = erlang_b_recursive(a, c);
        const T rnew = (lambda + r) * b;
        const T gap = num_abs(T(rnew - r));
        r = rnew;
        if (gap < tol) {
            converged = true;
            break;
        }
    }

    MmccRetrialFpResult<T> out;
    out.retrialRate = r;
    out.blockingProbability = erlang_b_recursive(T((lambda + r) / mu), c);
    out.iterations = it;
    out.converged = converged;
    return out;
}

/** qsys_mmcc_retrial_fp with the reference defaults tol = 1e-10, maxiter = 10000. */
template <class T>
MmccRetrialFpResult<T> qsys_mmcc_retrial_fp(const T& lambda, const T& mu, unsigned c) {
    return qsys_mmcc_retrial_fp(lambda, mu, c, T(num_traits<T>::from_double(1e-10)),
                                static_cast<std::size_t>(10000));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MMCC_RETRIAL_FP_H
