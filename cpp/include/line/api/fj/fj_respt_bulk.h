/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RESPT_BULK_H
#define LINE_API_FJ_RESPT_BULK_H

/**
 * Centralized splitting analysed as an M[K]/M/c bulk arrival system.
 *
 * Templated port of matlab/src/api/fj/fj_respt_bulk.m.
 *
 * A request forks into K tasks held in a single central queue and served by c
 * identical servers, so the same server may serve several tasks of the same
 * request: an M[K]/M/c queue with fixed batch size K. Its level chain is solved
 * by truncation, which is exact up to the tail mass discarded.
 *
 * The request response time is the completion of the LAST of the K tasks. By
 * PASTA the batch finds n tasks in system, its last task is the (n+K)-th in
 * line, and under first come first served with c exponential servers it starts
 * service after max(0, n+K-c) departures, each an exponential of rate c mu:
 *
 *   E[R_request] = sum_n p_n [ max(0, n+K-c)/(c mu) + 1/mu ].
 *
 * This lower bounds the distributed splitting fork-join system, because no task
 * is bound to a particular server.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fj {

/** [Rreq, Rtask, Q, p] of fj_respt_bulk. */
template <class T>
struct FJResptBulkResult {
    T Rreq;
    T Rtask;
    T Q;
    std::vector<T> p;
};

/**
 * @param K      batch size, that is the number of tasks per request
 * @param lambda arrival rate of requests
 * @param mu     per-server task service rate
 * @param c      number of servers
 * @param nmax   truncation level of the task-count chain, 0 for the default
 * @return       the request and task response times, the mean queue and the distribution
 */
template <class T>
FJResptBulkResult<T> fj_respt_bulk(unsigned K, const T& lambda, const T& mu, unsigned c,
                                   std::size_t nmax = 0) {
    detail::require_positive_K(K, "fj_respt_bulk");
    if (c < 1) throw InputError("fj_respt_bulk: c must be a positive integer");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (!(lambda > zero) || !(mu > zero))
        throw InputError("fj_respt_bulk: lambda and mu must be positive");

    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    const T ct = num_traits<T>::from_int(static_cast<long>(c));
    const T rho = lambda * Kt / (ct * mu);
    if (rho >= one)
        throw NumericError("fj_respt_bulk: unstable system, rho = lambda*K/(c*mu) >= 1");

    if (nmax == 0) {
        const double rd = num_traits<T>::to_double(rho);
        const double want = static_cast<double>(K) + 40.0 * static_cast<double>(c) / (1.0 - rd);
        nmax = static_cast<std::size_t>(want > 200.0 ? want : 200.0);
    }
    const std::size_t ns = nmax + 1;

    // Generator of the task-count chain: batch arrivals of K, service min(n,c) mu
    Matrix<T> Q(ns, ns, zero);
    for (std::size_t i = 0; i < ns; ++i) {
        if (i > 0) {
            const std::size_t busy = (i < c) ? i : c;
            const T srv = num_traits<T>::from_int(static_cast<long>(busy)) * mu;
            Q(i, i - 1) = Q(i, i - 1) + srv;
            Q(i, i) = Q(i, i) - srv;
        }
        const std::size_t j = i + K;
        if (j < ns) {
            Q(i, j) = Q(i, j) + lambda;
            Q(i, i) = Q(i, i) - lambda;
        }
    }

    FJResptBulkResult<T> out;
    out.p = mc::ctmc_solve(Q);

    out.Q = zero;
    for (std::size_t n = 0; n < ns; ++n)
        out.Q += num_traits<T>::from_int(static_cast<long>(n)) * out.p[n];
    out.Rtask = out.Q / (lambda * Kt);

    // Last of the K tasks of a tagged batch: the (n+K)-th in line on arrival
    out.Rreq = zero;
    for (std::size_t n = 0; n < ns; ++n) {
        const std::size_t ahead = n + K;
        const T wait = (ahead > c) ? num_traits<T>::from_int(static_cast<long>(ahead - c)) /
                                         (ct * mu)
                                   : zero;
        out.Rreq += out.p[n] * (wait + one / mu);
    }
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RESPT_BULK_H
