/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DQSYS_DQSYS_BERNOULLI1_H
#define LINE_API_DQSYS_DQSYS_BERNOULLI1_H

/**
 * State dependent Bernoulli server on a discrete time scale.
 *
 * Templated port of matlab/src/api/dqsys/dqsys_bernoulli1.m. Time advances in
 * slots. In the slot starting at t with n jobs present the job in service
 * departs with probability p(n) and an arrival occurs with probability b(n),
 * independently; both are recorded at the end of the slot with the departure
 * resolved first (Daduna's LA rule and D/A rule). The queue length at slot
 * boundaries is a discrete birth-death chain with
 *
 *   pi(n) = [prod_{m=0}^{n-1} b(m) / prod_{m=0}^{n} c(m)]
 *         * [prod_{m=1}^{n-1} q(m) / prod_{m=1}^{n} p(m)] / H,
 *
 * c = 1-b and q = 1-p, which is theorem 2.3 of Daduna (2001), and corollary 2.8
 * once b(n) = 0 above the capacity. For constant b and p it collapses to the
 * Geo/Geo/1 law of `dqsys_geogeo1` under the LAS_DA convention.
 *
 * The law seen by an arriving customer, with himself not counted, is theorem
 * 2.11 and is returned in `arrivalPmf`. It is not the time-stationary law:
 * discrete time has no PASTA analogue, and the two differ even when the arrival
 * stream is a state independent Bernoulli process. In that state independent
 * case pi_1 is exactly the EAS-convention queue length law of `dqsys_geogeo1`,
 * geometric with ratio r = b(1-p)/(p(1-b)).
 *
 * Both laws are built by their exact product recurrences rather than in log
 * space, so the exact instantiation returns them with no rounding.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace dqsys {

/** Steady-state quantities of a finite-buffer Bernoulli server. */
template <class T>
struct Bernoulli1Result {
    std::size_t capacity = 0;    ///< buffer capacity in jobs
    std::vector<T> arrivalProb;  ///< offered b(n), n = 0..L
    std::vector<T> serviceProb;  ///< p(n), indexed by n-1
    std::vector<T> pmf;          ///< time-stationary law, theorem 2.3
    std::vector<T> arrivalPmf;   ///< arrival law, theorem 2.11
    T emptyProb;
    T utilization;               ///< 1 - pi(0)
    T throughput;                ///< carried departures per slot
    T lossProb;                  ///< fraction of offered arrivals lost
    T meanQueueLength;
    T meanWaitingQueue;
    T meanSojournTime;           ///< in slots, by Little's law
    T meanWaitingTime;           ///< in slots
    T normConst;                 ///< H of theorem 2.3
};

/**
 * Finite buffer of L jobs. An arrival in a slot that finds L jobs present is
 * lost, which is the loss system of corollary 2.8.
 *
 * @param b offered arrival probabilities b(n) for n = 0..L, or one entry for a
 *          state independent stream
 * @param p service probabilities p(n) for n = 1..L, or one entry for a state
 *          independent server
 * @param L buffer capacity in jobs
 */
template <class T>
Bernoulli1Result<T> dqsys_bernoulli1(const std::vector<T>& b, const std::vector<T>& p,
                                     std::size_t L) {
    if (L < 1) {
        throw InputError("dqsys_bernoulli1: L must be a positive integer");
    }
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<T> boff(L + 1);
    if (b.size() == 1) {
        for (std::size_t n = 0; n <= L; ++n) boff[n] = b[0];
    } else if (b.size() == L + 1) {
        boff = b;
    } else {
        throw InputError("dqsys_bernoulli1: the arrival probability vector must have one entry "
                         "per state 0..L");
    }
    std::vector<T> pv(L);
    if (p.size() == 1) {
        for (std::size_t n = 0; n < L; ++n) pv[n] = p[0];
    } else if (p.size() == L) {
        pv = p;
    } else if (p.size() == L + 1) {
        // A vector of length L+1 is accepted with its first entry, which would
        // be p(0), ignored.
        pv.assign(p.begin() + 1, p.end());
    } else {
        throw InputError("dqsys_bernoulli1: the service probability vector must have L or L+1 "
                         "entries");
    }
    for (std::size_t n = 0; n <= L; ++n) {
        if (boff[n] < zero || boff[n] > one) {
            throw InputError("dqsys_bernoulli1: arrival probabilities must lie in [0,1]");
        }
    }
    for (std::size_t n = 0; n < L; ++n) {
        if (!(pv[n] > zero) || pv[n] > one) {
            throw InputError("dqsys_bernoulli1: service probabilities must lie in (0,1]");
        }
    }
    std::vector<T> badm = boff;
    badm[L] = zero;                       // an arrival finding L jobs is lost
    for (std::size_t n = 0; n < L; ++n) {
        if (!(badm[n] < one)) {
            // c(n)=0 makes the weight of theorem 2.3 diverge at n; p(n)=1 is
            // fine and truncates the chain instead, which is example 2.9.
            throw InputError("dqsys_bernoulli1: arrival probabilities below the capacity must be "
                             "strictly less than one");
        }
    }

    // Theorem 2.3 by its exact recurrence, u(0) = 1/c(0) and
    // u(n) = u(n-1) b(n-1) q(n-1) / (c(n) p(n)) with q(0) := 1.
    std::vector<T> u(L + 1, zero);
    u[0] = one / (one - badm[0]);
    for (std::size_t n = 1; n <= L; ++n) {
        const T q = (n >= 2) ? (one - pv[n - 2]) : one;
        u[n] = u[n - 1] * badm[n - 1] * q / ((one - badm[n]) * pv[n - 1]);
    }
    T H = zero;
    for (std::size_t n = 0; n <= L; ++n) H = H + u[n];

    Bernoulli1Result<T> r;
    r.capacity = L;
    r.arrivalProb = boff;
    r.serviceProb = pv;
    r.pmf.resize(L + 1);
    for (std::size_t n = 0; n <= L; ++n) r.pmf[n] = u[n] / H;

    // Theorem 2.11 by its exact recurrence, v(0) = b(0)/(c(0) c(1)) and
    // v(n) = v(n-1) b(n) q(n) / (c(n+1) p(n)).
    std::vector<T> v(L, zero);
    v[0] = badm[0] / ((one - badm[0]) * (one - badm[1]));
    for (std::size_t n = 1; n < L; ++n) {
        v[n] = v[n - 1] * badm[n] * (one - pv[n - 1]) / ((one - badm[n + 1]) * pv[n - 1]);
    }
    T Ha = zero;
    for (std::size_t n = 0; n < L; ++n) Ha = Ha + v[n];
    r.arrivalPmf.resize(L);
    if (Ha > zero) {
        for (std::size_t n = 0; n < L; ++n) r.arrivalPmf[n] = v[n] / Ha;
    }

    r.emptyProb = r.pmf[0];
    r.utilization = one - r.pmf[0];
    T q = zero, t = zero, offered = zero;
    for (std::size_t n = 0; n <= L; ++n) {
        q = q + r.pmf[n] * num_traits<T>::from_int(static_cast<long>(n));
        offered = offered + r.pmf[n] * boff[n];
        if (n >= 1) t = t + r.pmf[n] * pv[n - 1];
    }
    r.meanQueueLength = q;
    r.throughput = t;
    r.lossProb = (offered > zero) ? (r.pmf[L] * boff[L] / offered) : zero;
    r.meanWaitingQueue = q - r.utilization;
    r.meanSojournTime = (t > zero) ? (q / t) : zero;
    r.meanWaitingTime = (t > zero) ? (r.meanWaitingQueue / t) : zero;
    r.normConst = H;
    return r;
}

}  // namespace dqsys
}  // namespace line

#endif  // LINE_API_DQSYS_DQSYS_BERNOULLI1_H
