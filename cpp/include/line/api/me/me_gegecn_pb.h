/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_ME_ME_GEGECN_PB_H
#define LINE_API_ME_ME_GEGECN_PB_H

/**
 * Blocking probability seen by ONE arrival stream of a censored GE/GE/c/K;N
 * queue.
 *
 * Templated port of `matlab/src/api/me/me_gegecn_pb.m`, equation (4.3) of
 * Kouvatsos (1994), evaluated on the queue-length distribution `me_gegecn`
 * returns.
 *
 * WHY AN ARRIVAL CAN BE BLOCKED WITH ROOM TO SPARE. A GE arrival process is a
 * BATCH process, so a batch arriving to a queue holding n < N jobs can still
 * overflow the residual room. `(1-tau)^(N-n)` is the probability that it does,
 * and the first sum carries an extra factor for the servers still idle. With a
 * Poisson stream (Ca = 1, so tau = 1) every term but n = N vanishes and this
 * collapses to the PASTA value p(N) -- which is the cheapest check that the
 * formula is being evaluated correctly.
 *
 * THE STREAM SCV IS PER STREAM, NOT PER NODE, and that is the whole point of
 * the function existing separately. One node solution `p` yields a DIFFERENT
 * blocking probability for each flow merging into the queue -- the external
 * arrivals, the flow from each upstream station, and the flow released by each
 * holding node -- which is exactly how PBe_j, PB^i_j and PB^{h_ij}_j are
 * obtained in the transfer-blocking algorithm of Tahilramani, Manjunath and
 * Bose (1999).
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace me {

/**
 * Port of `me_gegecn_pb`.
 *
 * @param p  queue-length distribution, `p[idx] = Pr{n = K + idx}`
 * @param K  minimum number of jobs in the queue
 * @param N  buffer capacity in jobs
 * @param c  number of servers
 * @param Cs squared coefficient of variation of the service times
 * @param Ca squared coefficient of variation of the interarrival times OF THE
 *           STREAM whose blocking probability is requested
 * @return the probability that an arrival of this stream finds the queue full
 */
template <class T>
T me_gegecn_pb(const std::vector<T>& p, long K, long N, long c, const T& Cs, const T& Ca) {
    static_assert(num_traits<T>::has_transcendental,
                  "me_gegecn_pb requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (N <= K) throw InputError("me_gegecn_pb: the capacity must exceed the minimum occupancy");
    const std::size_t nn = static_cast<std::size_t>(N - K + 1);
    if (p.size() != nn)
        throw InputError("me_gegecn_pb: the distribution must hold N-K+1 entries");

    const T tau = T(two / (Ca + one));
    const T sigma = T(two / (Cs + one));
    const T omtau = T(one - tau);

    // w(idx) = (1-tau)^(N-n), the probability that a batch overflows the room.
    std::vector<T> w(nn, one);
    for (std::size_t idx = 0; idx < nn; ++idx) {
        const long n = K + static_cast<long>(idx);
        w[idx] = num_pow_int(omtau, static_cast<unsigned>(N - n));
    }

    T PB = zero;
    // Jobs arriving while some servers are still idle: n = K,...,c-1
    if (K < c) {
        const std::size_t last =
            static_cast<std::size_t>(std::min<long>(c - K, static_cast<long>(nn)));
        const T den = T(sigma * omtau + tau);
        for (std::size_t idx = 0; idx < last; ++idx) {
            const long n = K + static_cast<long>(idx);
            const T fac = num_pow_int(T(sigma / den), static_cast<unsigned>(c - n));
            PB += T(w[idx] * fac * p[idx]);
        }
    }
    // Jobs arriving with all servers busy: n = max(c,K),...,N
    const long lo = std::max(c, K);
    if (lo <= N)
        for (std::size_t idx = static_cast<std::size_t>(lo - K); idx < nn; ++idx)
            PB += T(w[idx] * p[idx]);
    return PB;
}

}  // namespace me
}  // namespace line

#endif  // LINE_API_ME_ME_GEGECN_PB_H
