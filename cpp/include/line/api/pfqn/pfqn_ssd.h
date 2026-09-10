/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_SSD_H
#define LINE_API_PFQN_PFQN_SSD_H

/**
 * Server-Station Disaggregation bounds for a multiserver closed network
 * (Dallery and Suri, SIGMETRICS 1986).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ssd.m. Single-class model: L is the
 * per-station demand vector, N the population, Z the think time.
 *
 * Theorem 5 eq. (6) carries no think time of its own. The objection this file
 * used to raise still stands as far as it went: inserting a bare +Z into the
 * Theorem 5 form is NOT a bound, because the BJB optimistic step needs
 * sum_k Q_k(N-1) = N-1, which fails once Z X(N-1) jobs sit at the terminal.
 * What that argument missed is that the queueing term can be corrected instead
 * of the delay being disaggregated separately. The reference now scales it by
 * the terminal-workload factor of Lazowska et al. 1984, Table 5.2:
 * (N-1) Y_l / (1 + Z/(N R_l)) on the lower bound and (N-1) Y_u / (1 + Z/R_u) on
 * the upper. Both reduce to the Z=0 forms exactly, so this generalises rather
 * than replaces eq. (6), and the Z>0 case is computed instead of refused.
 *
 * All operations stay in the field, so the bound is exact in rational
 * arithmetic: a bound computed exactly is worth having, since a bound violated
 * only by rounding is indistinguishable from a real violation.
 */

#include <algorithm>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
struct SsdBounds {
    T Xlo;
    T Xhi;
};

/**
 * @param L demands, @param N population, @param Z think time,
 * @param nservers per-station server counts (empty for all ones)
 */
template <class T>
SsdBounds<T> pfqn_ssd(const std::vector<T>& L, const T& N, const T& Z,
                      const std::vector<T>& nservers) {
    const std::size_t K = L.size();
    if (K == 0) throw InputError("pfqn_ssd: empty demand vector");
    if (!nservers.empty() && nservers.size() != K)
        throw InputError("pfqn_ssd: server-count vector has the wrong length");
    const T one = num_traits<T>::from_int(1);

    T Rl = num_traits<T>::from_int(0), Ru = num_traits<T>::from_int(0);
    T Yl = num_traits<T>::from_int(0);
    std::size_t b = 0;
    for (std::size_t i = 0; i < K; ++i) {
        const T c = nservers.empty() ? one : nservers[i];
        if (c == num_traits<T>::from_int(0)) throw InputError("pfqn_ssd: zero server count");
        Rl += L[i];
        const T lc = L[i] / c;
        Ru += lc;
        if (lc > Yl) {
            Yl = lc;
            b = i;
        }
    }
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    const T Yu = Ru / Kt;
    const T cb = nservers.empty() ? one : nservers[b];

    const T zero = num_traits<T>::from_int(0);
    // Lazowska Table 5.2 terminal-workload scaling; at Z=0 both factors are 1.
    // A zero total demand forces Yl = Yu = 0, so the term vanishes and the
    // division that would be undefined is never reached.
    const T ql = (Rl == zero) ? zero : (N - one) * Yl / (one + Z / (N * Rl));
    const T qu = (Ru == zero) ? zero : (N - one) * Yu / (one + Z / Ru);

    SsdBounds<T> r;
    r.Xlo = N / (Rl + Z + ql);  // Theorem 5 lower, think-time corrected
    T hi = N / (Ru + Z + qu);   // Theorem 5 upper, eq. (6), think-time corrected
    const T cap = cb / L[b];
    if (cap < hi) hi = cap;
    const T pop = N / (Rl + Z);
    if (pop < hi) hi = pop;
    r.Xhi = hi;
    return r;
}

template <class T>
SsdBounds<T> pfqn_ssd(const std::vector<T>& L, const T& N, const T& Z) {
    return pfqn_ssd(L, N, Z, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif
