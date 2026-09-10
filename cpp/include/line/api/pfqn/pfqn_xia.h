/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_XIA_H
#define LINE_API_PFQN_XIA_H

/**
 * Xia's asymptotic approximation of the normalizing constant of a
 * load-dependent (multiserver) closed network.
 *
 * Templated port of jar/src/main/java/jline/api/pfqn/ld/Pfqn_xia.java. MATLAB
 * has no counterpart.
 *
 * The demands are first rescaled so that the largest per-server utilization
 * rho_i = L_i / s_i is one. The stations that attain it are the BOTTLENECK SET
 * B; they saturate and contribute the M/M/s saturated term, while every other
 * station contributes its finite-capacity Erlang-like partial sum
 *
 *   F(u, k) = sum_{j < k} u^j / j!  +  (u^k / k!) / (1 - u/k),
 *
 * the closed form of the geometric tail beyond the k-th server. The result is
 *
 *   log G ~ -log((B-1)!) - N log(c) + sum_{b in B} [ s_b log L_b - log(s_b!) ]
 *                                   + sum_{k not in B} log F(L_k, s_k),
 *
 * c being the rescaling factor. Note the leading behaviour in N enters ONLY
 * through -N log c: this is the large-population limit, so the approximation
 * does not resolve the O(1) corrections that a finite population carries.
 *
 * A NON-BOTTLENECK STATION WITH u > k GIVES A NEGATIVE F, whose logarithm is
 * NaN and poisons the whole constant. The reference guards only against an
 * INFINITE F (u == k exactly), not a negative one, so it propagates the NaN;
 * that is reproduced rather than patched, because suppressing the term would
 * quietly return a plausible number for a model the expansion does not cover.
 * The condition cannot arise when every station has one server, but it can once
 * s varies. Note that at Real precision Boost raises on log of a negative where
 * double returns NaN, so the same input surfaces as an exception rather than a
 * NaN; both are refusals and neither is a silent wrong number.
 *
 * Arithmetic: TRANSCENDENTAL. Logs and factorials throughout.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** F(u,k) = sum_{j<k} u^j/j! + (u^k/k!)/(1 - u/k). */
template <class T>
T pfqn_xia_F(const T& u, const T& k) {
    const T one = num_traits<T>::from_int(1);
    const long kk = static_cast<long>(num_traits<T>::to_double(k));
    T ret = num_traits<T>::from_int(0);
    T upow = one;
    T fact = one;
    for (long j = 0; j < kk; ++j) {
        ret += T(upow / fact);
        upow *= u;
        fact *= num_traits<T>::from_int(j + 1);
    }
    // the tail term: u^k / k! / (1 - u/k)
    const T denom = T(one - T(u / k));
    return T(ret + T(T(upow / fact) / denom));
}

}  // namespace detail

/**
 * @param L (M) service demands
 * @param N population
 * @param s (M) server counts
 * @return log of the approximate normalizing constant
 */
template <class T>
T pfqn_xia(const std::vector<T>& L, int N, const std::vector<T>& s) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_xia is an asymptotic expansion in logs and needs transcendental "
                  "arithmetic");
    using std::log;
    const std::size_t M = L.size();
    if (M == 0) throw InputError("pfqn_xia requires at least one station");
    if (s.size() != M) throw InputError("pfqn_xia: L and s disagree on the station count");
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < M; ++i) {
        if (L[i] <= zero) throw InputError("pfqn_xia requires positive demands");
        if (s[i] <= zero) throw InputError("pfqn_xia requires positive server counts");
    }

    std::vector<T> rho(M);
    for (std::size_t i = 0; i < M; ++i) rho[i] = T(L[i] / s[i]);
    T rmax = rho[0];
    for (std::size_t i = 1; i < M; ++i)
        if (rho[i] > rmax) rmax = rho[i];
    const T scalefactor = T(num_traits<T>::from_int(1) / rmax);
    std::vector<T> Ls(M), rs(M);
    for (std::size_t i = 0; i < M; ++i) {
        Ls[i] = T(L[i] * scalefactor);
        rs[i] = T(rho[i] * scalefactor);
    }
    T rsmax = rs[0];
    for (std::size_t i = 1; i < M; ++i)
        if (rs[i] > rsmax) rsmax = rs[i];

    std::vector<std::size_t> bnk, nbnk;
    for (std::size_t i = 0; i < M; ++i) {
        if (rs[i] == rsmax)
            bnk.push_back(i);
        else
            nbnk.push_back(i);
    }
    const std::size_t B = bnk.size();

    T logGasy = T(-detail::num_factln<T>(num_traits<T>::from_int(static_cast<long>(B) - 1)) -
                  num_traits<T>::from_int(N) * log(scalefactor));
    for (std::size_t b = 0; b < B; ++b) {
        const std::size_t i = bnk[b];
        logGasy += s[i] * log(Ls[i]) - detail::num_factln<T>(s[i]);
    }
    for (std::size_t k = 0; k < nbnk.size(); ++k) {
        const std::size_t i = nbnk[k];
        const T f = detail::pfqn_xia_F(Ls[i], s[i]);
        const double fd = num_traits<T>::to_double(f);
        if (std::isfinite(fd)) logGasy += log(f);
    }
    return logGasy;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_XIA_H
