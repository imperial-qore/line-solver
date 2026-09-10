/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_RGF_H
#define LINE_API_PFQN_PFQN_RGF_H

/**
 * Recursion by Generating Functions (RGF) for the normalizing constant of a
 * SINGLE-CLASS closed product-form network with replicated stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_rgf.m, i.e. Property 1 of
 * J. Coury and P. G. Harrison, "Asymptotic properties of queuing networks",
 * IEE Proc.-Comput. Digit. Tech. 144(5):247-254, 1997.
 *
 * WHAT MAKES IT DIFFERENT FROM BUZEN. Convolution proceeds per GENERATING
 * FUNCTION rather than per station: a group of m stations sharing one demand p
 * collapses into the single negative-binomial sequence r(k) = C(k+m-1,k) p^k,
 * so the group costs one convolution pass instead of m. The delay contributes
 * the Poisson sequence r(k) = Z^k / k!. Convolving the G distinct sequences
 * yields g(0..N) exactly, and lG = log g(N).
 *
 * Cost O(G N^2) against Buzen's O(M N), so RGF is the cheaper route precisely
 * when the model is heavily replicated and the population moderate (G N < M).
 *
 * ARITHMETIC. The whole recursion runs in the LOG domain -- that is the point
 * of the routine, since no intermediate can then overflow or underflow -- so it
 * is gated on num_traits<T>::has_transcendental. A caller that wants the same
 * constant in exact arithmetic wants pfqn_ca, which computes it by the ordinary
 * convolution over the rationals.
 *
 * SINGLE CLASS BY CONSTRUCTION. Grouping stations by their demand only defines
 * a sequence when the demand is a scalar, so a multiclass argument is refused
 * here by name; `pfqn_nc`'s 'rgf' branch sends a multiclass model to pfqn_ca
 * instead, which is what the reference does.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_rgf, mirroring [G, lG, lg]. */
template <class T>
struct RgfResult {
    T G;                ///< normalizing constant
    T lG;               ///< its logarithm, i.e. lg[N]
    std::vector<T> lg;  ///< log g(0), log g(1), ..., log g(N)
};

namespace detail {

/** Log-domain linear convolution truncated at the common length (logconv.m). */
template <class T>
std::vector<T> rgf_logconv(const std::vector<T>& u, const std::vector<T>& v) {
    const std::size_t n = u.size();
    std::vector<T> c(n, num_traits<T>::from_double(-std::numeric_limits<double>::infinity()));
    for (std::size_t k = 0; k < n; ++k) {
        std::vector<T> t(k + 1);
        for (std::size_t i = 0; i <= k; ++i) t[i] = T(u[i] + v[k - i]);
        // logsumexp returns the maximum unchanged when every term is -inf,
        // which is the reference's `if isinf(vm), c(k) = vm` branch.
        c[k] = logsumexp(t);
    }
    return c;
}

}  // namespace detail

/**
 * @param L (M) service demands of the queueing stations
 * @param N population, a nonnegative integer
 * @param Z aggregate delay demand (think time); zero for none
 */
template <class T>
RgfResult<T> pfqn_rgf(const std::vector<T>& L, int N, const T& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_rgf requires transcendental arithmetic: the recursion is carried in the "
                  "log domain so that no intermediate can overflow. Use pfqn_ca for the same "
                  "constant in exact arithmetic");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T ninf = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());

    if (N < 0) throw InputError("pfqn_rgf: requires a nonnegative integer population");
    if (Z < zero) throw InputError("pfqn_rgf: requires a nonnegative think time");
    for (const T& v : L)
        if (v < zero) throw InputError("pfqn_rgf: requires nonnegative demands");

    const std::size_t Np = static_cast<std::size_t>(N);
    RgfResult<T> res;
    // g(k) = 1 for the empty network, 0 elsewhere before any node is folded in.
    res.lg.assign(Np + 1, ninf);
    res.lg[0] = zero;

    // Delay node: the Poisson sequence Z^k / k!.
    if (Z > zero) {
        std::vector<T> lr(Np + 1);
        for (std::size_t k = 0; k <= Np; ++k) {
            const T kT = num_traits<T>::from_int(static_cast<long>(k));
            lr[k] = T(kT * log(Z) - detail::num_factln<T>(kT));
        }
        res.lg = detail::rgf_logconv(res.lg, lr);
    }

    // Queueing stations grouped by identical demand: one sequence per group.
    // The reference groups through unique(L), which sorts; the multiplicity is
    // all that enters the sequence, so sorting a copy reproduces it exactly.
    std::vector<T> pos;
    for (const T& v : L)
        if (v > zero) pos.push_back(v);
    std::sort(pos.begin(), pos.end(), [](const T& a, const T& b) { return a < b; });
    std::size_t i = 0;
    while (i < pos.size()) {
        std::size_t j = i;
        while (j < pos.size() && pos[j] == pos[i]) ++j;
        const std::size_t m = j - i;
        const T p = pos[i];
        std::vector<T> lr(Np + 1);
        for (std::size_t k = 0; k <= Np; ++k) {
            const T kT = num_traits<T>::from_int(static_cast<long>(k));
            if (m == 1) {
                // 1 / (1 - p u)
                lr[k] = T(kT * log(p));
            } else {
                const T mT = num_traits<T>::from_int(static_cast<long>(m));
                lr[k] = T(detail::num_lgamma<T>(T(kT + mT)) - detail::num_factln<T>(kT) -
                          detail::num_lgamma<T>(mT) + kT * log(p));
            }
        }
        res.lg = detail::rgf_logconv(res.lg, lr);
        i = j;
    }

    res.lG = res.lg[Np];
    res.G = exp(res.lG);
    return res;
}

/** Overload without a delay. */
template <class T>
RgfResult<T> pfqn_rgf(const std::vector<T>& L, int N) {
    return pfqn_rgf(L, N, num_traits<T>::from_int(0));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_RGF_H
