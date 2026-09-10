/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_PBH_H
#define LINE_API_PFQN_PFQN_PBH_H

/**
 * Performance Bound Hierarchy (Eager and Sevcik 1983, ACM TOCS 1(2):99-115)
 * for single-class closed product-form networks, and the two iterative
 * families that are defined in terms of it.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_pbh.m, pfqn_pbk.m and
 * pfqn_bjbk.m. `level` MVA steps from an ABA-initialized residence give nested
 * optimistic and pessimistic bounds that converge to exact MVA as level -> N.
 *
 * MATLAB REDUNDANCY, reproduced rather than hidden: pfqn_pbk.m and
 * pfqn_bjbk.m are both one-line forwarders to pfqn_pbh with the same
 * arguments, so PB(k) and BJB(k) return identical numbers for every input.
 * The two names are kept because the surrounding solver code refers to both,
 * and collapsing them here would hide the fact in the port.
 *
 * ARITHMETIC. The recursion is a finite sequence of field operations -- no
 * root, no logarithm anywhere -- so these bounds are EXACT in rational
 * arithmetic and are deliberately left ungated. That is worth having: a bound
 * violated only by rounding cannot be told apart from a real violation.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_pbh, mirroring [Xlo, Xhi, Qlo, Qhi]. */
template <class T>
struct PbhBounds {
    T Xlo;
    T Xhi;
    std::vector<T> Qlo;
    std::vector<T> Qhi;
};

namespace detail {

/** Per-station residence vector of the level-`level` bound, one side. */
template <class T>
std::vector<T> pbh_residence(const std::vector<T>& L, int N, const T& Z, int level, bool optimistic) {
    const std::size_t K = L.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::size_t b = 0;
    for (std::size_t i = 1; i < K; ++i)
        if (L[i] > L[b]) b = i;
    const int lv = std::min(level, N);
    const int n0 = N - lv;
    T Lsum = zero;
    for (const T& x : L) Lsum += x;

    std::vector<T> Rk(K, zero);
    if (optimistic) {
        T asym = T(num_traits<T>::from_int(n0) * L[b] - Z);
        if (asym < Lsum) asym = Lsum;
        const T v = T(asym / num_traits<T>::from_int(static_cast<long>(K)));
        Rk.assign(K, v);
    } else {
        Rk.assign(K, zero);
        Rk[b] = num_traits<T>::from_int(n0);
    }
    if (n0 == 0) Rk.assign(K, zero);

    for (int n = n0 + 1; n <= N; ++n) {
        T Rtot = zero;
        for (const T& x : Rk) Rtot += x;
        if (n == 1 || T(Z + Rtot) == zero) {
            Rk = L;
        } else {
            const T f = T(num_traits<T>::from_int(n - 1) / T(Z + Rtot));
            for (std::size_t i = 0; i < K; ++i) Rk[i] = T(L[i] * T(one + f * Rk[i]));
        }
    }
    return Rk;
}

}  // namespace detail

/**
 * @param L     (M) per-station demands
 * @param N     population
 * @param Z     think time
 * @param level hierarchy level >= 0, clamped to N
 */
template <class T>
PbhBounds<T> pfqn_pbh(const std::vector<T>& L, int N, const T& Z, int level) {
    const std::size_t K = L.size();
    if (K == 0) throw InputError("pfqn_pbh: empty demand vector");
    if (N < 0) throw InputError("pfqn_pbh: negative population");
    if (level < 0) throw InputError("pfqn_pbh: negative hierarchy level");
    const T zero = num_traits<T>::from_int(0);
    T Lmax = L[0], Lsum = zero;
    for (const T& x : L) {
        Lsum += x;
        if (x > Lmax) Lmax = x;
    }

    const std::vector<T> Ro = detail::pbh_residence(L, N, Z, level, true);
    const std::vector<T> Rp = detail::pbh_residence(L, N, Z, level, false);

    T Rosum = zero, Rpsum = zero;
    for (const T& x : Ro) Rosum += x;
    for (const T& x : Rp) Rpsum += x;

    // Joint with the asymptotic residence lower bound R(N) >= max(N Lmax - Z, sum L).
    T RoC = Rosum;
    T asym = T(num_traits<T>::from_int(N) * Lmax - Z);
    if (asym < Lsum) asym = Lsum;
    if (RoC < asym) RoC = asym;

    PbhBounds<T> r;
    if (Lmax == zero) throw NumericError("pfqn_pbh: all demands are zero");
    r.Xhi = T(num_traits<T>::from_int(1) / Lmax);
    const T alt = T(num_traits<T>::from_int(N) / T(Z + RoC));
    if (alt < r.Xhi) r.Xhi = alt;
    r.Xlo = T(num_traits<T>::from_int(N) / T(Z + Rpsum));

    r.Qlo.resize(K);
    r.Qhi.resize(K);
    for (std::size_t i = 0; i < K; ++i) {
        r.Qlo[i] = T(r.Xlo * Ro[i]);
        r.Qhi[i] = T(r.Xhi * Rp[i]);
    }
    return r;
}

template <class T>
PbhBounds<T> pfqn_pbh(const std::vector<T>& L, int N, const T& Z) {
    return pfqn_pbh(L, N, Z, 1);
}

/** PB(k), the iterative Eager-Sevcik proportional bound. Forwards to pfqn_pbh. */
template <class T>
PbhBounds<T> pfqn_pbk(const std::vector<T>& L, int N, const T& Z, int k) {
    return pfqn_pbh(L, N, Z, k);
}

template <class T>
PbhBounds<T> pfqn_pbk(const std::vector<T>& L, int N, const T& Z) {
    return pfqn_pbh(L, N, Z, 1);
}

/** BJB(k), the iterative Balanced Job Bound. Forwards to pfqn_pbh. */
template <class T>
PbhBounds<T> pfqn_bjbk(const std::vector<T>& L, int N, const T& Z, int k) {
    return pfqn_pbh(L, N, Z, k);
}

template <class T>
PbhBounds<T> pfqn_bjbk(const std::vector<T>& L, int N, const T& Z) {
    return pfqn_pbh(L, N, Z, 1);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_PBH_H
