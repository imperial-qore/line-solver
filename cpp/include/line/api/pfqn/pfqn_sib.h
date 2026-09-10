/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_SIB_H
#define LINE_API_PFQN_PFQN_SIB_H

/**
 * Successively Improving Bounds (Srinivasan 1985/1987) on the cycle time and
 * throughput of a single-class closed product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sib.m, including its local
 * functions phi_u1, sigma and betaL. The hierarchy bounds
 * phi(K) = sum_m rho_m Q_m(K) from both sides, with level 1 the closed form of
 * Theorem 2.1 and higher levels the S_i power-sum forms of Theorems 3.5 and
 * 3.6, and reads the cycle time off W(N) = sum(L) (1 + phi(N-1)).
 *
 * REFERENCE DEFECT ABOVE LEVEL 1, and the port does NOT reproduce it.
 * pfqn_sib.m declares phi_u1, sigma and betaL as NESTED functions, which in
 * MATLAB share the parent's workspace rather than getting their own. phi_u1
 * assigns `eta = (K-1)/K` internally, and the parent has already set
 * `eta = (N-2)/(N-1)` for the Theorem 3.5 prefactor 0.5/eta. At level 1 sigma
 * returns before it ever calls phi_u1, so eta survives and the bound is
 * correct; from level 2 on, sigma calls phi_u1(N-3), which OVERWRITES the
 * parent's eta, and the prefactor 0.5/eta is then evaluated with the wrong
 * value. On L = [1/2, 1/3, 1/5], N = 5 the parent eta is 3/4 and the clobbered
 * one is 1/2, so phi_u_n comes out 1.5x too large -- 2.4919 instead of
 * 1.6613 -- which exceeds the Section-2 baseline 1.7798, so `min` discards it.
 * The net effect is that MATLAB's SIB hierarchy is INERT above level 1: levels
 * 2, 3 and 4 all return the level-0 baseline, and the "bound" reported at
 * level 2 (X in [1.7407, 1.9335]) is LOOSER than the one at level 1
 * ([1.7607, 1.9335]), which contradicts the monotonicity the method is for.
 * A C++ port has no such aliasing -- the three helpers are lambdas with their
 * own scope -- so this port computes the intended Theorem 3.5 value and its
 * bounds do tighten with the level: level 2 gives 1.8182, level 3 gives
 * 1.8369, both still below the exact 1.8558. Reproducing the MATLAB value here
 * would mean deliberately reintroducing a scoping accident, so the divergence
 * is documented and pinned in the tests instead.
 *
 * DELAY IS REJECTED, following the reference: the no-delay phi bounds do not
 * bracket the with-delay congestion, so pfqn_sib.m raises
 * pfqn_sib:delayUnsupported for Z > 0 rather than return an invalid bracket.
 * The port throws InputError in the same case; returning a bound that is not
 * one would be worse than refusing.
 *
 * ARITHMETIC. Both Theorem 3.5 and Theorem 2.1 solve a quadratic, so the
 * bounds carry a square root and the routine is gated on
 * num_traits<T>::has_transcendental. That is a real restriction rather than a
 * formality: unlike the CBH and PBH families, SIB cannot be evaluated exactly.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_sib, mirroring [Xlo, Xhi, Wlo, Whi]. */
template <class T>
struct SibBounds {
    T Xlo;
    T Xhi;
    T Wlo;
    T Whi;
};

/**
 * @param L     (M) fixed-rate demands; delay demand is NOT accepted
 * @param N     population, at least 2
 * @param Z     think time, must be zero (see the header note)
 * @param level bound level >= 1, default 3 at the convenience overload
 */
template <class T>
SibBounds<T> pfqn_sib(const std::vector<T>& L, int N, const T& Z, int level) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_sib requires transcendental arithmetic (the bounds solve a quadratic)");
    using std::sqrt;
    const std::size_t M = L.size();
    if (M == 0) throw InputError("pfqn_sib: empty demand vector");
    if (N < 2) throw InputError("pfqn_sib: population must be at least two");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (Z > zero)
        throw InputError(
            "pfqn_sib supports Z=0 only (delay needs the Section-3.2 demand substitution)");
    const int lv = std::max(1, level);

    T Lsum = zero;
    for (const T& x : L) Lsum += x;
    if (Lsum == zero) throw InputError("pfqn_sib: all demands are zero");
    std::vector<T> rho(M);
    T rho_u = zero;
    for (std::size_t i = 0; i < M; ++i) {
        rho[i] = T(L[i] / Lsum);
        if (rho[i] > rho_u) rho_u = rho[i];
    }
    const int imax = lv + 3;
    std::vector<T> S(static_cast<std::size_t>(imax) + 1, zero);  // S[i], i = 1..imax
    for (int i = 1; i <= imax; ++i) {
        T s = zero;
        for (std::size_t m = 0; m < M; ++m) s += num_pow_int(rho[m], static_cast<unsigned>(i));
        S[static_cast<std::size_t>(i)] = s;
    }
    const T S2 = S[2];

    // alpha_i, i = 0..level (eqs. 3.5-3.6).
    std::vector<T> alpha(static_cast<std::size_t>(lv) + 1, zero);
    alpha[0] = S2;
    for (int i = 1; i <= lv; ++i) {
        T acc = zero;
        for (int j = 0; j <= i - 1; ++j)
            acc += S[static_cast<std::size_t>(i + 1 - j)] * alpha[static_cast<std::size_t>(j)];
        alpha[static_cast<std::size_t>(i)] = T(S[static_cast<std::size_t>(i + 2)] - acc);
    }

    // Level-1 upper bound on phi(K), eq. 3.18.
    const auto phi_u1 = [&](int K) -> T {
        if (K <= 0) return zero;
        if (K == 1) return S2;
        const T KT = num_traits<T>::from_int(K);
        const T eta = T(T(KT - one) / KT);
        const T T1 = T(T(KT - one) * rho_u - one);
        const T disc = T(T1 * T1 + num_traits<T>::from_int(4) * T(KT - one) * S2);
        return T(T(num_traits<T>::from_rational(1, 2) / eta) * T(T1 + sqrt(disc)));
    };

    // eq. (3.22c). NN plays the role of N-1.
    const auto sigma = [&](int NN, int i) -> T {
        T s = zero;
        if (i <= 0) return s;
        const T Dbar = T(one + phi_u1(NN - 2));
        T pnum = one;
        for (int j = 1; j <= i; ++j) {
            pnum *= num_traits<T>::from_int(NN - 1 - (j - 1));
            s += T(T(rho_u * S[static_cast<std::size_t>(j + 1)] - S[static_cast<std::size_t>(j + 2)]) *
                   pnum / num_pow_int(Dbar, static_cast<unsigned>(j)));
        }
        return s;
    };

    // eq. (3.23c). NN plays the role of N-1.
    const auto betaL = [&](int NN, int i) -> T {
        T b = zero;
        if (i <= 0) return b;
        for (int j = 1; j <= i - 1; ++j) {
            T p = one;
            for (int m = 2; m <= j; ++m)
                p *= T(num_traits<T>::from_int(NN - m) / T(one + phi_u1(NN - m)));
            b += alpha[static_cast<std::size_t>(j)] * p;
        }
        T p = one;
        for (int m = 2; m <= i; ++m) p *= T(num_traits<T>::from_int(NN - m) / T(one + phi_u1(NN - m)));
        const T Nim2 = num_traits<T>::from_int(NN - i - 2);  // MATLAB Nim2 = NN-1-i-1
        if (alpha[static_cast<std::size_t>(i) - 1] == zero)
            throw NumericError("pfqn_sib: zero alpha coefficient in the level correction");
        const T corr =
            T(one + T(alpha[static_cast<std::size_t>(i)] / alpha[static_cast<std::size_t>(i) - 1]) *
                        Nim2 / T(one + Nim2 * alpha[0]));
        b += alpha[static_cast<std::size_t>(i)] * p * corr;
        return b;
    };

    const int NN = N - 1;
    T phi_lo = T(num_traits<T>::from_int(N - 1) * S2);
    const T T1s2 = T(num_traits<T>::from_int(N - 1) * rho_u - one);
    T phi_hi = T(num_traits<T>::from_rational(1, 2) *
                 T(T1s2 + sqrt(T(T1s2 * T1s2 + num_traits<T>::from_int(4 * (N - 1)) * S2))));

    if (N >= 3) {
        const T eta = T(num_traits<T>::from_int(N - 2) / num_traits<T>::from_int(N - 1));
        const T T1u = T(num_traits<T>::from_int(N - 2) * rho_u - one);
        const T su = sigma(NN, lv - 1);
        T d = T(T1u * T1u + num_traits<T>::from_int(4 * (N - 2)) * T(S2 - su));
        if (d < zero) d = zero;
        const T phi_u_n = T(T(num_traits<T>::from_rational(1, 2) / eta) * T(T1u + sqrt(d)));
        if (phi_u_n < phi_hi) phi_hi = phi_u_n;

        const T T1l = T(num_traits<T>::from_int(N - 2) * S2 - one);
        const T bl = betaL(NN, lv - 1);
        T dl = T(T1l * T1l + num_traits<T>::from_int(4 * (N - 2)) *
                                 T(S2 + num_traits<T>::from_int(N - 2) * bl));
        if (dl < zero) dl = zero;
        // eq (3.23) divides by 2*eta, the SAME constant eq (3.22) applies as
        // 0.5/eta above; the Greek eta on the scan was read as the level index n
        const T phi_l_n = T(T(T1l + sqrt(dl)) / T(num_traits<T>::from_int(2) * eta));
        if (phi_l_n > phi_lo) phi_lo = phi_l_n;
    }

    if (phi_lo < zero) phi_lo = zero;
    if (phi_hi < phi_lo) phi_hi = phi_lo;

    SibBounds<T> r;
    r.Wlo = T(Lsum * T(one + phi_lo) + Z);
    r.Whi = T(Lsum * T(one + phi_hi) + Z);
    r.Xlo = T(num_traits<T>::from_int(N) / r.Whi);
    r.Xhi = T(num_traits<T>::from_int(N) / r.Wlo);
    return r;
}

template <class T>
SibBounds<T> pfqn_sib(const std::vector<T>& L, int N, const T& Z) {
    return pfqn_sib(L, N, Z, 3);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_SIB_H
