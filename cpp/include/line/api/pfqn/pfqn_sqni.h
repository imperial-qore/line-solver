/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_SQNI_H
#define LINE_API_PFQN_PFQN_SQNI_H

/**
 * Square-root non-iterative (SQNI) approximation for a single queueing station
 * with per-class delay.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sqni.m. Each class throughput is
 * the admissible root of a quadratic assembled from the balanced-job estimate
 * B_r of the other classes' contribution:
 *
 *   X_r = (Z_r - sqrt(disc) - B_r + L_r Ntot) / (2 L_r Z_r)
 *   disc = B_r^2 - 2 B_r L_r Ntot - 2 B_r Z_r + L_r^2 Ntot^2 + 2 L_r Ntot Z_r
 *          - 4 N_r L_r Z_r + Z_r^2
 *
 * with the discriminant clamped at zero, as in the reference. Classes with
 * Z_r = 0 (self-looping classes, whose quadratic is degenerate) are handled by
 * the reference's two-pass structure: their queue length is pre-set to N_r,
 * and their throughput is filled in from the total queue length afterwards.
 *
 * ARITHMETIC. The square root is essential, so the routine is gated on
 * num_traits<T>::has_transcendental.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_sqni, mirroring [Q, U, X] for the single station. */
template <class T>
struct SqniResult {
    std::vector<T> Q;
    std::vector<T> U;
    std::vector<T> X;
};

/**
 * @param N (R) population, @param L (R) demand at the station,
 * @param Z (R) think times
 */
template <class T>
SqniResult<T> pfqn_sqni(const std::vector<T>& N, const std::vector<T>& L, const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_sqni requires transcendental arithmetic (the throughput solves a quadratic)");
    using std::sqrt;
    const std::size_t C = L.size();
    if (N.size() != C || Z.size() != C)
        throw InputError("pfqn_sqni: N, L and Z must have the same length");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2), four = num_traits<T>::from_int(4);

    SqniResult<T> r;
    r.Q.assign(C, zero);
    r.U.assign(C, zero);
    r.X.assign(C, zero);
    T Nt = zero;
    for (const T& v : N) Nt += v;
    if (Nt <= zero) return r;

    if (Nt == one) {
        for (std::size_t s = 0; s < C; ++s) {
            const T den = T(Z[s] + L[s]);
            if (den == zero) throw NumericError("pfqn_sqni: zero cycle time");
            r.X[s] = T(N[s] / den);
            r.U[s] = T(r.X[s] * L[s]);
            r.Q[s] = r.U[s];
        }
        return r;
    }

    for (std::size_t s = 0; s < C; ++s)
        if (Z[s] == zero) r.Q[s] = N[s];

    for (std::size_t s = 0; s < C; ++s) {
        if (Z[s] == zero) continue;  // handled after the main loop
        const T Nr = N[s], Lr = L[s], Zr = Z[s];
        // B_r: balanced-job estimate of the delay-resident population of the
        // other classes, at the population with one class-r job removed.
        T inner = zero;
        for (std::size_t t = 0; t < C; ++t) {
            T N1 = N[t];
            if (t == s) N1 -= one;
            const T d = T(Z[t] + L[t] + L[t] * T(Nt - two));
            if (d == zero) continue;
            inner += T(Z[t] * N1 / d);
        }
        T Brsum = zero;
        for (std::size_t t = 0; t < C; ++t) {
            if (t == s) continue;
            const T d = T(Z[t] + L[t] + L[t] * T(Nt - one - inner));
            if (d == zero) continue;
            Brsum += T(N[t] / d * Z[t]);
        }
        const T Br = T(Lr * Brsum);
        if (Lr == zero) {
            r.X[s] = T(Nr / Zr);
        } else {
            T disc = T(Br * Br - two * Br * Lr * Nt - two * Br * Zr + Lr * Lr * Nt * Nt +
                       two * Lr * Nt * Zr - four * Nr * Lr * Zr + Zr * Zr);
            if (disc < zero) disc = zero;
            r.X[s] = T(T(Zr - sqrt(disc) - Br + Lr * Nt) / T(two * Lr * Zr));
        }
        r.U[s] = T(r.X[s] * L[s]);
        r.Q[s] = T(N[s] - r.X[s] * Z[s]);
    }

    T Qtot = zero;
    for (const T& v : r.Q) Qtot += v;
    for (std::size_t s = 0; s < C; ++s) {
        if (Z[s] != zero) continue;
        if (L[s] == zero) throw NumericError("pfqn_sqni: a class has neither demand nor think time");
        r.X[s] = T(N[s] / T(L[s] * T(one + Qtot)));
        r.U[s] = T(r.X[s] * L[s]);
        r.Q[s] = T(N[s] - r.X[s] * Z[s]);
    }
    return r;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_SQNI_H
