/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_AGHQ_H
#define LINE_API_PFQN_PFQN_AGHQ_H

/**
 * Adaptive Gauss-Hermite quadrature of the simplex factor of the McKenna-Mitra integral.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_aghq.m. Rescaling by the logistic-expansion
 * mode and curvature, w = w* + A^(-1/2) z, and applying the q-node probabilists'
 * Gauss-Hermite rule in each of the M-1 simplex directions gives a convergent rule whose
 * q = 1 member is pfqn_le itself: one node at the mode with weight sqrt(2 pi). So LE is
 * the first term of a convergent quadrature rather than an approximation of unknown
 * accuracy. Cost q^(M-1), which is what confines the method to small M.
 *
 * A tensor rule is NOT invariant to the choice of square root of A: any B with B B' =
 * inv(A) is admissible and they place the nodes differently. The principal-axis frame is
 * used here, as in the reference results; where two curvatures are close to equal that
 * frame is close to arbitrary and two valid rules can part company well above their own
 * error, converging back together as q grows. Do NOT compare across codebases node by
 * node.
 *
 * With Z > 0 the radius is integrated numerically (pfqn_simplex.h) and the rule is applied
 * to the M-1 simplex directions, so every node costs one radial quadrature. Because the
 * radius is integrated rather than Laplaced, q = 1 there is the logistic expansion with an
 * exact radius, which is NOT pfqn_le's own Z > 0 branch.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental exactly as pfqn_le is.
 *
 * References:
 * J. McKenna, D. Mitra, "Integral Representations and Asymptotic Expansions for Closed
 * Markovian Queueing Networks: Normal Usage", Bell Syst. Tech. J. 61(5), 1982.
 * Q. Liu, D. A. Pierce, "A Note on Gauss-Hermite Quadrature", Biometrika 81(3), 1994.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_le.h"
#include "line/api/pfqn/pfqn_simplex.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
using AghqResult = LeResult<T>;

template <class T>
AghqResult<T> pfqn_aghq(const Matrix<T>& L, const std::vector<T>& N,
                        const std::vector<T>& Z, std::size_t q) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_aghq requires transcendental arithmetic (quadrature of an integral)");
    using std::exp;
    using std::log;
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    if (q == 0) throw InputError("pfqn_aghq: q must be at least 1");

    T Ntot = zero, Lsum = zero, Zsum = zero;
    for (const T& x : N) Ntot += x;
    for (const T& x : Z) Zsum += x;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Lsum += L(i, r);
    if (M == 0 || N.empty() || Ntot == zero || num_traits<T>::to_double(Lsum) < 1e-4) {
        return pfqn_le(L, N, Z);
    }

    const std::size_t d = M - 1;
    const T half = num_traits<T>::from_double(0.5);
    AghqResult<T> res;

    if (Z.empty() || num_traits<T>::to_double(Zsum) < lang::GlobalConstants::Zero) {
        const std::vector<T> umax = pfqn_le_fpi(L, N);
        // M = 1 leaves a 0x0 reduced Hessian, which pfqn_le_hessian refuses to form;
        // d = 0 then makes the rule a single node and ld = 0.
        const Matrix<T> A = (M < 2) ? Matrix<T>(0, 0) : pfqn_le_hessian(L, N, umax);
        const T ld = (M < 2) ? zero : detail::pfqn_logdet(A);
        T S = zero, sum_lu = zero;
        for (std::size_t r = 0; r < R; ++r) {
            T uL = zero;
            for (std::size_t i = 0; i < M; ++i) uL += umax[i] * L(i, r);
            S += N[r] * log(uL);
        }
        for (std::size_t i = 0; i < M; ++i) sum_lu += log(umax[i]);
        const T h0 = T(S + sum_lu);
        std::vector<T> w0(d);
        for (std::size_t i = 0; i < d; ++i) w0[i] = log(T(umax[i] / umax[M - 1]));
        const Matrix<T>* Lp = &L;
        const std::vector<T>* Np = &N;
        auto h = [Lp, Np, M, R, zero](const std::vector<T>& w) {
            using std::log;
            std::vector<T> x = simplex::softmax_gauge(w);
            T acc = zero;
            for (std::size_t r = 0; r < R; ++r) {
                T c = zero;
                for (std::size_t i = 0; i < M; ++i) c += x[i] * (*Lp)(i, r);
                acc += (*Np)[r] * log(c);
            }
            for (std::size_t i = 0; i < M; ++i) acc += log(x[i]);
            return acc;
        };
        T lacc = simplex::aghq_rule(h, w0, h0, A, q, d);
        T lG = detail::num_factln<T>(T(Ntot + num_traits<T>::from_int(static_cast<long>(M) - 1)));
        for (std::size_t r = 0; r < R; ++r) lG -= detail::num_factln<T>(N[r]);
        res.lG = T(lG + h0 + lacc - half * ld);
    } else {
        simplex::Mode<T> mode = simplex::simplex_mode(L, N, Z);
        std::vector<T> w0(d);
        for (std::size_t i = 0; i < d; ++i) w0[i] = log(T(mode.x[i] / mode.x[M - 1]));
        const Matrix<T>* Lp = &L;
        const std::vector<T>* Np = &N;
        const std::vector<T>* Zp = &Z;
        auto h = [Lp, Np, Zp, M, R, zero](const std::vector<T>& w) {
            using std::log;
            std::vector<T> x = simplex::softmax_gauge(w);
            std::vector<T> c(R, zero);
            for (std::size_t r = 0; r < R; ++r)
                for (std::size_t i = 0; i < M; ++i) c[r] += x[i] * (*Lp)(i, r);
            T acc = simplex::radial(c, *Np, *Zp, M).lJ;
            for (std::size_t i = 0; i < M; ++i) acc += log(x[i]);
            return acc;
        };
        T lacc = simplex::aghq_rule(h, w0, mode.h0, mode.A, q, d);
        T lG = zero;
        for (std::size_t r = 0; r < R; ++r) lG -= detail::num_factln<T>(N[r]);
        res.lG = T(lG + mode.h0 + lacc - half * mode.ld);
    }
    res.G = exp(res.lG);
    return res;
}

template <class T>
AghqResult<T> pfqn_aghq(const Matrix<T>& L, const std::vector<T>& N,
                        const std::vector<T>& Z) {
    return pfqn_aghq(L, N, Z, static_cast<std::size_t>(3));
}

template <class T>
AghqResult<T> pfqn_aghq(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_aghq(L, N, std::vector<T>(), static_cast<std::size_t>(3));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_AGHQ_H
