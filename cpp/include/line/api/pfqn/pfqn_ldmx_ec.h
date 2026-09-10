/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LDMX_EC_H
#define LINE_API_PFQN_LDMX_EC_H

/**
 * Bruell-Balbo-Afshari effective-capacity terms for a MIXED open/closed
 * network with limited load dependence.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ldmx_ec.m.
 *
 * Station i is limited load dependent: its rate lattice mu(i,.) is arbitrary
 * up to the saturation level b_i, the first k with mu(i,k) = mu(i,end), and
 * constant beyond it. Writing C(i,k) = 1/mu(i,k) and Lo_i = sum_r lambda_r
 * D(i,r) for the open load, the E-function
 *
 *   E_i(n) = sum_{n0 >= 0} C(n+n0, n0) Lo_i^{n0} prod_{j=n+1}^{n+n0} C(i,j)
 *
 * has the closed form E_i(n) = 1/(1 - Lo_i C(i,b_i))^{n+1} once n >= b_i, and
 * below saturation is assembled from three finite pieces E1 + E2 - E3, with
 * E1 the geometric tail, E2 the exact head, and E3 the head of the geometric
 * that E1 double counts. Eprime is the same construction one index up, and the
 * effective capacity is the ratio
 *
 *   EC_i(n) = C(i,n) E_i(n) / E_i(n-1),   n = 1, ..., Nt.
 *
 * Its reciprocal is the load-dependent rate that turns the mixed model into a
 * purely closed one, which is what pfqn_ncldmx and pfqn_mvaldmx consume.
 *
 * Arithmetic: EXACT-CAPABLE. Every term is a sum, product, quotient or integer
 * power in the field of the inputs; the reference itself uses no
 * transcendental function. The one thing the caller must respect is the
 * stability condition Lo_i C(i,b_i) < 1: at Lo_i C(i,b_i) = 1 the geometric
 * denominator vanishes and the routine throws, and above 1 the sum that E
 * represents diverges even though the closed form still evaluates. That check
 * is exact in rational arithmetic and is not a tolerance.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
struct LdmxEcResult {
    Matrix<T> EC;       ///< (M x Nt) effective capacity, EC(i,n) for n = 1..Nt
    Matrix<T> E;        ///< (M x Nt+1) E-function, column 1+n
    Matrix<T> Eprime;   ///< (M x Nt+1) Eprime-function, column 1+n
    std::vector<T> Lo;  ///< (M) open load per station
};

/**
 * @param lambda (R) arrival rates, zero on the closed classes
 * @param D      (M x R) service demands
 * @param mu     (M x Nt) load-dependent rate lattice
 */
template <class T>
LdmxEcResult<T> pfqn_ldmx_ec(const std::vector<T>& lambda, const Matrix<T>& D,
                             const Matrix<T>& mu) {
    const std::size_t M = mu.rows();
    const std::size_t Nt = mu.cols();
    if (M == 0 || Nt == 0) throw InputError("pfqn_ldmx_ec: the rate lattice is empty");
    if (D.rows() != M) throw InputError("pfqn_ldmx_ec: D and mu disagree on the station count");
    if (lambda.size() != D.cols())
        throw InputError("pfqn_ldmx_ec: lambda and D disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    LdmxEcResult<T> res;
    res.Lo.assign(M, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < lambda.size(); ++r) res.Lo[i] += lambda[r] * D(i, r);

    // b_i: the first index at which the rate reaches its saturation value.
    std::vector<std::size_t> b(M, 1);
    for (std::size_t i = 0; i < M; ++i) {
        const T& last = mu(i, Nt - 1);
        std::size_t k = 1;
        while (k < Nt && mu(i, k - 1) != last) ++k;
        b[i] = k;
    }
    std::size_t bmax = 1;
    for (std::size_t i = 0; i < M; ++i) bmax = b[i] > bmax ? b[i] : bmax;

    // C = 1/mu, extended past the lattice by the saturation rate, exactly as
    // the reference pads mu with max(b)+1 copies of its last column.
    const std::size_t Cw = Nt + bmax + 2;
    Matrix<T> C(M, Cw);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t k = 0; k < Cw; ++k) {
            const T& rate = k < Nt ? mu(i, k) : mu(i, Nt - 1);
            if (rate == zero) throw NumericError("pfqn_ldmx_ec: a load-dependent rate is zero");
            C(i, k) = one / rate;
        }
    }
    // Cq(i,k) is C(i,k) at the reference's 1-based rate index k >= 1.
    const auto Cq = [&](std::size_t i, std::size_t k) -> const T& { return C(i, k - 1); };

    res.EC = Matrix<T>(M, Nt, zero);
    res.E = Matrix<T>(M, Nt + 1, zero);
    res.Eprime = Matrix<T>(M, Nt + 1, zero);

    for (std::size_t i = 0; i < M; ++i) {
        const T Cb = Cq(i, b[i]);
        const T denom = one - res.Lo[i] * Cb;
        if (denom == zero)
            throw NumericError(
                "pfqn_ldmx_ec: the station is saturated by the open classes (Lo * C(b) = 1), the "
                "effective capacity is unbounded");
        const T geo = one / denom;
        // Number of head terms: n0 = 0 .. b_i - 2, empty when b_i = 1.
        const std::size_t nhead = b[i] >= 2 ? b[i] - 1 : 0;

        std::vector<T> E1(Nt + 1, zero);
        for (std::size_t n = 0; n <= Nt; ++n) {
            if (n >= b[i]) {
                res.E(i, n) = num_pow_int(geo, static_cast<unsigned>(n + 1));
                res.Eprime(i, n) = Cb * res.E(i, n);
                continue;
            }
            // ---- E1: the geometric tail --------------------------------------
            if (n == 0) {
                E1[0] = geo;
                for (std::size_t j = 1; j + 1 <= b[i]; ++j) E1[0] *= Cq(i, j) / Cb;
            } else {
                E1[n] = geo * Cb / Cq(i, n) * E1[n - 1];
            }
            // ---- E2 and E2prime: the exact head ------------------------------
            T E2 = zero, E2p = zero;
            {
                T F2 = one;             // n0 = 0
                T F2p = Cq(i, n + 1);   // n0 = 0
                for (std::size_t n0 = 0; n0 < nhead; ++n0) {
                    if (n0 > 0) {
                        const T w = num_traits<T>::from_int(static_cast<long>(n + n0)) /
                                    num_traits<T>::from_int(static_cast<long>(n0));
                        F2 = w * res.Lo[i] * Cq(i, n + n0) * F2;
                        F2p = w * res.Lo[i] * Cq(i, n + n0 + 1) * F2p;
                    }
                    E2 += F2;
                    E2p += F2p;
                }
            }
            // ---- E3: the head of the geometric that E1 double counts ----------
            T E3 = zero;
            {
                // F3(n, 0) obeys its own recursion in n, so it is rebuilt here
                // from the n = 0 seed rather than carried across iterations.
                T F3 = one;
                for (std::size_t j = 1; j + 1 <= b[i]; ++j) F3 *= Cq(i, j) / Cb;
                for (std::size_t k = 1; k <= n; ++k) F3 = Cb / Cq(i, k) * F3;
                for (std::size_t n0 = 0; n0 < nhead; ++n0) {
                    if (n0 > 0) {
                        const T w = num_traits<T>::from_int(static_cast<long>(n + n0)) /
                                    num_traits<T>::from_int(static_cast<long>(n0));
                        F3 = w * res.Lo[i] * Cb * F3;
                    }
                    E3 += F3;
                }
            }
            res.E(i, n) = E1[n] + E2 - E3;
            if (n + 1 < b[i])
                res.Eprime(i, n) = Cb * E1[n] + E2p - Cb * E3;
            else
                res.Eprime(i, n) = Cb * res.E(i, n);
        }

        for (std::size_t n = 1; n <= Nt; ++n) {
            if (res.E(i, n - 1) == zero) throw NumericError("pfqn_ldmx_ec: E vanishes");
            res.EC(i, n - 1) = Cq(i, n) * res.E(i, n) / res.E(i, n - 1);
        }
    }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LDMX_EC_H
