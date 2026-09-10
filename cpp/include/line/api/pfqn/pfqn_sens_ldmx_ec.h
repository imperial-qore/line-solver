/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SENS_LDMX_EC_H
#define LINE_API_PFQN_SENS_LDMX_EC_H

/**
 * Effective capacity terms of the mixed load-dependent MVA of
 * Bruell-Balbo-Afshari, together with their exact derivatives with respect to
 * the open-class load Lo(i) of each station.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sens_ldmx_ec.m. The primal terms
 * duplicate pfqn_ldmx_ec.h; they are recomputed here rather than reused
 * because every intermediate (E1, F2, F3, F2prime) has to be carried alongside
 * its derivative, and the reference does the same.
 *
 * Lo(i) = sum_r lambda(r) D(i,r) is the only channel through which a demand
 * enters E, Eprime and EC, since the load-dependent rates mu do not depend on
 * the demands. Station i's terms depend on Lo(i) alone, so one derivative per
 * station suffices and the chain rule then gives the derivative with respect
 * to any demand-scaling parameter. This is the factorization behind equations
 * (19), (21) and (24)-(31) of the reference, reproduced here term by term.
 *
 * Arithmetic. Every term is a rational function of Lo and of the reciprocal
 * rates, differentiated symbolically, so the routine stays in the field and
 * instantiates at line::Rational with no gate.
 *
 * Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
 * Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
 * Communications 39(6):828-832, 1991.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
struct SensLdmxEcResult {
    Matrix<T> EC;       ///< (M x Nt) effective capacity, EC(i,n-1) for n = 1..Nt
    Matrix<T> E;        ///< (M x Nt+1) E-function, column n
    Matrix<T> Eprime;   ///< (M x Nt+1) Eprime-function, column n
    std::vector<T> Lo;  ///< (M) open load per station
    Matrix<T> dEC;      ///< (M x Nt) dEC(i,n)/dLo(i)
    Matrix<T> dE;       ///< (M x Nt+1) dE(i,n)/dLo(i)
    Matrix<T> dEprime;  ///< (M x Nt+1) dEprime(i,n)/dLo(i)
};

/**
 * @param lambda (R) arrival rates, zero on the closed classes
 * @param D      (M x R) service demands
 * @param mu     (M x Nt) load-dependent rate lattice, limited load dependence
 */
template <class T>
SensLdmxEcResult<T> pfqn_sens_ldmx_ec(const std::vector<T>& lambda, const Matrix<T>& D,
                                      const Matrix<T>& mu) {
    const std::size_t M = mu.rows();
    const std::size_t Nt = mu.cols();
    if (M == 0 || Nt == 0) throw InputError("pfqn_sens_ldmx_ec: the rate lattice is empty");
    if (D.rows() != M)
        throw InputError("pfqn_sens_ldmx_ec: D and mu disagree on the station count");
    if (lambda.size() != D.cols())
        throw InputError("pfqn_sens_ldmx_ec: lambda and D disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    SensLdmxEcResult<T> res;
    res.Lo.assign(M, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < lambda.size(); ++r) res.Lo[i] += lambda[r] * D(i, r);

    std::vector<std::size_t> b(M, 1);
    for (std::size_t i = 0; i < M; ++i) {
        const T& last = mu(i, Nt - 1);
        std::size_t k = 1;
        while (k < Nt && mu(i, k - 1) != last) ++k;
        b[i] = k;
    }

    // padding-width rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    const std::size_t Cw = 2 * Nt + 2;
    Matrix<T> C(M, Cw, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < Cw; ++k) {
            const T& rate = k < Nt ? mu(i, k) : mu(i, Nt - 1);
            if (rate == zero) throw NumericError("pfqn_sens_ldmx_ec: a load-dependent rate is zero");
            C(i, k) = one / rate;
        }
    // 1-based rate index of the reference.
    const auto Cq = [&](std::size_t i, std::size_t k) -> const T& { return C(i, k - 1); };

    res.EC = Matrix<T>(M, Nt, zero);
    res.E = Matrix<T>(M, Nt + 1, zero);
    res.Eprime = Matrix<T>(M, Nt + 1, zero);
    res.dEC = Matrix<T>(M, Nt, zero);
    res.dE = Matrix<T>(M, Nt + 1, zero);
    res.dEprime = Matrix<T>(M, Nt + 1, zero);

    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t bi = b[i];
        const T Cb = Cq(i, bi);
        const T den = one - res.Lo[i] * Cb;
        if (den == zero)
            throw NumericError(
                "pfqn_sens_ldmx_ec: the station is saturated by the open classes (Lo * C(b) = 1)");
        const std::size_t nhead = bi >= 2 ? bi - 1 : 0;  // n0 = 0 .. bi-2

        std::vector<T> E1(Nt + 1, zero), dE1(Nt + 1, zero);
        std::vector<T> E2(Nt + 1, zero), dE2(Nt + 1, zero);
        std::vector<T> E3(Nt + 1, zero), dE3(Nt + 1, zero);
        std::vector<T> E2p(Nt + 1, zero), dE2p(Nt + 1, zero);
        Matrix<T> F2(Nt + 1, nhead + 1, zero), dF2(Nt + 1, nhead + 1, zero);
        Matrix<T> F3(Nt + 1, nhead + 1, zero), dF3(Nt + 1, nhead + 1, zero);
        Matrix<T> F2p(Nt + 1, nhead + 1, zero), dF2p(Nt + 1, nhead + 1, zero);

        for (std::size_t n = 0; n <= Nt; ++n) {
            if (n >= bi) {
                // E(n) = 1/den^(n+1)  =>  dE/dLo = (n+1) Cb / den^(n+2)
                res.E(i, n) = one / num_pow_int(den, static_cast<unsigned>(n + 1));
                res.dE(i, n) = num_traits<T>::from_int(static_cast<long>(n + 1)) * Cb /
                               num_pow_int(den, static_cast<unsigned>(n + 2));
                res.Eprime(i, n) = Cb * res.E(i, n);
                res.dEprime(i, n) = Cb * res.dE(i, n);
                continue;
            }
            // ---- E1, eq. (25)-(26) -----------------------------------------
            if (n == 0) {
                E1[0] = one / den;
                dE1[0] = Cb / (den * den);
                for (std::size_t j = 1; j + 1 <= bi; ++j) {
                    E1[0] *= Cq(i, j) / Cb;
                    dE1[0] *= Cq(i, j) / Cb;
                }
            } else {
                const T fac = Cb / Cq(i, n);
                E1[n] = one / den * fac * E1[n - 1];
                dE1[n] = Cb / (den * den) * fac * E1[n - 1] + one / den * fac * dE1[n - 1];
            }

            // ---- F2, eq. (27)-(28) ------------------------------------------
            for (std::size_t n0 = 0; n0 <= nhead; ++n0) {
                if (n0 == 0) {
                    F2(n, 0) = one;
                    dF2(n, 0) = zero;
                } else {
                    const T coef = num_traits<T>::from_int(static_cast<long>(n + n0)) /
                                   num_traits<T>::from_int(static_cast<long>(n0)) * Cq(i, n + n0);
                    F2(n, n0) = coef * res.Lo[i] * F2(n, n0 - 1);
                    dF2(n, n0) = coef * (F2(n, n0 - 1) + res.Lo[i] * dF2(n, n0 - 1));
                }
            }
            E2[n] = zero;
            dE2[n] = zero;
            for (std::size_t n0 = 0; n0 < nhead; ++n0) {  // n0 = 0 .. bi-2
                E2[n] += F2(n, n0);
                dE2[n] += dF2(n, n0);
            }

            // ---- F3, eq. (29)-(30) ------------------------------------------
            for (std::size_t n0 = 0; n0 <= nhead; ++n0) {
                if (n == 0 && n0 == 0) {
                    F3(0, 0) = one;
                    for (std::size_t j = 1; j + 1 <= bi; ++j) F3(0, 0) *= Cq(i, j) / Cb;
                    dF3(0, 0) = zero;
                } else if (n > 0 && n0 == 0) {
                    const T fac = Cb / Cq(i, n);
                    F3(n, 0) = fac * F3(n - 1, 0);
                    dF3(n, 0) = fac * dF3(n - 1, 0);
                } else {
                    const T coef = num_traits<T>::from_int(static_cast<long>(n + n0)) /
                                   num_traits<T>::from_int(static_cast<long>(n0)) * Cb;
                    F3(n, n0) = coef * res.Lo[i] * F3(n, n0 - 1);
                    dF3(n, n0) = coef * (F3(n, n0 - 1) + res.Lo[i] * dF3(n, n0 - 1));
                }
            }
            E3[n] = zero;
            dE3[n] = zero;
            for (std::size_t n0 = 0; n0 < nhead; ++n0) {
                E3[n] += F3(n, n0);
                dE3[n] += dF3(n, n0);
            }

            // ---- F2prime -----------------------------------------------------
            for (std::size_t n0 = 0; n0 <= nhead; ++n0) {
                if (n0 == 0) {
                    F2p(n, 0) = Cq(i, n + 1);
                    dF2p(n, 0) = zero;
                } else {
                    const T coef = num_traits<T>::from_int(static_cast<long>(n + n0)) /
                                   num_traits<T>::from_int(static_cast<long>(n0)) *
                                   Cq(i, n + n0 + 1);
                    F2p(n, n0) = coef * res.Lo[i] * F2p(n, n0 - 1);
                    dF2p(n, n0) = coef * (F2p(n, n0 - 1) + res.Lo[i] * dF2p(n, n0 - 1));
                }
            }
            E2p[n] = zero;
            dE2p[n] = zero;
            for (std::size_t n0 = 0; n0 < nhead; ++n0) {
                E2p[n] += F2p(n, n0);
                dE2p[n] += dF2p(n, n0);
            }

            // ---- E = E1 + E2 - E3, eq. (23)-(24) ------------------------------
            res.E(i, n) = E1[n] + E2[n] - E3[n];
            res.dE(i, n) = dE1[n] + dE2[n] - dE3[n];
            if (n + 1 < bi) {
                res.Eprime(i, n) = Cb * E1[n] + E2p[n] - Cb * E3[n];
                res.dEprime(i, n) = Cb * dE1[n] + dE2p[n] - Cb * dE3[n];
            } else {
                res.Eprime(i, n) = Cb * res.E(i, n);
                res.dEprime(i, n) = Cb * res.dE(i, n);
            }
        }

        // EC(n) = C(n) E(n)/E(n-1), eq. (19), and the quotient rule for dEC.
        for (std::size_t n = 1; n <= Nt; ++n) {
            if (res.E(i, n - 1) == zero) throw NumericError("pfqn_sens_ldmx_ec: E vanishes");
            res.EC(i, n - 1) = Cq(i, n) * res.E(i, n) / res.E(i, n - 1);
            res.dEC(i, n - 1) = Cq(i, n) *
                                (res.dE(i, n) * res.E(i, n - 1) - res.E(i, n) * res.dE(i, n - 1)) /
                                (res.E(i, n - 1) * res.E(i, n - 1));
        }
    }

    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SENS_LDMX_EC_H
