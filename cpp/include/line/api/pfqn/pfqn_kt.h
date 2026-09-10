/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_KT_H
#define LINE_API_PFQN_PFQN_KT_H

/**
 * Knessl-Tier asymptotic expansion of the normalizing constant.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_kt.m, including the two fixes the
 * MATLAB header documents (the think-time term, absent from stock pfqn_kt in
 * both the exponent and the Hessian, and the evaluation of the exponent at the
 * exact saddle point rather than at the AQL throughput). In LINE's convention
 * Cauchy extraction plus steepest descent on
 *
 *   F(u) = sum_r Z_r u_r - sum_k log(1 - U_k) - sum_r N_r log u_r,  U_k = L(k,:) u
 *
 * gives log G = F(u*) - sum_r log u*_r - (R/2) log(2 pi) - (1/2) log det H with
 *
 *   H_rs = delta_rs N_r/u_r^2 + sum_k L_kr L_ks/(1 - U_k)^2
 *
 * and u* the solution of N_r = u_r (Z_r + sum_k L_kr/(1 - U_k)), found by
 * damped Newton from the AMVA throughput.
 *
 * SELF-LOOPING CLASSES. A class that visits exactly one station and has no
 * think time would drive U_k to 1. Extracting u_r^N_r from 1/(1-U_ist) is exact
 * and leaves L(ist,r)^N_r with that station's factor raised to (1-V_ist)^-(1+N_r),
 * so the class is dropped, the station replicated N_r times and N_r log L(ist,r)
 * added to the exponent. Classes looping at the SAME station share one factor
 * of exponent 1+sum N_r and add the multinomial (sum N_r)!/prod N_r!.
 *
 * ARITHMETIC. Logarithms, an asymptotic expansion and a Newton iteration, so
 * gated on num_traits<T>::has_transcendental.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_aql.h"
#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_kt, mirroring [G, lG, X, Q]. */
template <class T>
struct KtResult {
    T G;
    T lG;
    std::vector<T> X;  ///< AMVA throughput used to seed the saddle point
    Matrix<T> Q;       ///< AMVA queue lengths
};

/**
 * @param L0 (M x R) demands
 * @param N0 (R) population
 * @param Z0 (R) think times
 */
template <class T>
KtResult<T> pfqn_kt(const Matrix<T>& L0, const std::vector<T>& N0, const std::vector<T>& Z0) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_kt requires transcendental arithmetic (steepest-descent expansion of log G)");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T fineTol = num_traits<T>::from_double(1e-8);  // GlobalConstants.FineTol

    KtResult<T> res;
    T Nt0 = zero;
    for (const T& v : N0) Nt0 += v;
    if (L0.empty() || N0.empty() || Nt0 == zero) {
        res.G = one;
        res.lG = zero;
        return res;
    }
    const std::size_t Rorig = L0.cols(), Morig = L0.rows();
    std::vector<T> Zin = Z0;
    if (Zin.empty()) Zin.assign(Rorig, zero);
    if (Zin.size() != Rorig || N0.size() != Rorig)
        throw InputError("pfqn_kt: L, N and Z disagree on the class count");

    // A class with no jobs contributes a factor of 1 to G, but its saddle point is
    // u_r -> 0, where N_r*log(u_r) and N_r/u_r^2 are indeterminate and lG comes back
    // NaN. Solve the reduced model, as the self-looping fold below already does.
    std::size_t nkeep = 0;
    for (std::size_t r = 0; r < Rorig; ++r)
        if (N0[r] > zero) ++nkeep;
    if (nkeep > 0 && nkeep < Rorig) {
        Matrix<T> Lk(Morig, nkeep, zero);
        std::vector<T> Nk, Zk;
        Nk.reserve(nkeep);
        Zk.reserve(nkeep);
        std::size_t c = 0;
        for (std::size_t r = 0; r < Rorig; ++r) {
            if (!(N0[r] > zero)) continue;
            for (std::size_t i = 0; i < Morig; ++i) Lk(i, c) = L0(i, r);
            Nk.push_back(N0[r]);
            Zk.push_back(Zin[r]);
            ++c;
        }
        return pfqn_kt(Lk, Nk, Zk);
    }

    // Fold self-looping classes into replicated stations.
    T slcdemandfactor = zero;
    std::vector<bool> isslc(Rorig, false);
    std::vector<std::vector<T>> rows;  // extra replicated station rows
    for (std::size_t i = 0; i < Morig; ++i) {
        std::vector<T> row(Rorig);
        for (std::size_t r = 0; r < Rorig; ++r) row[r] = L0(i, r);
        rows.push_back(row);
    }
    std::vector<std::size_t> slcstation(Rorig, 0);
    if (Rorig > 1) {
        for (std::size_t r = 0; r < Rorig; ++r) {
            std::size_t nnz = 0, ist = 0;
            for (std::size_t i = 0; i < Morig; ++i)
                if (L0(i, r) != zero) {
                    ++nnz;
                    ist = i;
                }
            if (nnz != 1 || Zin[r] != zero) continue;
            isslc[r] = true;
            slcstation[r] = ist;
        }
        // Classes looping at the SAME station share one (1-V)^-(1+sum N) factor
        // and contribute the multinomial (sum N)!/prod N_r!.
        std::vector<bool> done(Rorig, false);
        for (std::size_t r = 0; r < Rorig; ++r) {
            if (!isslc[r] || done[r]) continue;
            const std::size_t ist = slcstation[r];
            T ntot = zero;
            for (std::size_t s = r; s < Rorig; ++s) {
                if (!isslc[s] || slcstation[s] != ist) continue;
                done[s] = true;
                ntot += N0[s];
                slcdemandfactor += T(N0[s] * log(L0(ist, s))) - detail::num_factln<T>(N0[s]);
            }
            slcdemandfactor += detail::num_factln<T>(ntot);
            const long reps = static_cast<long>(num_traits<T>::to_double(ntot));
            for (long k = 0; k < reps; ++k) {
                std::vector<T> row(Rorig);
                for (std::size_t s = 0; s < Rorig; ++s) row[s] = L0(ist, s);
                rows.push_back(row);
            }
        }
    }
    std::vector<std::size_t> keep;
    for (std::size_t r = 0; r < Rorig; ++r)
        if (!isslc[r]) keep.push_back(r);
    if (keep.empty()) {  // every class self-loops: the demand factors are exact
        res.lG = slcdemandfactor;
        res.G = exp(slcdemandfactor);
        return res;
    }
    const std::size_t M = rows.size(), R = keep.size();
    Matrix<T> L(M, R);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) L(i, r) = rows[i][keep[r]];
    std::vector<T> N(R), Z(R);
    for (std::size_t r = 0; r < R; ++r) {
        N[r] = N0[keep[r]];
        Z[r] = Zin[keep[r]];
    }
    T Ntot = zero;
    for (const T& v : N) Ntot += v;
    if (Ntot == zero) {  // only self-looping classes carried jobs
        res.lG = slcdemandfactor;
        res.G = exp(slcdemandfactor);
        return res;
    }

    // AMVA seed.
    const AmvaResult<T> amva =
        (num_traits<T>::to_double(Ntot) <= 4.0) ? pfqn_bs(L, N, Z) : pfqn_aql(L, N, Z);
    res.X = amva.XN;
    res.Q = amva.QN;

    std::vector<T> u = amva.XN;
    // Keep the saddle point inside the domain U_k < 1.
    T Umax = zero;
    for (std::size_t k = 0; k < M; ++k) {
        T s = zero;
        for (std::size_t r = 0; r < R; ++r) s += L(k, r) * u[r];
        if (s > Umax) Umax = s;
    }
    if (Umax >= one) {
        const T f = T(T(one - num_traits<T>::from_double(1e-6)) / Umax);
        for (std::size_t r = 0; r < R; ++r) u[r] = T(u[r] * f);
    }

    bool converged = false;
    for (int it = 0; it < 200; ++it) {
        std::vector<T> D(M);
        for (std::size_t k = 0; k < M; ++k) {
            T s = zero;
            for (std::size_t r = 0; r < R; ++r) s += L(k, r) * u[r];
            D[k] = T(one / T(one - s));
        }
        std::vector<T> g(R);
        T gn = zero;
        for (std::size_t r = 0; r < R; ++r) {
            T s = zero;
            for (std::size_t k = 0; k < M; ++k) s += L(k, r) * D[k];
            g[r] = T(u[r] * T(Z[r] + s) - N[r]);
            gn += g[r] * g[r];
        }
        using std::sqrt;
        if (sqrt(gn) <= T(num_traits<T>::from_double(1e-12) * Ntot)) {
            converged = true;
            break;
        }
        Matrix<T> J(R, R, zero);
        for (std::size_t r = 0; r < R; ++r) {
            T s = zero;
            for (std::size_t k = 0; k < M; ++k) s += L(k, r) * D[k];
            J(r, r) = T(Z[r] + s);
            for (std::size_t sIdx = 0; sIdx < R; ++sIdx) {
                T acc = zero;
                for (std::size_t k = 0; k < M; ++k) acc += L(k, r) * T(D[k] * D[k]) * L(k, sIdx);
                J(r, sIdx) += u[r] * acc;
            }
        }
        std::vector<T> rhs(R);
        for (std::size_t r = 0; r < R; ++r) rhs[r] = T(-g[r]);
        const std::vector<T> du = solve(J, rhs);
        T alpha = one;
        for (int b = 0; b < 60; ++b) {
            bool ok = true;
            for (std::size_t r = 0; r < R && ok; ++r)
                if (T(u[r] + alpha * du[r]) <= zero) ok = false;
            if (ok) {
                for (std::size_t k = 0; k < M && ok; ++k) {
                    T s = zero;
                    for (std::size_t r = 0; r < R; ++r) s += L(k, r) * T(u[r] + alpha * du[r]);
                    if (s >= one) ok = false;
                }
            }
            if (ok) break;
            alpha = T(alpha / num_traits<T>::from_int(2));
        }
        if (num_traits<T>::to_double(alpha) < 1e-12) break;
        for (std::size_t r = 0; r < R; ++r) u[r] = T(u[r] + alpha * du[r]);
    }

    std::vector<T> us = converged ? u : amva.XN;
    if (converged) {
        // Second acceptance test of the reference, on the residual at u.
        std::vector<T> D(M);
        for (std::size_t k = 0; k < M; ++k) {
            T s = zero;
            for (std::size_t r = 0; r < R; ++r) s += L(k, r) * u[r];
            D[k] = T(one / T(one - s));
        }
        T gn = zero;
        for (std::size_t r = 0; r < R; ++r) {
            T s = zero;
            for (std::size_t k = 0; k < M; ++k) s += L(k, r) * D[k];
            const T gr = T(u[r] * T(Z[r] + s) - N[r]);
            gn += gr * gr;
        }
        using std::sqrt;
        if (!(sqrt(gn) <= T(num_traits<T>::from_double(1e-8) * Ntot))) us = amva.XN;
    }

    // Assemble the expansion at us.
    std::vector<T> Uk(M), D(M);
    for (std::size_t k = 0; k < M; ++k) {
        T s = zero;
        for (std::size_t r = 0; r < R; ++r) s += L(k, r) * us[r];
        Uk[k] = s;
        T den = T(one - s);
        if (den < fineTol) den = fineTol;
        D[k] = T(one / den);
    }
    Matrix<T> H(R, R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        if (us[r] == zero) throw NumericError("pfqn_kt: zero saddle-point coordinate");
        H(r, r) = T(N[r] / T(us[r] * us[r]));
        for (std::size_t s = 0; s < R; ++s) {
            T acc = zero;
            for (std::size_t k = 0; k < M; ++k) acc += L(k, r) * T(D[k] * D[k]) * L(k, s);
            H(r, s) += acc;
        }
    }
    T F = zero;
    for (std::size_t r = 0; r < R; ++r) F += Z[r] * us[r];
    for (std::size_t k = 0; k < M; ++k) {
        T den = T(one - Uk[k]);
        if (den < fineTol) den = fineTol;
        F -= log(den);
    }
    for (std::size_t r = 0; r < R; ++r) F -= N[r] * log(us[r]);

    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    T lG = F;
    for (std::size_t r = 0; r < R; ++r) lG -= log(us[r]);
    lG -= T(num_traits<T>::from_rational(static_cast<long>(R), 2) * log(twopi));
    lG -= T(num_traits<T>::from_rational(1, 2) * detail::pfqn_logdet(H));
    lG += slcdemandfactor;
    res.lG = lG;
    res.G = exp(lG);
    return res;
}

template <class T>
KtResult<T> pfqn_kt(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_kt(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_KT_H
