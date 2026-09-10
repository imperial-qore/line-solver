/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MVAOI_MARG_H
#define LINE_API_PFQN_MVAOI_MARG_H

/**
 * Exact marginal load-dependent MVA for a closed network of delay,
 * load-independent and ANY number of order-independent (OI) stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mvaoi_marg.m. This is the
 * marginal-distribution counterpart of pfqn_mvaoi (the mean-value CMVA form):
 * the same X and Q by a completely different route, which is what makes the
 * pair worth having.
 *
 * Each OI station carries its own joint COUNT-VECTOR marginal
 *
 *   pM_i(n | k) = (1/mu_i(n)) sum_r X_r(k) pM_i(n - e_r | k - e_r),   n != 0
 *   pM_i(0 | k) = 1 - sum_{n != 0} pM_i(n | k)
 *
 * driven by the common per-class throughput. The recursion is exact per station
 * because in product form pM_i(n|k) = Phi_i(n) G_{-i}(k-n)/G(k) with
 * X_r(k) = G(k-e_r)/G(k) and the balanced-fairness identity
 * Phi_i(n) = (1/mu_i(n)) sum_r Phi_i(n - e_r).
 *
 * Because the OI rate is class dependent, the mean-value response-time formula
 * is NOT exact, so X_r(k) is closed at every population level by population
 * conservation, X_r A_r + sum_i QM_ir(k;X) = k_r with A_r the non-OI residence
 * sum, and QM_ir read off the exact marginal. That closure is implicit in X and
 * the reference solves it by damped substitution (factor 1/2, tolerance 1e-13,
 * 2000 sweeps). This port keeps the same iteration; changing the damping or the
 * sweep count changes the last digits.
 *
 * RATE HANDLE CONVENTION, a genuine trap. pfqn_mvaoi's handles take the
 * per-class COUNT VECTOR n. This routine's take the MICROSTATE, the ordered
 * list of class indices with repetition (MATLAB's repelem(1:R, n)). The
 * reference keeps both conventions and converts between them in oi_rate; the
 * port keeps them too rather than silently unifying, so a handle written for
 * one routine is not accidentally accepted by the other. Microstate indices are
 * ZERO-based here, MATLAB's are one-based; a permutation-invariant rate, which
 * is what "order independent" means, cannot tell the difference, and any handle
 * that could is not an OI rate.
 *
 * Arithmetic: no transcendental, so it is left UNGATED and instantiates at
 * Rational. As with pfqn_momlin the exact instantiation is available rather
 * than advisable: the closure is a fixed point reached only in the limit, so
 * rational iterates grow without buying accuracy. Use double or Real.
 *
 * REFERENCE DEFECTS: none found.
 */

#include <algorithm>
#include <cstddef>
#include <functional>
#include <map>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_mvaoi_marg, mirroring [XN, QN]. */
template <class T>
struct MvaoiMargResult {
    std::vector<T> X;  ///< (R) per-class throughput
    Matrix<T> Q;       ///< (M x R) per-class queue length at every station
};

namespace detail {

/** All integer vectors 0 <= v <= bound, first component varying fastest. */
inline std::vector<std::vector<int>> enum_vecs(const std::vector<int>& bound) {
    std::vector<std::vector<int>> out;
    std::size_t total = 1;
    for (std::size_t r = 0; r < bound.size(); ++r)
        total *= static_cast<std::size_t>(bound[r]) + 1;
    out.reserve(total);
    std::vector<int> v(bound.size(), 0);
    for (std::size_t i = 0; i < total; ++i) {
        std::size_t li = i;
        for (std::size_t r = 0; r < bound.size(); ++r) {
            v[r] = static_cast<int>(li % (static_cast<std::size_t>(bound[r]) + 1));
            li /= static_cast<std::size_t>(bound[r]) + 1;
        }
        out.push_back(v);
    }
    return out;
}

/** Microstate of a count vector: class indices with repetition, ascending. */
inline std::vector<int> microstate(const std::vector<int>& n) {
    std::vector<int> mi;
    for (std::size_t r = 0; r < n.size(); ++r)
        for (int c = 0; c < n[r]; ++c) mi.push_back(static_cast<int>(r));
    return mi;
}

}  // namespace detail

/**
 * @param D       (M x R) per-class demands; the rows of OI stations are ignored
 * @param N       (R) closed populations
 * @param isDelay (M) true for infinite-server stations
 * @param mu      (M) rate handles; callable only at the OI stations, and taking
 *                the MICROSTATE (see the header note), not the count vector
 */
template <class T>
MvaoiMargResult<T> pfqn_mvaoi_marg(const Matrix<T>& D, const std::vector<int>& N,
                                   const std::vector<bool>& isDelay,
                                   const std::vector<std::function<T(const std::vector<int>&)>>& mu) {
    const std::size_t M = D.rows(), R = N.size();
    if (M == 0 || R == 0) throw InputError("pfqn_mvaoi_marg: empty model");
    if (D.cols() != R) throw InputError("pfqn_mvaoi_marg: D and N disagree on the class count");
    if (isDelay.size() != M || mu.size() != M)
        throw InputError("pfqn_mvaoi_marg: isDelay and mu must have one entry per station");
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] < 0) throw InputError("pfqn_mvaoi_marg: negative population");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T half = num_traits<T>::from_rational(1, 2);
    const T xtol = num_traits<T>::from_double(1e-13);

    std::vector<std::size_t> oi_list;
    std::vector<bool> isOI(M, false);
    for (std::size_t i = 0; i < M; ++i)
        if (mu[i]) {
            isOI[i] = true;
            oi_list.push_back(i);
        }
    const std::size_t nOI = oi_list.size();
    if (nOI == 0)
        throw InputError(
            "pfqn_mvaoi_marg: at least one order-independent station is required");

    // OI total rate at a count vector, through the microstate as the reference does
    std::vector<std::function<T(const std::vector<int>&)>> muM(nOI);
    for (std::size_t o = 0; o < nOI; ++o) {
        const std::function<T(const std::vector<int>&)> f = mu[oi_list[o]];
        muM[o] = [f, zero](const std::vector<int>& n) -> T {
            int tot = 0;
            for (int x : n) tot += x;
            if (tot == 0) return zero;
            return f(detail::microstate(n));
        };
    }

    std::map<std::vector<int>, std::vector<T>> X_cache;
    std::map<std::vector<int>, Matrix<T>> Q_cache;
    std::vector<std::map<std::vector<int>, std::map<std::vector<int>, T>>> pM(nOI);

    const std::vector<int> zeroK(R, 0);
    X_cache[zeroK] = std::vector<T>(R, zero);
    Q_cache[zeroK] = Matrix<T>(M, R, zero);
    for (std::size_t o = 0; o < nOI; ++o) pM[o][zeroK][zeroK] = one;

    std::vector<std::vector<int>> pops = detail::enum_vecs(N);
    std::stable_sort(pops.begin(), pops.end(),
                     [](const std::vector<int>& a, const std::vector<int>& b) {
                         int sa = 0, sb = 0;
                         for (int x : a) sa += x;
                         for (int x : b) sb += x;
                         if (sa != sb) return sa < sb;
                         return a < b;
                     });

    // marginal of one OI station at population k, given the throughput Xk
    const auto oi_marginal = [&](const std::vector<int>& k, const std::vector<T>& Xk,
                                 const std::function<T(const std::vector<int>&)>& mrate,
                                 const std::map<std::vector<int>, std::map<std::vector<int>, T>>& cache) {
        std::map<std::vector<int>, T> out;
        const std::vector<std::vector<int>> vecs = detail::enum_vecs(k);
        T psum = zero;
        bool haveZero = false;
        for (std::size_t idx = 0; idx < vecs.size(); ++idx) {
            const std::vector<int>& n = vecs[idx];
            int tot = 0;
            for (int x : n) tot += x;
            if (tot == 0) {
                haveZero = true;
                continue;
            }
            const T rate = mrate(n);
            if (!(rate > zero)) continue;
            T acc = zero;
            for (std::size_t r = 0; r < R; ++r) {
                if (n[r] < 1 || k[r] < 1) continue;
                std::vector<int> nr = n, kr = k;
                nr[r] -= 1;
                kr[r] -= 1;
                auto itk = cache.find(kr);
                if (itk == cache.end()) continue;
                auto itn = itk->second.find(nr);
                if (itn == itk->second.end()) continue;
                acc += Xk[r] * itn->second;
            }
            const T p = acc / rate;
            out[n] = p;
            psum += p;
        }
        if (haveZero) out[std::vector<int>(R, 0)] = one - psum;
        return out;
    };

    for (std::size_t pidx = 0; pidx < pops.size(); ++pidx) {
        const std::vector<int>& k = pops[pidx];
        int ktot = 0;
        for (int x : k) ktot += x;
        if (ktot == 0) continue;

        // non-OI response times by the arrival theorem at k - e_r
        Matrix<T> Rfix(M, R, zero);
        std::vector<T> A(R, zero);
        for (std::size_t r = 0; r < R; ++r) {
            if (k[r] == 0) continue;
            std::vector<int> kr = k;
            kr[r] -= 1;
            const Matrix<T>& Qkr = Q_cache.at(kr);
            for (std::size_t i = 0; i < M; ++i) {
                if (isOI[i]) continue;
                if (isDelay[i]) {
                    Rfix(i, r) = D(i, r);
                } else {
                    T s = zero;
                    for (std::size_t t = 0; t < R; ++t) s += Qkr(i, t);
                    Rfix(i, r) = D(i, r) * (one + s);
                }
                A[r] += Rfix(i, r);
            }
        }

        std::vector<T> Xk(R, zero);
        for (std::size_t r = 0; r < R; ++r)
            if (k[r] > 0) Xk[r] = num_traits<T>::from_int(k[r]) / (A[r] + one);

        std::vector<std::map<std::vector<int>, T>> marg(nOI);
        for (int it = 0; it < 2000; ++it) {
            std::vector<T> QMtot(R, zero);
            for (std::size_t o = 0; o < nOI; ++o) {
                marg[o] = oi_marginal(k, Xk, muM[o], pM[o]);
                for (auto itm = marg[o].begin(); itm != marg[o].end(); ++itm)
                    for (std::size_t r = 0; r < R; ++r)
                        QMtot[r] += num_traits<T>::from_int(itm->first[r]) * itm->second;
            }
            std::vector<T> Xnew(R, zero);
            for (std::size_t r = 0; r < R; ++r)
                if (k[r] > 0 && A[r] > zero) {
                    const T v = (num_traits<T>::from_int(k[r]) - QMtot[r]) / A[r];
                    Xnew[r] = (v > zero) ? v : zero;
                }
            T mx = zero;
            for (std::size_t r = 0; r < R; ++r) {
                const T d = num_abs(T(Xnew[r] - Xk[r]));
                if (d > mx) mx = d;
            }
            if (mx < xtol) {
                Xk = Xnew;
                break;
            }
            for (std::size_t r = 0; r < R; ++r) Xk[r] = half * Xk[r] + half * Xnew[r];
        }
        for (std::size_t o = 0; o < nOI; ++o) marg[o] = oi_marginal(k, Xk, muM[o], pM[o]);

        Matrix<T> Qk(M, R, zero);
        for (std::size_t o = 0; o < nOI; ++o)
            for (auto itm = marg[o].begin(); itm != marg[o].end(); ++itm)
                for (std::size_t r = 0; r < R; ++r)
                    Qk(oi_list[o], r) += num_traits<T>::from_int(itm->first[r]) * itm->second;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i)
                if (!isOI[i]) Qk(i, r) = Xk[r] * Rfix(i, r);

        X_cache[k] = Xk;
        Q_cache[k] = Qk;
        for (std::size_t o = 0; o < nOI; ++o) pM[o][k] = marg[o];
    }

    MvaoiMargResult<T> res;
    res.X = X_cache.at(N);
    res.Q = Q_cache.at(N);
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MVAOI_MARG_H
