/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_DAC_H
#define LINE_API_PFQN_PFQN_DAC_H

/**
 * Distribution Analysis by Chain (de Souza e Silva, UCLA CSD-870023, 1987):
 * the JOINT queue-length distribution of a closed product-form network with
 * single-server, infinite-server and queue-dependent centers.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_dac.m together with its local
 * functions dac_lattice, dac_compositions, dac_rank and dac_step. The
 * recursion runs over a related network in which every chain holds one
 * customer, a transformation that leaves the aggregate queue-length
 * distribution unchanged. Adding one customer with demands r to a network of
 * k customers gives
 *
 *   c_j      = sum_{n=1..k+1} (n/mu_j(n)) P_j^{k}(n-1)
 *   lambda   = 1 / sum_j r_j c_j
 *   P^{k+1}(n) = lambda sum_j r_j (n_j/mu_j(n_j)) P^{k}(n - e_j)
 *
 * so probability mass is conserved by construction and the recursion is
 * numerically stable, unlike a normalizing-constant route. Per-chain
 * throughputs and queue lengths come from re-running the tail of the recursion
 * with each chain placed last, sharing the common prefix.
 *
 * ARITHMETIC. Every step is an addition, a multiplication or a division of
 * field elements, so the whole joint distribution is EXACT in rational
 * arithmetic and the routine is deliberately left ungated. That is the point
 * of the algorithm here: an exact joint distribution is what availability
 * modelling needs, and it is the natural oracle for the marginal that MVA
 * returns.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_dac, mirroring [Pjoint, states, XN, QN, UN, CN, pi]. */
template <class T>
struct DacResult {
    std::vector<T> Pjoint;                 ///< one probability per row of `states`
    std::vector<std::vector<int>> states;  ///< aggregate states, J columns
    std::vector<T> XN;                     ///< (R) chain throughput
    Matrix<T> QN;                          ///< (M x R) mean queue length
    std::vector<T> UN;                     ///< (M) utilization, 1 - P(empty)
    std::vector<T> CN;                     ///< (R) cycle time, excluding think time
    Matrix<T> PI;                          ///< (M x (Nt+1)) marginal distribution
};

namespace detail {

/** All J-part compositions of k, in lexicographic order (dac_compositions). */
inline void dac_compositions(std::size_t J, int k, std::vector<std::vector<int>>& out) {
    out.clear();
    std::vector<int> n(J, 0);
    if (J == 0) return;
    n[J - 1] = k;
    // Lexicographic order on the leading parts: enumerate by decrementing the
    // rightmost position that can give a unit to a position on its left.
    while (true) {
        out.push_back(n);
        std::size_t i = J - 1;
        while (i > 0 && n[i] == 0) --i;
        if (i == 0) break;
        // Move one unit from position i to position i-1 and flush the rest
        // rightwards, which walks the compositions in lexicographic order.
        n[i - 1] += 1;
        const int rest = n[i] - 1;
        for (std::size_t d = i; d < J; ++d) n[d] = 0;
        n[J - 1] = rest;
    }
}

/** Lexicographic rank (1-based, as in dac_rank) of a composition. */
inline std::size_t dac_rank(const std::vector<int>& n, int k, std::size_t J,
                            const std::vector<std::vector<double>>& C) {
    std::size_t idx = 1;
    int rem = k;
    // 1-based vs 0-based binomial index rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    for (std::size_t j = 0; j + 1 < J; ++j) {
        const std::size_t col = J - j - 2;
        for (int v = 0; v < n[j]; ++v) {
            const int s = rem - v;
            idx += static_cast<std::size_t>(
                C[static_cast<std::size_t>(s + static_cast<int>(J - j) - 2)][col]);
        }
        rem -= n[j];
    }
    return idx;
}

}  // namespace detail

/**
 * @param L  (M x R) demands
 * @param N  (R) population
 * @param Z  (R) think times; a non-zero total appends an IS center, so the
 *           states then have M+1 columns
 * @param mu (M x Nt) load-dependent rates, empty for all ones
 */
template <class T>
DacResult<T> pfqn_dac(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                      const Matrix<T>& mu) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_dac: L and N disagree on the class count");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    int Nt = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_dac: the population vector must be non-negative");
        Nt += v;
    }
    const std::size_t mucols = static_cast<std::size_t>(Nt > 0 ? Nt : 1);
    if (!mu.empty()) {
        if (mu.rows() != M) throw InputError("pfqn_dac: mu must have one row per station");
        if (Nt > 0 && mu.cols() < static_cast<std::size_t>(Nt))
            throw InputError("pfqn_dac: mu must have at least sum(N) columns");
    }

    std::vector<T> Zv = Z;
    if (Zv.empty()) Zv.assign(R, zero);
    T Zsum = zero;
    for (const T& v : Zv) Zsum += v;
    const bool hasZ = Zsum > zero;
    const std::size_t J = M + (hasZ ? 1 : 0);

    Matrix<T> Lx(J, R), mux(J, mucols);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < R; ++r) Lx(i, r) = L(i, r);
        for (std::size_t k = 0; k < mucols; ++k)
            mux(i, k) = mu.empty() ? one : mu(i, k);
    }
    if (hasZ) {
        for (std::size_t r = 0; r < R; ++r) Lx(M, r) = Zv[r];
        for (std::size_t k = 0; k < mucols; ++k)
            mux(M, k) = num_traits<T>::from_int(static_cast<long>(k) + 1);
    }

    DacResult<T> res;
    if (Nt == 0) {
        res.states.assign(1, std::vector<int>(J, 0));
        res.Pjoint.assign(1, one);
        res.XN.assign(R, zero);
        res.QN = Matrix<T>(M, R, zero);
        res.UN.assign(M, zero);
        res.CN.assign(R, zero);
        res.PI = Matrix<T>(M, 1, one);
        return res;
    }

    std::vector<std::size_t> active;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > 0) active.push_back(r);
    for (std::size_t idx = 0; idx < active.size(); ++idx) {
        bool any = false;
        for (std::size_t j = 0; j < J; ++j)
            if (Lx(j, active[idx]) > zero) any = true;
        if (!any) throw InputError("pfqn_dac: a chain has null demand at every center");
    }

    // Lattice of aggregate states, level by level, with the successor map.
    std::vector<std::vector<std::vector<int>>> lv(static_cast<std::size_t>(Nt) + 1);
    for (int k = 0; k <= Nt; ++k) detail::dac_compositions(J, k, lv[static_cast<std::size_t>(k)]);
    std::vector<std::vector<double>> C(
        static_cast<std::size_t>(Nt) + J + 1, std::vector<double>(J + 2, 0.0));
    for (std::size_t a = 0; a + 1 <= static_cast<std::size_t>(Nt) + J; ++a)
        for (std::size_t b = 0; b <= (a < J ? a : J); ++b)
            C[a][b] = (b == 0) ? 1.0 : (C[a - 1][b - 1] + C[a - 1][b]);
    std::vector<std::vector<std::vector<std::size_t>>> succ(static_cast<std::size_t>(Nt) + 1);
    for (int k = 0; k < Nt; ++k) {
        const std::vector<std::vector<int>>& Sv = lv[static_cast<std::size_t>(k)];
        succ[static_cast<std::size_t>(k)].assign(Sv.size(), std::vector<std::size_t>(J, 0));
        for (std::size_t i = 0; i < Sv.size(); ++i)
            for (std::size_t j = 0; j < J; ++j) {
                std::vector<int> t = Sv[i];
                t[j] += 1;
                succ[static_cast<std::size_t>(k)][i][j] = detail::dac_rank(t, k + 1, J, C) - 1;
            }
    }

    // One recursion step: add one customer with demands r to a k-customer net.
    const auto dac_step = [&](const std::vector<T>& p, const std::vector<T>& r, int k,
                              std::vector<T>& pn, T& lam, std::vector<T>& Lq) {
        const std::vector<std::vector<int>>& Sv = lv[static_cast<std::size_t>(k)];
        Matrix<T> marg(J, static_cast<std::size_t>(k) + 1, zero);
        for (std::size_t i = 0; i < Sv.size(); ++i)
            for (std::size_t j = 0; j < J; ++j)
                marg(j, static_cast<std::size_t>(Sv[i][j])) += p[i];
        std::vector<T> c(J, zero);
        for (std::size_t j = 0; j < J; ++j)
            for (int n = 1; n <= k + 1; ++n) {
                const T mur = mux(j, static_cast<std::size_t>(n) - 1);
                if (mur == zero) throw NumericError("pfqn_dac: zero load-dependent rate");
                c[j] += T(num_traits<T>::from_int(n) / mur) *
                        marg(j, static_cast<std::size_t>(n) - 1);
            }
        T den = zero;
        for (std::size_t j = 0; j < J; ++j) den += r[j] * c[j];
        if (den == zero) throw NumericError("pfqn_dac: zero cycle time for the added customer");
        lam = T(one / den);
        Lq.assign(J, zero);
        for (std::size_t j = 0; j < J; ++j) Lq[j] = T(lam * r[j] * c[j]);

        pn.assign(lv[static_cast<std::size_t>(k) + 1].size(), zero);
        for (std::size_t j = 0; j < J; ++j) {
            if (!(r[j] > zero)) continue;
            for (std::size_t i = 0; i < Sv.size(); ++i) {
                const int nj = Sv[i][j] + 1;
                const T mur = mux(j, static_cast<std::size_t>(nj) - 1);
                const T w = T(lam * r[j] * T(num_traits<T>::from_int(nj) / mur) * p[i]);
                pn[succ[static_cast<std::size_t>(k)][i][j]] += w;
            }
        }
    };

    // Chain order: the D distinct chains last, sharing the common prefix.
    const std::size_t D = active.size();
    std::vector<std::size_t> prefix;
    for (std::size_t idx = 0; idx < D; ++idx)
        for (int t = 0; t + 1 < N[active[idx]]; ++t) prefix.push_back(active[idx]);

    std::vector<T> p(1, one);
    int k = 0;
    std::vector<T> pn, Lq;
    T lam = zero;
    const auto column = [&](std::size_t r) {
        std::vector<T> col(J);
        for (std::size_t j = 0; j < J; ++j) col[j] = Lx(j, r);
        return col;
    };
    for (std::size_t idx = 0; idx < prefix.size(); ++idx) {
        dac_step(p, column(prefix[idx]), k, pn, lam, Lq);
        p = pn;
        ++k;
    }

    std::vector<std::vector<T>> Sp(D);
    Sp[0] = p;
    std::vector<T> pb = p;
    int kb = k;
    for (std::size_t idx = 0; idx < D; ++idx) {
        dac_step(pb, column(active[idx]), kb, pn, lam, Lq);
        pb = pn;
        ++kb;
        if (idx + 1 < D) Sp[idx + 1] = pb;
    }
    res.Pjoint = pb;
    res.states = lv[static_cast<std::size_t>(Nt)];

    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    {
        const std::size_t r = active[D - 1];
        res.XN[r] = T(num_traits<T>::from_int(N[r]) * lam);
        for (std::size_t i = 0; i < M; ++i)
            res.QN(i, r) = T(num_traits<T>::from_int(N[r]) * Lq[i]);
    }
    for (std::size_t idx = 0; idx + 1 < D; ++idx) {
        std::vector<T> pc = Sp[idx];
        int kc = k + static_cast<int>(idx);
        std::vector<std::size_t> order;
        for (std::size_t t = idx + 1; t < D; ++t) order.push_back(active[t]);
        order.push_back(active[idx]);
        for (std::size_t t = 0; t < order.size(); ++t) {
            dac_step(pc, column(order[t]), kc, pn, lam, Lq);
            pc = pn;
            ++kc;
        }
        const std::size_t r = active[idx];
        res.XN[r] = T(num_traits<T>::from_int(N[r]) * lam);
        for (std::size_t i = 0; i < M; ++i)
            res.QN(i, r) = T(num_traits<T>::from_int(N[r]) * Lq[i]);
    }

    res.PI = Matrix<T>(M, static_cast<std::size_t>(Nt) + 1, zero);
    for (std::size_t i = 0; i < res.states.size(); ++i)
        for (std::size_t j = 0; j < M; ++j)
            res.PI(j, static_cast<std::size_t>(res.states[i][j])) += res.Pjoint[i];
    res.UN.assign(M, zero);
    for (std::size_t j = 0; j < M; ++j) res.UN[j] = T(one - res.PI(j, 0));
    res.CN.assign(R, zero);
    for (std::size_t idx = 0; idx < D; ++idx) {
        const std::size_t r = active[idx];
        if (res.XN[r] == zero) throw NumericError("pfqn_dac: zero chain throughput");
        res.CN[r] = T(num_traits<T>::from_int(N[r]) / res.XN[r] - Zv[r]);
    }
    return res;
}

template <class T>
DacResult<T> pfqn_dac(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z) {
    return pfqn_dac(L, N, Z, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_DAC_H
