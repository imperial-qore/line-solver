/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MVAC_H
#define LINE_API_PFQN_MVAC_H

/**
 * MVAC: exact mean value analysis BY CHAIN of a closed multichain product-form
 * network (Conway, de Souza e Silva and Lavenberg, IEEE Trans. Computers
 * 38(3):432-442, 1989).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mvac.m.
 *
 * Where the classic MVA recursion of pfqn_mva recurs on the POPULATION vector
 * and costs O(prod(N+1)), MVAC recurs on the CHAINS: each chain is reduced to
 * single-customer chains and the removed ones are replaced by self-looping
 * single-customer (SCSL) chains pinned at a service center. The multiplicity
 * vector v = (v_1,...,v_J), v_j the number of SCSL chains at center j, indexes
 * the recursion in place of the population. Writing L^k_j(v) for the mean
 * number at center j with the SCSL customers excluded,
 *
 *   lambda^k_k(v) = 1 / ( a_k + sum_{j SSFR} a_jk (L^{k-1}_j(v) + v_j) )     (10)
 *   L^k_{jk}(v)   = lambda^k_k(v) a_jk (1 + L^{k-1}_j(v) + v_j),  j SSFR     (9a)
 *   L^k_{jk}(v)   = lambda^k_k(v) a_jk,                           j IS       (9b)
 *   L^k_i(v)      = sum_j L^k_{jk}(v) L^{k-1}_i(v + 1_j) + L^k_{ik}(v)       (7)
 *   L^k_{il}(v)   = sum_j L^k_{jk}(v) L^{k-1}_{il}(v + 1_j),  l < k          (6)
 *
 * with L^0 = 0, read off at k = K, v = 0. Part 1 evaluates (10), (9) and (7);
 * part 2 evaluates (6) for the chains that visit at least one IS center, whose
 * throughput then follows from Little's law there; the chains that visit only
 * SSFR centers need a re-execution of part 1 with their label interchanged with
 * K, which is cheap because the levels below the interchanged label are
 * untouched and are reused.
 *
 * IDENTICAL CHAINS. Classes with N_r > 1, and classes with identical demand
 * columns, collapse into one subset of identical single-customer chains: only
 * the representative is analyzed and its per-chain measures are multiplied by
 * the class population. The cost therefore depends on the number of DISTINCT
 * chains, not on K. The subsets are found by MATLAB's
 * unique(...,'rows','stable'), whose first-appearance order the port
 * reproduces, because the chain labelling (representatives last, IS-visiting
 * ones before the rest) is built from it and part 3 interchanges labels by
 * position.
 *
 * NO NORMALIZING CONSTANT is formed anywhere, so MVAC does not suffer the
 * underflow and overflow that complicate RECAL and convolution.
 *
 * Arithmetic: EXACT-CAPABLE. Additions, multiplications and divisions in the
 * field of the inputs only, no logarithm and no tolerance. Instantiated at
 * Rational it returns the same throughputs and queue lengths as pfqn_mva and
 * pfqn_ca as exact fractions.
 *
 * REFERENCE DEFECTS: none found. Agreement with pfqn_mva is to the last ulp on
 * every model tried, including models with a delay, without a delay, with
 * repeated demand columns and with empty classes.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_comb_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_mvac, mirroring [XN, QN, UN, CN]. */
template <class T>
struct MvacResult {
    std::vector<T> X;  ///< (R) per-class throughput
    Matrix<T> Q;       ///< (M x R) per-class queue length at the SSFR queues
    Matrix<T> U;       ///< (M x R) per-class utilization
    Matrix<T> C;       ///< (M x R) per-class residence time
};

namespace detail {

/**
 * Chain subsets, labelling and multiplicity-vector lattice shared by pfqn_mvac
 * and pfqn_mvacld. Both build exactly the same objects from the same A matrix;
 * only the recursion that runs over them differs.
 */
template <class T>
struct MvacSetup {
    std::size_t J1 = 0, J = 0, D = 0, S = 0;
    std::size_t K = 0;
    Matrix<T> a;                            ///< (J x K) per-chain demands
    std::vector<std::size_t> posr;          ///< classes with N > 0
    std::vector<std::size_t> grpOfClass;    ///< subset of each populated class (0-based)
    std::vector<std::size_t> gorder;        ///< subset order: IS-visiting first
    std::vector<std::vector<int>> Vlist;    ///< multiplicity vectors, by increasing sum
    std::vector<int> vsum;
    std::vector<std::size_t> off, cnt;      ///< block offset and size per component sum
    Matrix<long> succ;                      ///< index of v + 1_j, or -1
};

template <class T>
MvacSetup<T> mvac_setup(const Matrix<T>& A, const std::vector<int>& N, std::size_t J1,
                        std::size_t J, std::size_t K) {
    const T zero = num_traits<T>::from_int(0);
    MvacSetup<T> st;
    st.J1 = J1;
    st.J = J;
    st.K = K;

    for (std::size_t r = 0; r < N.size(); ++r)
        if (N[r] > 0) st.posr.push_back(r);

    // distinct demand columns, first-appearance order (MATLAB unique 'stable')
    std::vector<std::vector<T>> Adist;
    st.grpOfClass.assign(st.posr.size(), 0);
    for (std::size_t p = 0; p < st.posr.size(); ++p) {
        std::vector<T> colv(J, zero);
        for (std::size_t j = 0; j < J; ++j) colv[j] = A(j, st.posr[p]);
        std::size_t g = Adist.size();
        for (std::size_t q = 0; q < Adist.size(); ++q)
            if (Adist[q] == colv) {
                g = q;
                break;
            }
        if (g == Adist.size()) Adist.push_back(colv);
        st.grpOfClass[p] = g;
    }
    const std::size_t Dall = Adist.size();

    std::vector<bool> visitsIS(Dall, false);
    for (std::size_t g = 0; g < Dall; ++g)
        for (std::size_t j = J1; j < J; ++j)
            if (Adist[g][j] > zero) {
                visitsIS[g] = true;
                break;
            }
    for (std::size_t g = 0; g < Dall; ++g)
        if (visitsIS[g]) st.gorder.push_back(g);
    for (std::size_t g = 0; g < Dall; ++g)
        if (!visitsIS[g]) st.gorder.push_back(g);
    st.D = Dall;
    st.S = 0;
    for (std::size_t g = 0; g < Dall; ++g)
        if (!visitsIS[g]) ++st.S;

    // chain labels: the D representatives take K-D+1..K, the rest fill 1..K-D
    std::vector<std::size_t> mult(st.D, 0);
    for (std::size_t g = 0; g < st.D; ++g)
        for (std::size_t p = 0; p < st.posr.size(); ++p)
            if (st.grpOfClass[p] == st.gorder[g])
                mult[g] += static_cast<std::size_t>(N[st.posr[p]]);
    std::vector<std::size_t> chainGroup(K, 0);
    for (std::size_t g = 0; g < st.D; ++g) chainGroup[K - st.D + g] = g;
    std::size_t p = 0;
    for (std::size_t g = 0; g < st.D; ++g)
        for (std::size_t cc = 0; cc + 1 < mult[g]; ++cc) chainGroup[p++] = g;

    st.a = Matrix<T>(J, K, zero);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t j = 0; j < J; ++j) st.a(j, k) = Adist[st.gorder[chainGroup[k]]][j];

    // multiplicity vectors, enumerated by increasing component sum
    st.off.assign(K + 1, 0);
    st.cnt.assign(K + 1, 0);
    for (std::size_t t = 0; t <= K; ++t) {
        const std::vector<std::vector<int>> Vt =
            multichoose_rows(static_cast<int>(J), static_cast<int>(t));
        st.off[t] = st.Vlist.size();
        st.cnt[t] = Vt.size();
        for (std::size_t i = 0; i < Vt.size(); ++i) {
            st.Vlist.push_back(Vt[i]);
            st.vsum.push_back(static_cast<int>(t));
        }
    }
    const std::size_t nv = st.Vlist.size();
    st.succ = Matrix<long>(nv, J, -1);
    for (std::size_t vi = 0; vi < nv; ++vi) {
        const int t = st.vsum[vi];
        if (static_cast<std::size_t>(t) + 1 > K) continue;
        for (std::size_t j = 0; j < J; ++j) {
            std::vector<int> w = st.Vlist[vi];
            w[j] += 1;
            for (std::size_t q = 0; q < st.cnt[t + 1]; ++q)
                if (st.Vlist[st.off[t + 1] + q] == w) {
                    st.succ(vi, j) = static_cast<long>(st.off[t + 1] + q);
                    break;
                }
        }
    }
    return st;
}

/** Part 1 of the basic step: (10), (9a)-(9b) and (7) for k = k0..K over sum(v) <= K-k. */
template <class T>
void mvac_part1(std::size_t k0, const MvacSetup<T>& st, std::vector<Matrix<T>>& Lall,
                std::vector<Matrix<T>>& Ljkall, std::vector<std::vector<T>>& lamall) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t J = st.J, J1 = st.J1, K = st.K, nv = st.Vlist.size();
    for (std::size_t k = k0; k <= K; ++k) {
        const Matrix<T>& Lp = Lall[k - 1];
        Matrix<T> Lk(J, nv, zero), Ljk(J, nv, zero);
        std::vector<T> lamv(nv, zero);
        for (std::size_t vi = 0; vi < nv; ++vi) {
            if (static_cast<std::size_t>(st.vsum[vi]) > K - k) continue;
            T den = zero;
            for (std::size_t j = 0; j < J; ++j) den += st.a(j, k - 1);  // a_k over ALL centers
            for (std::size_t j = 0; j < J1; ++j)
                den += (Lp(j, vi) + num_traits<T>::from_int(st.Vlist[vi][j])) * st.a(j, k - 1);
            if (den <= zero) throw NumericError("pfqn_mvac: a chain has zero total demand");
            const T lam = one / den;
            for (std::size_t j = 0; j < J1; ++j)
                Ljk(j, vi) = lam * (one + Lp(j, vi) + num_traits<T>::from_int(st.Vlist[vi][j])) *
                             st.a(j, k - 1);
            for (std::size_t j = J1; j < J; ++j) Ljk(j, vi) = lam * st.a(j, k - 1);
            lamv[vi] = lam;
            for (std::size_t i = 0; i < J; ++i) {
                T s = Ljk(i, vi);
                for (std::size_t j = 0; j < J; ++j) {
                    const long sj = st.succ(vi, j);
                    if (sj < 0) throw NumericError("pfqn_mvac: multiplicity vector out of range");
                    s += Ljk(j, vi) * Lp(i, static_cast<std::size_t>(sj));
                }
                Lk(i, vi) = s;
            }
        }
        Lall[k] = Lk;
        Ljkall[k] = Ljk;
        lamall[k] = lamv;
    }
}

}  // namespace detail

/**
 * @param L (M x R) demands of the single-server fixed-rate queues
 * @param N (R) closed populations
 * @param Z (Mz x R) demands of the infinite-server centers, one row per center
 */
template <class T>
MvacResult<T> pfqn_mvac(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    const std::size_t M = L.rows(), R = L.cols();
    if (M == 0 || R == 0) throw InputError("pfqn_mvac: empty demand matrix");
    if (N.size() != R) throw InputError("pfqn_mvac: L and N disagree on the class count");
    if (!Z.empty() && Z.cols() != R)
        throw InputError("pfqn_mvac: the think time matrix and the demand matrix disagree");
    const T zero = num_traits<T>::from_int(0);

    MvacResult<T> res;
    res.X.assign(R, zero);
    res.Q = Matrix<T>(M, R, zero);
    res.U = Matrix<T>(M, R, zero);
    res.C = Matrix<T>(M, R, zero);

    std::size_t K = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < 0) throw InputError("pfqn_mvac: the population vector must be nonnegative");
        K += static_cast<std::size_t>(N[r]);
    }
    if (K == 0) return res;

    // discard the centers no chain visits
    std::vector<std::size_t> ssfrIdx, isIdx;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (L(i, r) > zero) {
                ssfrIdx.push_back(i);
                break;
            }
    for (std::size_t i = 0; i < Z.rows(); ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (Z(i, r) > zero) {
                isIdx.push_back(i);
                break;
            }
    const std::size_t J1 = ssfrIdx.size(), J = J1 + isIdx.size();
    if (J == 0) throw InputError("pfqn_mvac: all service demands are zero, the throughput is unbounded");

    Matrix<T> A(J, R, zero);
    for (std::size_t j = 0; j < J1; ++j)
        for (std::size_t r = 0; r < R; ++r) A(j, r) = L(ssfrIdx[j], r);
    for (std::size_t j = J1; j < J; ++j)
        for (std::size_t r = 0; r < R; ++r) A(j, r) = Z(isIdx[j - J1], r);

    detail::MvacSetup<T> st = detail::mvac_setup(A, N, J1, J, K);
    const std::size_t nv = st.Vlist.size();
    const std::size_t z0 = 0;  // the zero multiplicity vector is enumerated first

    std::vector<Matrix<T>> Lall(K + 1), Ljkall(K + 1);
    std::vector<std::vector<T>> lamall(K + 1);
    Lall[0] = Matrix<T>(J, nv, zero);

    std::vector<T> lamChain(K, zero);
    Matrix<T> Lchain(J, K, zero);

    detail::mvac_part1(1, st, Lall, Ljkall, lamall);
    lamChain[K - 1] = lamall[K][z0];
    for (std::size_t j = 0; j < J; ++j) Lchain(j, K - 1) = Ljkall[K](j, z0);

    // ---- part 2: chains that visit at least one IS center ---------------------------
    const std::size_t D = st.D, S = st.S;
    const long lmaxK = std::min(static_cast<long>(K) - 1, static_cast<long>(K - S));
    if (D >= 2 && lmaxK >= static_cast<long>(K - D + 1)) {
        std::vector<Matrix<T>> L2prev(K + 1), L2cur(K + 1);
        for (std::size_t k = K - D + 2; k <= K; ++k) {
            L2cur.assign(K + 1, Matrix<T>());
            const long lhi = std::min(static_cast<long>(k) - 1, static_cast<long>(K - S));
            for (long l = static_cast<long>(K - D + 1); l <= lhi; ++l) {
                Matrix<T> acc(J, nv, zero);
                for (std::size_t vi = 0; vi < nv; ++vi) {
                    if (static_cast<std::size_t>(st.vsum[vi]) != K - k) continue;
                    for (std::size_t j = 0; j < J; ++j) {
                        const long sj = st.succ(vi, j);
                        if (sj < 0) throw NumericError("pfqn_mvac: multiplicity vector out of range");
                        const Matrix<T>& prev =
                            (l == static_cast<long>(k) - 1) ? Ljkall[k - 1] : L2prev[l];
                        for (std::size_t i = 0; i < J; ++i)
                            acc(i, vi) += Ljkall[k](j, vi) * prev(i, static_cast<std::size_t>(sj));
                    }
                }
                L2cur[l] = acc;
            }
            if (k == K) {
                for (long l = static_cast<long>(K - D + 1); l <= lmaxK; ++l) {
                    for (std::size_t j = 0; j < J; ++j)
                        Lchain(j, static_cast<std::size_t>(l) - 1) = L2cur[l](j, z0);
                    std::size_t jIS = J;
                    for (std::size_t j = J1; j < J; ++j)
                        if (st.a(j, static_cast<std::size_t>(l) - 1) > zero) {
                            jIS = j;
                            break;
                        }
                    if (jIS == J) throw NumericError("pfqn_mvac: chain has no IS center");
                    lamChain[static_cast<std::size_t>(l) - 1] =
                        Lchain(jIS, static_cast<std::size_t>(l) - 1) /
                        st.a(jIS, static_cast<std::size_t>(l) - 1);
                }
            }
            L2prev = L2cur;
        }
    }

    // ---- part 3: chains that visit only SSFR centers, by label interchange ----------
    std::vector<std::size_t> perm(K);
    for (std::size_t k = 0; k < K; ++k) perm[k] = k;
    for (std::size_t l = 1; l + 1 <= S && S >= 1 && l < S; ++l) {
        std::swap(perm[K - l - 1], perm[K - 1]);
        for (std::size_t j = 0; j < J; ++j) {
            const T tmp = st.a(j, K - l - 1);
            st.a(j, K - l - 1) = st.a(j, K - 1);
            st.a(j, K - 1) = tmp;
        }
        detail::mvac_part1(K - l, st, Lall, Ljkall, lamall);
        lamChain[perm[K - 1]] = lamall[K][z0];
        for (std::size_t j = 0; j < J; ++j) Lchain(j, perm[K - 1]) = Ljkall[K](j, z0);
    }

    // ---- expand the per-chain measures back to per-class ----------------------------
    for (std::size_t g = 0; g < D; ++g) {
        const std::size_t kg = K - D + g;
        for (std::size_t p = 0; p < st.posr.size(); ++p) {
            if (st.grpOfClass[p] != st.gorder[g]) continue;
            const std::size_t r = st.posr[p];
            const T nr = num_traits<T>::from_int(N[r]);
            res.X[r] = nr * lamChain[kg];
            for (std::size_t j = 0; j < J1; ++j) {
                res.Q(ssfrIdx[j], r) = nr * Lchain(j, kg);
                res.U(ssfrIdx[j], r) = res.X[r] * L(ssfrIdx[j], r);
                res.C(ssfrIdx[j], r) = res.Q(ssfrIdx[j], r) / res.X[r];
            }
        }
    }
    // an empty class reports its bare demand as residence time, as pfqn_mva does
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] == 0)
            for (std::size_t i = 0; i < M; ++i) res.C(i, r) = L(i, r);
    return res;
}

/** Overload with the zero think-time default. */
template <class T>
MvacResult<T> pfqn_mvac(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_mvac(L, N, Matrix<T>(1, L.cols(), num_traits<T>::from_int(0)));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MVAC_H
