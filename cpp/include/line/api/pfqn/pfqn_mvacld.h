/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MVACLD_H
#define LINE_API_PFQN_MVACLD_H

/**
 * MVAC for networks with queue-length dependent (QLD) service centers, the
 * Section V extension of Conway, de Souza e Silva and Lavenberg (1989).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mvacld.m. pfqn_mvac implements
 * Sections II-IV, which cover single-server fixed-rate and infinite-server
 * centers only.
 *
 * Where pfqn_mvac propagates the MEAN queue lengths through (7), the QLD
 * extension propagates the MARGINAL DISTRIBUTIONS P^k_j(n,v). That is forced
 * by load dependence, since the rate seen by a job depends on the whole
 * occupancy, but it also SIMPLIFIES the recursion: (21)-(25) read level k-1
 * only at the shifted vectors v + 1_i, so the basic step sweeps I_k alone where
 * pfqn_mvac must sweep I_k u ... u I_K. The marginals come out as a first-class
 * output for free. In the reference-station-free form used here,
 *
 *   c_i(k,v)      = sum_{n=0}^{k-1} P^{k-1}_i(n, v+1_i) mu_i(n+v_i+1)/(n+v_i+1)
 *   L^k_{jk}(v)   = (a_jk/c_j) / sum_m (a_mk/c_m)
 *   lambda^k_k(v) = 1 / sum_m (a_mk/c_m)
 *   P^k_j(n,v)    = L^k_{jk}(v) P^{k-1}_j(n-1, v+1_j)
 *                   + sum_{m != j} L^k_{mk}(v) P^{k-1}_j(n, v+1_m)
 *
 * with c_i = 1 identically at an IS center, which is (22). Equation (25) is
 * self-normalizing, sum_n P^k_j(n,v) = sum_m L^k_{mk}(v) = 1, so no normalizing
 * constant is formed and every quantity in the recursion is positive: unlike
 * the load-dependent MVA of pfqn_mvald it CANNOT produce negative
 * probabilities and needs no stabilization. That is the practical reason to
 * prefer it.
 *
 * OUTPUT CONVENTIONS differ from pfqn_mvac and follow the load-dependent
 * family (pfqn_mvald, pfqn_dac): U is PER-STATION, 1 - P_j(0), because for a
 * load-dependent center the per-class product X_r L_{jr} is not the
 * utilization; and C is the per-class CYCLE TIME exclusive of think time,
 * N_r/X_r - Z_r, not an (M x R) residence time.
 *
 * Parts 2 and 3 are unchanged from pfqn_mvac, since (6) holds verbatim in the
 * presence of QLD centers, and are driven from the same chain setup.
 *
 * Arithmetic: EXACT-CAPABLE, field operations only.
 *
 * REFERENCE DEFECTS: none found. With mu identically one the results agree with
 * pfqn_mvac, and with a genuine load-dependent rate they agree with
 * pfqn_mvald.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_mvac.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_mvacld, mirroring [XN, QN, UN, CN, pij]. */
template <class T>
struct MvacldResult {
    std::vector<T> X;   ///< (R) per-class throughput
    Matrix<T> Q;        ///< (M x R) per-class queue length
    std::vector<T> U;   ///< (M) per-station utilization, 1 - P_j(0)
    std::vector<T> C;   ///< (R) per-class cycle time exclusive of think time
    Matrix<T> pij;      ///< (M x sumN+1) marginal queue-length probabilities
};

namespace detail {

/**
 * Level-k marginals over I_k, held as a (J1 * nvk) x (k+1) matrix with row
 * j * nvk + vloc, so that a plain Matrix carries the reference's three
 * dimensions without a bespoke tensor type.
 */
template <class T>
inline const T& pall_at(const Matrix<T>& P, std::size_t j, std::size_t vloc, std::size_t n,
                        std::size_t nvk) {
    return P(j * nvk + vloc, n);
}

/** Part 1 of the QLD basic step: (21)-(25) for k = k0..K over v in I_k. */
template <class T>
void mvacld_part1(std::size_t k0, const MvacSetup<T>& st, const Matrix<T>& MU,
                  std::vector<Matrix<T>>& Pall, std::vector<Matrix<T>>& Ljkall,
                  std::vector<std::vector<T>>& lamall) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t J = st.J, J1 = st.J1, K = st.K;
    for (std::size_t k = k0; k <= K; ++k) {
        const std::size_t t = K - k;
        const std::size_t nvk = st.cnt[t];
        const std::size_t nvp = st.cnt[t + 1];
        const Matrix<T>& Pp = Pall[k - 1];
        Matrix<T> Lk(J, nvk, zero), Pk(J1 * nvk, k + 1, zero);
        std::vector<T> lam(nvk, zero);
        for (std::size_t vloc = 0; vloc < nvk; ++vloc) {
            const std::size_t vi = st.off[t] + vloc;
            const std::vector<int>& v = st.Vlist[vi];
            std::vector<std::size_t> sloc(J, 0);
            for (std::size_t j = 0; j < J; ++j) {
                const long sj = st.succ(vi, j);
                if (sj < 0) throw NumericError("pfqn_mvacld: multiplicity vector out of range");
                sloc[j] = static_cast<std::size_t>(sj) - st.off[t + 1];
            }
            // (21): mean rate at which an SCSL chain pinned at center i is served
            std::vector<T> c(J, one);
            for (std::size_t i = 0; i < J1; ++i) {
                T ci = zero;
                for (std::size_t n = 0; n + 1 <= k; ++n) {
                    const std::size_t occ = n + static_cast<std::size_t>(v[i]) + 1;
                    if (occ > K) throw NumericError("pfqn_mvacld: rate index beyond sum(N)");
                    ci += pall_at(Pp, i, sloc[i], n, nvp) * MU(i, occ - 1) /
                          num_traits<T>::from_int(static_cast<long>(occ));
                }
                c[i] = ci;
            }
            // (23)-(24) in reference-station-free form
            T sw = zero;
            std::vector<T> w(J, zero);
            for (std::size_t j = 0; j < J; ++j) {
                if (st.a(j, k - 1) > zero) {
                    if (c[j] <= zero) throw NumericError("pfqn_mvacld: nonpositive service rate");
                    w[j] = st.a(j, k - 1) / c[j];
                }
                sw += w[j];
            }
            if (sw <= zero) throw NumericError("pfqn_mvacld: a chain has zero total demand");
            lam[vloc] = one / sw;
            for (std::size_t j = 0; j < J; ++j) Lk(j, vloc) = w[j] / sw;
            // (25): condition on the center holding the single chain-k customer
            for (std::size_t j = 0; j < J1; ++j) {
                for (std::size_t n = 0; n <= k; ++n) {
                    T s = zero;
                    if (n >= 1) s = Lk(j, vloc) * pall_at(Pp, j, sloc[j], n - 1, nvp);
                    if (n + 1 <= k)
                        for (std::size_t mm = 0; mm < J; ++mm)
                            if (mm != j) s += Lk(mm, vloc) * pall_at(Pp, j, sloc[mm], n, nvp);
                    Pk(j * nvk + vloc, n) = s;
                }
            }
        }
        Pall[k] = Pk;
        Ljkall[k] = Lk;
        lamall[k] = lam;
    }
}

}  // namespace detail

/**
 * @param L  (M x R) demands of the queue-length dependent centers
 * @param N  (R) closed populations
 * @param Z  (Mz x R) demands of the infinite-server centers
 * @param mu (M x n) load-dependent rates, mu(j, k-1) the total rate of center j
 *           with k jobs present; empty for the fixed-rate default
 */
template <class T>
MvacldResult<T> pfqn_mvacld(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                            const Matrix<T>& mu) {
    const std::size_t M = L.rows(), R = L.cols();
    if (M == 0 || R == 0) throw InputError("pfqn_mvacld: empty demand matrix");
    if (N.size() != R) throw InputError("pfqn_mvacld: L and N disagree on the class count");
    if (!Z.empty() && Z.cols() != R)
        throw InputError("pfqn_mvacld: the think time matrix and the demand matrix disagree");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    std::size_t K = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < 0) throw InputError("pfqn_mvacld: the population vector must be nonnegative");
        K += static_cast<std::size_t>(N[r]);
    }

    MvacldResult<T> res;
    res.X.assign(R, zero);
    res.Q = Matrix<T>(M, R, zero);
    res.U.assign(M, zero);
    res.C.assign(R, zero);
    res.pij = Matrix<T>(M, K + 1, zero);
    for (std::size_t i = 0; i < M; ++i) res.pij(i, 0) = one;  // an unvisited center holds no jobs
    if (K == 0) return res;

    if (!mu.empty()) {
        if (mu.rows() != M)
            throw InputError("pfqn_mvacld: the rate matrix and the demand matrix disagree on centers");
        if (mu.cols() < K)
            throw InputError("pfqn_mvacld: the rate matrix must supply a rate for every population up to sum(N)");
    }

    std::vector<std::size_t> ldIdx, isIdx;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (L(i, r) > zero) {
                ldIdx.push_back(i);
                break;
            }
    for (std::size_t i = 0; i < Z.rows(); ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (Z(i, r) > zero) {
                isIdx.push_back(i);
                break;
            }
    const std::size_t J1 = ldIdx.size(), J = J1 + isIdx.size();
    if (J == 0)
        throw InputError("pfqn_mvacld: all service demands are zero, the throughput is unbounded");

    Matrix<T> A(J, R, zero);
    for (std::size_t j = 0; j < J1; ++j)
        for (std::size_t r = 0; r < R; ++r) A(j, r) = L(ldIdx[j], r);
    for (std::size_t j = J1; j < J; ++j)
        for (std::size_t r = 0; r < R; ++r) A(j, r) = Z(isIdx[j - J1], r);

    Matrix<T> MU(J1, K, one);
    if (!mu.empty())
        for (std::size_t j = 0; j < J1; ++j)
            for (std::size_t n = 0; n < K; ++n) {
                if (mu(ldIdx[j], n) <= zero)
                    throw InputError(
                        "pfqn_mvacld: the service rates must be strictly positive for every "
                        "population up to sum(N)");
                MU(j, n) = mu(ldIdx[j], n);
            }

    detail::MvacSetup<T> st = detail::mvac_setup(A, N, J1, J, K);
    const std::size_t D = st.D, S = st.S;

    std::vector<Matrix<T>> Pall(K + 1), Ljkall(K + 1);
    std::vector<std::vector<T>> lamall(K + 1);
    Pall[0] = Matrix<T>(J1 * st.cnt[K], 1, one);  // P^0_j(0,v) = 1 over I_0

    std::vector<T> lamChain(K, zero);
    Matrix<T> Lchain(J, K, zero);

    detail::mvacld_part1(1, st, MU, Pall, Ljkall, lamall);
    lamChain[K - 1] = lamall[K][0];
    for (std::size_t j = 0; j < J; ++j) Lchain(j, K - 1) = Ljkall[K](j, 0);
    // The marginals of the ORIGINAL network live at k = K, v = 0; the label
    // interchanges of part 3 overwrite that level, so capture them now.
    for (std::size_t j = 0; j < J1; ++j)
        for (std::size_t n = 0; n <= K; ++n) res.pij(ldIdx[j], n) = Pall[K](j, n);

    // ---- part 2: chains that visit at least one IS center ---------------------------
    const long lmaxK = std::min(static_cast<long>(K) - 1, static_cast<long>(K - S));
    if (D >= 2 && lmaxK >= static_cast<long>(K - D + 1)) {
        std::vector<Matrix<T>> L2prev(K + 1), L2cur(K + 1);
        for (std::size_t k = K - D + 2; k <= K; ++k) {
            const std::size_t t = K - k;
            L2cur.assign(K + 1, Matrix<T>());
            const long lhi = std::min(static_cast<long>(k) - 1, static_cast<long>(K - S));
            for (long l = static_cast<long>(K - D + 1); l <= lhi; ++l) {
                Matrix<T> acc(J, st.cnt[t], zero);
                for (std::size_t vloc = 0; vloc < st.cnt[t]; ++vloc) {
                    const std::size_t vi = st.off[t] + vloc;
                    for (std::size_t j = 0; j < J; ++j) {
                        const long sj = st.succ(vi, j);
                        if (sj < 0)
                            throw NumericError("pfqn_mvacld: multiplicity vector out of range");
                        const std::size_t sloc = static_cast<std::size_t>(sj) - st.off[t + 1];
                        const Matrix<T>& prev =
                            (l == static_cast<long>(k) - 1) ? Ljkall[k - 1] : L2prev[l];
                        for (std::size_t i = 0; i < J; ++i)
                            acc(i, vloc) += Ljkall[k](j, vloc) * prev(i, sloc);
                    }
                }
                L2cur[l] = acc;
            }
            if (k == K) {
                for (long l = static_cast<long>(K - D + 1); l <= lmaxK; ++l) {
                    const std::size_t l0 = static_cast<std::size_t>(l) - 1;
                    for (std::size_t j = 0; j < J; ++j) Lchain(j, l0) = L2cur[l](j, 0);
                    std::size_t jIS = J;
                    for (std::size_t j = J1; j < J; ++j)
                        if (st.a(j, l0) > zero) {
                            jIS = j;
                            break;
                        }
                    if (jIS == J) throw NumericError("pfqn_mvacld: chain has no IS center");
                    lamChain[l0] = Lchain(jIS, l0) / st.a(jIS, l0);
                }
            }
            L2prev = L2cur;
        }
    }

    // ---- part 3: chains that visit no IS center, by label interchange ---------------
    std::vector<std::size_t> perm(K);
    for (std::size_t k = 0; k < K; ++k) perm[k] = k;
    for (std::size_t l = 1; l + 1 <= S; ++l) {
        std::swap(perm[K - l - 1], perm[K - 1]);
        for (std::size_t j = 0; j < J; ++j) {
            const T tmp = st.a(j, K - l - 1);
            st.a(j, K - l - 1) = st.a(j, K - 1);
            st.a(j, K - 1) = tmp;
        }
        detail::mvacld_part1(K - l, st, MU, Pall, Ljkall, lamall);
        lamChain[perm[K - 1]] = lamall[K][0];
        for (std::size_t j = 0; j < J; ++j) Lchain(j, perm[K - 1]) = Ljkall[K](j, 0);
    }

    // ---- expand the per-chain measures back to per-class ----------------------------
    for (std::size_t g = 0; g < D; ++g) {
        const std::size_t kg = K - D + g;
        for (std::size_t p = 0; p < st.posr.size(); ++p) {
            if (st.grpOfClass[p] != st.gorder[g]) continue;
            const std::size_t r = st.posr[p];
            const T nr = num_traits<T>::from_int(N[r]);
            res.X[r] = nr * lamChain[kg];
            for (std::size_t j = 0; j < J1; ++j) res.Q(ldIdx[j], r) = nr * Lchain(j, kg);
        }
    }
    for (std::size_t i = 0; i < M; ++i) res.U[i] = one - res.pij(i, 0);
    for (std::size_t p = 0; p < st.posr.size(); ++p) {
        const std::size_t r = st.posr[p];
        if (res.X[r] <= zero) throw NumericError("pfqn_mvacld: nonpositive throughput");
        T zr = zero;
        for (std::size_t i = 0; i < Z.rows(); ++i) zr += Z(i, r);
        res.C[r] = num_traits<T>::from_int(N[r]) / res.X[r] - zr;
    }
    return res;
}

/** Overload with the fixed-rate default. */
template <class T>
MvacldResult<T> pfqn_mvacld(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    return pfqn_mvacld(L, N, Z, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MVACLD_H
