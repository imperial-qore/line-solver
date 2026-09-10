/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_BMAPPHNN_RETRIAL_H
#define LINE_API_QSYS_QSYS_BMAPPHNN_RETRIAL_H

/**
 * The BMAP/PH/N/N bufferless retrial queue with flexible retrial admission
 * control.
 *
 * Templated port of matlab/src/api/qsys/qsys_bmapphnn_retrial.m, which
 * implements Dudin et al., "Analysis of BMAP/PH/N-Type Queueing System with
 * Flexible Retrials Admission Control", Mathematics 2025, 13(9), 1434.
 *
 * THE MODEL. N servers and no waiting room. Arrivals come in batches from a
 * BMAP (D0, D1, ..., DK) on V states; a batch that finds fewer than its size
 * free servers fills what it can and the excess either joins the orbit, with
 * probability 1 - p, or is lost, with probability p. Service is PH (beta, S)
 * on M phases, so the server-side state with n busy servers is the multiset of
 * their phases, of which there are T_n = C(n+M-1, M-1); the level of the
 * generator is the orbit population and each level carries V * sum_n T_n
 * states. Each orbiting customer retries at rate alpha and abandons at rate
 * gamma, and a retrial SUCCEEDS only while n <= R(nu), which is the admission
 * control: above the threshold the retrial finds the system closed and the
 * customer stays in orbit.
 *
 * THE STATE SPACE IS A LEVEL-DEPENDENT QBD and the reference does not solve it
 * as one: it truncates the orbit at a finite level, assembles the whole
 * generator densely, and solves pi Q = 0 with the last column replaced by the
 * normalization. This port does the same, so that the two agree, including on
 * the truncation. See the two reference defects below, which are both about
 * that truncation.
 *
 * WHAT THE PORT CHANGES, AND WHY IT DOES NOT CHANGE THE ANSWER. The reference
 * rebuilds B, B_bar, Gamma and B_tilde inside buildGeneratorLevel, i.e. once
 * per (i, j) block pair, although none of them depends on i or j; on the
 * default truncation that is over a hundred rebuilds of the same matrices.
 * They are hoisted here. The generator entries are identical.
 *
 * REFERENCE DEFECTS in qsys_bmapphnn_retrial.m:
 *
 *  1. AN UNSTABLE MODEL RETURNS A FINITE, PLAUSIBLE-LOOKING ANSWER. The
 *     truncated chain is always positive recurrent, and line 202 renormalizes
 *     it (pi = pi / sum(pi)), so divergence appears as a number rather than as
 *     an error. MATLAB reproduction, from matlab/:
 *         for L = [100 200 400 800 1600]
 *           r = qsys_bmapphnn_retrial({-2,2}, 1, -1, 3, 0.5, 0, 0, 0, 'MaxLevel', L);
 *           fprintf('%d %g %g\n', L, r.L_orbit, sum(r.pi(end,:)));
 *         end
 *     gives L_orbit 94.5, 194.0, 393.7, 793.6, 1593.5 -- growing linearly with
 *     the truncation level -- while the mass sitting at the TOP level stays at
 *     about 0.18 instead of decaying, and the reported throughput drifts from
 *     1.8607 to 1.8741. With R = 0 a retrial succeeds only when the system is
 *     completely empty, so the orbit is not stable even though the offered
 *     load rho = 2/3 is well below one. Nothing warns.
 *  2. THE DEFAULT TRUNCATION LEVEL IGNORES EVERY ORBIT PARAMETER. Line 127
 *     sets truncLevel = max(100, ceil(50/(1 - min(rho, 0.99)))) with
 *     rho = lambda b1 / N, so it depends on neither alpha nor gamma nor the
 *     blocking probability, which are exactly what govern the decay of the
 *     orbit tail. MATLAB reproduction: with alpha = 0.02 the default level is
 *     150, and L_orbit converges to 57.73566168 (reached by MaxLevel 400),
 *     while MaxLevel 100 gives 57.5935191, an error of 0.25 per cent with no
 *     indication that anything was truncated.
 *
 *     Neither is worked around here. The port returns truncLevel and, as the
 *     diagnostic the reference lacks, topLevelMass: the probability mass at
 *     the highest retained level. A caller can test it (it should be
 *     negligible, of the order of 1e-16 on a converged instance) and raise
 *     maxLevel when it is not. The default level and the returned means are
 *     the reference's, unchanged.
 *
 *  3. THE DOCUMENTED MEANING OF R DOES NOT MATCH THE CODE. The header of the
 *     .m file says "When n > R(nu), arriving customers go to orbit", but R
 *     enters the generator only through computeGamma and computeBbar, both of
 *     which act on RETRIALS from the orbit. A fresh arrival always takes a
 *     free server if there is one, whatever R is. The code is consistent with
 *     the paper's "flexible retrials admission control"; the docstring is not.
 *
 *  4. computeC's first branch builds a (T_n x 1) zero column, a shape that
 *     cannot be a valid block, and the caller silently discards it through a
 *     size test (size(C_nk, 2) == T(N+1)). It also assigns the same zero
 *     matrix twice under an if that cannot change it. Harmless, and the port
 *     simply does not create the malformed shape.
 *
 * ARITHMETIC. Nothing here needs a transcendental function: the generator is a
 * rational expression in the inputs and the solution is one linear solve, so
 * the header is instantiable at T = Rational and the exact instantiation
 * returns the stationary law of the TRUNCATED chain exactly. That is worth
 * having for a small truncation, where it makes the generator itself
 * checkable, and impractical for the default one, where the dense solve is
 * hundreds of dimensions wide.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Options of qsys_bmapphnn_retrial. */
struct BmapPhNnRetrialOptions {
    /** Orbit truncation level; 0 selects the reference's automatic choice. */
    std::size_t maxLevel = 0;
};

/** Result of qsys_bmapphnn_retrial. */
template <class T>
struct BmapPhNnRetrialResult {
    T L_orbit;         ///< mean number of customers in orbit
    T N_server;        ///< mean number of busy servers
    T L_system;        ///< L_orbit + N_server
    T utilization;     ///< N_server / N
    T throughput;      ///< N_server / b1
    T P_idle;          ///< probability that every server is idle
    T P_empty_orbit;   ///< probability that the orbit is empty
    T P_empty_system;  ///< probability that the system is empty
    Matrix<T> pi;      ///< stationary law, (truncLevel + 1) x (V d)
    std::size_t truncLevel;  ///< truncation level used
    T topLevelMass;    ///< mass at the highest retained level; see defect 2
    bool clipped;      ///< a negative probability had to be clipped to zero
};

namespace detail {

/**
 * Weak compositions of n into M parts in the reference's reverse
 * lexicographic order (generateCompositions in the .m file). The order matters:
 * it fixes the index of every server-phase multiset.
 */
inline std::vector<std::vector<int>> retrial_compositions(int n, int M) {
    std::vector<std::vector<int>> out;
    if (M <= 0) return out;
    if (M == 1) {
        out.push_back(std::vector<int>(1, n));
        return out;
    }
    for (int m1 = n; m1 >= 0; --m1) {
        const std::vector<std::vector<int>> sub = retrial_compositions(n - m1, M - 1);
        for (std::size_t s = 0; s < sub.size(); ++s) {
            std::vector<int> row;
            row.reserve(static_cast<std::size_t>(M));
            row.push_back(m1);
            row.insert(row.end(), sub[s].begin(), sub[s].end());
            out.push_back(row);
        }
    }
    return out;
}

/** Index of a composition in the reverse-lexicographic list, or -1. */
inline int retrial_find(const std::vector<std::vector<int>>& comps, const std::vector<int>& key) {
    for (std::size_t j = 0; j < comps.size(); ++j)
        if (comps[j] == key) return static_cast<int>(j);
    return -1;
}

/** Everything the block builders need, assembled once. */
template <class T>
struct RetrialCtx {
    std::vector<Matrix<T>> D;  ///< D0 ... DK
    std::vector<T> beta;
    Matrix<T> S;
    std::vector<T> S0;
    int M = 0, N = 0, V = 0, K = 0;
    std::size_t d = 0;
    std::vector<std::size_t> Tn;   ///< T_n, n = 0..N
    std::vector<long> R;           ///< admission threshold per BMAP state
    T alpha, gamma, p;
    std::vector<std::vector<std::vector<int>>> comps;  ///< comps[n]
};

/** Starting index of the states with n busy servers, 0-based. */
template <class T>
std::size_t retrial_offset(const RetrialCtx<T>& c, int n) {
    std::size_t off = 0;
    for (int i = 0; i < n; ++i) off += c.Tn[static_cast<std::size_t>(i)];
    return off;
}

/** L_n: service completions, T_n x T_{n-1}. */
template <class T>
Matrix<T> retrial_L(const RetrialCtx<T>& c, int n) {
    const T zero = num_traits<T>::from_int(0);
    if (n == 0) return Matrix<T>();
    Matrix<T> L(c.Tn[static_cast<std::size_t>(n)], c.Tn[static_cast<std::size_t>(n - 1)], zero);
    const std::vector<std::vector<int>>& cn = c.comps[static_cast<std::size_t>(n)];
    const std::vector<std::vector<int>>& cm = c.comps[static_cast<std::size_t>(n - 1)];
    for (std::size_t i = 0; i < cn.size(); ++i)
        for (int l = 0; l < c.M; ++l) {
            if (cn[i][static_cast<std::size_t>(l)] <= 0) continue;
            std::vector<int> mp = cn[i];
            --mp[static_cast<std::size_t>(l)];
            const int j = retrial_find(cm, mp);
            if (j >= 0)
                L(i, static_cast<std::size_t>(j)) +=
                    num_traits<T>::from_int(cn[i][static_cast<std::size_t>(l)]) *
                    c.S0[static_cast<std::size_t>(l)];
        }
    return L;
}

/** A_n: service phase changes, T_n x T_n. */
template <class T>
Matrix<T> retrial_A(const RetrialCtx<T>& c, int n) {
    const T zero = num_traits<T>::from_int(0);
    if (n == 0) return Matrix<T>(1, 1, zero);
    Matrix<T> A(c.Tn[static_cast<std::size_t>(n)], c.Tn[static_cast<std::size_t>(n)], zero);
    const std::vector<std::vector<int>>& cn = c.comps[static_cast<std::size_t>(n)];
    for (std::size_t i = 0; i < cn.size(); ++i)
        for (int l = 0; l < c.M; ++l) {
            if (cn[i][static_cast<std::size_t>(l)] <= 0) continue;
            for (int lp = 0; lp < c.M; ++lp) {
                if (lp == l) continue;
                if (!(c.S(static_cast<std::size_t>(l), static_cast<std::size_t>(lp)) > zero))
                    continue;
                std::vector<int> mp = cn[i];
                --mp[static_cast<std::size_t>(l)];
                ++mp[static_cast<std::size_t>(lp)];
                const int j = retrial_find(cn, mp);
                if (j >= 0)
                    A(i, static_cast<std::size_t>(j)) +=
                        num_traits<T>::from_int(cn[i][static_cast<std::size_t>(l)]) *
                        c.S(static_cast<std::size_t>(l), static_cast<std::size_t>(lp));
            }
        }
    return A;
}

/** P_n: a new customer enters service, T_n x T_{n+1}. */
template <class T>
Matrix<T> retrial_P(const RetrialCtx<T>& c, int n) {
    const T zero = num_traits<T>::from_int(0);
    if (n >= c.N) return Matrix<T>();
    Matrix<T> P(c.Tn[static_cast<std::size_t>(n)], c.Tn[static_cast<std::size_t>(n + 1)], zero);
    const std::vector<std::vector<int>>& cn = c.comps[static_cast<std::size_t>(n)];
    const std::vector<std::vector<int>>& cp = c.comps[static_cast<std::size_t>(n + 1)];
    for (std::size_t i = 0; i < cn.size(); ++i)
        for (int l = 0; l < c.M; ++l) {
            if (!(c.beta[static_cast<std::size_t>(l)] > zero)) continue;
            std::vector<int> mp = cn[i];
            ++mp[static_cast<std::size_t>(l)];
            const int j = retrial_find(cp, mp);
            if (j >= 0)
                P(i, static_cast<std::size_t>(j)) += c.beta[static_cast<std::size_t>(l)];
        }
    return P;
}

/** Delta_n: total exit rate of each server-phase multiset, diagonal T_n. */
template <class T>
Matrix<T> retrial_Delta(const RetrialCtx<T>& c, int n) {
    const T zero = num_traits<T>::from_int(0);
    if (n == 0) return Matrix<T>(1, 1, zero);
    Matrix<T> D(c.Tn[static_cast<std::size_t>(n)], c.Tn[static_cast<std::size_t>(n)], zero);
    const std::vector<std::vector<int>>& cn = c.comps[static_cast<std::size_t>(n)];
    for (std::size_t i = 0; i < cn.size(); ++i) {
        T total = zero;
        for (int l = 0; l < c.M; ++l)
            total += num_traits<T>::from_int(cn[i][static_cast<std::size_t>(l)]) *
                     T(-c.S(static_cast<std::size_t>(l), static_cast<std::size_t>(l)));
        D(i, i) = total;
    }
    return D;
}

/** Gamma^(nu): the indicator of n > R(nu), diagonal d. */
template <class T>
Matrix<T> retrial_Gamma(const RetrialCtx<T>& c, int nu) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> G(c.d, c.d, zero);
    std::size_t off = 0;
    for (int n = 0; n <= c.N; ++n) {
        const std::size_t width = c.Tn[static_cast<std::size_t>(n)];
        if (static_cast<long>(n) > c.R[static_cast<std::size_t>(nu)])
            for (std::size_t t = 0; t < width; ++t) G(off + t, off + t) = one;
        off += width;
    }
    return G;
}

/** The scalar of G_{n,n}^{(nu,nu')}: batch losses when the batch cannot fit. */
template <class T>
T retrial_G_scalar(const RetrialCtx<T>& c, int n, int nu, int nuPrime) {
    const T zero = num_traits<T>::from_int(0);
    if (n <= c.N - c.K) return zero;
    T total = zero;
    for (int k = c.N - n + 1; k <= c.K; ++k)
        if (k >= 1)
            total += c.D[static_cast<std::size_t>(k)](static_cast<std::size_t>(nu),
                                                      static_cast<std::size_t>(nuPrime));
    return T(c.p * total);
}

/**
 * The product P_n P_{n+1} ... P_{n+k-1}, which places k arriving customers
 * into service starting from n busy servers.
 */
template <class T>
Matrix<T> retrial_Pprod(const RetrialCtx<T>& c, const std::vector<Matrix<T>>& P, int n, int upto) {
    Matrix<T> prod = eye<T>(c.Tn[static_cast<std::size_t>(n)]);
    for (int j = n; j < upto; ++j)
        if (j < c.N) prod = matmul(prod, P[static_cast<std::size_t>(j)]);
    return prod;
}

/** B^(nu): the within-level block for a BMAP state that does not change. */
template <class T>
Matrix<T> retrial_B(const RetrialCtx<T>& c, const std::vector<Matrix<T>>& L,
                    const std::vector<Matrix<T>>& A, const std::vector<Matrix<T>>& P,
                    const std::vector<Matrix<T>>& Delta, int nu) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> B(c.d, c.d, zero);
    for (int n = 0; n <= c.N; ++n) {
        const std::size_t rs = retrial_offset(c, n);
        const std::size_t w = c.Tn[static_cast<std::size_t>(n)];
        const T g = retrial_G_scalar(c, n, nu, nu);
        if (n == 0) {
            B(rs, rs) = g;
        } else {
            for (std::size_t a = 0; a < w; ++a)
                for (std::size_t b = 0; b < w; ++b)
                    B(rs + a, rs + b) = A[static_cast<std::size_t>(n)](a, b) +
                                        Delta[static_cast<std::size_t>(n)](a, b) +
                                        (a == b ? g : zero);
        }
        if (n >= 1) {
            const std::size_t cs = retrial_offset(c, n - 1);
            for (std::size_t a = 0; a < w; ++a)
                for (std::size_t b = 0; b < c.Tn[static_cast<std::size_t>(n - 1)]; ++b)
                    B(rs + a, cs + b) = L[static_cast<std::size_t>(n)](a, b);
        }
        for (int k = 1; k <= c.K; ++k) {
            if (n + k > c.N) continue;
            const std::size_t cs = retrial_offset(c, n + k);
            const T dk = c.D[static_cast<std::size_t>(k)](static_cast<std::size_t>(nu),
                                                          static_cast<std::size_t>(nu));
            const Matrix<T> pp = retrial_Pprod(c, P, n, n + k);
            for (std::size_t a = 0; a < w; ++a)
                for (std::size_t b = 0; b < c.Tn[static_cast<std::size_t>(n + k)]; ++b)
                    B(rs + a, cs + b) = dk * pp(a, b);
        }
    }
    return B;
}

/** B_bar^(nu): a successful retrial, valid only while n <= R(nu). */
template <class T>
Matrix<T> retrial_Bbar(const RetrialCtx<T>& c, const std::vector<Matrix<T>>& P, int nu) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> B(c.d, c.d, zero);
    const long cap = c.R[static_cast<std::size_t>(nu)] < static_cast<long>(c.N - 1)
                         ? c.R[static_cast<std::size_t>(nu)]
                         : static_cast<long>(c.N - 1);
    for (long n = 0; n <= cap; ++n) {
        const std::size_t rs = retrial_offset(c, static_cast<int>(n));
        const std::size_t cs = retrial_offset(c, static_cast<int>(n) + 1);
        const Matrix<T>& Pn = P[static_cast<std::size_t>(n)];
        for (std::size_t a = 0; a < Pn.rows(); ++a)
            for (std::size_t b = 0; b < Pn.cols(); ++b) B(rs + a, cs + b) = Pn(a, b);
    }
    return B;
}

/** B_tilde^(nu,nu'): the within-level block when the BMAP state changes. */
template <class T>
Matrix<T> retrial_Btilde(const RetrialCtx<T>& c, const std::vector<Matrix<T>>& P, int nu,
                         int nuPrime) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> B(c.d, c.d, zero);
    for (int n = 0; n <= c.N; ++n) {
        const std::size_t rs = retrial_offset(c, n);
        const std::size_t w = c.Tn[static_cast<std::size_t>(n)];
        const T g = retrial_G_scalar(c, n, nu, nuPrime);
        for (std::size_t a = 0; a < w; ++a) B(rs + a, rs + a) = g;
        for (int k = 1; k <= c.K; ++k) {
            if (n + k > c.N) continue;
            const std::size_t cs = retrial_offset(c, n + k);
            const T dk = c.D[static_cast<std::size_t>(k)](static_cast<std::size_t>(nu),
                                                          static_cast<std::size_t>(nuPrime));
            const Matrix<T> pp = retrial_Pprod(c, P, n, n + k);
            for (std::size_t a = 0; a < w; ++a)
                for (std::size_t b = 0; b < c.Tn[static_cast<std::size_t>(n + k)]; ++b)
                    B(rs + a, cs + b) = dk * pp(a, b);
        }
    }
    return B;
}

/**
 * C_{n,k}^(nu,nu'): a batch that overflows the free servers sends k customers
 * to the orbit and leaves all N servers busy. Returns an empty matrix when the
 * block does not exist, in place of the reference's malformed zero column.
 */
template <class T>
Matrix<T> retrial_C(const RetrialCtx<T>& c, const std::vector<Matrix<T>>& P, int n, int k, int nu,
                    int nuPrime) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (n < c.N - c.K + k) return Matrix<T>();
    if (n < c.N) {
        const int batch = c.N - n + k;
        if (batch < 1 || batch > c.K) return Matrix<T>();
        const T db = c.D[static_cast<std::size_t>(batch)](static_cast<std::size_t>(nu),
                                                          static_cast<std::size_t>(nuPrime));
        Matrix<T> pp = retrial_Pprod(c, P, n, c.N);
        for (std::size_t a = 0; a < pp.rows(); ++a)
            for (std::size_t b = 0; b < pp.cols(); ++b) pp(a, b) = T(one - c.p) * db * pp(a, b);
        return pp;
    }
    if (k < 1 || k > c.K) return Matrix<T>();
    const T dk = c.D[static_cast<std::size_t>(k)](static_cast<std::size_t>(nu),
                                                  static_cast<std::size_t>(nuPrime));
    Matrix<T> out(c.Tn[static_cast<std::size_t>(c.N)], c.Tn[static_cast<std::size_t>(c.N)], zero);
    for (std::size_t a = 0; a < out.rows(); ++a) out(a, a) = T(one - c.p) * dk;
    return out;
}

}  // namespace detail

/**
 * The BMAP/PH/N/N bufferless retrial queue.
 *
 * @param D     the BMAP as {D0, D1, ..., DK}, each V x V
 * @param beta  PH service entry vector, length M
 * @param S     PH service sub-generator, M x M
 * @param N     number of servers, which is also the capacity
 * @param alpha retrial rate per orbiting customer
 * @param gamma abandonment rate per orbiting customer
 * @param p     probability that an overflowing batch is lost rather than
 *              joining the orbit
 * @param R     admission threshold per BMAP state; a single entry is broadcast
 * @param opt   truncation level
 */
template <class T>
BmapPhNnRetrialResult<T> qsys_bmapphnn_retrial(const std::vector<Matrix<T>>& D,
                                               const std::vector<T>& beta, const Matrix<T>& S,
                                               int N, const T& alpha, const T& gamma, const T& p,
                                               const std::vector<long>& R,
                                               const BmapPhNnRetrialOptions& opt) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (D.size() < 2) throw InputError("qsys_bmapphnn_retrial: the BMAP needs D0 and at least D1");
    if (N < 1) throw InputError("qsys_bmapphnn_retrial: at least one server is required");
    const int V = static_cast<int>(D[0].rows());
    if (V < 1) throw InputError("qsys_bmapphnn_retrial: empty BMAP");
    for (std::size_t k = 0; k < D.size(); ++k)
        if (D[k].rows() != static_cast<std::size_t>(V) || D[k].cols() != static_cast<std::size_t>(V))
            throw InputError("qsys_bmapphnn_retrial: every BMAP matrix must be V x V");
    const int M = static_cast<int>(S.rows());
    if (M < 1 || S.cols() != static_cast<std::size_t>(M) ||
        beta.size() != static_cast<std::size_t>(M))
        throw InputError("qsys_bmapphnn_retrial: the PH service is malformed");
    if (R.empty()) throw InputError("qsys_bmapphnn_retrial: no admission threshold given");

    detail::RetrialCtx<T> c;
    c.D = D;
    c.beta = beta;
    c.S = S;
    c.M = M;
    c.N = N;
    c.V = V;
    c.K = static_cast<int>(D.size()) - 1;
    c.alpha = alpha;
    c.gamma = gamma;
    c.p = p;
    c.R.assign(static_cast<std::size_t>(V), R[0]);
    if (R.size() > 1) {
        if (R.size() != static_cast<std::size_t>(V))
            throw InputError("qsys_bmapphnn_retrial: R must be a scalar or one entry per BMAP "
                             "state");
        c.R = R;
    }

    // S0 = -S e, the absorption rates
    c.S0.assign(static_cast<std::size_t>(M), zero);
    for (int i = 0; i < M; ++i) {
        T s = zero;
        for (int j = 0; j < M; ++j)
            s += S(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
        c.S0[static_cast<std::size_t>(i)] = -s;
    }

    // T_n = C(n + M - 1, M - 1), and the compositions that index those states
    c.Tn.assign(static_cast<std::size_t>(N) + 1, 0);
    c.comps.resize(static_cast<std::size_t>(N) + 1);
    c.d = 0;
    for (int n = 0; n <= N; ++n) {
        c.comps[static_cast<std::size_t>(n)] = detail::retrial_compositions(n, M);
        c.Tn[static_cast<std::size_t>(n)] = c.comps[static_cast<std::size_t>(n)].size();
        c.d += c.Tn[static_cast<std::size_t>(n)];
    }

    // mean arrival rate and mean service time, for rho and the throughput
    Matrix<T> Dsum(static_cast<std::size_t>(V), static_cast<std::size_t>(V), zero);
    Matrix<T> DkSum(static_cast<std::size_t>(V), static_cast<std::size_t>(V), zero);
    for (std::size_t k = 0; k < D.size(); ++k)
        for (std::size_t i = 0; i < Dsum.rows(); ++i)
            for (std::size_t j = 0; j < Dsum.cols(); ++j) {
                Dsum(i, j) += D[k](i, j);
                if (k >= 1)
                    DkSum(i, j) += num_traits<T>::from_int(static_cast<long>(k)) * D[k](i, j);
            }
    // theta Dsum = 0, sum theta = 1, by replacing the last equation
    Matrix<T> A = Dsum.transpose();
    for (std::size_t j = 0; j < A.cols(); ++j) A(A.rows() - 1, j) = one;
    std::vector<T> rhs(A.rows(), zero);
    rhs[A.rows() - 1] = one;
    const std::vector<T> theta = solve(A, rhs);
    T lambda = zero;
    for (std::size_t i = 0; i < theta.size(); ++i)
        for (std::size_t j = 0; j < DkSum.cols(); ++j) lambda += theta[i] * DkSum(i, j);

    Matrix<T> negS = S;
    for (std::size_t i = 0; i < negS.rows(); ++i)
        for (std::size_t j = 0; j < negS.cols(); ++j) negS(i, j) = -negS(i, j);
    const std::vector<T> tau = mulvec(inverse(negS), ones<T>(static_cast<std::size_t>(M)));
    T b1 = zero;
    for (int i = 0; i < M; ++i) b1 += beta[static_cast<std::size_t>(i)] * tau[static_cast<std::size_t>(i)];
    if (!(b1 > zero)) throw InputError("qsys_bmapphnn_retrial: the mean service time must be positive");

    // truncation level, the reference's formula
    std::size_t truncLevel = opt.maxLevel;
    if (truncLevel == 0) {
        const double rho = num_traits<T>::to_double(T(lambda * b1)) / static_cast<double>(N);
        const double capped = rho < 0.99 ? rho : 0.99;
        const double lv = std::ceil(50.0 / (1.0 - capped));
        truncLevel = lv > 100.0 ? static_cast<std::size_t>(lv) : static_cast<std::size_t>(100);
    }

    // Blocks that do not depend on the level, built once (the reference
    // rebuilds them per block pair).
    std::vector<Matrix<T>> L(static_cast<std::size_t>(N) + 1), Am(static_cast<std::size_t>(N) + 1),
        P(static_cast<std::size_t>(N) + 1), Delta(static_cast<std::size_t>(N) + 1);
    for (int n = 0; n <= N; ++n) {
        L[static_cast<std::size_t>(n)] = detail::retrial_L(c, n);
        Am[static_cast<std::size_t>(n)] = detail::retrial_A(c, n);
        P[static_cast<std::size_t>(n)] = detail::retrial_P(c, n);
        Delta[static_cast<std::size_t>(n)] = detail::retrial_Delta(c, n);
    }
    std::vector<Matrix<T>> B(static_cast<std::size_t>(V)), Bbar(static_cast<std::size_t>(V)),
        Gam(static_cast<std::size_t>(V));
    for (int nu = 0; nu < V; ++nu) {
        B[static_cast<std::size_t>(nu)] = detail::retrial_B(c, L, Am, P, Delta, nu);
        Bbar[static_cast<std::size_t>(nu)] = detail::retrial_Bbar(c, P, nu);
        Gam[static_cast<std::size_t>(nu)] = detail::retrial_Gamma(c, nu);
    }
    std::vector<std::vector<Matrix<T>>> Btilde(
        static_cast<std::size_t>(V), std::vector<Matrix<T>>(static_cast<std::size_t>(V)));
    for (int nu = 0; nu < V; ++nu)
        for (int nup = 0; nup < V; ++nup)
            if (nu != nup) Btilde[static_cast<std::size_t>(nu)][static_cast<std::size_t>(nup)] =
                    detail::retrial_Btilde(c, P, nu, nup);

    const std::size_t Vd = static_cast<std::size_t>(V) * c.d;
    const std::size_t total = (truncLevel + 1) * Vd;
    Matrix<T> Q(total, total, zero);

    for (std::size_t i = 0; i <= truncLevel; ++i) {
        const T iT = num_traits<T>::from_int(static_cast<long>(i));
        // diagonal block
        for (int nu = 0; nu < V; ++nu) {
            const std::size_t rs = i * Vd + static_cast<std::size_t>(nu) * c.d;
            for (int nup = 0; nup < V; ++nup) {
                const std::size_t cs = i * Vd + static_cast<std::size_t>(nup) * c.d;
                const T d0 = D[0](static_cast<std::size_t>(nu), static_cast<std::size_t>(nup));
                if (nu == nup) {
                    for (std::size_t a = 0; a < c.d; ++a)
                        for (std::size_t b = 0; b < c.d; ++b) {
                            T v = B[static_cast<std::size_t>(nu)](a, b) +
                                  iT * alpha * Gam[static_cast<std::size_t>(nu)](a, b);
                            if (a == b) v += d0 - iT * T(gamma + alpha);
                            Q(rs + a, cs + b) = v;
                        }
                } else {
                    const Matrix<T>& bt =
                        Btilde[static_cast<std::size_t>(nu)][static_cast<std::size_t>(nup)];
                    for (std::size_t a = 0; a < c.d; ++a)
                        for (std::size_t b = 0; b < c.d; ++b)
                            Q(rs + a, cs + b) = bt(a, b) + (a == b ? d0 : zero);
                }
            }
        }
        // subdiagonal block: abandonment and successful retrials
        if (i >= 1) {
            for (int nu = 0; nu < V; ++nu) {
                const std::size_t rs = i * Vd + static_cast<std::size_t>(nu) * c.d;
                const std::size_t cs = (i - 1) * Vd + static_cast<std::size_t>(nu) * c.d;
                for (std::size_t a = 0; a < c.d; ++a)
                    for (std::size_t b = 0; b < c.d; ++b) {
                        T v = iT * alpha * Bbar[static_cast<std::size_t>(nu)](a, b);
                        if (a == b) v += iT * gamma;
                        Q(rs + a, cs + b) = v;
                    }
            }
        }
        // superdiagonal blocks: batches that overflow into the orbit
        for (int k = 1; k <= c.K; ++k) {
            const std::size_t j = i + static_cast<std::size_t>(k);
            if (j > truncLevel) break;
            for (int nu = 0; nu < V; ++nu) {
                const std::size_t rs = i * Vd + static_cast<std::size_t>(nu) * c.d;
                for (int nup = 0; nup < V; ++nup) {
                    const std::size_t cs = j * Vd + static_cast<std::size_t>(nup) * c.d;
                    const std::size_t ncol = detail::retrial_offset(c, N);
                    for (int n = 0; n <= N; ++n) {
                        const Matrix<T> Cnk = detail::retrial_C(c, P, n, k, nu, nup);
                        if (Cnk.rows() == 0) continue;
                        if (Cnk.cols() != c.Tn[static_cast<std::size_t>(N)]) continue;
                        const std::size_t nrow = detail::retrial_offset(c, n);
                        for (std::size_t a = 0; a < Cnk.rows(); ++a)
                            for (std::size_t b = 0; b < Cnk.cols(); ++b)
                                Q(rs + nrow + a, cs + ncol + b) += Cnk(a, b);
                    }
                }
            }
        }
    }

    // conservative diagonal, as the reference does after assembly
    for (std::size_t i = 0; i < total; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < total; ++j) s += Q(i, j);
        Q(i, i) = Q(i, i) - s;
    }

    // pi Q = 0 with the last column carrying the normalization
    for (std::size_t i = 0; i < total; ++i) Q(i, total - 1) = one;
    Matrix<T> QT = Q.transpose();
    std::vector<T> rhs2(total, zero);
    rhs2[total - 1] = one;
    std::vector<T> pi = solve(QT, rhs2);

    BmapPhNnRetrialResult<T> r;
    r.clipped = false;
    T mass = zero;
    for (std::size_t i = 0; i < total; ++i) {
        if (pi[i] < zero) {
            pi[i] = zero;
            r.clipped = true;
        }
        mass += pi[i];
    }
    if (!(mass > zero)) throw NumericError("qsys_bmapphnn_retrial: the solution has no mass");
    for (std::size_t i = 0; i < total; ++i) pi[i] = pi[i] / mass;

    r.pi = Matrix<T>(truncLevel + 1, Vd);
    for (std::size_t i = 0; i <= truncLevel; ++i)
        for (std::size_t j = 0; j < Vd; ++j) r.pi(i, j) = pi[i * Vd + j];

    r.L_orbit = zero;
    for (std::size_t i = 1; i <= truncLevel; ++i) {
        T lv = zero;
        for (std::size_t j = 0; j < Vd; ++j) lv += r.pi(i, j);
        r.L_orbit += num_traits<T>::from_int(static_cast<long>(i)) * lv;
    }
    r.N_server = zero;
    r.P_idle = zero;
    for (std::size_t i = 0; i <= truncLevel; ++i)
        for (int nu = 0; nu < V; ++nu) {
            const std::size_t base = static_cast<std::size_t>(nu) * c.d;
            for (int n = 0; n <= N; ++n) {
                const std::size_t off = base + detail::retrial_offset(c, n);
                for (std::size_t t = 0; t < c.Tn[static_cast<std::size_t>(n)]; ++t)
                    r.N_server += num_traits<T>::from_int(n) * r.pi(i, off + t);
            }
            r.P_idle += r.pi(i, base);
        }
    r.P_empty_orbit = zero;
    for (std::size_t j = 0; j < Vd; ++j) r.P_empty_orbit += r.pi(0, j);
    r.P_empty_system = zero;
    for (int nu = 0; nu < V; ++nu) r.P_empty_system += r.pi(0, static_cast<std::size_t>(nu) * c.d);
    r.topLevelMass = zero;
    for (std::size_t j = 0; j < Vd; ++j) r.topLevelMass += r.pi(truncLevel, j);

    r.L_system = T(r.L_orbit + r.N_server);
    r.utilization = T(r.N_server / num_traits<T>::from_int(N));
    r.throughput = T(r.N_server / b1);
    r.truncLevel = truncLevel;
    return r;
}

/** qsys_bmapphnn_retrial with the reference's automatic truncation level. */
template <class T>
BmapPhNnRetrialResult<T> qsys_bmapphnn_retrial(const std::vector<Matrix<T>>& D,
                                               const std::vector<T>& beta, const Matrix<T>& S,
                                               int N, const T& alpha, const T& gamma, const T& p,
                                               long R) {
    return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, std::vector<long>(1, R),
                                 BmapPhNnRetrialOptions());
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_BMAPPHNN_RETRIAL_H
