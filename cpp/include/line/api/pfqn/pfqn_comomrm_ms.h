/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_COMOMRM_MS_H
#define LINE_API_PFQN_COMOMRM_MS_H

/**
 * CoMoM for the MULTISERVER repairman model: one queueing station with S
 * servers (optionally replicated m times), plus a delay.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_comomrm_ms.m, and the home of the
 * bidiagonal transfer-matrix recursion shared with pfqn_comomrm_ld.
 *
 * The basis is the vector of normalizing constants indexed by the queue
 * occupancy, h_k = G(k jobs at the queue), and adding one class-r job applies
 * the bidiagonal transfer matrix
 *
 *   T_r = Z_r I + superdiag_k ( L_r (Nt+1-k) / mu(Nt+1-k) ),
 *   h  <- T_r h / n_r,
 *
 * once per job, for n_r = 1, ..., N_r. G(N) is the sum of the resulting vector
 * and the queue-length marginal is its reversal, normalized.
 *
 * SCALING. The reference renormalizes h to unit 1-norm after every step and
 * accumulates the discarded factors in log space, because in IEEE double the
 * unscaled vector underflows. That renormalization is a pure change of
 * representation: the discarded factors multiply back to exactly the sum of
 * the unscaled vector, so this port drops it and returns
 *
 *   G = Gremaind * sum_k h_k
 *
 * with h_k the UNSCALED basis. In an exact field the two agree identically; in
 * double they agree to rounding. The marginal is unaffected either way, since
 * it is a ratio within one vector.
 *
 * Arithmetic: EXACT-CAPABLE. Only additions, multiplications and divisions in
 * the field of the inputs. The multiserver rate lattice comes from pfqn_mu_ms
 * (m > 1) or is min(k, S) (m = 1), both exact.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_mu_ms.h"
#include "line/api/pfqn/pfqn_nc_sanitize.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
struct ComomRmResult {
    T G;                 ///< normalizing constant
    double lG;           ///< its logarithm
    std::vector<T> prob;  ///< (Nt+1) queue-length marginal, prob[k] = P(n = k)
};

namespace detail {

/**
 * Bidiagonal CoMoM recursion for a single queueing station with rate lattice
 * mu(1..Nt) and think times Z, over the sanitized classes.
 *
 * @param L  (R) per-class demand at the queueing station
 * @param N  (R) per-class population
 * @param Z  (R) per-class think time
 * @param mu (Nt) rate lattice, mu[k-1] the rate with k jobs present
 * @return the unnormalized basis h of length Nt+1, h[j] indexed as in the
 *         reference (h[Nt] is the empty-network seed)
 */
template <class T>
std::vector<T> comomrm_bidiag(const std::vector<T>& L, const std::vector<int>& N,
                              const std::vector<T>& Z, const std::vector<T>& mu) {
    const std::size_t R = N.size();
    int Nt = 0;
    for (int v : N) Nt += v;
    if (static_cast<int>(mu.size()) < Nt)
        throw InputError("comomrm_bidiag: the rate lattice is shorter than the total population");

    const T zero = num_traits<T>::from_int(0);
    std::vector<T> h(static_cast<std::size_t>(Nt) + 1, zero);
    h[static_cast<std::size_t>(Nt)] = num_traits<T>::from_int(1);

    for (std::size_t r = 0; r < R; ++r) {
        for (int nr = 1; nr <= N[r]; ++nr) {
            // 0-based row indexing rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
            std::vector<T> hn(static_cast<std::size_t>(Nt) + 1, zero);
            const T inv = num_traits<T>::from_int(1) / num_traits<T>::from_int(nr);
            for (std::size_t k = 0; k <= static_cast<std::size_t>(Nt); ++k) {
                T acc = Z[r] * h[k];
                if (k < static_cast<std::size_t>(Nt)) {
                    const std::size_t j = static_cast<std::size_t>(Nt) - k;  // 1 .. Nt
                    if (mu[j - 1] == zero)
                        throw NumericError("comomrm_bidiag: a load-dependent rate is zero");
                    acc += L[r] * num_traits<T>::from_int(static_cast<long>(j)) / mu[j - 1] *
                           h[k + 1];
                }
                hn[k] = acc * inv;
            }
            h.swap(hn);
        }
    }
    return h;
}

/** Assemble G and the marginal from the unscaled basis and the sanitize factor. */
template <class T>
ComomRmResult<T> comomrm_finish(const std::vector<T>& h, const T& Gremaind) {
    const T zero = num_traits<T>::from_int(0);
    T s = zero;
    for (const T& x : h) s += x;
    ComomRmResult<T> res;
    res.G = Gremaind * s;
    res.lG = num_traits<T>::log_as_double(res.G);
    res.prob.assign(h.size(), zero);
    if (s != zero)
        for (std::size_t k = 0; k < h.size(); ++k) res.prob[k] = h[h.size() - 1 - k] / s;
    return res;
}

}  // namespace detail

/**
 * @param L (1 x R) demands at the single queueing station
 * @param N (R) populations
 * @param Z (1 x R) think times
 * @param m replication factor of the queueing station
 * @param S number of servers per replica
 */
template <class T>
ComomRmResult<T> pfqn_comomrm_ms(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                                 int m, int S) {
    if (!L.empty() && L.rows() != 1)
        throw InputError("pfqn_comomrm_ms: the solver accepts at most a single queueing station");
    if (m < 1) throw InputError("pfqn_comomrm_ms: the replication factor must be at least one");
    if (S < 1) throw InputError("pfqn_comomrm_ms: the server count must be at least one");

    const NcSanitizeResult<T> san = pfqn_nc_sanitize(L, N, Z);
    int Nt = 0;
    for (int v : san.N) Nt += v;
    if (Nt == 0) {
        ComomRmResult<T> res;
        res.G = san.Gremaind;
        res.lG = san.lGremaind;
        res.prob.assign(1, num_traits<T>::from_int(1));
        return res;
    }

    std::vector<T> mu;
    if (m > 1) {
        mu = pfqn_mu_ms<T>(Nt, m, S);
    } else {
        mu.assign(static_cast<std::size_t>(Nt), num_traits<T>::from_int(1));
        for (int k = 1; k <= Nt; ++k)
            mu[static_cast<std::size_t>(k - 1)] = num_traits<T>::from_int(k < S ? k : S);
    }

    const std::size_t Rk = san.N.size();
    std::vector<T> Lv(Rk, num_traits<T>::from_int(0)), Zv(Rk, num_traits<T>::from_int(0));
    for (std::size_t r = 0; r < Rk; ++r) {
        if (!san.L.empty()) Lv[r] = san.L(0, r);
        for (std::size_t k = 0; k < san.Z.rows(); ++k) Zv[r] += san.Z(k, r);
    }
    return detail::comomrm_finish(detail::comomrm_bidiag(Lv, san.N, Zv, mu), san.Gremaind);
}

/** Overload with the single-replica default. */
template <class T>
ComomRmResult<T> pfqn_comomrm_ms(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                                 int S) {
    return pfqn_comomrm_ms(L, N, Z, 1, S);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_COMOMRM_MS_H
