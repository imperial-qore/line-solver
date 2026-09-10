/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_LRUM_MAP_LEVELSTATS_H
#define LINE_API_CACHE_LRUM_MAP_LEVELSTATS_H

/**
 * Level statistics of one item's embedded (list, phase) chain in the
 * LRU(m)-MAP TTL approximation.
 *
 * Templated port of matlab/src/api/cache/cache_lrum_map_levelstats.m,
 * cross-checked against the levelStats method of
 * jar/src/main/java/jline/api/cache/Cache_t_lrum_map.java. Evaluates eqs.
 * (5)-(9) of Gast and Van Houdt, Performance Evaluation 2017, for an item
 * whose request process is the MAP (D0, D1) and whose lists have
 * characteristic times T(1..h):
 *
 *   E_l    = exp(D0 T_l),        N_l = (I - E_l)(-D0)^-1,     A_l = N_l D1,
 *   A_0    = (-D0)^-1 D1,        N_0 = (-D0)^-1,
 *   R_h    = A_{h-1} (I - A_h)^-1,
 *   R_l    = A_{l-1} (I - R_{l+1} E_{l+1})^-1,   l < h  (A_0 for l = 1),
 *   pi_0   = pi_0 R_1 E_1,       pi_l = pi_{l-1} R_l,
 *
 * with the level-0 boundary vector pi_0 normalized to sum one. The reported
 * probabilities weight each level by its mean holding time pi_l N_l e, and the
 * hit fractions divide the hit throughput of each list by the item's
 * stationary request rate pi D1 e.
 *
 * PERRON VECTOR WITHOUT AN EIGENDECOMPOSITION: MATLAB computes pi_0 with
 * eig(M') and picks the eigenvector of largest real part; the JAR runs power
 * iteration. Neither is needed. M = R_1 E_1 has Perron root exactly one -- it
 * is the transition matrix of the embedded chain of returns to level 0, so it
 * is stochastic in the (list, phase) sense -- and this was verified in MATLAB
 * on the reference instances, where max(real(eig(M'))) = 1 to 1e-15 at
 * characteristic times both at and far from the capacity fixed point. So pi_0
 * is the stationary vector of a stochastic matrix, i.e. dtmc_solve(M), which
 * the port already solves EXACTLY by LU. The assumption is not taken on faith:
 * the residual ||pi_0 (M - I)||_inf is checked and a NumericError is raised if
 * the Perron root is not one, rather than silently returning the stationary
 * vector of a matrix that does not have one.
 *
 * ARITHMETIC: exp(D0 T_l) is a tolerance-controlled approximation, so this
 * requires transcendental arithmetic. Everything else (the inverses, pi_0, the
 * stationary phase vector) is a linear solve and is exact.
 *
 * DIVERGENCE, MATLAB vs JAR: MATLAB uses one MAP (D0, D1) per item, shared by
 * every list. The JAR signature takes a MatrixCell per item and reads a
 * DIFFERENT (D0, D1) per list, D0c.get(l) and D1c.get(l), which is a strictly
 * larger model and is not what cache_ttl_lrum_map.m builds or what the paper
 * states. This port follows MATLAB.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Per-item level statistics, MATLAB's [prob, occ, hitfrac]. */
template <class T>
struct CacheLrumMapLevelStats {
    std::vector<T> prob;     ///< (h+1) time-stationary probability of level 0..h
    std::vector<T> occ;      ///< (h) occupancy of lists 1..h, prob(2:end)
    std::vector<T> hitfrac;  ///< (h) fraction of the item's requests hitting in list l
};

/**
 * @param D0 (d x d) hidden-transition matrix of the item's request MAP
 * @param D1 (d x d) arrival matrix of the item's request MAP
 * @param Tv (h) characteristic times, all positive
 */
template <class T>
CacheLrumMapLevelStats<T> cache_lrum_map_levelstats(const Matrix<T>& D0, const Matrix<T>& D1,
                                                    const std::vector<T>& Tv) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_lrum_map_levelstats requires transcendental arithmetic");
    const std::size_t d = D0.rows();
    const std::size_t h = Tv.size();
    if (d == 0 || D0.cols() != d || D1.rows() != d || D1.cols() != d)
        throw InputError("cache_lrum_map_levelstats: D0 and D1 must be square and of equal size");
    if (h == 0) throw InputError("cache_lrum_map_levelstats: no characteristic times");
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t l = 0; l < h; ++l)
        if (!(Tv[l] > zero))
            throw InputError("cache_lrum_map_levelstats: characteristic times must be positive");

    Matrix<T> negD0 = D0;
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t j = 0; j < d; ++j) negD0(i, j) = -D0(i, j);
    const Matrix<T> iD0 = inverse(negD0);
    const Matrix<T> I = eye<T>(d);

    std::vector<Matrix<T>> E(h), N(h), A(h);
    for (std::size_t l = 0; l < h; ++l) {
        E[l] = expm(D0, Tv[l]);
        Matrix<T> IE(d, d);
        for (std::size_t i = 0; i < d; ++i)
            for (std::size_t j = 0; j < d; ++j) IE(i, j) = I(i, j) - E[l](i, j);
        N[l] = matmul(IE, iD0);
        A[l] = matmul(N[l], D1);
    }
    const Matrix<T> A0 = matmul(iD0, D1);
    const Matrix<T>& N0 = iD0;

    // R recursion, eqs. (6)-(7); R[l] is the paper's R_{l+1}.
    std::vector<Matrix<T>> R(h);
    for (std::size_t li = h; li-- > 0;) {
        Matrix<T> lhs;   // the numerator, A_{l-1} or A_0
        Matrix<T> inner(d, d);
        if (li + 1 == h) {
            lhs = (h == 1) ? A0 : A[li - 1];
            for (std::size_t i = 0; i < d; ++i)
                for (std::size_t j = 0; j < d; ++j) inner(i, j) = I(i, j) - A[li](i, j);
        } else {
            lhs = (li == 0) ? A0 : A[li - 1];
            const Matrix<T> RE = matmul(R[li + 1], E[li + 1]);
            for (std::size_t i = 0; i < d; ++i)
                for (std::size_t j = 0; j < d; ++j) inner(i, j) = I(i, j) - RE(i, j);
        }
        R[li] = matmul(lhs, inverse(inner));
    }

    // pi_0: stationary vector of the stochastic matrix M = R_1 E_1.
    const Matrix<T> M = matmul(R[0], E[0]);
    const std::vector<T> pi0 = mc::dtmc_solve(M);
    {
        const std::vector<T> pM = vecmul(pi0, M);
        T res = zero;
        for (std::size_t i = 0; i < d; ++i) {
            const T e = num_abs(T(pM[i] - pi0[i]));
            if (e > res) res = e;
        }
        if (num_traits<T>::to_double(res) > 1e-8)
            throw NumericError(
                "cache_lrum_map_levelstats: R_1 exp(D0 T_1) has no unit Perron root, so the "
                "level-0 balance equation has no stochastic solution");
    }

    std::vector<std::vector<T>> pih(h);
    pih[0] = vecmul(pi0, R[0]);
    for (std::size_t l = 1; l < h; ++l) pih[l] = vecmul(pih[l - 1], R[l]);

    std::vector<T> holding(h + 1, zero);
    {
        const std::vector<T> v = vecmul(pi0, N0);
        for (std::size_t i = 0; i < d; ++i) holding[0] += v[i];
    }
    for (std::size_t l = 0; l < h; ++l) {
        const std::vector<T> v = vecmul(pih[l], N[l]);
        for (std::size_t i = 0; i < d; ++i) holding[1 + l] += v[i];
    }
    T denom = zero;
    for (std::size_t l = 0; l <= h; ++l) denom += holding[l];
    if (denom == zero) throw NumericError("cache_lrum_map_levelstats: zero total holding time");

    CacheLrumMapLevelStats<T> out;
    out.prob.resize(h + 1);
    for (std::size_t l = 0; l <= h; ++l) out.prob[l] = holding[l] / denom;
    out.occ.assign(out.prob.begin() + 1, out.prob.end());

    // Request-weighted hit fractions: hit throughput of list l over the item's
    // stationary request rate.
    const mam::Map<T> item{D0, D1};
    const T lam = mam::map_lambda(item);
    out.hitfrac.assign(h, zero);
    if (lam > zero) {
        for (std::size_t l = 0; l < h; ++l) {
            const std::vector<T> v = vecmul(vecmul(pih[l], N[l]), D1);
            T s = zero;
            for (std::size_t i = 0; i < d; ++i) s += v[i];
            out.hitfrac[l] = s / denom / lam;
        }
    }
    return out;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_LRUM_MAP_LEVELSTATS_H
