/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MM1_DPS_H
#define LINE_API_QSYS_QSYS_MM1_DPS_H

/**
 * Multiclass M/M/1 under DPS (discriminatory processor sharing), solved
 * numerically on the truncated population chain.
 *
 * Templated port of matlab/src/api/qsys/qsys_mm1_dps.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mm1_dps.java.
 *
 * The state is the per-class population vector (n_1,...,n_K); class k arrives
 * at rate lambda_k and completes at rate
 *
 *   mu_k n_k w_k / sum_j n_j w_j,
 *
 * so the server capacity is split in proportion to the weighted populations.
 * The chain is truncated at a total population level taken from the geometric
 * tail bound, N = max(16, ceil(log(tol)/log(rho))), and the level is doubled
 * until the per-class mean counts move by less than tol or the hard cutoff is
 * reached. Response times follow by Little's law, T_k = E[N_k]/lambda_k.
 *
 * The truncation blocks arrivals at the top level rather than dropping them,
 * which is what makes the answer a lower bound that converges from below as
 * the cutoff grows; that is the reference behaviour and is kept.
 *
 * ARITHMETIC. The cutoff is chosen through a logarithm and the doubling loop
 * is driven to a tolerance, so the function is gated on transcendental
 * arithmetic. The inner solve is exact field arithmetic given the truncation,
 * but the truncation itself is the approximation.
 *
 * COST. The truncated state space has C(N+K,K) states and the stationary
 * solve is dense here, so the MATLAB default maxCutoff of 2048 is not usable
 * at K >= 2; callers should pass a cutoff matched to the load. The defaults
 * are kept as the MATLAB ones so that a like-for-like comparison is possible
 * on small instances.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

template <class T>
struct Mm1DpsResult {
    std::vector<T> T_;  ///< per-class mean response time
    T rho;              ///< total utilization sum_k lambda_k/mu_k
};

namespace detail {

/** Recursive helper of dps_states. */
inline void dps_states_rec(std::size_t K, unsigned rem, std::vector<unsigned>& cur,
                           std::vector<std::vector<unsigned>>& out) {
    if (cur.size() == K) {
        out.push_back(cur);
        return;
    }
    for (unsigned v = 0; v <= rem; ++v) {
        cur.push_back(v);
        dps_states_rec(K, rem - v, cur, out);
        cur.pop_back();
    }
}

/** All population vectors of K classes with total at most Ncut. */
inline std::vector<std::vector<unsigned>> dps_states(std::size_t K, unsigned Ncut) {
    std::vector<std::vector<unsigned>> out;
    std::vector<unsigned> cur;
    cur.reserve(K);
    dps_states_rec(K, Ncut, cur, out);
    return out;
}

/** Per-class mean populations of the chain truncated at total level Ncut. */
template <class T>
std::vector<T> dps_solve_trunc(const std::vector<T>& lambda, const std::vector<T>& mu,
                               const std::vector<T>& w, unsigned Ncut) {
    const std::size_t K = lambda.size();
    const std::vector<std::vector<unsigned>> S = dps_states(K, Ncut);
    const std::size_t n = S.size();
    std::vector<std::size_t> stride(K, 1);
    for (std::size_t k = 1; k < K; ++k) stride[k] = stride[k - 1] * (Ncut + 1);
    std::vector<std::size_t> index(stride[K - 1] * (Ncut + 1), n);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t key = 0;
        for (std::size_t k = 0; k < K; ++k) key += S[i][k] * stride[k];
        index[key] = i;
    }

    const T zero = num_traits<T>::from_int(0);
    Matrix<T> Q(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        unsigned tot = 0;
        for (std::size_t k = 0; k < K; ++k) tot += S[i][k];
        T den = zero;
        for (std::size_t k = 0; k < K; ++k)
            den += num_traits<T>::from_int(static_cast<long>(S[i][k])) * w[k];
        std::size_t key = 0;
        for (std::size_t k = 0; k < K; ++k) key += S[i][k] * stride[k];
        for (std::size_t k = 0; k < K; ++k) {
            if (tot < Ncut) Q(i, index[key + stride[k]]) += lambda[k];
            if (S[i][k] > 0)
                Q(i, index[key - stride[k]]) +=
                    mu[k] * num_traits<T>::from_int(static_cast<long>(S[i][k])) * w[k] / den;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < n; ++j)
            if (j != i) s += Q(i, j);
        Q(i, i) = -s;
    }
    const std::vector<T> pi = mc::ctmc_solve(Q);
    std::vector<T> EN(K, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < K; ++k)
            EN[k] += pi[i] * num_traits<T>::from_int(static_cast<long>(S[i][k]));
    return EN;
}

}  // namespace detail

/**
 * @param lambda    per-class Poisson arrival rates
 * @param mu        per-class exponential service rates
 * @param w         per-class DPS weights, positive
 * @param tol       convergence tolerance on the mean counts (MATLAB 1e-10)
 * @param maxCutoff hard bound on the truncation level (MATLAB 2048)
 */
template <class T>
Mm1DpsResult<T> qsys_mm1_dps(const std::vector<T>& lambda, const std::vector<T>& mu,
                             const std::vector<T>& w, const T& tol, unsigned maxCutoff) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mm1_dps requires transcendental arithmetic");
    const std::size_t K = lambda.size();
    if (mu.size() != K || w.size() != K)
        throw InputError("qsys_mm1_dps: lambda, mu and w must have the same length");
    if (K == 0) throw InputError("qsys_mm1_dps: at least one class is required");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    for (std::size_t k = 0; k < K; ++k)
        if (lambda[k] <= zero || mu[k] <= zero || w[k] <= zero)
            throw InputError("qsys_mm1_dps: lambda, mu and w must all be positive");
    T rho = zero;
    for (std::size_t k = 0; k < K; ++k) rho += lambda[k] / mu[k];
    if (rho >= one) throw InputError("qsys_mm1_dps: system is unstable, rho >= 1");

    const double lt = num_traits<T>::log_as_double(tol);
    const double lr = num_traits<T>::log_as_double(rho);
    unsigned N = static_cast<unsigned>(std::max(16.0, std::ceil(lt / lr)));
    if (N > maxCutoff) N = maxCutoff;

    std::vector<T> ENprev = detail::dps_solve_trunc(lambda, mu, w, N);
    while (N < maxCutoff) {
        const unsigned N2 = std::min(2u * N, maxCutoff);
        const std::vector<T> EN = detail::dps_solve_trunc(lambda, mu, w, N2);
        T gap = zero;
        for (std::size_t k = 0; k < K; ++k) {
            const T d = num_abs(T(EN[k] - ENprev[k]));
            if (d > gap) gap = d;
        }
        ENprev = EN;
        if (gap < tol) break;
        N = N2;
        if (N2 == maxCutoff) break;
    }

    Mm1DpsResult<T> r;
    r.rho = rho;
    r.T_.resize(K);
    for (std::size_t k = 0; k < K; ++k) r.T_[k] = ENprev[k] / lambda[k];
    return r;
}

/** MATLAB defaults: tol = 1e-10, maxCutoff = 2048. */
template <class T>
Mm1DpsResult<T> qsys_mm1_dps(const std::vector<T>& lambda, const std::vector<T>& mu,
                             const std::vector<T>& w) {
    return qsys_mm1_dps(lambda, mu, w, T(num_traits<T>::from_double(1e-10)), 2048u);
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MM1_DPS_H
