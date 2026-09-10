/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_QUORUM_MOMENTS_H
#define LINE_API_FJ_QUORUM_MOMENTS_H

/**
 * Mean and variance of a k-of-n (quorum) join completion time, from the mean
 * and variance of each branch.
 *
 * Templated port of matlab/src/api/fj/fj_quorum_moments.m, cross-checked
 * against jar/src/main/java/jline/api/fj/FJ_quorum.java (identical
 * three-point fit, merged grid and Poisson-binomial recurrence).
 *
 * Each branch is expanded into a three-point discrete step CDF matching its
 * first two moments. At any time the number of completed branches is
 * Poisson-binomial, so its distribution is built by the recurrence
 * q_j <- q_{j-1} F_i + q_j (1 - F_i) over branches, and the k-th order
 * statistic is the upper tail sum_{j>=k} q_j. At k = n this reduces to the
 * product of the branch CDFs (the ordinary AND-join) and at k = 1 to
 * 1 - prod(1 - F_i) (the minimum), two identities the port must reproduce.
 *
 * The recurrence is the reason for the whole construction: it computes the
 * same quantity as the inclusion-exclusion identity used by LQNS, but every
 * term is a probability in [0,1] and none is subtracted, so it avoids the
 * catastrophic cancellation the alternating binomial sum incurs as n grows.
 *
 * static_assert(num_traits<T>::has_transcendental) -- the three-point fit
 * takes the standard deviation, so a single sqrt puts the whole function
 * outside the field. Everything downstream of that sqrt is rational.
 *
 * Follows Omari, Franks, Woodside and Pan, as implemented in LQNS 6.x
 * (randomvar.cc).
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fj {

namespace detail {

/** A step CDF as abscissae t with cumulative masses A. */
template <class T>
struct StepCdf {
    std::vector<T> t;
    std::vector<T> A;
};

/** Two-moment fit of a branch to a three-point discrete distribution. */
template <class T>
StepCdf<T> three_point_fit(const T& mu, const T& var) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (mu < zero || var < zero)
        throw InputError("fj_quorum_moments: branch mean and variance must be non-negative");
    StepCdf<T> f;
    if (mu == zero) return f;  // no mass: caller replaces it with an instant completion
    const T sd = num_sqrt(var);
    if (sd == zero) {
        f.t.push_back(mu);
        f.A.push_back(one);
        return f;
    }
    const T t1 = (mu > sd) ? T(mu - sd) : zero;
    const T t2 = mu;
    const T t3 = (sd >= mu) ? T(mu + two * var / mu) : T(mu + two * sd);
    const T delta = t1 * t1 * (t3 - t2) + t2 * t2 * (t1 - t3) + t3 * t3 * (t2 - t1);
    if (delta == zero) {  // abscissae not distinct, the fit is undetermined
        f.t.push_back(mu);
        f.A.push_back(one);
        return f;
    }
    const T temp = var + mu * mu;
    const T a1 = (temp * (t3 - t2) + t2 * t2 * (mu - t3) + t3 * t3 * (t2 - mu)) / delta;
    const T a3 = (t1 * t1 * (mu - t2) + t2 * t2 * (t1 - mu) + temp * (t2 - t1)) / delta;
    f.t = {t1, t2, t3};
    f.A = {a1, T(one - a3), one};
    return f;
}

}  // namespace detail

/** Branch counts above this are refused: evaluation is cubic in n. */
constexpr std::size_t FJ_QUORUM_MAX_BRANCHES = 512;

/**
 * @param branchMeans mean completion time of each branch
 * @param branchVars  variance of each branch, same length
 * @param k           quorum size, 1 <= k <= n
 * @return            [m, v], the mean and variance of the k-th smallest
 */
template <class T>
FJQuorumMomentsResult<T> fj_quorum_moments(const std::vector<T>& branchMeans,
                                           const std::vector<T>& branchVars, std::size_t k) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_quorum_moments requires transcendental arithmetic");
    const std::size_t n = branchMeans.size();
    if (n != branchVars.size())
        throw InputError("fj_quorum_moments: branchMeans and branchVars must have the same length");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (n == 0) return {zero, zero};
    if (k < 1 || k > n) throw InputError("fj_quorum_moments: k must satisfy 1 <= k <= n");
    if (n > FJ_QUORUM_MAX_BRANCHES)
        throw InputError(
            "fj_quorum_moments: too many branches; evaluation is cubic in the branch count");

    std::vector<detail::StepCdf<T>> branches(n);
    for (std::size_t i = 0; i < n; ++i) {
        branches[i] = detail::three_point_fit(branchMeans[i], branchVars[i]);
        if (branches[i].t.empty()) {  // a massless branch completes instantly
            branches[i].t.push_back(zero);
            branches[i].A.push_back(one);
        }
    }

    std::vector<T> grid;
    for (std::size_t i = 0; i < n; ++i) grid.insert(grid.end(), branches[i].t.begin(), branches[i].t.end());
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
    const std::size_t ngrid = grid.size();

    // Tabulate each branch along the merged grid by a single monotone walk.
    Matrix<T> fv(n, ngrid, zero);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t p = 0;
        T cur = zero;
        for (std::size_t g = 0; g < ngrid; ++g) {
            while (p < branches[i].t.size() && branches[i].t[p] <= grid[g]) {
                cur = branches[i].A[p];
                ++p;
            }
            fv(i, g) = cur;
        }
    }

    std::vector<T> A(ngrid, zero);
    for (std::size_t g = 0; g < ngrid; ++g) {
        std::vector<T> q(n + 1, zero);
        q[0] = one;
        for (std::size_t i = 0; i < n; ++i) {
            const T f = fv(i, g);
            for (std::size_t j = std::min(i + 1, n); j >= 1; --j) q[j] = q[j - 1] * f + q[j] * (one - f);
            q[0] *= (one - f);
        }
        T tail = zero;
        for (std::size_t j = k; j <= n; ++j) tail += q[j];
        A[g] = tail;
    }

    T m = zero, prev = zero;
    for (std::size_t g = 0; g < ngrid; ++g) {
        m += (A[g] - prev) * grid[g];
        prev = A[g];
    }
    T v = zero;
    prev = zero;
    for (std::size_t g = 0; g < ngrid; ++g) {
        const T d = grid[g] - m;
        v += (A[g] - prev) * d * d;
        prev = A[g];
    }
    if (v < zero) v = zero;
    return {m, v};
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_QUORUM_MOMENTS_H
