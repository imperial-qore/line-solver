/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_DTMC_TRANSIENT_H
#define LINE_API_MC_DTMC_TRANSIENT_H

/**
 * Discrete-time transient distributions, hitting times and uniformization.
 *
 * Templated port of matlab/src/api/mc/dtmc_transient.m,
 * matlab/src/api/mc/dtmc_hitting_time.m and
 * matlab/lib/kpctoolbox/mc/dtmc_uniformization.m.
 *
 * dtmc_transient returns the WHOLE trajectory, steps+1 rows with row 0 the
 * initial law, not just the distribution at the last step.
 *
 * dtmc_hitting_time solves (I - P_TT) h = 1 on the non-target states only. A
 * state that cannot reach the target set makes that system singular, and the
 * right answer there is an INFINITE hitting time, not the least-squares
 * solution of the singular system: the code detects the singularity and reports
 * infinity for the unreachable block rather than a finite fabricated number.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_uniformization.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** Trajectory of the law over steps transitions, row k holding pi0 P^k. */
template <class T>
Matrix<T> dtmc_transient(const Matrix<T>& P, const std::vector<T>& pi0, std::size_t steps) {
    const std::size_t n = P.rows();
    if (P.cols() != n) throw InputError("dtmc_transient: P is not square");
    std::vector<T> pik = pi0;
    if (pik.empty())
        pik.assign(n, num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<int>(n)));
    if (pik.size() != n) throw InputError("dtmc_transient: pi0 has the wrong length");
    Matrix<T> out(steps + 1, n, num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < n; ++j) out(0, j) = pik[j];
    for (std::size_t k = 1; k <= steps; ++k) {
        std::vector<T> next(n, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) next[j] += pik[i] * P(i, j);
        pik = next;
        for (std::size_t j = 0; j < n; ++j) out(k, j) = pik[j];
    }
    return out;
}

/** Uniform initial law and one step, matching the one-argument MATLAB call. */
template <class T>
Matrix<T> dtmc_transient(const Matrix<T>& P) {
    return dtmc_transient(P, std::vector<T>(), 1);
}

/** Expected number of steps to reach the target set, zero on the target set itself. */
template <class T>
std::vector<T> dtmc_hitting_time(const Matrix<T>& P, const std::vector<std::size_t>& target) {
    static_assert(num_traits<T>::has_transcendental,
                  "dtmc_hitting_time requires a backend with an infinity, since a state that "
                  "cannot reach the target set has an infinite hitting time");
    const std::size_t n = P.rows();
    if (P.cols() != n) throw InputError("dtmc_hitting_time: P is not square");
    std::vector<bool> is_target(n, false);
    for (std::size_t k = 0; k < target.size(); ++k) {
        if (target[k] >= n) throw InputError("dtmc_hitting_time: target index out of range");
        is_target[target[k]] = true;
    }
    std::vector<std::size_t> nt;
    for (std::size_t i = 0; i < n; ++i)
        if (!is_target[i]) nt.push_back(i);
    std::vector<T> h(n, num_traits<T>::from_int(0));
    if (nt.empty()) return h;
    const std::size_t m = nt.size();
    Matrix<T> A(m, m, num_traits<T>::from_int(0));
    for (std::size_t a = 0; a < m; ++a)
        for (std::size_t b = 0; b < m; ++b)
            A(a, b) = (a == b ? num_traits<T>::from_int(1) : num_traits<T>::from_int(0)) -
                      P(nt[a], nt[b]);
    std::vector<T> b(m, num_traits<T>::from_int(1));
    // A state from which the target set is unreachable makes I - P_TT singular
    // there, and its hitting time is infinite rather than any finite solution.
    std::vector<bool> reaches(m, false);
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t a = 0; a < m; ++a) {
            if (reaches[a]) continue;
            for (std::size_t j = 0; j < n; ++j) {
                if (P(nt[a], j) == num_traits<T>::from_int(0)) continue;
                bool ok = is_target[j];
                if (!ok)
                    for (std::size_t c = 0; c < m; ++c)
                        if (nt[c] == j && reaches[c]) ok = true;
                if (ok) {
                    reaches[a] = true;
                    changed = true;
                    break;
                }
            }
        }
    }
    std::vector<std::size_t> keep;
    for (std::size_t a = 0; a < m; ++a) {
        if (reaches[a])
            keep.push_back(a);
        else
            h[nt[a]] = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    }
    if (keep.empty()) return h;
    Matrix<T> Ak(keep.size(), keep.size(), num_traits<T>::from_int(0));
    std::vector<T> bk(keep.size(), num_traits<T>::from_int(1));
    for (std::size_t a = 0; a < keep.size(); ++a)
        for (std::size_t c = 0; c < keep.size(); ++c) Ak(a, c) = A(keep[a], keep[c]);
    const std::vector<T> hk = solve(Ak, bk);
    for (std::size_t a = 0; a < keep.size(); ++a) h[nt[keep[a]]] = hk[a];
    return h;
}

/** Transient law of a DTMC through the uniformized generator of P. */
template <class T>
UniformizationResult<T> dtmc_uniformization(const std::vector<T>& pi0, const Matrix<T>& P,
                                            const T& t, double tol = 1e-12, long maxiter = -1) {
    return ctmc_uniformization(pi0, ctmc_makeinfgen(P), t, tol, maxiter);
}

}  // namespace mc
}  // namespace line

#endif
