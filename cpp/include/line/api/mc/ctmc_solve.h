/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_SOLVE_H
#define LINE_API_MC_CTMC_SOLVE_H

/**
 * Steady-state distribution of a continuous-time Markov chain.
 *
 * Templated port of matlab/lib/kpctoolbox/mc/ctmc_solve.m (the numeric branch)
 * and jar/src/main/java/jline/api/mc/Ctmc_solve.java. Solves pi Q = 0 with
 * sum(pi) = 1 by replacing the last column of the generator with ones and
 * solving Q' x = e_n, exactly as MATLAB's `Qnnz' \ bnnz` does.
 *
 * Every step is a field operation, so the exact instantiation returns pi as a
 * vector of rationals: for a generator with rational rates that is the true
 * stationary distribution with no rounding whatsoever, which is what makes the
 * exact backend worth its cost here (the double path of the JAR needs 50.9 s
 * at n=800 where this LU needs tens of ms, see _kb/14-cpp-multiprecision).
 *
 * Reducible generators are handled as MATLAB does: the weakly connected
 * components are solved separately and renormalized. States that are ISOLATED
 * (no transition in and none out) are trimmed first -- an absorbing state is NOT
 * isolated and is kept, being where the stationary mass ends up -- and a
 * generator that trims to nothing is an error rather than a plausible-looking
 * uniform vector.
 *
 * ABOVE GMRES_MIN_STATES THE LU IS ABANDONED FOR THE KRYLOV PATH, as in the
 * other three codebases: restarted GMRES first, BiCGSTAB if it reports a nonzero
 * flag, and the direct solve only if both do. The gate is compile-time as well
 * as size-based -- an exact instantiation always takes the LU, because the
 * iteration stops on a residual tolerance and normalizes by a Euclidean norm, so
 * there is nothing exact for it to converge to.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_bicgstab.h"
#include "line/api/mc/ctmc_gmres.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * Set the diagonal so that every row sums to zero (ctmc_makeinfgen).
 * Any pre-existing diagonal entry is discarded, as in MATLAB.
 */
template <class T>
Matrix<T> ctmc_makeinfgen(const Matrix<T>& Q) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_makeinfgen: generator is not square");
    Matrix<T> R = Q;
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        R(i, i) = zero;
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += R(i, j);
        R(i, i) = -s;
    }
    return R;
}

namespace detail {

/** Weakly connected components of the graph of |Q + Q'| > 0. */
template <class T>
std::vector<std::vector<std::size_t>> weak_components(const Matrix<T>& Q) {
    const std::size_t n = Q.rows();
    const T zero = num_traits<T>::from_int(0);
    std::vector<int> comp(n, -1);
    std::vector<std::vector<std::size_t>> out;
    for (std::size_t s = 0; s < n; ++s) {
        if (comp[s] >= 0) continue;
        std::vector<std::size_t> stack{s}, members;
        comp[s] = static_cast<int>(out.size());
        while (!stack.empty()) {
            const std::size_t u = stack.back();
            stack.pop_back();
            members.push_back(u);
            for (std::size_t v = 0; v < n; ++v) {
                if (v == u || comp[v] >= 0) continue;
                if (Q(u, v) != zero || Q(v, u) != zero) {
                    comp[v] = comp[u];
                    stack.push_back(v);
                }
            }
        }
        std::sort(members.begin(), members.end());
        out.push_back(members);
    }
    return out;
}

/** Submatrix on the given index set, re-normalized as a generator. */
template <class T>
Matrix<T> submatrix(const Matrix<T>& Q, const std::vector<std::size_t>& idx) {
    Matrix<T> S(idx.size(), idx.size());
    for (std::size_t a = 0; a < idx.size(); ++a)
        for (std::size_t b = 0; b < idx.size(); ++b) S(a, b) = Q(idx[a], idx[b]);
    return S;
}

}  // namespace detail

/**
 * @param Qin generator; the diagonal is recomputed, so an off-diagonal rate
 *            matrix is accepted directly
 * @return stationary distribution as a row vector of length n, summing to one
 */
template <class T>
std::vector<T> ctmc_solve(const Matrix<T>& Qin) {
    const std::size_t n = Qin.rows();
    if (Qin.cols() != n) throw InputError("ctmc_solve: generator is not square");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    if (n == 0) throw InputError("ctmc_solve: empty generator");
    if (n == 1) return std::vector<T>{one};

    const Matrix<T> Q = ctmc_makeinfgen(Qin);

    bool allZero = true;
    for (std::size_t i = 0; i < n && allZero; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (Q(i, j) != zero) {
                allZero = false;
                break;
            }
    if (allZero) {
        // No transitions at all: every distribution satisfies pi Q = 0, so the
        // stationary distribution is not unique and uniform is as good as any.
        return std::vector<T>(n, one / num_traits<T>::from_int(static_cast<long>(n)));
    }

    // Reducible: solve each weakly connected component and renormalize.
    const std::vector<std::vector<std::size_t>> comps = detail::weak_components(Q);
    if (comps.size() > 1) {
        std::vector<T> pi(n, zero);
        for (const std::vector<std::size_t>& c : comps) {
            const std::vector<T> pc = ctmc_solve(ctmc_makeinfgen(detail::submatrix(Q, c)));
            for (std::size_t k = 0; k < c.size(); ++k) pi[c[k]] = pc[k];
        }
        T s = zero;
        for (const T& v : pi) s += v;
        if (s == zero) throw NumericError("ctmc_solve: components sum to zero");
        for (T& v : pi) v /= s;
        return pi;
    }

    // Trim ISOLATED states -- no transition in and none out -- repeatedly.
    //
    // AN ABSORBING STATE IS NOT ISOLATED AND IS KEPT. Its column carries the flow that
    // reaches it, and it is where the stationary mass ends up; the trim used to also
    // require a nonzero ROW, which dropped exactly that state, left its feeders with
    // nothing to flow into, and cascaded through them until the generator was empty and
    // this function refused a chain whose distribution is unique (Q = [0 0; 1 -1] has
    // pi = [1 0]). The MATLAB and JAR twins carried the same rule and are fixed with
    // this one; see _kb/06-solver-catalog.md, where it cost SolverMAM a host-dependent
    // answer. The test here scans stored entries rather than a summed row, so it stays
    // exact and needs no tolerance -- which also keeps it correct at T = Rational.
    std::vector<std::size_t> keep(n);
    for (std::size_t i = 0; i < n; ++i) keep[i] = i;
    Matrix<T> Qk = Q;
    for (;;) {
        const std::size_t m = Qk.rows();
        std::vector<std::size_t> active;
        for (std::size_t i = 0; i < m; ++i) {
            bool colNz = false;
            for (std::size_t j = 0; j < m; ++j) {
                if (Qk(j, i) != zero) { colNz = true; break; }
            }
            if (colNz) active.push_back(i);
        }
        if (active.empty())
            throw NumericError(
                "ctmc_solve: the generator has no connected state, every state was eliminated as "
                "isolated; it admits no unique stationary distribution");
        if (active.size() == m) break;
        std::vector<std::size_t> keep2(active.size());
        for (std::size_t k = 0; k < active.size(); ++k) keep2[k] = keep[active[k]];
        keep = keep2;
        Qk = ctmc_makeinfgen(detail::submatrix(Qk, active));
    }

    // Replace the last column with ones and solve the transposed system.
    const std::size_t m = Qk.rows();
    Matrix<T> A(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) A(i, j) = (i == m - 1) ? one : Qk(j, i);
    std::vector<T> b(m, zero);
    b[m - 1] = one;

    // The iterative path. Its accuracy is a residual tolerance of 1e-12, tighter
    // than any fixed-point tolerance a caller sets, so switching to it cannot
    // move a reported metric; its failure modes are reported through the flag
    // rather than thrown, which is what makes the fallback chain possible.
    std::vector<T> x;
    if constexpr (num_traits<T>::has_transcendental) {
        if (m > GMRES_MIN_STATES) {
            const GmresResult<T> g = ctmc_gmres(A, b);
            if (g.flag == 0) {
                x = g.x;
            } else {
                const BicgstabResult<T> bs = ctmc_bicgstab(A, b);
                if (bs.flag == 0) x = bs.x;
            }
        }
    }
    if (x.empty()) x = solve(A, b);

    std::vector<T> pi(n, zero);
    for (std::size_t k = 0; k < m; ++k) pi[keep[k]] = x[k];
    return pi;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_SOLVE_H
