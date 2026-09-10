/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_DTMC_SOLVE_H
#define LINE_API_MC_DTMC_SOLVE_H

/**
 * Equilibrium distribution of a discrete-time Markov chain, and stochastic
 * complementation.
 *
 * dtmc_solve is the port of matlab/lib/kpctoolbox/mc/dtmc_solve.m: the
 * stationary vector of P is the stationary vector of the generator P - I, so
 * the whole implementation delegates to ctmc_solve. Keeping that delegation
 * literal matters for parity, since every reducibility and trimming rule then
 * lives in exactly one place.
 *
 * ctmc_stochcomp is the port of matlab/src/api/mc/ctmc_stochcomp.m: the
 * stochastic complement of the state subset I,
 *   S = Q11 + Q12 (-Q22)^-1 Q21,
 * which is itself a generator on I with the same stationary distribution up to
 * renormalization. The MATLAB version switches to GMRES above 6000 states;
 * there is no iterative path here yet, and none is needed for the exact
 * arithmetic, where an iterative method has no meaning.
 */

#include <cstddef>
#include <cstring>
#include <mutex>
#include <type_traits>
#include <utility>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

namespace detail {

/**
 * Exact sameness of two matrix entries, for the dtmc_solve memo below.
 *
 * For a hardware float this is a BIT comparison, not `==`: `==` says -0.0 equals
 * +0.0 and says NaN equals nothing, and the memo must not conflate two inputs
 * that could produce different output bits, nor be clever about NaN. Bits make
 * both cases a miss, which is always safe. Other arithmetics (the multiprecision
 * reals, the exact rationals) have no signed zero and no NaN, so value equality
 * is exact sameness there.
 */
template <class T>
inline typename std::enable_if<std::is_floating_point<T>::value, bool>::type dtmc_same(
    const T& a, const T& b) {
    return std::memcmp(&a, &b, sizeof(T)) == 0;
}

template <class T>
inline typename std::enable_if<!std::is_floating_point<T>::value, bool>::type dtmc_same(
    const T& a, const T& b) {
    return a == b;
}

template <class T>
bool dtmc_same_matrix(const Matrix<T>& a, const Matrix<T>& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j)
            if (!dtmc_same(a(i, j), b(i, j))) return false;
    return true;
}

/**
 * Bounded memo, one instantiation per arithmetic type.
 *
 * The layered fixed point asks dtmc_solve the same question over and over: a
 * layer's routing does not change between SolverLN iterations, only its rates
 * do, and visits do not depend on rates, so the per-chain body of
 * sn_refresh_visits (NetworkStruct::refresh_chains) re-solves one identical
 * chain per layer per iteration. Native Python measured 2410 calls carrying ONE
 * distinct matrix on lqn_ofbiz. The twins are Dtmc_solve.java,
 * python/line_solver/api/mc/dtmc.py and matlab/lib/kpctoolbox/mc/dtmc_solve.m;
 * keep the four in step.
 */
template <class T>
struct DtmcSolveMemo {
    static constexpr std::size_t kMax = 32;
    std::vector<std::pair<Matrix<T>, std::vector<T>>> entries;
    std::mutex mu;

    static DtmcSolveMemo& instance() {
        static DtmcSolveMemo m;
        return m;
    }
};

}  // namespace detail

/** Stationary distribution of a stochastic matrix P. */
template <class T>
std::vector<T> dtmc_solve(const Matrix<T>& P) {
    const std::size_t n = P.rows();
    if (P.cols() != n) throw InputError("dtmc_solve: transition matrix is not square");

    detail::DtmcSolveMemo<T>& memo = detail::DtmcSolveMemo<T>::instance();
    {
        std::lock_guard<std::mutex> lk(memo.mu);
        for (std::size_t k = memo.entries.size(); k-- > 0;) {
            if (detail::dtmc_same_matrix(memo.entries[k].first, P)) {
                // Returned BY VALUE, so a caller that writes into the result cannot
                // reach the stored copy.
                return memo.entries[k].second;
            }
        }
    }

    Matrix<T> Q = P;
    const T one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < n; ++i) Q(i, i) -= one;
    std::vector<T> pi = ctmc_solve(Q);

    {
        std::lock_guard<std::mutex> lk(memo.mu);
        memo.entries.emplace_back(P, pi);
        if (memo.entries.size() > detail::DtmcSolveMemo<T>::kMax) memo.entries.erase(memo.entries.begin());
    }
    return pi;
}

template <class T>
struct StochCompResult {
    Matrix<T> S;    ///< stochastic complement on the selected states
    Matrix<T> Q11;  ///< the four blocks, as MATLAB returns them
    Matrix<T> Q12;
    Matrix<T> Q21;
    Matrix<T> Q22;
    Matrix<T> T12;  ///< Q12 (-Q22)^-1 Q21, the correction term
};

/**
 * @param Q generator
 * @param I state subset to keep; defaults to the first ceil(n/2) states
 */
template <class T>
StochCompResult<T> ctmc_stochcomp(const Matrix<T>& Q, const std::vector<std::size_t>& I) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_stochcomp: generator is not square");

    std::vector<bool> selected(n, false);
    for (std::size_t k : I) {
        if (k >= n) throw InputError("ctmc_stochcomp: state index out of range");
        selected[k] = true;
    }
    std::vector<std::size_t> Ic;
    for (std::size_t i = 0; i < n; ++i)
        if (!selected[i]) Ic.push_back(i);
    if (I.empty()) throw InputError("ctmc_stochcomp: empty state subset");

    StochCompResult<T> r;
    r.Q11 = detail::submatrix(Q, I);
    r.Q22 = detail::submatrix(Q, Ic);
    r.Q12 = Matrix<T>(I.size(), Ic.size());
    r.Q21 = Matrix<T>(Ic.size(), I.size());
    for (std::size_t a = 0; a < I.size(); ++a)
        for (std::size_t b = 0; b < Ic.size(); ++b) r.Q12(a, b) = Q(I[a], Ic[b]);
    for (std::size_t a = 0; a < Ic.size(); ++a)
        for (std::size_t b = 0; b < I.size(); ++b) r.Q21(a, b) = Q(Ic[a], I[b]);

    if (Ic.empty()) {
        r.T12 = Matrix<T>(I.size(), I.size(), num_traits<T>::from_int(0));
        r.S = r.Q11;
        return r;
    }

    // T = (-Q22) \ Q21, solved column by column with one shared factorization.
    Matrix<T> A = r.Q22;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) A(i, j) = -A(i, j);
    Matrix<T> LU = A;
    const std::vector<std::size_t> piv = lu_factor(LU);

    Matrix<T> Tm(Ic.size(), I.size());
    for (std::size_t c = 0; c < I.size(); ++c) {
        std::vector<T> rhs(Ic.size());
        for (std::size_t i = 0; i < Ic.size(); ++i) rhs[i] = r.Q21(i, c);
        lu_solve(LU, piv, rhs);
        for (std::size_t i = 0; i < Ic.size(); ++i) Tm(i, c) = rhs[i];
    }

    // T = Q12 * T; S = Q11 + T.
    r.T12 = Matrix<T>(I.size(), I.size(), num_traits<T>::from_int(0));
    for (std::size_t a = 0; a < I.size(); ++a)
        for (std::size_t c = 0; c < I.size(); ++c) {
            T s = num_traits<T>::from_int(0);
            for (std::size_t k = 0; k < Ic.size(); ++k) s += r.Q12(a, k) * Tm(k, c);
            r.T12(a, c) = s;
        }
    r.S = Matrix<T>(I.size(), I.size());
    for (std::size_t a = 0; a < I.size(); ++a)
        for (std::size_t c = 0; c < I.size(); ++c) r.S(a, c) = r.Q11(a, c) + r.T12(a, c);
    return r;
}

/** Default subset: the first ceil(n/2) states, as in MATLAB. */
template <class T>
StochCompResult<T> ctmc_stochcomp(const Matrix<T>& Q) {
    const std::size_t half = (Q.rows() + 1) / 2;
    std::vector<std::size_t> I(half);
    for (std::size_t i = 0; i < half; ++i) I[i] = i;
    return ctmc_stochcomp(Q, I);
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_DTMC_SOLVE_H
