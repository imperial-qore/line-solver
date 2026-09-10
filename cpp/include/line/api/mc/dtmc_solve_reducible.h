/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_DTMC_SOLVE_REDUCIBLE_H
#define LINE_API_MC_DTMC_SOLVE_REDUCIBLE_H

/**
 * Limiting distribution of a discrete-time Markov chain whose transition
 * matrix may be reducible.
 *
 * Templated port of matlab/src/api/mc/dtmc_solve_reducible.m and
 * jar/src/main/java/jline/api/mc/Dtmc_solve_reducible.java. The states are
 * partitioned into strongly connected components; the components are lumped
 * into a chain Pl on which the recurrent ones are made absorbing; the limiting
 * matrix of Pl distributes the initial mass over the recurrent components; and
 * within each component the conditional limiting vector is the stationary
 * vector of the restricted chain.
 *
 * EXACT, DELIBERATELY, AND THIS IS WHERE THE PORT BEATS THE REFERENCE. MATLAB
 * computes the limiting matrix of the lumped chain by spectral decomposition
 * (spectd), with a power iteration capped at 1000 steps and tolerance 1e-10 as
 * the fallback whenever the eigenvector matrix has condition number above
 * 1e10 -- which is exactly the situation a lumped chain with repeated unit
 * eigenvalues produces, so the fallback is the common case rather than the
 * rare one, and it converges only linearly in the subdominant eigenvalue.
 * Here the same matrix is obtained in closed form: the recurrent lumped states
 * are absorbing, so the limit is the absorption probability
 *   PI(t, r) = [(I - Pl_TT)^-1 Pl_TR](t, r),  PI(r, r') = delta,
 * with I - Pl_TT non-singular because every transient component reaches a
 * recurrent one. That is one finite linear solve, no iteration and no
 * tolerance, so the whole routine instantiates at Rational and returns the
 * true limiting distribution of a rational chain.
 *
 * The only tolerance left is the one MATLAB uses to decide that a state has no
 * incoming mass (column sum below 1e-12) when no initial vector is supplied;
 * it is a structural test on the input, not a convergence criterion, and it is
 * exposed as a parameter.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/dtmc_makestochastic.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/api/mc/stronglyconncomp.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct ReducibleResult {
    std::vector<T> pi;               ///< limiting distribution, length N
    Matrix<T> pis;                   ///< numSCC x N, limiting vector per starting component
    Matrix<T> pi0;                   ///< numSCC x numSCC, the lumped starting vectors (empty if irreducible)
    std::vector<std::size_t> scc;    ///< component index of each state, 1-based
    std::vector<bool> isrec;         ///< recurrence flag per component
    Matrix<T> Pl;                    ///< lumped chain
    Matrix<T> pil;                   ///< numSCC x numSCC, limiting vector of the lumped chain
};

namespace detail {

/**
 * Limiting matrix lim_k Pl^k of a lumped chain whose recurrent states are
 * absorbing, in closed form. See the header note: this replaces MATLAB's
 * spectd / power-iteration pair and is exact.
 */
template <class T>
Matrix<T> lumped_limiting_matrix(const Matrix<T>& Pl, const std::vector<bool>& isrec) {
    const std::size_t m = Pl.rows();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<std::size_t> tr, rec;
    for (std::size_t i = 0; i < m; ++i) (isrec[i] ? rec : tr).push_back(i);

    Matrix<T> PI(m, m, zero);
    for (std::size_t i : rec) PI(i, i) = one;
    if (tr.empty()) return PI;
    if (rec.empty())
        throw NumericError(
            "dtmc_solve_reducible: the chain has no recurrent component, so no limiting "
            "distribution exists");

    // (I - Pl_TT) X = Pl_TR, solved with one factorization for all columns.
    Matrix<T> A(tr.size(), tr.size(), zero);
    for (std::size_t a = 0; a < tr.size(); ++a)
        for (std::size_t b = 0; b < tr.size(); ++b)
            A(a, b) = (a == b ? one : zero) - Pl(tr[a], tr[b]);
    Matrix<T> LU = A;
    const std::vector<std::size_t> piv = lu_factor(LU);
    for (std::size_t c = 0; c < rec.size(); ++c) {
        std::vector<T> rhs(tr.size());
        for (std::size_t a = 0; a < tr.size(); ++a) rhs[a] = Pl(tr[a], rec[c]);
        lu_solve(LU, piv, rhs);
        for (std::size_t a = 0; a < tr.size(); ++a) PI(tr[a], rec[c]) = rhs[a];
    }
    return PI;
}

}  // namespace detail

/**
 * @param P   transition matrix, possibly reducible
 * @param pin initial distribution; empty to let the routine pick one
 * @param zeroColTol column-sum threshold below which a state is treated as
 *                   having no incoming mass (MATLAB 1e-12)
 */
template <class T>
ReducibleResult<T> dtmc_solve_reducible(const Matrix<T>& P, const std::vector<T>& pin,
                                        double zeroColTol = 1e-12) {
    const std::size_t N = P.rows();
    if (P.cols() != N) throw InputError("dtmc_solve_reducible: transition matrix is not square");
    if (!pin.empty() && pin.size() != N)
        throw InputError("dtmc_solve_reducible: initial vector has the wrong length");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    const SccResult s = stronglyconncomp(P);
    const std::size_t numSCC = s.numSCC();

    ReducibleResult<T> r;
    r.scc = s.scc;
    r.isrec = s.recurrent;

    if (numSCC == 1) {
        r.pi = dtmc_solve(P);
        r.pis = Matrix<T>(1, N);
        for (std::size_t j = 0; j < N; ++j) r.pis(0, j) = r.pi[j];
        r.pi0 = Matrix<T>();
        r.Pl = P;
        r.pil = r.pis;
        return r;
    }

    // Lumped chain: mass flowing between distinct components, row-normalized,
    // then recurrent components made absorbing.
    Matrix<T> Pl(numSCC, numSCC, zero);
    for (std::size_t i = 0; i < numSCC; ++i)
        for (std::size_t j = 0; j < numSCC; ++j) {
            if (i == j) continue;
            T acc = zero;
            for (std::size_t a : s.members[i])
                for (std::size_t b : s.members[j]) acc += P(a, b);
            Pl(i, j) = acc;
        }
    Pl = dtmc_makestochastic(Pl);
    for (std::size_t i = 0; i < numSCC; ++i)
        if (s.recurrent[i]) {
            for (std::size_t j = 0; j < numSCC; ++j) Pl(i, j) = zero;
            Pl(i, i) = one;
        }
    r.Pl = Pl;

    // Probability of starting in each component.
    std::vector<T> pinl(numSCC, zero);
    if (pin.empty()) {
        const T tol = num_traits<T>::from_double(zeroColTol);
        for (std::size_t i = 0; i < numSCC; ++i) pinl[i] = one;
        for (std::size_t j = 0; j < N; ++j) {
            T cs = zero;
            for (std::size_t i = 0; i < N; ++i) cs += P(i, j);
            if (cs < tol) pinl[s.scc[j] - 1] = zero;
        }
        T tot = zero;
        for (const T& v : pinl) tot += v;
        if (tot == zero) {
            // empty-component uniform-weighting rationale: see _kb/03-api-layer.md (cpp port notes: mc)
            for (std::size_t i = 0; i < numSCC; ++i)
                pinl[i] = one / num_traits<T>::from_int(static_cast<long>(numSCC));
        } else {
            for (T& v : pinl) v /= tot;
        }
    } else {
        for (std::size_t i = 0; i < numSCC; ++i) {
            T acc = zero;
            for (std::size_t a : s.members[i]) acc += pin[a];
            pinl[i] = acc;
        }
    }

    const Matrix<T> PI = detail::lumped_limiting_matrix(Pl, s.recurrent);

    // Conditional limiting vector inside each component, computed once. It is
    // computed for EVERY component, not only the ones some starting component
    // reaches with positive weight: the `pis` rows below are addressed BY SCC
    // INDEX by the single-transient-component branch at the end, so a row left
    // unfilled is not an absent row, it is a row of zeros masquerading as a
    // distribution. `dtmc_solve_reducible.m:153-158` carries the same note.
    std::vector<std::vector<T>> within(numSCC);
    for (std::size_t j = 0; j < numSCC; ++j)
        within[j] = dtmc_solve(detail::submatrix(P, s.members[j]));

    r.pi0 = Matrix<T>(numSCC, numSCC, zero);
    r.pil = Matrix<T>(numSCC, numSCC, zero);
    r.pis = Matrix<T>(numSCC, N, zero);
    r.pi.assign(N, zero);
    for (std::size_t i = 0; i < numSCC; ++i) {
        r.pi0(i, i) = one;
        for (std::size_t j = 0; j < numSCC; ++j) r.pil(i, j) = PI(i, j);
        for (std::size_t j = 0; j < numSCC; ++j) {
            if (r.pil(i, j) == zero) continue;
            for (std::size_t k = 0; k < s.members[j].size(); ++k)
                r.pis(i, s.members[j][k]) = r.pil(i, j) * within[j][k];
        }
        // Only a component that CAN be started in enters the mixture; every
        // row is filled above regardless, for the reason given there.
        if (!(pinl[i] > zero)) continue;
        for (std::size_t k = 0; k < N; ++k) r.pi[k] += r.pis(i, k) * pinl[i];
    }

    // A single transient component and no explicit start: that component IS the
    // starting state, so its row is the answer rather than the weighted mean.
    std::size_t nTrans = 0, transIdx = 0;
    for (std::size_t i = 0; i < numSCC; ++i)
        if (!s.recurrent[i]) {
            ++nTrans;
            transIdx = i;
        }
    if (nTrans == 1 && pin.empty())
        for (std::size_t k = 0; k < N; ++k) r.pi[k] = r.pis(transIdx, k);

    return r;
}

/** Overload without an initial vector. */
template <class T>
ReducibleResult<T> dtmc_solve_reducible(const Matrix<T>& P, double zeroColTol = 1e-12) {
    return dtmc_solve_reducible(P, std::vector<T>(), zeroColTol);
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_DTMC_SOLVE_REDUCIBLE_H
