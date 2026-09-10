/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_SOLVE_REDUCIBLE_BLKDECOMP_H
#define LINE_API_MC_CTMC_SOLVE_REDUCIBLE_BLKDECOMP_H

/**
 * Limiting distribution of a reducible CTMC by direct block decomposition of
 * the generator.
 *
 * Templated port of matlab/src/api/mc/ctmc_solve_reducible_blkdecomp.m and
 * jar/src/main/java/jline/api/mc/Ctmc_solve_reducible_blkdecomp.java. Unlike
 * ctmc_solve_reducible it never uniformizes: the states are split into
 * transient and recurrent classes by strong connectivity, the expected sojourn
 * of the transient part is obtained from sojourn Q_tt = -p0_t (Q_tt is Hurwitz,
 * so this is a plain non-singular solve), the absorption probabilities follow
 * as hit = sojourn Q_ta + p0_r, and each recurrent class contributes its own
 * stationary vector scaled by the probability of reaching it.
 *
 * EXACT, AND THE EXACTNESS IS THE POINT. Every step is a finite linear solve
 * over the field of the rates: strong connectivity is combinatorial, Q_tt is
 * inverted once per starting class, and each recurrent class goes through
 * ctmc_solve. At Rational the absorption probabilities are exact rationals,
 * where the reference computes them in double precision on a matrix that is
 * ill-conditioned precisely when the transient class is nearly closed -- the
 * regime the routine exists to handle.
 *
 * The two thresholds MATLAB uses are structural tests on the input, not
 * convergence criteria, and are exposed as parameters: a reachability
 * probability below 1e-15 is treated as unreachable, and a state whose column
 * of |Q| sums below 1e-12 is treated as having no incoming rate.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/stronglyconncomp.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct BlkDecompResult {
    std::vector<T> pi;             ///< limiting distribution, length N
    Matrix<T> pis;                 ///< numSCC x N, one row per starting component
    Matrix<T> pi0;                 ///< numSCC x N, the uniform-in-component starting vectors
    std::vector<std::size_t> scc;  ///< component index of each state, 1-based
    std::vector<bool> isrec;       ///< recurrence flag per component
};

/**
 * @param Qin generator; the diagonal is recomputed
 * @param pin initial distribution; empty if not available
 * @param reachTol probability below which a recurrent class counts as unreached
 * @param zeroColTol column-sum threshold for a state with no incoming rate
 */
template <class T>
BlkDecompResult<T> ctmc_solve_reducible_blkdecomp(const Matrix<T>& Qin, const std::vector<T>& pin,
                                                  double reachTol = 1e-15,
                                                  double zeroColTol = 1e-12) {
    const std::size_t N = Qin.rows();
    if (Qin.cols() != N) throw InputError("ctmc_solve_reducible_blkdecomp: generator is not square");
    if (!pin.empty() && pin.size() != N)
        throw InputError("ctmc_solve_reducible_blkdecomp: initial vector has the wrong length");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    const Matrix<T> Q = ctmc_makeinfgen(Qin);

    // Adjacency from the off-diagonal entries, which are non-negative in a
    // generator, so "non-zero" and MATLAB's "> 0" coincide.
    Matrix<T> Adj = Q;
    for (std::size_t i = 0; i < N; ++i) Adj(i, i) = zero;
    const SccResult s = stronglyconncomp(Adj);
    const std::size_t numSCC = s.numSCC();

    BlkDecompResult<T> r;
    r.scc = s.scc;
    r.isrec = s.recurrent;

    if (numSCC == 1) {
        r.pi = ctmc_solve(Q);
        r.pis = Matrix<T>(1, N);
        for (std::size_t j = 0; j < N; ++j) r.pis(0, j) = r.pi[j];
        r.pi0 = Matrix<T>();
        return r;
    }

    std::vector<std::size_t> transStates, recStates, transSccIds, recSccIds;
    for (std::size_t c = 0; c < numSCC; ++c) {
        if (s.recurrent[c]) {
            recSccIds.push_back(c);
            recStates.insert(recStates.end(), s.members[c].begin(), s.members[c].end());
        } else {
            transSccIds.push_back(c);
            transStates.insert(transStates.end(), s.members[c].begin(), s.members[c].end());
        }
    }
    std::sort(transStates.begin(), transStates.end());
    std::sort(recStates.begin(), recStates.end());
    const std::size_t nt = transStates.size(), nr = recStates.size();
    if (nr == 0)
        throw NumericError(
            "ctmc_solve_reducible_blkdecomp: no recurrent class, the chain admits no limiting "
            "distribution");

    // Q_tt (transposed, ready for the solve) and Q_ta.
    Matrix<T> Q_ttT(nt, nt), Q_ta(nt, nr);
    std::vector<std::size_t> lupiv;
    Matrix<T> LU;
    if (nt > 0) {
        for (std::size_t a = 0; a < nt; ++a) {
            for (std::size_t b = 0; b < nt; ++b) Q_ttT(b, a) = Q(transStates[a], transStates[b]);
            for (std::size_t b = 0; b < nr; ++b) Q_ta(a, b) = Q(transStates[a], recStates[b]);
        }
        LU = Q_ttT;
        lupiv = lu_factor(LU);
    }

    const T rtol = num_traits<T>::from_double(reachTol);
    r.pis = Matrix<T>(numSCC, N, zero);
    r.pi0 = Matrix<T>(numSCC, N, zero);
    // Position of each state within recStates, for the per-class gather below.
    std::vector<std::size_t> recPos(N, static_cast<std::size_t>(-1));
    for (std::size_t k = 0; k < nr; ++k) recPos[recStates[k]] = k;

    for (std::size_t c = 0; c < numSCC; ++c) {
        std::vector<T> p0(N, zero);
        const T w = one / num_traits<T>::from_int(static_cast<long>(s.members[c].size()));
        for (std::size_t a : s.members[c]) p0[a] = w;
        for (std::size_t j = 0; j < N; ++j) r.pi0(c, j) = p0[j];

        std::vector<T> hit(nr, zero);
        if (nt > 0) {
            bool anyT = false;
            std::vector<T> rhs(nt);
            for (std::size_t a = 0; a < nt; ++a) {
                rhs[a] = -p0[transStates[a]];
                if (rhs[a] != zero) anyT = true;
            }
            if (anyT) {
                // Solve sojourn * Q_tt = -p0_t, i.e. Q_tt' * sojourn' = -p0_t'.
                std::vector<T> sojourn = rhs;
                lu_solve(LU, lupiv, sojourn);
                for (std::size_t b = 0; b < nr; ++b) {
                    T acc = zero;
                    for (std::size_t a = 0; a < nt; ++a) acc += sojourn[a] * Q_ta(a, b);
                    hit[b] = acc;
                }
            }
        }
        for (std::size_t k = 0; k < nr; ++k) hit[k] += p0[recStates[k]];

        for (std::size_t cr : recSccIds) {
            const std::vector<std::size_t>& idx = s.members[cr];
            T reach = zero;
            for (std::size_t a : idx) reach += hit[recPos[a]];
            if (reach < rtol) continue;
            if (idx.size() == 1) {
                r.pis(c, idx[0]) = reach;
            } else {
                const std::vector<T> pi_c = ctmc_solve(detail::submatrix(Q, idx));
                for (std::size_t k = 0; k < idx.size(); ++k) r.pis(c, idx[k]) = pi_c[k] * reach;
            }
        }
    }

    // Probability of starting in each component.
    std::vector<T> pinl(numSCC, zero);
    if (pin.empty()) {
        const T ztol = num_traits<T>::from_double(zeroColTol);
        for (std::size_t c = 0; c < numSCC; ++c) pinl[c] = one;
        for (std::size_t j = 0; j < N; ++j) {
            T cs = zero;
            for (std::size_t i = 0; i < N; ++i) cs += num_abs(T(Q(i, j)));
            if (cs < ztol) pinl[s.scc[j] - 1] = zero;
        }
        T tot = zero;
        for (const T& v : pinl) tot += v;
        if (tot > zero) {
            for (T& v : pinl) v /= tot;
        } else {
            for (std::size_t c = 0; c < numSCC; ++c)
                pinl[c] = one / num_traits<T>::from_int(static_cast<long>(numSCC));
        }
    } else {
        for (std::size_t c = 0; c < numSCC; ++c) {
            T acc = zero;
            for (std::size_t a : s.members[c]) acc += pin[a];
            pinl[c] = acc;
        }
    }

    r.pi.assign(N, zero);
    for (std::size_t c = 0; c < numSCC; ++c) {
        if (!(pinl[c] > zero)) continue;
        for (std::size_t j = 0; j < N; ++j) r.pi[j] += r.pis(c, j) * pinl[c];
    }
    if (transSccIds.size() == 1 && pin.empty())
        for (std::size_t j = 0; j < N; ++j) r.pi[j] = r.pis(transSccIds[0], j);

    T tot = zero;
    for (const T& v : r.pi) tot += v;
    if (tot > zero)
        for (T& v : r.pi) v /= tot;
    return r;
}

/** Overload without an initial vector. */
template <class T>
BlkDecompResult<T> ctmc_solve_reducible_blkdecomp(const Matrix<T>& Q, double reachTol = 1e-15,
                                                  double zeroColTol = 1e-12) {
    return ctmc_solve_reducible_blkdecomp(Q, std::vector<T>(), reachTol, zeroColTol);
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_SOLVE_REDUCIBLE_BLKDECOMP_H
