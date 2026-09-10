/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_JOINTMARG_H
#define LINE_API_PFQN_JOINTMARG_H

/**
 * Joint probability of the per-station TOTAL queue lengths.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_jointmarg.m, of
 * jline.api.pfqn.Pfqn_jointmarg and of python pfqn_jointmarg.
 *
 * Joint probability that station i holds n(i) jobs IN TOTAL, all classes summed
 * out, in a closed multiclass product-form network:
 *
 *   P(n_1,...,n_M) = perm(A) / ( prod_r N_r! * prod_{j in INFSET} n_j! * G(N) )
 *
 * with A the demand matrix whose column r is repeated N_r times and whose row i
 * is repeated n_i times, so A is square of order sum(N).
 *
 * HOW THIS DIFFERS FROM pfqn_joint_total, which is the same identity with the
 * delay taken as ONE aggregated row: here every infinite-server station keeps
 * its own row and contributes its own 1/n_j!. The queueing stations contribute
 * the n_i! that the permanent identity supplies; the infinite servers do not.
 * Dividing once, as a single aggregated delay row would, leaves the law
 * unnormalized as soon as the model has two delays.
 *
 * The identity holds for load-independent single-server queues plus infinite
 * servers. Multiserver and load-dependent stations break the n_i! factor and
 * are the caller's responsibility to exclude (see solver_nc_jointmarg).
 *
 * ZERO ELEMENTS are safe under the exact engine and only under it: a station
 * holding no jobs contributes no row, a class with no jobs contributes no
 * column, a zero demand is an ordinary zero entry of A, and the permanent of
 * the empty matrix is 1. The approximate engines are REFUSED on a matrix with a
 * structural zero rather than having it floored at eps: Sinkhorn scaling needs
 * full support, and the Bethe gap is a state-dependent lower bound that does
 * not cancel when the estimates are normalized against each other.
 *
 * Arithmetic: EXACT-CAPABLE on the "exact" engine, which reaches pfqn_perm and
 * uses additions, multiplications and exact binomials only. The four
 * approximate engines are double-precision Monte Carlo or message passing, so
 * they collapse an exact T to double before running.
 *
 * Reference:
 *   H. J. Ryser, "Combinatorial Mathematics", Carus Mathematical Monographs 14,
 *   Mathematical Association of America, 1963.
 */

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/perm/perm_approx.h"
#include "line/api/perm/perm_sampling.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_perm.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace jointmargdetail {

/** Lowercase copy, so the engine name is matched case-insensitively. */
inline std::string lower(const std::string& s) {
    std::string out = s;
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(out[i])));
    return out;
}

/**
 * Row i of L repeated n(i) times, kept in the (rows x R) form pfqn_perm takes,
 * whose column r stands for N_r identical columns.
 *
 * A station holding no jobs drops out here, which is what makes a zero entry of
 * the occupancy vector free of any special case: the expanded matrix stays
 * square of order sum(N).
 */
template <class T>
Matrix<T> replicate_rows(const Matrix<T>& L, const std::vector<int>& n) {
    const std::size_t M = L.rows(), R = L.cols();
    std::vector<std::size_t> rowOf;
    for (std::size_t i = 0; i < M; ++i)
        for (int c = 0; c < n[i]; ++c) rowOf.push_back(i);
    Matrix<T> A(rowOf.size(), R);
    for (std::size_t a = 0; a < rowOf.size(); ++a)
        for (std::size_t r = 0; r < R; ++r) A(a, r) = L(rowOf[a], r);
    return A;
}

/** The fully expanded square matrix, which the approximate engines need. */
template <class T>
Matrix<double> expand_to_double(const Matrix<T>& A, const std::vector<int>& N) {
    std::vector<std::size_t> colOf;
    for (std::size_t r = 0; r < N.size(); ++r)
        for (int c = 0; c < N[r]; ++c) colOf.push_back(r);
    Matrix<double> out(A.rows(), colOf.size(), 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t b = 0; b < colOf.size(); ++b)
            out(i, b) = num_traits<T>::to_double(A(i, colOf[b]));
    return out;
}

/**
 * The row-replicated matrix in double, keeping only the classes that hold jobs,
 * with their populations returned as the column multiplicities.
 *
 * This is what `expand_to_double` expands, one step earlier: perm(Ar, mult)
 * equals perm of the expansion, and the saddle point wants the unexpanded form
 * because its expansion is asymptotic in mult. A class with no jobs is dropped
 * rather than passed with multiplicity zero, so a zero demand in such a column
 * cannot trip the full-support check.
 */
template <class T>
Matrix<double> rows_to_double(const Matrix<T>& A, const std::vector<int>& N,
                              std::vector<std::size_t>* mult) {
    std::vector<std::size_t> keep;
    mult->clear();
    for (std::size_t r = 0; r < N.size(); ++r)
        if (N[r] > 0) {
            keep.push_back(r);
            mult->push_back(static_cast<std::size_t>(N[r]));
        }
    Matrix<double> out(A.rows(), keep.size(), 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t l = 0; l < keep.size(); ++l)
            out(i, l) = num_traits<T>::to_double(A(i, keep[l]));
    return out;
}

/**
 * First (station, class) whose zero demand actually reaches the replicated
 * matrix, 1-based, or (0,0) when there is none.
 *
 * A class with no jobs or a station with no jobs contributes nothing, so its
 * zeros are irrelevant.
 */
template <class T>
std::pair<std::size_t, std::size_t> first_zero(const Matrix<T>& L, const std::vector<int>& N,
                                               const std::vector<int>& n) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < L.rows(); ++i) {
        if (n[i] == 0) continue;
        for (std::size_t r = 0; r < L.cols(); ++r)
            if (N[r] > 0 && !(L(i, r) > zero)) return std::make_pair(i + 1, r + 1);
    }
    return std::make_pair(static_cast<std::size_t>(0), static_cast<std::size_t>(0));
}

}  // namespace jointmargdetail

/**
 * @param n       (M) per-station total queue lengths, infinite servers included
 * @param L       (M x R) demand matrix, infinite-server rows included
 * @param N       (R) per-class populations
 * @param infset  rows of L that are infinite-server stations, 0-based
 * @param G       the normalizing constant G(N)
 * @param engine  "exact" (default), "spm", "bethe", "heur", "huberlaw" or
 *                "adapart". "spm" is the only engine that does NOT expand the
 *                matrix to order sum(N): it takes the row-replicated matrix with
 *                the class populations as column multiplicities, which is the
 *                regime its saddle-point expansion is asymptotically exact in, so
 *                its cost does not grow with the population and its relative error
 *                is O((R-1)/min(N)). Measured on a 3-station 2-class model, 12.8%
 *                at N = (1,1), 4.2% at (3,3), 2.1% at (6,6); it degrades the other
 *                way round, when the CLASS COUNT grows at fixed population (2.7%
 *                at R = 2, 21% at R = 7, both at N_r = 3), because R-1 is the
 *                dimension being expanded in. The bias is nearly constant across
 *                the lattice, so a caller that renormalizes a full sweep keeps far
 *                less of it: total variation distance 5.0e-3 at N = (1,1), 8.4e-4
 *                at (3,3), 4.3e-4 at (5,5), better than "bethe" and "heur" at
 *                every population measured.
 * @param seed    seed of the two sampling engines
 */
template <class T>
T pfqn_jointmarg(const std::vector<int>& n, const Matrix<T>& L, const std::vector<int>& N,
                 const std::vector<std::size_t>& infset, const T& G,
                 const std::string& engine = "exact", std::uint64_t seed = 0) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_jointmarg: L and N disagree on the class count");
    if (n.size() != M) throw InputError("pfqn_jointmarg: the occupancy vector has the wrong length");
    const T zero = num_traits<T>::from_int(0);
    if (G == zero) throw NumericError("pfqn_jointmarg: the normalizing constant is zero");
    for (std::size_t k = 0; k < infset.size(); ++k)
        if (infset[k] >= M) throw InputError("pfqn_jointmarg: infset indexes a station outside L");

    long Ntot = 0, ntot = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_jointmarg: negative population");
        Ntot += v;
    }
    for (int v : n) {
        if (v < 0) throw InputError("pfqn_jointmarg: negative occupancy");
        ntot += v;
    }
    // Infeasible occupancies are not an error: the caller sweeps a lattice.
    if (ntot != Ntot) return zero;
    if (Ntot == 0) return num_traits<T>::from_int(1) / G;

    const std::string eng = jointmargdetail::lower(engine);
    const Matrix<T> A = jointmargdetail::replicate_rows(L, n);

    T F;
    if (eng == "exact") {
        F = pfqn_perm(A, N);
    } else {
        const std::pair<std::size_t, std::size_t> z = jointmargdetail::first_zero(L, N, n);
        if (z.first != 0)
            throw InputError("pfqn_jointmarg: the '" + eng +
                             "' permanent engine cannot be applied: the demand of class " +
                             std::to_string(z.second) + " at station " + std::to_string(z.first) +
                             " is zero, so the replicated matrix has no full support. "
                             "Use engine 'exact'.");
        double f;
        if (eng == "spm") {
            // Never the expanded matrix: the saddle point is asymptotic in the
            // column multiplicities, which are the class populations themselves.
            std::vector<std::size_t> mult;
            const Matrix<double> Ar = jointmargdetail::rows_to_double(A, N, &mult);
            f = perm::perm_spm(Ar, mult);
        } else {
            const Matrix<double> Ad = jointmargdetail::expand_to_double(A, N);
            if (eng == "bethe") {
                f = perm::perm_bethe(Ad);
            } else if (eng == "heur") {
                f = perm::perm_heur(Ad);
            } else if (eng == "huberlaw") {
                f = perm::perm_huberlaw(Ad, seed);
            } else if (eng == "adapart") {
                f = perm::perm_adapart(Ad, seed);
            } else {
                throw InputError("pfqn_jointmarg: unrecognized permanent engine '" + engine +
                                 "'. Use exact, spm, bethe, heur, huberlaw or adapart.");
            }
        }
        F = num_traits<T>::from_double(f);
    }

    for (int v : N) F /= num_factorial<T>(static_cast<unsigned>(v));
    // Every infinite server keeps its own row, so every one of them divides by
    // its own n_j!; the queueing stations do not.
    for (std::size_t k = 0; k < infset.size(); ++k)
        F /= num_factorial<T>(static_cast<unsigned>(n[infset[k]]));
    return F / G;
}

/** Overload computing G with pfqn_ca first, matching the reference's default. */
template <class T>
T pfqn_jointmarg(const std::vector<int>& n, const Matrix<T>& L, const std::vector<int>& N,
                 const std::vector<std::size_t>& infset, const std::string& engine = "exact",
                 std::uint64_t seed = 0) {
    const std::size_t M = L.rows(), R = L.cols();
    // G does not depend on how the delay stations are split: they aggregate by
    // the multinomial theorem, so the constant may be taken with the
    // infinite-server rows summed into the think time.
    std::vector<bool> isinf(M, false);
    for (std::size_t k = 0; k < infset.size(); ++k)
        if (infset[k] < M) isinf[infset[k]] = true;
    std::size_t nq = 0;
    for (std::size_t i = 0; i < M; ++i)
        if (!isinf[i]) ++nq;
    Matrix<T> Lq(nq, R);
    std::size_t a = 0;
    for (std::size_t i = 0; i < M; ++i) {
        if (isinf[i]) continue;
        for (std::size_t r = 0; r < R; ++r) Lq(a, r) = L(i, r);
        ++a;
    }
    Matrix<T> Z;
    if (!infset.empty()) {
        Z = Matrix<T>(1, R, num_traits<T>::from_int(0));
        for (std::size_t k = 0; k < infset.size(); ++k)
            for (std::size_t r = 0; r < R; ++r) Z(0, r) += L(infset[k], r);
    }
    return pfqn_jointmarg(n, L, N, infset, pfqn_ca(Lq, N, Z).G, engine, seed);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_JOINTMARG_H
