/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_KMS_H
#define LINE_API_MC_CTMC_KMS_H

/**
 * Koury-McAllister-Stewart aggregation-disaggregation for a nearly completely
 * decomposable CTMC.
 *
 * Templated port of matlab/src/api/mc/ctmc_kms.m and
 * jar/src/main/java/jline/api/mc/Ctmc_kms.java. Starting from the Courtois
 * approximation, each sweep conditions the current iterate within every
 * macro-state, aggregates those conditional vectors into a macro-state chain,
 * solves it for the macro-state weights, and disaggregates by one block
 * Gauss-Seidel sweep on the uniformized matrix:
 *
 *   pn (D - U) = zn L,   D = blkdiag(I - P), L and U the strictly block lower
 *                        and upper parts of P, zn the reweighted conditional.
 *
 * Everything is done in the PERMUTED, macro-state-major index space, because
 * that is the space P and the block offsets live in; ctmc_courtois returns its
 * vector already mapped back to the original ordering, so the initial iterate
 * is permuted on entry and the result is unpermuted on exit. Mixing the two
 * spaces is the bug this port is written to avoid: with contiguous macro-states
 * the two orderings coincide and the error is invisible.
 *
 * As in the reference there is no convergence test: the caller asks for a fixed
 * number of sweeps. Above 6000 states the block solve is attempted by GMRES
 * first, exactly as MATLAB does, since that is where the direct factorization
 * stops fitting in memory; a non-zero GMRES flag falls back to the direct
 * solve, so the result never depends on whether the iteration converged.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC: it is seeded by ctmc_courtois, whose
 * epsMAX is an eigenvalue modulus, and its own large-system path is GMRES.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_courtois.h"
#include "line/api/mc/ctmc_bicgstab.h"
#include "line/api/mc/ctmc_gmres.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct KmsResult {
    std::vector<T> p;       ///< estimate after numSteps sweeps, ORIGINAL ordering
    std::vector<T> p_1;     ///< the previous iterate, ORIGINAL ordering
    std::vector<T> pcourt;  ///< the Courtois starting point, ORIGINAL ordering
    Matrix<T> Qperm;        ///< Q reordered by macro-state
    T eps;                  ///< NCD index, as ctmc_courtois defines it
    T epsMAX;               ///< maximum admissible NCD index
};

namespace detail {

/**
 * Solve x' M = rhs' for the row vector x, i.e. M' x = rhs. Uses GMRES first
 * above the reference's 6000-state dispatch threshold, and the direct solve
 * whenever GMRES does not converge.
 */
template <class T>
std::vector<T> aggregation_block_solve(const Matrix<T>& M, const std::vector<T>& rhs) {
    const std::size_t n = M.rows();
    Matrix<T> Mt(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Mt(i, j) = M(j, i);
    if (n > GMRES_MIN_STATES) {
        const GmresResult<T> g = ctmc_gmres(Mt, rhs);
        if (g.flag == 0) return g.x;
        // Short-recurrence retry before the cubic factorization, as in ctmc_solve.
        const BicgstabResult<T> bs = ctmc_bicgstab(Mt, rhs);
        if (bs.flag == 0) return bs.x;
    }
    return solve(Mt, rhs);
}

}  // namespace detail

/**
 * @param Q        generator
 * @param MS       macro-states partitioning 0..n-1
 * @param numSteps number of aggregation-disaggregation sweeps
 */
template <class T>
KmsResult<T> ctmc_kms(const Matrix<T>& Q, const std::vector<std::vector<std::size_t>>& MS,
                      std::size_t numSteps) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_kms requires transcendental arithmetic: it is seeded by ctmc_courtois, "
                  "whose epsMAX is an eigenvalue modulus, and its large-system block solve is "
                  "GMRES, which stops on a residual tolerance");
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_kms: generator is not square");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t nMacro = MS.size();

    const detail::CourtoisCore<T> c = detail::courtois_core(Q, MS, detail::courtois_default_rate(Q, MS));
    const std::vector<T> pMacro = dtmc_solve(c.G);

    // Block offsets in the permuted ordering.
    std::vector<std::size_t> off(nMacro + 1, 0);
    for (std::size_t i = 0; i < nMacro; ++i) off[i + 1] = off[i] + MS[i].size();

    std::vector<T> pn(n, zero);
    for (std::size_t i = 0; i < nMacro; ++i)
        for (std::size_t a = off[i]; a < off[i + 1]; ++a) pn[a] = pMacro[i] * c.pmicro[a];

    KmsResult<T> r;
    r.pcourt = detail::unpermute_states(pn, c.v);
    r.Qperm = c.Qperm;
    r.eps = c.eps;
    r.epsMAX = c.epsMAX;

    std::vector<T> pn_1 = pn;
    for (std::size_t step = 0; step < numSteps; ++step) {
        pn_1 = pn;

        // Aggregation: condition within each macro-state, then lump.
        std::vector<T> pcond = pn_1;
        for (std::size_t I = 0; I < nMacro; ++I) {
            T s = zero;
            for (std::size_t a = off[I]; a < off[I + 1]; ++a) s += pn_1[a];
            if (s > zero)
                for (std::size_t a = off[I]; a < off[I + 1]; ++a) pcond[a] /= s;
        }

        Matrix<T> G(nMacro, nMacro, zero);
        for (std::size_t I = 0; I < nMacro; ++I)
            for (std::size_t J = 0; J < nMacro; ++J) {
                T acc = zero;
                for (std::size_t b = off[J]; b < off[J + 1]; ++b) {
                    T s = zero;
                    for (std::size_t a = off[I]; a < off[I + 1]; ++a) s += c.P(b, a);
                    acc += pcond[b] * s;
                }
                G(I, J) = acc;
            }
        Matrix<T> Gt(nMacro, nMacro);
        for (std::size_t i = 0; i < nMacro; ++i)
            for (std::size_t j = 0; j < nMacro; ++j) Gt(i, j) = G(j, i);
        const std::vector<T> w = dtmc_solve(Gt);

        // Disaggregation: one block Gauss-Seidel sweep, pn (D - U) = zn L.
        std::vector<T> zn(n, zero);
        for (std::size_t I = 0; I < nMacro; ++I)
            for (std::size_t a = off[I]; a < off[I + 1]; ++a) zn[a] = w[I] * pcond[a];

        Matrix<T> M(n, n, zero);
        std::vector<T> rhs(n, zero);
        for (std::size_t I = 0; I < nMacro; ++I)
            for (std::size_t J = 0; J < nMacro; ++J)
                for (std::size_t a = off[I]; a < off[I + 1]; ++a)
                    for (std::size_t b = off[J]; b < off[J + 1]; ++b) {
                        if (I > J) {
                            rhs[b] += zn[a] * c.P(a, b);  // (zn L)_b
                        } else if (I == J) {
                            M(a, b) = (a == b ? one : zero) - c.P(a, b);
                        } else {
                            M(a, b) = -c.P(a, b);  // D - U
                        }
                    }

        pn = detail::aggregation_block_solve(M, rhs);
        T tot = zero;
        for (const T& x : pn) tot += x;
        if (tot == zero) throw NumericError("ctmc_kms: the disaggregation sweep returned a null vector");
        for (T& x : pn) x /= tot;
    }

    r.p = detail::unpermute_states(pn, c.v);
    r.p_1 = detail::unpermute_states(pn_1, c.v);
    return r;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_KMS_H
