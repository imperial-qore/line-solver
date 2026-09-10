/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_TAKAHASHI_H
#define LINE_API_MC_CTMC_TAKAHASHI_H

/**
 * Takahashi's aggregation-disaggregation for a nearly completely decomposable
 * CTMC.
 *
 * Templated port of matlab/src/api/mc/ctmc_takahashi.m and
 * jar/src/main/java/jline/api/mc/Ctmc_takahashi.java. Like KMS it starts from
 * the Courtois approximation and alternates aggregation with disaggregation,
 * but the disaggregation step is a per-macro-state fixed point rather than a
 * Gauss-Seidel sweep: for macro-state I,
 *
 *   (I - P_II') x_I = b_I,   b_I(i) = sum_{K != I} gamma_K GI(K, i),
 *
 * with gamma the aggregated macro-state distribution and GI(K, .) the
 * conditional one-step flow out of macro-state K.
 *
 * Everything here is in the ORIGINAL state ordering: P is the uniformization
 * of the UNPERMUTED Q and each macro-state is addressed by the indices MS[I]
 * carries. Contiguous block offsets would silently solve for a different
 * partition of the same block sizes whenever a macro-state is not a contiguous
 * range, which is the defect the MATLAB version records in its own comments.
 *
 * The uniformization rate is passed explicitly as (21/20) max|Q|, matching the
 * rate ctmc_courtois derives. The reference used to leave it to the
 * ctmc_randomization default, max|Q| + rand, which made the whole iteration
 * irreproducible for no benefit: the aggregation and disaggregation equations
 * are homogeneous in P - I = Q/q, so the fixed point does not depend on q.
 *
 * REFERENCE DEFECT, repaired here. MATLAB guards the aggregation against an
 * empty macro-state (S > 1e-14) but then divides by that same S unguarded when
 * building GI, so a macro-state carrying no mass poisons the whole iterate with
 * NaN. Here a macro-state whose current mass is below the threshold contributes
 * nothing to either G or GI, which is the limit of the expression as S -> 0.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC: seeded by ctmc_courtois, whose epsMAX is
 * an eigenvalue modulus, and its large-system block solve is GMRES.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_bicgstab.h"
#include "line/api/mc/ctmc_courtois.h"
#include "line/api/mc/ctmc_kms.h"
#include "line/api/mc/ctmc_randomization.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct TakahashiResult {
    std::vector<T> p;       ///< estimate after numSteps sweeps
    std::vector<T> p_1;     ///< the previous iterate
    std::vector<T> pcourt;  ///< the Courtois starting point
    Matrix<T> Qperm;        ///< Q reordered by macro-state
    T eps;                  ///< NCD index, as ctmc_courtois defines it
    T epsMAX;               ///< maximum admissible NCD index
};

/**
 * @param Q        generator
 * @param MS       macro-states partitioning 0..n-1
 * @param numSteps number of aggregation-disaggregation sweeps
 * @param massTol  macro-state mass below which its contribution is dropped
 */
template <class T>
TakahashiResult<T> ctmc_takahashi(const Matrix<T>& Q, const std::vector<std::vector<std::size_t>>& MS,
                                  std::size_t numSteps, double massTol = 1e-14) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_takahashi requires transcendental arithmetic: it is seeded by "
                  "ctmc_courtois, whose epsMAX is an eigenvalue modulus, and its large-system "
                  "block solve is GMRES, which stops on a residual tolerance");
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_takahashi: generator is not square");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t nMacro = MS.size();
    const T tol = num_traits<T>::from_double(massTol);

    const detail::CourtoisCore<T> c = detail::courtois_core(Q, MS, detail::courtois_default_rate(Q, MS));
    const std::vector<T> pMacro = dtmc_solve(c.G);

    std::vector<std::size_t> off(nMacro + 1, 0);
    for (std::size_t i = 0; i < nMacro; ++i) off[i + 1] = off[i] + MS[i].size();
    std::vector<T> pperm(n, zero);
    for (std::size_t i = 0; i < nMacro; ++i)
        for (std::size_t a = off[i]; a < off[i + 1]; ++a) pperm[a] = pMacro[i] * c.pmicro[a];

    TakahashiResult<T> r;
    r.pcourt = detail::unpermute_states(pperm, c.v);
    r.Qperm = c.Qperm;
    r.eps = c.eps;
    r.epsMAX = c.epsMAX;

    // Uniformization of the UNPERMUTED generator, at the rate Courtois derives.
    const T qmax = ctmc_maxabs(Q);
    if (qmax == zero) throw InputError("ctmc_takahashi: the generator has no transitions");
    const Matrix<T> P = ctmc_randomization(Q, T(qmax * num_traits<T>::from_rational(21, 20))).P;

    std::vector<T> pn = r.pcourt, pn_1 = r.pcourt;
    for (std::size_t step = 0; step < numSteps; ++step) {
        pn_1 = pn;

        // Aggregation.
        Matrix<T> G(nMacro, nMacro, zero);
        std::vector<T> S(nMacro, zero);
        for (std::size_t I = 0; I < nMacro; ++I)
            for (std::size_t i : MS[I]) S[I] += pn_1[i];
        for (std::size_t I = 0; I < nMacro; ++I) {
            if (!(S[I] > tol)) continue;
            for (std::size_t J = 0; J < nMacro; ++J) {
                if (I == J) continue;
                T acc = zero;
                for (std::size_t i : MS[I])
                    for (std::size_t j : MS[J]) acc += P(i, j) * pn_1[i] / S[I];
                G(I, J) = acc;
            }
        }
        for (std::size_t I = 0; I < nMacro; ++I) {
            T rs = zero;
            for (std::size_t J = 0; J < nMacro; ++J)
                if (J != I) rs += G(I, J);
            G(I, I) = one - rs;
        }
        const std::vector<T> gamma = dtmc_solve(G);

        // Conditional one-step flow out of each macro-state.
        Matrix<T> GI(nMacro, n, zero);
        for (std::size_t I = 0; I < nMacro; ++I) {
            if (!(S[I] > tol)) continue;
            for (std::size_t j = 0; j < n; ++j) {
                T acc = zero;
                for (std::size_t i : MS[I]) acc += P(i, j) * pn_1[i];
                GI(I, j) = acc / S[I];
            }
        }

        // Disaggregation, one macro-state at a time.
        for (std::size_t I = 0; I < nMacro; ++I) {
            const std::size_t sz = MS[I].size();
            Matrix<T> A(sz, sz, zero);
            std::vector<T> b(sz, zero);
            for (std::size_t i = 0; i < sz; ++i) {
                for (std::size_t j = 0; j < sz; ++j)
                    A(i, j) = (i == j ? one : zero) - P(MS[I][j], MS[I][i]);
                for (std::size_t K = 0; K < nMacro; ++K)
                    if (K != I) b[i] += gamma[K] * GI(K, MS[I][i]);
            }
            std::vector<T> xI;
            if (sz > GMRES_MIN_STATES) {
                const GmresResult<T> g = ctmc_gmres(A, b);
                if (g.flag == 0) {
                    xI = g.x;
                } else {
                    // Short-recurrence retry before the cubic factorization, as in ctmc_solve.
                    const BicgstabResult<T> bs = ctmc_bicgstab(A, b);
                    if (bs.flag == 0) xI = bs.x;
                }
            }
            if (xI.empty()) xI = solve(A, b);
            for (std::size_t i = 0; i < sz; ++i) pn[MS[I][i]] = xI[i];
        }

        T tot = zero;
        for (const T& x : pn) tot += x;
        if (tot == zero) throw NumericError("ctmc_takahashi: the disaggregation step returned a null vector");
        for (T& x : pn) x /= tot;
    }

    r.p = pn;
    r.p_1 = pn_1;
    return r;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_TAKAHASHI_H
