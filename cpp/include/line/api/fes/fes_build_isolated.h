/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_BUILD_ISOLATED_H
#define LINE_API_FES_BUILD_ISOLATED_H

/**
 * Service demands and visit ratios of an isolated subnetwork, from the
 * stochastic complement of its routing matrix.
 *
 * Templated port of matlab/src/api/fes/fes_build_isolated.m. Given the
 * stochastic complement S of the routing chain restricted to a subset of
 * stations, the per-class visit ratios are the stationary distribution of the
 * embedded DTMC of that class,
 *
 *     v_k P_k = v_k,   sum_i v_k(i) = 1,
 *
 * and the demands follow as L(i,k) = v_k(i) / rate(i,k), renormalized so that
 * the first station of the subset is visited once, which is the MVA
 * convention the flow-equivalent-server construction expects.
 *
 * THE sn DEPENDENCY IS LIFTED INTO THE SIGNATURE. The MATLAB entry point takes
 * `(sn, subsetIndices, stochCompS)` and opens with a loop that plucks five
 * fields out of `sn` -- `nclasses`, `stationToNode`, `nodetype`, `nservers`,
 * `rates` -- to produce the per-station server counts, the delay flags and the
 * service rates of the subset. Nothing after that loop reads `sn`, and the
 * routing already arrives as a plain matrix argument. The extraction is
 * therefore not part of the algorithm: this port takes the extracted
 * quantities directly, in the same spirit as npfqn_sqd.
 *
 * Note in particular that `mi` and `isDelay` are PURE PASS-THROUGH in the
 * reference: they are assembled from `sn` and returned, and no later line
 * reads them. They are not arguments here for that reason -- a caller that
 * wants them holds them already. The numeric core needs only the service
 * rates and the stochastic complement.
 *
 * INDEXING. stochCompS is indexed (station-1)*K + class over the SUBSET, so
 * it is (M_sub K x M_sub K); the reference forms it that way and indexes it
 * with the subset position i, not the original station id. Its class-k block
 * is read at rows and columns (i-1)*K + k.
 *
 * ARITHMETIC. Exact at Rational. The stationary distribution is obtained from
 * the singular system (I - P_k' + e e'/M) x = e/M, one linear solve per class,
 * which is a finite sequence of field operations; there is no iteration and no
 * tolerance in the solve itself. The two guards that DO carry a tolerance are
 * reproduced from the reference and are structural rather than numerical: a
 * row of P_k whose sum is below FineTol = 1e-8 is treated as "no routing
 * defined" and made a self-loop, and a row whose sum differs from one by more
 * than FineTol is renormalized. At Rational a caller supplying an exactly
 * stochastic complement never reaches either branch.
 *
 * DIVERGENCE, stated because it changes which inputs take the fallback path.
 * The reference guards the solve with `rank(A) == M_sub` and substitutes the
 * uniform distribution when the rank is deficient. A rank test is a singular
 * value decomposition, which is double-only in this tree and has no exact
 * instantiation (util/eig.h documents why). The port instead attempts the
 * solve and takes the same uniform fallback when the factorization finds an
 * exactly zero pivot. The two agree on every non-degenerate input: A is
 * nonsingular precisely when the class-k chain has a unique stationary
 * distribution, which is when the solve succeeds. They can differ only for a
 * matrix that is numerically rank deficient yet still factorizable, where the
 * reference falls back and the port returns the (ill-conditioned) solve; that
 * is the direction that preserves information rather than discarding it.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace fes {

/** Demands and visit ratios of the isolated subnetwork. */
template <class T>
struct FesIsolated {
    Matrix<T> L;       ///< (M_sub x K) service demands
    Matrix<T> visits;  ///< (M_sub x K) visit ratios, each class summing to one
};

/**
 * Build the isolated subnetwork's demands and visit ratios.
 *
 * @param rates      (M_sub x K) service rates of the subset stations; an entry
 *                   that is zero, negative or non-finite marks a disabled
 *                   service and yields a zero demand, as in the reference
 * @param stochCompS (M_sub K x M_sub K) stochastic complement of the routing
 *                   chain over the subset, indexed (i-1)*K + k
 */
template <class T>
FesIsolated<T> fes_build_isolated(const Matrix<T>& rates, const Matrix<T>& stochCompS) {
    const std::size_t M = rates.rows();
    const std::size_t K = rates.cols();
    if (M == 0) throw InputError("fes_build_isolated: empty station subset");
    if (K == 0) throw InputError("fes_build_isolated: no classes");
    if (stochCompS.rows() < M * K || stochCompS.cols() < M * K)
        throw InputError(
            "fes_build_isolated: the stochastic complement is too small for the subset; it must be "
            "at least (M_sub*K) square, indexed (station-1)*K + class over the SUBSET");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T fineTol = num_traits<T>::from_double(1e-8);
    const T Mt = num_traits<T>::from_int(static_cast<long>(M));

    FesIsolated<T> out;
    out.L = Matrix<T>(M, K, zero);
    out.visits = Matrix<T>(M, K, zero);

    for (std::size_t k = 0; k < K; ++k) {
        // Class-k routing block of the stochastic complement.
        Matrix<T> P(M, M, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) P(i, j) = stochCompS(i * K + k, j * K + k);

        // Row repair, exactly as the reference: a row with no routing becomes
        // a self-loop, a row that is off by more than FineTol is renormalized.
        for (std::size_t i = 0; i < M; ++i) {
            T rowsum = zero;
            for (std::size_t j = 0; j < M; ++j) rowsum += P(i, j);
            if (rowsum > fineTol && num_abs(T(rowsum - one)) > fineTol) {
                for (std::size_t j = 0; j < M; ++j) P(i, j) /= rowsum;
            } else if (rowsum < fineTol) {
                for (std::size_t j = 0; j < M; ++j) P(i, j) = zero;
                P(i, i) = one;
            }
        }

        // Stationary distribution from (I - P' + e e'/M) x = e/M.
        Matrix<T> A(M, M, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) {
                const T id = (i == j) ? one : zero;
                A(i, j) = id - P(j, i) + one / Mt;
            }
        std::vector<T> rhs(M, T(one / Mt));
        std::vector<T> pi;
        bool ok = true;
        try {
            pi = solve(A, rhs);
        } catch (const NumericError&) {
            ok = false;  // exactly singular: the reference's rank-deficient branch
        }
        if (!ok) {
            pi.assign(M, T(one / Mt));
        } else {
            // Clamp and renormalize, as the reference does.
            for (std::size_t i = 0; i < M; ++i)
                if (pi[i] < zero) pi[i] = zero;
            T s = zero;
            for (std::size_t i = 0; i < M; ++i) s += pi[i];
            if (s > zero) {
                for (std::size_t i = 0; i < M; ++i) pi[i] /= s;
            } else {
                pi.assign(M, T(one / Mt));
            }
        }
        for (std::size_t i = 0; i < M; ++i) out.visits(i, k) = pi[i];
    }

    // Demands L(i,k) = visits(i,k) / rate(i,k), zero where service is disabled.
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const T r = rates(i, k);
            if (!(r > zero) || !std::isfinite(num_traits<T>::to_double(r)))
                out.L(i, k) = zero;
            else
                out.L(i, k) = out.visits(i, k) / r;
        }

    // Renormalize so the first station of the subset is visited once, the MVA
    // convention. A class that never visits station 0 is left as it is.
    for (std::size_t k = 0; k < K; ++k) {
        const T v0 = out.visits(0, k);
        if (v0 > zero)
            for (std::size_t i = 0; i < M; ++i) out.L(i, k) /= v0;
    }
    return out;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_BUILD_ISOLATED_H
