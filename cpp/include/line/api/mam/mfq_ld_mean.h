/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MFQ_LD_MEAN_H
#define LINE_API_MAM_MFQ_LD_MEAN_H

/**
 * Stationary mean fluid level E[X] of a first- or second-order level-dependent
 * (multi-regime) Markovian fluid queue, in closed form from the
 * matrix-exponential building blocks.
 *
 * Port of matlab/src/api/mam/mfq_ld_mean.m and the BUTools
 * LevelDependentFluidStationaryMean it wraps. Its inputs are exactly the
 * outputs of mfq_ld_solve, and it calls nothing else, so it is complete on its
 * own terms: a caller holding the blocks from any source can use it.
 *
 * THE DENSITY. Over regime k, that is over the level interval
 * [T(k), T(k+1)], the stationary density is the sum of a forward and a backward
 * matrix exponential,
 *
 *     pi_k(x) = iniF_k exp(KF_k (x - T(k))) cloF_k
 *             + iniB_k exp(KB_k (T(k+1) - x)) cloB_k,
 *
 * one anchored at each end of the regime, plus point masses at the K+1
 * thresholds. E[X] is then the mass contribution sum_j T(j) masses_j e plus
 * the integral of x pi_k(x) over each regime.
 *
 * THE INTEGRALS. Both regime integrals reduce to
 *   J0 = int_0^L exp(M u) du and J1 = int_0^L u exp(M u) du,
 * which are read off ONE matrix exponential of the nilpotent augmentation
 *
 *     A = [ M  I  0 ;  0  0  I ;  0  0  0 ],   W = exp(A L),
 *     J0 = W(1:n, n+1:2n),   J1 = L J0 - W(1:n, 2n+1:3n).
 *
 * That form is used rather than the obvious M^-1 (exp(M L) - I) because KF and
 * KB are SINGULAR whenever a regime has a zero-drift direction, which is the
 * normal case and not an edge case; the augmentation never inverts anything and
 * is exact for a singular M. The forward integral is anchored at T(k), giving
 * the weight T(k) J0 + J1, and the backward one runs the other way, giving
 * T(k+1) J0 - J1.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental because it calls expm,
 * which is a tolerance-controlled Pade approximation in any arithmetic. The
 * rest is finite exact linear algebra.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * The matrix-exponential building blocks of a multi-regime fluid queue, the
 * output of mfq_ld_solve. With K regimes there are K+1 point-mass vectors, at
 * the levels 0 = T(0), T(1), ..., T(K).
 */
template <class T>
struct LevelDependentFluidBlocks {
    std::vector<std::vector<T>> masses;  ///< K+1 point-mass vectors of length N
    std::vector<std::vector<T>> iniF;    ///< K forward initial vectors
    std::vector<Matrix<T>> KF;           ///< K forward matrix exponents
    std::vector<Matrix<T>> cloF;         ///< K forward closing matrices
    std::vector<std::vector<T>> iniB;    ///< K backward initial vectors
    std::vector<Matrix<T>> KB;           ///< K backward matrix exponents
    std::vector<Matrix<T>> cloB;         ///< K backward closing matrices
    std::vector<T> Thr;                  ///< K regime thresholds T(1)..T(K)
};

namespace mfq_ld_detail {

/**
 * J0 = int_0^L exp(M u) du and J1 = int_0^L u exp(M u) du, from one matrix
 * exponential of a nilpotent block augmentation. Valid for a singular M.
 */
template <class T>
void exp_int_moments(const Matrix<T>& M, const T& L, Matrix<T>& J0, Matrix<T>& J1) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t n = M.rows();
    Matrix<T> A(3 * n, 3 * n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) A(i, j) = M(i, j);
        A(i, n + i) = one;
        A(n + i, 2 * n + i) = one;
    }
    const Matrix<T> W = expm(A, L);
    J0 = Matrix<T>(n, n, zero);
    J1 = Matrix<T>(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            J0(i, j) = W(i, n + j);
            J1(i, j) = L * W(i, n + j) - W(i, 2 * n + j);
        }
}

}  // namespace mfq_ld_detail

/**
 * Stationary mean fluid level of a level-dependent fluid queue.
 *
 * @param b the building blocks returned by mfq_ld_solve
 */
template <class T>
T mfq_ld_mean(const LevelDependentFluidBlocks<T>& b) {
    static_assert(num_traits<T>::has_transcendental,
                  "mfq_ld_mean evaluates matrix exponentials");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = b.Thr.size();
    if (K == 0) throw InputError("mfq_ld_mean: at least one regime is required");
    if (b.masses.size() != K + 1)
        throw InputError("mfq_ld_mean: expected K+1 point-mass vectors");
    if (b.iniF.size() != K || b.KF.size() != K || b.cloF.size() != K || b.iniB.size() != K ||
        b.KB.size() != K || b.cloB.size() != K)
        throw InputError("mfq_ld_mean: expected K forward and K backward blocks");
    const std::size_t N = b.masses[0].size();

    // Thresholds with the implicit zero at the front: T(0) = 0.
    std::vector<T> Tv(K + 1, zero);
    for (std::size_t k = 0; k < K; ++k) Tv[k + 1] = b.Thr[k];

    T res = zero;
    // Contribution of the point masses, located at the levels T(0..K).
    for (std::size_t j = 0; j <= K; ++j) {
        if (b.masses[j].size() != N)
            throw InputError("mfq_ld_mean: the point-mass vectors have different lengths");
        T s = zero;
        for (const T& v : b.masses[j]) s += v;
        res += Tv[j] * s;
    }
    // Contribution of the continuous density in each regime.
    const std::vector<T> h = ones<T>(N);
    for (std::size_t k = 0; k < K; ++k) {
        const T L = Tv[k + 1] - Tv[k];
        Matrix<T> J0F, J1F, J0B, J1B;
        mfq_ld_detail::exp_int_moments(b.KF[k], L, J0F, J1F);
        mfq_ld_detail::exp_int_moments(b.KB[k], L, J0B, J1B);
        Matrix<T> WF = J1F, WB = J1B;
        for (std::size_t i = 0; i < WF.rows(); ++i)
            for (std::size_t j = 0; j < WF.cols(); ++j) WF(i, j) += Tv[k] * J0F(i, j);
        for (std::size_t i = 0; i < WB.rows(); ++i)
            for (std::size_t j = 0; j < WB.cols(); ++j)
                WB(i, j) = Tv[k + 1] * J0B(i, j) - WB(i, j);
        {
            const std::vector<T> t = vecmul(vecmul(b.iniF[k], WF), b.cloF[k]);
            for (std::size_t j = 0; j < N; ++j) res += t[j] * h[j];
        }
        {
            const std::vector<T> t = vecmul(vecmul(b.iniB[k], WB), b.cloB[k]);
            for (std::size_t j = 0; j < N; ++j) res += t[j] * h[j];
        }
    }
    return res;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MFQ_LD_MEAN_H
