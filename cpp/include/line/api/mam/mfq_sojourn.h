/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MFQ_SOJOURN_H
#define LINE_API_MAM_MFQ_SOJOURN_H

/**
 * Sojourn-time distribution of a Markov-modulated fluid queue, as a
 * matrix-exponential or phase-type representation (alpha, A).
 *
 * Port of matlab/src/api/mam/mfq_sojourn.m and the BUTools FluidQueueSTD it
 * wraps. A background chain with generator Q modulates an input fluid rate
 * matrix Rin and an output (service) rate matrix Rout, both diagonal. The
 * quantity returned is the law of the time a fluid DROP spends in the queue.
 *
 * CONSTRUCTION. The net drift is Rin - Rout, so mfq_general_solve returns the
 * stationary level law as a point mass mass0 at zero plus ini exp(K x) clo. The
 * arrival rate of fluid is
 *
 *     lambda = sum( mass0 Rin + ini (-K)^-1 clo Rin ),
 *
 * the mass at zero and the density integrated over the level each weighted by
 * the input rate. A drop that arrives when the level is x and the background
 * state is j leaves after the queue has drained x at the state-dependent rate,
 * so its sojourn time is the first passage of a level-dependent process whose
 * generator on the product space (background state, K-phase) is
 *
 *     A = B ( kron(Q', I) + kron(Rout, K) ) B^-1        (ME form)
 *     A = kron(Rout, Delta^-1 K' Delta) + kron(Q, I)    (PH form)
 *
 * of order N * Np, with N the background order and Np the number of up-drift
 * states. The two forms are the same operator in two different bases: the ME
 * form uses the similarity that maps the closing vector to a vector of ones,
 * the PH form the diagonal similarity Delta = diag(ini (-K)^-1)/lambda, which
 * makes the representation a genuine phase type (non-negative generator, so it
 * can be sampled and fed to any PH consumer) at the cost of requiring that
 * diagonal to be positive.
 *
 * WHICH FORM TO ASK FOR. The ME form always exists; the PH form exists only
 * when Delta is invertible, and neither the reference nor this port can promise
 * that in advance. Ask for PH when the result must be sampled or handed to a
 * PH-only routine, ME otherwise. The two give the SAME distribution and the
 * tests assert exactly that, by comparing their first three moments and their
 * CDFs -- which is a real check, since the two are computed through disjoint
 * code paths.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental through
 * mfq_general_solve, whose Riccati iteration terminates on a tolerance.
 * Everything on top of it is finite exact linear algebra.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/mfq_solve.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** A matrix-exponential or phase-type representation (alpha, A). */
template <class T>
struct MeRepresentation {
    std::vector<T> alpha;  ///< initial row vector
    Matrix<T> A;           ///< generator-like matrix exponent
};

/**
 * Sojourn time of a drop in a Markov-modulated fluid queue.
 *
 * @param Q         generator of the background chain, N x N
 * @param Rin       diagonal input fluid rate matrix, N x N
 * @param Rout      diagonal output (service) fluid rate matrix, N x N
 * @param Q0        level-zero generator; an empty matrix means Q0 = Q
 * @param transToPH true for a phase-type representation, false for ME
 * @param prec      tolerance handed to mfq_general_solve
 */
template <class T>
MeRepresentation<T> mfq_sojourn(const Matrix<T>& Q, const Matrix<T>& Rin, const Matrix<T>& Rout,
                                const Matrix<T>& Q0, bool transToPH, const T& prec) {
    static_assert(num_traits<T>::has_transcendental,
                  "mfq_sojourn requires transcendental arithmetic");
    using namespace mfq_detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t N = Q.rows();
    if (Q.cols() != N || Rin.rows() != N || Rin.cols() != N || Rout.rows() != N ||
        Rout.cols() != N)
        throw InputError("mfq_sojourn: Q, Rin and Rout must be square and of equal order");

    const Matrix<T> Rnet = sub(Rin, Rout);
    const GeneralFluidSolution<T> s = mfq_general_solve(Q, Rnet, Q0, prec);
    const std::size_t Np = s.ini.size();
    if (Np == 0) throw NumericError("mfq_sojourn: the fluid model has no up-drift states");

    const Matrix<T> negKinv = inverse(scale(s.K, T(-one)));
    const std::vector<T> iniKi = vecmul(s.ini, negKinv);  // ini (-K)^-1

    // lambda = sum( mass0 Rin + iniKi clo Rin ).
    T lambda = zero;
    {
        const std::vector<T> a = vecmul(s.mass0, Rin);
        const std::vector<T> b = vecmul(vecmul(iniKi, s.clo), Rin);
        for (std::size_t j = 0; j < N; ++j) lambda += a[j] + b[j];
    }
    if (lambda <= zero) throw NumericError("mfq_sojourn: non-positive fluid arrival rate");

    MeRepresentation<T> out;
    const std::size_t n = N * Np;
    if (transToPH) {
        std::vector<T> delta(Np);
        for (std::size_t i = 0; i < Np; ++i) {
            delta[i] = iniKi[i] / lambda;
            if (delta[i] == zero)
                throw NumericError(
                    "mfq_sojourn: the PH similarity is singular on this model; ask for the ME "
                    "representation instead");
        }
        // alpha = reshape(clo Rin, 1, N Np) * kron(I_N, Delta), column-major:
        // index (j, i) -> j Np + i.
        const Matrix<T> cloRin = matmul(s.clo, Rin);
        out.alpha.assign(n, zero);
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t i = 0; i < Np; ++i) out.alpha[j * Np + i] = cloRin(i, j) * delta[i];
        // A = kron(Rout, Delta^-1 K' Delta) + kron(Q, I_Np).
        Matrix<T> Mi(Np, Np, zero);
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t k = 0; k < Np; ++k) Mi(i, k) = s.K(k, i) * delta[k] / delta[i];
        out.A = Matrix<T>(n, n, zero);
        for (std::size_t j1 = 0; j1 < N; ++j1)
            for (std::size_t j2 = 0; j2 < N; ++j2)
                for (std::size_t i1 = 0; i1 < Np; ++i1)
                    for (std::size_t i2 = 0; i2 < Np; ++i2) {
                        T v = Rout(j1, j2) * Mi(i1, i2);
                        if (i1 == i2) v += Q(j1, j2);
                        out.A(j1 * Np + i1, j2 * Np + i2) = v;
                    }
    } else {
        // B maps the closing vector reshape((-K)^-1 clo Rin, N Np, 1) to ones.
        const Matrix<T> W = matmul(matmul(negKinv, s.clo), Rin);  // Np x N
        std::vector<T> clovec(n, zero);
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t i = 0; i < Np; ++i) clovec[j * Np + i] = W(i, j);
        const Matrix<T> B = mfq_transform_to_ones(clovec);
        const Matrix<T> Bi = inverse(B);
        // alpha = kron(ones(1,N), ini/lambda) * B^-1.
        std::vector<T> rep(n, zero);
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t i = 0; i < Np; ++i) rep[j * Np + i] = s.ini[i] / lambda;
        out.alpha = vecmul(rep, Bi);
        // A = B ( kron(Q', I_Np) + kron(Rout, K) ) B^-1.
        Matrix<T> Mid(n, n, zero);
        for (std::size_t j1 = 0; j1 < N; ++j1)
            for (std::size_t j2 = 0; j2 < N; ++j2)
                for (std::size_t i1 = 0; i1 < Np; ++i1)
                    for (std::size_t i2 = 0; i2 < Np; ++i2) {
                        T v = Rout(j1, j2) * s.K(i1, i2);
                        if (i1 == i2) v += Q(j2, j1);
                        Mid(j1 * Np + i1, j2 * Np + i2) = v;
                    }
        out.A = matmul(matmul(B, Mid), Bi);
    }
    return out;
}

/** mfq_sojourn with the regular boundary, an ME representation and prec = 1e-14. */
template <class T>
MeRepresentation<T> mfq_sojourn(const Matrix<T>& Q, const Matrix<T>& Rin, const Matrix<T>& Rout) {
    return mfq_sojourn(Q, Rin, Rout, Matrix<T>(0, 0), false,
                       T(num_traits<T>::from_double(1e-14)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MFQ_SOJOURN_H
