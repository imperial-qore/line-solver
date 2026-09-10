/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MFQ_FLUFLU_SOJOURN_H
#define LINE_API_MAM_MFQ_FLUFLU_SOJOURN_H

/**
 * Sojourn-time distribution of a fluid queue whose SERVICE is itself a
 * Markov-modulated fluid flow, as an ME or PH representation (alpha, A).
 *
 * Port of matlab/src/api/mam/mfq_fluflu_sojourn.m and the BUTools FluFluSTD it
 * wraps. Arrivals are a fluid flow modulated by (Qin, Rin) and the server
 * drains fluid at a rate modulated by an independent (Qout, Rout).
 *
 * CONSTRUCTION. The two modulating chains are combined on the product space,
 * but NOT as a plain Kronecker sum: because the two flows run on different
 * clocks, the combined model is the TIME-CHANGED fluid queue
 *
 *     Rh = kron(Rin, I) - kron(I, Rout),
 *     Qh = kron(Qin, Rout) + kron(Rin, Qout),
 *
 * in which each chain's generator is weighted by the OTHER flow's rate. Qh is
 * therefore not a generator, and mfq_general_solve does not require it to be:
 * its algebra is generic. Solving that fluid model gives the level law of the
 * queue, and the sojourn-time representation follows from it exactly as in
 * mfq_sojourn, with the closing vector weighted by kron(Rin, I) when service
 * continues while the server's own fluid level is at zero, and by
 * kron(Rin, Rout)/mu when service stops there. That is the srv0stop flag, and
 * it is the only place where the two conventions differ.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental through
 * mfq_general_solve.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/mfq_solve.h"
#include "line/api/mam/mfq_sojourn.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Sojourn time of a drop in a fluid queue with fluid-modulated service.
 *
 * @param Qin       generator of the arrival-modulating chain, Na x Na
 * @param Rin       diagonal arrival fluid rate matrix, Na x Na
 * @param Qout      generator of the service-modulating chain, Ns x Ns
 * @param Rout      diagonal service fluid rate matrix, Ns x Ns
 * @param srv0stop  true if service stops when the server fluid level hits zero
 * @param transToPH true for a phase-type representation, false for ME
 * @param prec      tolerance handed to mfq_general_solve
 */
template <class T>
MeRepresentation<T> mfq_fluflu_sojourn(const Matrix<T>& Qin, const Matrix<T>& Rin,
                                       const Matrix<T>& Qout, const Matrix<T>& Rout,
                                       bool srv0stop, bool transToPH, const T& prec) {
    static_assert(num_traits<T>::has_transcendental,
                  "mfq_fluflu_sojourn requires transcendental arithmetic");
    using namespace mfq_detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t Na = Qin.rows(), Ns = Qout.rows();
    if (Qin.cols() != Na || Rin.rows() != Na || Rin.cols() != Na)
        throw InputError("mfq_fluflu_sojourn: Qin and Rin must be square and of equal order");
    if (Qout.cols() != Ns || Rout.rows() != Ns || Rout.cols() != Ns)
        throw InputError("mfq_fluflu_sojourn: Qout and Rout must be square and of equal order");

    const Matrix<T> Iin = eye<T>(Na);
    const Matrix<T> Iout = eye<T>(Ns);
    const Matrix<T> KRinI = kron(Rin, Iout);
    const Matrix<T> Rh = sub(KRinI, kron(Iin, Rout));
    const Matrix<T> Qh = add(kron(Qin, Rout), kron(Rin, Qout));

    const GeneralFluidSolution<T> s = mfq_general_solve(Qh, Rh, Matrix<T>(0, 0), prec);
    const std::size_t Np = s.ini.size();
    if (Np == 0) throw NumericError("mfq_fluflu_sojourn: the combined model has no up-drift state");

    // lambda and mu, the mean fluid rates of the two flows.
    T lambda = zero, mu = zero;
    {
        const std::vector<T> pin = mc::ctmc_solve(Qin);
        for (std::size_t i = 0; i < Na; ++i) lambda += pin[i] * Rin(i, i);
        const std::vector<T> pout = mc::ctmc_solve(Qout);
        for (std::size_t i = 0; i < Ns; ++i) mu += pout[i] * Rout(i, i);
    }
    if (lambda <= zero) throw NumericError("mfq_fluflu_sojourn: non-positive arrival fluid rate");
    if (srv0stop && mu <= zero)
        throw NumericError("mfq_fluflu_sojourn: non-positive service fluid rate");

    // The closing weight: kron(Rin, I)/lambda, or kron(Rin, Rout)/(lambda mu)
    // when service stops at a zero server level.
    Matrix<T> Wt = srv0stop ? kron(Rin, Rout) : KRinI;
    const T denom = srv0stop ? T(lambda * mu) : lambda;
    Wt = scale(Wt, T(one / denom));
    const Matrix<T> cloW = matmul(s.clo, Wt);  // Np x (Na Ns)

    const Matrix<T> negKinv = inverse(scale(s.K, T(-one)));
    MeRepresentation<T> out;
    if (transToPH) {
        const std::vector<T> delta = vecmul(s.ini, negKinv);
        for (std::size_t i = 0; i < Np; ++i)
            if (delta[i] == zero)
                throw NumericError(
                    "mfq_fluflu_sojourn: the PH similarity is singular on this model; ask for the "
                    "ME representation instead");
        // A = Delta^-1 Kh' Delta.
        out.A = Matrix<T>(Np, Np, zero);
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t k = 0; k < Np; ++k) out.A(i, k) = s.K(k, i) * delta[k] / delta[i];
        // alpha = row sums of Delta cloW.
        out.alpha.assign(Np, zero);
        for (std::size_t i = 0; i < Np; ++i) {
            T acc = zero;
            for (std::size_t j = 0; j < cloW.cols(); ++j) acc += delta[i] * cloW(i, j);
            out.alpha[i] = acc;
        }
    } else {
        std::vector<T> clovec(Np, zero);
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t j = 0; j < cloW.cols(); ++j) clovec[i] += cloW(i, j);
        const Matrix<T> B = mfq_transform_to_ones(clovec);
        const Matrix<T> Bi = inverse(B);
        out.A = matmul(matmul(B, s.K), Bi);
        out.alpha = vecmul(vecmul(s.ini, negKinv), Bi);
    }
    return out;
}

/** mfq_fluflu_sojourn with an ME representation and prec = 1e-14. */
template <class T>
MeRepresentation<T> mfq_fluflu_sojourn(const Matrix<T>& Qin, const Matrix<T>& Rin,
                                       const Matrix<T>& Qout, const Matrix<T>& Rout,
                                       bool srv0stop) {
    return mfq_fluflu_sojourn(Qin, Rin, Qout, Rout, srv0stop, false,
                              T(num_traits<T>::from_double(1e-14)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MFQ_FLUFLU_SOJOURN_H
