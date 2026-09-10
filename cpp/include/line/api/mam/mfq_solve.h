/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MFQ_SOLVE_H
#define LINE_API_MAM_MFQ_SOLVE_H

/**
 * Core of the Markovian fluid queue: the fundamental matrices Psi, K, U and the
 * matrix-exponential stationary solution of a general fluid model.
 *
 * These are the two BUTools routines that every mfq_* entry point of
 * matlab/src/api/mam/ ultimately calls (FluidFundamentalMatrices and
 * GeneralFluidSolve). They carry no LINE-level name of their own, so they live
 * here as shared infrastructure rather than as a ported API function, and
 * mfq_sojourn.h and mfq_fluflu_sojourn.h are the thin entry points over them.
 *
 * THE MODEL. A background chain with generator Q modulates a fluid level whose
 * drift in state i is R(i,i), of any sign, zero included. Writing the states in
 * the order (zero drift, positive drift, negative drift), censoring out the
 * zero-drift states and rescaling time by |1/R| turns the model into a fluid
 * queue with drifts +-1, whose first passage matrix Psi solves the Riccati
 * equation
 *
 *     Fpm + Fpp Psi + Psi Fmm + Psi Fmp Psi = 0.
 *
 * Psi(i,j) is the probability that, starting in an up state i at a level, the
 * process returns to that level in the down state j. From it, K = Fpp + Psi Fmp
 * generates the stationary density and U = Fmm + Fmp Psi is the generator of
 * the chain seen at level zero. The stationary law is then a point mass at
 * level zero plus the density pi(x) = ini exp(K x) clo.
 *
 * THE RICCATI SOLVER. ADDA (Wang, Wang and Li 2011), the BUTools default: an
 * alternating-directional doubling iteration whose per-step cost is four
 * inverses and whose error SQUARES each step, so the 1e-14 default is reached
 * in a few tens of iterations rather than the hundreds a linearly convergent
 * fixed point would need. SDA is the same iteration with a common shift and is
 * offered because it is the form quoted in most of the literature; the two
 * differ only in the shift and in the scaling of E and F. The reference's third
 * option, cyclic reduction, is NOT ported: it is an alternative with the same
 * convergence order and no accuracy advantage on this problem, and it would
 * duplicate the whole solver for nothing. Callers asking for it get an
 * UnsupportedError naming it, rather than a silent substitution.
 *
 * WHERE THE PORT IS STRICTER THAN THE REFERENCE. GeneralFluidSolve censors the
 * zero-drift states through pinv(-Qv00), a PSEUDO-inverse. That matrix is
 * singular exactly when the zero-drift states form a closed set, i.e. when the
 * fluid can be trapped at a constant level forever and the model has no
 * stationary fluid law to report. The reference then returns whatever the
 * pseudo-inverse yields, without comment. The port uses the ordinary inverse
 * and raises NumericError naming the condition, because a pseudo-inverse of a
 * singular censoring operator is not the answer to a different question, it is
 * an answer to no question. On every non-degenerate model the two agree
 * exactly, the pseudo-inverse of a non-singular matrix being its inverse.
 *
 * The two overdetermined normalizations that the reference solves with MATLAB's
 * backslash (an (n+1) x n system that is consistent by construction) are solved
 * here through the normal equations, which return the same vector for a
 * consistent full-column-rank system and need no least-squares factorization.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental: ADDA terminates on a
 * tolerance and its scaling step takes a square root. Everything else is finite
 * exact linear algebra.
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
namespace mam {

/** Which doubling iteration to run for Psi. */
enum class RiccatiMethod { ADDA, SDA };

/** Psi, K and U of a fluid queue with drifts normalized to +-1. */
template <class T>
struct FluidFundamental {
    Matrix<T> Psi;     ///< first return matrix, up states to down states
    Matrix<T> K;       ///< Fpp + Psi Fmp, the density generator
    Matrix<T> U;       ///< Fmm + Fmp Psi, the level-zero generator
    unsigned iterations;
    bool converged;
};

namespace mfq_detail {

/** Sub-block A(r0..r0+nr-1, c0..c0+nc-1). */
template <class T>
Matrix<T> block(const Matrix<T>& A, std::size_t r0, std::size_t c0, std::size_t nr,
                std::size_t nc) {
    Matrix<T> B(nr, nc, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < nr; ++i)
        for (std::size_t j = 0; j < nc; ++j) B(i, j) = A(r0 + i, c0 + j);
    return B;
}

/** Elementwise A + B. */
template <class T>
Matrix<T> add(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) += B(i, j);
    return C;
}

/** Elementwise A - B. */
template <class T>
Matrix<T> sub(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) -= B(i, j);
    return C;
}

/** A scaled by s. */
template <class T>
Matrix<T> scale(const Matrix<T>& A, const T& s) {
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) *= s;
    return C;
}

/** MATLAB norm(A, 1): the largest absolute column sum. */
template <class T>
T norm1(const Matrix<T>& A) {
    T best = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < A.cols(); ++j) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < A.rows(); ++i) s += num_abs(A(i, j));
        if (s > best) best = s;
    }
    return best;
}

/**
 * Solution of the CONSISTENT overdetermined system M x = b through the normal
 * equations M^T M x = M^T b. Used for the two normalizations that the reference
 * writes as a MATLAB backslash on an (n+1) x n system.
 *
 * Deliberately NOT line::lstsq from util/lstsq.h, and deliberately not named
 * lstsq either. Two reasons, in order of importance. The systems here are
 * consistent and of full column rank BY CONSTRUCTION -- they are a square
 * balance system plus one normalization row -- so the normal equations are
 * exact for them and the rank-revealing SVD path buys nothing; and squaring the
 * condition number, the usual objection to normal equations, is not a concern
 * on a system whose extra row is a normalization. Switching solvers would also
 * perturb values that currently agree with MATLAB digit for digit, for no
 * demonstrated gain. The distinct name additionally keeps argument-dependent
 * lookup from finding both this and line::lstsq on a Matrix<T> argument, which
 * is an ambiguity rather than an overload.
 */
template <class T>
std::vector<T> normal_equations_solve(const Matrix<T>& M, const std::vector<T>& b) {
    const std::size_t m = M.rows(), n = M.cols();
    if (b.size() != m)
        throw InputError("mfq normal_equations_solve: right-hand side length mismatch");
    Matrix<T> N(n, n, num_traits<T>::from_int(0));
    std::vector<T> r(n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            T s = num_traits<T>::from_int(0);
            for (std::size_t k = 0; k < m; ++k) s += M(k, i) * M(k, j);
            N(i, j) = s;
        }
        T s = num_traits<T>::from_int(0);
        for (std::size_t k = 0; k < m; ++k) s += M(k, i) * b[k];
        r[i] = s;
    }
    return solve(N, r);
}

}  // namespace mfq_detail

/**
 * Psi, K and U of a fluid queue whose drifts have been normalized to +-1.
 *
 * @param Fpp up-to-up block, Fpm up-to-down, Fmp down-to-up, Fmm down-to-down
 * @param precision stopping tolerance on the doubling residual
 * @param maxNumIt  iteration cap
 * @param method    ADDA (the BUTools default) or SDA
 * @param Fpm generator block from the up-phases to the down-phases
 * @param Fmp generator block from the down-phases to the up-phases
 * @param Fmm generator block within the down-phases
 */
template <class T>
FluidFundamental<T> mfq_fundamental(const Matrix<T>& Fpp, const Matrix<T>& Fpm,
                                    const Matrix<T>& Fmp, const Matrix<T>& Fmm,
                                    const T& precision, unsigned maxNumIt,
                                    RiccatiMethod method) {
    static_assert(num_traits<T>::has_transcendental,
                  "mfq_fundamental runs a tolerance-terminated doubling iteration");
    using namespace mfq_detail;
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t sA = Fpp.rows();
    const std::size_t sD = Fmm.rows();
    if (Fpp.cols() != sA || Fmm.cols() != sD || Fpm.rows() != sA || Fpm.cols() != sD ||
        Fmp.rows() != sD || Fmp.cols() != sA)
        throw InputError("mfq_fundamental: the four blocks are not conformable");

    FluidFundamental<T> out;
    out.iterations = 0;
    out.converged = true;
    if (sA == 0) {
        out.Psi = Matrix<T>(0, sD, zero);
        out.K = Matrix<T>(0, 0, zero);
        out.U = Fmm;
        return out;
    }

    // ADDA / SDA (Wang, Wang, Li 2011).
    Matrix<T> A = scale(Fpp, T(-num_traits<T>::from_int(1)));
    const Matrix<T>& B = Fpm;
    const Matrix<T>& C = Fmp;
    Matrix<T> D = scale(Fmm, T(-num_traits<T>::from_int(1)));
    T gamma1 = A(0, 0), gamma2 = (sD > 0) ? D(0, 0) : zero;
    for (std::size_t i = 0; i < sA; ++i)
        if (A(i, i) > gamma1) gamma1 = A(i, i);
    for (std::size_t i = 0; i < sD; ++i)
        if (D(i, i) > gamma2) gamma2 = D(i, i);
    if (method == RiccatiMethod::SDA) {
        if (gamma2 > gamma1) gamma1 = gamma2;
        gamma2 = gamma1;
    }
    const Matrix<T> IA = eye<T>(sA);
    const Matrix<T> ID = eye<T>(sD);
    for (std::size_t i = 0; i < sA; ++i) A(i, i) += gamma2;
    for (std::size_t i = 0; i < sD; ++i) D(i, i) += gamma1;
    const T g = gamma1 + gamma2;

    const Matrix<T> Dginv0 = inverse(D);
    Matrix<T> Vginv = inverse(sub(D, matmul(matmul(C, inverse(A)), B)));
    Matrix<T> Wginv = inverse(sub(A, matmul(matmul(B, Dginv0), C)));
    Matrix<T> Eg = sub(ID, scale(Vginv, g));
    Matrix<T> Fg = sub(IA, scale(Wginv, g));
    Matrix<T> Gg = scale(matmul(matmul(Dginv0, C), Wginv), g);
    Matrix<T> Hg = scale(matmul(matmul(Wginv, B), Dginv0), g);

    T diff = num_traits<T>::from_int(1);
    unsigned numit = 0;
    while (diff > precision && numit < maxNumIt) {
        Vginv = matmul(Eg, inverse(sub(ID, matmul(Gg, Hg))));
        Wginv = matmul(Fg, inverse(sub(IA, matmul(Hg, Gg))));
        Gg = add(Gg, matmul(matmul(Vginv, Gg), Fg));
        Hg = add(Hg, matmul(matmul(Wginv, Hg), Eg));
        Eg = matmul(Vginv, Eg);
        Fg = matmul(Wginv, Fg);
        const T neg = norm1(Eg);
        const T nfg = norm1(Fg);
        if (method == RiccatiMethod::ADDA) {
            const T eta = sqrt(nfg / neg);
            Eg = scale(Eg, eta);
            Fg = scale(Fg, T(num_traits<T>::from_int(1) / eta));
            diff = neg * nfg;
        } else {
            diff = (neg < nfg) ? neg : nfg;
        }
        ++numit;
    }
    out.iterations = numit;
    out.converged = (numit < maxNumIt);
    out.Psi = Hg;
    out.K = add(Fpp, matmul(out.Psi, Fmp));
    out.U = add(Fmm, matmul(Fmp, out.Psi));
    return out;
}

/** Stationary matrix-exponential solution of a general Markovian fluid model. */
template <class T>
struct GeneralFluidSolution {
    std::vector<T> mass0;  ///< P(level 0, state j), length N
    std::vector<T> ini;    ///< initial vector of the density, length Np
    Matrix<T> K;           ///< matrix exponent of the density, Np x Np
    Matrix<T> clo;         ///< closing matrix of the density, Np x N
};

/**
 * Stationary law of a general Markovian fluid model, pi(x) = ini exp(K x) clo
 * above level zero plus the point mass mass0 at zero.
 *
 * @param Q  generator of the background chain, N x N
 * @param R  diagonal drift matrix, N x N, entries of any sign
 * @param Q0 boundary generator at level zero; pass an empty matrix for the
 *           regular boundary behaviour Q0 = Q
 * @param prec tolerance, used both to classify a drift as zero and to stop the
 *             Riccati iteration, exactly as in the reference
 */
template <class T>
GeneralFluidSolution<T> mfq_general_solve(const Matrix<T>& Q, const Matrix<T>& R,
                                          const Matrix<T>& Q0, const T& prec) {
    static_assert(num_traits<T>::has_transcendental,
                  "mfq_general_solve calls a tolerance-terminated Riccati solver");
    using namespace mfq_detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t N = Q.rows();
    if (Q.cols() != N || R.rows() != N || R.cols() != N)
        throw InputError("mfq_general_solve: Q and R must be square and of equal order");

    // Partition the state space by the sign of the drift: zero, positive,
    // negative, in that order.
    std::vector<std::size_t> ixz, ixp, ixn;
    for (std::size_t i = 0; i < N; ++i) {
        const T d = R(i, i);
        if (num_abs(d) <= prec)
            ixz.push_back(i);
        else if (d > prec)
            ixp.push_back(i);
        else
            ixn.push_back(i);
    }
    const std::size_t Nz = ixz.size(), Np = ixp.size(), Nn = ixn.size();
    // Nn==0 vs Np==0 rationale: see _kb/03-api-layer.md (cpp port notes: mam) and _kb/14-cpp-multiprecision.md
    if (Nn == 0)
        throw InputError(
            "mfq_general_solve: every state has an up drift, so the fluid level grows without "
            "bound and has no stationary law");

    std::vector<std::size_t> perm;
    perm.insert(perm.end(), ixz.begin(), ixz.end());
    perm.insert(perm.end(), ixp.begin(), ixp.end());
    perm.insert(perm.end(), ixn.begin(), ixn.end());
    Matrix<T> P(N, N, zero);
    for (std::size_t i = 0; i < N; ++i) P(i, perm[i]) = one;

    // Qv = P Q P^T, Rv = P R P^T: P is a permutation, so its inverse is P^T.
    Matrix<T> Qv(N, N, zero), Rv(N, N, zero);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) {
            Qv(i, j) = Q(perm[i], perm[j]);
            Rv(i, j) = R(perm[i], perm[j]);
        }

    // Censor the zero-drift states out.
    Matrix<T> iQv00(Nz, Nz, zero);
    if (Nz > 0) {
        Matrix<T> negQ00(Nz, Nz, zero);
        for (std::size_t i = 0; i < Nz; ++i)
            for (std::size_t j = 0; j < Nz; ++j) negQ00(i, j) = -Qv(i, j);
        try {
            iQv00 = inverse(negQ00);
        } catch (const NumericError&) {
            throw NumericError(
                "mfq_general_solve: the zero-drift states form a closed set, so the fluid can be "
                "trapped at a constant level and the model has no stationary fluid law; the "
                "MATLAB reference silently substitutes a pseudo-inverse here");
        }
    }
    const std::size_t Npn = Np + Nn;
    Matrix<T> Qbar = block(Qv, Nz, Nz, Npn, Npn);
    if (Nz > 0) {
        const Matrix<T> L = block(Qv, Nz, 0, Npn, Nz);
        const Matrix<T> Rt = block(Qv, 0, Nz, Nz, Npn);
        Qbar = add(Qbar, matmul(matmul(L, iQv00), Rt));
    }
    Matrix<T> absRi(Npn, Npn, zero);
    for (std::size_t i = 0; i < Npn; ++i) absRi(i, i) = one / num_abs(Rv(Nz + i, Nz + i));
    const Matrix<T> Qz = matmul(absRi, Qbar);

    const FluidFundamental<T> ff =
        mfq_fundamental(block(Qz, 0, 0, Np, Np), block(Qz, 0, Np, Np, Nn),
                        block(Qz, Np, 0, Nn, Np), block(Qz, Np, Np, Nn, Nn), prec, 150u,
                        RiccatiMethod::ADDA);
    const Matrix<T>& Psi = ff.Psi;
    const Matrix<T>& K = ff.K;
    const Matrix<T>& U = ff.U;

    // Pm = [I_Np, Psi], iCn and iCp the down and up blocks of absRi.
    Matrix<T> Pm(Np, Npn, zero);
    for (std::size_t i = 0; i < Np; ++i) {
        Pm(i, i) = one;
        for (std::size_t j = 0; j < Nn; ++j) Pm(i, Np + j) = Psi(i, j);
    }
    const Matrix<T> iCp = block(absRi, 0, 0, Np, Np);
    const Matrix<T> iCn = block(absRi, Np, Np, Nn, Nn);

    // clo = [ (iCp Qv(p,z) + Psi iCn Qv(n,z)) iQv00 , Pm absRi ], Np x N.
    Matrix<T> clo(Np, N, zero);
    if (Nz > 0) {
        const Matrix<T> Qpz = block(Qv, Nz, 0, Np, Nz);
        const Matrix<T> Qnz = block(Qv, Nz + Np, 0, Nn, Nz);
        const Matrix<T> lhs =
            matmul(add(matmul(iCp, Qpz), matmul(matmul(Psi, iCn), Qnz)), iQv00);
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t j = 0; j < Nz; ++j) clo(i, j) = lhs(i, j);
    }
    {
        const Matrix<T> rhs = matmul(Pm, absRi);
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t j = 0; j < Npn; ++j) clo(i, Nz + j) = rhs(i, j);
    }

    GeneralFluidSolution<T> out;
    out.K = K;
    const std::vector<T> eN = ones<T>(N);
    const Matrix<T> negKinv = inverse(scale(K, T(-one)));

    if (Q0.rows() == 0) {
        // Regular boundary behaviour, Q0 = Q.
        // clo goes back to the original state ordering: clo * P.
        Matrix<T> cloP(Np, N, zero);
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t j = 0; j < N; ++j) cloP(i, perm[j]) = clo(i, j);

        std::vector<T> Ua(Nn, zero);
        if (Nz > 0) {
            const Matrix<T> Qnz = block(Qv, Nz + Np, 0, Nn, Nz);
            const std::vector<T> t = mulvec(matmul(matmul(iCn, Qnz), iQv00), ones<T>(Nz));
            for (std::size_t i = 0; i < Nn; ++i) Ua[i] += t[i];
        }
        {
            const std::vector<T> t = mulvec(iCn, ones<T>(Nn));
            for (std::size_t i = 0; i < Nn; ++i) Ua[i] += t[i];
        }
        {
            const Matrix<T> Qnp = block(Qz, Np, 0, Nn, Np);
            const std::vector<T> t = mulvec(matmul(matmul(Qnp, negKinv), cloP), eN);
            for (std::size_t i = 0; i < Nn; ++i) Ua[i] += t[i];
        }
        // pm solves pm [U, Ua] = [0 ... 0, 1], an (Nn+1)-equation system in Nn
        // unknowns that is consistent by construction.
        Matrix<T> Msys(Nn + 1, Nn, zero);
        for (std::size_t j = 0; j < Nn; ++j)
            for (std::size_t i = 0; i < Nn; ++i) Msys(j, i) = U(i, j);
        for (std::size_t i = 0; i < Nn; ++i) Msys(Nn, i) = Ua[i];
        std::vector<T> rhs(Nn + 1, zero);
        rhs[Nn] = one;
        const std::vector<T> pm = normal_equations_solve(Msys, rhs);

        std::vector<T> m0(N, zero);
        {
            const std::vector<T> pmiCn = vecmul(pm, iCn);
            if (Nz > 0) {
                const Matrix<T> Qnz = block(Qv, Nz + Np, 0, Nn, Nz);
                const std::vector<T> t = vecmul(vecmul(pmiCn, Qnz), iQv00);
                for (std::size_t j = 0; j < Nz; ++j) m0[j] = t[j];
            }
            for (std::size_t j = 0; j < Nn; ++j) m0[Nz + Np + j] = pmiCn[j];
        }
        out.mass0.assign(N, zero);
        for (std::size_t j = 0; j < N; ++j) out.mass0[perm[j]] = m0[j];
        out.ini = vecmul(pm, block(Qz, Np, 0, Nn, Np));
        out.clo = cloP;
    } else {
        if (Q0.rows() != N || Q0.cols() != N)
            throw InputError("mfq_general_solve: Q0 must be square and of the order of Q");
        Matrix<T> Q0v(N, N, zero);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) Q0v(i, j) = Q0(perm[i], perm[j]);

        // M = [-clo Rv ; Q0v(n, :) ; Q0v(z, :)], N x N, and Ma the normalizer.
        Matrix<T> Msys(N + 1, N, zero);
        const Matrix<T> cloRv = matmul(clo, Rv);
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t j = 0; j < N; ++j) Msys(j, i) = -cloRv(i, j);
        for (std::size_t i = 0; i < Nn; ++i)
            for (std::size_t j = 0; j < N; ++j) Msys(j, Np + i) = Q0v(Nz + Np + i, j);
        for (std::size_t i = 0; i < Nz; ++i)
            for (std::size_t j = 0; j < N; ++j) Msys(j, Np + Nn + i) = Q0v(i, j);
        {
            const std::vector<T> s = mulvec(matmul(negKinv, clo), eN);
            for (std::size_t i = 0; i < Np; ++i) Msys(N, i) = s[i];
            for (std::size_t i = Np; i < N; ++i) Msys(N, i) = one;
        }
        std::vector<T> rhs(N + 1, zero);
        rhs[N] = one;
        const std::vector<T> sol = normal_equations_solve(Msys, rhs);

        out.ini.assign(sol.begin(), sol.begin() + Np);
        Matrix<T> cloP(Np, N, zero);
        for (std::size_t i = 0; i < Np; ++i)
            for (std::size_t j = 0; j < N; ++j) cloP(i, perm[j]) = clo(i, j);
        out.clo = cloP;
        std::vector<T> m0(N, zero);
        for (std::size_t j = 0; j < Nz; ++j) m0[j] = sol[Np + Nn + j];
        for (std::size_t j = 0; j < Nn; ++j) m0[Nz + Np + j] = sol[Np + j];
        out.mass0.assign(N, zero);
        for (std::size_t j = 0; j < N; ++j) out.mass0[perm[j]] = m0[j];
    }
    return out;
}

/** mfq_general_solve with the regular boundary and the BUTools default prec = 1e-14. */
template <class T>
GeneralFluidSolution<T> mfq_general_solve(const Matrix<T>& Q, const Matrix<T>& R) {
    return mfq_general_solve(Q, R, Matrix<T>(0, 0), T(num_traits<T>::from_double(1e-14)));
}

/**
 * The similarity transformation B with B v = e, for a non-negative column
 * vector v. Port of BUTools TransformToOnes: sort v decreasing so that a
 * non-zero entry leads, then take the lower-triangular matrix of reciprocal
 * partial sums. It works even when v has zero entries, which is why the naive
 * diag(1/v) is not used.
 */
template <class T>
Matrix<T> mfq_transform_to_ones(const std::vector<T>& v) {
    const std::size_t m = v.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (m == 0) throw InputError("mfq_transform_to_ones: empty vector");
    std::vector<std::size_t> ix(m);
    for (std::size_t i = 0; i < m; ++i) ix[i] = i;
    // Stable sort on -v, i.e. decreasing v, matching MATLAB's sort(-clovec).
    for (std::size_t i = 1; i < m; ++i)
        for (std::size_t j = i; j > 0 && v[ix[j]] > v[ix[j - 1]]; --j) std::swap(ix[j], ix[j - 1]);
    std::vector<T> cp(m);
    for (std::size_t i = 0; i < m; ++i) cp[i] = v[ix[i]];
    Matrix<T> Bt(m, m, zero);
    T acc = zero;
    for (std::size_t i = 0; i < m; ++i) {
        acc += cp[i];
        if (acc == zero)
            throw NumericError("mfq_transform_to_ones: the closing vector sums to zero");
        for (std::size_t j = 0; j <= i; ++j) Bt(i, j) = one / acc;
    }
    // B = Bt * P, with P(i, ix(i)) = 1.
    Matrix<T> B(m, m, zero);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) B(i, ix[j]) = Bt(i, j);
    return B;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MFQ_SOLVE_H
