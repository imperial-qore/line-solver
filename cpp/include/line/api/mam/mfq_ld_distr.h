/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MFQ_LD_DISTR_H
#define LINE_API_MAM_MFQ_LD_DISTR_H

/**
 * Stationary density and distribution of a first- or second-order
 * level-dependent (multi-regime) Markovian fluid queue, evaluated at requested
 * fluid levels.
 *
 * Port of matlab/src/api/mam/mfq_ld_distr.m and the BUTools
 * LevelDependentFluidStationaryDistr it wraps. Its inputs are the
 * matrix-exponential building blocks that mfq_ld_solve returns, the same
 * LevelDependentFluidBlocks that mfq_ld_mean consumes, so the two are
 * interchangeable consumers of one solve.
 *
 * THE DENSITY. Over regime k, that is over the level interval
 * [T(k), T(k+1)], the stationary density is anchored at both ends,
 *
 *     pi_k(x) = iniF_k exp(KF_k (x - T(k))) cloF_k
 *             + iniB_k exp(KB_k (T(k+1) - x)) cloB_k,
 *
 * with point masses at the K+1 thresholds. Four quantities are offered:
 *   Pdf  the per-state density
 *   Pdfd its derivative
 *   Cdf  P(X < p), so a point mass sitting exactly at p is EXCLUDED
 *   Cdfm P(X <= p), so that mass is INCLUDED
 * The Cdf/Cdfm distinction is not cosmetic here: a level-dependent fluid queue
 * puts genuine atoms at its thresholds, so the two differ by a finite amount at
 * every threshold and by nothing anywhere else.
 *
 * INTEGRATING A SINGULAR EXPONENT. The cumulative forms need
 * int_0^L exp(M u) du. When M is non-singular that is
 * (-M)^-1 (I - exp(M L)); when M has a zero eigenvalue -- which is the normal
 * case for one of the two directions, not an edge case -- the inverse does not
 * exist and the reference deflates it instead. With l and r the left and right
 * null vectors of M normalized so that l r = 1,
 *
 *     int_0^L exp(M u) du = (-(M - r l))^-1 (I - exp((M - r l) L))
 *                         + r l (L + exp(-L) - 1),
 *
 * the rank-one shift moving the zero eigenvalue to -1 and its contribution
 * being added back in closed form. The null vectors come from BUTools CRPSolve,
 * which is a LINEAR SOLVE and not an eigenvector computation: it replaces the
 * first column of M by ones and solves, which is legitimate because M has zero
 * row sums by construction. That is worth stating because it is the reason this
 * function does NOT need an eigendecomposition for its arithmetic.
 *
 * WHERE THE EIGENVALUES DO ENTER, AND WHY Real IS STILL HONEST. The reference
 * chooses which of KF, KB to deflate by comparing min |eig(KF)| against
 * min |eig(KB)|, i.e. by asking which is closer to singular. That comparison is
 * the ONLY use of an eigendecomposition in the whole function, and its result
 * is a BRANCH SELECTION, a discrete choice between two formulas for the same
 * integral. No eigenvalue flows into any returned number. The port therefore
 * converts the two matrices to double for the comparison alone, exactly as
 * util/eig.h instructs its callers to do, and carries out the integral itself
 * in T. At Real50 the returned values are genuinely Real50-accurate; what is
 * computed in double is which of two algebraically equivalent routes to take.
 * If the build has no LAPACK, eig_values refuses by name, and this function
 * refuses with it rather than guessing the branch.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental: expm is a
 * tolerance-controlled Pade approximation in any arithmetic, and the deflation
 * formula evaluates exp(-L). See the paragraph above for why the LAPACK
 * dependency does not reduce this to double.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "line/api/mam/mfq_ld_mean.h"
#include "line/num/number.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Which functional of the stationary level law to evaluate. */
enum class FluidDistrKind {
    Pdf,   ///< per-state density
    Pdfd,  ///< derivative of the density
    Cdf,   ///< P(X < p), excluding an atom sitting exactly at p
    Cdfm   ///< P(X <= p), including it
};

namespace mfq_ld_detail {

/**
 * BUTools CRPSolve: the stationary vector of a continuous-time RATIONAL
 * process, pi M = 0 with sum(pi) = 1. M needs zero row sums but, unlike a
 * generator, is not required to have non-negative off-diagonals. Replacing the
 * first column by ones turns the singular system into a non-singular one whose
 * unique solution is the normalized null vector, so this is a linear solve and
 * NOT an eigenvector computation.
 */
template <class T>
std::vector<T> crp_solve(const Matrix<T>& M) {
    const std::size_t n = M.rows();
    if (M.cols() != n) throw InputError("crp_solve: matrix is not square");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    // pi M' = m with M' = M, first column replaced by ones, and m = e_1.
    // Transposed for the row-vector solve.
    Matrix<T> A(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(j, i) = (j == 0) ? one : M(i, j);
    std::vector<T> rhs(n, zero);
    rhs[0] = one;
    return solve(A, rhs);
}

/** Smallest eigenvalue modulus of M, computed in double: a branch test only. */
template <class T>
double min_abs_eig(const Matrix<T>& M) {
    const std::size_t n = M.rows();
    Matrix<double> Md(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Md(i, j) = num_traits<T>::to_double(M(i, j));
    const std::vector<std::complex<double>> ev = eig_values(Md);
    double best = std::abs(ev[0]);
    for (const std::complex<double>& v : ev) {
        const double a = std::abs(v);
        if (a < best) best = a;
    }
    return best;
}

/**
 * int_0^L exp(M u) du for a SINGULAR M, by the rank-one deflation of the
 * reference: shift the zero eigenvalue to -1 with M - r l, integrate the
 * non-singular shift, and add the deflated direction's contribution
 * r l (L + exp(-L) - 1) back in closed form.
 */
template <class T>
Matrix<T> integ_exp_singular(const Matrix<T>& M, const T& L) {
    using std::exp;
    const std::size_t n = M.rows();
    const T one = num_traits<T>::from_int(1);
    Matrix<T> Mt(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Mt(i, j) = M(j, i);
    std::vector<T> l = crp_solve(M);   // left null vector
    std::vector<T> r = crp_solve(Mt);  // right null vector, as a row of M'
    T lr = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) lr += l[i] * r[i];
    if (lr == num_traits<T>::from_int(0))
        throw NumericError("mfq_ld_distr: the null vectors of a regime exponent are orthogonal");
    for (T& v : l) v /= lr;

    Matrix<T> rl(n, n);  // the rank-one r l
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) rl(i, j) = r[i] * l[j];
    Matrix<T> S = M;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) S(i, j) -= rl(i, j);

    Matrix<T> negS(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negS(i, j) = -S(i, j);
    const Matrix<T> E = expm(S, L);
    Matrix<T> ImE = eye<T>(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) ImE(i, j) -= E(i, j);
    Matrix<T> out = matmul(inverse(negS), ImE);
    const T w = L + exp(-L) - one;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) out(i, j) += rl(i, j) * w;
    return out;
}

/** int_0^L exp(M u) du for a NON-singular M. */
template <class T>
Matrix<T> integ_exp_regular(const Matrix<T>& M, const T& L) {
    const std::size_t n = M.rows();
    Matrix<T> negM(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negM(i, j) = -M(i, j);
    const Matrix<T> E = expm(M, L);
    Matrix<T> ImE = eye<T>(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) ImE(i, j) -= E(i, j);
    return matmul(inverse(negM), ImE);
}

/**
 * The pair of integrals for one regime. The direction whose exponent is closer
 * to singular is deflated, the other inverted directly; see the header note on
 * why deciding that in double costs the result nothing.
 */
template <class T>
void integ_exp_pair(const Matrix<T>& KA, const Matrix<T>& KB, const T& L, Matrix<T>& JA,
                    Matrix<T>& JB) {
    if (min_abs_eig(KA) > min_abs_eig(KB)) {
        JA = integ_exp_regular(KA, L);
        JB = integ_exp_singular(KB, L);
    } else {
        JA = integ_exp_singular(KA, L);
        JB = integ_exp_regular(KB, L);
    }
}

}  // namespace mfq_ld_detail

/**
 * Stationary density or distribution of a level-dependent fluid queue.
 *
 * @param b      the building blocks returned by mfq_ld_solve
 * @param what   which functional to evaluate
 * @param points the fluid levels at which to evaluate it
 * @return       one row of N per-state values per requested point
 */
template <class T>
std::vector<std::vector<T>> mfq_ld_distr(const LevelDependentFluidBlocks<T>& b,
                                         FluidDistrKind what, const std::vector<T>& points) {
    static_assert(num_traits<T>::has_transcendental,
                  "mfq_ld_distr evaluates matrix exponentials");
    using namespace mfq_ld_detail;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = b.Thr.size();
    if (K == 0) throw InputError("mfq_ld_distr: at least one regime is required");
    if (b.masses.size() != K + 1)
        throw InputError("mfq_ld_distr: expected K+1 point-mass vectors");
    if (b.iniF.size() != K || b.KF.size() != K || b.cloF.size() != K || b.iniB.size() != K ||
        b.KB.size() != K || b.cloB.size() != K)
        throw InputError("mfq_ld_distr: expected K forward and K backward blocks");
    const std::size_t N = b.masses[0].size();

    std::vector<T> Tv(K + 1, zero);
    for (std::size_t k = 0; k < K; ++k) Tv[k + 1] = b.Thr[k];
    const bool cumulate = (what == FluidDistrKind::Cdf || what == FluidDistrKind::Cdfm);

    std::vector<std::vector<T>> res;
    res.reserve(points.size());
    for (const T& p : points) {
        if (p < zero) throw InputError("mfq_ld_distr: the evaluation points must be non-negative");
        std::vector<T> pres(N, zero);
        // cumulative-form threshold walk rationale: see _kb/03-api-layer.md (cpp port notes: mam)
        std::size_t k = 0;
        while (k < K && p >= Tv[k]) {
            if (cumulate) {
                if (k > 0) {
                    Matrix<T> sF, sB;
                    integ_exp_pair(b.KF[k - 1], b.KB[k - 1], T(Tv[k] - Tv[k - 1]), sF, sB);
                    const std::vector<T> vF =
                        vecmul(vecmul(b.iniF[k - 1], sF), b.cloF[k - 1]);
                    const std::vector<T> vB =
                        vecmul(vecmul(b.iniB[k - 1], sB), b.cloB[k - 1]);
                    for (std::size_t j = 0; j < N; ++j) pres[j] += vF[j] + vB[j];
                }
                // Cdf excludes an atom sitting exactly at p, Cdfm includes it.
                if (p > Tv[k] || what == FluidDistrKind::Cdfm)
                    for (std::size_t j = 0; j < N; ++j) pres[j] += b.masses[k][j];
            }
            ++k;
        }
        if (k == K && p == Tv[K] && what == FluidDistrKind::Cdfm)
            for (std::size_t j = 0; j < N; ++j) pres[j] += b.masses[K][j];
        // The regime holding p is the one below threshold k.
        const std::size_t kk = k - 1;
        const T prem = p - Tv[kk];
        const T Tk = Tv[kk + 1] - Tv[kk];

        if (what == FluidDistrKind::Pdf) {
            const std::vector<T> vF =
                vecmul(vecmul(b.iniF[kk], expm(b.KF[kk], prem)), b.cloF[kk]);
            const std::vector<T> vB =
                vecmul(vecmul(b.iniB[kk], expm(b.KB[kk], T(Tk - prem))), b.cloB[kk]);
            for (std::size_t j = 0; j < N; ++j) pres[j] = vF[j] + vB[j];
        } else if (what == FluidDistrKind::Pdfd) {
            const std::vector<T> vF = vecmul(
                vecmul(vecmul(b.iniF[kk], b.KF[kk]), expm(b.KF[kk], prem)), b.cloF[kk]);
            const std::vector<T> vB =
                vecmul(vecmul(vecmul(b.iniB[kk], b.KB[kk]), expm(b.KB[kk], T(Tk - prem))),
                       b.cloB[kk]);
            for (std::size_t j = 0; j < N; ++j) pres[j] = vF[j] - vB[j];
        } else {
            Matrix<T> sF, sB;
            integ_exp_pair(b.KF[kk], b.KB[kk], prem, sF, sB);
            const std::vector<T> vF = vecmul(vecmul(b.iniF[kk], sF), b.cloF[kk]);
            const std::vector<T> vB = vecmul(
                vecmul(vecmul(b.iniB[kk], expm(b.KB[kk], T(Tk - prem))), sB), b.cloB[kk]);
            for (std::size_t j = 0; j < N; ++j) pres[j] += vF[j] + vB[j];
        }
        res.push_back(pres);
    }
    return res;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MFQ_LD_DISTR_H
