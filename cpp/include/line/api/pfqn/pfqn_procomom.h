/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PROCOMOM_H
#define LINE_API_PFQN_PROCOMOM_H

/**
 * ProCoMoM: marginal queue-length probabilities of a closed multiclass
 * product-form network by the class-oriented method of moments.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_procomom.m and
 * matlab/src/api/pfqn/pfqn_procomom2.m.
 *
 * pfqn_procomom carries, instead of the plain moments of pfqn_comom, the
 * generating coefficients of ONE station's marginal distribution: pk(:, j+1)
 * holds the basis coefficients of P(n_station = j). The basis is the same
 * Dn = multichoose(R,M) layout as pfqn_comom, with the reference's UP shifts
 * (Dn + e_s) where pfqn_comom uses down shifts, and the station of interest is
 * rotated into the last row of L so the two extra couplings, DC (same level,
 * one job less at the station) and DD (previous level), always refer to row M.
 * Every station is solved in turn and its distribution normalized to sum 1.
 *
 * pfqn_procomom2 is the two-node special case, a queue plus a delay, where the
 * whole recursion collapses to a product of bidiagonal transfer matrices,
 *
 *   F = prod_r T_r^{N_r} / N_r!,   T_r(row,row) = Z_r,
 *                                  T_r(row,row+1) = (n + m - 1) L_r / mu(n),
 *
 * and the unnormalized marginal is F e_last. No linear system is solved at all.
 *
 * SYSTEM SHAPE. The ProCoMoM matrix is NOT square: countrows(r) contributes
 * M rows per basis vector in the propagation branch but M + r - 1 (or 1) in the
 * equation branch, so for M = R = 2 the class-1 step is 6 x 6 and the class-2
 * step is 7 x 6. The reference solves it with an economy QR when full rank and
 * with a truncated SVD when not. Both are done here by util/lstsq.h in exact
 * arithmetic: the normal equations give the identical least-squares vector, and
 * a full-rank factorization gives the identical minimum-norm vector.
 *
 * NOT PORTED, deliberately: the reference's rank-deficiency remedy, which
 * re-solves the whole problem with the demands randomly perturbed by
 * 1e-10..1e-4 of their scale (`rng(23000,'twister')`) and keeps whichever
 * attempt looks best behaved. It is a floating-point workaround for a
 * floating-point rank test, it changes the answer, and reproducing it would
 * require MATLAB's Mersenne-Twister stream bit for bit. This port returns the
 * exact pseudoinverse solution, which is what the perturbation approximates,
 * and reports the rank deficiency in `rankdef` instead of hiding it.
 *
 * Arithmetic: EXACT-CAPABLE. Field operations and least-squares solves only.
 *
 * REFERENCE DEFECTS in pfqn_procomom2.m, all reproducible, none corrected in
 * MATLAB by this port (it does not edit MATLAB):
 *
 *  1. `pfqn_procomom2(L,N,Z)` raises "Not enough input arguments". The
 *     nargin < 4 branch evaluates `ones(m, sum(N)+1)`, but `m` is assigned only
 *     by the nargin < 5 branch BELOW it. Reproduce with
 *     pfqn_procomom2([0.4 0.6],[2 1],[1 2]).
 *  2. Debug scaffolding runs on EVERY call: `QN=double(Q)` and a full
 *     `pfqn_mvald(...)` are executed without semicolons, so the routine prints
 *     QN, XNMVA, QNMVA and pik to the console and pays for an O(prod(N+1))
 *     load-dependent MVA it does not use. Visible in the output of the call
 *     above.
 *  3. Dead branch: `if any(~isfinite(pk(1,1))) || ~isfinite(G) ... elseif
 *     ~isfinite(G)` -- the elseif is subsumed by the first test and can never
 *     be taken.
 *
 * This port implements the transfer-matrix method only, with no printing and
 * no MVA call, and accepts the three-argument form.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_comb_common.h"
#include "line/api/pfqn/pfqn_comom.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lstsq.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_procomom, mirroring [Pr, Q]. */
template <class T>
struct ProcomomResult {
    Matrix<T> Pr;       ///< (M x sumN+1) marginals, Pr(k, j) = P(n_k = j)
    std::vector<T> Q;   ///< (M) mean queue lengths
    bool rankdef;       ///< true when some ProCoMoM system was rank deficient
};

/**
 * Marginal queue-length distributions of every station.
 *
 * @param L    (M x R) demand matrix
 * @param N    (R) populations
 * @param Z    (R) think times
 * @param atol tolerance below which a class demand counts as zero
 */
template <class T>
ProcomomResult<T> pfqn_procomom(const Matrix<T>& L, const std::vector<int>& N,
                                const std::vector<T>& Z, const T& atol) {
    const std::size_t M = L.rows(), R = L.cols();
    if (M == 0 || R == 0) throw InputError("pfqn_procomom: empty demand matrix");
    if (N.size() != R) throw InputError("pfqn_procomom: L and N disagree on the class count");
    if (Z.size() != R) throw InputError("pfqn_procomom: L and Z disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    int sumN = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < 0) throw InputError("pfqn_procomom: negative population");
        sumN += N[r];
    }
    const std::size_t NP = static_cast<std::size_t>(sumN) + 1;

    // ---- per-class rescaling; normalized marginals are invariant to it -------------
    Matrix<T> Ls(M, R, zero);
    std::vector<T> Zs(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        T mx = L(0, r);
        for (std::size_t i = 1; i < M; ++i)
            if (L(i, r) > mx) mx = L(i, r);
        if (mx < atol) mx = one;
        for (std::size_t i = 0; i < M; ++i) Ls(i, r) = L(i, r) / mx;
        Zs[r] = Z[r] / mx;
    }

    const std::vector<std::vector<int>> Dn =
        detail::comom_basis(static_cast<int>(R), static_cast<int>(M));
    const std::size_t numDn = Dn.size();
    const std::size_t basisSize = numDn * M;
    const std::vector<int> zeroDn(R, 0);

    // zero-based basis column of (dn, i) with i = 1..M
    const auto phash = [&](const std::vector<int>& dn, std::size_t i) -> long {
        const int pos = matchrow(Dn, dn);
        if (pos < 0) return -1;
        return static_cast<long>(pos) * static_cast<long>(M) + static_cast<long>(i) - 1;
    };

    ProcomomResult<T> res;
    res.Pr = Matrix<T>(M, NP, zero);
    res.Q.assign(M, zero);
    res.rankdef = false;

    for (std::size_t station = 0; station < M; ++station) {
        // rotate the station of interest into the last row
        Matrix<T> Lr = Ls;
        for (std::size_t r = 0; r < R; ++r) {
            const T tmp = Lr(station, r);
            Lr(station, r) = Lr(M - 1, r);
            Lr(M - 1, r) = tmp;
        }

        Matrix<T> pk(basisSize, NP, zero);
        for (std::size_t kk = 1; kk <= M; ++kk) pk(static_cast<std::size_t>(phash(zeroDn, kk)), 0) = one;

        std::vector<int> Ncur(R, 0);
        for (std::size_t r = 0; r < R; ++r) {
            for (int Nr = 1; Nr <= N[r]; ++Nr) {
                Ncur[r] = Nr;
                const Matrix<T> pklast = pk;
                pk = Matrix<T>(basisSize, NP, zero);

                // ---- genpmatrix ------------------------------------------------
                std::size_t numRows = 0;
                for (std::size_t d = 0; d < numDn; ++d) {
                    int s1 = 0;
                    if (r + 1 <= R - 1)
                        for (std::size_t c = r; c + 1 < R; ++c) s1 += Dn[d][c];
                    if (r + 1 <= R - 1 && s1 > 0) {
                        numRows += M;
                    } else {
                        int s3 = 0;
                        for (std::size_t c = 0; c <= r; ++c) s3 += Dn[d][c];
                        numRows += (s3 < static_cast<int>(M)) ? (M + r) : 1;
                    }
                }
                Matrix<T> Ag(numRows, basisSize, zero), Bg(numRows, basisSize, zero),
                    DCg(numRows, basisSize, zero), DDg(numRows, basisSize, zero);
                std::size_t row = 0;
                for (std::size_t d = 0; d < numDn; ++d) {
                    const std::vector<int>& dv = Dn[d];
                    int s1 = 0;
                    if (r + 1 <= R - 1)
                        for (std::size_t c = r; c + 1 < R; ++c) s1 += dv[c];
                    if (r + 1 <= R - 1 && s1 > 0) {
                        int s2 = 0;
                        for (std::size_t c = r + 1; c + 1 < R; ++c) s2 += dv[c];
                        for (std::size_t k = 1; k <= M; ++k) {
                            const long cA = phash(dv, k);
                            Ag(row, static_cast<std::size_t>(cA)) = one;
                            if (s2 > 0) {
                                Bg(row, static_cast<std::size_t>(cA)) = one;
                            } else {
                                std::vector<int> sh = dv;
                                sh[r] -= 1;
                                const long cB = phash(sh, k);
                                if (cB >= 0) Bg(row, static_cast<std::size_t>(cB)) = one;
                            }
                            ++row;
                        }
                    } else {
                        int s3 = 0;
                        for (std::size_t c = 0; c <= r; ++c) s3 += dv[c];
                        if (s3 < static_cast<int>(M)) {
                            for (std::size_t k = 1; k + 1 <= M; ++k) {
                                Ag(row, static_cast<std::size_t>(phash(dv, k + 1))) = one;
                                Ag(row, static_cast<std::size_t>(phash(dv, 1))) = -one;
                                for (std::size_t s = 0; s < r; ++s) {
                                    std::vector<int> sh = dv;
                                    sh[s] += 1;  // UP shift
                                    const long c = phash(sh, k + 1);
                                    if (c >= 0) Ag(row, static_cast<std::size_t>(c)) = -Lr(k - 1, s);
                                }
                                Bg(row, static_cast<std::size_t>(phash(dv, k + 1))) = Lr(k - 1, r);
                                ++row;
                            }
                            for (std::size_t s = 0; s < r; ++s) {
                                Ag(row, static_cast<std::size_t>(phash(dv, 1))) =
                                    num_traits<T>::from_int(Ncur[s] - dv[s]);
                                std::vector<int> sh = dv;
                                sh[s] += 1;  // UP shift
                                const long cb = phash(sh, 1);
                                if (cb >= 0) {
                                    Ag(row, static_cast<std::size_t>(cb)) = -Zs[s];
                                    DCg(row, static_cast<std::size_t>(cb)) = Lr(M - 1, s);
                                }
                                for (std::size_t k = 1; k + 1 <= M; ++k) {
                                    const long c = phash(sh, k + 1);
                                    if (c >= 0) Ag(row, static_cast<std::size_t>(c)) = -Lr(k - 1, s);
                                }
                                ++row;
                            }
                        }
                        // extra population constraint of class r, always present
                        Ag(row, static_cast<std::size_t>(phash(dv, 1))) =
                            num_traits<T>::from_int(Ncur[r] - dv[r]);
                        Bg(row, static_cast<std::size_t>(phash(dv, 1))) = Zs[r];
                        for (std::size_t k = 1; k + 1 <= M; ++k)
                            Bg(row, static_cast<std::size_t>(phash(dv, k + 1))) = Lr(k - 1, r);
                        DDg(row, static_cast<std::size_t>(phash(dv, 1))) = Lr(M - 1, r);
                        ++row;
                    }
                }
                if (row != numRows) throw NumericError("pfqn_procomom: row count mismatch");

                // ---- level recursion -------------------------------------------
                int sumNcur = 0;
                for (int x : Ncur) sumNcur += x;
                const T tol = line::detail::lstsq_tolerance(Ag);
                for (int n = 0; n <= sumNcur; ++n) {
                    std::vector<T> rhs(numRows, zero);
                    for (std::size_t i = 0; i < numRows; ++i) {
                        T s = zero;
                        for (std::size_t j = 0; j < basisSize; ++j)
                            s += Bg(i, j) * pklast(j, static_cast<std::size_t>(n));
                        if (n >= 1) {
                            const T nn = num_traits<T>::from_int(n);
                            T sc = zero, sd = zero;
                            for (std::size_t j = 0; j < basisSize; ++j) {
                                sc += DCg(i, j) * pk(j, static_cast<std::size_t>(n) - 1);
                                sd += DDg(i, j) * pklast(j, static_cast<std::size_t>(n) - 1);
                            }
                            s += nn * sc + nn * sd;
                        }
                        rhs[i] = s;
                    }
                    const LstsqResult<T> sol = lstsq(Ag, rhs, tol);
                    if (sol.rankdef) res.rankdef = true;
                    for (std::size_t j = 0; j < basisSize; ++j)
                        pk(j, static_cast<std::size_t>(n)) = sol.x[j];
                }
            }
        }

        const std::size_t outRow = static_cast<std::size_t>(phash(zeroDn, 1));
        T total = zero;
        for (std::size_t j = 0; j < NP; ++j) total += pk(outRow, j);
        if (total != zero)
            for (std::size_t j = 0; j < NP; ++j) res.Pr(station, j) = pk(outRow, j) / total;
    }

    for (std::size_t i = 0; i < M; ++i) {
        T q = zero;
        for (std::size_t j = 0; j < NP; ++j) q += num_traits<T>::from_int(static_cast<long>(j)) * res.Pr(i, j);
        res.Q[i] = q;
    }
    return res;
}

/** Overload with the reference's default tolerance. */
template <class T>
ProcomomResult<T> pfqn_procomom(const Matrix<T>& L, const std::vector<int>& N,
                                const std::vector<T>& Z) {
    return pfqn_procomom(L, N, Z, num_traits<T>::from_double(1e-14));
}

/** Return value of pfqn_procomom2, mirroring [pk, lG, G, T, F, B]. */
template <class T>
struct Procomom2Result {
    std::vector<T> pk;             ///< (sumN+1) marginal, pk[n] = P(n jobs at the queue)
    T G;                           ///< normalizing constant
    double lG;                     ///< its logarithm
    std::vector<Matrix<T>> Tr;     ///< per-class transfer matrices
    Matrix<T> F;                   ///< prod_r Tr^{N_r} / N_r!
    Matrix<T> B;                   ///< prod_r Tr
};

/**
 * Queue-plus-delay marginal by the transfer-matrix form of ProCoMoM.
 *
 * @param L  (R) demands at the single queueing station
 * @param N  (R) populations
 * @param Z  (R) think times
 * @param mu (sumN) load-dependent rates of the queue, mu[n-1] with n jobs;
 *           empty for the load-independent default
 * @param m  multiplicity of the queueing station
 */
template <class T>
Procomom2Result<T> pfqn_procomom2(const std::vector<T>& L, const std::vector<int>& N,
                                  const std::vector<T>& Z, const std::vector<T>& mu, int m) {
    const std::size_t R = L.size();
    if (R == 0) throw InputError("pfqn_procomom2: empty demand vector");
    if (N.size() != R) throw InputError("pfqn_procomom2: L and N disagree on the class count");
    if (Z.size() != R) throw InputError("pfqn_procomom2: L and Z disagree on the class count");
    if (m < 1) throw InputError("pfqn_procomom2: the multiplicity must be at least one");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    int sumN = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < 0) throw InputError("pfqn_procomom2: negative population");
        sumN += N[r];
    }
    const std::size_t dim = static_cast<std::size_t>(sumN) + 1;

    // mu_[n] is the reference's mu(1+n): mu_[0] = 1, mu_[n] the rate with n jobs.
    std::vector<T> mu_(dim, one);
    if (!mu.empty()) {
        if (mu.size() < static_cast<std::size_t>(sumN))
            throw InputError("pfqn_procomom2: mu must supply a rate for every population up to sum(N)");
        for (std::size_t n = 1; n < dim; ++n) {
            if (mu[n - 1] <= zero)
                throw InputError("pfqn_procomom2: the load-dependent rates must be positive");
            mu_[n] = mu[n - 1];
        }
    }

    Procomom2Result<T> res;
    res.Tr.reserve(R);
    for (std::size_t r = 0; r < R; ++r) {
        Matrix<T> Tm(dim, dim, zero);
        for (int n = sumN; n >= 1; --n) {
            const std::size_t row = static_cast<std::size_t>(sumN - n);
            Tm(row, row) = Z[r];
            Tm(row, row + 1) =
                num_traits<T>::from_int(n + m - 1) * L[r] / mu_[static_cast<std::size_t>(n)];
        }
        Tm(dim - 1, dim - 1) = Z[r];
        res.Tr.push_back(Tm);
    }

    const auto matmul = [&](const Matrix<T>& A, const Matrix<T>& Bm) {
        Matrix<T> C(dim, dim, zero);
        for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t k = 0; k < dim; ++k) {
                if (A(i, k) == zero) continue;
                for (std::size_t j = 0; j < dim; ++j) C(i, j) += A(i, k) * Bm(k, j);
            }
        return C;
    };

    Matrix<T> F(dim, dim, zero), B(dim, dim, zero);
    for (std::size_t i = 0; i < dim; ++i) {
        F(i, i) = one;
        B(i, i) = one;
    }
    for (std::size_t r = 0; r < R; ++r) {
        Matrix<T> P(dim, dim, zero);
        for (std::size_t i = 0; i < dim; ++i) P(i, i) = one;
        for (int e = 0; e < N[r]; ++e) P = matmul(P, res.Tr[r]);
        const T fac = num_factorial<T>(static_cast<unsigned>(N[r]));
        for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t j = 0; j < dim; ++j) P(i, j) /= fac;
        F = matmul(F, P);
        B = matmul(B, res.Tr[r]);
    }
    res.F = F;
    res.B = B;

    // pk = (F e_last)', reversed so that pk[n] = P(n jobs at the queue)
    std::vector<T> v(dim, zero);
    for (std::size_t i = 0; i < dim; ++i) v[i] = F(i, dim - 1);
    T G = zero;
    for (std::size_t i = 0; i < dim; ++i) G += v[i];
    res.G = G;
    res.lG = num_traits<T>::log_as_double(G);
    res.pk.assign(dim, zero);
    if (G != zero)
        for (std::size_t i = 0; i < dim; ++i) res.pk[dim - 1 - i] = v[i] / G;
    return res;
}

/** Overload with the load-independent, unit-multiplicity defaults. */
template <class T>
Procomom2Result<T> pfqn_procomom2(const std::vector<T>& L, const std::vector<int>& N,
                                  const std::vector<T>& Z) {
    return pfqn_procomom2(L, N, Z, std::vector<T>(), 1);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PROCOMOM_H
