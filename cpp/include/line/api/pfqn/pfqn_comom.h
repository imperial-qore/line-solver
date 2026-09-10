/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_COMOM_H
#define LINE_API_PFQN_COMOM_H

/**
 * CoMoM (class-oriented method of moments), the general basis formulation, and
 * the original repairman-model implementation that preceded pfqn_comomrm.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_comom.m and
 * matlab/src/api/pfqn/pfqn_comomrm_orig.m. Both accept at most ONE queueing
 * station (plus a delay); the reference raises an error for M > 1 and so does
 * this port.
 *
 * WHERE THIS DIFFERS FROM pfqn_comomrm. pfqn_comomrm (already ported) writes
 * C^{-1} out in closed form and never forms a linear system. pfqn_comom builds
 * the full (numDn (M+1)) x (numDn (M+1)) matrices A, B and DA of the CoMoM
 * basis from scratch on the first job of each class, updates A by DA on every
 * later job, and SOLVES A h = B h n_r / (sum(n) + M - 1) at each step. It is
 * the slower but structurally general form, and it is a genuinely independent
 * route to the same normalizing constant, which is what makes it worth
 * carrying as a cross-check.
 *
 * BASIS LAYOUT. Dn = multichoose(R, M) with its LAST column zeroed and the
 * rows reordered by sort_by_nnz_pos. For M = 1 that is the R-row set
 * {0, e_1, ..., e_{R-1}} with the zero row first. The basis position of the
 * moment G(n - d) with "level" index i (i = 1 is the plain constant, i > 1 the
 * per-station moment) is
 *
 *   col(d, 1)   = numDn M + pos(d),      col(d, i>1) = pos(d) M + i - 2
 *
 * in zero-based columns, which is the reference's hash() with its 1-based
 * indexing removed. Call sites index into that layout by POSITION, so the
 * enumeration order of multichoose and the bubble sort that follows it are
 * reproduced verbatim in pfqn_comb_common.h rather than replaced.
 *
 * SCALING. The reference renormalizes h by |sum(h)| after every job and folds
 * the discarded factors into a log accumulator, then also takes abs(h). The
 * basis entries are normalizing constants and are positive, so the abs() is a
 * no-op and the scale factors cancel in the final lG. This port carries the
 * UNSCALED basis, exactly as pfqn_comomrm.h does, which is what makes
 *
 *   G(N) = h[matrixDim - R] * (sum(N) + M - 1)! / prod_r N_r! * prod_r Lmax_r^{N_r}
 *
 * an exact identity in the field of the inputs.
 *
 * SORTING IN pfqn_comomrm_orig. The reference sorts the classes by ascending
 * think time AFTER calling pfqn_nc_sanitize, but permutes only L and Z, leaving
 * N in place. That would desynchronize N from its demands, except that
 * pfqn_nc_sanitize has ALREADY ordered the classes by ascending think time, so
 * the second sort is the identity on every input. This port permutes L, Z and N
 * together, which agrees with the reference wherever the reference is defined
 * and is correct if the ordering guarantee were ever to change. Verified
 * identical on L = [0.5 0.25 0.125], N = [1 2 1], Z = [4 0.25 0.5], whose
 * post-sanitize think times [8 1 4] are NOT sorted before the second pass:
 * MATLAB returns lG = -0.439231970578982 against the exact -0.43923197057898189.
 *
 * Arithmetic: EXACT-CAPABLE. Additions, multiplications, divisions and one
 * linear solve per job, all in the field of the inputs. The reference's logs
 * are only the scale bookkeeping and the factln seeding, both of which are
 * ratios of factorials formed directly here.
 *
 * REFERENCE DEFECTS: none found in either routine. Both reproduce pfqn_ca to
 * the last few ulps on every model tried, including the non-identity think-time
 * ordering above.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_comb_common.h"
#include "line/api/pfqn/pfqn_comomrm.h"
#include "line/api/pfqn/pfqn_nc_sanitize.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** Dn basis of pfqn_comom / pfqn_procomom: multichoose(R,M), last column zeroed, sorted. */
inline std::vector<std::vector<int>> comom_basis(int R, int M) {
    std::vector<std::vector<int>> Dn = multichoose_rows(R, M);
    for (std::size_t d = 0; d < Dn.size(); ++d) Dn[d][static_cast<std::size_t>(R) - 1] = 0;
    sort_by_nnz_pos(Dn);
    return Dn;
}

}  // namespace detail

/**
 * CoMoM on the general basis (matlab pfqn_comom.m).
 *
 * @param L    (1 x R) demands at the single queueing station
 * @param N    (R) populations
 * @param Z    (R) think times
 * @param atol tolerance below which a demand counts as zero
 */
template <class T>
ComomResult<T> pfqn_comom(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                          const T& atol) {
    const std::size_t M = L.rows();
    const std::size_t R = L.cols();
    if (M > 1)
        throw InputError(
            "pfqn_comom: at most one queueing station is supported (repairman models with a "
            "delay); use pfqn_ca or pfqn_recal for M > 1");
    if (M == 0 || R == 0) throw InputError("pfqn_comom: empty demand matrix");
    if (N.size() != R) throw InputError("pfqn_comom: L and N disagree on the class count");
    if (Z.size() != R) throw InputError("pfqn_comom: L and Z disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    // ---- per-class rescaling (reference: Lmax with the sub-atol entries taken from Z) --
    std::vector<T> Lmax(R, one);
    for (std::size_t r = 0; r < R; ++r) {
        T mx = (L(0, r) < atol) ? Z[r] : L(0, r);
        for (std::size_t i = 1; i < M; ++i) {
            const T v = (L(i, r) < atol) ? Z[r] : L(i, r);
            if (v > mx) mx = v;
        }
        if (mx <= zero) {
            if (N[r] > 0)
                throw InputError(
                    "pfqn_comom: a populated class has neither service demand nor think time");
            mx = one;
        }
        Lmax[r] = mx;
    }
    Matrix<T> Ls(M, R, zero);
    std::vector<T> Zs(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        for (std::size_t i = 0; i < M; ++i) Ls(i, r) = L(i, r) / Lmax[r];
        Zs[r] = Z[r] / Lmax[r];
    }

    // ---- basis layout ---------------------------------------------------------------
    const std::vector<std::vector<int>> Dn =
        detail::comom_basis(static_cast<int>(R), static_cast<int>(M));
    const std::size_t numDn = Dn.size();
    const std::size_t dim = numDn * (M + 1);

    // zero-based column of moment (d, level i), i = 1..M+1 in the reference's numbering
    const auto col = [&](const std::vector<int>& d, std::size_t i) {
        const int pos = matchrow(Dn, d);
        if (pos < 0) throw NumericError("pfqn_comom: basis vector outside the Dn set");
        const std::size_t p = static_cast<std::size_t>(pos);
        return (i == 1) ? (numDn * M + p) : (p * M + i - 2);
    };

    std::vector<int> zeroRow(R, 0);
    std::vector<T> h(dim, zero);
    for (std::size_t i = 0; i <= M; ++i) h[col(zeroRow, i + 1)] = one;

    // ---- iterate one job at a time ---------------------------------------------------
    std::vector<int> nvec(R, 0);
    Matrix<T> A(dim, dim, zero), B(dim, dim, zero), DA(dim, dim, zero);
    for (std::size_t r = 0; r < R; ++r) {
        for (int Nr = 1; Nr <= N[r]; ++Nr) {
            nvec[r] += 1;
            if (Nr == 1) {
                A = Matrix<T>(dim, dim, zero);
                B = Matrix<T>(dim, dim, zero);
                DA = Matrix<T>(dim, dim, zero);
                std::size_t row = 0;
                for (std::size_t d = 0; d < numDn; ++d) {
                    const std::vector<int>& dv = Dn[d];
                    int s1 = 0;  // sum of dv over the reference's columns r..R-1 (1-based)
                    for (std::size_t c = r; c + 1 < R; ++c) s1 += dv[c];
                    if (s1 > 0) {
                        int s2 = 0;
                        for (std::size_t c = r + 1; c + 1 < R; ++c) s2 += dv[c];
                        for (std::size_t k = 0; k <= M; ++k) {
                            const std::size_t c1 = col(dv, k + 1);
                            A(row, c1) = one;
                            if (s2 > 0) {
                                B(row, c1) = one;
                            } else {
                                std::vector<int> dm = dv;
                                dm[r] -= 1;
                                B(row, col(dm, k + 1)) = one;
                            }
                            ++row;
                        }
                    } else {
                        int s3 = 0;
                        for (std::size_t c = 0; c <= r; ++c) s3 += dv[c];
                        if (s3 < static_cast<int>(M)) {
                            for (std::size_t k = 1; k <= M; ++k) {
                                A(row, col(dv, k + 1)) = one;
                                A(row, col(dv, 1)) = -one;
                                for (std::size_t s = 0; s < r; ++s) {
                                    std::vector<int> dp = dv;
                                    dp[s] += 1;
                                    A(row, col(dp, k + 1)) = -Ls(k - 1, s);
                                }
                                B(row, col(dv, k + 1)) = Ls(k - 1, r);
                                ++row;
                            }
                            for (std::size_t s = 0; s < r; ++s) {
                                A(row, col(dv, 1)) = num_traits<T>::from_int(nvec[s] - dv[s]);
                                std::vector<int> dp = dv;
                                dp[s] += 1;
                                A(row, col(dp, 1)) = -Zs[s];
                                for (std::size_t k = 1; k <= M; ++k)
                                    A(row, col(dp, k + 1)) = -Ls(k - 1, s);
                                ++row;
                            }
                        }
                    }
                }
                // population constraint of the class being filled
                for (std::size_t d = 0; d < numDn; ++d) {
                    const std::vector<int>& dv = Dn[d];
                    int s1 = 0;
                    for (std::size_t c = r; c + 1 < R; ++c) s1 += dv[c];
                    if (s1 > 0) continue;
                    const std::size_t c0 = col(dv, 1);
                    A(row, c0) = num_traits<T>::from_int(nvec[r] - dv[r]);
                    DA(row, c0) = one;
                    B(row, c0) = Zs[r];
                    for (std::size_t k = 1; k <= M; ++k) B(row, col(dv, k + 1)) = Ls(k - 1, r);
                    ++row;
                }
                if (row != dim)
                    throw NumericError("pfqn_comom: the CoMoM system is not square");
            } else {
                for (std::size_t i = 0; i < dim; ++i)
                    for (std::size_t j = 0; j < dim; ++j) A(i, j) += DA(i, j);
            }
            int nt = 0;
            for (int x : nvec) nt += x;
            const T fac = num_traits<T>::from_int(nvec[r]) /
                          num_traits<T>::from_int(nt + static_cast<long>(M) - 1);
            std::vector<T> b(dim, zero);
            for (std::size_t i = 0; i < dim; ++i) {
                T s = zero;
                for (std::size_t j = 0; j < dim; ++j) s += B(i, j) * h[j];
                b[i] = s * fac;
            }
            h = solve(A, b);
        }
    }

    // ---- undo the scaling ------------------------------------------------------------
    int Ntot = 0;
    for (int x : N) Ntot += x;
    T fact = num_factorial<T>(static_cast<unsigned>(Ntot + static_cast<long>(M) - 1));
    for (std::size_t r = 0; r < R; ++r) {
        fact /= num_factorial<T>(static_cast<unsigned>(N[r]));
        fact *= num_pow_int(Lmax[r], static_cast<unsigned>(N[r]));
    }

    ComomResult<T> res;
    res.basis = h;
    res.G = fact * h[dim - R];
    res.lG = num_traits<T>::log_as_double(res.G);
    return res;
}

/** Overload with the reference's default tolerance. */
template <class T>
ComomResult<T> pfqn_comom(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z) {
    return pfqn_comom(L, N, Z, num_traits<T>::from_double(1e-14));
}

/**
 * Original CoMoM for the finite repairman model (matlab pfqn_comomrm_orig.m).
 *
 * @param L    (1 x R) demands at the single queueing station
 * @param N    (R) populations
 * @param Z    (K x R) think times
 * @param atol tolerance passed to pfqn_nc_sanitize
 */
template <class T>
ComomResult<T> pfqn_comomrm_orig(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                                 const T& atol) {
    if (!L.empty() && L.rows() != 1)
        throw InputError("pfqn_comomrm_orig: the solver accepts at most a single queueing station");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    // basis multiplicity rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    const T mT = one;

    const NcSanitizeResult<T> san = pfqn_nc_sanitize(std::vector<T>(), L, N, Z, atol);
    const std::size_t R = san.N.size();

    ComomResult<T> res;
    if (R == 0) {
        res.G = san.Gremaind;
        res.lG = san.lGremaind;
        res.basis.assign(1, one);
        return res;
    }

    // ---- second rescaling, an identity after pfqn_nc_sanitize but reproduced ---------
    std::vector<T> Lv(R, zero), Zv(R, zero), Lmax(R, one);
    for (std::size_t r = 0; r < R; ++r) {
        Lv[r] = san.L.empty() ? zero : san.L(0, r);
        for (std::size_t k = 0; k < san.Z.rows(); ++k) Zv[r] += san.Z(k, r);
        T mx = (Lv[r] < atol) ? Zv[r] : Lv[r];
        if (mx <= zero) {
            if (san.N[r] > 0)
                throw InputError(
                    "pfqn_comomrm_orig: a populated class has neither demand nor think time");
            mx = one;
        }
        Lmax[r] = mx;
    }
    std::vector<int> Nv(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        Lv[r] /= Lmax[r];
        Zv[r] /= Lmax[r];
        Nv[r] = san.N[r];
    }

    // ---- ascending think time; identity given pfqn_nc_sanitize's ordering -----------
    std::vector<std::size_t> ord(R);
    for (std::size_t r = 0; r < R; ++r) ord[r] = r;
    std::stable_sort(ord.begin(), ord.end(),
                     [&](std::size_t a, std::size_t b) { return Zv[a] < Zv[b]; });
    std::vector<T> Lp(R, zero), Zp(R, zero);
    std::vector<int> Np(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        Lp[r] = Lv[ord[r]];
        Zp[r] = Zv[ord[r]];
        Np[r] = Nv[ord[r]];
    }

    // ---- the transfer-matrix recursion on the unscaled basis -------------------------
    std::vector<int> nvec(R, 0);
    std::vector<T> h(2, one), hprev(2, one);
    Matrix<T> F1, F2;
    for (std::size_t r = 0; r < R; ++r) {
        const std::size_t rr = r + 1;  // the reference's 1-based class index
        for (int Nr = 1; Nr <= Np[r]; ++Nr) {
            nvec[r] += 1;
            int nt = 0;
            for (int x : nvec) nt += x;
            if (Nr == 1) {
                if (rr > 1) {
                    // basis expansion rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                    const std::size_t p = rr - 1;
                    const T w = num_traits<T>::from_int(Np[r - 1]) /
                                num_traits<T>::from_int(nt - 1);
                    std::vector<T> hn(2 * rr, zero);
                    for (std::size_t i = 0; i < p; ++i) hn[i] = h[i];
                    for (std::size_t i = 0; i < p; ++i) hn[rr + i] = h[p + i];
                    hn[p] = hprev[0] * w;
                    hn[2 * rr - 1] = hprev[p] * w;
                    h.swap(hn);
                }
                // C, A12, B2r of the reference, then F1r and F2r.
                Matrix<T> C(rr, rr, zero), A12(rr, rr, zero), B1r(rr, 2 * rr, zero),
                    B2r(rr, 2 * rr, zero);
                C(0, 0) = one;
                for (std::size_t s = 0; s + 1 < rr; ++s) C(0, 1 + s) = -Lp[s];
                A12(0, 0) = -one;
                B1r(0, 0) = Lp[r];
                for (std::size_t s = 0; s + 1 < rr; ++s) {
                    A12(1 + s, 0) = num_traits<T>::from_int(Np[s]);
                    A12(1 + s, 1 + s) = -Zp[s];
                    C(1 + s, 1 + s) = -mT * Lp[s];
                }
                for (std::size_t i = 0; i < rr; ++i) {
                    B2r(i, i) = mT * Lp[r];
                    B2r(i, rr + i) = Zp[r];
                }
                // F1r = [ C^{-1} B1r ; 0 ],  F2r = [ -C^{-1} A12 B2r ; B2r ]
                Matrix<T> LU = C;
                const std::vector<std::size_t> piv = lu_factor(LU);
                Matrix<T> iCB1(rr, 2 * rr, zero), iCA12(rr, rr, zero);
                for (std::size_t j = 0; j < 2 * rr; ++j) {
                    std::vector<T> rhs(rr, zero);
                    for (std::size_t i = 0; i < rr; ++i) rhs[i] = B1r(i, j);
                    lu_solve(LU, piv, rhs);
                    for (std::size_t i = 0; i < rr; ++i) iCB1(i, j) = rhs[i];
                }
                for (std::size_t j = 0; j < rr; ++j) {
                    std::vector<T> rhs(rr, zero);
                    for (std::size_t i = 0; i < rr; ++i) rhs[i] = A12(i, j);
                    lu_solve(LU, piv, rhs);
                    for (std::size_t i = 0; i < rr; ++i) iCA12(i, j) = rhs[i];
                }
                F1 = Matrix<T>(2 * rr, 2 * rr, zero);
                F2 = Matrix<T>(2 * rr, 2 * rr, zero);
                for (std::size_t i = 0; i < rr; ++i)
                    for (std::size_t j = 0; j < 2 * rr; ++j) {
                        F1(i, j) = iCB1(i, j);
                        T s = zero;
                        for (std::size_t k = 0; k < rr; ++k) s += iCA12(i, k) * B2r(k, j);
                        F2(i, j) = -s;
                        F2(rr + i, j) = B2r(i, j);
                    }
            }
            hprev = h;
            const T nr = num_traits<T>::from_int(nvec[r]);
            const T den = num_traits<T>::from_int(nt);  // sum(nvec) + M - 1 with M = 1
            std::vector<T> hn(2 * rr, zero);
            for (std::size_t i = 0; i < 2 * rr; ++i) {
                T s = zero;
                for (std::size_t j = 0; j < 2 * rr; ++j) s += (nr * F1(i, j) + F2(i, j)) * hprev[j];
                hn[i] = s / den;
            }
            h.swap(hn);
        }
    }

    int Ntot = 0;
    for (std::size_t r = 0; r < R; ++r) Ntot += Np[r];
    T fact = num_factorial<T>(static_cast<unsigned>(Ntot));  // (sum N + M - 1)! with M = 1
    for (std::size_t r = 0; r < R; ++r) {
        fact /= num_factorial<T>(static_cast<unsigned>(Np[r]));
        fact *= num_pow_int(Lmax[r], static_cast<unsigned>(Nv[r]));
    }

    res.basis = h;
    res.G = san.Gremaind * fact * h[h.size() - R];
    res.lG = num_traits<T>::log_as_double(res.G);
    return res;
}

/** Overload with the exact (zero-tolerance) tests. */
template <class T>
ComomResult<T> pfqn_comomrm_orig(const Matrix<T>& L, const std::vector<int>& N,
                                 const Matrix<T>& Z) {
    return pfqn_comomrm_orig(L, N, Z, num_traits<T>::from_int(0));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_COMOM_H
