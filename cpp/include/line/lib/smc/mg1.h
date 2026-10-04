/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LIB_SMC_MG1_H
#define LINE_LIB_SMC_MG1_H

/**
 * @file
 * @ingroup api_mam
 * The M/G/1-type and GI/M/1-type fundamental-matrix solvers of MAMSolver /
 * SMCSolver, ported from `matlab/lib/thirdparty/MG1files`:
 * `stat.m`, `MG1_EG.m`, `MG1_Decay.m`, `GIM1_Caudal.m`, `MG1_Shifts.m`,
 * `MG1_CR.m`, `MG1_FI.m`, `MG1_NI.m` (with `solveSylvPowersDirectSum.m` and
 * `solveSylvPowersRealSchur_FW.m`), `MG1_RR.m` (with `MG1_RR_Btemp.m` and `MG1_RR_tempB.m`),
 * `MG1_IS.m` and `GIM1_R.m`.
 *
 * These are the third-party numerics the ETAQA aggregation sits on: `MG1_CR`
 * returns the minimal nonnegative G of an M/G/1-type chain and `GIM1_R` the
 * minimal nonnegative R of a GI/M/1-type one, and without them
 * `solver_mam_bmap_map_1` and `solver_mam_map_bmap_1` have nothing to
 * aggregate. The port follows the MATLAB line by line, including the block
 * index arithmetic, the stopping tests and the constants, so that a divergence
 * against the reference is a bug here and not a design difference.
 *
 * BLOCK LAYOUT. The reference passes a block sequence as one wide matrix
 * `A = [A0 A1 A2 ... Amax]`, m rows by m*(max+1) columns. This port carries the
 * same sequence as a `std::vector<Matrix<double>>` of m x m blocks, which is
 * the identical object with the index arithmetic done once, in `blocks_of` /
 * `hcat`, instead of at every use. The GI/M/1 side stacks its blocks
 * VERTICALLY in the reference; `GIM1_R_ETAQA` is what transposes that stack
 * into the horizontal one, so everything below is horizontal.
 *
 * DOUBLE ONLY, and not by preference. `MG1_Decay` and `GIM1_Caudal` bisect on
 * the Perron-Frobenius eigenvalue of A(z), which needs LAPACK; `MG1_CR`
 * evaluates its polynomials at complex roots of unity through an FFT, whose
 * twiddle factors are cos/sin; `MG1_pi_ETAQA` drops a column chosen by a
 * numerical rank test, which needs an SVD. None of the three has a
 * multiprecision or exact counterpart in this tree, so the whole family is
 * declared on `Matrix<double>` and the solvers that call it refuse at any
 * other arithmetic rather than down-converting behind the caller's back.
 *
 * `MG1_Shifts` IS PORTED IN FULL, all three ShiftTypes at either drift. Its
 * reference writes the last block of a row shift as
 *     rowhatA(1,maxd*i:end) = uT*A(:,maxd*i:end)
 * with `i` the value 1 left over from the `beta` loop, so the line addresses
 * column `maxd`, not block `maxd`. It is nevertheless CORRECT: both sides take
 * the same columns, so columns `maxd*m+1:end` receive exactly `uT*A_maxd`, and
 * the stray columns before them are overwritten by the loop that follows. An
 * earlier version of this header read that line as defective and refused the
 * branches; the transient-chain Newton iteration under GIM1_R's Ramaswami dual
 * is what reaches them. With maxd = 1 the `beta` loop never runs, `i` is the
 * imaginary unit and the reference errors; this port computes the last block
 * the line intends.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/num/complex_number.h"
#include "line/num/number.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/fft.h"
#include "line/util/linalg.h"
#include "line/util/lstsq.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"
#include "line/util/svd.h"

namespace line {
namespace smc {

using Blocks = std::vector<Matrix<double>>;

// ---------------------------------------------------------------------------
// Small matrix helpers, named after the MATLAB they stand for
// ---------------------------------------------------------------------------

/** A + B. Templated so the point-wise CR step can add complex blocks. */
template <class T>
Matrix<T> madd(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols())
        throw InputError("smc: matrix addition with mismatched shapes");
    Matrix<T> C(A.rows(), A.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) + B(i, j);
    return C;
}

/** A - B. */
template <class T>
Matrix<T> msub(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols())
        throw InputError("smc: matrix subtraction with mismatched shapes");
    Matrix<T> C(A.rows(), A.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) - B(i, j);
    return C;
}

/** c * A. */
inline Matrix<double> mscale(const Matrix<double>& A, double c) {
    Matrix<double> C(A.rows(), A.cols(), 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = c * A(i, j);
    return C;
}

/** sum(A,2), the row sums, as a column held in a vector. */
inline std::vector<double> rowsums(const Matrix<double>& A) {
    std::vector<double> s(A.rows(), 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) s[i] += A(i, j);
    return s;
}

/** norm(A,inf), the largest absolute row sum. */
inline double inf_norm(const Matrix<double>& A) {
    double best = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < A.cols(); ++j) s += std::fabs(A(i, j));
        if (s > best) best = s;
    }
    return best;
}

/** norm(A,inf) over a whole block sequence stacked vertically. */
inline double inf_norm(const Blocks& blk, std::size_t from) {
    double best = 0.0;
    for (std::size_t i = from; i < blk.size(); ++i) best = std::max(best, inf_norm(blk[i]));
    return best;
}

/** max(max(abs(A-B))). */
inline double max_abs_diff(const Matrix<double>& A, const Matrix<double>& B) {
    double best = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            best = std::max(best, std::fabs(A(i, j) - B(i, j)));
    return best;
}

/** max(sum(A)), the largest column sum WITHOUT absolute values, as in MATLAB. */
inline double max_col_sum(const Matrix<double>& A) {
    double best = -std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < A.cols(); ++j) {
        double s = 0.0;
        for (std::size_t i = 0; i < A.rows(); ++i) s += A(i, j);
        if (s > best) best = s;
    }
    return best;
}

/** Splits the wide `[A0 A1 ... Amax]` into its m x m blocks. */
inline Blocks blocks_of(const Matrix<double>& A, std::size_t m) {
    if (m == 0 || A.cols() % m != 0)
        throw InputError("smc: the block sequence has an incorrect number of columns");
    const std::size_t nb = A.cols() / m;
    Blocks out(nb, Matrix<double>(A.rows(), m, 0.0));
    for (std::size_t b = 0; b < nb; ++b)
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < m; ++j) out[b](i, j) = A(i, b * m + j);
    return out;
}

/** Re-assembles a block sequence into the wide `[A0 A1 ... Amax]`. */
inline Matrix<double> hcat(const Blocks& blk) {
    if (blk.empty()) return Matrix<double>();
    const std::size_t r = blk[0].rows(), c = blk[0].cols();
    Matrix<double> A(r, c * blk.size(), 0.0);
    for (std::size_t b = 0; b < blk.size(); ++b)
        for (std::size_t i = 0; i < r; ++i)
            for (std::size_t j = 0; j < c; ++j) A(i, b * c + j) = blk[b](i, j);
    return A;
}

/** Stacks a block sequence vertically, `[A0; A1; ...; Amax]`. */
inline Matrix<double> vcat(const Blocks& blk) {
    if (blk.empty()) return Matrix<double>();
    const std::size_t r = blk[0].rows(), c = blk[0].cols();
    Matrix<double> A(r * blk.size(), c, 0.0);
    for (std::size_t b = 0; b < blk.size(); ++b)
        for (std::size_t i = 0; i < r; ++i)
            for (std::size_t j = 0; j < c; ++j) A(b * r + i, j) = blk[b](i, j);
    return A;
}

/** Splits a vertical stack into its blocks of `r` rows. */
inline Blocks vblocks_of(const Matrix<double>& A, std::size_t r) {
    if (r == 0 || A.rows() % r != 0)
        throw InputError("smc: the stacked block sequence has an incorrect number of rows");
    const std::size_t nb = A.rows() / r;
    Blocks out(nb, Matrix<double>(r, A.cols(), 0.0));
    for (std::size_t b = 0; b < nb; ++b)
        for (std::size_t i = 0; i < r; ++i)
            for (std::size_t j = 0; j < A.cols(); ++j) out[b](i, j) = A(b * r + i, j);
    return out;
}

// ---------------------------------------------------------------------------
// stat.m
// ---------------------------------------------------------------------------

/**
 * Stationary distribution of a stochastic matrix: the left eigenvector for
 * eigenvalue 1, nonnegative and summing to one.
 *
 * Port of `stat.m`, including its shape: `[A - I, e]` is S x (S+1), so the
 * reference's `y / B` is a least-squares solve of a consistent overdetermined
 * system and not a square solve. The normalization is IN the system (the
 * appended column of ones against the appended 1 on the right), which is why
 * the result needs no rescaling afterwards.
 */
inline std::vector<double> stat(const Matrix<double>& A) {
    const std::size_t S = A.rows();
    if (A.cols() != S) throw InputError("stat: matrix is not square");
    // x [A - I, e] = [0 ... 0 1]  <=>  [A - I, e]^T x^T = [0 ... 0 1]^T
    Matrix<double> Bt(S + 1, S, 0.0);
    for (std::size_t i = 0; i < S; ++i)
        for (std::size_t j = 0; j < S; ++j) Bt(j, i) = A(i, j) - (i == j ? 1.0 : 0.0);
    for (std::size_t i = 0; i < S; ++i) Bt(S, i) = 1.0;
    std::vector<double> y(S + 1, 0.0);
    y[S] = 1.0;
    return lstsq(Bt, y, detail::lstsq_tolerance(Bt)).x;
}

/** theta A, the row vector times matrix product used throughout. */
inline std::vector<double> rowvec_times(const std::vector<double>& v, const Matrix<double>& A) {
    return vecmul(v, A);
}

/** The inner product of a row vector with a column held as a vector. */
inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw InputError("smc: inner product with mismatched lengths");
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

/** The drift of an M/G/1-type sequence, and the invariant vector it uses. */
struct Drift {
    double value = 0.0;
    std::vector<double> theta;  ///< stat(A0 + A1 + ... + Amax)
};

/**
 * `drift = theta * beta` with `beta = (Amax)e + (Amax+Amax-1)e + ...`, the
 * expected level increment per transition of the phase process. Repeated
 * verbatim in MG1_EG, MG1_Shifts, MG1_pi_ETAQA and GIM1_R, so it lives here.
 */
inline Drift mg1_drift(const Blocks& A) {
    const std::size_t dega = A.size() - 1;
    Matrix<double> sumA = A[dega];
    std::vector<double> beta = rowsums(sumA);
    for (std::size_t i = dega; i-- > 1;) {
        sumA = madd(sumA, A[i]);
        const std::vector<double> rs = rowsums(sumA);
        for (std::size_t k = 0; k < beta.size(); ++k) beta[k] += rs[k];
    }
    sumA = madd(sumA, A[0]);
    Drift d;
    d.theta = stat(sumA);
    d.value = dot(d.theta, beta);
    return d;
}

// ---------------------------------------------------------------------------
// MG1_Decay.m and GIM1_Caudal.m
// ---------------------------------------------------------------------------

/**
 * `max(eig(M))` with MATLAB's semantics on a complex spectrum: the element of
 * largest modulus, ties broken by the larger phase angle. For the nonnegative
 * A(z) both callers evaluate, this is the Perron-Frobenius eigenvalue and is
 * real; the comparisons below then take its real part, which is what MATLAB's
 * relational operators do on a complex value.
 */
inline std::complex<double> max_eig(const Matrix<double>& M) {
    const std::vector<std::complex<double>> ev = eig_values(M);
    if (ev.empty()) throw NumericError("smc: empty spectrum");
    std::complex<double> best = ev[0];
    for (const std::complex<double>& z : ev) {
        const double mz = std::abs(z), mb = std::abs(best);
        if (mz > mb || (mz == mb && std::arg(z) > std::arg(best))) best = z;
    }
    return best;
}

/** A(z) = A0 + A1 z + ... + Amax z^max, by Horner as the reference writes it. */
inline Matrix<double> poly_at(const Blocks& A, double z) {
    Matrix<double> temp = A.back();
    for (std::size_t i = A.size() - 1; i-- > 0;) temp = madd(mscale(temp, z), A[i]);
    return temp;
}

/**
 * The Perron-Frobenius eigenvector of M, right (`left == false`) or left, scaled to unit
 * sum: the null vector of M - PF(M) I (or its transpose), read off the SVD. The reference
 * picks the column of `eig` at the largest eigenvalue; both callers normalize it by its sum,
 * so the scale and sign `eig` happens to return do not matter.
 */
inline std::vector<double> pf_vector(const Matrix<double>& M, bool left) {
    const std::size_t m = M.rows();
    const double lambda = max_eig(M).real();
    Matrix<double> K(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) K(i, j) = left ? M(j, i) : M(i, j);
    for (std::size_t i = 0; i < m; ++i) K(i, i) -= lambda;
    const SvdFactors sv = svd_full(K);
    std::vector<double> v(m);
    double s = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        v[i] = sv.Vt(m - 1, i);
        s += v[i];
    }
    for (double& x : v) x /= s;
    return v;
}

/**
 * Decay rate of a recurrent M/G/1-type chain: the unique z > 1 with
 * PF(A(z)) = z. Port of `MG1_Decay.m`. When `uT` is given it receives the
 * left PF eigenvector of A(z) at the last bisection point, as the reference's
 * second output, scaled to unit sum.
 */
inline double mg1_decay(const Blocks& A, std::vector<double>* uT = nullptr) {
    double eta = 1.0, new_eta = 0.0;
    Matrix<double> temp;
    while (new_eta - eta < 0.0) {
        eta += 1.0;
        temp = poly_at(A, eta);
        new_eta = max_eig(temp).real();
    }
    double eta_min = eta - 1.0, eta_max = eta;
    eta = eta_min + 0.5;
    while (eta_max - eta_min > 1e-15) {
        temp = poly_at(A, eta);
        new_eta = max_eig(temp).real();
        if (new_eta < eta) {
            eta_min = eta;
        } else {
            eta_max = eta;
        }
        eta = (eta_min + eta_max) / 2.0;
    }
    if (uT) *uT = pf_vector(temp, true);
    return eta;
}

/**
 * Caudal characteristic of a GI/M/1-type chain: the spectral radius of R, the
 * unique z in (0,1) with PF(A(z)) = z. Port of `GIM1_Caudal.m`. When `v` is
 * given it receives the right PF eigenvector of A(z) at the last bisection
 * point, scaled to unit sum.
 */
inline double gim1_caudal(const Blocks& A, std::vector<double>* v = nullptr) {
    double eta_min = 0.0, eta_max = 1.0, eta = 0.5;
    Matrix<double> temp;
    while (eta_max - eta_min > 1e-15) {
        temp = poly_at(A, eta);
        const double new_eta = max_eig(temp).real();
        if (new_eta > eta) {
            eta_min = eta;
        } else {
            eta_max = eta;
        }
        eta = (eta_min + eta_max) / 2.0;
    }
    if (v) *v = pf_vector(temp, false);
    return eta;
}

// ---------------------------------------------------------------------------
// MG1_Shifts.m
// ---------------------------------------------------------------------------

/** What `MG1_Shifts` returns: the shifted sequence and the drift it measured. */
struct ShiftResult {
    Blocks hatA;
    double drift = 0.0;
    double tau = 1.0;
    std::vector<double> v;
};

/**
 * Shift technique for the M/G/1-type sequence. Port of `MG1_Shifts.m`, ShiftType
 * 'one' (the default), 'tau' and 'dbl'.
 *
 * For a positive recurrent chain (drift < 1) 'one' shifts the eigenvalue 1 of A(z) to zero
 * by subtracting `(A0+...+Ai)e u^T` from block i with `u^T = e^T/m`, and 'tau' shifts the
 * decay rate tau to infinity by subtracting `e rowhatA_i`; for a transient one (drift >= 1)
 * 'one' shifts 1 to infinity through the row `theta(Amax+...+Ai)` and 'tau' shifts the
 * caudal value to zero through the column built on its right eigenvector v. 'dbl' applies
 * both, in the reference's order. The solver converges in the next root instead of stalling
 * on the removed one and puts the removed rank-one term back on G (`mg1_unshift`).
 */
inline ShiftResult mg1_shifts(const Blocks& Ain, const std::string& shift_type) {
    if (shift_type != "one" && shift_type != "tau" && shift_type != "dbl")
        throw InputError("MG1_Shifts: ShiftType '" + shift_type + "' is not one of one, tau, dbl");
    Blocks A = Ain;
    const std::size_t m = A[0].rows();
    const std::size_t maxd = A.size() - 1;
    const Drift d = mg1_drift(A);
    ShiftResult out;
    out.drift = d.value;
    out.tau = 1.0;
    out.v.assign(m, 0.0);
    Blocks hatA(A.size(), Matrix<double>(m, m, 0.0));
    auto minus_I1 = [&](Blocks& X, double s) {
        for (std::size_t i = 0; i < m; ++i) X[1](i, i) += s;
    };
    // hatA_i = A_i - e rowhat_i, rowhat_maxd = u A_maxd, rowhat_i = f*rowhat_{i+1} + u A_i
    auto row_shift = [&](const Blocks& X, const std::vector<double>& u, double f) {
        Blocks R(X.size());
        std::vector<std::vector<double>> row(X.size());
        row[maxd] = vecmul(u, X[maxd]);
        for (std::size_t i = maxd; i-- > 0;) {
            row[i] = vecmul(u, X[i]);
            for (std::size_t j = 0; j < m; ++j) row[i][j] += f * row[i + 1][j];
        }
        for (std::size_t b = 0; b < X.size(); ++b) {
            R[b] = X[b];
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j) R[b](i, j) -= row[b][j];
        }
        return R;
    };

    if (d.value < 1.0) {
        if (shift_type == "tau" || shift_type == "dbl") {  // shift tau to infinity
            std::vector<double> uT;
            out.tau = mg1_decay(A, &uT);
            minus_I1(A, -1.0);
            hatA = row_shift(A, uT, out.tau);
        }
        if (shift_type == "dbl") A = hatA;
        if (shift_type == "one") minus_I1(A, -1.0);
        if (shift_type == "one" || shift_type == "dbl") {  // shift one to zero
            std::vector<double> col(m, 0.0);
            for (std::size_t b = 0; b < A.size(); ++b) {
                const std::vector<double> rs = rowsums(A[b]);
                for (std::size_t i = 0; i < m; ++i) col[i] += rs[i];
                hatA[b] = A[b];
                for (std::size_t i = 0; i < m; ++i)
                    for (std::size_t j = 0; j < m; ++j)
                        hatA[b](i, j) -= col[i] / static_cast<double>(m);
            }
        }
    } else {
        if (shift_type == "one" || shift_type == "dbl") {  // shift one to infinity
            minus_I1(A, -1.0);
            hatA = row_shift(A, d.theta, 1.0);
        }
        if (shift_type == "dbl") {
            A = hatA;
            minus_I1(A, 1.0);
        }
        if (shift_type == "tau" || shift_type == "dbl") {  // shift tau to zero
            std::vector<double> v;
            out.tau = gim1_caudal(A, &v);
            minus_I1(A, -1.0);
            out.v = v;
            std::vector<double> col = mulvec(A[0], v);
            for (std::size_t b = 0; b < A.size(); ++b) {
                if (b > 0) {
                    const std::vector<double> Av = mulvec(A[b], v);
                    for (std::size_t i = 0; i < m; ++i) col[i] = col[i] / out.tau + Av[i];
                }
                hatA[b] = A[b];
                for (std::size_t i = 0; i < m; ++i)
                    for (std::size_t j = 0; j < m; ++j) hatA[b](i, j) -= col[i];
            }
        }
    }
    for (std::size_t i = 0; i < m; ++i) hatA[1](i, i) += 1.0;
    out.hatA = hatA;
    return out;
}

/**
 * Put back on G the rank-one term a shift removed: `ones/m` for a 'one' shift at
 * drift < 1, `tau v e^T` for a 'tau' shift at drift > 1, both for 'dbl'. The
 * reference repeats this switch verbatim in MG1_CR, MG1_FI and MG1_NI.
 */
inline void mg1_unshift(Matrix<double>& G, const ShiftResult& sh, const std::string& shift_type) {
    const std::size_t m = G.rows();
    const bool one = shift_type == "one" || shift_type == "dbl";
    const bool tau = shift_type == "tau" || shift_type == "dbl";
    if (one && sh.drift < 1.0)
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) += 1.0 / static_cast<double>(m);
    if (tau && sh.drift > 1.0)
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) += sh.tau * sh.v[i];
}

// ---------------------------------------------------------------------------
// MG1_EG.m
// ---------------------------------------------------------------------------

/**
 * G in closed form when A0 has rank one. Port of `MG1_EG.m`.
 *
 * `found` is false when the shortcut does not apply, which is the reference's
 * empty return. A rank-one A0 means every down-transition forgets the phase it
 * came from, so G is the same rank-one matrix `e beta` in the recurrent case;
 * this is not an approximation and it is why an M/M/1-shaped input never
 * enters cyclic reduction at all.
 */
inline Matrix<double> mg1_eg(const Blocks& Ain, bool& found) {
    found = false;
    Blocks A = Ain;
    const std::size_t m = A[0].rows();
    const std::size_t dega = A.size() - 1;
    const Drift d = mg1_drift(A);

    if (matrix_rank(A[0]) != 1) return Matrix<double>();

    if (d.value < 1.0) {
        // A0 = alpha beta: G = e beta with beta the normalized first nonzero row.
        const std::vector<double> rs = rowsums(A[0]);
        std::size_t first = m;
        for (std::size_t i = 0; i < m; ++i)
            if (rs[i] > 0.0) {
                first = i;
                break;
            }
        if (first == m) return Matrix<double>();
        Matrix<double> G(m, m, 0.0);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) = A[0](first, j) / rs[first];
        found = true;
        return G;
    }
    if (d.value > 1.0) {
        // Transient chain: G through the Ramaswami dual and its caudal value.
        Blocks At(A.size(), Matrix<double>(m, m, 0.0));
        for (std::size_t b = 0; b < A.size(); ++b)
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j)
                    At[b](i, j) = A[b](j, i) * d.theta[j] / d.theta[i];
        const double etahat = gim1_caudal(At);
        Matrix<double> temp = At[dega];
        for (std::size_t i = dega; i-- > 1;) temp = madd(mscale(temp, etahat), At[i]);
        Matrix<double> M = matmul(At[0], inverse(msub(eye<double>(m), temp)));
        Matrix<double> G(m, m, 0.0);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) = M(j, i) * d.theta[j] / d.theta[i];
        found = true;
        return G;
    }
    return Matrix<double>();
}

// ---------------------------------------------------------------------------
// MG1_CR.m
// ---------------------------------------------------------------------------

/** Options of `MG1_CR`, with the reference's defaults. */
struct Mg1CrOptions {
    std::string mode = "ShiftPWCR";  ///< 'ShiftPWCR' or 'PWCR'
    std::string shift_type = "one";
    std::size_t max_num_it = 50;
    std::size_t max_num_root = 2048;
    double epsilon = 1e-16;
};

namespace cr_detail {

/** Blockwise DFT along the block index, MATLAB's fft over the block sequence. */
inline std::vector<Matrix<std::complex<double>>> block_dft(const Blocks& blk, std::size_t n,
                                                           std::size_t use, bool inverse) {
    const std::size_t m = blk[0].rows();
    std::vector<Matrix<std::complex<double>>> out(
        n, Matrix<std::complex<double>>(m, m, std::complex<double>(0.0, 0.0)));
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            for (std::size_t k = 0; k < n; ++k)
                buf[k] = (k < use && k < blk.size()) ? std::complex<double>(blk[k](i, j), 0.0)
                                                     : std::complex<double>(0.0, 0.0);
            dft(buf, inverse);
            for (std::size_t k = 0; k < n; ++k) out[k](i, j) = buf[k];
        }
    return out;
}

/** The inverse of the above, keeping the real part as `real(ifft(...))` does. */
inline Blocks block_idft_real(const std::vector<Matrix<std::complex<double>>>& F) {
    const std::size_t n = F.size(), m = F[0].rows();
    Blocks out(n, Matrix<double>(m, m, 0.0));
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            for (std::size_t k = 0; k < n; ++k) buf[k] = F[k](i, j);
            dft(buf, true);
            for (std::size_t k = 0; k < n; ++k) out[k](i, j) = buf[k].real();
        }
    return out;
}

inline Matrix<std::complex<double>> cmul(const Matrix<std::complex<double>>& A,
                                         const Matrix<std::complex<double>>& B) {
    return matmul(A, B);
}

inline Matrix<std::complex<double>> cinv_i_minus(const Matrix<std::complex<double>>& A) {
    const std::size_t m = A.rows();
    Matrix<std::complex<double>> M(m, m, std::complex<double>(0.0, 0.0));
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j)
            M(i, j) = (i == j ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0)) -
                      A(i, j);
    return inverse(M);
}

/** The tail norm the reference measures, over blocks deg/2 .. deg-1. */
inline double tail_norm(const Blocks& blk) {
    const std::size_t deg = blk.size();
    // MATLAB's `for i=deg/2:deg-1` starts at a HALF-INTEGER when deg is odd,
    // and a half-integer index never matches, so the loop body runs from
    // ceil(deg/2). Reproduced, or a degree-1 sequence would be measured here
    // where the reference measures nothing.
    const std::size_t start = (deg + 1) / 2;
    double best = 0.0;
    for (std::size_t i = start; i + 1 <= deg && i < deg; ++i) best = std::max(best, inf_norm(blk[i]));
    return best;
}

/** The even-subscript blocks A0, A2, A4, ... of a sequence. */
inline Blocks even_blocks(const Blocks& b) {
    Blocks out;
    for (std::size_t i = 0; i < b.size(); i += 2) out.push_back(b[i]);
    return out;
}

/** The odd-subscript blocks A1, A3, A5, ... of a sequence. */
inline Blocks odd_blocks(const Blocks& b) {
    Blocks out;
    for (std::size_t i = 1; i < b.size(); i += 2) out.push_back(b[i]);
    return out;
}

}  // namespace cr_detail

/**
 * Cyclic reduction for M/G/1-type Markov chains [Bini, Meini]. Port of
 * `MG1_CR.m`, default mode 'ShiftPWCR' with ShiftType 'one'.
 *
 * WHAT THE ALGORITHM DOES, since the transcription is otherwise opaque. One
 * step of cyclic reduction eliminates every odd level of the chain and leaves
 * a chain of the same M/G/1-type shape on the even ones, so the level
 * distance halves per iteration and the iteration converges quadratically.
 * Doing that on the block sequences directly is a polynomial composition; the
 * reference instead evaluates the four sequences at the (nj+1)-th roots of
 * unity, does the composition POINT-WISE (a small dense inverse per root), and
 * interpolates back with an inverse transform -- the "point-wise" in PWCR. The
 * number of roots doubles until the interpolated tail is below (nj+1) eps,
 * which is the reference's own accuracy control and the reason MaxNumRoot
 * exists.
 *
 * Everything runs on the TRANSPOSED blocks, as the reference does after
 * `D=D'`, and the final G is transposed back.
 */
inline Matrix<double> mg1_cr(const Blocks& Ain, const Mg1CrOptions& opts = Mg1CrOptions()) {
    if (Ain.empty()) throw InputError("MG1_CR: empty block sequence");
    const std::size_t m = Ain[0].rows();
    if (opts.mode != "ShiftPWCR" && opts.mode != "PWCR")
        throw UnsupportedError("MG1_CR: Mode '" + opts.mode +
                               "' is not supported; the reference offers 'PWCR' and 'ShiftPWCR'");

    bool eg_found = false;
    const Matrix<double> Geg = mg1_eg(Ain, eg_found);
    if (eg_found) return Geg;

    Blocks A = Ain;
    ShiftResult sh;
    if (opts.mode == "ShiftPWCR") {
        sh = mg1_shifts(A, opts.shift_type);
        A = sh.hatA;
    }

    // D = A', padded with zero blocks to 2^(1+floor(log2(maxd)))+1 of them.
    const std::size_t maxd = A.size() - 1;
    if (maxd == 0) throw InputError("MG1_CR: the sequence needs at least two blocks");
    std::size_t target = 1;
    while (target < maxd) target <<= 1;  // 2^ceil(log2(maxd))
    if (target == maxd) target <<= 1;    // 2^(1+floor(log2(maxd))) when maxd is a power of two
    target += 1;
    Blocks D(target, Matrix<double>(m, m, 0.0));
    for (std::size_t b = 0; b <= maxd; ++b) D[b] = A[b].transpose();

    Blocks Aeven = cr_detail::even_blocks(D);
    Blocks Aodd = cr_detail::odd_blocks(D);
    Blocks Ahatodd(Aeven.begin() + 1, Aeven.end());
    Ahatodd.push_back(D.back());
    Blocks Ahateven = Aodd;

    Matrix<double> Rj = D[1];
    for (std::size_t i = 2; i < D.size(); ++i) Rj = madd(Rj, D[i]);
    Rj = matmul(D[0], inverse(msub(eye<double>(m), Rj)));

    Matrix<double> G(m, m, 0.0);
    Blocks Anew, Ahatnew;
    std::size_t numit = 0;
    while (numit < opts.max_num_it) {
        ++numit;
        std::size_t nj = Aodd.size() - 1;
        double nAnew = 0.0, nAhatnew = 0.0;

        if (nj > 0) {
            const std::size_t n = nj + 1;
            const std::vector<Matrix<std::complex<double>>> T1 =
                cr_detail::block_dft(Aodd, n, n, false);
            const std::vector<Matrix<std::complex<double>>> T2 =
                cr_detail::block_dft(Aeven, n, n, false);
            const std::vector<Matrix<std::complex<double>>> T3 =
                cr_detail::block_dft(Ahatodd, n, n, false);
            const std::vector<Matrix<std::complex<double>>> T4 =
                cr_detail::block_dft(Ahateven, n, n, false);
            std::vector<Matrix<std::complex<double>>> Ah(n), An(n);
            const double pi = 3.14159265358979323846;
            for (std::size_t c = 0; c < n; ++c) {
                const Matrix<std::complex<double>> W = cr_detail::cinv_i_minus(T1[c]);
                Ah[c] = madd(T4[c], cr_detail::cmul(cr_detail::cmul(T2[c], W), T3[c]));
                const double ang = -2.0 * pi * static_cast<double>(c) / static_cast<double>(n);
                const std::complex<double> w(std::cos(ang), std::sin(ang));
                Matrix<std::complex<double>> first(m, m, std::complex<double>(0.0, 0.0));
                for (std::size_t i = 0; i < m; ++i)
                    for (std::size_t j = 0; j < m; ++j) first(i, j) = w * T1[c](i, j);
                An[c] = madd(first, cr_detail::cmul(cr_detail::cmul(T2[c], W), T2[c]));
            }
            Ahatnew = cr_detail::block_idft_real(Ah);
            Anew = cr_detail::block_idft_real(An);
        } else {
            const Matrix<double> temp =
                matmul(Aeven[0], inverse(msub(eye<double>(m), Aodd[0])));
            Ahatnew.assign(1, madd(Ahateven[0], matmul(temp, Ahatodd[0])));
            Anew.clear();
            Anew.push_back(matmul(temp, Aeven[0]));
            Anew.push_back(Aodd[0]);
        }

        nAnew = cr_detail::tail_norm(Anew);
        nAhatnew = cr_detail::tail_norm(Ahatnew);

        // Double the number of roots until the interpolated tail is negligible.
        while ((nAnew > static_cast<double>(nj + 1) * opts.epsilon ||
                nAhatnew > static_cast<double>(nj + 1) * opts.epsilon) &&
               nj + 1 < opts.max_num_root) {
            nj = 2 * (nj + 1) - 1;
            const std::size_t n = nj + 1;
            const std::size_t stopv = std::min(n, Aodd.size());
            const std::vector<Matrix<std::complex<double>>> T1 =
                cr_detail::block_dft(Aodd, n, stopv, false);
            const std::vector<Matrix<std::complex<double>>> T2 =
                cr_detail::block_dft(Aeven, n, stopv, false);
            const std::vector<Matrix<std::complex<double>>> T3 =
                cr_detail::block_dft(Ahatodd, n, stopv, false);
            const std::vector<Matrix<std::complex<double>>> T4 =
                cr_detail::block_dft(Ahateven, n, stopv, false);
            std::vector<Matrix<std::complex<double>>> Ah(n), An(n);
            const double pi = 3.14159265358979323846;
            for (std::size_t c = 0; c < n; ++c) {
                const Matrix<std::complex<double>> W = cr_detail::cinv_i_minus(T1[c]);
                Ah[c] = madd(T4[c], cr_detail::cmul(cr_detail::cmul(T2[c], W), T3[c]));
                const double ang = -2.0 * pi * static_cast<double>(c) / static_cast<double>(n);
                const std::complex<double> w(std::cos(ang), std::sin(ang));
                Matrix<std::complex<double>> first(m, m, std::complex<double>(0.0, 0.0));
                for (std::size_t i = 0; i < m; ++i)
                    for (std::size_t j = 0; j < m; ++j) first(i, j) = w * T1[c](i, j);
                An[c] = madd(first, cr_detail::cmul(cr_detail::cmul(T2[c], W), T2[c]));
            }
            Ahatnew = cr_detail::block_idft_real(Ah);
            Anew = cr_detail::block_idft_real(An);
            nAnew = cr_detail::tail_norm(Anew);
            nAhatnew = cr_detail::tail_norm(Ahatnew);
        }

        if (nj > 1) {
            const std::size_t keep = (nj + 1) / 2;
            Anew.resize(std::min(keep, Anew.size()));
            Ahatnew.resize(std::min(keep, Ahatnew.size()));
        }

        Aeven = cr_detail::even_blocks(Anew);
        Aodd = cr_detail::odd_blocks(Anew);
        Ahateven = cr_detail::even_blocks(Ahatnew);
        Ahatodd = cr_detail::odd_blocks(Ahatnew);

        if (opts.mode == "PWCR") {
            Matrix<double> Rnewj = Anew.size() > 1 ? Anew[1] : Matrix<double>(m, m, 0.0);
            for (std::size_t i = 2; i < Anew.size(); ++i) Rnewj = madd(Rnewj, Anew[i]);
            Rnewj = matmul(Anew[0], inverse(msub(eye<double>(m), Rnewj)));
            const Matrix<double> U =
                Anew.size() > 1
                    ? msub(eye<double>(m),
                           matmul(Anew[0], inverse(msub(eye<double>(m), Anew[1]))))
                    : eye<double>(m);
            if (max_abs_diff(Rj, Rnewj) < opts.epsilon || max_col_sum(U) < opts.epsilon) {
                G = Ahatnew[0];
                for (std::size_t i = 1; i < Ahatnew.size(); ++i)
                    G = madd(G, matmul(Rnewj, Ahatnew[i]));
                G = matmul(D[0], inverse(msub(eye<double>(m), G)));
                break;
            }
            Rj = Rnewj;
            double tail_sum = 0.0;
            for (std::size_t i = 1; i < Ahatnew.size(); ++i) tail_sum += Ahatnew[i].sum();
            const double sv = svd_values(Anew[0]).empty() ? 0.0 : svd_values(Anew[0])[0];
            const Matrix<double> V =
                msub(eye<double>(m), matmul(D[0], inverse(msub(eye<double>(m), Ahatnew[0]))));
            if (sv < opts.epsilon || tail_sum < opts.epsilon || max_col_sum(V) < opts.epsilon) {
                G = matmul(D[0], inverse(msub(eye<double>(m), Ahatnew[0])));
                break;
            }
        } else {
            const Matrix<double> Gold = G;
            G = matmul(D[0], inverse(msub(eye<double>(m), Ahatnew[0])));
            if (inf_norm(msub(G, Gold)) < opts.epsilon || inf_norm(Ahatnew, 1) < opts.epsilon)
                break;
        }
    }
    if (numit == opts.max_num_it && !Ahatnew.empty())
        G = matmul(D[0], inverse(msub(eye<double>(m), Ahatnew[0])));

    G = G.transpose();

    // Undo the shift: put back the rank-one term it removed from G.
    if (opts.mode == "ShiftPWCR") mg1_unshift(G, sh, opts.shift_type);
    return G;
}

// ---------------------------------------------------------------------------
// MG1_FI.m
// ---------------------------------------------------------------------------

/** Options of `MG1_FI`, with the reference's defaults. */
struct Mg1FiOptions {
    std::string mode = "U-Based";  ///< 'Natural', 'Traditional', 'U-Based', or 'Shift<Mode>'
    std::string shift_type = "one";
    std::size_t max_num_it = 10000;
    double tol = 1e-14;
};

/**
 * Functional iterations for M/G/1-type Markov chains [Neuts]. Port of
 * `MG1_FI.m` for the three modes and the shift variants; the reference's
 * `NonZeroBlocks` option is not exposed, because it changes only which
 * products are skipped when some A_i vanish and converges to the same G.
 *
 * 'U-Based' is the default and the one `GIM1_R(...,'FI')` uses: it solves
 * G = (I - sum_{j>=1} A_j G^{j-1})^{-1} A0, which is the U-based iteration and
 * converges monotonically from below to the minimal nonnegative solution.
 */
inline Matrix<double> mg1_fi(const Blocks& Ain, const Mg1FiOptions& opts = Mg1FiOptions()) {
    const std::size_t m = Ain[0].rows();
    const std::size_t maxd = Ain.size() - 1;

    bool eg_found = false;
    const Matrix<double> Geg = mg1_eg(Ain, eg_found);
    if (eg_found) return Geg;

    Blocks A = Ain;
    const bool shifted = opts.mode.find("Shift") != std::string::npos;
    ShiftResult sh;
    if (shifted) {
        sh = mg1_shifts(A, opts.shift_type);
        A = sh.hatA;
    }

    const bool natural = opts.mode.find("Natural") != std::string::npos;
    const bool traditional = opts.mode.find("Traditional") != std::string::npos;
    const bool ubased = opts.mode.find("U-Based") != std::string::npos;
    if (!natural && !traditional && !ubased)
        throw UnsupportedError("MG1_FI: Mode '" + opts.mode + "' is not supported");

    Matrix<double> G(m, m, 0.0);
    double check = 1.0;
    std::size_t numit = 0;
    while (check > opts.tol && numit < opts.max_num_it) {
        const Matrix<double> Gold = G;
        if (natural) {
            G = A[maxd];
            for (std::size_t j = maxd; j-- > 0;) G = madd(A[j], matmul(G, Gold));
        } else if (traditional) {
            G = A[maxd];
            for (std::size_t j = maxd; j-- > 2;) G = madd(A[j], matmul(G, Gold));
            G = madd(A[0], matmul(G, matmul(Gold, Gold)));
            G = matmul(inverse(msub(eye<double>(m), A[1])), G);
        } else {
            G = A[maxd];
            for (std::size_t j = maxd; j-- > 1;) G = madd(A[j], matmul(G, Gold));
            G = matmul(inverse(msub(eye<double>(m), G)), A[0]);
        }
        check = inf_norm(msub(G, Gold));
        ++numit;
    }

    if (shifted) mg1_unshift(G, sh, opts.shift_type);
    return G;
}

// ---------------------------------------------------------------------------
// Dense block helpers for MG1_NI / MG1_RR / MG1_IS
// ---------------------------------------------------------------------------

namespace mg1x_detail {

/** A(r0:r0+nr, c0:c0+nc), half-open, 0-based. */
inline Matrix<double> sub(const Matrix<double>& A, std::size_t r0, std::size_t nr, std::size_t c0,
                          std::size_t nc) {
    Matrix<double> S(nr, nc, 0.0);
    for (std::size_t i = 0; i < nr; ++i)
        for (std::size_t j = 0; j < nc; ++j) S(i, j) = A(r0 + i, c0 + j);
    return S;
}

/** A(r0:, c0:) = S. */
inline void put(Matrix<double>& A, std::size_t r0, std::size_t c0, const Matrix<double>& S) {
    for (std::size_t i = 0; i < S.rows(); ++i)
        for (std::size_t j = 0; j < S.cols(); ++j) A(r0 + i, c0 + j) = S(i, j);
}

inline Matrix<double> trans(const Matrix<double>& A) {
    Matrix<double> T(A.cols(), A.rows(), 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) T(j, i) = A(i, j);
    return T;
}

/** [A B], horizontally. */
inline Matrix<double> hjoin(const Matrix<double>& A, const Matrix<double>& B) {
    if (A.rows() != B.rows()) throw InputError("smc: horizontal join with mismatched rows");
    Matrix<double> C(A.rows(), A.cols() + B.cols(), 0.0);
    put(C, 0, 0, A);
    put(C, 0, A.cols(), B);
    return C;
}

/** [A; B], vertically. */
inline Matrix<double> vjoin(const Matrix<double>& A, const Matrix<double>& B) {
    if (A.cols() != B.cols()) throw InputError("smc: vertical join with mismatched columns");
    Matrix<double> C(A.rows() + B.rows(), A.cols(), 0.0);
    put(C, 0, 0, A);
    put(C, A.rows(), 0, B);
    return C;
}

/** Z \ b for a square Z. */
inline std::vector<double> lsolve(const Matrix<double>& Z, const std::vector<double>& b) {
    return solve(Z, b);
}

/**
 * Economy QR, A = Q R with Q (r x k) orthonormal columns and R (k x c), k = min(r, c), by
 * Householder reflections. The sign convention differs from LAPACK's, which MG1_RR cannot see:
 * its generators enter B only through a product invariant under an orthogonal change of basis.
 */
inline void qr_econ(const Matrix<double>& A, Matrix<double>& Q, Matrix<double>& R) {
    const std::size_t r = A.rows(), c = A.cols(), k = std::min(r, c);
    Matrix<double> W = A;
    std::vector<std::vector<double>> vs;
    for (std::size_t j = 0; j < k; ++j) {
        double nrm = 0.0;
        for (std::size_t i = j; i < r; ++i) nrm += W(i, j) * W(i, j);
        nrm = std::sqrt(nrm);
        std::vector<double> v(r, 0.0);
        if (nrm > 0.0) {
            const double alpha = W(j, j) > 0.0 ? -nrm : nrm;
            for (std::size_t i = j; i < r; ++i) v[i] = W(i, j);
            v[j] -= alpha;
            double vn = 0.0;
            for (std::size_t i = j; i < r; ++i) vn += v[i] * v[i];
            if (vn > 0.0) {
                for (std::size_t col = j; col < c; ++col) {
                    double s = 0.0;
                    for (std::size_t i = j; i < r; ++i) s += v[i] * W(i, col);
                    s = 2.0 * s / vn;
                    for (std::size_t i = j; i < r; ++i) W(i, col) -= s * v[i];
                }
                const double vnr = std::sqrt(vn);
                for (std::size_t i = j; i < r; ++i) v[i] /= vnr;
            } else {
                std::fill(v.begin(), v.end(), 0.0);
            }
        }
        vs.push_back(v);
    }
    R = Matrix<double>(k, c, 0.0);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = i; j < c; ++j) R(i, j) = W(i, j);
    // Q = H_0 H_1 ... H_{k-1} applied to the first k columns of the identity.
    Q = Matrix<double>(r, k, 0.0);
    for (std::size_t i = 0; i < k; ++i) Q(i, i) = 1.0;
    for (std::size_t j = k; j-- > 0;) {
        const std::vector<double>& v = vs[j];
        for (std::size_t col = 0; col < k; ++col) {
            double s = 0.0;
            for (std::size_t i = j; i < r; ++i) s += v[i] * Q(i, col);
            s *= 2.0;
            for (std::size_t i = j; i < r; ++i) Q(i, col) -= s * v[i];
        }
    }
}

}  // namespace mg1x_detail

// ---------------------------------------------------------------------------
// solveSylvPowersDirectSum.m, solveSylvPowersRealSchur_FW.m
// ---------------------------------------------------------------------------

/**
 * Solve sum_{j=1}^N B_j Y A^{j-1} = C directly, through the Kronecker form of vec(Y).
 * Port of `solveSylvPowersDirectSum.m`; `B` holds the N blocks B_1..B_N (n x n), A is m x m.
 */
inline Matrix<double> sylv_powers_direct(const Matrix<double>& A, const Blocks& B,
                                         const Matrix<double>& C) {
    const std::size_t m = A.rows(), n = C.rows(), N = B.size();
    const Matrix<double> At = mg1x_detail::trans(A);
    Matrix<double> P = eye<double>(m);  // (A')^(j-1)
    Matrix<double> Z(m * n, m * n, 0.0);
    for (std::size_t j = 0; j < N; ++j) {
        for (std::size_t a = 0; a < m; ++a)
            for (std::size_t b = 0; b < m; ++b) {
                const double p = P(a, b);
                if (p == 0.0) continue;
                for (std::size_t k = 0; k < n; ++k)
                    for (std::size_t l = 0; l < n; ++l) Z(a * n + k, b * n + l) += p * B[j](k, l);
            }
        P = matmul(P, At);
    }
    std::vector<double> c(m * n);
    for (std::size_t col = 0; col < m; ++col)
        for (std::size_t k = 0; k < n; ++k) c[col * n + k] = C(k, col);
    const std::vector<double> y = mg1x_detail::lsolve(Z, c);
    Matrix<double> Y(n, m, 0.0);
    for (std::size_t col = 0; col < m; ++col)
        for (std::size_t k = 0; k < n; ++k) Y(k, col) = y[col * n + k];
    return Y;
}

/**
 * Solve sum_{j=1}^N B_j Y A^{j-1} = C by a real Schur form A = U T U': with Y = X U' the
 * system becomes sum_j B_j X T^{j-1} = C U, and T quasi-triangular makes it a forward
 * substitution over the columns of X, one n x n solve per 1 x 1 block and one 2n x 2n solve
 * per 2 x 2 block. Port of `solveSylvPowersRealSchur_FW.m`, including its 1e-13 test on the
 * subdiagonal that tells the two block shapes apart.
 */
inline Matrix<double> sylv_powers_real_schur(const Matrix<double>& A, const Blocks& B,
                                             const Matrix<double>& C) {
    const std::size_t m = A.rows(), n = C.rows(), N = B.size();
    const RealSchur sf = schur_decomposition(A);
    const Matrix<double>& T = sf.T;
    const Matrix<double> E = matmul(C, sf.Z);
    Blocks Tp(N, eye<double>(m));
    for (std::size_t j = 1; j < N; ++j) Tp[j] = matmul(Tp[j - 1], T);
    // S(a,b) = sum_j B_j (T^{j-1})(a,b)
    auto S = [&](std::size_t a, std::size_t b) {
        Matrix<double> out(n, n, 0.0);
        for (std::size_t j = 0; j < N; ++j) {
            const double t = Tp[j](a, b);
            if (t == 0.0) continue;
            for (std::size_t k = 0; k < n; ++k)
                for (std::size_t l = 0; l < n; ++l) out(k, l) += t * B[j](k, l);
        }
        return out;
    };
    // sum_{l<k} S(l,col) X(:,l)
    Matrix<double> X(n, m, 0.0);
    auto W = [&](std::size_t k, std::size_t col) {
        std::vector<double> w(n, 0.0);
        for (std::size_t l = 0; l < k; ++l) {
            const Matrix<double> Sl = S(l, col);
            for (std::size_t a = 0; a < n; ++a)
                for (std::size_t b = 0; b < n; ++b) w[a] += Sl(a, b) * X(b, l);
        }
        return w;
    };
    const double epsilon = 10e-14;
    std::size_t k = 0;
    while (k < m) {
        if (k + 1 == m || std::abs(T(k + 1, k)) < epsilon) {
            const Matrix<double> Z = S(k, k);
            const std::vector<double> w = W(k, k);
            std::vector<double> rhs(n);
            for (std::size_t a = 0; a < n; ++a) rhs[a] = E(a, k) - w[a];
            const std::vector<double> x = mg1x_detail::lsolve(Z, rhs);
            for (std::size_t a = 0; a < n; ++a) X(a, k) = x[a];
            k += 1;
        } else {
            Matrix<double> Z(2 * n, 2 * n, 0.0);
            mg1x_detail::put(Z, 0, 0, S(k, k));
            mg1x_detail::put(Z, 0, n, S(k + 1, k));
            mg1x_detail::put(Z, n, 0, S(k, k + 1));
            mg1x_detail::put(Z, n, n, S(k + 1, k + 1));
            const std::vector<double> w0 = W(k, k), w1 = W(k, k + 1);
            std::vector<double> rhs(2 * n);
            for (std::size_t a = 0; a < n; ++a) {
                rhs[a] = E(a, k) - w0[a];
                rhs[n + a] = E(a, k + 1) - w1[a];
            }
            const std::vector<double> x = mg1x_detail::lsolve(Z, rhs);
            for (std::size_t a = 0; a < n; ++a) {
                X(a, k) = x[a];
                X(a, k + 1) = x[n + a];
            }
            k += 2;
        }
    }
    return matmul(X, mg1x_detail::trans(sf.Z));
}

// ---------------------------------------------------------------------------
// MG1_NI.m
// ---------------------------------------------------------------------------

/** Options of `MG1_NI`, with the reference's defaults. */
struct Mg1NiOptions {
    /// 'DirectSum', 'RealSchur' or 'ComplexSchur', each optionally with the 'Shift' suffix
    std::string mode = "RealSchurShift";
    std::size_t max_num_it = 50;
    std::string shift_type = "one";
    double epsilon = 1e-14;
};

/**
 * Newton iteration for M/G/1-type Markov chains. Port of `MG1_NI.m`.
 *
 * Each step linearizes G = sum_i A_i G^i at the current G and solves the resulting
 * Sylvester-power equation sum_j B_j Y G^{j-1} = G - B_0 for the correction Y, so the
 * iteration converges quadratically; the shift variants run it on the sequence with the
 * unit root moved to zero (`mg1_shifts`) and put the removed rank-one term back at the end.
 *
 * 'ComplexSchur' asks the reference for a complex Schur form of G. It solves the SAME linear
 * system as 'RealSchur', whose real quasi-triangular form this port uses for both, so the two
 * modes differ only in rounding: there is no complex Schur factorization in this tree.
 */
inline Matrix<double> mg1_ni(const Blocks& Ain, const Mg1NiOptions& opts = Mg1NiOptions()) {
    if (Ain.size() < 2) throw InputError("MG1_NI: the sequence needs at least two blocks");
    std::string base = opts.mode;
    const bool shifted = base.size() > 5 && base.compare(base.size() - 5, 5, "Shift") == 0;
    if (shifted) base = base.substr(0, base.size() - 5);
    if (base != "DirectSum" && base != "RealSchur" && base != "ComplexSchur")
        throw InputError("MG1_NI: Mode '" + opts.mode + "' is not one of DirectSum, RealSchur, "
                         "ComplexSchur, each optionally with the Shift suffix");

    bool eg_found = false;
    const Matrix<double> Geg = mg1_eg(Ain, eg_found);
    if (eg_found) return Geg;

    Blocks A = Ain;
    ShiftResult sh;
    if (shifted) {
        sh = mg1_shifts(A, opts.shift_type);
        A = sh.hatA;
    }

    const std::size_t m = A[0].rows();
    const std::size_t N = A.size() - 1;
    Matrix<double> G(m, m, 0.0);
    double check = 1.0;
    std::size_t numit = 0;
    while (check > opts.epsilon && numit < opts.max_num_it) {
        const Matrix<double> Gold = G;
        Blocks Bk(N + 1);
        Bk[N] = A[N];
        for (std::size_t i = N; i-- > 0;) Bk[i] = madd(A[i], matmul(Bk[i + 1], G));
        const Matrix<double> C = msub(G, Bk[0]);
        Blocks Bs(Bk.begin() + 1, Bk.end());
        for (std::size_t i = 0; i < m; ++i) Bs[0](i, i) -= 1.0;
        const Matrix<double> Y =
            base == "DirectSum" ? sylv_powers_direct(G, Bs, C) : sylv_powers_real_schur(G, Bs, C);
        G = madd(G, Y);
        check = inf_norm(msub(G, Gold));
        ++numit;
    }

    if (shifted) mg1_unshift(G, sh, opts.shift_type);
    return G;
}

// ---------------------------------------------------------------------------
// MG1_RR.m (with MG1_RR_Btemp.m and MG1_RR_tempB.m)
// ---------------------------------------------------------------------------

/** Options of `MG1_RR`, with the reference's defaults. */
struct Mg1RrOptions {
    std::string mode = "Direct";  ///< 'Direct', 'DispStruct' or 'DispStructFFT'
    std::size_t max_num_it = 50;
};

namespace rr_detail {

/**
 * The displacement representation B = L(b) + L(c1) L(Z r1)' + L(c2) L(Z r2)' of the
 * Ramaswami-reduction matrix, with L(x) the block lower-triangular Toeplitz matrix whose
 * first block column is x (N*m x m) and Z the block down-shift. Only the generators are
 * stored; `left` and `right` are `MG1_RR_tempB.m` and `MG1_RR_Btemp.m`.
 */
struct Disp {
    std::size_t N, m;
    Matrix<double> b, c1, c2, r1, r2;

    /** Block k (0-based) of a stacked N*m x m generator. */
    Matrix<double> blk(const Matrix<double>& x, std::size_t k) const {
        return mg1x_detail::sub(x, k * m, m, 0, m);
    }

    /** temp * B, `MG1_RR_tempB.m`. */
    Matrix<double> left(const Matrix<double>& temp) const {
        const std::size_t hd = temp.rows();
        Matrix<double> out(hd, N * m, 0.0);
        // temp * L(x): block column l is sum_{i>=l} temp_i x_{i-l}
        auto times_L = [&](const Matrix<double>& x) {
            Matrix<double> r(hd, N * m, 0.0);
            for (std::size_t l = 0; l < N; ++l) {
                Matrix<double> acc(hd, m, 0.0);
                for (std::size_t i = l; i < N; ++i)
                    acc = madd(acc, matmul(mg1x_detail::sub(temp, 0, hd, i * m, m), blk(x, i - l)));
                mg1x_detail::put(r, 0, l * m, acc);
            }
            return r;
        };
        // (temp L(c)) L(Z r)': L(Z r)' is block upper triangular with block (i,l) =
        // (Zr)_{l-i}' = r_{l-i-1}' for l > i, zero otherwise.
        auto times_LZrT = [&](const Matrix<double>& P, const Matrix<double>& r) {
            Matrix<double> q(hd, N * m, 0.0);
            for (std::size_t l = 1; l < N; ++l) {
                Matrix<double> acc(hd, m, 0.0);
                for (std::size_t i = 0; i < l; ++i)
                    acc = madd(acc, matmul(mg1x_detail::sub(P, 0, hd, i * m, m),
                                           mg1x_detail::trans(blk(r, l - i - 1))));
                mg1x_detail::put(q, 0, l * m, acc);
            }
            return q;
        };
        out = madd(out, times_LZrT(times_L(c1), r1));
        out = madd(out, times_LZrT(times_L(c2), r2));
        out = madd(out, times_L(b));
        return out;
    }

    /** B * temp, `MG1_RR_Btemp.m`. */
    Matrix<double> right(const Matrix<double>& temp) const {
        const std::size_t hd = temp.cols();
        // L(Z r)' temp: block row i is sum_{l>i} r_{l-i-1}' temp_l
        auto LZrT_times = [&](const Matrix<double>& r) {
            Matrix<double> q(N * m, hd, 0.0);
            for (std::size_t i = 0; i + 1 < N; ++i) {
                Matrix<double> acc(m, hd, 0.0);
                for (std::size_t l = i + 1; l < N; ++l)
                    acc = madd(acc, matmul(mg1x_detail::trans(blk(r, l - i - 1)),
                                           mg1x_detail::sub(temp, l * m, m, 0, hd)));
                mg1x_detail::put(q, i * m, 0, acc);
            }
            return q;
        };
        // L(x) P: block row i is sum_{l<=i} x_{i-l} P_l
        auto L_times = [&](const Matrix<double>& x, const Matrix<double>& P) {
            Matrix<double> q(N * m, hd, 0.0);
            for (std::size_t i = 0; i < N; ++i) {
                Matrix<double> acc(m, hd, 0.0);
                for (std::size_t l = 0; l <= i; ++l)
                    acc = madd(acc, matmul(blk(x, i - l), mg1x_detail::sub(P, l * m, m, 0, hd)));
                mg1x_detail::put(q, i * m, 0, acc);
            }
            return q;
        };
        Matrix<double> out = L_times(c1, LZrT_times(r1));
        out = madd(out, L_times(c2, LZrT_times(r2)));
        out = madd(out, L_times(b, temp));
        return out;
    }
};

inline double sum_all(const Matrix<double>& A) {
    double s = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) s += A(i, j);
    return s;
}

}  // namespace rr_detail

/**
 * Ramaswami reduction for M/G/1-type Markov chains [Bini, Meini, Ramaswami]. Port of
 * `MG1_RR.m`.
 *
 * The chain is reduced to a QBD whose level is N blocks wide, and cyclic reduction on that
 * QBD is run in the Sherman-Morrison-Woodbury form of the reference, carrying only the
 * rank-m corrections `uhat` and `vT`; G is then read off by the formula at the end of
 * Section 3 of the paper, G = (I - vT uhat)^{-1} A0.
 *
 * 'Direct' (the default) stores the N*m x N*m matrix B outright. 'DispStruct' stores it
 * through its displacement generators (`rr_detail::Disp`), compressed to rank 2m by a QR of
 * each side and an SVD of the small core at every step. 'DispStructFFT' is the reference's
 * FFT evaluation of the same block-Toeplitz products; this port evaluates those products
 * directly, so it returns the 'DispStruct' answer, which the FFT reproduces up to rounding.
 */
inline Matrix<double> mg1_rr(const Blocks& D, const Mg1RrOptions& opts = Mg1RrOptions()) {
    using namespace mg1x_detail;
    if (D.size() < 2) throw InputError("MG1_RR: the sequence needs at least two blocks");
    if (opts.mode != "Direct" && opts.mode != "DispStruct" && opts.mode != "DispStructFFT")
        throw InputError("MG1_RR: Mode '" + opts.mode +
                         "' is not one of Direct, DispStruct, DispStructFFT");

    bool eg_found = false;
    const Matrix<double> Geg = mg1_eg(D, eg_found);
    if (eg_found) return Geg;

    const std::size_t m = D[0].rows();
    const std::size_t N = D.size() - 1;
    const Matrix<double> I = eye<double>(m);
    Matrix<double> D0 = D[0];
    Matrix<double> vT(m, N * m, 0.0);
    for (std::size_t k = 1; k <= N; ++k) put(vT, 0, (k - 1) * m, D[k]);
    const Matrix<double> vT_orig = vT;
    Matrix<double> uhat(N * m, m, 0.0);
    put(uhat, 0, 0, I);

    // [0; uhat(m+1:end,:)]
    auto ZZTu = [&](const Matrix<double>& u) {
        Matrix<double> z = u;
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) z(i, j) = 0.0;
        return z;
    };

    if (opts.mode == "Direct") {
        Matrix<double> B(N * m, N * m, 0.0);
        for (std::size_t i = 0; i + m < N * m; ++i) B(m + i, i) = 1.0;
        double check = std::min(inf_norm(B), inf_norm(D0));
        std::size_t numit = 0;
        while (check > 1e-15 && numit < opts.max_num_it) {
            ++numit;
            const Matrix<double> S = inverse(msub(I, matmul(vT, uhat)));
            const Matrix<double> T = matmul(S, matmul(vT, ZZTu(uhat)));
            const Matrix<double> mid = madd(madd(matmul(S, sub(vT, 0, m, 0, m)), T), I);
            const Matrix<double> newD0 = matmul(matmul(D0, mid), D0);
            Matrix<double> newB = matmul(uhat, matmul(S, matmul(vT, B)));
            newB = matmul(B, madd(B, newB));
            const Matrix<double> newuhat = madd(uhat, matmul(B, matmul(matmul(uhat, mid), D0)));
            const Matrix<double> newvT = madd(vT, matmul(matmul(matmul(D0, S), vT), B));
            B = newB;
            vT = newvT;
            uhat = newuhat;
            D0 = newD0;
            check = std::min(inf_norm(B), inf_norm(D0));
        }
    } else {
        if (N < 2)
            throw InputError("MG1_RR: Mode '" + opts.mode +
                             "' needs at least three blocks A0, A1, A2; use Mode 'Direct'");
        rr_detail::Disp B{N, m, Matrix<double>(N * m, m, 0.0), Matrix<double>(N * m, m, 0.0),
                          Matrix<double>(N * m, m, 0.0), Matrix<double>(N * m, m, 0.0),
                          Matrix<double>(N * m, m, 0.0)};
        if (N >= 2) put(B.b, m, 0, I);
        bool first = true;
        double check = std::min(rr_detail::sum_all(B.b), inf_norm(D0));
        std::size_t numit = 0;
        while ((check > 1e-15 || first) && numit < opts.max_num_it) {
            ++numit;
            first = false;
            // step 1
            const Matrix<double> S = inverse(msub(I, matmul(vT, uhat)));
            const Matrix<double> T = matmul(S, matmul(vT, ZZTu(uhat)));
            const Matrix<double> mid = madd(madd(matmul(S, sub(vT, 0, m, 0, m)), T), I);
            // step 2
            const Matrix<double> newD0 = matmul(matmul(D0, mid), D0);
            // step 3
            const Matrix<double> newuhat = madd(uhat, B.right(matmul(matmul(uhat, mid), D0)));
            const Matrix<double> newvT = madd(vT, B.left(matmul(matmul(D0, S), vT)));
            // [S T; S I+T]
            Matrix<double> ST(2 * m, 2 * m, 0.0);
            put(ST, 0, 0, S);
            put(ST, 0, m, T);
            put(ST, m, 0, S);
            put(ST, m, m, madd(I, T));
            // (I - A)^{-1} X through its structure, for X with N*m rows
            auto inv_struct = [&](const Matrix<double>& X) {
                const std::size_t w = X.cols();
                const Matrix<double> t = matmul(ST, vjoin(matmul(vT, X), sub(X, 0, m, 0, w)));
                Matrix<double> out(N * m, w, 0.0);
                put(out, 0, 0, madd(sub(X, 0, m, 0, w), sub(t, 0, m, 0, w)));
                if (N > 1)
                    put(out, m, 0,
                        madd(sub(X, m, (N - 1) * m, 0, w),
                             matmul(sub(uhat, m, (N - 1) * m, 0, m), sub(t, m, m, 0, w))));
                return out;
            };
            // step 4
            const Matrix<double> newb = B.right(inv_struct(B.b));

            // step 5, part 1: H = [e1pZu, e2pZ2u*S, e2pZ2u*T + Z2u]
            Matrix<double> e1pZu(N * m, m, 0.0), Z2u(N * m, m, 0.0), e2pZ2u(N * m, m, 0.0);
            put(e1pZu, 0, 0, I);
            if (N > 1) put(e1pZu, m, 0, sub(uhat, m, (N - 1) * m, 0, m));
            if (N > 1) put(e2pZ2u, m, 0, I);
            if (N > 2) {
                put(Z2u, 2 * m, 0, sub(uhat, m, (N - 2) * m, 0, m));
                put(e2pZ2u, 2 * m, 0, sub(uhat, m, (N - 2) * m, 0, m));
            }
            const Matrix<double> H =
                hjoin(hjoin(e1pZu, matmul(e2pZ2u, S)), madd(matmul(e2pZ2u, T), Z2u));
            // KT = [S*[vT(:,m+1:end) 0]; -vT; [-I 0]]
            Matrix<double> vshift(m, N * m, 0.0);
            if (N > 1) put(vshift, 0, 0, sub(vT, 0, m, m, (N - 1) * m));
            Matrix<double> negI(m, N * m, 0.0);
            for (std::size_t i = 0; i < m; ++i) negI(i, i) = -1.0;
            const Matrix<double> KT =
                vjoin(vjoin(matmul(S, vshift), mscale(vT, -1.0)), negI);
            // W = [c1 c2 B*H B*inv_struct([c1 c2])]
            const Matrix<double> C12 = hjoin(B.c1, B.c2);
            const Matrix<double> W =
                hjoin(hjoin(C12, B.right(H)), B.right(inv_struct(C12)));
            // step 6, part 1
            Matrix<double> Q1, R1;
            qr_econ(W, Q1, R1);
            // YT = [ (row block through Zu) * B ; KT * B ; [r1 r2]' ]
            const Matrix<double> R12t = trans(hjoin(B.r1, B.r2));  // 2m x N*m
            const Matrix<double> Zu = ZZTu(uhat);
            const Matrix<double> tq =
                matmul(hjoin(sub(R12t, 0, 2 * m, 0, m), matmul(R12t, Zu)), ST);  // 2m x 2m
            Matrix<double> temp2 = matmul(sub(tq, 0, 2 * m, 0, m), vT);
            put(temp2, 0, 0, madd(sub(temp2, 0, 2 * m, 0, m), sub(tq, 0, 2 * m, m, m)));
            const Matrix<double> YT =
                vjoin(vjoin(B.left(madd(temp2, R12t)), B.left(KT)), R12t);  // 7m x N*m
            // step 6, part 2
            Matrix<double> Q2, R2;
            qr_econ(trans(YT), Q2, R2);
            // step 7
            const SvdFactors sv = svd_full(matmul(R1, trans(R2)));
            // step 8
            Matrix<double> US(sv.U.rows(), 2 * m, 0.0);
            for (std::size_t i = 0; i < sv.U.rows(); ++i)
                for (std::size_t j = 0; j < 2 * m; ++j) US(i, j) = sv.U(i, j) * sv.s[j];
            const Matrix<double> cc = matmul(Q1, US);
            Matrix<double> V2(sv.Vt.cols(), 2 * m, 0.0);
            for (std::size_t i = 0; i < sv.Vt.cols(); ++i)
                for (std::size_t j = 0; j < 2 * m; ++j) V2(i, j) = sv.Vt(j, i);
            const Matrix<double> rr = matmul(Q2, V2);
            B.c1 = sub(cc, 0, N * m, 0, m);
            B.c2 = sub(cc, 0, N * m, m, m);
            B.r1 = sub(rr, 0, N * m, 0, m);
            B.r2 = sub(rr, 0, N * m, m, m);

            B.b = newb;
            vT = newvT;
            uhat = newuhat;
            D0 = newD0;
            check = std::min(rr_detail::sum_all(B.b), inf_norm(D0));
        }
    }

    // G = (I - [A1 ... AN] uhat)^{-1} A0
    return matmul(inverse(msub(I, matmul(vT_orig, uhat))), D[0]);
}

// ---------------------------------------------------------------------------
// MG1_IS.m
// ---------------------------------------------------------------------------

/** Options of `MG1_IS`, with the reference's defaults. */
struct Mg1IsOptions {
    std::string mode = "MSignBalzer";  ///< 'MSignStandard', 'MSignBalzer' or 'Schur'
    std::size_t max_num_it = 50;
};

/**
 * Invariant subspace method for M/G/1-type Markov chains [Akar, Sohraby]. Port of `MG1_IS.m`.
 *
 * F(z) = zI - A(z) is mapped by the Moebius transform z = (1+s)/(1-s) onto a matrix
 * polynomial H(s) whose companion matrix, after a rank-one correction that moves the unit
 * root off the imaginary axis, has exactly m eigenvalues in the open left half plane. A basis
 * T of that invariant subspace gives G = (T1 + T2)(T1 - T2)^{-1}, invariant to the choice of
 * basis. The matrix-sign modes find the subspace as the range of sign(Z) - I; 'Schur' orders
 * a real Schur form with the left half plane first, as the reference's ordschur 'lhp' does.
 */
inline Matrix<double> mg1_is(const Blocks& D, const Mg1IsOptions& opts = Mg1IsOptions()) {
    using namespace mg1x_detail;
    if (opts.mode != "MSignStandard" && opts.mode != "MSignBalzer" && opts.mode != "Schur")
        throw InputError("MG1_IS: Mode '" + opts.mode +
                         "' is not one of MSignStandard, MSignBalzer, Schur");
    if (D.size() < 3)
        throw InputError("MG1_IS: the sequence needs at least three blocks A0, A1, A2; the "
                         "reference reads rows m+1:2m of an m*max-row subspace basis");

    bool eg_found = false;
    const Matrix<double> Geg = mg1_eg(D, eg_found);
    if (eg_found) return Geg;

    const double epsilon = 1e-14;
    const std::size_t m = D[0].rows();
    const std::size_t f = D.size() - 1;
    const double drift = mg1_drift(D).value - 1.0;
    const double sgn = drift > 0.0 ? 1.0 : (drift < 0.0 ? -1.0 : 0.0);

    // Step 1, F(z) = z - A(z)
    Blocks F(f + 1);
    for (std::size_t i = 0; i <= f; ++i) F[i] = mscale(D[i], -1.0);
    for (std::size_t i = 0; i < m; ++i) F[1](i, i) += 1.0;

    // Step 2, H(s) = sum_i F_i (1-s)^(f-i) (1+s)^i = sum_j H_j s^j
    Blocks H(f + 1, Matrix<double>(m, m, 0.0));
    for (std::size_t i = 0; i <= f; ++i) {
        std::vector<double> c{1.0};
        auto conv = [](const std::vector<double>& a, double s1) {  // a * [1 s1]
            std::vector<double> r(a.size() + 1, 0.0);
            for (std::size_t k = 0; k < a.size(); ++k) {
                r[k] += a[k];
                r[k + 1] += s1 * a[k];
            }
            return r;
        };
        for (std::size_t j = 0; j < f - i; ++j) c = conv(c, -1.0);
        for (std::size_t j = 0; j < i; ++j) c = conv(c, 1.0);
        for (std::size_t j = 0; j <= f; ++j) H[j] = madd(H[j], mscale(F[i], c[j]));
    }

    // Step 3, hatH_i = H_f^{-1} H_i
    const Matrix<double> Hfinv = inverse(H[f]);
    Blocks hatH(f);
    for (std::size_t i = 0; i < f; ++i) hatH[i] = matmul(Hfinv, H[i]);

    // Step 4, y and xT; x0T = [0 1] / [hatH0 e], a least-squares solve as in MATLAB
    const std::size_t mf = m * f;
    std::vector<double> y(mf, 0.0);
    for (std::size_t i = 0; i < m; ++i) y[i] = 1.0;
    Matrix<double> Mt(m + 1, m, 0.0);  // [hatH0 e]'
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) Mt(j, i) = hatH[0](i, j);
        Mt(m, i) = 1.0;
    }
    std::vector<double> rhs(m + 1, 0.0);
    rhs[m] = 1.0;
    const std::vector<double> x0 = lstsq(Mt, rhs).x;
    std::vector<double> xT(mf, 0.0);
    for (std::size_t i = 1; i < f; ++i) {
        const std::vector<double> xi = vecmul(x0, hatH[i]);
        for (std::size_t j = 0; j < m; ++j) xT[(i - 1) * m + j] = xi[j];
    }
    for (std::size_t j = 0; j < m; ++j) xT[(f - 1) * m + j] = x0[j];

    // Step 5, the companion matrix plus the rank-one correction
    Matrix<double> Z(mf, mf, 0.0);
    for (std::size_t i = 1; i < f; ++i) put(Z, (i - 1) * m, i * m, eye<double>(m));
    for (std::size_t i = 0; i < f; ++i) put(Z, (f - 1) * m, i * m, mscale(hatH[i], -1.0));
    double xy = 0.0;
    for (std::size_t k = 0; k < mf; ++k) xy += xT[k] * y[k];
    for (std::size_t a = 0; a < mf; ++a)
        for (std::size_t b = 0; b < mf; ++b) Z(a, b) += sgn * (y[a] / xy) * xT[b];

    Matrix<double> T;
    if (opts.mode != "Schur") {
        // Step 6, matrix sign function, optionally with Balzer's determinantal scaling
        Matrix<double> Zold = Z;
        Matrix<double> Znew = mscale(madd(Zold, inverse(Zold)), 0.5);
        std::size_t numit = 0;
        double check = 1.0;
        while (check > epsilon && numit < opts.max_num_it) {
            ++numit;
            Zold = Znew;
            const double determ =
                opts.mode == "MSignStandard"
                    ? 0.5
                    : 1.0 / (1.0 + std::pow(std::abs(lu_det(Zold)), 1.0 / static_cast<double>(mf)));
            Znew = madd(mscale(Zold, determ), mscale(inverse(Zold), 1.0 - determ));
            check = max_col_sum(msub(Znew, Zold)) / max_col_sum(Zold);
        }
        // Step 7, T = orth(Znew - I)
        const Matrix<double> K = msub(Znew, eye<double>(mf));
        const SvdFactors sv = svd_full(K);
        // MATLAB orth: rank tolerance max(size) * eps(max(s))
        const double tol = static_cast<double>(mf) *
                           (std::nextafter(sv.s[0], std::numeric_limits<double>::infinity()) - sv.s[0]);
        std::size_t r = 0;
        for (double s : sv.s)
            if (s > tol) ++r;
        if (r != m)
            throw NumericError("MG1_IS: the sign iteration left an invariant subspace of dimension " +
                               std::to_string(r) + ", not m = " + std::to_string(m));
        T = sub(sv.U, 0, mf, 0, m);
    } else {
        const RealSchur sf = schur_decomposition(Z);
        std::vector<double> key(mf, 0.0);
        for (std::size_t k = 0; k < mf; ++k) key[k] = sf.T(k, k) < 0.0 ? 1.0 : 0.0;
        const RealSchur ord = schur_reorder(sf, key);
        T = sub(ord.Z, 0, mf, 0, m);
    }

    // Step 8
    const Matrix<double> T1 = sub(T, 0, m, 0, m), T2 = sub(T, m, m, 0, m);
    return matmul(madd(T1, T2), inverse(msub(T1, T2)));
}

// ---------------------------------------------------------------------------
// GIM1_R.m
// ---------------------------------------------------------------------------

/**
 * R of a GI/M/1-type Markov chain, through the G of its dual. Port of
 * `GIM1_R.m` for Dual 'A', 'R' and 'B' and every Algor: 'FI', 'CR', 'NI', 'RR', 'IS'.
 *
 * THE DUAL IS THE WHOLE IDEA. There is no cyclic reduction for R directly, so
 * the chain is transposed into an M/G/1-type one whose G carries the same
 * information: the Ramaswami dual `diag(theta)^-1 A_i' diag(theta)` for a
 * transient chain, and the Bright dual, which additionally rescales block i by
 * eta^(i-1) with eta the caudal characteristic, for a positive recurrent one.
 * 'A' picks between them by the drift, which is what makes it the fastest
 * default. R is then read back off G by the inverse similarity, times eta in
 * the Bright case.
 */
inline Matrix<double> gim1_r(const Blocks& Ain, const std::string& dual,
                             const std::string& algor) {
    const std::size_t m = Ain[0].rows();
    const std::size_t dega = Ain.size() - 1;
    Blocks A = Ain;

    // drift > 1: positive recurrent GI/M/1; drift < 1: transient.
    const Drift d = mg1_drift(A);
    std::vector<double> theta = d.theta;
    const bool ram = (dual == "R") || (dual == "A" && d.value <= 1.0);
    double eta = 1.0;

    if (ram) {
        for (std::size_t b = 0; b <= dega; ++b) {
            Matrix<double> Bb(m, m, 0.0);
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j) Bb(i, j) = A[b](j, i) * theta[j] / theta[i];
            A[b] = Bb;
        }
    } else if (dual == "B" || dual == "A") {
        eta = (d.value > 1.0) ? gim1_caudal(Ain) : mg1_decay(Ain);
        // theta of A0 + A1 eta + ... + Amax eta^max, shifted to be stochastic.
        Matrix<double> sumAeta = mscale(Ain[dega], std::pow(eta, static_cast<double>(dega)));
        for (std::size_t i = dega; i-- > 0;)
            sumAeta = madd(sumAeta, mscale(Ain[i], std::pow(eta, static_cast<double>(i))));
        Matrix<double> shifted = sumAeta;
        for (std::size_t i = 0; i < m; ++i) shifted(i, i) += (1.0 - eta);
        theta = stat(shifted);
        for (std::size_t b = 0; b <= dega; ++b) {
            Matrix<double> Bb(m, m, 0.0);
            const double s = std::pow(eta, static_cast<double>(b) - 1.0);
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j)
                    Bb(i, j) = s * A[b](j, i) * theta[j] / theta[i];
            A[b] = Bb;
        }
    } else {
        throw InputError("GIM1_R: Dual '" + dual + "' is not one of 'A', 'B', 'R'");
    }

    Matrix<double> G;
    if (algor == "FI") {
        G = mg1_fi(A);
    } else if (algor == "CR") {
        G = mg1_cr(A);
    } else if (algor == "NI") {
        G = mg1_ni(A);
    } else if (algor == "RR") {
        G = mg1_rr(A);
    } else if (algor == "IS") {
        G = mg1_is(A);
    } else {
        throw InputError("GIM1_R: Algorithm '" + algor + "' is not supported");
    }

    Matrix<double> R(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) R(i, j) = G(j, i) * theta[j] / theta[i];
    if (!ram)
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) R(i, j) *= eta;
    return R;
}

}  // namespace smc
}  // namespace line

#endif  // LINE_LIB_SMC_MG1_H
