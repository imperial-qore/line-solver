/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LIB_SMC_MG1_H
#define LINE_LIB_SMC_MG1_H

/**
 * The M/G/1-type and GI/M/1-type fundamental-matrix solvers of MAMSolver /
 * SMCSolver, ported from `matlab/lib/thirdparty/MG1files`:
 * `stat.m`, `MG1_EG.m`, `MG1_Decay.m`, `GIM1_Caudal.m`, `MG1_Shifts.m`,
 * `MG1_CR.m`, `MG1_FI.m` and `GIM1_R.m`.
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
 * TWO REFERENCE BRANCHES ARE REFUSED BY NAME rather than transcribed, because
 * they cannot run in MATLAB either. `MG1_Shifts` writes
 *     rowhatA(1,maxd*i:end) = uT*A(:,maxd*i:end)
 * in three places (the drift < 1 'tau' branch and the drift > 1 'one' branch),
 * where `i` is not a loop variable at that point: it is either MATLAB's
 * imaginary unit, which makes the index complex and errors, or the leftover
 * value 1 from the `beta` loop above, which addresses column `maxd` of a matrix
 * whose blocks start every m columns. Either way the line does not compute
 * "the last block", which is what the surrounding code needs. ETAQA reaches
 * neither branch: it shifts a positive recurrent chain (drift < 1) with the
 * default ShiftType 'one', which is the branch that is correct and is ported.
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
 * Decay rate of a recurrent M/G/1-type chain: the unique z > 1 with
 * PF(A(z)) = z. Port of `MG1_Decay.m`; the eigenvector output the reference
 * offers is not returned, because the only caller that wants it is the
 * 'tau' shift, which is refused (see the header note).
 */
inline double mg1_decay(const Blocks& A) {
    double eta = 1.0, new_eta = 0.0;
    while (new_eta - eta < 0.0) {
        eta += 1.0;
        new_eta = max_eig(poly_at(A, eta)).real();
    }
    double eta_min = eta - 1.0, eta_max = eta;
    eta = eta_min + 0.5;
    while (eta_max - eta_min > 1e-15) {
        new_eta = max_eig(poly_at(A, eta)).real();
        if (new_eta < eta) {
            eta_min = eta;
        } else {
            eta_max = eta;
        }
        eta = (eta_min + eta_max) / 2.0;
    }
    return eta;
}

/**
 * Caudal characteristic of a GI/M/1-type chain: the spectral radius of R, the
 * unique z in (0,1) with PF(A(z)) = z. Port of `GIM1_Caudal.m`.
 */
inline double gim1_caudal(const Blocks& A) {
    double eta_min = 0.0, eta_max = 1.0, eta = 0.5;
    while (eta_max - eta_min > 1e-15) {
        const double new_eta = max_eig(poly_at(A, eta)).real();
        if (new_eta > eta) {
            eta_min = eta;
        } else {
            eta_max = eta;
        }
        eta = (eta_min + eta_max) / 2.0;
    }
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
 * Shift technique for the M/G/1-type sequence. Port of `MG1_Shifts.m`,
 * ShiftType 'one', which is the default and the only type ETAQA uses.
 *
 * For a positive recurrent chain (drift < 1) the eigenvalue 1 of A(z) is
 * shifted to zero by subtracting `(A0+...+Ai)e u^T` from block i with
 * `u^T = e^T/m`; cyclic reduction then converges linearly in the SECOND
 * largest root instead of stalling on the unit one, which is the entire point
 * of running CR on the shifted sequence and undoing the shift on G afterwards.
 *
 * 'tau' and 'dbl', and 'one' at drift > 1, are refused: see the header.
 */
inline ShiftResult mg1_shifts(const Blocks& Ain, const std::string& shift_type) {
    if (shift_type != "one")
        throw UnsupportedError(
            "MG1_Shifts: ShiftType '" + shift_type +
            "' is not ported. Its reference branch writes rowhatA(1,maxd*i:end) with `i` "
            "undefined at that point, so it addresses column maxd of a sequence whose blocks "
            "start every m columns and does not compute the last block it needs; the port "
            "refuses rather than transcribe a line that cannot run. ETAQA uses ShiftType 'one'");

    Blocks A = Ain;
    const std::size_t m = A[0].rows();
    const Drift d = mg1_drift(A);

    if (!(d.value < 1.0))
        throw UnsupportedError(
            "MG1_Shifts: the drift > 1 branch of ShiftType 'one' is not ported, for the same "
            "reason as ShiftType 'tau': it shifts one to infinity through the same defective "
            "rowhatA(1,maxd*i:end) line. A transient M/G/1-type chain is outside what ETAQA "
            "solves here");

    // Shift one to zero: A1 <- A1 - I, then hatA_i = A_i - (A0+...+Ai) e u^T.
    for (std::size_t i = 0; i < m; ++i) A[1](i, i) -= 1.0;
    std::vector<double> col(m, 0.0);
    Blocks hatA(A.size(), Matrix<double>(m, m, 0.0));
    for (std::size_t b = 0; b < A.size(); ++b) {
        const std::vector<double> rs = rowsums(A[b]);
        for (std::size_t i = 0; i < m; ++i) col[i] += rs[i];
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j)
                hatA[b](i, j) = A[b](i, j) - col[i] / static_cast<double>(m);
    }
    for (std::size_t i = 0; i < m; ++i) hatA[1](i, i) += 1.0;

    ShiftResult out;
    out.hatA = hatA;
    out.drift = d.value;
    out.v.assign(m, 0.0);
    return out;
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
    double drift = 0.0;
    if (opts.mode == "ShiftPWCR") {
        const ShiftResult sh = mg1_shifts(A, opts.shift_type);
        A = sh.hatA;
        drift = sh.drift;
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

    // Undo the shift: shifting one to zero removed the rank-one term e u^T
    // from G, so it goes back on.
    if (opts.mode == "ShiftPWCR" && drift < 1.0)
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) += 1.0 / static_cast<double>(m);
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
    double drift = 0.0;
    if (shifted) {
        const ShiftResult sh = mg1_shifts(A, opts.shift_type);
        A = sh.hatA;
        drift = sh.drift;
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

    if (shifted && drift < 1.0)
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) += 1.0 / static_cast<double>(m);
    return G;
}

// ---------------------------------------------------------------------------
// GIM1_R.m
// ---------------------------------------------------------------------------

/**
 * R of a GI/M/1-type Markov chain, through the G of its dual. Port of
 * `GIM1_R.m` for Dual 'A', 'R' and 'B' and Algor 'FI' and 'CR'.
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
    } else {
        throw UnsupportedError(
            "GIM1_R: Algor '" + algor +
            "' is not ported; MG1_NI, MG1_RR and MG1_IS have no C++ counterpart. "
            "'FI' (the one GIM1_R_ETAQA asks for) and 'CR' are available");
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
