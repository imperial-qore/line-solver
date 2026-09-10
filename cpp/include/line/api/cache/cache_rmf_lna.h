/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_CACHE_RMF_LNA_H
#define LINE_API_CACHE_CACHE_RMF_LNA_H

/**
 * Stationary covariance of a RANDOM(m) cache occupancy, under the LNA.
 *
 * Port of `cache_rmf_lna` in python/line_solver/api/cache/rmf.py, itself a twin
 * of the MATLAB local `rmf_lna_covariance`. It solves the Lyapunov equation
 *
 *     F'(x) W + W F'(x)' + Q(x) = 0
 *
 * with exactly the `jacobian` and `noise_matrix` that `cache_miss_rmf.h`
 * already uses for the 1/N mean correction, so the mean and the covariance
 * linearise about the IDENTICAL drift. Reusing those two rather than restating
 * them is deliberate: the reference's `rmf_jacobian` is documented there as NOT
 * being the derivative of `rmf_drift` at every entry, and a second transcription
 * would silently pick the other one.
 *
 * THE SUBSPACE IS THE POINT, AND THE REASON THE SOLVE IS POSSIBLE AT ALL. The
 * Jacobian is singular twice over, because a RANDOM(m) cache conserves two
 * things: every item is in exactly one list (`sum_k x[i,k] = 1`) and every list
 * holds exactly its capacity (`sum_i x[i,k] = m[k]`). Every jump is a SWAP,
 * `(e_i - e_j) tensor (e_{k+1} - e_k)`, so the fluctuation lives on the tensor
 * product of the zero-sum ITEM space with the zero-sum LIST space -- the
 * double-centred subspace, of dimension `(n-1) * h`. Restricting to an
 * orthonormal basis of it is EXACT, not a regularization, and it is what makes
 * the covariance of a deterministic total come out as zero rather than as
 * whatever a pseudo-inverse would have produced.
 *
 * A fixed point that is not exponentially stable ON THAT SUBSPACE has no
 * stationary covariance, and is refused rather than answered.
 *
 * ARITHMETIC: double. The eigenvalue test and the Lyapunov solve are both
 * floating point.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/cache/cache_miss_rmf.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/sylvester.h"

namespace line {
namespace cache {

namespace lnadetail {

/**
 * Orthonormal basis of `{u in R^n : sum(u) = 0}`, as an n-by-(n-1) matrix.
 *
 * Built by Gram-Schmidt on the centring projector's columns, which is what the
 * reference's QR of `I - ones/n` produces. Any orthonormal basis of that
 * subspace serves: the restricted solve is basis-independent because the answer
 * is mapped back by `V W_r V'`.
 */
inline Matrix<double> centered_basis(std::size_t n) {
    if (n <= 1) return Matrix<double>(n, 0, 0.0);
    std::vector<std::vector<double>> cols;
    for (std::size_t c = 0; c < n && cols.size() + 1 < n; ++c) {
        std::vector<double> v(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            v[i] = (i == c ? 1.0 : 0.0) - 1.0 / static_cast<double>(n);
        for (std::size_t b = 0; b < cols.size(); ++b) {
            double d = 0.0;
            for (std::size_t i = 0; i < n; ++i) d += v[i] * cols[b][i];
            for (std::size_t i = 0; i < n; ++i) v[i] -= d * cols[b][i];
        }
        double nv = 0.0;
        for (std::size_t i = 0; i < n; ++i) nv += v[i] * v[i];
        nv = std::sqrt(nv);
        if (nv <= 1e-12) continue;
        for (std::size_t i = 0; i < n; ++i) v[i] /= nv;
        cols.push_back(v);
    }
    Matrix<double> U(n, cols.size(), 0.0);
    for (std::size_t c = 0; c < cols.size(); ++c)
        for (std::size_t i = 0; i < n; ++i) U(i, c) = cols[c][i];
    return U;
}

/** `kron(A, B)` in MATLAB's ordering. */
inline Matrix<double> kron_d(const Matrix<double>& A, const Matrix<double>& B) {
    Matrix<double> K(A.rows() * B.rows(), A.cols() * B.cols(), 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            for (std::size_t k = 0; k < B.rows(); ++k)
                for (std::size_t l = 0; l < B.cols(); ++l)
                    K(i * B.rows() + k, j * B.cols() + l) = A(i, j) * B(k, l);
    return K;
}

}  // namespace lnadetail

/**
 * @param x   the fluid fixed point, flattened item-major (`i + k*n`)
 * @param p   per-item popularities
 * @param m   per-list capacities
 * @param n   number of items
 * @param h   number of cache lists
 * @return    (dim x dim) stationary covariance, symmetric
 */
inline Matrix<double> cache_rmf_lna(const std::vector<double>& x, const std::vector<double>& p,
                                    const std::vector<double>& m, std::size_t n, std::size_t h,
                                    std::size_t dim) {
    if (x.size() != dim) throw InputError("cache_rmf_lna: the fixed point has the wrong length");

    const Matrix<double> Fp = rmf_detail::jacobian<double>(x, p, m, n, h);
    Matrix<double> Q = rmf_detail::noise_matrix<double>(x, p, m, n, h);
    // Symmetrize: the noise intensity is symmetric by construction and the
    // asymmetry that survives assembly is round-off.
    for (std::size_t i = 0; i < Q.rows(); ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const double v = 0.5 * (Q(i, j) + Q(j, i));
            Q(i, j) = v;
            Q(j, i) = v;
        }

    // The double-centred subspace: zero-sum over items, zero-sum over lists.
    const Matrix<double> Ui = lnadetail::centered_basis(n);
    const Matrix<double> Ul = lnadetail::centered_basis(h + 1);
    const Matrix<double> V = lnadetail::kron_d(Ul, Ui);  // item-major, i + k*n
    if (V.cols() == 0) return Matrix<double>(dim, dim, 0.0);
    if (V.rows() != dim)
        throw InputError(
            "cache_rmf_lna: the centred basis does not span the state dimension; n, h and dim "
            "disagree");

    const std::size_t d = V.cols();
    Matrix<double> Ar(d, d, 0.0), Qr(d, d, 0.0);
    for (std::size_t a = 0; a < d; ++a)
        for (std::size_t b = 0; b < d; ++b) {
            double sa = 0.0, sq = 0.0;
            for (std::size_t i = 0; i < dim; ++i)
                for (std::size_t j = 0; j < dim; ++j) {
                    sa += V(i, a) * Fp(i, j) * V(j, b);
                    sq += V(i, a) * Q(i, j) * V(j, b);
                }
            Ar(a, b) = sa;
            Qr(a, b) = sq;
        }
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const double v = 0.5 * (Qr(i, j) + Qr(j, i));
            Qr(i, j) = v;
            Qr(j, i) = v;
        }

    // Exponential stability ON THE SUBSPACE. Off it the Jacobian is singular by
    // conservation, so testing the full spectrum would refuse every model.
    const std::vector<std::complex<double>> ev = eig_values(Ar);
    double worst = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < ev.size(); ++i) worst = std::max(worst, ev[i].real());
    if (worst >= -std::sqrt(2.220446049250313e-16))
        throw InputError(
            "cache_rmf_lna: the cache fluid fixed point is not exponentially stable on the "
            "reachable subspace, so the occupancy process has no stationary covariance");

    // A W + W A' + Q = 0 with A = Ar; `lyap_solve(A, B, C)` solves A X + X B + C = 0.
    Matrix<double> ArT(d, d, 0.0);
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t j = 0; j < d; ++j) ArT(i, j) = Ar(j, i);
    Matrix<double> Wr = lyap_solve(Ar, ArT, Qr);
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const double v = 0.5 * (Wr(i, j) + Wr(j, i));
            Wr(i, j) = v;
            Wr(j, i) = v;
        }

    Matrix<double> W(dim, dim, 0.0);
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim; ++j) {
            double s = 0.0;
            for (std::size_t a = 0; a < d; ++a)
                for (std::size_t b = 0; b < d; ++b) s += V(i, a) * Wr(a, b) * V(j, b);
            W(i, j) = s;
        }
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const double v = 0.5 * (W(i, j) + W(j, i));
            W(i, j) = v;
            W(j, i) = v;
        }
    return W;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_CACHE_RMF_LNA_H
