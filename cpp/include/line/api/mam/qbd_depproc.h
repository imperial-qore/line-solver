/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_QBD_DEPPROC_H
#define LINE_API_MAM_QBD_DEPPROC_H

/**
 * Departure process of a MAP/MAP/1 queue: the ETAQA-truncated MAP descriptor
 * under FCFS and under PS, and the joint moments of consecutive
 * inter-departure times.
 *
 * Templated port of matlab/src/api/mam/qbd_depproc_etaqa.m,
 * qbd_depproc_etaqa_ps.m and qbd_depproc_jointmom.m.
 *
 * The queue is the usual MAP/MAP/1 QBD, level = number in system, phase =
 * (arrival phase, service phase) with the arrival phase major:
 *
 *     F  = D1^a (x) I_ns       L  = D0^a (+) D0^s
 *     B  = I_na (x) D1^s       L0 = D0^a (x) I_ns
 *
 * ETAQA keeps levels 0..n-1 explicitly and lumps every level from n upwards
 * into a single aggregate block, using G to describe how the aggregate returns
 * to level n-1: Bbar = B + F G and Bhat = F G, with Lhat = F + L. The
 * descriptor has (n+1) blocks of order na*ns; a transition that carries a
 * departure goes into D1 and everything else into D0.
 *
 * ARITHMETIC. All three entry points call qbd_fundmat for R and G, so they are
 * gated on num_traits<T>::has_transcendental; see qbd_r.h for why cyclic
 * reduction cannot be exact. Everything downstream of R and G is a finite
 * sequence of matrix products and inverses.
 *
 * REFERENCE DEFECTS (reproduced here verbatim, see the tests).
 *
 * 1. qbd_depproc_etaqa returns a pair (D0, D1) that is NOT a conservative
 *    generator: (D0 + D1) e is nonzero in the last two block rows, so the
 *    result is not a MAP and map_lambda / map_prob applied to it are
 *    meaningless. Two independent causes, both an off-by-one in the block
 *    index arithmetic:
 *      a) the line
 *           D0(((n-1)*lvlsz+1):n*lvlsz, ((n-1)*lvlsz+1):n*lvlsz) = Lhat
 *         addresses the FULL padded matrix, whose block rows are 0..n after
 *         the two zero paddings, so it writes Lhat = F + L into block
 *         (n-1, n-1), the last EXPLICIT level, and not into block (n, n), the
 *         aggregate. Block (n-1, n-1) already has the up-block F sitting at
 *         (n-1, n), so that row acquires F twice and its row sum becomes F e
 *         instead of 0, while the aggregate keeps a bare L.
 *      b) the aggregate row carries both Bbar = B + F G at (n, n-1) and
 *         Bhat = F G at (n, n), so F G is counted twice there; the row sum is
 *         F G e = F e rather than 0.
 *    Measured on the M/M/1 instance of the test (lambda = 0.6, mu = 1, n = 4):
 *    ||(D0 + D1) e||_inf = 0.6 = lambda, concentrated in exactly those two
 *    rows and zero everywhere else.
 *
 * 2. qbd_depproc_etaqa_ps additionally puts the FULL down-block Bbar at
 *    (n, n-1) into BOTH D0 and D1, instead of splitting it (1 - 1/n) / (1/n)
 *    the way it splits B at every explicit level and Bhat at the aggregate.
 *    The aggregate therefore fires that transition at twice its rate, once as
 *    a departure and once silently, and (D0 + D1) e picks up a further Bbar e.
 *    It also divides by j at level j starting from j = 1, so the level-1 row
 *    gets B * (1 - 1/1) = 0 for the non-departure part, which is right, but
 *    the loop stops at n-1 and never treats the boundary the same way.
 *
 * The port does NOT repair either: these functions are the reference for the
 * JAR and Python ports and a silent divergence would be worse than a
 * reproduced defect. qbd_depproc_etaqa_residual below measures ||(D0+D1)e||_inf
 * so a caller can see it, and the tests pin the measured value.
 */

#include <cstddef>
#include <utility>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/qbd_mapmap1.h"
#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace depproc_detail {

/** Copy src into dest at block position (br, bc) of blocks of size sz. */
template <class T>
void put_block(Matrix<T>& dest, std::size_t br, std::size_t bc, const Matrix<T>& src) {
    const std::size_t r0 = br * src.rows(), c0 = bc * src.cols();
    if (r0 + src.rows() > dest.rows() || c0 + src.cols() > dest.cols())
        throw InputError("qbd_depproc: block write out of range");
    for (std::size_t i = 0; i < src.rows(); ++i)
        for (std::size_t j = 0; j < src.cols(); ++j) dest(r0 + i, c0 + j) = src(i, j);
}

/** Add src into dest at block position (br, bc). */
template <class T>
void add_block(Matrix<T>& dest, std::size_t br, std::size_t bc, const Matrix<T>& src) {
    const std::size_t r0 = br * src.rows(), c0 = bc * src.cols();
    if (r0 + src.rows() > dest.rows() || c0 + src.cols() > dest.cols())
        throw InputError("qbd_depproc: block write out of range");
    for (std::size_t i = 0; i < src.rows(); ++i)
        for (std::size_t j = 0; j < src.cols(); ++j) dest(r0 + i, c0 + j) += src(i, j);
}

/** The ETAQA pieces shared by the FCFS and the PS construction. */
template <class T>
struct EtaqaPieces {
    Matrix<T> F, L, B, L0;
    Matrix<T> R, G;
    Matrix<T> Lhat;  ///< F + L
    Matrix<T> Bbar;  ///< B + F G
    Matrix<T> Bhat;  ///< F G
    std::size_t lvlsz;
};

template <class T>
EtaqaPieces<T> etaqa_pieces(const Map<T>& arrival, const Map<T>& service) {
    using namespace qbd_detail;
    const QbdMapMap1Blocks<T> blk = qbd_mapmap1_blocks(arrival, service);
    EtaqaPieces<T> p;
    p.F = blk.F;
    p.L = blk.L;
    p.B = blk.B;
    p.L0 = blk.Lbar;
    p.lvlsz = blk.L.rows();
    const QbdFundMat<T> fm = qbd_fundmat(p.B, p.L, p.F);
    p.R = fm.R;
    // G recomputation rationale: see _kb/03-api-layer.md (cpp port notes: mam)
    p.G = matmul(inverse(msub(mscale(p.L, T(num_traits<T>::from_int(-1))), matmul(p.R, p.B))), p.B);
    p.Lhat = madd(p.F, p.L);
    p.Bbar = madd(p.B, matmul(p.F, p.G));
    p.Bhat = matmul(p.F, p.G);
    return p;
}

}  // namespace depproc_detail

/**
 * MAP descriptor of the departure process of a MAP/MAP/1-FCFS queue,
 * ETAQA-truncated at level n (qbd_depproc_etaqa.m).
 *
 * @param n number of explicitly represented levels; the descriptor has n+1
 *          blocks of order na*ns, the last one being the aggregate
 * @param arrival the arrival MAP
 * @param service the service MAP
 */
template <class T>
Map<T> qbd_depproc_etaqa(const Map<T>& arrival, const Map<T>& service, std::size_t n) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_depproc_etaqa requires transcendental arithmetic");
    if (n < 2) throw InputError("qbd_depproc_etaqa: the truncation level must be at least 2");
    using namespace depproc_detail;
    const EtaqaPieces<T> p = etaqa_pieces(arrival, service);
    const std::size_t m = p.lvlsz;
    const std::size_t dim = (n + 1) * m;
    const T zero = num_traits<T>::from_int(0);

    Matrix<T> D0(dim, dim, zero);
    // Diagonal L on blocks 1..n and superdiagonal F on blocks (1,2)..(n-1,n),
    // which is the padded image of kron(I_n, L) + kron(superdiag, F).
    for (std::size_t r = 1; r <= n; ++r) put_block(D0, r, r, p.L);
    for (std::size_t r = 1; r + 1 <= n; ++r) put_block(D0, r, r + 1, p.F);
    // Boundary row: [L0, F] over the first two block columns.
    put_block(D0, 0, 0, p.L0);
    put_block(D0, 0, 1, p.F);
    // Reference defect 1a: Lhat lands on block (n-1, n-1), not on the aggregate.
    put_block(D0, n - 1, n - 1, p.Lhat);

    Matrix<T> D1(dim, dim, zero);
    put_block(D1, n, n - 1, p.Bbar);
    put_block(D1, n, n, p.Bhat);
    // The leading n x n block region is overwritten with the plain subdiagonal
    // B, which leaves the aggregate row and column untouched.
    for (std::size_t r = 1; r + 1 <= n; ++r) put_block(D1, r, r - 1, p.B);

    Map<T> out;
    out.D0 = D0;
    out.D1 = D1;
    return out;
}

/**
 * MAP descriptor of the departure process of a MAP/MAP/1-PS queue,
 * ETAQA-truncated at level n (qbd_depproc_etaqa_ps.m).
 *
 * With j jobs sharing the server a completion is a departure of the tagged
 * job with probability 1/j, so the down-block at level j splits into
 * B * (1/j) into D1 and B * (1 - 1/j) into D0.
 */
template <class T>
Map<T> qbd_depproc_etaqa_ps(const Map<T>& arrival, const Map<T>& service, std::size_t n) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_depproc_etaqa_ps requires transcendental arithmetic");
    if (n < 2) throw InputError("qbd_depproc_etaqa_ps: the truncation level must be at least 2");
    using namespace depproc_detail;
    using namespace qbd_detail;
    const EtaqaPieces<T> p = etaqa_pieces(arrival, service);
    const std::size_t m = p.lvlsz;
    const std::size_t dim = (n + 1) * m;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T nT = num_traits<T>::from_int(static_cast<long>(n));

    Matrix<T> D0(dim, dim, zero);
    for (std::size_t r = 1; r <= n; ++r) put_block(D0, r, r, p.L);
    for (std::size_t r = 1; r + 1 <= n; ++r) put_block(D0, r, r + 1, p.F);
    put_block(D0, 0, 0, p.L0);
    put_block(D0, 0, 1, p.F);
    put_block(D0, n - 1, n - 1, p.Lhat);
    // Reference defect 2: the full Bbar goes into D0 as well as into D1.
    add_block(D0, n, n - 1, p.Bbar);
    add_block(D0, n, n, mscale(p.Bhat, T((nT - one) / nT)));

    Matrix<T> D1(dim, dim, zero);
    put_block(D1, n, n - 1, p.Bbar);
    put_block(D1, n, n, mscale(p.Bhat, T(one / nT)));

    for (std::size_t j = 1; j + 1 <= n; ++j) {
        const T jT = num_traits<T>::from_int(static_cast<long>(j));
        put_block(D0, j, j - 1, mscale(p.B, T(one - one / jT)));
        put_block(D1, j, j - 1, mscale(p.B, T(one / jT)));
    }

    Map<T> out;
    out.D0 = D0;
    out.D1 = D1;
    return out;
}

/**
 * ||(D0 + D1) e||_inf, zero for a genuine MAP. Exposed so a caller can see
 * reference defects 1 and 2 rather than discovering them downstream.
 */
template <class T>
T qbd_depproc_residual(const Map<T>& m) {
    const std::size_t n = m.D0.rows();
    T worst = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < m.D0.cols(); ++j) s += m.D0(i, j) + m.D1(i, j);
        const T a = num_abs(T(s));
        if (a > worst) worst = a;
    }
    return worst;
}

/**
 * Joint moments E[X_0^i X_1^j] of consecutive inter-departure times of a
 * MAP/MAP/1-FCFS queue (qbd_depproc_jointmom.m).
 *
 * The initial vector is built from the level-0 vector of the QBD and the
 * first three terms of the geometric tail,
 *
 *     z = [v0 R F, v0 R^2 F, v0 R^3 (I-R)^-1 F] / lambda_s,
 *
 * normalized to a probability vector, and the moments follow from the
 * three-level block operators
 *
 *     M0 = [L0 F 0; 0 L F; 0 0 L+F],   M1 = [0 0 0; B 0 0; 0 B 0]
 *
 * as JM = z i! (-M0)^{-i-1} M1 j! (-M0)^{-j} e.
 *
 * @param iset one (i, j) pair per requested moment
 * @param arrival the arrival MAP
 * @param service the service MAP
 */
template <class T>
std::vector<T> qbd_depproc_jointmom(const Map<T>& arrival, const Map<T>& service,
                                    const std::vector<std::pair<unsigned, unsigned>>& iset) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_depproc_jointmom requires transcendental arithmetic");
    using namespace depproc_detail;
    using namespace qbd_detail;
    const EtaqaPieces<T> p = etaqa_pieces(arrival, service);
    const std::size_t m = p.lvlsz;
    const T zero = num_traits<T>::from_int(0);

    const Matrix<T> pi = qbd_pi(p.B, p.L0, p.R);
    std::vector<T> v0(m);
    for (std::size_t j = 0; j < m; ++j) v0[j] = pi(0, j);

    const T lambdaS = map_lambda(service);
    if (lambdaS == zero) throw NumericError("qbd_depproc_jointmom: zero service rate");
    const T inv = num_traits<T>::from_int(1) / lambdaS;

    const std::vector<T> v0R = vecmul(v0, p.R);
    const std::vector<T> v0R2 = vecmul(v0R, p.R);
    const std::vector<T> v0R3 = vecmul(v0R2, p.R);
    const Matrix<T> ImRinv = inverse(msub(eye<T>(m), p.R));

    // Departure epochs are the B transitions, so the embedded vector weighs the
    // level probabilities by B and not by the arrival matrix F.
    const std::vector<T> v0D = vecmul(v0R, p.B);
    const std::vector<T> v1D = vecmul(v0R2, p.B);
    const std::vector<T> v2Dp = vecmul(vecmul(v0R3, ImRinv), p.B);

    std::vector<T> z(3 * m);
    for (std::size_t j = 0; j < m; ++j) {
        z[j] = inv * v0D[j];
        z[m + j] = inv * v1D[j];
        z[2 * m + j] = inv * v2Dp[j];
    }
    T zs = zero;
    for (const T& v : z) zs += v;
    if (zs == zero) throw NumericError("qbd_depproc_jointmom: degenerate initial vector");
    for (T& v : z) v /= zs;

    Matrix<T> M0(3 * m, 3 * m, zero), M1(3 * m, 3 * m, zero);
    put_block(M0, 0, 0, p.L0);
    put_block(M0, 0, 1, p.F);
    put_block(M0, 1, 1, p.L);
    put_block(M0, 1, 2, p.F);
    put_block(M0, 2, 2, p.Lhat);
    put_block(M1, 1, 0, p.B);
    put_block(M1, 2, 1, p.B);

    const Matrix<T> negM0inv = inverse(mscale(M0, T(num_traits<T>::from_int(-1))));
    const std::vector<T> e = ones<T>(3 * m);

    std::vector<T> out;
    out.reserve(iset.size());
    for (std::size_t k = 0; k < iset.size(); ++k) {
        const unsigned i = iset[k].first, j = iset[k].second;
        std::vector<T> row = vecmul(z, matpow(negM0inv, i + 1));
        for (T& v : row) v *= num_factorial<T>(i);
        row = vecmul(row, M1);
        row = vecmul(row, matpow(negM0inv, j));
        for (T& v : row) v *= num_factorial<T>(j);
        T s = zero;
        for (std::size_t q = 0; q < row.size(); ++q) s += row[q] * e[q];
        out.push_back(s);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_QBD_DEPPROC_H
