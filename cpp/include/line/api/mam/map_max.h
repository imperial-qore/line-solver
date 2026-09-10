/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_MAX_H
#define LINE_API_MAM_MAP_MAX_H

/**
 * Maximum of two independent MAPs, and its marked k-stage generalization.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_max.m and
 * matlab/lib/m3a/m3a/mmap/mmap_max.m. These are the synchronization primitives
 * behind a fork-join: the joined interval is the maximum of the two branch
 * intervals, not their sum, so the phase space carries a RACE followed by the
 * residual of whichever branch is still running.
 *
 * map_max therefore has order na*nb + na + nb: the product block while both
 * branches are alive, then one absorbing-residual block per branch. Reading the
 * order as na*nb, as for a superposition, drops exactly the residual phases
 * that make the maximum different from the minimum.
 *
 * mmap_max keeps k rounds of the race, so its order is na*nb*(1+2k), with the
 * marks of both branches carried through unchanged.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace detail {

/** Copies src into the (br, bc) block of dst, blocks being rb x cb. */
template <class T>
void set_block(Matrix<T>& dst, std::size_t br, std::size_t bc, const Matrix<T>& src) {
    for (std::size_t i = 0; i < src.rows(); ++i)
        for (std::size_t j = 0; j < src.cols(); ++j) dst(br + i, bc + j) = src(i, j);
}

}  // namespace detail

/** MAP of the maximum of two independent MAPs. */
template <class T>
Map<T> map_max(const Map<T>& A, const Map<T>& B) {
    const std::size_t na = A.D0.rows(), nb = B.D0.rows();
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> a(na, zero), b(nb, zero);
    for (std::size_t i = 0; i < na; ++i)
        for (std::size_t j = 0; j < na; ++j) a[i] -= A.D0(i, j);
    for (std::size_t i = 0; i < nb; ++i)
        for (std::size_t j = 0; j < nb; ++j) b[i] -= B.D0(i, j);
    const std::size_t nab = na * nb, N = nab + nb + na;
    Matrix<T> M0(N, N, zero);
    const Matrix<T> S = krons(A.D0, B.D0);
    detail::set_block(M0, 0, 0, S);
    // kron(a, I_nb) sits to the right of the product block, then kron(I_na, b)
    for (std::size_t i = 0; i < na; ++i)
        for (std::size_t j = 0; j < nb; ++j) M0(i * nb + j, nab + j) = a[i];
    for (std::size_t i = 0; i < na; ++i)
        for (std::size_t j = 0; j < nb; ++j) M0(i * nb + j, nab + nb + i) = b[j];
    detail::set_block(M0, nab, nab, B.D0);
    detail::set_block(M0, nab + nb, nab + nb, A.D0);
    const std::vector<T> pa = map_pie(A), pb = map_pie(B);
    std::vector<T> pie(N, zero);
    for (std::size_t i = 0; i < na; ++i)
        for (std::size_t j = 0; j < nb; ++j) pie[i * nb + j] = pa[i] * pb[j];
    std::vector<T> d(N, zero);
    for (std::size_t j = 0; j < nb; ++j) d[nab + j] = b[j];
    for (std::size_t i = 0; i < na; ++i) d[nab + nb + i] = a[i];
    Matrix<T> M1(N, N, zero);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) M1(i, j) = d[i] * pie[j];
    return Map<T>{M0, M1};
}

/** MMAP of the maximum over k synchronization rounds of two independent MMAPs. */
template <class T>
Mmap<T> mmap_max(const Mmap<T>& a, const Mmap<T>& b, unsigned k) {
    if (k == 0) throw InputError("mmap_max: k must be positive");
    if (a.classes() != b.classes())
        throw InputError("mmap_max: the two MMAPs must carry the same number of classes");
    const std::size_t na = a.order(), nb = b.order(), n = na * nb;
    const std::size_t nblk = 1 + 2 * static_cast<std::size_t>(k);
    const std::size_t N = n * nblk;
    const T zero = num_traits<T>::from_int(0);
    const Matrix<T> Ia = eye<T>(na), Ib = eye<T>(nb);
    const Matrix<T> A0B0 = krons(a.D0, b.D0);
    const Matrix<T> A1IB = kron(a.D1, Ib);
    const Matrix<T> IAB1 = kron(Ia, b.D1);
    const Matrix<T> IAB0 = kron(Ia, b.D0);
    const Matrix<T> A0IB = kron(a.D0, Ib);
    Matrix<T> M0(N, N, zero);
    detail::set_block(M0, 0, 0, A0B0);
    detail::set_block(M0, 0, n, A1IB);
    detail::set_block(M0, 0, 2 * n, IAB1);
    for (std::size_t bi = 1; bi + 2 <= nblk - 1; ++bi)
        detail::set_block(M0, bi * n, bi * n, A0B0);
    detail::set_block(M0, (nblk - 2) * n, (nblk - 2) * n, IAB0);
    detail::set_block(M0, (nblk - 1) * n, (nblk - 1) * n, A0IB);
    for (unsigned i = 2; i <= k; ++i) {
        const std::size_t r = (1 + 2 * (i - 2)) * n, c = (3 + 2 * (i - 2)) * n;
        detail::set_block(M0, r, c, A1IB);
        detail::set_block(M0, r + n, c + n, IAB1);
    }
    Mmap<T> out;
    out.D0 = M0;
    out.D1 = Matrix<T>(N, N, zero);
    detail::set_block(out.D1, n, 0, IAB1);
    detail::set_block(out.D1, 2 * n, 0, A1IB);
    for (unsigned i = 2; i <= k; ++i) {
        const std::size_t r = (1 + 2 * (i - 1)) * n, c = (1 + 2 * (i - 2)) * n;
        detail::set_block(out.D1, r, c, IAB1);
        detail::set_block(out.D1, r + n, c + n, A1IB);
    }
    for (std::size_t cls = 0; cls < a.classes(); ++cls) {
        Matrix<T> Mc(N, N, zero);
        const Matrix<T> IABc = kron(Ia, b.Dc[cls]);
        const Matrix<T> AcIB = kron(a.Dc[cls], Ib);
        detail::set_block(Mc, n, 0, IABc);
        detail::set_block(Mc, 2 * n, 0, AcIB);
        for (unsigned i = 2; i <= k; ++i) {
            const std::size_t r = (1 + 2 * (i - 1)) * n, c = (1 + 2 * (i - 2)) * n;
            detail::set_block(Mc, r, c, IABc);
            detail::set_block(Mc, r + n, c + n, AcIB);
        }
        out.Dc.push_back(Mc);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif
