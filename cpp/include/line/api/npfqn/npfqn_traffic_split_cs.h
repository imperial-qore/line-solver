/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_TRAFFIC_SPLIT_CS_H
#define LINE_API_NPFQN_TRAFFIC_SPLIT_CS_H

/**
 * Splitting of a marked MAP departure flow with class switching.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_traffic_split_cs.m,
 * cross-checked against
 * jar/src/main/java/jline/api/npfqn/Npfqn_traffic_split_cs.java (identical).
 *
 * An MMAP is the cell {D0, D1, D1^(1), ..., D1^(R)}: the hidden generator D0,
 * the aggregate arrival matrix D1 = sum_r D1^(r), and one marking matrix per
 * class. P(r, (j-1)R + s) is the probability that a class-r departure flows to
 * destination j in class s. For each destination j the port builds
 *
 *   D1^(s)_j = sum_r D1^(r) P(r, (j-1)R + s)
 *   D1_j     = sum_s D1^(s)_j
 *   D0_j     = D0 + D1 - D1_j
 *
 * i.e. everything not routed to j is folded back into the hidden part, which
 * is exactly the MATLAB accumulation written out.
 *
 * Arithmetic. Only additions and multiplications by routing probabilities,
 * plus the max(.,0) clipping of mmap_normalize, so the algorithm is exact at
 * T = Rational and needs no transcendental function.
 *
 * mmap_normalize (matlab/lib/m3a/m3a/mmap/mmap_normalize.m, JAR
 * jline.api.mam.Mmap_normalize) is inlined here as a detail helper because the
 * mam/M3A domain is not part of this port. It is a verbatim transcription: the
 * NaN test mirrors MATLAB's "if isnan(X)" on a matrix, which is true only when
 * every entry is NaN, and is vacuously false in an exact field.
 */

#include <cstddef>
#include <vector>

#include "line/api/npfqn/npfqn_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

/** An MMAP as the MATLAB cell {D0, D1, D1^(1), ..., D1^(R)}. */
template <class T>
using Mmap = std::vector<Matrix<T>>;

namespace detail {

/**
 * Port of mmap_normalize: clip the off-diagonal of D0 and every marking matrix
 * at zero, rebuild D1 as the sum of the markings, and reset the diagonal of D0
 * so that D0 + D1 is a generator.
 */
template <class T>
void mmap_normalize(Mmap<T>& M) {
    if (M.empty()) return;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = M[0].rows();
    const std::size_t C = M.size() - 2;

    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j)
            if (i != j && M[0](i, j) < zero) M[0](i, j) = zero;

    M[1] = Matrix<T>(M[0].rows(), M[0].cols(), zero);
    for (std::size_t c = 0; c < C; ++c) {
        Matrix<T>& Dc = M[2 + c];
        bool allNaN = Dc.rows() > 0 && Dc.cols() > 0;
        for (std::size_t i = 0; i < Dc.rows(); ++i)
            for (std::size_t j = 0; j < Dc.cols(); ++j) {
                if (!num_isnan(Dc(i, j))) allNaN = false;
                if (Dc(i, j) < zero) Dc(i, j) = zero;
            }
        if (allNaN) Dc = Matrix<T>(Dc.rows(), Dc.cols(), zero);
        for (std::size_t i = 0; i < Dc.rows(); ++i)
            for (std::size_t j = 0; j < Dc.cols(); ++j) M[1](i, j) += Dc(i, j);
    }

    for (std::size_t k = 0; k < K; ++k) {
        M[0](k, k) = zero;
        T s = zero;
        for (std::size_t j = 0; j < M[0].cols(); ++j) s += M[0](k, j);
        for (std::size_t j = 0; j < M[1].cols(); ++j) s += M[1](k, j);
        M[0](k, k) = -s;
    }
}

}  // namespace detail

/**
 * @param MMAP the departure MMAP {D0, D1, D1^(1), ..., D1^(R)}
 * @param P    (R x J) with J = M*R, P(r, (j-1)R + s) in MATLAB 1-based terms
 * @return one MMAP per destination station, in destination order
 */
template <class T>
std::vector<Mmap<T>> npfqn_traffic_split_cs(const Mmap<T>& MMAP, const Matrix<T>& P) {
    if (MMAP.size() < 3) throw InputError("npfqn_traffic_split_cs: MMAP has no marking matrices");
    const std::size_t R = P.rows();
    const std::size_t J = P.cols();
    if (R == 0 || J % R != 0)
        throw InputError("npfqn_traffic_split_cs: the class-switching matrix is not R x (M R)");
    if (MMAP.size() - 2 != R)
        throw InputError("npfqn_traffic_split_cs: MMAP and P disagree on the class count");
    const std::size_t M = J / R;

    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = MMAP[0].rows();

    std::vector<Mmap<T>> SMMAP(M);
    for (std::size_t jst = 0; jst < M; ++jst) {
        Mmap<T>& S = SMMAP[jst];
        S.assign(2 + R, Matrix<T>(n, MMAP[0].cols(), zero));
        // D0 starts as D0 + D1, with every arrival folded into the hidden part
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = 0; b < MMAP[0].cols(); ++b)
                S[0](a, b) = MMAP[0](a, b) + MMAP[1](a, b);
        for (std::size_t s = 0; s < R; ++s) {
            for (std::size_t r = 0; r < R; ++r) {
                const T& p = P(r, jst * R + s);
                for (std::size_t a = 0; a < n; ++a)
                    for (std::size_t b = 0; b < MMAP[2 + r].cols(); ++b) {
                        const T v = MMAP[2 + r](a, b) * p;
                        S[2 + s](a, b) += v;
                        S[1](a, b) += v;
                        S[0](a, b) -= v;
                    }
            }
        }
        detail::mmap_normalize(S);
    }
    return SMMAP;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_TRAFFIC_SPLIT_CS_H
