/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_QBD_BMAPBMAP1_H
#define LINE_API_MAM_QBD_BMAPBMAP1_H

/**
 * Level blocks of a BMAP/MAP/1 queue: batch Markovian arrivals against a
 * single-departure MAP service process.
 *
 * Templated port of matlab/src/api/mam/qbd_bmapbmap1.m. The level is the
 * number in system and the phase is the pair (arrival phase, service phase)
 * with the arrival phase major, the same layout as qbd_mapmap1. A batch of
 * size b raises the level by b, so the process is an M/G/1-type chain rather
 * than a QBD, with one up-block per batch size:
 *
 *     A1[b] = (D1^a p_b) (x) I_ns     an arrival batch of size b, level up by b
 *     A0    = D0^a (+) D0^s           no event
 *     A_1   = I_na (x) D1^s           a service completion, level down by one
 *     A0bar = D0^a (x) I_ns           level zero, no server busy
 *     B0    = D0^a (+) I_ns           level-zero local block
 *     B1[b] = A1[b]                   level-zero up-blocks
 *
 * REFERENCE DEFECT. The MATLAB function declares NO output argument and
 * assembles the blocks into local variables that are discarded on return; the
 * commented-out tail of the file shows the intended Q matrix. Calling it
 * therefore computes the blocks and throws them away, and in a context that
 * expects a value MATLAB raises "Output argument not assigned". The port
 * returns the blocks in a struct, which is the only content the function
 * actually has; nothing else can be inferred from it, in particular no
 * solution of the chain, because none is written.
 *
 * Two further observations on the reference, reproduced rather than repaired:
 *   - despite the name, the SERVICE process is an ordinary MAP: A_1 uses
 *     D1^s with no batch-size distribution, so only the arrivals are batched.
 *   - B0 = krons(D0^a, I_ns) is the Kronecker SUM of D0^a with the identity,
 *     i.e. D0^a (x) I + I (x) I, which adds an unconditional unit rate on the
 *     diagonal of every phase; the analogous boundary block in qbd_mapmap1 is
 *     the Kronecker PRODUCT kron(D0^a, I_ns), which is what A0bar holds here.
 *     Both are returned under their own names so the discrepancy is visible.
 *
 * ARITHMETIC. Kronecker products and sums only, so exact at Rational.
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

/** The level blocks assembled by qbd_bmapbmap1. */
template <class T>
struct QbdBmapBmap1Blocks {
    std::vector<Matrix<T>> A1;  ///< A1[b-1] is the up-block for a batch of size b
    Matrix<T> A0;               ///< local block
    Matrix<T> Am1;              ///< down-block (A_1 in the reference)
    Matrix<T> A0bar;            ///< kron(D0^a, I_ns)
    Matrix<T> B0;               ///< boundary local block, krons(D0^a, I_ns)
    std::vector<Matrix<T>> B1;  ///< boundary up-blocks, equal to A1
};

/**
 * Level blocks of a BMAP/MAP/1 queue.
 *
 * @param arrival arrival MAP, whose D1 is split by the batch-size law
 * @param pbatch  batch-size probabilities, pbatch[b-1] = P(batch = b)
 * @param service service MAP
 */
template <class T>
QbdBmapBmap1Blocks<T> qbd_bmapbmap1(const Map<T>& arrival, const std::vector<T>& pbatch,
                                    const Map<T>& service) {
    const std::size_t na = arrival.order();
    const std::size_t ns = service.order();
    if (na == 0 || ns == 0) throw InputError("qbd_bmapbmap1: empty MAP");
    if (pbatch.empty()) throw InputError("qbd_bmapbmap1: empty batch-size distribution");

    QbdBmapBmap1Blocks<T> out;
    const Matrix<T> Ins = eye<T>(ns);
    for (std::size_t b = 0; b < pbatch.size(); ++b) {
        Matrix<T> scaled(na, na);
        for (std::size_t i = 0; i < na; ++i)
            for (std::size_t j = 0; j < na; ++j) scaled(i, j) = arrival.D1(i, j) * pbatch[b];
        out.A1.push_back(kron(scaled, Ins));
    }
    out.A0 = krons(arrival.D0, service.D0);
    out.Am1 = kron(eye<T>(na), service.D1);
    out.A0bar = kron(arrival.D0, Ins);
    out.B0 = krons(arrival.D0, Ins);
    out.B1 = out.A1;
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_QBD_BMAPBMAP1_H
