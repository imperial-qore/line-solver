/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PAS_PLACEMENT_H
#define LINE_API_PFQN_PAS_PLACEMENT_H

/**
 * Placement order of a pass-and-swap (P&S) order-independent network.
 *
 * Templated port of matlab/src/api/pfqn/pas_placement.m. An ordering
 * c = (c_1, ..., c_l) is feasible iff class a never precedes class b whenever
 * H(b, a) is nonzero (Comte and Dorsman, 2021, arXiv:2009.12299), so H is read
 * as "row must be placed before column" and its transitive closure P is the
 * full precedence relation. The closure is taken by repeated Boolean squaring
 * against H until it stops growing, which is the reference's own iteration and
 * terminates in at most R steps.
 *
 * placeable(x) returns the classes that may be placed next given the vector x
 * of remaining per-class counts: class j is placeable iff x(j) > 0 and no
 * still-present class must precede it, sum_i x(i) P(i, j) == 0.
 *
 * ARITHMETIC. The closure is Boolean and the placeable test is a sum of
 * counts, so nothing here rounds and the header is instantiated at Rational as
 * well as double and Real. P is returned as a 0/1 matrix of T, matching
 * MATLAB's double(P), so it can be multiplied straight into the caller's
 * arithmetic.
 *
 * EMPTY H. MATLAB returns P = [] and a placeable that admits every present
 * class. Reproduced: an empty H yields an empty P, and placeable then ignores
 * the precedence test.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Result of pas_placement. */
template <class T>
struct PasPlacement {
    Matrix<T> P;  ///< precedence closure, P(i, j) = 1 iff i must precede j

    /**
     * Classes that may be placed next, given the remaining per-class counts.
     * @return 0-based class indices, in increasing order (MATLAB's find order)
     */
    std::vector<std::size_t> placeable(const std::vector<T>& x) const {
        const T zero = num_traits<T>::from_int(0);
        std::vector<std::size_t> idx;
        if (P.rows() == 0) {
            for (std::size_t j = 0; j < x.size(); ++j)
                if (x[j] > zero) idx.push_back(j);
            return idx;
        }
        if (x.size() != P.rows()) throw InputError("pas_placement: x has the wrong class count");
        for (std::size_t j = 0; j < P.cols(); ++j) {
            if (!(x[j] > zero)) continue;
            T s = zero;
            for (std::size_t i = 0; i < P.rows(); ++i) s += x[i] * P(i, j);
            if (s == zero) idx.push_back(j);
        }
        return idx;
    }
};

/**
 * Precedence closure of a swap graph.
 *
 * @param H (R x R) swap graph; H(b, a) nonzero forces b before a
 */
template <class T>
PasPlacement<T> pas_placement(const Matrix<T>& H) {
    PasPlacement<T> out;
    if (H.rows() == 0 || H.cols() == 0) return out;
    if (H.rows() != H.cols()) throw InputError("pas_placement: the swap graph is not square");
    const std::size_t R = H.rows();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    std::vector<char> P(R * R, 0), Hb(R * R, 0);
    for (std::size_t i = 0; i < R; ++i)
        for (std::size_t j = 0; j < R; ++j) {
            const char b = H(i, j) != zero ? 1 : 0;
            Hb[i * R + j] = b;
            P[i * R + j] = b;
        }
    for (std::size_t it = 0; it < R; ++it) {
        std::vector<char> next(P);
        for (std::size_t i = 0; i < R; ++i)
            for (std::size_t k = 0; k < R; ++k) {
                if (!P[i * R + k]) continue;
                for (std::size_t j = 0; j < R; ++j)
                    if (Hb[k * R + j]) next[i * R + j] = 1;
            }
        if (next == P) break;
        P.swap(next);
    }

    out.P = Matrix<T>(R, R, zero);
    for (std::size_t i = 0; i < R; ++i)
        for (std::size_t j = 0; j < R; ++j)
            if (P[i * R + j]) out.P(i, j) = one;
    return out;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PAS_PLACEMENT_H
