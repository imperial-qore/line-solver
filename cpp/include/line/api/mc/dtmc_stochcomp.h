/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_DTMC_STOCHCOMP_H
#define LINE_API_MC_DTMC_STOCHCOMP_H

/**
 * Stochastic complement of a DTMC partition, a port of
 * matlab/lib/kpctoolbox/mc/dtmc_stochcomp.m.
 *
 * For a row-stochastic P partitioned into a retained set I and its complement
 * Ic, the stochastic complement over I is
 *   S = P11 + P12 (Id - P22)^-1 P21,
 * the routing seen by an observer who watches only the states in I, censoring
 * the excursions through Ic. It is the same construction `NetworkStruct`'s
 * `station_routing` performs over the stateful nodes, exposed here as a free
 * function so the cacheqn driver can complement a rtnodes matrix it has
 * rewritten in place (relabelling a Cache as a class switch) without going back
 * through `route_eff`.
 *
 * ARITHMETIC: field. The only operation is the linear solve (Id - P22) X = P21
 * by Gaussian elimination with partial pivoting, so it is exact under Rational.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * @param P    (n x n) row-stochastic transition matrix
 * @param keep the 0-based indices of the states to retain (the set I)
 * @return the (|I| x |I|) stochastic complement S over the retained states, in
 *         the order given by `keep`
 */
template <class T>
Matrix<T> dtmc_stochcomp(const Matrix<T>& P, const std::vector<std::size_t>& keep) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t n = P.rows();
    if (P.cols() != n) throw InputError("dtmc_stochcomp: the matrix is not square");

    std::vector<bool> kept(n, false);
    for (std::size_t i : keep) {
        if (i >= n) throw InputError("dtmc_stochcomp: a retained index is out of range");
        kept[i] = true;
    }
    std::vector<std::size_t> drop;
    for (std::size_t i = 0; i < n; ++i)
        if (!kept[i]) drop.push_back(i);

    const std::size_t nk = keep.size(), nd = drop.size();
    Matrix<T> P11(nk, nk, zero);
    for (std::size_t a = 0; a < nk; ++a)
        for (std::size_t b = 0; b < nk; ++b) P11(a, b) = P(keep[a], keep[b]);
    if (nd == 0) return P11;

    Matrix<T> P12(nk, nd, zero), P21(nd, nk, zero), A(nd, nd, zero);
    for (std::size_t a = 0; a < nk; ++a)
        for (std::size_t b = 0; b < nd; ++b) P12(a, b) = P(keep[a], drop[b]);
    for (std::size_t a = 0; a < nd; ++a) {
        for (std::size_t b = 0; b < nk; ++b) P21(a, b) = P(drop[a], keep[b]);
        for (std::size_t b = 0; b < nd; ++b) A(a, b) = T((a == b ? one : zero) - P(drop[a], drop[b]));
    }

    // X = (Id - P22)^-1 P21 by Gaussian elimination with partial pivoting.
    Matrix<T> X = P21;
    for (std::size_t col = 0; col < nd; ++col) {
        std::size_t best = col;
        double bv = std::fabs(num_traits<T>::to_double(A(col, col)));
        for (std::size_t r = col + 1; r < nd; ++r) {
            const double v = std::fabs(num_traits<T>::to_double(A(r, col)));
            if (v > bv) { bv = v; best = r; }
        }
        if (best != col) {
            for (std::size_t b = 0; b < nd; ++b) std::swap(A(col, b), A(best, b));
            for (std::size_t b = 0; b < nk; ++b) std::swap(X(col, b), X(best, b));
        }
        if (A(col, col) == zero)
            throw NumericError("dtmc_stochcomp: the complement block is singular");
        for (std::size_t r = 0; r < nd; ++r) {
            if (r == col) continue;
            const T f = T(A(r, col) / A(col, col));
            if (f == zero) continue;
            for (std::size_t b = 0; b < nd; ++b) A(r, b) = T(A(r, b) - f * A(col, b));
            for (std::size_t b = 0; b < nk; ++b) X(r, b) = T(X(r, b) - f * X(col, b));
        }
    }
    for (std::size_t r = 0; r < nd; ++r)
        for (std::size_t b = 0; b < nk; ++b) X(r, b) = T(X(r, b) / A(r, r));

    Matrix<T> S = P11;
    for (std::size_t a = 0; a < nk; ++a)
        for (std::size_t b = 0; b < nk; ++b) {
            T acc = zero;
            for (std::size_t d = 0; d < nd; ++d) acc = T(acc + P12(a, d) * X(d, b));
            S(a, b) = T(S(a, b) + acc);
        }
    return S;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_DTMC_STOCHCOMP_H
