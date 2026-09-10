/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MUSHIFT_H
#define LINE_API_PFQN_MUSHIFT_H

/**
 * Shift the load-dependent service-rate lattice of selected stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mushift.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/ld/Pfqn_mushift.java.
 *
 * The load-dependent convolution and MVA recursions need mu(i, n+1) when a job
 * is added at station i. This drops the first column of that station's rate
 * row and the last column of every other row, giving an (M x N-1) lattice in
 * which station i is "one job ahead":
 *
 *   mushifted(i, j) = mu(i, j+1)   for the shifted stations
 *   mushifted(m, j) = mu(m, j)     otherwise,          j = 1 .. N-1
 *
 * Arithmetic: EXACT-CAPABLE. The routine only copies entries, so it is exact
 * in every arithmetic and carries no transcendental gate.
 *
 * Note on the MATLAB loop. The reference recomputes the whole matrix inside a
 * loop over iset, so with more than one index only the LAST one is actually
 * shifted; every earlier one is overwritten. That is a defect, not a
 * convention -- the routine's own name and its single caller (one station at a
 * time) say each listed station should be shifted -- so this port shifts every
 * station in iset. With the single-element iset that the reference is called
 * with, the two agree exactly.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param mu   (M x N) rate lattice, N >= 1
 * @param iset 0-based station indices to shift
 * @return (M x N-1) shifted lattice
 */
template <class T>
Matrix<T> pfqn_mushift(const Matrix<T>& mu, const std::vector<std::size_t>& iset) {
    const std::size_t M = mu.rows();
    const std::size_t N = mu.cols();
    if (N < 1) throw InputError("pfqn_mushift: the rate lattice is empty");
    for (std::size_t i : iset)
        if (i >= M) throw InputError("pfqn_mushift: station index out of range");

    std::vector<bool> shift(M, false);
    for (std::size_t i : iset) shift[i] = true;

    Matrix<T> out(M, N - 1);
    for (std::size_t m = 0; m < M; ++m)
        for (std::size_t j = 0; j + 1 < N; ++j) out(m, j) = shift[m] ? mu(m, j + 1) : mu(m, j);
    return out;
}

/** Single-station overload, the form the reference is actually called with. */
template <class T>
Matrix<T> pfqn_mushift(const Matrix<T>& mu, std::size_t i) {
    return pfqn_mushift(mu, std::vector<std::size_t>{i});
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MUSHIFT_H
