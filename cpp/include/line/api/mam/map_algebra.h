/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_ALGEBRA_H
#define LINE_API_MAM_MAP_ALGEBRA_H

/**
 * MAP algebra missing from the moment and transform headers: time reversal, the
 * Kronecker product composition, the subdominant eigenvalue of the embedded
 * chain and the large-order threshold.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_timereverse.m, map_kpc.m,
 * map_gamma2.m and map_largemap.m.
 *
 * map_kpc composes two MAPs into one of order na*nb whose autocorrelation
 * decays with the PRODUCT of the two decay rates; note the sign, D0 of the
 * composition is MINUS the Kronecker product of the two D0 blocks, because
 * kron of two matrices with negative diagonals has a positive one.
 *
 * map_gamma2 is the second largest eigenvalue in modulus of the embedded
 * chain, which is the geometric decay rate of the autocorrelation of a
 * second-order MAP and its leading term in general. It is genuinely COMPLEX for
 * a MAP with oscillating correlation, so it is returned as such rather than as
 * its modulus.
 */

#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Order above which a MAP counts as large for the fitting heuristics. */
inline std::size_t map_largemap() { return 100; }

/** Time-reversed MAP, diag(pi)^-1 M' diag(pi) applied to D0 and D1. */
template <class T>
Map<T> map_timereverse(const Map<T>& m) {
    const std::size_t n = m.D0.rows();
    const std::vector<T> piq = map_prob(m);
    Map<T> out;
    out.D0 = Matrix<T>(n, n, num_traits<T>::from_int(0));
    out.D1 = Matrix<T>(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            out.D0(i, j) = m.D0(j, i) * piq[j] / piq[i];
            out.D1(i, j) = m.D1(j, i) * piq[j] / piq[i];
        }
    return out;
}

/** Kronecker product composition of two MAPs. */
template <class T>
Map<T> map_kpc(const Map<T>& a, const Map<T>& b) {
    Matrix<T> D0 = kron(a.D0, b.D0);
    for (std::size_t i = 0; i < D0.rows(); ++i)
        for (std::size_t j = 0; j < D0.cols(); ++j) D0(i, j) = -D0(i, j);
    return Map<T>{D0, kron(a.D1, b.D1)};
}

/** Left-folded composition of a whole list of MAPs. */
template <class T>
Map<T> map_kpc(const std::vector<Map<T>>& maps) {
    if (maps.size() < 2) throw InputError("map_kpc: at least two MAPs are required");
    Map<T> out = map_kpc(maps[0], maps[1]);
    for (std::size_t k = 2; k < maps.size(); ++k) out = map_kpc(out, maps[k]);
    return out;
}

/** Subdominant eigenvalue of the embedded chain, the leading ACF decay rate. */
inline std::complex<double> map_gamma2(const Map<double>& m) {
    const Matrix<double> P = map_embedded(m);
    std::vector<std::complex<double>> ev = eig_values(P);
    if (ev.size() < 2) throw InputError("map_gamma2: the MAP must have order at least 2");
    std::sort(ev.begin(), ev.end(),
              [](const std::complex<double>& x, const std::complex<double>& y) {
                  return std::abs(x) > std::abs(y);
              });
    return ev[1];
}

}  // namespace mam
}  // namespace line

#endif
