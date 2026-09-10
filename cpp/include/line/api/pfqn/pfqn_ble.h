/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_BLE_H
#define LINE_API_PFQN_PFQN_BLE_H

/**
 * Logistic expansion with the eps->0 bias correction (BLE).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ble.m. Cas17 Theorem 4.1 holds for
 * eps >= eps_N > 0; the K(1+eps*N) self-looping populations are what make the
 * integrand concentrate. Evaluated at eps->0, as pfqn_le does, the curvature at
 * the saddle tends to 1 rather than growing with N, so Laplace's method has no
 * asymptotic regime there and carries an O(1) relative bias of e/sqrt(2 pi) PER
 * LAPLACED DIRECTION. The count is the exponent on sqrt(2 pi) in the branch
 * taken: M-1 with Z = 0, where the radial integral is exact as Gamma(N+M), and M
 * with Z > 0, where the radius is Laplaced too. Measured over the 1562 models of
 * the Cas17 dataset (Zenodo 546873, sec5.3.1, sigma = 100) the Z > 0 deficit is M
 * to within 0.01 units. The published expansion is NOT in error and the
 * correction is EMPIRICAL, not part of Cas17; see _kb/03-api-layer.md.
 *
 * ARITHMETIC. Inherited from pfqn_le: the correction is a logarithm, so the
 * routine is gated on num_traits<T>::has_transcendental for the same reason.
 *
 * DEGENERATE BRANCH. When there is no queueing station to expand around, pfqn_le
 * returns the exact delay term and there is no Laplace step to correct, so
 * pfqn_ble returns it unchanged; this is the MATLAB branching.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_le.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_ble, mirroring [Gn, lGn]. */
template <class T>
using BleResult = LeResult<T>;

/**
 * Logistic expansion estimate of the normalizing constant, bias-corrected.
 *
 * @param L (M x R) demands, @param N (R) population, @param Z (R) think times
 *          (pass an empty vector or all zeros for the Z = 0 branch)
 */
template <class T>
BleResult<T> pfqn_ble(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_ble requires transcendental arithmetic (Laplace approximation of an integral)");
    using std::exp;
    using std::log;
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0);

    BleResult<T> res = pfqn_le(L, N, Z);

    T Ntot = zero, Lsum = zero, Zsum = zero;
    for (const T& x : N) Ntot += x;
    for (const T& x : Z) Zsum += x;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Lsum += L(i, r);
    if (M == 0 || N.empty() || Ntot == zero || num_traits<T>::to_double(Lsum) < 1e-4) {
        return res;
    }

    // Same predicate pfqn_le branches on, so the count always matches the branch.
    const bool no_delay =
        Z.empty() || num_traits<T>::to_double(Zsum) < lang::GlobalConstants::Zero;
    const long n_gauss = static_cast<long>(M) - (no_delay ? 1 : 0);
    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    res.lG += num_traits<T>::from_int(n_gauss) *
              T(num_traits<T>::from_int(1) - log(twopi) / num_traits<T>::from_int(2));
    res.G = exp(res.lG);
    return res;
}

template <class T>
BleResult<T> pfqn_ble(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_ble(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_BLE_H
