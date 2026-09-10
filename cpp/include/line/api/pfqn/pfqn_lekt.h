/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_LEKT_H
#define LINE_API_PFQN_PFQN_LEKT_H

/**
 * The common corrected asymptotic expansion (LE-KT), computed on the cheaper side.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_lekt.m. The corrected logistic
 * expansion (pfqn_ble) and the corrected Knessl-Tier expansion (pfqn_bkt) are
 * ONE estimator, evaluated in M-1 and in R dimensions. With a think time their
 * stationary points are one point in dual coordinates,
 *
 *   xi_r = N_r / (Z_r + v u'L_r)   (the class throughputs of the LE fixed point)
 *   v u_k = 1 / (1 - U_k)          (the M/M/1 factor of the KT saddle)
 *
 * and Sylvester's identity det(I_R + C'C) = det(I_M + CC') exchanges the R x R
 * Hessian determinant for the M x M one, after which every 2 pi cancels on both
 * sides; the two agree to the accuracy of the saddle-point solvers. Without a
 * think time the LE branch integrates the radius exactly as Gamma(N+M) while KT
 * Laplaces it, so they differ by the constant (1 - log(2 pi)/2) - r(N+M), r the
 * Stirling remainder of a Gamma direction; the common estimator is defined as
 * the KT value, and the LE side here carries M (1 - log(2 pi)/2) - r(N+M) in
 * place of pfqn_ble's (M-1)(1 - log(2 pi)/2). See _kb/03-api-layer.md.
 *
 * ROUTE. The KT side is an R-dimensional convex solve and an R x R
 * determinant, the LE side an M-dimensional fixed point and an (M-1) x (M-1)
 * one, so the KT side is taken when R <= M, and whenever a class self-loops
 * (one nonzero demand and no think time), which pfqn_kt extracts exactly.
 *
 * ARITHMETIC. Inherited from both sides: gated on has_transcendental.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_bkt.h"
#include "line/api/pfqn/pfqn_ble.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_lekt, mirroring [Gn, lGn, route]. */
template <class T>
struct LektResult {
    T G;
    T lG;
    std::string route;  ///< "kt" or "le"
};

/** "kt" when R <= M or a class self-loops, "le" otherwise. */
template <class T>
std::string pfqn_lekt_route(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.rows(), R = L.cols();
    bool selfloop = false;
    if (R > 1) {
        for (std::size_t r = 0; r < R && !selfloop; ++r) {
            std::size_t nnz = 0;
            for (std::size_t i = 0; i < M; ++i)
                if (L(i, r) != zero) ++nnz;
            const T z = Z.empty() ? zero : Z[r];
            if (nnz == 1 && z == zero) selfloop = true;
        }
    }
    return (R <= M || selfloop) ? "kt" : "le";
}

/**
 * @param L (M x R) demands, @param N (R) population, @param Z (R) think times
 *          (pass an empty vector or all zeros for the Z = 0 branch)
 */
template <class T>
LektResult<T> pfqn_lekt(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_lekt requires transcendental arithmetic (Laplace approximation of an integral)");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.rows(), R = L.cols();
    std::vector<T> Zc = Z;
    if (Zc.empty()) Zc.assign(R, zero);

    LektResult<T> res;
    res.route = pfqn_lekt_route(L, N, Zc);
    if (res.route == "kt") {
        BktResult<T> k = pfqn_bkt(L, N, Zc);
        res.G = k.G;
        res.lG = k.lG;
        return res;
    }
    BleResult<T> b = pfqn_ble(L, N, Zc);
    res.G = b.G;
    res.lG = b.lG;

    T Ntot = zero, Lsum = zero, Zsum = zero;
    for (const T& x : N) Ntot += x;
    for (const T& x : Zc) Zsum += x;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Lsum += L(i, r);
    if (M == 0 || N.empty() || Ntot == zero || num_traits<T>::to_double(Lsum) < 1e-4) {
        return res;  // pfqn_ble's degenerate branch: the delay term is exact
    }
    const bool no_delay = num_traits<T>::to_double(Zsum) < lang::GlobalConstants::Zero;
    if (!no_delay) return res;
    // the Z = 0 branch of pfqn_ble counts M-1 directions; the common estimator
    // carries M kappa - r(N+M)
    const T half = num_traits<T>::from_double(0.5);
    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    const T eta = T(Ntot + num_traits<T>::from_int(static_cast<long>(M)));
    const T kappa = T(num_traits<T>::from_int(1) - log(twopi) / num_traits<T>::from_int(2));
    const T r = T(detail::num_lgamma<T>(eta) - (eta - half) * log(eta) + eta - half * log(twopi));
    res.lG += kappa - r;
    res.G = exp(res.lG);
    return res;
}

template <class T>
LektResult<T> pfqn_lekt(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_lekt(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_LEKT_H
