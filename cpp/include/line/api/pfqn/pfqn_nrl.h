/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_NRL_H
#define LINE_API_PFQN_PFQN_NRL_H

/**
 * Norlund-Rice inversion of the normalizing constant, in its logistic (NRL)
 * and probit (NRP) substitutions.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_nrl.m and pfqn_nrp.m. Both scale
 * the demands column-wise into [0,1], apply the Laplace approximation of
 * laplaceapprox.m to the corresponding infradius integrand at the origin, and
 * undo the scaling with sum_r N_r log Lmax_r. The two differ only in the
 * change of variable, so they are ported together.
 *
 * THE CURVATURE TERM IS HALF THE LOG-DETERMINANT, and it is the single thing
 * most worth checking when a ported value disagrees with MATLAB here. Both
 * routines consume laplaceapprox's logI, which is log(I) and therefore carries
 * -(1/2) log det(-H); using the full log-determinant instead shifts every
 * value by (1/2) log det(-H) with no other symptom. On the single-class model
 * L = [1/2, 1/3, 1/5], N = 5, Z = 0 that shift is 0.347 nats, turning the
 * correct -2.1159 into -1.7694, both of which look plausible next to the exact
 * -1.9951.
 *
 * DELAY. A non-zero think time is appended as one more station whose rate row
 * is 1, 2, ..., Ntot, i.e. an infinite server, exactly as in the reference.
 *
 * ARITHMETIC. A Laplace approximation of a contour integral, so gated on
 * num_traits<T>::has_transcendental. pfqn_nrp additionally carries the
 * double-precision ceiling of the normal CDF; see infradius_h.h.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/pfqn/infradius_h.h"
#include "line/api/pfqn/laplaceapprox.h"
#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** Shared body of pfqn_nrl and pfqn_nrp; `probit` selects the substitution. */
template <class T>
T pfqn_nr_impl(const Matrix<T>& L0, const std::vector<T>& N, const std::vector<T>& Z,
               const Matrix<T>& alpha0, bool probit) {
    using std::log;
    const std::size_t M = L0.rows(), R = L0.cols();
    if (N.size() != R) throw InputError("pfqn_nr: L and N disagree on the class count");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T Nt = zero, Zsum = zero;
    for (const T& v : N) Nt += v;
    for (const T& v : Z) Zsum += v;
    if (Nt < zero) throw InputError("pfqn_nr: negative population");
    if (Nt == zero) return zero;
    const std::size_t Ntot = static_cast<std::size_t>(num_traits<T>::to_double(Nt));

    // Append the delay as an infinite-server station.
    const std::size_t Mx = M + (Zsum > zero ? 1 : 0);
    Matrix<T> L(Mx, R);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) L(i, r) = L0(i, r);
    Matrix<T> alpha(Mx, Ntot);
    for (std::size_t i = 0; i < M && i < alpha0.rows(); ++i)
        for (std::size_t k = 0; k < Ntot; ++k)
            alpha(i, k) = k < alpha0.cols() ? alpha0(i, k) : one;
    if (Zsum > zero) {
        if (Z.size() != R) throw InputError("pfqn_nr: Z has the wrong length");
        for (std::size_t r = 0; r < R; ++r) L(M, r) = Z[r];
        for (std::size_t k = 0; k < Ntot; ++k)
            alpha(M, k) = num_traits<T>::from_int(static_cast<long>(k) + 1);
    }

    // A single queueing station with no delay has a closed form.
    if (M == 1 && Zsum == zero) {
        T lG = detail::num_factln<T>(Nt);
        for (std::size_t r = 0; r < R; ++r) {
            lG -= detail::num_factln<T>(N[r]);
            if (L(0, r) > zero) lG += N[r] * log(L(0, r));
        }
        for (std::size_t k = 0; k < Ntot; ++k) lG -= log(alpha(0, k));
        return lG;
    }

    // Scale each class column by its largest demand.
    std::vector<T> Lmax(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        T m = L(0, r);
        for (std::size_t i = 1; i < Mx; ++i)
            if (L(i, r) > m) m = L(i, r);
        if (m <= zero) throw InputError("pfqn_nr: a class has no positive demand");
        Lmax[r] = m;
    }
    Matrix<T> Ls(Mx, R);
    for (std::size_t i = 0; i < Mx; ++i)
        for (std::size_t r = 0; r < R; ++r) Ls(i, r) = T(L(i, r) / Lmax[r]);

    const std::function<T(const std::vector<T>&)> h = [&](const std::vector<T>& x) {
        return probit ? infradius_hnorm(x, Ls, N, alpha) : infradius_h(x, Ls, N, alpha);
    };
    const std::vector<T> x0(R, zero);
    const LaplaceResult<T> la = laplaceapprox<T>(h, x0);
    T lG = la.logI;
    for (std::size_t r = 0; r < R; ++r) lG += N[r] * log(Lmax[r]);
    return lG;
}

}  // namespace detail

/**
 * Norlund-Rice logit approximation of log G.
 *
 * @param L     (M x R) demands
 * @param N     (R) population
 * @param Z     (R) think times, empty for zero
 * @param alpha (M x Ntot) load-dependent rates, empty for all ones
 */
template <class T>
T pfqn_nrl(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
           const Matrix<T>& alpha) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_nrl requires transcendental arithmetic (Laplace approximation of a contour "
                  "integral)");
    return detail::pfqn_nr_impl(L, N, Z, alpha, false);
}

/** Norlund-Rice probit approximation of log G. */
template <class T>
T pfqn_nrp(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
           const Matrix<T>& alpha) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_nrp requires transcendental arithmetic (Laplace approximation of a contour "
                  "integral)");
    return detail::pfqn_nr_impl(L, N, Z, alpha, true);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_NRL_H
