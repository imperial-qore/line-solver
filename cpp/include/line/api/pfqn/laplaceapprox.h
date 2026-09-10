/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LAPLACEAPPROX_H
#define LINE_API_PFQN_LAPLACEAPPROX_H

/**
 * Laplace approximation of a multidimensional integral around a given point.
 *
 * Templated port of matlab/src/api/pfqn/laplaceapprox.m together with the two
 * thirdparty helpers it calls, matlab/lib/thirdparty/num_hess.m and num_grad.m
 * (central differences at step h, the Hessian built as a difference of
 * gradients). Those helpers have no pfqn counterpart of their own, so they are
 * ported here as num_grad / num_hess in the pfqn::detail namespace rather than
 * given a header of their own.
 *
 *   I = h(x0) sqrt((2 pi)^d / det(-H)),  H = Hessian of log h at x0.
 *
 * MATLAB's retry ladder on a negative det(-H) -- widen the differentiation
 * step from 1e-5 to 1e-4 to 1e-3, then warn and carry on -- is reproduced
 * exactly, including the fact that the final value is still returned when the
 * determinant stays negative. Callers (pfqn_nrl, pfqn_nrp) take the real part
 * of the logarithm afterwards, which is what makes that path survivable.
 *
 * logI IS the logarithm of I. laplaceapprox.m returns
 * logI = log(h(x0)) + (d/2) log(2 pi) - (1/2) log(det(-H)), i.e. exactly
 * log(I): the curvature term carries HALF the log-determinant, because I
 * carries sqrt(det(-H)) in its denominator. That is worth stating explicitly
 * because the two are easy to desynchronize -- the reference itself carried a
 * full log(detnH) at one point, which shifted every pfqn_nrl and pfqn_nrp
 * value by (1/2) log det(-H), and pfqn_nrl / pfqn_nrp read logI, not I.
 *
 * ARITHMETIC. Both the (2 pi)^{d/2} factor and the finite-difference Hessian
 * are inexact, so the routine is gated on num_traits<T>::has_transcendental.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** Central-difference gradient at step h (matlab/lib/thirdparty/num_grad.m). */
template <class T>
std::vector<T> num_grad(const std::function<T(const std::vector<T>&)>& f, const std::vector<T>& X,
                        const T& h) {
    std::vector<T> df(X.size());
    for (std::size_t i = 0; i < X.size(); ++i) {
        std::vector<T> x1 = X, x2 = X;
        x1[i] = T(X[i] - h);
        x2[i] = T(X[i] + h);
        df[i] = T(T(f(x2) - f(x1)) / T(num_traits<T>::from_int(2) * h));
    }
    return df;
}

/** Central-difference Hessian at step h (matlab/lib/thirdparty/num_hess.m). */
template <class T>
Matrix<T> num_hess(const std::function<T(const std::vector<T>&)>& f, const std::vector<T>& X,
                   const T& h) {
    const std::size_t n = X.size();
    Matrix<T> H(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<T> x1 = X, x2 = X;
        x1[i] = T(X[i] - h);
        x2[i] = T(X[i] + h);
        const std::vector<T> df1 = num_grad<T>(f, x1, h);
        const std::vector<T> df2 = num_grad<T>(f, x2, h);
        for (std::size_t j = 0; j < n; ++j)
            H(i, j) = T(T(df2[j] - df1[j]) / T(num_traits<T>::from_int(2) * h));
    }
    return H;
}

}  // namespace detail

/** Return value of laplaceapprox, mirroring [I, H, logI]. */
template <class T>
struct LaplaceResult {
    T I;
    Matrix<T> H;
    T logI;
    bool detNegative;  ///< true when det(-H) stayed negative (MATLAB warns)
};

/**
 * @param h  integrand, evaluated as a callable on a d-vector
 * @param x0 expansion point
 */
template <class T>
LaplaceResult<T> laplaceapprox(const std::function<T(const std::vector<T>&)>& h,
                               const std::vector<T>& x0) {
    static_assert(num_traits<T>::has_transcendental,
                  "laplaceapprox requires transcendental arithmetic (log, sqrt, (2 pi)^{d/2})");
    using std::log;
    using std::sqrt;
    const std::size_t d = x0.size();
    if (d == 0) throw InputError("laplaceapprox: empty expansion point");
    const std::function<T(const std::vector<T>&)> logh = [&h](const std::vector<T>& x) {
        using std::log;
        return T(log(h(x)));
    };

    const double steps[3] = {1e-5, 1e-4, 1e-3};
    Matrix<T> H;
    T detnH = num_traits<T>::from_int(0);
    bool neg = true;
    for (int k = 0; k < 3; ++k) {
        H = detail::num_hess<T>(logh, x0, num_traits<T>::from_double(steps[k]));
        Matrix<T> nH(d, d);
        for (std::size_t i = 0; i < d; ++i)
            for (std::size_t j = 0; j < d; ++j) nH(i, j) = T(-H(i, j));
        detnH = detail::pfqn_det(nH);
        if (detnH >= num_traits<T>::from_int(0)) {
            neg = false;
            break;
        }
    }

    LaplaceResult<T> r;
    r.H = H;
    r.detNegative = neg;
    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    const T h0 = h(x0);
    // negative det(-H) rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    const T adet = num_abs(detnH);
    if (adet == num_traits<T>::from_int(0)) throw NumericError("laplaceapprox: singular Hessian");
    r.I = neg ? num_traits<T>::from_int(0)
              : T(h0 * sqrt(T(num_pow_int(twopi, static_cast<unsigned>(d)) / detnH)));
    r.logI = T(log(h0) + T(num_traits<T>::from_rational(static_cast<long>(d), 2) * log(twopi)) -
               T(num_traits<T>::from_rational(1, 2) * log(adet)));
    return r;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LAPLACEAPPROX_H
