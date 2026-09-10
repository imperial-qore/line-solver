/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_LQN_JACOBIAN_H
#define LINE_API_INFER_INFER_LQN_JACOBIAN_H

/**
 * Forward finite-difference sensitivity matrix of an observation map.
 *
 * Templated port of matlab/src/api/infer/infer_lqn_jacobian.m. The JAR carries
 * the same computation inline inside jline/api/infer/InferLqn.java rather than
 * as a separate entry point.
 *
 * H = dh/da at the parameter vector a, column by column,
 *   H(:,i) = (h(a + d_i e_i) - h(a)) / d_i,  d_i = fd_step max(|a_i|, fd_floor)
 * which is the approximate sensitivity matrix H_k of the EKF update of Zheng,
 * Yang, Woodside, Litoiu, Iszlai, "Tracking Time-Varying Parameters in
 * Software Systems with Extended Kalman Filters", CASCON 2005. h is evaluated
 * numel(a)+1 times, and h0 = h(a) is returned alongside because the caller
 * (the EKF innovation) always needs it too.
 *
 * The observation map is a std::function, so the rest of infer_lqn -- which
 * needs the LayeredNetwork model layer to evaluate h -- stays out of this
 * header: any caller that can evaluate h can use this.
 *
 * ARITHMETIC: additions, one division per column and a comparison, so a finite
 * field computation and no transcendental gate. Note what that does and does
 * not buy: in the exact instantiation the returned matrix is the exact
 * DIFFERENCE QUOTIENT of h, not the derivative -- the truncation error of the
 * forward difference is a property of the formula, not of the arithmetic.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/** Mirrors the [H, h0] return list of the MATLAB function. */
template <class T>
struct JacobianResult {
    Matrix<T> H;      ///< (no x np) sensitivity matrix
    std::vector<T> h0;  ///< (no) observation at the base point
};

/**
 * @param hfun     observation map, a parameter vector to an observation vector
 * @param a        (np) point at which the sensitivity is taken
 * @param fd_step  relative perturbation, MATLAB's default 1e-3
 * @param fd_floor minimum absolute perturbation scale, MATLAB's default 1e-6
 */
template <class T>
JacobianResult<T> infer_lqn_jacobian(
    const std::function<std::vector<T>(const std::vector<T>&)>& hfun, const std::vector<T>& a,
    double fd_step = 1e-3, double fd_floor = 1e-6) {
    const std::size_t np = a.size();
    if (np == 0) throw InputError("infer_lqn_jacobian: empty parameter vector");
    if (!(fd_step > 0.0)) throw InputError("infer_lqn_jacobian: the step must be positive");

    JacobianResult<T> out;
    out.h0 = hfun(a);
    const std::size_t no = out.h0.size();
    out.H = Matrix<T>(no, np, num_traits<T>::from_int(0));

    const T step = num_traits<T>::from_double(fd_step);
    const T floor = num_traits<T>::from_double(fd_floor);
    for (std::size_t i = 0; i < np; ++i) {
        const T mag = num_abs(T(a[i]));
        const T d = step * (mag > floor ? mag : floor);
        std::vector<T> ap = a;
        ap[i] = ap[i] + d;
        const std::vector<T> hi = hfun(ap);
        if (hi.size() != no)
            throw InputError("infer_lqn_jacobian: the observation map changed its output length");
        for (std::size_t k = 0; k < no; ++k) out.H(k, i) = (hi[k] - out.h0[k]) / d;
    }
    return out;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_LQN_JACOBIAN_H
