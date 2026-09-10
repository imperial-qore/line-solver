/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_MAP_MOMENTS_H
#define LINE_API_FES_MAP_MOMENTS_H

/**
 * Moments and index of dispersion of an inter-departure MAP.
 *
 * Templated port of matlab/src/api/fes/fes_map_moments.m and fes_map_euler.m,
 * mirrored by the JAR and native Python. Evaluates equations (4), (5) and (7) of
 * Casale, Mi, Cherkasova and Smirni, IEEE Trans. Soft. Eng. 37(5), 2011.
 *
 * The inverse (-T0)^-1 is dense even when T0 is sparse, so it is never formed:
 * the moments follow from the vector recursion v_{k+1} = v_k (-T0)^-1, each step
 * being one linear solve. Method "euler" replaces the solve by the quadrature
 * v (-T0)^-1 = v int_0^inf exp(T0 t) dt of Section 5.2.2, integrated by the
 * trapezoid rule with the Euler propagator exp(T0 dt) ~ I + T0 dt and a step
 * below the inverse of the largest diagonal element in absolute value, as in the
 * uniformization method. Method "ssolve" is the default because it is exact and
 * faster; "euler" reproduces the reference implementation of the paper and is
 * first order in the step.
 *
 * ARITHMETIC: "ssolve" uses field operations only and is exact at T = Rational;
 * "euler" is an approximation at every arithmetic.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace fes {

/** Descriptors a MAP(2) is fitted against. */
template <class T>
struct FesMapMoments {
    T e1;
    T e2;
    T e3;
    T e11;
    T idc;
};

namespace detail {

/** Solve the row system x A = b through the transpose. */
template <class T>
std::vector<T> solve_row(const Matrix<T>& A, const std::vector<T>& b) {
    Matrix<T> At(A.cols(), A.rows());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) At(j, i) = A(i, j);
    return solve(At, b);
}

}  // namespace detail

/**
 * Approximate v (-T0)^-1 by the trapezoid rule with the Euler propagator.
 *
 * @param v        row vector to be multiplied by (-T0)^-1
 * @param T0       hidden transitions of the MAP, a stable matrix
 * @param dt       integration step, below 1/max(abs(diag(T0)))
 * @param tol      relative mass left when the integration stops
 * @param iter_max maximum number of integration steps
 */
template <class T>
std::vector<T> fes_map_euler(const std::vector<T>& v, const Matrix<T>& T0, const T& dt, double tol,
                             std::size_t iter_max) {
    const T zero = num_traits<T>::from_int(0);
    const T half = num_traits<T>::from_rational(1, 2);
    std::vector<T> y(v.size(), zero);
    std::vector<T> z = v;
    double nrm0 = 0;
    for (const T& x : v) nrm0 += std::abs(num_traits<T>::to_double(x));

    for (std::size_t it = 0; it < iter_max; ++it) {
        const std::vector<T> zT0 = vecmul(z, T0);
        std::vector<T> znext(z.size());
        for (std::size_t i = 0; i < z.size(); ++i) znext[i] = z[i] + dt * zT0[i];
        for (std::size_t i = 0; i < y.size(); ++i) y[i] += dt * half * (z[i] + znext[i]);
        z = znext;
        double nrm = 0;
        for (const T& x : z) nrm += std::abs(num_traits<T>::to_double(x));
        if (nrm <= tol * nrm0) break;
    }
    return y;
}

/**
 * @param map         the pair (T0, T1) of the inter-departure MAP
 * @param method      "ssolve" for the linear solve, "euler" for the quadrature
 * @param step_safety fraction of the uniformization bound used as Euler step
 * @param tol         relative mass left when the Euler quadrature stops
 * @param iter_max    maximum number of Euler integration steps
 */
template <class T>
FesMapMoments<T> fes_map_moments(const mam::Map<T>& map, const std::string& method = "ssolve",
                                 double step_safety = 0.1, double tol = 1e-12,
                                 std::size_t iter_max = 1000000) {
    const Matrix<T>& T0 = map.D0;
    const Matrix<T>& T1 = map.D1;
    const std::size_t dim = T0.rows();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    const Matrix<T> Q = mam::map_infgen(map);
    const std::vector<T> phi = mc::ctmc_solve(Q);
    std::vector<T> pie = vecmul(phi, T1);
    T lambda = zero;
    for (const T& x : pie) lambda += x;
    for (std::size_t i = 0; i < pie.size(); ++i) pie[i] = pie[i] / lambda;

    Matrix<T> negT0(dim, dim, zero);
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim; ++j) negT0(i, j) = -T0(i, j);

    const bool euler = (method == "euler");
    if (!euler && method != "ssolve")
        throw InputError("fes_map_moments: unknown method, use ssolve or euler");

    T dt = zero;
    if (euler) {
        double dmax = 0;
        for (std::size_t i = 0; i < dim; ++i) {
            const double d = std::abs(num_traits<T>::to_double(T0(i, i)));
            if (d > dmax) dmax = d;
        }
        dt = num_traits<T>::from_double(step_safety / dmax);
    }

    const std::vector<T> v1 = euler ? fes_map_euler(pie, T0, dt, tol, iter_max)
                                    : detail::solve_row(negT0, pie);
    const std::vector<T> v2 = euler ? fes_map_euler(v1, T0, dt, tol, iter_max)
                                    : detail::solve_row(negT0, v1);
    const std::vector<T> v3 = euler ? fes_map_euler(v2, T0, dt, tol, iter_max)
                                    : detail::solve_row(negT0, v2);
    const std::vector<T> v2T1 = vecmul(v2, T1);
    const std::vector<T> v4 = euler ? fes_map_euler(v2T1, T0, dt, tol, iter_max)
                                    : detail::solve_row(negT0, v2T1);

    FesMapMoments<T> out;
    out.e1 = zero;
    out.e2 = zero;
    out.e3 = zero;
    out.e11 = zero;
    for (std::size_t i = 0; i < dim; ++i) {
        out.e1 += v1[i];
        out.e2 += v2[i];
        out.e3 += v3[i];
        out.e11 += v4[i];
    }
    out.e2 = num_traits<T>::from_int(2) * out.e2;
    out.e3 = num_traits<T>::from_int(6) * out.e3;

    // equation (7), with pie inv(Q + e phi) obtained from the rank-one update
    // y Q = pie - phi under the normalization y e = 1
    Matrix<T> A = Q;
    for (std::size_t i = 0; i < dim; ++i) A(i, dim - 1) = one;
    std::vector<T> rhs(dim);
    for (std::size_t i = 0; i < dim; ++i) rhs[i] = pie[i] - phi[i];
    rhs[dim - 1] = one;
    const std::vector<T> y = detail::solve_row(A, rhs);
    const std::vector<T> yT1 = vecmul(y, T1);
    T s = zero;
    for (const T& x : yT1) s += x;
    out.idc = one + num_traits<T>::from_int(2) * (lambda - s);
    return out;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_MAP_MOMENTS_H
