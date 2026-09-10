/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_TRANSIENT_SENS_H
#define LINE_API_MC_CTMC_TRANSIENT_SENS_H

/**
 * Sensitivity of the transient distribution of a CTMC to a scalar parameter.
 *
 * Templated port of matlab/src/api/mc/ctmc_transient_sens.m. Differentiating
 * the forward equations with respect to theta, with an initial vector that does
 * not depend on theta, gives Trivedi and Bobbio (2017), Eq. (9.82),
 *
 *   d/dt (dpi/dtheta) = (dpi/dtheta) Q + pi (dQ/dtheta),   dpi(0)/dtheta = 0.
 *
 * The sensitivity equation is driven by pi(t), so the two cannot be advanced
 * separately: state and sensitivity are integrated as ONE augmented system of
 * size 2n on a single adaptive grid, which is also what keeps the two
 * consistent at every returned time point.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC, for the same reason as ctmc_transient:
 * it is the same ode23 controller on a larger system.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_transient.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct TransientSensResult {
    std::vector<T> t;  ///< accepted time points, the first being t0
    Matrix<T> pi;      ///< distribution at each time point
    Matrix<T> dpi;     ///< sensitivity at each time point
};

/**
 * @param Q   generator
 * @param dQ  derivative of the generator with respect to theta, same size
 * @param pi0 initial distribution (row vector)
 * @param t0  initial time
 * @param t1  final time
 * @param rtol relative tolerance of the ODE integrator
 * @param atol absolute tolerance of the ODE integrator
 */
template <class T>
TransientSensResult<T> ctmc_transient_sens(const Matrix<T>& Q, const Matrix<T>& dQ,
                                           const std::vector<T>& pi0, const T& t0, const T& t1,
                                           double rtol = 1e-3, double atol = 1e-6) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_transient_sens requires transcendental arithmetic: it is the ode23 "
                  "controller of ctmc_transient applied to the augmented 2n system");
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_transient_sens: generator is not square");
    if (dQ.rows() != n || dQ.cols() != n)
        throw InputError("ctmc_transient_sens: dQ must have the same size as Q");
    if (pi0.size() != n) throw InputError("ctmc_transient_sens: pi0 has the wrong length");
    const T zero = num_traits<T>::from_int(0);

    std::vector<T> v0(2 * n, zero);
    for (std::size_t i = 0; i < n; ++i) v0[i] = pi0[i];

    std::vector<T> tv;
    std::vector<std::vector<T>> yv;
    detail::ode23<T>(
        [&Q, &dQ, n, &zero](const std::vector<T>& v, std::vector<T>& dv) {
            for (std::size_t j = 0; j < n; ++j) {
                T a = zero, b = zero;
                for (std::size_t i = 0; i < n; ++i) {
                    a += v[i] * Q(i, j);
                    b += v[n + i] * Q(i, j) + v[i] * dQ(i, j);
                }
                dv[j] = a;
                dv[n + j] = b;
            }
        },
        t0, t1, v0, rtol, atol, tv, yv);

    TransientSensResult<T> r;
    r.t = tv;
    r.pi = Matrix<T>(yv.size(), n);
    r.dpi = Matrix<T>(yv.size(), n);
    for (std::size_t k = 0; k < yv.size(); ++k)
        for (std::size_t j = 0; j < n; ++j) {
            r.pi(k, j) = yv[k][j];
            r.dpi(k, j) = yv[k][n + j];
        }
    return r;
}

/** Overload starting from the uniform distribution, as MATLAB's short forms do. */
template <class T>
TransientSensResult<T> ctmc_transient_sens(const Matrix<T>& Q, const Matrix<T>& dQ, const T& t0,
                                           const T& t1, double rtol = 1e-3, double atol = 1e-6) {
    const std::size_t n = Q.rows();
    if (n == 0) throw InputError("ctmc_transient_sens: empty generator");
    const std::vector<T> pi0(n, num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(n)));
    return ctmc_transient_sens(Q, dQ, pi0, t0, t1, rtol, atol);
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_TRANSIENT_SENS_H
