/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MXM1_H
#define LINE_API_QSYS_QSYS_MXM1_H

/**
 * M^X/M/1: the batch-arrival queue with exponential service.
 *
 * Templated port of matlab/src/api/qsys/qsys_mxm1.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mxm1.java.
 *
 *   lambda = lambda_batch E[X],   rho = lambda/mu
 *   Wq = rho/(mu(1-rho)) + (E[X^2]-E[X]) / (2 mu E[X] (1-rho))
 *   W  = Wq + 1/mu,   Q = lambda W
 *
 * The first term is the M/M/1 delay of the batch stream and the second is the
 * delay a job suffers behind its own batch mates. Only the first two moments
 * of the batch size enter, and only rationally, so this is exact for
 * T = Rational.
 *
 * MATLAB dispatches on the argument shapes to accept the batch law as
 * (E[X], E[X^2]), as (E[X], Var[X]) with a 'variance' flag, or as a support
 * with a pmf. C++ overloading cannot see MATLAB's shape test, so the three
 * forms are separate named entry points: qsys_mxm1, qsys_mxm1_variance and
 * qsys_mxm1_pmf. All three funnel into the same closed form.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

template <class T>
struct MxM1Result {
    T W;   ///< mean time in system
    T Wq;  ///< mean waiting time in queue
    T U;   ///< server utilization rho
    T Q;   ///< mean number in system, lambda W
};

/**
 * @param lambda_batch batch arrival rate
 * @param mu           service rate
 * @param E_X          mean batch size
 * @param E_X2         second raw moment of the batch size
 */
template <class T>
MxM1Result<T> qsys_mxm1(const T& lambda_batch, const T& mu, const T& E_X, const T& E_X2) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T lambda = lambda_batch * E_X;
    const T rho = lambda / mu;
    if (rho >= one) throw InputError("qsys_mxm1: system is unstable, rho >= 1");
    MxM1Result<T> r;
    r.Wq = rho / (mu * (one - rho)) + (E_X2 - E_X) / (two * mu * E_X * (one - rho));
    r.W = r.Wq + one / mu;
    r.U = rho;
    r.Q = lambda * r.W;
    return r;
}

/** Variance form: MATLAB's qsys_mxm1(..., Var_X, 'variance'). */
template <class T>
MxM1Result<T> qsys_mxm1_variance(const T& lambda_batch, const T& mu, const T& E_X,
                                 const T& Var_X) {
    return qsys_mxm1(lambda_batch, mu, E_X, T(Var_X + E_X * E_X));
}

/**
 * Support-and-pmf form: MATLAB's qsys_mxm1(..., batch_sizes, pmf). The pmf is
 * renormalized, as in MATLAB, so an unnormalized weight vector is accepted.
 */
template <class T>
MxM1Result<T> qsys_mxm1_pmf(const T& lambda_batch, const T& mu,
                            const std::vector<T>& batch_sizes, const std::vector<T>& pmf) {
    if (batch_sizes.size() != pmf.size())
        throw InputError("qsys_mxm1: batch sizes and pmf must have the same length");
    if (batch_sizes.empty()) throw InputError("qsys_mxm1: empty batch-size support");
    const T zero = num_traits<T>::from_int(0);
    T tot = zero;
    for (const T& v : pmf) tot += v;
    if (tot == zero) throw InputError("qsys_mxm1: pmf sums to zero");
    T E_X = zero, E_X2 = zero;
    for (std::size_t i = 0; i < pmf.size(); ++i) {
        const T w = pmf[i] / tot;
        E_X += batch_sizes[i] * w;
        E_X2 += batch_sizes[i] * batch_sizes[i] * w;
    }
    return qsys_mxm1(lambda_batch, mu, E_X, E_X2);
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MXM1_H
