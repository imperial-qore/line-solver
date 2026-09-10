/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MM1K_LOSS_H
#define LINE_API_QSYS_QSYS_MM1K_LOSS_H

/**
 * Blocking probability of the M/M/1/K queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_mm1k_loss.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_mm1k_loss.java (identical).
 *
 *   rho  = lambda/mu
 *   Ploss = (1-rho)/(1-rho^(K+1)) * rho^K
 *
 * Only integer powers appear, so this is exact for T = Rational. It is the
 * closed form the Niu-Cooper transform-free M/G/1/K analysis collapses onto
 * when the service is exponential, and qsys_mg1k_loss must reproduce it.
 *
 * The formula has a removable singularity at rho = 1, where the true value is
 * 1/(K+1). MATLAB divides by zero and returns NaN there; the port raises
 * instead, since in an exact field the quotient has no value at all.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

template <class T>
struct Mm1kLossResult {
    T lossProbability;
    T utilization;  ///< rho = lambda/mu, the offered load (not the carried one)
};

/**
 * @param lambda arrival rate
 * @param mu     service rate
 * @param K      system capacity, jobs in service included, K >= 1
 */
template <class T>
Mm1kLossResult<T> qsys_mm1k_loss(const T& lambda, const T& mu, unsigned K) {
    if (K == 0) throw InputError("qsys_mm1k_loss: K must be at least 1");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    if (one - num_pow_int(rho, K + 1) == num_traits<T>::from_int(0))
        throw InputError("qsys_mm1k_loss: rho == 1, the closed form is a removable singularity");
    Mm1kLossResult<T> r;
    r.utilization = rho;
    r.lossProbability = (one - rho) / (one - num_pow_int(rho, K + 1)) * num_pow_int(rho, K);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MM1K_LOSS_H
