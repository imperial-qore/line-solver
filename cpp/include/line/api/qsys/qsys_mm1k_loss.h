/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MM1K_LOSS_H
#define LINE_API_QSYS_QSYS_MM1K_LOSS_H

/**
 * @file
 * @ingroup api_qsys
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
 * 1/(K+1): the stationary law is uniform over 0..K, so the full state carries
 * 1/(K+1) like every other. This port RAISED there until 2026-09-13, while
 * MATLAB, the JAR and python all returned the limit, so the one caller that
 * reaches a saturated M/M/1/K got an exception from C++ and a number from the
 * other three. The limit is a rational, not an approximation, so it is exact
 * for T = Rational as well.
 */

#include <cmath>

#include "line/api/qsys/qsys_types.h"
#include "line/lang/lang_types.h"
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
 * @brief Blocking probability of the M/M/1/K queue.
 *
 * @param lambda arrival rate
 * @param mu     service rate
 * @param K      system capacity, jobs in service included, K >= 1
 */
template <class T>
Mm1kLossResult<T> qsys_mm1k_loss(const T& lambda, const T& mu, unsigned K) {
    if (K == 0) throw InputError("qsys_mm1k_loss: K must be at least 1");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    Mm1kLossResult<T> r;
    r.utilization = rho;
    // The removable singularity, taken as the reference takes it. The tolerance
    // matches qsys_mm1k_loss.m; on an exact T the subtraction is exact and the
    // test fires only at rho = 1 itself.
    if (std::fabs(num_traits<T>::to_double(rho) - 1.0) < lang::GlobalConstants::FineTol) {
        r.lossProbability = one / num_traits<T>::from_int(static_cast<int>(K) + 1);
        return r;
    }
    r.lossProbability = (one - rho) / (one - num_pow_int(rho, K + 1)) * num_pow_int(rho, K);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MM1K_LOSS_H
