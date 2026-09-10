/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_POLLING_POLLING_QSYS_1LIMITED_H
#define LINE_API_POLLING_POLLING_QSYS_1LIMITED_H

/**
 * Mean waiting times in polling systems: 1-limited and decrementing service.
 *
 * Templated port of matlab/src/api/polling/polling_qsys_1limited.m and
 * polling_qsys_decrementing.m. The MATLAB versions take MAP descriptors and
 * reduce them to the first two moments of the arrival, service and switchover
 * processes; this port takes those moments directly, so the MAP reduction stays
 * in line::mam (map_lambda, map_mean, map_moment, map_var) and the waiting-time
 * formula is a pure rational function of the moments. Callers holding MAPs
 * compose the two.
 *
 * Both formulas stay in the field, so a polling system with rational parameters
 * has an exactly representable mean waiting time. Worth having: the denominator
 * 1 - rho - lambda R vanishes at the stability boundary, and near it a rounded
 * evaluation can return a finite but meaningless number where the exact one
 * shows the pole.
 */

#include <cstddef>
#include <limits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace polling {

/**
 * Per-queue first two moments of the arrival, service and switchover
 * processes, the reduced form both formulas consume.
 */
template <class T>
struct PollingMoments {
    std::vector<T> lambda;  ///< arrival rate per queue
    std::vector<T> b;       ///< mean service time per queue
    std::vector<T> b2;      ///< second raw moment of the service time
    std::vector<T> r;       ///< mean switchover time per queue
    std::vector<T> delta2;  ///< variance of the switchover time per queue

    std::size_t size() const { return lambda.size(); }
    void validate(const char* who) const {
        if (lambda.empty()) throw InputError(std::string(who) + ": empty system");
        if (b.size() != lambda.size() || b2.size() != lambda.size() ||
            r.size() != lambda.size() || delta2.size() != lambda.size())
            throw InputError(std::string(who) + ": moment vectors have different lengths");
    }
};

/**
 * 1-limited polling: one job served per visit.
 *
 * W_i = (1-rho+rho_i)/(1-rho-lambda_i R) * (1-rho)/((1-rho)rho + sum rho_j^2)
 *       * (rho/(2(1-rho)) sum_j lambda_j b2_j + rho sum_j delta2_j/(2R)
 *          + R/(2(1-rho)) rho_i (1+rho_i))
 */
template <class T>
std::vector<T> polling_qsys_1limited(const PollingMoments<T>& m) {
    m.validate("polling_qsys_1limited");
    const std::size_t n = m.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    std::vector<T> rho1(n);
    T rho = zero, R = zero, sumLb2 = zero, sumDelta2 = zero, sumRho2 = zero;
    for (std::size_t i = 0; i < n; ++i) {
        rho1[i] = m.lambda[i] * m.b[i];
        rho += rho1[i];
        R += m.r[i];
        sumLb2 += m.lambda[i] * m.b2[i];
        sumDelta2 += m.delta2[i];
        sumRho2 += rho1[i] * rho1[i];
    }
    if (rho >= one) throw NumericError("polling_qsys_1limited: unstable system, rho >= 1");
    if (R == zero) throw InputError("polling_qsys_1limited: zero total switchover time");

    const T den2 = (one - rho) * rho + sumRho2;
    if (den2 == zero) throw NumericError("polling_qsys_1limited: degenerate load configuration");

    std::vector<T> W(n);
    for (std::size_t i = 0; i < n; ++i) {
        const T den1 = one - rho - m.lambda[i] * R;
        // AT THE POLE THE ANSWER IS INFINITE, NOT UNAVAILABLE. Queue i is stable
        // under 1-limited service while rho + lambda_i R < 1; on the boundary its
        // mean waiting time diverges, and polling_qsys_1limited.m divides through
        // and returns Inf rather than erroring. Refusing here turned a model the
        // reference answers into a failed solve (polling_klimited, whose second
        // class sits exactly on the boundary). An exact field has no infinity, so
        // there the pole is still reported rather than misrepresented as a number.
        if (den1 == zero) {
            if constexpr (num_traits<T>::is_exact) {
                throw NumericError(
                    "polling_qsys_1limited: queue at the stability boundary, "
                    "rho + lambda_i R = 1, and the exact field has no infinity");
            } else {
                W[i] = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
                continue;
            }
        }
        T w = (one - rho + rho1[i]) / den1;
        w *= (one - rho) / den2;
        w *= rho / (two * (one - rho)) * sumLb2 + rho * sumDelta2 / (two * R) +
             R / (two * (one - rho)) * (rho1[i] * (one + rho1[i]));
        W[i] = w;
    }
    return W;
}

/**
 * Decrementing service, symmetric systems only (the MATLAB version rejects
 * asymmetric parameters with a 1e-6 relative tolerance; here the check is
 * exact, which is the right test in an exact field and a stricter one in
 * double).
 *
 * W = delta2/(2r) + (N lambda b2 (1 - lambda r) + (r + lambda delta2)(N - rho))
 *                   / (2 (1 - rho - lambda r (N - rho)))
 */
template <class T>
std::vector<T> polling_qsys_decrementing(const PollingMoments<T>& m) {
    m.validate("polling_qsys_decrementing");
    const std::size_t n = m.size();
    for (std::size_t i = 1; i < n; ++i)
        if (m.lambda[i] != m.lambda[0] || m.b[i] != m.b[0] || m.b2[i] != m.b2[0] ||
            m.r[i] != m.r[0] || m.delta2[i] != m.delta2[0])
            throw InputError(
                "polling_qsys_decrementing: only symmetric systems are supported (identical "
                "arrival, service and switchover parameters across all queues)");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T N = num_traits<T>::from_int(static_cast<long>(n));
    const T lam = m.lambda[0], b2s = m.b2[0], r = m.r[0], d2 = m.delta2[0];
    const T rho = N * lam * m.b[0];

    const T denom = two * (one - rho - lam * r * (N - rho));
    if (denom <= zero)
        throw NumericError(
            "polling_qsys_decrementing: unstable system, rho + lambda r (N - rho) >= 1");

    const T residualSwitchover = (r > zero) ? d2 / (two * r) : zero;
    const T w = residualSwitchover +
                (N * lam * b2s * (one - lam * r) + (r + lam * d2) * (N - rho)) / denom;
    return std::vector<T>(n, w);
}

}  // namespace polling
}  // namespace line

#endif  // LINE_API_POLLING_POLLING_QSYS_1LIMITED_H
