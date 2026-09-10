/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_POLLING_POLLING_QSYS_EXHAUSTIVE_H
#define LINE_API_POLLING_POLLING_QSYS_EXHAUSTIVE_H

/**
 * Mean waiting times of a polling system under exhaustive service.
 *
 * Templated port of matlab/src/api/polling/polling_qsys_exhaustive.m, which is
 * the station-time method of Ferguson and Aminetzah (1985) as reported by
 * Takagi, ACM Computing Surveys 20(1), 1988, eq. (15). The MATLAB version takes
 * MAP descriptors and immediately reduces them to the first two moments of the
 * arrival, service and switchover processes; this port takes those moments
 * directly through line::polling::PollingMoments, exactly as the 1-limited and
 * decrementing ports do, so the MAP reduction stays in line::mam.
 *
 * The method solves an n^2 x n^2 linear system for the station times r_ij and
 * then reads the waiting times off them. Everything is a rational function of
 * the input moments plus one exact linear solve, so the whole computation stays
 * in the field: a polling system with rational parameters has an exactly
 * representable mean waiting time and the port is instantiable at Rational.
 * That is worth having here because the denominators 1 - rho and 1 - rho_i both
 * vanish at a stability boundary and appear cubed, so a rounded evaluation near
 * one can return a finite but meaningless number.
 *
 * ORACLES USED IN THE TESTS.
 *  - n = 1 collapse. The formula reduces to the M/G/1 queue with multiple
 *    vacations, W = lambda b2/(2(1-rho)) + E[R^2]/(2 E[R]), which the port
 *    reproduces exactly (in Rational, digit for digit).
 *  - Pseudo-conservation law of Boxma and Groenendijk (1987),
 *      sum_i rho_i W_i = rho/(2(1-rho)) sum_i lambda_i b2_i
 *                      + rho (delta2tot + R^2)/(2R)
 *                      + R (rho^2 - sum_i rho_i^2)/(2(1-rho))
 *                      + sum_i E[M_i],
 *    with E[M_i] = 0 under exhaustive service and E[M_i] = rho_i^2 R/(1-rho)
 *    under gated service. The port satisfies it as an identity, so exhaustive
 *    and gated are cross-checked against one common invariant.
 *
 * A simulator is deliberately NOT used as an oracle: LINE's LDES polling server
 * parks at the last visited queue rather than roving, which puts a systematic
 * negative offset of a few percent on its waiting times relative to Takagi.
 */

#include <cstddef>
#include <vector>

#include "line/api/polling/polling_qsys_1limited.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace polling {

namespace detail {

/** Derived aggregates shared by the exhaustive and gated formulas. */
template <class T>
struct PollingAggregates {
    std::vector<T> rho1;  ///< per-queue load lambda_i b_i
    T rho;                ///< total load
    T R;                  ///< total switchover time
};

template <class T>
PollingAggregates<T> polling_aggregates(const PollingMoments<T>& m, const char* who) {
    const std::size_t n = m.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    PollingAggregates<T> a;
    a.rho1.resize(n);
    a.rho = zero;
    a.R = zero;
    for (std::size_t i = 0; i < n; ++i) {
        a.rho1[i] = m.lambda[i] * m.b[i];
        if (a.rho1[i] <= zero)
            throw InputError(std::string(who) + ": every queue must carry a positive load");
        a.rho += a.rho1[i];
        a.R += m.r[i];
    }
    if (a.rho >= one) throw NumericError(std::string(who) + ": unstable system, rho >= 1");
    if (a.R <= zero) throw InputError(std::string(who) + ": zero total switchover time");
    return a;
}

}  // namespace detail

/**
 * Exhaustive service: the server empties a queue completely before switching.
 *
 * Takagi (1988) eq. (15) in the station-time form. The unknowns are the n^2
 * station times r_ij collected row-major as index (i-1) n + j, 1-based in the
 * reference and shifted by one here.
 *
 * @param m per-queue first two moments of arrivals, service and switchover
 * @return  the n mean waiting times
 */
template <class T>
std::vector<T> polling_qsys_exhaustive(const PollingMoments<T>& m) {
    m.validate("polling_qsys_exhaustive");
    const std::size_t n = m.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const detail::PollingAggregates<T> a =
        detail::polling_aggregates(m, "polling_qsys_exhaustive");
    const std::vector<T>& rho1 = a.rho1;
    const T rho = a.rho, r = a.R;

    const std::size_t nn = n * n;
    Matrix<T> A(nn, nn, zero);
    std::vector<T> rhs(nn, zero);
    std::size_t row = 0;
    for (std::size_t i = 1; i <= n; ++i) {
        for (std::size_t j = 1; j <= n; ++j, ++row) {
            if (i > j) {
                for (std::size_t k = i + 1; k <= n; ++k) A(row, (j - 1) * n + k - 1) -= one;
                for (std::size_t k = 1; k + 1 <= j; ++k) A(row, (j - 1) * n + k - 1) -= one;
                for (std::size_t k = j; k + 1 <= i; ++k) A(row, (k - 1) * n + j - 1) -= one;
                A(row, (i - 1) * n + j - 1) += (one - rho1[i - 1]) / rho1[i - 1];
            } else if (j > i) {
                for (std::size_t k = i + 1; k + 1 <= j; ++k) A(row, (j - 1) * n + k - 1) -= one;
                for (std::size_t k = j; k <= n; ++k) A(row, (k - 1) * n + j - 1) -= one;
                for (std::size_t k = 1; k + 1 <= i; ++k) A(row, (k - 1) * n + j - 1) -= one;
                A(row, (i - 1) * n + j - 1) += (one - rho1[i - 1]) / rho1[i - 1];
            } else {
                A(row, (i - 1) * n + i - 1) += one;
                for (std::size_t k = 1; k <= n; ++k)
                    if (k != i) A(row, (i - 1) * n + k - 1) -= rho1[i - 1] / (one - rho1[i - 1]);
                // The switchover variance charged to queue i is the one of the
                // switchover that PRECEDES it, index i-1 with wraparound to n.
                const T dprev = (i > 1) ? m.delta2[i - 2] : m.delta2[n - 1];
                const T omr = one - rho1[i - 1];
                rhs[row] = dprev / (omr * omr) +
                           m.lambda[i - 1] * m.b2[i - 1] * r * omr / ((one - rho) * omr * omr * omr);
            }
        }
    }

    const std::vector<T> f = solve(A, rhs);

    std::vector<T> W(n);
    for (std::size_t i = 1; i <= n; ++i) {
        const T omr = one - rho1[i - 1];
        T w = m.lambda[i - 1] * m.b2[i - 1] / (two * omr);
        w += r * omr / (two * (one - rho));
        T s = zero;
        for (std::size_t j = 1; j <= n; ++j)
            if (j != i) s += f[(i - 1) * n + j - 1];
        s *= omr / rho1[i - 1];
        s += (i > 1) ? m.delta2[i - 2] : m.delta2[n - 1];
        s /= r * omr * two / (one - rho);
        W[i - 1] = w + s;
    }
    return W;
}

/**
 * Gated service: only the jobs found at the polling instant are served.
 *
 * Takagi (1988) eq. (20), same station-time unknowns and the same layout. The
 * switchover variance enters at index i here, not i-1 as in the exhaustive
 * form; that asymmetry is in the reference and is reproduced, and both forms
 * satisfy the pseudo-conservation law, which is what pins them.
 */
template <class T>
std::vector<T> polling_qsys_gated(const PollingMoments<T>& m) {
    m.validate("polling_qsys_gated");
    const std::size_t n = m.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const detail::PollingAggregates<T> a = detail::polling_aggregates(m, "polling_qsys_gated");
    const std::vector<T>& rho1 = a.rho1;
    const T rho = a.rho, r = a.R;

    const std::size_t nn = n * n;
    Matrix<T> A(nn, nn, zero);
    std::vector<T> rhs(nn, zero);
    std::size_t row = 0;
    for (std::size_t i = 1; i <= n; ++i) {
        for (std::size_t j = 1; j <= n; ++j, ++row) {
            if (i > j) {
                for (std::size_t k = i; k <= n; ++k) A(row, (j - 1) * n + k - 1) = -one;
                for (std::size_t k = 1; k + 1 <= j; ++k) A(row, (j - 1) * n + k - 1) = -one;
                for (std::size_t k = j; k + 1 <= i; ++k) A(row, (k - 1) * n + j - 1) = -one;
                A(row, (i - 1) * n + j - 1) = one / rho1[i - 1];
            } else if (j > i) {
                for (std::size_t k = i; k + 1 <= j; ++k) A(row, (j - 1) * n + k - 1) = -one;
                for (std::size_t k = j; k <= n; ++k) A(row, (k - 1) * n + j - 1) = -one;
                for (std::size_t k = 1; k + 1 <= i; ++k) A(row, (k - 1) * n + j - 1) = -one;
                A(row, (i - 1) * n + j - 1) = one / rho1[i - 1];
            } else {
                A(row, (i - 1) * n + i - 1) += one;
                for (std::size_t k = 1; k <= n; ++k)
                    if (k != i) A(row, (i - 1) * n + k - 1) = -rho1[i - 1];
                for (std::size_t k = 1; k <= n; ++k)
                    A(row, (k - 1) * n + i - 1) -= rho1[i - 1] * rho1[i - 1];
                rhs[row] = m.delta2[i - 1] + m.lambda[i - 1] * m.b2[i - 1] * r / (one - rho);
            }
        }
    }

    const std::vector<T> f = solve(A, rhs);

    std::vector<T> W(n);
    for (std::size_t i = 1; i <= n; ++i) {
        T w = (one + rho1[i - 1]) * r / (two * (one - rho));
        T s = zero;
        for (std::size_t j = 1; j <= n; ++j)
            if (j != i) s += f[(i - 1) * n + j - 1];
        s /= rho1[i - 1];
        for (std::size_t j = 1; j <= n; ++j) s += f[(j - 1) * n + i - 1];
        w += (one - rho) * (one + rho1[i - 1]) * s / (two * r);
        W[i - 1] = w;
    }
    return W;
}

}  // namespace polling
}  // namespace line

#endif  // LINE_API_POLLING_POLLING_QSYS_EXHAUSTIVE_H
