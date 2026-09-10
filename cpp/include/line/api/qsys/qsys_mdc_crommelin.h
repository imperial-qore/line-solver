/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MDC_CROMMELIN_H
#define LINE_API_QSYS_QSYS_MDC_CROMMELIN_H

/**
 * M/D/c by Crommelin's embedded chain.
 *
 * Templated port of jar/src/main/java/jline/api/qsys/Qsys_mdc_crommelin.java and
 * the native Python `qsys_mdc_crommelin`; MATLAB carries the method under
 * `qsys_dmc.m`'s "see also" rather than as its own file.
 *
 * The chain is embedded at multiples of the deterministic service time s:
 *   X_{n+1} = max(0, X_n - c) + A_n,   A_n ~ Poisson(lambda s),
 * because in one service period exactly min(X_n, c) jobs complete and the
 * arrivals in that period are Poisson. The embedded epochs are Poisson arrival
 * epochs, so PASTA makes the embedded law the time-average law, and the result
 * is EXACT for M/D/c under FCFS, not an approximation.
 *
 * THE TRUNCATION IS THE ONLY ERROR. The default level is
 * max(200, min(2500, 10/(1-rho) + 200)), which the references settled on
 * empirically at six digits for moderate c; the cap keeps the dense LU
 * affordable, since the transition matrix is triangular-banded but not sparse.
 * A caller comparing against another codebase must pass the SAME truncation to
 * both, as the defaults are the only free parameter.
 *
 * The Poisson weights are formed in logs and exponentiated once, which is what
 * keeps lambda s in the hundreds from overflowing the factorial.
 *
 * ARITHMETIC: transcendental, for the Poisson weights.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

template <class T>
struct MDcCrommelinResult {
    T meanQueueLength;   ///< E[N], jobs in system
    T meanWaitingQueue;  ///< Lq = E[(N-c)+]
    T meanWaitingTime;   ///< Wq = Lq/lambda
    T meanSojournTime;   ///< W = Wq + s
    T utilization;       ///< rho = lambda s / c
};

/**
 * @param lambda_arr Poisson arrival rate
 * @param s          deterministic service time
 * @param c          number of servers
 * @param truncation state-space cap; <= 0 selects the reference's automatic level
 */
template <class T>
MDcCrommelinResult<T> qsys_mdc_crommelin(const T& lambda_arr, const T& s, unsigned c,
                                         long truncation = -1) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mdc_crommelin forms Poisson weights in logs");
    const double lam = num_traits<T>::to_double(lambda_arr);
    const double sv = num_traits<T>::to_double(s);
    if (!(lam > 0.0)) throw InputError("qsys_mdc_crommelin: the arrival rate must be positive");
    if (!(sv > 0.0)) throw InputError("qsys_mdc_crommelin: the service time must be positive");
    if (c < 1) throw InputError("qsys_mdc_crommelin: the number of servers must be at least one");

    const double a = lam * sv;
    const double rho = a / static_cast<double>(c);
    if (!(rho < 1.0 - 1e-12))
        throw InputError("qsys_mdc_crommelin: the load must be strictly below one");

    const long autoN = std::max(
        200L, std::min(2500L, static_cast<long>(10.0 / (1.0 - rho)) + 200L));
    const std::size_t nMax = static_cast<std::size_t>(truncation > 0 ? truncation : autoN);
    const std::size_t n = nMax + 1;

    // Poisson(a) weights in logs, so a in the hundreds does not overflow.
    std::vector<double> pmf(n);
    const double logA = std::log(a);
    double logFact = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        pmf[k] = std::exp(-a + static_cast<double>(k) * logA - logFact);
        logFact += std::log(static_cast<double>(k + 1));
    }

    // (P' - I) with the last row replaced by the normalization sum(pi) = 1.
    // P(i,j) = pmf[j - max(0, i - c)] for j >= max(0, i - c), zero otherwise.
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> A(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t base = (i <= c) ? 0 : i - c;
        for (std::size_t j = base; j < n; ++j)
            A(j, i) = num_traits<T>::from_double(pmf[j - base]);
    }
    for (std::size_t i = 0; i < n; ++i) A(i, i) = T(A(i, i) - one);
    for (std::size_t j = 0; j < n; ++j) A(n - 1, j) = one;

    std::vector<T> rhs(n, zero);
    rhs[n - 1] = one;
    const std::vector<T> pi = solve(A, rhs);

    T meanN = zero, Lq = zero;
    for (std::size_t i = 0; i < n; ++i) {
        meanN += num_traits<T>::from_int(static_cast<long>(i)) * pi[i];
        if (i > c) Lq += num_traits<T>::from_int(static_cast<long>(i - c)) * pi[i];
    }

    MDcCrommelinResult<T> r;
    r.meanQueueLength = meanN;
    r.meanWaitingQueue = Lq;
    r.meanWaitingTime = T(Lq / lambda_arr);
    r.meanSojournTime = T(r.meanWaitingTime + s);
    r.utilization = num_traits<T>::from_double(rho);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MDC_CROMMELIN_H
