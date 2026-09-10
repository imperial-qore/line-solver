/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MTGS0_MOL_H
#define LINE_API_QSYS_MTGS0_MOL_H

/**
 * Modified-offered-load and pointwise-stationary approximations for a
 * time-varying multiserver system.
 *
 * Templated port of matlab/src/api/qsys/qsys_mtgs0_mol.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mtgs0_mol.java.
 *
 * THE ONE IDEA. A stationary loss system with offered load a blocks with
 * probability B(s,a). In a time-varying system the question is WHICH LOAD goes
 * into that formula. PSA uses the instantaneous one, lambda(t)E[S]. MOL uses the
 * offered load of the corresponding INFINITE-SERVER system,
 *
 *   m(t) = E[S] E[lambda(t - Se)] = int_0^Inf lambda(t-x) P(S>x) dx,
 *
 * which is EXACT there and therefore carries the time lag and the smoothing the
 * finite-server system also has. MOL is then B(s,m(t)). The difference between
 * the two is precisely the lag: PSA peaks when the arrival rate peaks, MOL peaks
 * later, and the real system peaks later too.
 *
 * WHAT TO EXPECT. Against the exact time-varying birth-death chain on a
 * sinusoidal rate, MOL cuts the mean RELATIVE error roughly threefold (0.13
 * against 0.44 at s = 100) because it gets the phase right; it does not always
 * win on ABSOLUTE error, which is dominated by the peak of the cycle. Under
 * constant input MOL is exact.
 *
 * ARITHMETIC. The offered load is a quadrature, so transcendental only. The
 * Erlang recursions themselves are exact and are exposed separately.
 *
 * Reference: W. A. Massey, W. Whitt (1994). An analysis of the modified offered
 * load approximation for the nonstationary Erlang loss model. Annals of Applied
 * Probability 4(4), 1145-1160; W. Whitt (1991). Management Science 37(3),
 * 307-314.
 */

#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "line/api/qsys/qsys_mtginf.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** MOL and PSA measures of a time-varying multiserver system. */
template <class T>
struct QsysMolResult {
    std::vector<T> times;         ///< the evaluation times
    std::vector<T> offeredLoad;   ///< m(t), the infinite-server load
    std::vector<T> instantLoad;   ///< lambda(t)E[S]
    std::vector<T> probBlockMOL;  ///< B(s,m(t)) or C(s,m(t))
    std::vector<T> probBlockPSA;  ///< the same at the instantaneous load
    std::vector<T> meanBusyMOL;   ///< carried load m(t)(1-B), or min(m,s) for the delay model
    std::vector<T> arrivalRate;   ///< lambda(t)
};

/**
 * Erlang B by the recursion B_j = a B_{j-1}/(j + a B_{j-1}), which never forms
 * a^s/s! and so never overflows.
 *
 * @param s number of servers
 * @param a offered load in erlangs
 */
template <class T>
T qsys_erlang_b(unsigned s, const T& a) {
    T b = num_traits<T>::from_int(1);
    for (unsigned j = 1; j <= s; ++j)
        b = a * b / (num_traits<T>::from_int(static_cast<long>(j)) + a * b);
    return b;
}

/**
 * Erlang C from the same recursion; 1 when the load saturates the servers.
 *
 * @param s number of servers
 * @param a offered load in erlangs
 */
template <class T>
T qsys_erlang_c(unsigned s, const T& a) {
    const T one = num_traits<T>::from_int(1);
    const T sT = num_traits<T>::from_int(static_cast<long>(s));
    if (a >= sT) return one;
    const T b = qsys_erlang_b(s, a);
    const T rho = a / sT;
    return b / (one - rho * (one - b));
}

/**
 * @param lambdaFun   the arrival rate
 * @param serviceCcdf G^c(x) = P(S > x)
 * @param ES          the mean service time
 * @param s           number of servers
 * @param tvals       times at which to evaluate
 * @param startTime   time the system started empty; -Inf assumes an infinite past
 * @param delay       use Erlang C rather than Erlang B
 */
template <class T>
QsysMolResult<T> qsys_mtgs0_mol(const std::function<T(const T&)>& lambdaFun,
                                const std::function<T(const T&)>& serviceCcdf, const T& ES,
                                unsigned s, const std::vector<T>& tvals,
                                double startTime = -std::numeric_limits<double>::infinity(),
                                bool delay = false) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mtgs0_mol integrates the offered load, so it needs inexact arithmetic");
    if (s < 1) throw InputError("qsys_mtgs0_mol: the number of servers s must be at least 1");
    const QsysMtginfResult<T> inf =
        qsys_mtginf<T>(lambdaFun, serviceCcdf, ES, tvals, startTime);
    QsysMolResult<T> r;
    r.times = inf.times;
    r.offeredLoad = inf.meanNumber;
    r.instantLoad = inf.offeredLoadPSA;
    r.arrivalRate = inf.arrivalRate;
    const std::size_t n = inf.times.size();
    r.probBlockMOL.resize(n);
    r.probBlockPSA.resize(n);
    r.meanBusyMOL.resize(n);
    const T sT = num_traits<T>::from_int(static_cast<long>(s));
    for (std::size_t i = 0; i < n; ++i) {
        r.probBlockMOL[i] = delay ? qsys_erlang_c(s, r.offeredLoad[i])
                                  : qsys_erlang_b(s, r.offeredLoad[i]);
        r.probBlockPSA[i] = delay ? qsys_erlang_c(s, r.instantLoad[i])
                                  : qsys_erlang_b(s, r.instantLoad[i]);
        r.meanBusyMOL[i] = delay ? detail::num_min(r.offeredLoad[i], sT)
                                 : T(r.offeredLoad[i] * (num_traits<T>::from_int(1) - r.probBlockMOL[i]));
    }
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MTGS0_MOL_H
