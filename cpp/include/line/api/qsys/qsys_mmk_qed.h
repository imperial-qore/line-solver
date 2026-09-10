/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MMK_QED_H
#define LINE_API_QSYS_MMK_QED_H

/**
 * Halfin-Whitt QED approximation for the M/M/s queue, and the square-root
 * staffing rule that inverts it.
 *
 * Templated port of matlab/src/api/qsys/qsys_mmk_qed.m, qsys_mmk_qed_alpha.m and
 * qsys_mmk_qed_staffing.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mmk_qed.java.
 *
 * Let s grow with the offered load a = lambda/mu so that the SERVER SLACK stays
 * of order sqrt(s), i.e. beta = (1-rho) sqrt(s) = (s-a)/sqrt(s) is held fixed.
 * The delay probability then has the non-degenerate limit
 *
 *   alpha(beta) = [ 1 + beta Phi(beta)/phi(beta) ]^(-1)
 *
 * with phi and Phi the standard normal density and cdf. Servers are busy a
 * fraction 1 - beta/sqrt(s) of the time, so efficiency tends to 1, and yet the
 * delay probability tends to a constant strictly between 0 and 1, so quality
 * does not collapse. Staffing inverts that: s = ceil(a + beta sqrt(a)).
 *
 * ARITHMETIC. erfc and exp put this in the transcendental family; the staffing
 * search is a bisection against a tolerance, so there is nothing exact to
 * preserve.
 *
 * NUMERICS. alpha is evaluated as phi/(phi + beta Phi) rather than as the
 * reciprocal of 1 + beta Phi/phi: the two agree, but the quotient overflows once
 * phi underflows (beta beyond about 38), whereas this form degrades to
 * 0/(0+beta) = 0, which is the correct limit. The exact Erlang C used by the
 * refinement goes through the Erlang B recursion for the same reason, a^j/j!
 * being infinite well below the s these rules propose.
 *
 * Reference: S. Halfin, W. Whitt (1981). Heavy-traffic limits for queues with
 * many exponential servers. Operations Research 29(3), 567-588.
 */

#include <cmath>
#include <cstddef>
#include <string>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** QED measures of the M/M/s queue. */
template <class T>
struct QsysQedResult {
    T offeredLoad;       ///< a = lambda/mu, in erlangs
    T trafficIntensity;  ///< rho = a/s
    T beta;              ///< the QED server slack (s-a)/sqrt(s)
    T probDelay;         ///< alpha(beta)
    T meanWaitDelayed;   ///< E[W | W > 0] = 1/(s mu - lambda), exact for M/M/s
    T meanWait;          ///< E[W] = alpha(beta)/(s mu - lambda)
    T meanQueueLength;   ///< E[Q] = lambda E[W]
    T meanNumber;        ///< E[N] = a + E[Q]
    T utilization;       ///< rho
};

/** Outcome of the square-root staffing rule. */
template <class T>
struct QsysQedStaffingResult {
    unsigned numServers;  ///< the recommended s
    T beta;               ///< the slack achieved, (s-a)/sqrt(s)
    T betaTarget;         ///< the slack the target asks for, before rounding s up
    T offeredLoad;        ///< a = lambda/mu
    T probDelay;          ///< the QED delay probability at the recommended s
    T meanWait;           ///< the QED mean wait at the recommended s
    T serviceLevel;       ///< P(W <= deadline), for the service-level criterion
    bool exactUsed;       ///< whether the exact Erlang C refinement was applied
};

/** Which target the staffing rule is asked to meet. */
enum class QedCriterion { Delay, MeanWait, ServiceLevel };

/**
 * The Halfin-Whitt delay-probability function alpha(beta); 1 at beta <= 0.
 *
 * @param beta the QED server-slack parameter
 */
template <class T>
T qsys_mmk_qed_alpha(const T& beta) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mmk_qed_alpha needs erfc and exp");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (beta <= zero) return one;
    const T two = num_traits<T>::from_int(2);
    using std::erfc;
    using std::exp;
    using std::sqrt;
    const T phi = exp(-beta * beta / two) / sqrt(two * num_traits<T>::from_double(M_PI));
    const T Phi = erfc(T(-beta / sqrt(two))) / two;
    return phi / (phi + beta * Phi);
}

/**
 * Erlang C by the Erlang B recursion B_j = a B_{j-1}/(j + a B_{j-1}), which
 * never forms a^j/j! and so never overflows.
 *
 * @param s      number of servers
 * @param lambda arrival rate
 * @param mu     service rate of one server
 */
template <class T>
T qsys_mmk_qed_erlangc(unsigned s, const T& lambda, const T& mu) {
    const T one = num_traits<T>::from_int(1);
    const T a = lambda / mu;
    T b = one;
    for (unsigned j = 1; j <= s; ++j)
        b = a * b / (num_traits<T>::from_int(static_cast<long>(j)) + a * b);
    const T rho = a / num_traits<T>::from_int(static_cast<long>(s));
    return rho >= one ? one : T(b / (one - rho * (one - b)));
}

/**
 * @param lambda arrival rate
 * @param mu     service rate of one server
 * @param s      number of servers, s >= 1
 */
template <class T>
QsysQedResult<T> qsys_mmk_qed(const T& lambda, const T& mu, unsigned s) {
    static_assert(num_traits<T>::has_transcendental, "qsys_mmk_qed needs erfc and exp");
    const T zero = num_traits<T>::from_int(0);
    if (lambda <= zero) throw InputError("qsys_mmk_qed: the arrival rate lambda must be positive");
    if (mu <= zero) throw InputError("qsys_mmk_qed: the service rate mu must be positive");
    if (s < 1) throw InputError("qsys_mmk_qed: the number of servers s must be at least 1");
    using std::sqrt;
    const T sT = num_traits<T>::from_int(static_cast<long>(s));
    QsysQedResult<T> r;
    r.offeredLoad = lambda / mu;
    r.trafficIntensity = r.offeredLoad / sT;
    r.beta = (sT - r.offeredLoad) / sqrt(sT);
    r.utilization = r.trafficIntensity;
    if (r.beta <= zero) {
        // Not a QED model: every arrival is delayed and there is no steady state.
        r.probDelay = num_traits<T>::from_int(1);
        const T inf = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
        r.meanWaitDelayed = r.meanWait = r.meanQueueLength = r.meanNumber = inf;
        return r;
    }
    r.probDelay = qsys_mmk_qed_alpha(r.beta);
    r.meanWaitDelayed = num_traits<T>::from_int(1) / (sT * mu - lambda);
    r.meanWait = r.probDelay * r.meanWaitDelayed;
    r.meanQueueLength = lambda * r.meanWait;
    r.meanNumber = r.offeredLoad + r.meanQueueLength;
    return r;
}

namespace detail {

/** Bisection for a root of an increasing f on (0, hi]; the bracket grows. */
template <class T, class Fn>
T qed_solve(Fn&& f) {
    const T lo0 = num_traits<T>::from_double(1e-9);
    T lo = lo0, hi = num_traits<T>::from_int(1);
    if (f(lo) > num_traits<T>::from_int(0)) return lo;
    while (f(hi) < num_traits<T>::from_int(0)) {
        hi *= num_traits<T>::from_int(2);
        if (hi > num_traits<T>::from_double(1e6))
            throw InputError("qsys_mmk_qed_staffing: no server slack meets the target");
    }
    const T two = num_traits<T>::from_int(2);
    for (int i = 0; i < 200; ++i) {
        const T mid = (lo + hi) / two;
        if (f(mid) < num_traits<T>::from_int(0)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return (lo + hi) / two;
}

/** The exact M/M/s measure against the target. */
template <class T>
bool qed_meets(const T& lambda, const T& mu, unsigned s, const T& target, QedCriterion crit,
               const T& deadline, const T& level) {
    const T sT = num_traits<T>::from_int(static_cast<long>(s));
    if (sT * mu <= lambda) return false;
    const T c = qsys_mmk_qed_erlangc(s, lambda, mu);
    const T wq = c / (sT * mu - lambda);
    switch (crit) {
        case QedCriterion::Delay: return c <= target;
        case QedCriterion::MeanWait: return wq <= target;
        case QedCriterion::ServiceLevel: {
            using std::exp;
            return (num_traits<T>::from_int(1) - c * exp(-(sT * mu - lambda) * deadline)) >= level;
        }
    }
    return false;
}

}  // namespace detail

/**
 * Square-root staffing of the M/M/s queue.
 *
 * @param lambda   arrival rate
 * @param mu       service rate of one server
 * @param target   the largest acceptable P(W>0) for Delay, the largest
 *                 acceptable E[W] for MeanWait, unused for ServiceLevel
 * @param crit     which target to meet
 * @param deadline the deadline of the service-level criterion
 * @param level    the probability that deadline must be met with
 * @param exact    walk s until the EXACT Erlang C measure meets the target
 */
template <class T>
QsysQedStaffingResult<T> qsys_mmk_qed_staffing(
    const T& lambda, const T& mu, const T& target, QedCriterion crit = QedCriterion::Delay,
    const T& deadline = num_traits<T>::from_int(0), const T& level = num_traits<T>::from_int(0),
    bool exact = false) {
    static_assert(num_traits<T>::has_transcendental, "qsys_mmk_qed_staffing needs erfc and exp");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (lambda <= zero)
        throw InputError("qsys_mmk_qed_staffing: the arrival rate lambda must be positive");
    if (mu <= zero) throw InputError("qsys_mmk_qed_staffing: the service rate mu must be positive");
    using std::ceil;
    using std::exp;
    using std::floor;
    using std::sqrt;
    const T a = lambda / mu;
    T betaTarget;
    switch (crit) {
        case QedCriterion::Delay:
            if (!(target > zero && target < one))
                throw InputError(
                    "qsys_mmk_qed_staffing: for the delay criterion the target must be in (0,1)");
            betaTarget = detail::qed_solve<T>([&](const T& b) { return target - qsys_mmk_qed_alpha(b); });
            break;
        case QedCriterion::MeanWait:
            if (!(target > zero))
                throw InputError(
                    "qsys_mmk_qed_staffing: for the meanwait criterion the target must be positive");
            // The residual is written target - E[W] so that it increases in beta.
            betaTarget = detail::qed_solve<T>(
                [&](const T& b) { return target - qsys_mmk_qed_alpha(b) / (mu * b * sqrt(a)); });
            break;
        case QedCriterion::ServiceLevel:
            if (!(level > zero && level < one) || !(deadline > zero))
                throw InputError("qsys_mmk_qed_staffing: the service level must be in (0,1) and "
                                 "the deadline positive");
            betaTarget = detail::qed_solve<T>([&](const T& b) {
                const T sApprox = a + b * sqrt(a);
                return (one - qsys_mmk_qed_alpha(b) * exp(-mu * b * sqrt(sApprox) * deadline)) - level;
            });
            break;
    }

    long sl = static_cast<long>(num_traits<T>::to_double(T(ceil(a + betaTarget * sqrt(a)))));
    if (sl < 1) sl = 1;
    unsigned s = static_cast<unsigned>(sl);
    if (num_traits<T>::from_int(static_cast<long>(s)) * mu <= lambda)
        s = static_cast<unsigned>(num_traits<T>::to_double(T(floor(a)))) + 1;
    if (exact) {
        while (!detail::qed_meets(lambda, mu, s, target, crit, deadline, level)) {
            ++s;
            if (s > 10000000u)
                throw InputError("qsys_mmk_qed_staffing: the exact refinement passed 10^7 servers "
                                 "without meeting the target");
        }
        while (s > 1 && detail::qed_meets(lambda, mu, s - 1, target, crit, deadline, level)) --s;
    }

    const QsysQedResult<T> qed = qsys_mmk_qed(lambda, mu, s);
    QsysQedStaffingResult<T> res;
    res.numServers = s;
    res.beta = qed.beta;
    res.betaTarget = betaTarget;
    res.offeredLoad = a;
    res.probDelay = qed.probDelay;
    res.meanWait = qed.meanWait;
    res.exactUsed = exact;
    res.serviceLevel = zero;
    if (crit == QedCriterion::ServiceLevel)
        res.serviceLevel =
            one - qed.probDelay *
                      exp(-(num_traits<T>::from_int(static_cast<long>(s)) * mu - lambda) * deadline);
    return res;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MMK_QED_H
