/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_FEEDBACK_ELIM_H
#define LINE_API_NPFQN_FEEDBACK_ELIM_H

/**
 * Near-immediate feedback elimination for the robust queueing network analyzer.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_feedback_elim.m, cross-checked
 * against jar/src/main/java/jline/api/npfqn/Npfqn_feedback_elim.java.
 *
 * WHY FEEDBACK BREAKS DECOMPOSITION. A parametric decomposition treats the
 * arrival stream at each station as renewal. Feedback destroys that badly: a
 * customer that leaves a busy station and comes straight back arrives exactly
 * when the station is busy, so the flow is strongly correlated with the queue it
 * feeds. The fix is not to model the correlation but to REMOVE the feedback, by
 * folding the repeated visits into the service time:
 *
 *   effective mean service  E[S]/(1-p)
 *   effective service SCV   p + (1-p)cs^2                            (37)
 *   fresh arrival rate      lambda(1-p)
 *   per-visit waiting time  (1-p) times the wait in the modified system
 *
 * The modified system has the SAME heavy-traffic limits for queue length,
 * workload, waiting time and external departures, so this is asymptotically
 * exact rather than merely plausible.
 *
 * NEAR-IMMEDIATE, NOT JUST IMMEDIATE. What matters is whether the customer
 * returns WITHOUT PASSING A BUSIER STATION: a detour through a station of lower
 * traffic intensity is fast on the time scale of the busy station. The
 * probability computed here is therefore the probability of returning to station
 * i through stations of strictly smaller rho only, obtained from the absorbing
 * chain restricted to those stations.
 *
 * ARITHMETIC. Only a linear solve, so this instantiates at T = Rational too.
 *
 * Reference: W. Whitt, W. You (2022). A robust queueing network analyzer based
 * on indices of dispersion. Naval Research Logistics 69(1), 36-56, Section 4.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

/** Outcome of the feedback elimination. */
template <class T>
struct FeedbackElimResult {
    std::vector<T> feedbackProb;    ///< p-hat per station
    std::vector<T> visitInflation;  ///< 1/(1-p), the mean visits per customer
    std::vector<T> modifiedScv;     ///< p + (1-p)cs^2, empty when no SCV was given
    std::vector<T> modifiedRates;   ///< lambda(1-p), empty when no rate was given
    Matrix<T> modifiedRouting;      ///< the immediate-feedback reduction of P
    bool reductionExact = false;    ///< whether that reduction describes this network
};

/**
 * @param P             routing matrix, substochastic
 * @param rho           traffic intensity of each station
 * @param cs2           service SCV of each station, empty to skip modifiedScv
 * @param lambda        arrival rate of each station, empty to skip modifiedRates
 * @param immediateOnly keep only the self-loops, i.e. Section 4.1 feedback
 */
template <class T>
FeedbackElimResult<T> npfqn_feedback_elim(const Matrix<T>& P, const std::vector<T>& rho,
                                          const std::vector<T>& cs2 = std::vector<T>(),
                                          const std::vector<T>& lambda = std::vector<T>(),
                                          bool immediateOnly = false) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t m = P.rows();
    if (P.cols() != m) throw InputError("npfqn_feedback_elim: the routing matrix must be square");
    for (std::size_t i = 0; i < m; ++i) {
        T row = zero;
        for (std::size_t j = 0; j < m; ++j) {
            if (P(i, j) < -num_traits<T>::from_double(1e-12))
                throw InputError("npfqn_feedback_elim: the routing matrix must be non-negative");
            row += P(i, j);
        }
        if (row > one + num_traits<T>::from_double(1e-9))
            throw InputError("npfqn_feedback_elim: the routing matrix must be substochastic");
    }
    if (rho.size() != m)
        throw InputError("npfqn_feedback_elim: one traffic intensity per station is required");

    FeedbackElimResult<T> res;
    res.feedbackProb.assign(m, zero);
    res.visitInflation.assign(m, zero);
    for (std::size_t i = 0; i < m; ++i) {
        T ret = P(i, i);
        if (!immediateOnly) {
            // Stations a customer may pass through on a near-immediate return:
            // those NOT MORE loaded than i. A detour through a busier station is
            // not fast on the time scale of station i, so it is not
            // near-immediate; one through a station of equal load is, which is
            // why the test is <= and not <. This is the cloud of eqs.
            // (3.8)-(3.9) with H = {i}, and the same one solver_rqna applies --
            // the two must not drift.
            const T slack = num_traits<T>::from_double(1e-9);
            std::vector<std::size_t> idx;
            for (std::size_t j = 0; j < m; ++j)
                if (j != i && rho[j] <= T(rho[i] + slack)) idx.push_back(j);
            if (!idx.empty()) {
                const std::size_t k = idx.size();
                Matrix<T> A(k, k);
                std::vector<T> b(k, zero);
                for (std::size_t a = 0; a < k; ++a) {
                    for (std::size_t c = 0; c < k; ++c)
                        A(a, c) = (a == c ? one : zero) - P(idx[a], idx[c]);
                    b[a] = P(idx[a], i);
                }
                // (I-Q)^-1 r: the probability of eventually reaching i from each
                // allowed station without leaving the allowed set.
                const std::vector<T> reach = solve(A, b);
                for (std::size_t a = 0; a < k; ++a) ret += P(i, idx[a]) * reach[a];
            }
        }
        if (ret < zero) ret = zero;
        const T cap = one - num_traits<T>::from_double(1e-12);
        if (ret > cap) ret = cap;
        res.feedbackProb[i] = ret;
        res.visitInflation[i] = one / (one - ret);
    }

    if (!cs2.empty()) {
        if (cs2.size() != m)
            throw InputError("npfqn_feedback_elim: one service SCV per station is required");
        res.modifiedScv.resize(m);
        for (std::size_t i = 0; i < m; ++i)                      // eq. (37)
            res.modifiedScv[i] = res.feedbackProb[i] + (one - res.feedbackProb[i]) * cs2[i];
    }
    if (!lambda.empty()) {
        if (lambda.size() != m)
            throw InputError("npfqn_feedback_elim: one arrival rate per station is required");
        res.modifiedRates.resize(m);
        for (std::size_t i = 0; i < m; ++i)
            res.modifiedRates[i] = lambda[i] * (one - res.feedbackProb[i]);
    }

    // The immediate-feedback reduction: drop the self-loop and renormalize the
    // rest of the row. For near-immediate feedback the return path runs through
    // other stations, so no row-local reduction exists and the elimination
    // applies to the service description instead.
    res.modifiedRouting = P;
    bool exact = true;
    for (std::size_t i = 0; i < m; ++i) {
        if (num_abs(T(res.feedbackProb[i] - P(i, i))) > num_traits<T>::from_double(1e-12))
            exact = false;
        const T loop = P(i, i);
        if (loop <= zero) continue;
        res.modifiedRouting(i, i) = zero;
        T rest = zero, total = zero;
        for (std::size_t j = 0; j < m; ++j) {
            rest += res.modifiedRouting(i, j);
            total += P(i, j);
        }
        if (rest > zero)
            for (std::size_t j = 0; j < m; ++j)
                res.modifiedRouting(i, j) = res.modifiedRouting(i, j) * (total - loop) / rest;
    }
    res.reductionExact = immediateOnly || exact;
    return res;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_FEEDBACK_ELIM_H
