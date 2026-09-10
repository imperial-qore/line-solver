/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_COMPAT_RATE_H
#define LINE_API_SN_SN_COMPAT_RATE_H

/**
 * Total service rate of a station served by heterogeneous server pools with a
 * class-compatibility graph, and the peak that normalizes its utilization.
 *
 * Port of matlab/src/api/sn/sn_compat_rate.m. A pool t holds counts(t) identical
 * servers, each running at rates(t), and may serve operand j when
 * compat(t, j) is nonzero. The rate the station clears in state n is
 *
 *     mu(n) = sum_t counts(t) * rates(t) * min(1, sum_{j: compat(t,j) != 0} n(j))
 *
 * the ACTIVATED-SERVER law: a pool contributes its full rate as soon as it is
 * compatible with at least one operand PRESENT. This is the order-independent
 * reading of a compatibility structure -- at an INTEGER state mu depends on n
 * only through its SUPPORT, so it is invariant to the arrival order and to any
 * permutation of the microstate, which is exactly the condition an OI station
 * has to meet (Dorsman & Gardner, Queueing Systems 107:205-256, 2024, Fig. 1).
 * It is also what pas_compatibility_5class.m encodes for a flat Network, so the
 * layered and flat readings of one compatibility matrix agree.
 *
 * WHY min(1, .) AND NOT AN INDICATOR. At every integer state the two agree
 * exactly -- a pool with at least one compatible job present is fully active,
 * one with none is idle -- so nothing about the OI law on the real state
 * lattice changes. They part company only at a FRACTIONAL argument, which is
 * what a mean-value solver hands this function: AMVA evaluates the rate at a
 * mean population, and under a hard indicator any operand with a mean above
 * zero, however small, activates every pool it touches. A compatibility
 * structure would then be invisible to AMVA whenever every operand is a little
 * bit busy -- which is nearly always. Scaling linearly below one job keeps the
 * structure visible at the evaluation point while leaving the integer-state law
 * untouched; it is the ordinary continuous relaxation of a step function, and
 * the CTMC and simulation paths, which only ever evaluate at integer states,
 * cannot tell the difference.
 *
 * IT IS NOT A MATCHING. A pool of two servers compatible with a class holding
 * ONE job contributes both servers here, which over-counts against a
 * non-redundant system where one server serves one job. That is deliberate:
 * the matching size depends on the counts and not only on the support, so it is
 * NOT order-independent and would take the station outside the product form the
 * OI closure is built on. A model that means the matching wants a different
 * station, not a different reading of this one.
 *
 * WHY THE PEAK IS SEPARATE. Utilization at a rate-scaled station is reported as
 * U = T*S/peak, and the peak is the rate with every pool active, sum_t
 * counts(t)*rates(t). It is a property of the DECLARATION, not of a state, so
 * it is computed once and handed to the solver beside the handle rather than
 * recovered from mu at a guessed state.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

/**
 * Rate cleared by the pools when the operands in `n` are present.
 *
 * @param compat  (npools x noperands), nonzero where the pool may serve
 * @param counts  (npools) servers held by each pool
 * @param rates   (npools) per-server rate of each pool
 * @param n       (noperands) per-operand population, integer or fractional
 */
template <class T>
T sn_compat_rate(const Matrix<T>& compat, const std::vector<double>& counts,
                 const std::vector<T>& rates, const std::vector<T>& n) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t npools = counts.size();
    if (rates.size() != npools)
        throw InputError("sn_compat_rate: one rate per pool is required");
    if (compat.rows() != npools)
        throw InputError("sn_compat_rate: compat must have one row per pool");
    if (n.size() != compat.cols())
        throw InputError("sn_compat_rate: n must have one entry per operand");
    const T one = num_traits<T>::from_int(1);
    T mu = zero;
    for (std::size_t t = 0; t < npools; ++t) {
        // The pool is activated ONCE by the jobs it can reach, not once per
        // operand: its weight is the compatible load, capped at one job.
        T load = zero;
        for (std::size_t j = 0; j < n.size(); ++j)
            if (compat(t, j) != zero && n[j] > zero) load = T(load + n[j]);
        if (load > one) load = one;
        mu = mu + num_traits<T>::from_double(counts[t]) * rates[t] * load;
    }
    return mu;
}

/** Rate with every pool active: sum_t counts(t)*rates(t). Normalizes U = T*S/peak. */
template <class T>
T sn_compat_peak(const std::vector<double>& counts, const std::vector<T>& rates) {
    const T zero = num_traits<T>::from_int(0);
    if (rates.size() != counts.size())
        throw InputError("sn_compat_peak: one rate per pool is required");
    T peak = zero;
    for (std::size_t t = 0; t < counts.size(); ++t)
        peak = peak + num_traits<T>::from_double(counts[t]) * rates[t];
    return peak;
}

/**
 * Rate scaling eta(n) a compatibility declaration imposes on its station.
 *
 * This is what SolverLN carries onto the layer station, and it is NOT
 * `sn_compat_rate / sn_compat_peak`. The denominator is the rate the SAME
 * population would obtain under FULL compatibility,
 *
 *     eta(n) = mu(n) / (peak * min(1, sum_j n(j)))
 *
 * so eta isolates the effect of the compatibility GRAPH and nothing else. The
 * denominator DAMPS BY OCCUPANCY RELATIVE TO THE SERVER COUNT, min(1, N/S),
 * because that is what the solver's own multiserver term contributes: it
 * applies min(N,S) servers at the average server rate peak/S, so
 *
 *     min(N,S) * (peak/S) * eta(n) = mu(n)
 *
 * and the station clears the activated-server rate exactly, at every state.
 * Damping by min(1, N) instead -- which this did until 2026-08-28 -- left the
 * effective law at min(N,S)/S * mu(n) and cancelled the REDUNDANCY SPEED-UP the
 * activated-server law exists to express: a pool of S servers facing one
 * compatible job clears S, not 1, because every one works on it and the first
 * to finish cancels the rest. eta is therefore above one at low occupancy,
 * which is the speed-up of servers that would otherwise be idle.
 */
template <class T>
T sn_compat_scaling(const Matrix<T>& compat, const std::vector<double>& counts,
                    const std::vector<T>& rates, const std::vector<T>& n) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    T total = zero;
    for (std::size_t j = 0; j < n.size(); ++j)
        if (n[j] > zero) total = T(total + n[j]);
    if (!(total > zero)) return one;  // an empty station: nothing to scale
    T nservers = zero;
    for (std::size_t t = 0; t < counts.size(); ++t)
        nservers = T(nservers + num_traits<T>::from_double(counts[t]));
    if (!(nservers > zero)) return one;
    const T occ = T(total / nservers);
    T ref = sn_compat_peak(counts, rates) * (occ > one ? one : occ);
    return T(sn_compat_rate(compat, counts, rates, n) / ref);
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_COMPAT_RATE_H
