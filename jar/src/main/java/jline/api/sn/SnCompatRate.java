/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import jline.util.matrix.Matrix;

/**
 * Total service rate of a station served by heterogeneous server pools with a
 * class-compatibility graph, and the peak that normalizes its utilization.
 *
 * <p>Port of {@code matlab/src/api/sn/sn_compat_rate.m}. A pool t holds
 * {@code counts[t]} identical servers, each running at {@code rates[t]}, and
 * may serve operand j when {@code compat(t, j)} is nonzero. The rate the
 * station clears in state n is
 *
 * <pre>
 *   mu(n) = sum_t counts[t]*rates[t]*min(1, sum_{j: compat(t,j) != 0} n(j))
 * </pre>
 *
 * <p>the ACTIVATED-SERVER law: a pool contributes its full rate as soon as it
 * is compatible with at least one operand PRESENT. This is the
 * order-independent reading of a compatibility structure -- at an INTEGER state
 * mu depends on n only through its SUPPORT, so it is invariant to the arrival
 * order and to any permutation of the microstate, which is exactly the condition
 * an OI station has to meet (Dorsman and Gardner, Queueing Systems 107:205-256,
 * 2024, Fig. 1). It is also what {@code pas_compatibility_5class.m} encodes for
 * a flat Network, so the layered and flat readings of one compatibility matrix
 * agree.
 *
 * <p>WHY min(1, .) AND NOT AN INDICATOR. At every integer state the two agree
 * exactly -- a pool with at least one compatible job present is fully active,
 * one with none is idle -- so nothing about the OI law on the real state lattice
 * changes. They part company only at a FRACTIONAL argument, which is what a
 * mean-value solver hands this function: AMVA evaluates the rate at a MEAN
 * population, and under a hard indicator any operand with a mean above zero,
 * however small, activates every pool it touches. A compatibility structure
 * would then be invisible to AMVA whenever every operand is a little bit busy --
 * which is nearly always. Scaling linearly below one job keeps the structure
 * visible at the evaluation point while leaving the integer-state law untouched;
 * it is the ordinary continuous relaxation of a step function, and the CTMC and
 * simulation paths, which only ever evaluate at integer states, cannot tell the
 * difference.
 *
 * <p>IT IS NOT A MATCHING. A pool of two servers compatible with a class
 * holding ONE job contributes both servers here, which over-counts against a
 * non-redundant system where one server serves one job. That is deliberate: the
 * matching size depends on the counts and not only on the support, so it is NOT
 * order independent and would take the station outside the product form the OI
 * closure is built on. A model that means the matching wants a different
 * station, not a different reading of this one.
 */
public class SnCompatRate {

    /**
     * Rate cleared by the pools when the operands in {@code n} are present.
     *
     * @param compat (npools x noperands), nonzero where the pool may serve
     * @param counts (npools) servers held by each pool
     * @param rates  (npools) per-server rate of each pool
     * @param n      (noperands) per-operand population, integer or fractional
     * @return the total service rate mu(n)
     */
    public static double snCompatRate(Matrix compat, double[] counts, double[] rates, double[] n) {
        int npools = counts.length;
        if (rates.length != npools) {
            throw new IllegalArgumentException("sn_compat_rate: one rate per pool is required");
        }
        if (compat.getNumRows() != npools) {
            throw new IllegalArgumentException("sn_compat_rate: compat must have one row per pool");
        }
        if (n.length != compat.getNumCols()) {
            throw new IllegalArgumentException("sn_compat_rate: n must have one entry per operand");
        }
        double mu = 0.0;
        for (int t = 0; t < npools; t++) {
            // The pool is activated ONCE by the jobs it can reach, not once per
            // operand: its weight is the compatible load, capped at one job.
            double load = 0.0;
            for (int j = 0; j < n.length; j++) {
                if (compat.get(t, j) != 0 && n[j] > 0) {
                    load += n[j];
                }
            }
            if (load > 1.0) {
                load = 1.0;
            }
            mu += counts[t] * rates[t] * load;
        }
        return mu;
    }

    /**
     * Rate with every pool active, {@code sum_t counts[t]*rates[t]}.
     *
     * <p>Utilization at a rate-scaled station is reported as U = T*S/peak, and
     * the peak is a property of the DECLARATION rather than of a state, so it is
     * computed once and handed to the solver beside the rate handle rather than
     * recovered from {@link #snCompatRate} at a guessed state.
     */
    /**
     * Rate scaling eta(n) a compatibility declaration imposes on its station.
     *
     * <p>This is what SolverLN carries onto the layer station, and it is NOT
     * {@code snCompatRate / snCompatPeak}. The denominator is the rate the SAME
     * population would obtain under FULL compatibility,
     *
     * <pre>
     *   eta(n) = mu(n) / (peak * min(1, sum_j n(j)))
     * </pre>
     *
     * <p>so eta isolates the effect of the compatibility GRAPH and nothing else.
     * The denominator DAMPS BY OCCUPANCY RELATIVE TO THE SERVER COUNT,
     * {@code min(1, N/S)}, because that is precisely what the solver's own
     * multiserver term contributes: it applies {@code min(N, S)} servers at the
     * average server rate {@code peak/S}, so
     *
     * <pre>
     *   min(N,S) * (peak/S) * eta(n) = mu(n)
     * </pre>
     *
     * and the station clears the activated-server rate exactly, at every state.
     *
     * <p>DAMPING BY {@code min(1, N)} INSTEAD -- which this did until
     * 2026-08-28 -- leaves the effective law at {@code min(N,S)/S * mu(n)},
     * which cancels the REDUNDANCY SPEED-UP the activated-server law exists to
     * express: a pool of S servers facing one compatible job clears S, not 1,
     * because every one of them works on it and the first to finish cancels the
     * rest. Under the old normalisation a fully-compatible pool reduced to the
     * plain multiserver, so the OI machinery did no work in the homogeneous
     * case and LDES, which simulates mu(n) directly, disagreed with it by that
     * factor -- measured 1.36858 against 1.15976 on {@code lqn_server_pools}.
     *
     * <p>eta is therefore ABOVE ONE at low occupancy, which is not a defect: it
     * is the speed-up carried by servers that would otherwise be idle.
     */
    public static double snCompatScaling(Matrix compat, double[] counts, double[] rates,
                                         double[] n) {
        double total = 0.0;
        for (int j = 0; j < n.length; j++) {
            if (n[j] > 0) {
                total += n[j];
            }
        }
        if (!(total > 0)) {
            return 1.0; // an empty station: nothing to scale
        }
        double nservers = 0.0;
        for (int t = 0; t < counts.length; t++) {
            nservers += counts[t];
        }
        if (!(nservers > 0)) {
            return 1.0;
        }
        double ref = snCompatPeak(counts, rates) * Math.min(1.0, total / nservers);
        return snCompatRate(compat, counts, rates, n) / ref;
    }

    public static double snCompatPeak(double[] counts, double[] rates) {
        if (rates.length != counts.length) {
            throw new IllegalArgumentException("sn_compat_peak: one rate per pool is required");
        }
        double peak = 0.0;
        for (int t = 0; t < counts.length; t++) {
            peak += counts[t] * rates[t];
        }
        return peak;
    }
}
