/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.dqsys;

/**
 * Steady-state quantities of a state dependent Bernoulli server.
 *
 * <p>Every time is a number of slots and every rate a per-slot probability.</p>
 */
public final class Bernoulli1Result {

    /** Buffer capacity in jobs, {@link Integer#MAX_VALUE} when unbounded. */
    public final int capacity;

    /** Offered per-slot arrival probability b(n), indexed by n = 0..L. */
    public final double[] arrivalProb;

    /** Per-slot service completion probability p(n), indexed by n-1. */
    public final double[] serviceProb;

    /** Time-stationary queue length law of theorem 2.3, {@code pmf[n]}. */
    public final double[] pmf;

    /** Arrival queue length law of theorem 2.11, {@code arrivalPmf[n]}. */
    public final double[] arrivalPmf;

    /** Probability of an empty system, pi(0). */
    public final double emptyProb;

    /** Fraction of slots with the server busy, 1 - pi(0). */
    public final double utilization;

    /** Carried departures per slot. */
    public final double throughput;

    /** Fraction of offered arrivals lost, 0 on an unbounded buffer. */
    public final double lossProb;

    /** Mean number of jobs in the system. */
    public final double meanQueueLength;

    /** Mean number of jobs waiting, i.e. not in service. */
    public final double meanWaitingQueue;

    /** Mean sojourn time in slots, by Little's law. */
    public final double meanSojournTime;

    /** Mean waiting time in slots. */
    public final double meanWaitingTime;

    /** Normalizing constant H of theorem 2.3. */
    public final double normConst;

    public Bernoulli1Result(int capacity, double[] arrivalProb, double[] serviceProb,
                            double[] pmf, double[] arrivalPmf, double emptyProb,
                            double utilization, double throughput, double lossProb,
                            double meanQueueLength, double meanWaitingQueue,
                            double meanSojournTime, double meanWaitingTime,
                            double normConst) {
        this.capacity = capacity;
        this.arrivalProb = arrivalProb;
        this.serviceProb = serviceProb;
        this.pmf = pmf;
        this.arrivalPmf = arrivalPmf;
        this.emptyProb = emptyProb;
        this.utilization = utilization;
        this.throughput = throughput;
        this.lossProb = lossProb;
        this.meanQueueLength = meanQueueLength;
        this.meanWaitingQueue = meanWaitingQueue;
        this.meanSojournTime = meanSojournTime;
        this.meanWaitingTime = meanWaitingTime;
        this.normConst = normConst;
    }
}
