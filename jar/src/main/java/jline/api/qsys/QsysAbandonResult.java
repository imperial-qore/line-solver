/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

/**
 * Steady-state measures of a multiserver queue with customer abandonment, as
 * produced by {@link Qsys_mgisrgi_whitt} and {@link Qsys_erlanga}.
 *
 * <p>Port of the struct returned by MATLAB qsys_mgisrgi_whitt.m.
 *
 * @since LINE 3.1.0
 */
public class QsysAbandonResult {
    /** P(N = k) for k = 0..s+r, N the number in system. */
    public final double[] queueLengthDist;
    /** P(an arrival is blocked) = p_{s+r}; zero when the waiting room is infinite. */
    public final double probLoss;
    /** P(W = 0) among entering customers. */
    public final double probNoWait;
    /** P(S), an entering customer is eventually served. */
    public final double probServed;
    /** P(A) = 1 - P(S), an entering customer eventually abandons. */
    public final double probAbandon;
    /** E[N], the mean number in system. */
    public final double meanNumber;
    /** Var[N]. */
    public final double varNumber;
    /** E[Q] with Q = (N-s)^+, the mean number waiting. */
    public final double meanQueueLength;
    /** Var[Q]. */
    public final double varQueueLength;
    /** E[min(N,s)]/s, the fraction of servers busy. */
    public final double utilization;
    /** Rate of served customers, lambda(1-P_loss)P(S). */
    public final double throughput;
    /** Rate of abandoning customers, lambda(1-P_loss)P(A). */
    public final double abandonRate;
    /** E[W|S], the mean wait of a customer who is served. */
    public final double meanWaitServed;
    /** Var[W|S]. */
    public final double varWaitServed;
    /** E[W|A], the mean time in queue of a customer who abandons. */
    public final double meanWaitAbandon;
    /** Var[W|A]. */
    public final double varWaitAbandon;
    /** E[W] over entering customers, the zero waits included. */
    public final double meanWait;
    /** E[W^2] over entering customers. */
    public final double secondMomentWait;
    /** delta_j, the abandonment rate of the customer jth from the end of the queue. */
    public final double[] abandonRates;
    /** Delta_k, the total abandonment rate with k waiting; index k, so Delta_0 = 0. */
    public final double[] totalAbandonRates;
    /** Waiting spaces actually used, which is the truncation level when r is infinite. */
    public final int numWaitingSpaces;
    /** Whether the patience law was exponential, in which case the answer is exact. */
    public final boolean exponentialPatience;
    /** The exponential patience rate, NaN for a general patience law. */
    public final double patienceRate;
    /** The times at which the waiting-time cdfs were evaluated, null when not requested. */
    public final double[] waitPoints;
    /** P(W &lt;= t | S) at those times, null when not requested. */
    public final double[] cdfWaitServed;
    /** P(W &lt;= t | A) at those times, null when not requested. */
    public final double[] cdfWaitAbandon;
    /** P(W &lt;= t) over entering customers, null when not requested. */
    public final double[] cdfWait;

    public QsysAbandonResult(double[] queueLengthDist, double probLoss, double probNoWait,
                             double probServed, double probAbandon, double meanNumber,
                             double varNumber, double meanQueueLength, double varQueueLength,
                             double utilization, double throughput, double abandonRate,
                             double meanWaitServed, double varWaitServed, double meanWaitAbandon,
                             double varWaitAbandon, double meanWait, double secondMomentWait,
                             double[] abandonRates, double[] totalAbandonRates,
                             int numWaitingSpaces, boolean exponentialPatience,
                             double patienceRate, double[] waitPoints, double[] cdfWaitServed,
                             double[] cdfWaitAbandon, double[] cdfWait) {
        this.queueLengthDist = queueLengthDist;
        this.probLoss = probLoss;
        this.probNoWait = probNoWait;
        this.probServed = probServed;
        this.probAbandon = probAbandon;
        this.meanNumber = meanNumber;
        this.varNumber = varNumber;
        this.meanQueueLength = meanQueueLength;
        this.varQueueLength = varQueueLength;
        this.utilization = utilization;
        this.throughput = throughput;
        this.abandonRate = abandonRate;
        this.meanWaitServed = meanWaitServed;
        this.varWaitServed = varWaitServed;
        this.meanWaitAbandon = meanWaitAbandon;
        this.varWaitAbandon = varWaitAbandon;
        this.meanWait = meanWait;
        this.secondMomentWait = secondMomentWait;
        this.abandonRates = abandonRates;
        this.totalAbandonRates = totalAbandonRates;
        this.numWaitingSpaces = numWaitingSpaces;
        this.exponentialPatience = exponentialPatience;
        this.patienceRate = patienceRate;
        this.waitPoints = waitPoints;
        this.cdfWaitServed = cdfWaitServed;
        this.cdfWaitAbandon = cdfWaitAbandon;
        this.cdfWait = cdfWait;
    }

    @Override
    public String toString() {
        return "QsysAbandonResult(probNoWait=" + probNoWait
                + ", probAbandon=" + probAbandon
                + ", meanQueueLength=" + meanQueueLength
                + ", meanWaitServed=" + meanWaitServed
                + ", meanWaitAbandon=" + meanWaitAbandon
                + ", utilization=" + utilization + ")";
    }
}
