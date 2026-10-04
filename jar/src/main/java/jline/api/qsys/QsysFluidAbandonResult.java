/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

/**
 * Steady state of the G/GI/s+GI fluid model, as produced by
 * {@link Qsys_ggisgi_fluid}.
 *
 * <p>Port of the struct returned by MATLAB qsys_ggisgi_fluid.m.
 *
 * @since LINE 3.1.0
 */
public class QsysFluidAbandonResult {
    /** "underloaded", "balanced" or "overloaded". */
    public final String regime;
    /** rho = lambda/(s*mu). */
    public final double trafficIntensity;
    /** w, the wait of every customer who is served; 0 unless overloaded. */
    public final double offeredWait;
    /** E[W] over all customers, int_0^w F^c(t)dt = m_a F_e(w). */
    public final double meanWait;
    /** w again, the fluid wait of served customers being deterministic. */
    public final double meanWaitServed;
    /** E[patience | patience &lt;= w]. */
    public final double meanWaitAbandon;
    /** 1 - 1/rho when overloaded, 0 otherwise. */
    public final double probAbandon;
    /** Q = lambda * meanWait, in customers. */
    public final double meanQueueLength;
    /** B = min(lambda/mu, s), in customers. */
    public final double meanNumberInService;
    /** B + Q. */
    public final double meanNumber;
    /** min(rho,1). */
    public final double utilization;
    /** min(lambda, s*mu). */
    public final double throughput;
    /** lambda - throughput. */
    public final double abandonRate;
    /** Ages at which the densities were evaluated, null when not requested. */
    public final double[] agePoints;
    /** b(x) per server at those ages, null when not requested. */
    public final double[] serviceAgeDensity;
    /** q(x) per server at those ages, null when not requested. */
    public final double[] queueAgeDensity;

    public QsysFluidAbandonResult(String regime, double trafficIntensity, double offeredWait,
                                  double meanWait, double meanWaitServed, double meanWaitAbandon,
                                  double probAbandon, double meanQueueLength,
                                  double meanNumberInService, double meanNumber,
                                  double utilization, double throughput, double abandonRate,
                                  double[] agePoints, double[] serviceAgeDensity,
                                  double[] queueAgeDensity) {
        this.regime = regime;
        this.trafficIntensity = trafficIntensity;
        this.offeredWait = offeredWait;
        this.meanWait = meanWait;
        this.meanWaitServed = meanWaitServed;
        this.meanWaitAbandon = meanWaitAbandon;
        this.probAbandon = probAbandon;
        this.meanQueueLength = meanQueueLength;
        this.meanNumberInService = meanNumberInService;
        this.meanNumber = meanNumber;
        this.utilization = utilization;
        this.throughput = throughput;
        this.abandonRate = abandonRate;
        this.agePoints = agePoints;
        this.serviceAgeDensity = serviceAgeDensity;
        this.queueAgeDensity = queueAgeDensity;
    }

    @Override
    public String toString() {
        return "QsysFluidAbandonResult(regime=" + regime
                + ", trafficIntensity=" + trafficIntensity
                + ", offeredWait=" + offeredWait
                + ", probAbandon=" + probAbandon
                + ", meanQueueLength=" + meanQueueLength + ")";
    }
}
