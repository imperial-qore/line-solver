/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

/**
 * Trajectory of the Gt/Mt/st+GI many-server fluid queue, as produced by
 * {@link Qsys_gtmtst_fluid}.
 *
 * <p>Port of the struct returned by MATLAB qsys_gtmtst_fluid.m. Every array is
 * indexed by the time grid.
 *
 * @since LINE 3.1.0
 */
public class QsysTvFluidResult {
    /** The time grid. */
    public final double[] times;
    /** 1 while the queue is overloaded, 0 while it is underloaded. */
    public final int[] regime;
    /** Fluid in service. */
    public final double[] B;
    /** Fluid in queue. */
    public final double[] Q;
    /** Total fluid, B + Q. */
    public final double[] X;
    /** Boundary waiting time, the age of the oldest fluid still waiting. */
    public final double[] w;
    /** Potential waiting time of a quantum arriving now with infinite patience. */
    public final double[] v;
    /** Service completion rate mu(t)B(t). */
    public final double[] sigma;
    /** Abandonment rate. */
    public final double[] alpha;
    /** B/s, the fraction of capacity in use. */
    public final double[] utilization;
    /** The arrival rate on the grid. */
    public final double[] arrivalRate;
    /** The staffing on the grid. */
    public final double[] staffing;
    /** Gamma(t) = s'(t) + s(t)mu(t), the rate at which capacity frees up. */
    public final double[] capacityRate;

    public QsysTvFluidResult(double[] times, int[] regime, double[] B, double[] Q, double[] X,
                             double[] w, double[] v, double[] sigma, double[] alpha,
                             double[] utilization, double[] arrivalRate, double[] staffing,
                             double[] capacityRate) {
        this.times = times;
        this.regime = regime;
        this.B = B;
        this.Q = Q;
        this.X = X;
        this.w = w;
        this.v = v;
        this.sigma = sigma;
        this.alpha = alpha;
        this.utilization = utilization;
        this.arrivalRate = arrivalRate;
        this.staffing = staffing;
        this.capacityRate = capacityRate;
    }

    @Override
    public String toString() {
        return "QsysTvFluidResult(points=" + (times == null ? 0 : times.length)
                + ", finalB=" + (B == null || B.length == 0 ? Double.NaN : B[B.length - 1])
                + ", finalQ=" + (Q == null || Q.length == 0 ? Double.NaN : Q[Q.length - 1]) + ")";
    }
}
