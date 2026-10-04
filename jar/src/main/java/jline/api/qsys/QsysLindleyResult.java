/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Conditional waiting-time moments of one Lindley step.
 *
 * <p>Returned by {@link Qsys_mm1_lindley} and {@link Qsys_hh1_lindley}. Every
 * array is indexed by the supplied current waiting times, so
 * {@code moments[i][m-1]} is the m-th conditional raw moment at {@code Wn[i]}.
 *
 * <p>Port of the MATLAB struct returned by qsys_mm1_lindley.m.
 *
 * @since LINE 3.1.0
 */
public class QsysLindleyResult {
    /** Conditional mean, one entry per supplied current waiting time. */
    public final double[] mean;
    /** Conditional variance, one entry per supplied current waiting time. */
    public final double[] var;
    /** Conditional raw moments, {@code moments[i][m-1]} for order m. */
    public final double[][] moments;
    /** Highest moment order computed. */
    public final int mmax;

    public QsysLindleyResult(double[] mean, double[] var, double[][] moments, int mmax) {
        this.mean = mean;
        this.var = var;
        this.moments = moments;
        this.mmax = mmax;
    }
}
