/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

/**
 * Standardized time series statistics of the batched quantile process.
 *
 * <p>Produced by {@link Sim_sts_quantile_areas} and consumed by {@link
 * Sim_fquest} and {@link Sim_firquest}. The three variance-parameter estimators
 * all target {@code sigma_p^2 = lim n Var(ytilde_p(n))}; see {@link
 * Sim_sts_quantile_areas} for their definitions.
 *
 * @since LINE 3.1.0
 */
public class StsQuantileStats {
    /** Signed STS areas, one per batch. */
    public final double[] areas;
    /** Batched quantile estimators, one per batch. */
    public final double[] bqe;
    /** Full-sample empirical p-quantile over all n observations. */
    public final double quantile;
    /** Batched STS area estimator of the variance parameter. */
    public final double Ap;
    /** Nonoverlapping batched quantile estimator of the variance parameter, NaN when b &lt; 2. */
    public final double Np;
    /** Combined estimator of the variance parameter, NaN when b &lt; 2. */
    public final double Vp;
    /** Batch count. */
    public final int b;
    /** Batch size. */
    public final int m;
    /** Observations used, b*m. */
    public final int n;

    public StsQuantileStats(double[] areas, double[] bqe, double quantile,
                            double Ap, double Np, double Vp, int b, int m, int n) {
        this.areas = areas;
        this.bqe = bqe;
        this.quantile = quantile;
        this.Ap = Ap;
        this.Np = Np;
        this.Vp = Vp;
        this.b = b;
        this.m = m;
        this.n = n;
    }
}
