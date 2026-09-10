/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Point estimate and confidence interval for a steady-state quantile.
 *
 * <p>Returned by {@link Sim_fquest} and {@link Sim_firquest}. When {@link
 * #heuristic} is true a stage test failed, the sample was too small for the
 * asymptotic justification, and the interval is the conservative fallback, which
 * may be asymmetric about the estimate.
 *
 * @since LINE 3.1.0
 */
public class QuantileCIResult {
    /** Full-sample empirical p-quantile of the truncated data. */
    public final double estimate;
    /** Lower confidence limit. */
    public final double lower;
    /** Upper confidence limit. */
    public final double upper;
    /** Half of the interval width; the asymmetric fallback attains it only on average. */
    public final double halfwidth;
    /** Final batch count, per replication for Sim_firquest. */
    public final int b;
    /** Final batch size. */
    public final int m;
    /** Observations used, b*m, or R*b*m for Sim_firquest. */
    public final int n;
    /** Number of replications, 1 for Sim_fquest. */
    public final int R;
    /** Observations deleted from the front of each sample path. */
    public final int truncated;
    /** Batched STS area estimator of the variance parameter. */
    public final double Ap;
    /** Nonoverlapping batched quantile estimator of the variance parameter. */
    public final double Np;
    /** Combined estimator of the variance parameter. */
    public final double Vp;
    /** True when a stage test failed and the interval is not asymptotically justified. */
    public final boolean heuristic;
    /** Diagnostic messages, empty on a clean run. */
    public final List<String> warnings;

    public QuantileCIResult(double estimate, double lower, double upper, double halfwidth,
                            int b, int m, int n, int R, int truncated,
                            double Ap, double Np, double Vp, boolean heuristic,
                            List<String> warnings) {
        this.estimate = estimate;
        this.lower = lower;
        this.upper = upper;
        this.halfwidth = halfwidth;
        this.b = b;
        this.m = m;
        this.n = n;
        this.R = R;
        this.truncated = truncated;
        this.Ap = Ap;
        this.Np = Np;
        this.Vp = Vp;
        this.heuristic = heuristic;
        this.warnings = Collections.unmodifiableList(
                warnings == null ? new ArrayList<String>() : new ArrayList<String>(warnings));
    }
}
