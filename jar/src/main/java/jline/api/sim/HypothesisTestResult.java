/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

/**
 * Outcome of one of the output-analysis hypothesis tests.
 *
 * <p>Shared by {@link Sim_vonneumann} and {@link Sim_shapirowilk}; {@link
 * #statistic} is the von Neumann ratio in the first case and the Shapiro-Wilk W
 * in the second.
 *
 * @since LINE 3.1.0
 */
public class HypothesisTestResult {
    /** The test statistic. */
    public final double statistic;
    /** The standardized statistic, NaN when the null law is used exactly. */
    public final double zscore;
    /** The p-value. */
    public final double pvalue;
    /** True when the null hypothesis is rejected at the requested level. */
    public final boolean reject;
    /** Number of observations the test used. */
    public final int nobs;

    public HypothesisTestResult(double statistic, double zscore, double pvalue,
                                boolean reject, int nobs) {
        this.statistic = statistic;
        this.zscore = zscore;
        this.pvalue = pvalue;
        this.reject = reject;
        this.nobs = nobs;
    }
}
