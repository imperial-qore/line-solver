/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

/**
 * Von Neumann ratio test for randomness of a sequence.
 *
 * <p>The statistic is the ratio of the mean square successive difference to the
 * variance,
 *
 * <pre>
 *   ratio = sum_{i=1}^{b-1} (x_{i+1}-x_i)^2 / sum_{i=1}^{b} (x_i - xbar)^2.
 * </pre>
 *
 * <p>Under the null hypothesis that the sequence is i.i.d. normal the ratio has
 * mean 2 and variance {@code 4(b-2)/((b-1)(b+1))}, and the standardized
 * statistic is asymptotically normal, giving a two-sided p-value. Serial
 * correlation of either sign moves the ratio away from 2: positive correlation
 * shrinks the successive differences and pushes the ratio below 2, negative
 * correlation pushes it above. Those null moments were confirmed by Monte Carlo
 * over b = 10, 16, 24, 32, 50 to within 0.3%.
 *
 * <p>Port of MATLAB sim_vonneumann.m.
 *
 * <p>Reference: J. von Neumann, "Distribution of the Ratio of the Mean Square
 * Successive Difference to the Variance", Ann. Math. Statist. 12(4), 1941;
 * L. C. Young, "Randomness in Ordered Sequences", Ann. Math. Statist. 12, 1941.
 *
 * @since LINE 3.1.0
 */
public final class Sim_vonneumann {
    private Sim_vonneumann() {}

    /**
     * Tests a sequence for randomness at the 5% level.
     *
     * @param x the sequence, at least 3 finite and not all equal values
     * @return the test outcome
     */
    public static HypothesisTestResult sim_vonneumann(double[] x) {
        return sim_vonneumann(x, 0.05);
    }

    /**
     * Tests a sequence for randomness.
     *
     * @param x     the sequence, at least 3 finite and not all equal values
     * @param alpha significance level in (0,1)
     * @return the test outcome
     */
    public static HypothesisTestResult sim_vonneumann(double[] x, double alpha) {
        if (x == null) {
            throw new IllegalArgumentException("The sequence must not be null");
        }
        if (!(alpha > 0.0) || !(alpha < 1.0)) {
            throw new IllegalArgumentException("alpha=" + alpha + " must lie in (0,1)");
        }
        int b = x.length;
        if (b < 3) {
            throw new IllegalArgumentException(
                    "At least 3 observations are required, got " + b);
        }

        double mean = 0.0;
        for (int i = 0; i < b; i++) {
            if (!Double.isFinite(x[i])) {
                throw new IllegalArgumentException("The sequence must be finite");
            }
            mean += x[i];
        }
        mean /= b;

        double den = 0.0;
        for (int i = 0; i < b; i++) {
            double d = x[i] - mean;
            den += d * d;
        }
        if (!(den > 0.0)) {
            throw new IllegalArgumentException(
                    "The sequence is constant, the ratio is undefined");
        }

        double num = 0.0;
        for (int i = 0; i < b - 1; i++) {
            double d = x[i + 1] - x[i];
            num += d * d;
        }

        double ratio = num / den;
        double sd = Math.sqrt(4.0 * (b - 2) / ((double) (b - 1) * (b + 1)));
        double zscore = (ratio - 2.0) / sd;
        double pvalue = 2.0 * (1.0 - SimDist.normcdf(Math.abs(zscore)));

        return new HypothesisTestResult(ratio, zscore, pvalue, pvalue < alpha, b);
    }
}
