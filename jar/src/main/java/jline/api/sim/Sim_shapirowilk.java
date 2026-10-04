/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

import java.util.Arrays;

/**
 * Shapiro-Wilk test for univariate normality, Royston's AS R94 algorithm.
 *
 * <p>The statistic is
 *
 * <pre>
 *   W = (sum_i a_i x_(i))^2 / sum_i (x_i - xbar)^2,
 * </pre>
 *
 * <p>with {@code x_(i)} the order statistics and {@code a} the antisymmetric
 * weight vector obtained by correcting the normalized expected normal order
 * statistics {@code m_i = Phi^{-1}((i-3/8)/(n+1/4))} in their two extreme
 * components. Small W means departure from normality, so the test is one-sided
 * in W and the p-value is an upper normal tail after Royston's normalizing
 * transform, which has three branches: n = 3 exact, 4 &lt;= n &lt;= 11, and
 * n &gt;= 12. Valid for 3 &lt;= n &lt;= 5000.
 *
 * <p>Port of MATLAB sim_shapirowilk.m. W and the p-value agree with
 * {@code scipy.stats.shapiro} to 5e-10 and 1.5e-7 respectively over n up to 2000.
 *
 * <p>Reference: J. P. Royston, "Approximating the Shapiro-Wilk W-test for
 * Non-normality", Statistics and Computing 2, 1992; J. P. Royston, "Remark
 * AS R94", Applied Statistics 44(4), 1995.
 *
 * @since LINE 3.1.0
 */
public final class Sim_shapirowilk {
    private static final double[] C1 =
            {0.0, 0.221157, -0.147981, -2.071190, 4.434685, -2.706056};
    private static final double[] C2 =
            {0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633};

    private Sim_shapirowilk() {}

    /**
     * Tests a sample for normality at the 5% level.
     *
     * @param x the sample, 3 to 5000 finite and not all equal values
     * @return the test outcome
     */
    public static HypothesisTestResult sim_shapirowilk(double[] x) {
        return sim_shapirowilk(x, 0.05);
    }

    /**
     * Tests a sample for normality.
     *
     * @param x     the sample, 3 to 5000 finite and not all equal values
     * @param alpha significance level in (0,1)
     * @return the test outcome
     */
    public static HypothesisTestResult sim_shapirowilk(double[] x, double alpha) {
        if (x == null) {
            throw new IllegalArgumentException("The sample must not be null");
        }
        if (!(alpha > 0.0) || !(alpha < 1.0)) {
            throw new IllegalArgumentException("alpha=" + alpha + " must lie in (0,1)");
        }
        int n = x.length;
        if (n < 3) {
            throw new IllegalArgumentException(
                    "At least 3 observations are required, got " + n);
        }
        if (n > 5000) {
            throw new IllegalArgumentException(
                    "The AS R94 approximation is valid up to n = 5000, got " + n);
        }

        double[] sorted = new double[n];
        double mean = 0.0;
        for (int i = 0; i < n; i++) {
            if (!Double.isFinite(x[i])) {
                throw new IllegalArgumentException("The sample must be finite");
            }
            sorted[i] = x[i];
            mean += x[i];
        }
        Arrays.sort(sorted);
        mean /= n;

        double ssd = 0.0;
        for (int i = 0; i < n; i++) {
            double d = sorted[i] - mean;
            ssd += d * d;
        }
        if (!(ssd > 0.0)) {
            throw new IllegalArgumentException("The sample is constant, W is undefined");
        }

        double[] a = weights(n);
        double dot = 0.0;
        for (int i = 0; i < n; i++) {
            dot += a[i] * sorted[i];
        }
        double W = Math.min(dot * dot / ssd, 1.0);

        double zscore;
        double pvalue;
        if (n == 3) {
            // exact null law, W is supported on [3/4, 1]
            pvalue = 6.0 / Math.PI * (Math.asin(Math.sqrt(W)) - Math.asin(Math.sqrt(0.75)));
            pvalue = Math.min(Math.max(pvalue, 0.0), 1.0);
            zscore = Double.NaN;
        } else {
            double w;
            double mu;
            double sigma;
            if (n <= 11) {
                double g = -2.273 + 0.459 * n;
                w = -Math.log(g - Math.log(1.0 - W));
                mu = 0.5440 - 0.39978 * n + 0.025054 * n * n - 0.0006714 * n * n * n;
                sigma = Math.exp(1.3822 - 0.77857 * n + 0.062767 * n * n
                        - 0.0020322 * n * n * n);
            } else {
                double ln = Math.log(n);
                w = Math.log(1.0 - W);
                mu = -1.5861 - 0.31082 * ln - 0.083751 * ln * ln + 0.0038915 * ln * ln * ln;
                sigma = Math.exp(-0.4803 - 0.082676 * ln + 0.0030302 * ln * ln);
            }
            zscore = (w - mu) / sigma;
            pvalue = 1.0 - SimDist.normcdf(zscore);
        }

        return new HypothesisTestResult(W, zscore, pvalue, pvalue < alpha, n);
    }

    /**
     * Royston AS R94 antisymmetric weight vector, {@code a[n-1-i] = -a[i]}.
     *
     * @param n sample size, at least 3
     * @return the weight vector, ascending with the order statistics
     */
    public static double[] weights(int n) {
        if (n < 3) {
            throw new IllegalArgumentException("At least 3 observations are required, got " + n);
        }
        double[] a = new double[n];
        if (n == 3) {
            a[0] = -Math.sqrt(0.5);
            a[1] = 0.0;
            a[2] = Math.sqrt(0.5);
            return a;
        }

        double[] m = new double[n];
        double mm = 0.0;
        for (int i = 1; i <= n; i++) {
            m[i - 1] = SimDist.norminv((i - 0.375) / (n + 0.25));
            mm += m[i - 1] * m[i - 1];
        }
        double root = Math.sqrt(mm);
        double u = 1.0 / Math.sqrt(n);

        System.arraycopy(m, 0, a, 0, n);
        double an = m[n - 1] / root + poly(C1, u);
        if (n > 5) {
            double anm1 = m[n - 2] / root + poly(C2, u);
            double phi = (mm - 2.0 * m[n - 1] * m[n - 1] - 2.0 * m[n - 2] * m[n - 2])
                    / (1.0 - 2.0 * an * an - 2.0 * anm1 * anm1);
            double sphi = Math.sqrt(phi);
            for (int i = 2; i <= n - 3; i++) {
                a[i] = m[i] / sphi;
            }
            a[n - 1] = an;
            a[n - 2] = anm1;
            a[0] = -an;
            a[1] = -anm1;
        } else {
            double phi = (mm - 2.0 * m[n - 1] * m[n - 1]) / (1.0 - 2.0 * an * an);
            double sphi = Math.sqrt(phi);
            for (int i = 1; i <= n - 2; i++) {
                a[i] = m[i] / sphi;
            }
            a[n - 1] = an;
            a[0] = -an;
        }
        return a;
    }

    private static double poly(double[] c, double x) {
        double v = 0.0;
        double p = 1.0;
        for (int i = 0; i < c.length; i++) {
            v += c[i] * p;
            p *= x;
        }
        return v;
    }
}
