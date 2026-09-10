/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.junit.jupiter.api.Test;

import jline.util.Maths;

import static jline.TestTools.MID_TOL;
import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Validates the Choudhury-Leung-Whitt (CLW) normalizing-constant algorithm
 * against Example 8.4 of Choudhury, Leung and Whitt, "Calculating Normalization
 * Constants of Closed Queuing Networks by Numerically Inverting Their
 * Generating Functions," J. ACM 42(5):935-970, 1995 (Table IV, Case 1) --
 * the paper's largest and most challenging numerical example.
 *
 * <p>The model has p=11 closed chains and q'=10 distinct single-server queues
 * (star structure: chain 1 is a hub visiting every queue, chains 2..11 each
 * have a private queue group), with multiplicities m_i=100 and populations
 * K_j=200 for j=2..11. Direct application of {@code Pfqn_clw.pfqn_clw} would
 * require an 11-dimensional nested lattice-Poisson inversion, which is
 * computationally intractable at this scale (cost is exponential in the
 * number of chains) -- exactly why the paper introduces a model-specific
 * dimension-reduction closed form (Section 5.4, eqs. 5.31-5.38) that inverts
 * only the hub chain's dimension, with the other 10 chains' contributions
 * folded in analytically.
 *
 * <p>This test implements that closed form directly (not via
 * {@code Pfqn_clw}, which does not include this model-specific reduction)
 * and evaluates it in the log domain, since g(K) reaches magnitudes up to
 * ~1e8575 -- far beyond double range. All quantities in the reduction are
 * non-negative, so log-domain accumulation via log-sum-exp is exact up to
 * floating-point rounding, and this is fast and requires no complex
 * arithmetic or numerical contour inversion.
 */
public class PfqnClwExample84Test {

    private static final int P = 10;
    private static final double RHO10 = 50.0;
    private static final double[] RHO_G0 = new double[P];
    private static final double[] RHO_GG = new double[P];
    private static final double[] RHO_1G = new double[P];
    private static final int MG = 100;
    private static final int KG = 200;

    static {
        for (int g = 1; g <= P; g++) {
            RHO_G0[g - 1] = 5.0 * (g + 1) - 10.0;
            RHO_GG[g - 1] = 0.1 * g;
            RHO_1G[g - 1] = 1.0 + 0.1 * g;
        }
    }

    /**
     * log([z1^j] S_g(z1)) for j=0..K1, where S_g is group g's closed-form
     * marginal generating function (eq. 5.35, single group), obtained via the
     * binomial-series expansion of each (1-rho_1g z1)^{-(mg+Kg-k)} term.
     */
    private static double[] groupLogSeries(int gIdx, int K1) {
        double r0 = RHO_G0[gIdx];
        double rg = RHO_GG[gIdx];
        double rho1 = RHO_1G[gIdx];
        double logRho1 = FastMath.log(rho1);

        double[] acc = new double[K1 + 1];
        java.util.Arrays.fill(acc, Double.NEGATIVE_INFINITY);

        int kmax = (r0 > 0) ? KG : 0;  // r0=0 => only the k=0 term is finite (0^k=0 for k>0)
        for (int k = 0; k <= kmax; k++) {
            int n = MG + KG - k;
            double logCk = (k == 0 ? 0.0 : k * FastMath.log(r0))
                    - Gamma.logGamma(k + 1)
                    + Maths.logBinomial(MG + KG - k - 1, KG - k)
                    + (KG - k) * FastMath.log(rg);
            for (int j = 0; j <= K1; j++) {
                double logTerm = logCk + Maths.logBinomial(n + j - 1, j) + j * logRho1;
                acc[j] = logAdd(acc[j], logTerm);
            }
        }
        return acc;
    }

    /** log(rho10^i / i!) for i=0..K1, the exp(rho10*z1) prefactor series. */
    private static double[] logExpSeries(int K1) {
        double[] s = new double[K1 + 1];
        double logRho10 = FastMath.log(RHO10);
        for (int i = 0; i <= K1; i++) {
            s[i] = i * logRho10 - Gamma.logGamma(i + 1);
        }
        return s;
    }

    private static double logAdd(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b;
        if (b == Double.NEGATIVE_INFINITY) return a;
        double m = Math.max(a, b);
        return m + FastMath.log(FastMath.exp(a - m) + FastMath.exp(b - m));
    }

    /** Log-domain convolution of two power series, truncated to degree K1. */
    private static double[] logConvolve(double[] a, double[] b, int K1) {
        double[] out = new double[K1 + 1];
        for (int m = 0; m <= K1; m++) {
            double max = Double.NEGATIVE_INFINITY;
            for (int idx = 0; idx <= m; idx++) {
                double v = a[idx] + b[m - idx];
                if (v > max) max = v;
            }
            if (max == Double.NEGATIVE_INFINITY) {
                out[m] = Double.NEGATIVE_INFINITY;
                continue;
            }
            double sum = 0;
            for (int idx = 0; idx <= m; idx++) {
                sum += FastMath.exp(a[idx] + b[m - idx] - max);
            }
            out[m] = max + FastMath.log(sum);
        }
        return out;
    }

    /** log10(g(K1, 200, 200, ..., 200)) via the closed-form dimension reduction. */
    private static double log10G(int K1) {
        double[] conv = groupLogSeries(0, K1);
        for (int g = 1; g < P; g++) {
            conv = logConvolve(conv, groupLogSeries(g, K1), K1);
        }
        double[] full = logConvolve(conv, logExpSeries(K1), K1);
        return full[K1] / FastMath.log(10);
    }

    private static void assertMatchesPaper(int K1, double mantissa, int exponent) {
        double expectedLog10 = FastMath.log10(mantissa) + exponent;
        double actualLog10 = log10G(K1);
        assertEquals(expectedLog10, actualLog10, MID_TOL,
                String.format("Example 8.4 Case 1, K1=%d: expected log10(g)=%.6f, got=%.6f",
                        K1, expectedLog10, actualLog10));
    }

    @Test
    public void testExample84Case1_K1_20() {
        assertMatchesPaper(20, 1.232036, 278);
    }

    @Test
    public void testExample84Case1_K1_200() {
        assertMatchesPaper(200, 2.941740, 579);
    }

    @Test
    public void testExample84Case1_K1_2000() {
        assertMatchesPaper(2000, 3.399948, 2037);
    }

    /**
     * The paper's largest numerical example: g(K) ~ 9.07177e8575.
     * Runtime is O(P * K1^2) in the log-domain convolution (~15-20s at K1=20000).
     */
    @Test
    public void testExample84Case1_K1_20000() {
        assertMatchesPaper(20000, 9.07177, 8575);
    }
}
