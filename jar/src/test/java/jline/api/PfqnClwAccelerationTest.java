/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.nc.Pfqn_clw;
import jline.io.Ret;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for the two speed-ups of Choudhury, Leung and Whitt,
 * "Calculating normalization constants of closed queueing networks by
 * numerically inverting their generating functions", J. ACM 42(5):935-970,
 * 1995, as applied by {@code Pfqn_clw}: Euler summation of the inner sums
 * (Section 2.4) and dimension reduction by decomposition (Section 3).
 *
 * <p>The expected values are the paper's own Tables I and III, which pin the
 * algorithm rather than the port; they agree across MATLAB, JAR, Python native
 * and C++. Table III is out of reach without the reduction: p = 11 chains cost
 * prod_j 2 l_j K_j contour points, and the reduction takes the inversion to
 * dimension 2 by inverting the hub chain and then each leaf separately.
 *
 * <p>Mirrors line-test.git/test/testsAPI/test_pfqn_clw_acceleration.m.
 */
public class PfqnClwAccelerationTest {

    private static final double LOG10 = Math.log(10.0);

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int j = 0; j < v.length; j++) m.set(0, j, v[j]);
        return m;
    }

    private static Matrix col(double... v) {
        Matrix m = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) m.set(i, 0, v[i]);
        return m;
    }

    /** Example 8.1: p = 1, q' = 10 distinct queues of multiplicity 5. */
    private static Matrix example81L() {
        double[] v = new double[10];
        for (int i = 0; i < 10; i++) v[i] = 0.1 * (i + 1);
        return col(v);
    }

    private static Matrix mult5() {
        double[] v = new double[10];
        java.util.Arrays.fill(v, 5.0);
        return col(v);
    }

    /**
     * Examples 8.3 and 8.4: chain 1 is a hub visiting all ten queues, chain g+1
     * has queue g to itself. The interdependence graph is a star, so D = {1}
     * leaves ten single-variable components (the paper's Figure 1).
     */
    private static Matrix starL() {
        Matrix L = new Matrix(10, 11);
        for (int g = 1; g <= 10; g++) {
            L.set(g - 1, 0, 1 + 0.1 * g);
            L.set(g - 1, g, 0.1 * g);
        }
        return L;
    }

    private static Matrix starZ() {
        double[] z = new double[11];
        z[0] = 50.0;
        for (int g = 1; g <= 10; g++) z[g] = 5.0 * (g + 1) - 10.0;
        return row(z);
    }

    @Test
    public void testEulerReproducesTableI() {
        // Table I: every population from 2 to 2e7. Without Euler summation the
        // last row alone would cost 4e7 contour points.
        int[] K = {2, 20, 200, 2000, 20000, 200000, 2000000, 20000000};
        double[][] ref = {{5.377500, 2}, {1.906584, 13}, {1.381312, 26}, {1.284918, 31},
                          {1.541538, 35}, {1.569301, 39}, {1.572100, 43}, {1.572380, 47}};
        for (int t = 0; t < K.length; t++) {
            Ret.pfqnNc r = Pfqn_clw.pfqn_clw(example81L(), row(K[t]), row(5.0), mult5());
            assertEquals(Math.log10(ref[t][0]) + ref[t][1], r.lG / LOG10, 1e-6,
                    "Table I at K1 = " + K[t]);
        }
    }

    @Test
    public void testEulerAgreesWithTheExactSum() {
        // The acceleration is adaptive: the Euler order is doubled until the
        // paper's own estimate |E(m,n) - E(m,n+1)| settles, so switching it off
        // must not move the answer. K1 = 40 is past the n+m = 31 threshold.
        Pfqn_clw.Options off = new Pfqn_clw.Options();
        off.euler = false;
        for (int K : new int[]{40, 100, 250}) {
            Ret.pfqnNc a = Pfqn_clw.pfqn_clw(example81L(), row(K), row(5.0), mult5());
            Ret.pfqnNc b = Pfqn_clw.pfqn_clw(example81L(), row(K), row(5.0), mult5(),
                    null, null, off);
            assertEquals(b.lG, a.lG, 1e-9 * Math.abs(b.lG),
                    "Euler summation moved the answer at K1 = " + K);
        }
    }

    @Test
    public void testDimredReproducesTableIII() {
        // Table III rows 1-4: eleven chains, which only the reduction reaches.
        // Mantissa and exponent are kept apart: 1.937826e683 is not a double.
        int[] K1 = {2, 20, 200, 2000};
        double[][] ref = {{1.235628, 25}, {7.503087, 45}, {5.970503, 129}, {1.937826, 683}};
        for (int t = 0; t < K1.length; t++) {
            double[] N = new double[11];
            java.util.Arrays.fill(N, 2.0);
            N[0] = K1[t];
            Ret.pfqnNc r = Pfqn_clw.pfqn_clw(starL(), row(N), starZ(), mult5());
            assertEquals(Math.log10(ref[t][0]) + ref[t][1], r.lG / LOG10, 1e-5,
                    "Table III at K1 = " + K1[t]);
        }
    }

    @Test
    public void testDimredRows5To8UnderThePaperScaleTuning() {
        // Table III rows 5-8 need the manual tuning of page 956 on the hub
        // chain, beta in [0.8, 1.2], which the paper prescribes for its largest
        // examples.
        int[] K1 = {2, 20, 200, 2000};
        double[] b1 = {0.8, 0.8, 0.8, 0.95};
        double[][] ref = {{3.004462, 107}, {1.677866, 133}, {8.032122, 260}, {1.617153, 926}};
        for (int t = 0; t < K1.length; t++) {
            double[] N = new double[11];
            N[0] = K1[t];
            for (int g = 1; g <= 10; g++) N[g] = 5.0 * g;
            Pfqn_clw.Options opt = new Pfqn_clw.Options();
            opt.beta = new double[11];
            java.util.Arrays.fill(opt.beta, 1.0);
            opt.beta[0] = b1[t];
            Ret.pfqnNc r = Pfqn_clw.pfqn_clw(starL(), row(N), starZ(), mult5(), null, null, opt);
            assertEquals(Math.log10(ref[t][0]) + ref[t][1], r.lG / LOG10, 1e-3,
                    "Table III row for K1 = " + K1[t]);
        }
    }

    @Test
    public void testDimredAgreesWithTheFullInversion() {
        // Reduction reorders the inversion and splits the factors; it must not
        // move the answer. Three leaves keep the undecomposed inversion
        // affordable, and exact convolution adjudicates both.
        Matrix L = new Matrix(3, 4);
        double[] z = new double[4];
        z[0] = 5.0;
        for (int g = 1; g <= 3; g++) {
            L.set(g - 1, 0, 1 + 0.1 * g);
            L.set(g - 1, g, 0.1 * g);
            z[g] = 5.0 * g - 5.0;
        }
        Matrix N = row(4, 3, 3, 3);
        Pfqn_clw.Options off = new Pfqn_clw.Options();
        off.dimred = false;
        double lca = jline.api.pfqn.nc.Pfqn_ca.pfqn_ca(L, N, row(z)).lG;
        double lred = Pfqn_clw.pfqn_clw(L, N, row(z)).lG;
        double lfull = Pfqn_clw.pfqn_clw(L, N, row(z), null, null, null, off).lG;
        assertEquals(lca, lred, 1e-7 * Math.abs(lca));
        assertEquals(lca, lfull, 1e-6 * Math.abs(lca));
    }

    @Test
    public void testDimredIsInertOnACoupledModel() {
        // Every chain visits every queue, so the interdependence graph is
        // complete, no subset D reduces the dimension, and the classical path
        // must be taken unchanged -- to the last bit, not merely to a tolerance.
        Matrix L = new Matrix(2, 3);
        L.set(0, 0, 0.1); L.set(0, 1, 0.2); L.set(0, 2, 0.15);
        L.set(1, 0, 0.3); L.set(1, 1, 0.05); L.set(1, 2, 0.1);
        Matrix N = row(2, 2, 1);
        Matrix Z = row(1.0, 0.5, 0.2);
        Pfqn_clw.Options off = new Pfqn_clw.Options();
        off.dimred = false;
        double a = Pfqn_clw.pfqn_clw(L, N, Z).lG;
        double b = Pfqn_clw.pfqn_clw(L, N, Z, null, null, null, off).lG;
        assertTrue(a == b, "dimension reduction changed a model it cannot reduce");
        // the MATLAB value, at the spread this model already has across the
        // codebases: p = 3 doubles the nesting depth and the cancellation with
        // it, and the ports sit ~4e-7 apart on it (see test_pfqn_mvac_oi_clw.cpp)
        assertEquals(-1.441982188151596, a, 1e-6);
    }

    @Test
    public void testDisconnectedModelNeedsNoCommittedVariable() {
        // Two chains that share no queue: the graph is already disconnected, so
        // the reduction is exact with D empty and the constant is the product of
        // the two single-chain constants.
        Matrix L = new Matrix(4, 2);
        L.set(0, 0, 0.4); L.set(1, 0, 0.2); L.set(2, 1, 0.5); L.set(3, 1, 0.3);
        double both = Pfqn_clw.pfqn_clw(L, row(3, 4), row(1.0, 2.0)).lG;
        double one = Pfqn_clw.pfqn_clw(col(0.4, 0.2), row(3), row(1.0)).lG;
        double two = Pfqn_clw.pfqn_clw(col(0.5, 0.3), row(4), row(2.0)).lG;
        assertEquals(one + two, both, 1e-9 * Math.abs(one + two));
    }
}
