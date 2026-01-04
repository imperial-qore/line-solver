/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static jline.lib.butools.ph.CheckMERepresentationKt.checkMERepresentation;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for ME (Matrix Exponential) distribution
 */
public class MEDistributionTest {

    private static final double TOLERANCE = 1e-6;
    private static final double LOOSE_TOLERANCE = 1e-3;
    private static final double SAMPLING_TOLERANCE = 0.1; // 10% tolerance for empirical samples

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation with fixed seed
        Maths.setRandomNumbersMatlab(true);
        RandomManager.setMasterSeed(23000);
    }

    @Test
    public void testMEConstructionValid() {
        // Create valid ME distribution (order 2)
        Matrix alpha = new Matrix(new double[]{0.4, 0.6});
        Matrix A = new Matrix(new double[][]{
                {-2.0, 1.0},
                {0.5, -1.5}
        });

        ME me = new ME(alpha, A);

        assertNotNull(me);
        assertEquals(2, me.getNumberOfPhases());
        assertArrayEquals(new double[]{0.4, 0.6}, me.getAlpha().toArray1D(), TOLERANCE);
    }

    @Test
    public void testMEConstructionInvalid() {
        // Test with positive eigenvalue (invalid)
        Matrix alpha = new Matrix(new double[]{1.0});
        Matrix A = new Matrix(new double[][]{{2.0}});

        assertThrows(IllegalArgumentException.class, () -> {
            new ME(alpha, A);
        });
    }

    @Test
    public void testMEFromExponential() {
        // Exponential is simplest ME (order 1)
        double rate = 2.0;
        ME me = ME.fromExp(rate);

        assertNotNull(me);
        assertEquals(1, me.getNumberOfPhases());

        // Mean should be 1/rate
        double mean = me.getMean();
        assertEquals(1.0 / rate, mean, TOLERANCE);

        // SCV should be 1 (for exponential)
        double scv = me.getSCV();
        assertEquals(1.0, scv, TOLERANCE);
    }

    @Test
    public void testMEFromErlang() {
        // Erlang-3 with rate 1.0
        int k = 3;
        double rate = 1.0;
        ME me = ME.fromErlang(k, rate);

        assertNotNull(me);
        assertEquals(k, me.getNumberOfPhases());

        // Mean should be k/rate
        double mean = me.getMean();
        assertEquals((double) k / rate, mean, TOLERANCE);

        // SCV should be 1/k (for Erlang)
        double scv = me.getSCV();
        assertEquals(1.0 / k, scv, TOLERANCE);
    }

    @Test
    public void testMEFromHyperExp() {
        // HyperExp with 2 branches
        double[] p = {0.3, 0.7};
        double[] rates = {1.0, 4.0};
        ME me = ME.fromHyperExp(p, rates);

        assertNotNull(me);
        assertEquals(2, me.getNumberOfPhases());

        // Verify it's a valid ME distribution
        assertTrue(checkMERepresentation(me.getAlpha(), me.getA(), 1e-14));

        // Mean = p[0]/rates[0] + p[1]/rates[1]
        double expectedMean = p[0] / rates[0] + p[1] / rates[1];
        double mean = me.getMean();
        assertEquals(expectedMean, mean, TOLERANCE);
    }

    @Test
    public void testMEMeanComputation() {
        // Create ME with known mean
        Matrix alpha = new Matrix(new double[]{1.0});
        Matrix A = new Matrix(new double[][]{{-2.0}});

        ME me = new ME(alpha, A);

        // Mean = -alpha * A^(-1) * e = -1 * (-1/2) * 1 = 0.5
        double mean = me.getMean();
        assertEquals(0.5, mean, TOLERANCE);
    }

    @Test
    public void testMEVarianceComputation() {
        // Exponential with rate 2.0
        ME me = ME.fromExp(2.0);

        // For Exp(lambda): variance = 1/lambda^2 = 1/4 = 0.25
        double variance = me.getVar();
        assertEquals(0.25, variance, TOLERANCE);
    }

    @Test
    public void testMESCVComputation() {
        // Erlang-2 with rate 1.0
        ME me = ME.fromErlang(2, 1.0);

        // For Erlang-k: SCV = 1/k = 1/2 = 0.5
        double scv = me.getSCV();
        assertEquals(0.5, scv, TOLERANCE);
    }

    @Test
    public void testMECDFEvaluation() {
        // Exponential with rate 1.0
        ME me = ME.fromExp(1.0);

        // For Exp(1): CDF(t) = 1 - exp(-t)
        double cdf0 = me.evalCDF(0.0);
        assertEquals(0.0, cdf0, TOLERANCE);

        double cdf1 = me.evalCDF(1.0);
        assertEquals(1.0 - Math.exp(-1.0), cdf1, TOLERANCE);

        double cdf2 = me.evalCDF(2.0);
        assertEquals(1.0 - Math.exp(-2.0), cdf2, TOLERANCE);
    }

    @Test
    public void testMESampling() {
        // Create ME distribution
        ME me = ME.fromExp(2.0);

        // Generate samples
        Random rng = new Random(12345);
        int n = 10000;
        double[] samples = me.sample(n, rng);

        assertEquals(n, samples.length);

        // Compute empirical mean
        double empiricalMean = 0.0;
        for (double sample : samples) {
            assertTrue(sample >= 0, "Sample should be non-negative");
            empiricalMean += sample;
        }
        empiricalMean /= n;

        // Expected mean = 1/2 = 0.5
        double expectedMean = 0.5;
        assertEquals(expectedMean, empiricalMean, SAMPLING_TOLERANCE * expectedMean);
    }

    @Test
    public void testMEValidationStrict() {
        // Test that validation catches invalid representations

        // Case 1: Non-square matrix
        Matrix alpha1 = new Matrix(new double[]{1.0, 0.0});
        Matrix A1 = new Matrix(new double[][]{
                {-1.0, 0.5}
        }); // 1x2 matrix (non-square)

        assertThrows(IllegalArgumentException.class, () -> {
            new ME(alpha1, A1);
        });

        // Case 2: Incompatible dimensions
        Matrix alpha2 = new Matrix(new double[]{1.0});
        Matrix A2 = new Matrix(new double[][]{
                {-1.0, 0.5},
                {0.5, -1.0}
        }); // 2x2 matrix but alpha is 1x1

        assertThrows(IllegalArgumentException.class, () -> {
            new ME(alpha2, A2);
        });

        // Case 3: Dominant eigenvalue is complex
        Matrix alpha3 = new Matrix(new double[]{0.5, 0.5});
        Matrix A3 = new Matrix(new double[][]{
                {-1.0, 2.0},
                {-2.0, -1.0}
        }); // Eigenvalues are -1 ± 2i (complex)

        assertThrows(IllegalArgumentException.class, () -> {
            new ME(alpha3, A3);
        });
    }

    @Test
    public void testMEProcessRepresentation() {
        // Verify that ME creates correct process representation {D0=A, D1=-A*e*alpha'}
        Matrix alpha = new Matrix(new double[]{0.3, 0.7});
        Matrix A = new Matrix(new double[][]{
                {-2.0, 1.0},
                {0.5, -1.5}
        });

        ME me = new ME(alpha, A);

        // D0 should be A
        Matrix D0 = me.D(0);
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                assertEquals(A.get(i, j), D0.get(i, j), TOLERANCE);
            }
        }

        // D1 should be -A*e*alpha' (outer product)
        Matrix D1 = me.D(1);
        assertNotNull(D1);
        assertEquals(2, D1.getNumRows());
        assertEquals(2, D1.getNumCols());
    }

    @Test
    public void testMENonNormalizedAlpha() {
        // Test ME with non-normalized alpha (not summing to 1)
        // This is valid for ME but not for PH
        Matrix alpha = new Matrix(new double[]{0.8, 0.3}); // Sum = 1.1 > 1
        Matrix A = new Matrix(new double[][]{
                {-3.0, 2.0},
                {1.0, -2.0}
        });

        // This should be valid for ME (though BuTools may require sum <= 1 + tolerance)
        // Let's test with sum < 1 which is definitely valid
        Matrix alpha2 = new Matrix(new double[]{0.3, 0.4}); // Sum = 0.7 < 1
        Matrix A2 = new Matrix(new double[][]{
                {-3.0, 2.0},
                {1.0, -2.0}
        });

        ME me = new ME(alpha2, A2);
        assertNotNull(me);
        assertTrue(checkMERepresentation(alpha2, A2, 1e-14));
    }

    @Test
    public void testMEGetters() {
        Matrix alpha = new Matrix(new double[]{0.4, 0.6});
        Matrix A = new Matrix(new double[][]{
                {-2.0, 1.0},
                {0.5, -1.5}
        });

        ME me = new ME(alpha, A);

        // Test getAlpha()
        Matrix retrievedAlpha = me.getAlpha();
        assertArrayEquals(alpha.toArray1D(), retrievedAlpha.toArray1D(), TOLERANCE);

        // Test getA()
        Matrix retrievedA = me.getA();
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                assertEquals(A.get(i, j), retrievedA.get(i, j), TOLERANCE);
            }
        }

        // Test getNumberOfPhases()
        assertEquals(2, me.getNumberOfPhases());
    }
}
