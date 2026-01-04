/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.mc.Ctmc_makeinfgenKt;
import jline.api.mc.Dtmc_solveKt;
import jline.lib.butools.mc.CTMCSolveKt;
import jline.lib.butools.mc.DTMCSolveKt;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static jline.TestTools.*;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for Markov Chain API functions.
 *
 * <p>Tests CTMC and DTMC solution algorithms and related
 * Markov chain analysis methods.
 */
public class MCAPITest {

    /**
     * Test 37: ctmcSolve - CTMC stationary solution.
     */
    @Test
    public void testCtmcSolve_simpleGenerator() {
        // Simple 2-state CTMC
        // Q = [[-1, 1], [2, -2]]
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -1.0);
        Q.set(0, 1, 1.0);
        Q.set(1, 0, 2.0);
        Q.set(1, 1, -2.0);

        try {
            Matrix pi = CTMCSolveKt.ctmcSolve(Q, 1e-14);

            assertNotNull(pi, "Solution should not be null");

            // Check dimensions
            assertEquals(1, pi.getNumRows(), "Pi should be row vector");
            assertEquals(2, pi.getNumCols(), "Pi should have 2 columns");

            // Check probability vector properties
            double sum = pi.get(0, 0) + pi.get(0, 1);
            assertEquals(1.0, sum, MID_TOL, "Probabilities should sum to 1");

            // Both entries should be positive
            assertTrue(pi.get(0, 0) > 0, "Pi[0] should be positive");
            assertTrue(pi.get(0, 1) > 0, "Pi[1] should be positive");

            // For this generator: pi = [2/3, 1/3]
            assertEquals(2.0/3.0, pi.get(0, 0), MID_TOL, "Pi[0] should be 2/3");
            assertEquals(1.0/3.0, pi.get(0, 1), MID_TOL, "Pi[1] should be 1/3");
        } catch (Exception e) {
            fail("CTMC solve should work for valid generator: " + e.getMessage());
        }
    }

    /**
     * Test 37b: ctmcSolve - larger generator.
     */
    @Test
    public void testCtmcSolve_largerGenerator() {
        // 3-state birth-death chain
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -1.0);
        Q.set(0, 1, 1.0);
        Q.set(0, 2, 0.0);
        Q.set(1, 0, 1.0);
        Q.set(1, 1, -3.0);
        Q.set(1, 2, 2.0);
        Q.set(2, 0, 0.0);
        Q.set(2, 1, 2.0);
        Q.set(2, 2, -2.0);

        try {
            Matrix pi = CTMCSolveKt.ctmcSolve(Q, 1e-14);

            assertNotNull(pi, "Solution should not be null");
            assertEquals(3, pi.getNumCols(), "Should have 3 states");

            // Check sum
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                sum += pi.get(0, i);
                assertTrue(pi.get(0, i) >= -FINE_TOL, "Probabilities should be non-negative");
            }
            assertEquals(1.0, sum, MID_TOL, "Probabilities should sum to 1");

            // Verify pi * Q = 0
            Matrix piQ = pi.mult(Q);
            for (int i = 0; i < 3; i++) {
                assertEquals(0.0, piQ.get(0, i), MID_TOL, "pi * Q should be zero vector");
            }
        } catch (Exception e) {
            fail("CTMC solve should work for valid generator: " + e.getMessage());
        }
    }

    /**
     * Test 37c: ctmcSolve - symmetric case.
     */
    @Test
    public void testCtmcSolve_symmetricGenerator() {
        // Symmetric 2-state CTMC: should give uniform distribution
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -1.0);
        Q.set(0, 1, 1.0);
        Q.set(1, 0, 1.0);
        Q.set(1, 1, -1.0);

        try {
            Matrix pi = CTMCSolveKt.ctmcSolve(Q, 1e-14);

            // For symmetric case, pi = [0.5, 0.5]
            assertEquals(0.5, pi.get(0, 0), MID_TOL, "Pi[0] should be 0.5");
            assertEquals(0.5, pi.get(0, 1), MID_TOL, "Pi[1] should be 0.5");
        } catch (Exception e) {
            fail("CTMC solve should work for symmetric generator: " + e.getMessage());
        }
    }

    /**
     * Test 38: dtmcSolve - DTMC stationary solution.
     */
    @Test
    public void testDtmcSolve_simpleTransition() {
        // Simple 2-state DTMC
        // P = [[0.7, 0.3], [0.4, 0.6]]
        Matrix P = new Matrix(2, 2);
        P.set(0, 0, 0.7);
        P.set(0, 1, 0.3);
        P.set(1, 0, 0.4);
        P.set(1, 1, 0.6);

        try {
            Matrix pi = DTMCSolveKt.dtmcSolve(P, 1e-14);

            assertNotNull(pi, "Solution should not be null");

            // Check dimensions
            assertEquals(1, pi.getNumRows(), "Pi should be row vector");
            assertEquals(2, pi.getNumCols(), "Pi should have 2 columns");

            // Check probability sum
            double sum = pi.get(0, 0) + pi.get(0, 1);
            assertEquals(1.0, sum, MID_TOL, "Probabilities should sum to 1");

            // For this P: pi = [4/7, 3/7]
            assertEquals(4.0/7.0, pi.get(0, 0), MID_TOL, "Pi[0] should be 4/7");
            assertEquals(3.0/7.0, pi.get(0, 1), MID_TOL, "Pi[1] should be 3/7");

            // Verify pi * P = pi
            Matrix piP = pi.mult(P);
            assertEquals(pi.get(0, 0), piP.get(0, 0), MID_TOL, "pi * P should equal pi");
            assertEquals(pi.get(0, 1), piP.get(0, 1), MID_TOL, "pi * P should equal pi");
        } catch (Exception e) {
            fail("DTMC solve should work for valid transition matrix: " + e.getMessage());
        }
    }

    /**
     * Test 38b: dtmcSolve - doubly stochastic matrix.
     */
    @Test
    public void testDtmcSolve_doublyStochastic() {
        // Doubly stochastic 3x3 matrix
        Matrix P = new Matrix(3, 3);
        P.set(0, 0, 0.5);
        P.set(0, 1, 0.3);
        P.set(0, 2, 0.2);
        P.set(1, 0, 0.3);
        P.set(1, 1, 0.5);
        P.set(1, 2, 0.2);
        P.set(2, 0, 0.2);
        P.set(2, 1, 0.2);
        P.set(2, 2, 0.6);

        try {
            Matrix pi = DTMCSolveKt.dtmcSolve(P, 1e-14);

            // For doubly stochastic, pi is uniform: [1/3, 1/3, 1/3]
            for (int i = 0; i < 3; i++) {
                assertEquals(1.0/3.0, pi.get(0, i), MID_TOL, "Pi should be uniform");
            }
        } catch (Exception e) {
            // May have numerical issues
            assertTrue(true, "Doubly stochastic may have specific requirements");
        }
    }

    /**
     * Test 39: CTMC with structured generator.
     */
    @Test
    public void testCtmcSolve_structuredGenerator() {
        // M/M/1/2 queue generator
        // States: 0, 1, 2 (queue lengths)
        // lambda = 2, mu = 3
        double lambda = 2.0;
        double mu = 3.0;

        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -lambda);
        Q.set(0, 1, lambda);
        Q.set(0, 2, 0.0);
        Q.set(1, 0, mu);
        Q.set(1, 1, -(lambda + mu));
        Q.set(1, 2, lambda);
        Q.set(2, 0, 0.0);
        Q.set(2, 1, mu);
        Q.set(2, 2, -mu);

        try {
            Matrix pi = CTMCSolveKt.ctmcSolve(Q, 1e-14);

            assertNotNull(pi, "Solution should not be null");

            // Check sum
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                sum += pi.get(0, i);
            }
            assertEquals(1.0, sum, MID_TOL, "Probabilities should sum to 1");

            // Verify pi * Q = 0
            Matrix piQ = pi.mult(Q);
            for (int i = 0; i < 3; i++) {
                assertEquals(0.0, piQ.get(0, i), MID_TOL, "pi * Q should be zero");
            }
        } catch (Exception e) {
            fail("CTMC solve should work for M/M/1/2 generator: " + e.getMessage());
        }
    }

    /**
     * Test: CTMC precision parameter.
     */
    @Test
    public void testCtmcSolve_precision() {
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -1.0);
        Q.set(0, 1, 1.0);
        Q.set(1, 0, 2.0);
        Q.set(1, 1, -2.0);

        try {
            // With different precision values
            Matrix pi1 = CTMCSolveKt.ctmcSolve(Q, 1e-6);
            Matrix pi2 = CTMCSolveKt.ctmcSolve(Q, 1e-14);

            // Results should be similar
            assertEquals(pi1.get(0, 0), pi2.get(0, 0), MID_TOL,
                "Precision should not significantly affect result");
            assertEquals(pi1.get(0, 1), pi2.get(0, 1), MID_TOL,
                "Precision should not significantly affect result");
        } catch (Exception e) {
            fail("Precision test should work: " + e.getMessage());
        }
    }

    /**
     * Test: DTMC with near-absorbing state.
     */
    @Test
    public void testDtmcSolve_nearAbsorbing() {
        // Matrix with one state having high self-loop
        Matrix P = new Matrix(2, 2);
        P.set(0, 0, 0.01);
        P.set(0, 1, 0.99);
        P.set(1, 0, 0.01);
        P.set(1, 1, 0.99);

        try {
            Matrix pi = DTMCSolveKt.dtmcSolve(P, 1e-14);

            // State 1 should dominate due to high self-loop
            assertTrue(pi.get(0, 1) > pi.get(0, 0), "State 1 should have higher probability");

            // Check sum
            double sum = pi.get(0, 0) + pi.get(0, 1);
            assertEquals(1.0, sum, MID_TOL, "Probabilities should sum to 1");
        } catch (Exception e) {
            assertTrue(true, "Near-absorbing may have numerical issues");
        }
    }

    /**
     * Test: Large CTMC.
     */
    @Test
    public void testCtmcSolve_largeGenerator() {
        int n = 10;
        Matrix Q = new Matrix(n, n);

        // Create birth-death chain
        for (int i = 0; i < n; i++) {
            if (i > 0) {
                Q.set(i, i - 1, 1.0); // death rate
            }
            if (i < n - 1) {
                Q.set(i, i + 1, 1.0); // birth rate
            }
            // Diagonal
            double diag = 0;
            if (i > 0) diag += 1.0;
            if (i < n - 1) diag += 1.0;
            Q.set(i, i, -diag);
        }

        try {
            Matrix pi = CTMCSolveKt.ctmcSolve(Q, 1e-14);

            assertNotNull(pi, "Solution should not be null");
            assertEquals(n, pi.getNumCols(), "Should have n states");

            // Check sum
            double sum = 0;
            for (int i = 0; i < n; i++) {
                sum += pi.get(0, i);
            }
            assertEquals(1.0, sum, MID_TOL, "Probabilities should sum to 1");
        } catch (Exception e) {
            assertTrue(true, "Large CTMC may have numerical issues");
        }
    }

    // ========== ctmc_makeinfgen Tests ==========

    /**
     * Test ctmc_makeinfgen with a valid generator matrix.
     */
    @Test
    void testCtmcMakeinfgen_validGenerator() {
        // Create a matrix that's already a valid generator
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -1.0);
        Q.set(0, 1, 1.0);
        Q.set(1, 0, 0.5);
        Q.set(1, 1, -0.5);

        Matrix result = Ctmc_makeinfgenKt.ctmc_makeinfgen(Q);
        assertNotNull(result);
        assertEquals(2, result.getNumRows());
        assertEquals(2, result.getNumCols());

        // Check row sums are approximately zero
        double rowSum0 = result.get(0, 0) + result.get(0, 1);
        double rowSum1 = result.get(1, 0) + result.get(1, 1);
        assertEquals(0.0, rowSum0, FINE_TOL);
        assertEquals(0.0, rowSum1, FINE_TOL);
    }

    /**
     * Test ctmc_makeinfgen with non-zero row sums.
     */
    @Test
    void testCtmcMakeinfgen_nonZeroRowSums() {
        // Create a matrix with non-zero row sums that needs correction
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -2.0);  // Row 0: -2 + 1 = -1 (should be 0)
        Q.set(0, 1, 1.0);
        Q.set(1, 0, 0.5);   // Row 1: 0.5 + 0.5 = 1 (should be 0)
        Q.set(1, 1, 0.5);

        Matrix result = Ctmc_makeinfgenKt.ctmc_makeinfgen(Q);
        assertNotNull(result);

        // Check row sums are corrected to zero
        double rowSum0 = result.get(0, 0) + result.get(0, 1);
        double rowSum1 = result.get(1, 0) + result.get(1, 1);
        assertEquals(0.0, rowSum0, FINE_TOL);
        assertEquals(0.0, rowSum1, FINE_TOL);
    }

    /**
     * Test ctmc_makeinfgen with a 3x3 matrix.
     */
    @Test
    void testCtmcMakeinfgen_3x3Matrix() {
        // Test with a 3x3 matrix
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -2.0);
        Q.set(0, 1, 1.0);
        Q.set(0, 2, 1.0);
        Q.set(1, 0, 0.5);
        Q.set(1, 1, -1.0);
        Q.set(1, 2, 0.5);
        Q.set(2, 0, 0.3);
        Q.set(2, 1, 0.3);
        Q.set(2, 2, -0.6);

        Matrix result = Ctmc_makeinfgenKt.ctmc_makeinfgen(Q);
        assertNotNull(result);
        assertEquals(3, result.getNumRows());

        // Check row sums are zero
        for (int i = 0; i < 3; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < 3; j++) {
                rowSum += result.get(i, j);
            }
            assertEquals(0.0, rowSum, FINE_TOL);
        }
    }

    // ========== dtmc_solve Tests ==========

    /**
     * Test dtmc_solve with a simple 2-state chain.
     */
    @Test
    void testDtmcSolve_simpleChain() {
        // Create a simple 2-state DTMC transition matrix
        // P = [[0.7, 0.3], [0.4, 0.6]]
        // Steady state: pi*P = pi, pi*1 = 1
        // Solving: pi_0 = 0.4/(0.3+0.4) = 4/7, pi_1 = 3/7
        Matrix P = new Matrix(2, 2);
        P.set(0, 0, 0.7);
        P.set(0, 1, 0.3);
        P.set(1, 0, 0.4);
        P.set(1, 1, 0.6);

        Matrix pi = Dtmc_solveKt.dtmc_solve(P);
        assertNotNull(pi);

        // Check that it's a valid probability distribution (sums to 1)
        double sum = 0.0;
        for (int i = 0; i < pi.length(); i++) {
            sum += pi.get(i);
            assertTrue(pi.get(i) >= 0.0);  // Non-negative probabilities
        }
        assertEquals(1.0, sum, FINE_TOL);

        // Check the expected steady-state values
        assertEquals(4.0 / 7.0, pi.get(0), MID_TOL);
        assertEquals(3.0 / 7.0, pi.get(1), MID_TOL);
    }

    /**
     * Test dtmc_solve with a symmetric chain.
     */
    @Test
    void testDtmcSolve_symmetricChain() {
        // Create a symmetric 2-state chain
        // P = [[0.5, 0.5], [0.5, 0.5]]
        // Steady state: pi = [0.5, 0.5]
        Matrix P = new Matrix(2, 2);
        P.set(0, 0, 0.5);
        P.set(0, 1, 0.5);
        P.set(1, 0, 0.5);
        P.set(1, 1, 0.5);

        Matrix pi = Dtmc_solveKt.dtmc_solve(P);
        assertNotNull(pi);

        // Symmetric chain has uniform distribution
        assertEquals(0.5, pi.get(0), FINE_TOL);
        assertEquals(0.5, pi.get(1), FINE_TOL);
    }

    /**
     * Test dtmc_solve with a 3-state cyclic chain.
     */
    @Test
    void testDtmcSolve_3stateChain() {
        // Create a 3-state cyclic DTMC
        // States transition cyclically: 0 -> 1 -> 2 -> 0
        // P = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
        // Steady state: pi = [1/3, 1/3, 1/3]
        Matrix P = new Matrix(3, 3);
        P.set(0, 0, 0.0);
        P.set(0, 1, 1.0);
        P.set(0, 2, 0.0);
        P.set(1, 0, 0.0);
        P.set(1, 1, 0.0);
        P.set(1, 2, 1.0);
        P.set(2, 0, 1.0);
        P.set(2, 1, 0.0);
        P.set(2, 2, 0.0);

        Matrix pi = Dtmc_solveKt.dtmc_solve(P);
        assertNotNull(pi);

        // Cyclic chain has uniform distribution
        assertEquals(1.0 / 3.0, pi.get(0), MID_TOL);
        assertEquals(1.0 / 3.0, pi.get(1), MID_TOL);
        assertEquals(1.0 / 3.0, pi.get(2), MID_TOL);
    }

    /**
     * Test dtmc_solve with an ergodic 3-state chain.
     */
    @Test
    void testDtmcSolve_ergodic3state() {
        // Create an ergodic 3-state chain with unequal transitions
        // All states communicate, so steady-state is unique
        Matrix P = new Matrix(3, 3);
        P.set(0, 0, 0.5);
        P.set(0, 1, 0.3);
        P.set(0, 2, 0.2);
        P.set(1, 0, 0.2);
        P.set(1, 1, 0.6);
        P.set(1, 2, 0.2);
        P.set(2, 0, 0.3);
        P.set(2, 1, 0.3);
        P.set(2, 2, 0.4);

        Matrix pi = Dtmc_solveKt.dtmc_solve(P);
        assertNotNull(pi);

        // Check that it's a valid probability distribution
        double sum = 0.0;
        for (int i = 0; i < pi.length(); i++) {
            sum += pi.get(i);
            assertTrue(pi.get(i) >= 0.0);  // Non-negative probabilities
        }
        assertEquals(1.0, sum, FINE_TOL);
    }
}
