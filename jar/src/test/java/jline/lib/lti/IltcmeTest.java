/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import java.util.function.UnaryOperator;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

/**
 * Tests for the unified iltcme inverse Laplace transform.
 * Validates CME, Euler, and Gaver methods against known analytical inversions.
 */
public class IltcmeTest {

    // F(s) = 1/(1+s) -> f(t) = exp(-t)
    private static final UnaryOperator<Complex> EXP_LT = s -> Complex.ONE.divide(Complex.ONE.add(s));

    // F(s) = 1/(1+s^2) -> f(t) = sin(t)
    private static final UnaryOperator<Complex> SIN_LT = s -> Complex.ONE.divide(Complex.ONE.add(s.multiply(s)));

    private static final double[] T_POINTS = {0.5, 1.0, 2.0, 5.0};
    private static final int MAX_FN_EVALS = 25;

    @Test
    public void testCmeExponential() {
        double[] result = iltcme.ilt(EXP_LT, T_POINTS, MAX_FN_EVALS, "cme");
        assertEquals(T_POINTS.length, result.length);
        for (int i = 0; i < T_POINTS.length; i++) {
            double expected = Math.exp(-T_POINTS[i]);
            assertEquals(expected, result[i], LOOSE_MID_TOL,
                    String.format("CME exp(-t) at t=%.1f: expected=%.10f, got=%.10f", T_POINTS[i], expected, result[i]));
        }
    }

    @Test
    public void testCmeSine() {
        double[] result = iltcme.ilt(SIN_LT, T_POINTS, MAX_FN_EVALS, "cme");
        assertEquals(T_POINTS.length, result.length);
        for (int i = 0; i < T_POINTS.length; i++) {
            double expected = Math.sin(T_POINTS[i]);
            assertEquals(expected, result[i], COARSE_TOL,
                    String.format("CME sin(t) at t=%.1f: expected=%.10f, got=%.10f", T_POINTS[i], expected, result[i]));
        }
    }

    @Test
    public void testEulerExponential() {
        double[] result = iltcme.ilt(EXP_LT, T_POINTS, MAX_FN_EVALS, "euler");
        assertEquals(T_POINTS.length, result.length);
        for (int i = 0; i < T_POINTS.length; i++) {
            double expected = Math.exp(-T_POINTS[i]);
            assertEquals(expected, result[i], LOOSE_FINE_TOL,
                    String.format("Euler exp(-t) at t=%.1f: expected=%.10f, got=%.10f", T_POINTS[i], expected, result[i]));
        }
    }

    @Test
    public void testEulerSine() {
        double[] result = iltcme.ilt(SIN_LT, T_POINTS, MAX_FN_EVALS, "euler");
        assertEquals(T_POINTS.length, result.length);
        for (int i = 0; i < T_POINTS.length; i++) {
            double expected = Math.sin(T_POINTS[i]);
            assertEquals(expected, result[i], LOOSE_FINE_TOL,
                    String.format("Euler sin(t) at t=%.1f: expected=%.10f, got=%.10f", T_POINTS[i], expected, result[i]));
        }
    }

    @Test
    public void testGaverExponential() {
        double[] result = iltcme.ilt(EXP_LT, T_POINTS, MAX_FN_EVALS, "gaver");
        assertEquals(T_POINTS.length, result.length);
        for (int i = 0; i < T_POINTS.length; i++) {
            double expected = Math.exp(-T_POINTS[i]);
            assertEquals(expected, result[i], 0.15,
                    String.format("Gaver exp(-t) at t=%.1f: expected=%.10f, got=%.10f", T_POINTS[i], expected, result[i]));
        }
    }

    @Test
    public void testGaverSine() {
        double[] result = iltcme.ilt(SIN_LT, T_POINTS, MAX_FN_EVALS, "gaver");
        assertEquals(T_POINTS.length, result.length);
        for (int i = 0; i < T_POINTS.length; i++) {
            double expected = Math.sin(T_POINTS[i]);
            assertEquals(expected, result[i], LOOSE_COARSE_TOL,
                    String.format("Gaver sin(t) at t=%.1f: expected=%.10f, got=%.10f", T_POINTS[i], expected, result[i]));
        }
    }

    @Test
    public void testCmeParameterLoading() {
        // Verify that CME parameters load correctly and produce reasonable results
        // With maxFnEvals=3, should use n=1 or n=2 (smallest that fits)
        double[] result = iltcme.ilt(EXP_LT, new double[]{1.0}, 3, "cme");
        assertEquals(1, result.length);
        assertTrue(Math.abs(result[0] - Math.exp(-1.0)) < 0.05,
                "CME with maxFnEvals=3 should still give reasonable result for exp(-t)");
    }

    @Test
    public void testMaxFnEvalsConstraint() {
        // With very large maxFnEvals, CME should select a high-order (low cv2) entry
        double[] resultLow = iltcme.ilt(EXP_LT, new double[]{1.0}, 5, "cme");
        double[] resultHigh = iltcme.ilt(EXP_LT, new double[]{1.0}, 100, "cme");
        // Higher order should give better or equal accuracy
        double errLow = Math.abs(resultLow[0] - Math.exp(-1.0));
        double errHigh = Math.abs(resultHigh[0] - Math.exp(-1.0));
        assertTrue(errHigh <= errLow + 1e-10,
                String.format("Higher maxFnEvals should give better accuracy: err5=%.2e, err100=%.2e", errLow, errHigh));
    }

    @Test
    public void testDefaultMethodIsCme() {
        // Calling without method argument should use CME
        double[] result = iltcme.ilt(EXP_LT, T_POINTS, MAX_FN_EVALS);
        double[] resultCme = iltcme.ilt(EXP_LT, T_POINTS, MAX_FN_EVALS, "cme");
        assertArrayEquals(resultCme, result, ZERO_TOL, "Default method should be CME");
    }

    @Test
    public void testInvalidMethod() {
        assertThrows(IllegalArgumentException.class, () -> {
            iltcme.ilt(EXP_LT, T_POINTS, MAX_FN_EVALS, "invalid");
        });
    }

}
