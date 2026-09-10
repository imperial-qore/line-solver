/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.function.Executable;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for {@link DistributionScaling#scaleRate(Distribution, double)},
 * the perturbation primitive of the finite-difference branch of
 * {@link jline.solvers.NetworkSolver#getSensitivityTable(String, double, String)}.
 *
 * <p>The contract the finite-difference branch relies on is that the scaling is a
 * pure time scaling X -> X/factor: the mean is divided by the factor and every
 * dimensionless shape statistic, the SCV in particular, is left alone. If the SCV
 * moved, the difference quotient would measure the derivative along a mixed
 * rate-and-variability direction and would not be d(.)/d(rate) at all, so the SCV
 * assertion here is the load-bearing one, not the mean.</p>
 *
 * <p>Java counterpart of MATLAB's {@code dist_scale_rate.m}.</p>
 */
public class DistributionScalingTest {

    /** Relative agreement required; the scaling is exact, so this is near machine eps. */
    private static final double REL_TOL = 1e-9;

    /**
     * Asserts the two defining properties on one distribution: the mean divides by the
     * factor and the SCV is untouched.
     */
    private static void assertPureRateScaling(Distribution d, double factor) {
        double mean0 = d.getMean();
        double scv0 = d.getSCV();
        Distribution s = DistributionScaling.scaleRate(d, factor);
        assertNotNull(s, "scaleRate must return a distribution");
        assertEquals(d.getClass(), s.getClass(), "the scaled copy must stay in the same family");
        String tag = d.getClass().getSimpleName();
        assertRel(mean0 / factor, s.getMean(), tag + " mean");
        assertRel(scv0, s.getSCV(), tag + " SCV");
        // The original must be left untouched: the sensitivity loop restores it after
        // each perturbation and would otherwise restore a corrupted process.
        assertRel(mean0, d.getMean(), tag + " original mean");
        assertRel(scv0, d.getSCV(), tag + " original SCV");
    }

    private static void assertRel(double expected, double actual, String msg) {
        // Relative to 1e-9, with a machine-scale absolute floor so that a statistic
        // that is exactly zero (the SCV of a Det process) is still testable.
        assertTrue(Math.abs(actual - expected) <= REL_TOL * Math.abs(expected) + 1e-12,
                msg + ": expected " + expected + " but got " + actual
                        + " (err " + Math.abs(actual - expected) + ")");
    }

    // ---------- per-family scaling ----------------------------------------

    @Test
    @DisplayName("Exp, Erlang and HyperExp scale as pure rate changes")
    public void testExponentialFamilies() {
        assertPureRateScaling(new Exp(2.5), 1.0001);
        assertPureRateScaling(new Exp(2.5), 0.7);
        assertPureRateScaling(new Erlang(3.0, 4), 1.3);
        assertPureRateScaling(new HyperExp(0.3, 1.5, 6.0), 1.7);
        assertPureRateScaling(new HyperExp(new double[]{0.2, 0.5, 0.3},
                new double[]{1.0, 4.0, 9.0}), 0.55);
    }

    @Test
    @DisplayName("Coxian and Cox2 scale the phase rates and keep the completion probabilities")
    public void testCoxianFamilies() {
        assertPureRateScaling(Coxian.fitMeanAndSCV(0.4, 2.5), 1.9);
        assertPureRateScaling(new Cox2(3.0, 7.0, 0.25), 0.6);
        // The completion probabilities are dimensionless and must be carried over.
        Cox2 c = new Cox2(3.0, 7.0, 0.25);
        Cox2 s = (Cox2) DistributionScaling.scaleRate(c, 2.0);
        assertRel(c.getPhi().get(0), s.getPhi().get(0), "Cox2 phi(1)");
        assertRel(2.0 * c.getMu().get(0), s.getMu().get(0), "Cox2 mu(1)");
        assertRel(2.0 * c.getMu().get(1), s.getMu().get(1), "Cox2 mu(2)");
    }

    @Test
    @DisplayName("APH, PH, MAP and MMPP2 scale their Markovian representations")
    public void testMarkovianFamilies() {
        assertPureRateScaling(APH.fitMeanAndSCV(0.8, 3.0), 1.45);
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 0.35);
        alpha.set(0, 1, 0.65);
        Matrix T = new Matrix(2, 2);
        T.set(0, 0, -5.0);
        T.set(0, 1, 2.0);
        T.set(1, 0, 1.0);
        T.set(1, 1, -4.0);
        assertPureRateScaling(new PH(alpha, T), 1.25);
        assertPureRateScaling(MAP.rand(3), 0.8);
        assertPureRateScaling(new MMPP2(4.0, 1.0, 0.5, 0.3), 1.6);
    }

    @Test
    @DisplayName("the non-Markovian families scale through their time parameters")
    public void testContinuousFamilies() {
        assertPureRateScaling(new Det(0.75), 2.0);
        assertPureRateScaling(new Uniform(0.5, 2.5), 1.4);
        assertPureRateScaling(new Gamma(2.5, 0.4), 1.15);
        assertPureRateScaling(new Pareto(3.0, 0.5), 0.65);
        assertPureRateScaling(new Weibull(1.7, 0.9), 1.35);
        assertPureRateScaling(new Lognormal(-0.5, 0.7), 1.8);
    }

    @Test
    @DisplayName("a unit factor is the identity on every statistic")
    public void testUnitFactorIsIdentity() {
        assertPureRateScaling(new Exp(2.0), 1.0);
        assertPureRateScaling(new Erlang(3.0, 2), 1.0);
        assertPureRateScaling(new Gamma(2.0, 0.5), 1.0);
    }

    // ---------- error contract --------------------------------------------

    @Test
    @DisplayName("an unsupported family is rejected, and the message names it")
    public void testUnsupportedFamilyThrows() {
        RuntimeException e = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                DistributionScaling.scaleRate(new Disabled(), 1.5);
            }
        });
        assertTrue(e.getMessage().contains("Disabled"),
                "the message must name the offending class: " + e.getMessage());
        assertTrue(e.getMessage().contains("Rate scaling is not defined"),
                "unexpected message: " + e.getMessage());
    }

    @Test
    @DisplayName("a non-positive or non-finite factor is rejected")
    public void testInvalidFactorThrows() {
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                DistributionScaling.scaleRate(new Exp(1.0), 0.0);
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                DistributionScaling.scaleRate(new Exp(1.0), -1.0);
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                DistributionScaling.scaleRate(new Exp(1.0), Double.NaN);
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                DistributionScaling.scaleRate(new Exp(1.0), Double.POSITIVE_INFINITY);
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                DistributionScaling.scaleRate(null, 1.5);
            }
        });
    }

    @Test
    @DisplayName("an Immediate process is scale-invariant")
    public void testImmediateIsScaleInvariant() {
        Distribution s = DistributionScaling.scaleRate(new Immediate(), 3.0);
        assertEquals(Immediate.class, s.getClass());
        assertEquals(0.0, s.getMean(), 1e-12, "an immediate process has zero mean at any scale");
    }

    @Test
    @DisplayName("an NHPP schedule is time-scaled, rates up and breakpoints compressed")
    public void testNhppScheduleIsTimeScaled() {
        double factor = 2.5;
        NHPP nhpp = new NHPP(new double[]{0.0, 1.0, 2.0}, new double[]{2.0, 4.0}, true);
        double mean0 = nhpp.getMean();
        Distribution scaled = DistributionScaling.scaleRate(nhpp, factor);
        assertEquals(NHPP.class, scaled.getClass());
        NHPP s = (NHPP) scaled;
        assertRel(mean0 / factor, s.getMean(), "NHPP mean");
        double[] breakpoints = s.getBreakpoints();
        double[] rates = s.getRates();
        assertRel(0.0 / factor, breakpoints[0], "NHPP breakpoint 0");
        assertRel(1.0 / factor, breakpoints[1], "NHPP breakpoint 1");
        assertRel(2.0 / factor, breakpoints[2], "NHPP breakpoint 2");
        assertRel(2.0 * factor, rates[0], "NHPP rate 0");
        assertRel(4.0 * factor, rates[1], "NHPP rate 1");
        assertTrue(s.isCyclic(), "the cyclic flag must carry over");
    }

    @Test
    @DisplayName("a Replayer trace is scaled sample by sample")
    public void testReplayerTraceIsScaled() {
        double factor = 2.0;
        Replayer replayer = new Replayer(new double[]{1.0, 2.0, 3.0, 6.0});
        double mean0 = replayer.getMean();
        double scv0 = replayer.getSCV();
        Distribution scaled = DistributionScaling.scaleRate(replayer, factor);
        assertEquals(Replayer.class, scaled.getClass());
        double[] data = ((Replayer) scaled).getData();
        assertRel(1.0 / factor, data[0], "sample 0");
        assertRel(6.0 / factor, data[3], "sample 3");
        assertRel(mean0 / factor, scaled.getMean(), "Replayer mean");
        assertRel(scv0, scaled.getSCV(), "Replayer SCV");
    }

    @Test
    @DisplayName("a file-backed Replayer is refused rather than rewritten")
    public void testFileBackedReplayerIsRefused() throws java.io.IOException {
        java.io.File trace = java.io.File.createTempFile("line_scaling_", ".trace");
        trace.deleteOnExit();
        java.io.PrintWriter out = new java.io.PrintWriter(trace);
        try {
            out.println("1.0");
            out.println("2.0");
            out.println("3.0");
        } finally {
            out.close();
        }
        final Replayer replayer = new Replayer(trace.getAbsolutePath());
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                DistributionScaling.scaleRate(replayer, 2.0);
            }
        });
    }

    @Test
    @DisplayName("scaling twice composes, and the inverse factor recovers the original")
    public void testCompositionAndInverse() {
        Exp d = new Exp(2.0);
        Distribution once = DistributionScaling.scaleRate(d, 1.5);
        Distribution twice = DistributionScaling.scaleRate(once, 1.0 / 1.5);
        assertRel(d.getMean(), twice.getMean(), "round-trip mean");
        assertRel(d.getSCV(), twice.getSCV(), "round-trip SCV");
        assertSame(Exp.class, twice.getClass());
    }
}
