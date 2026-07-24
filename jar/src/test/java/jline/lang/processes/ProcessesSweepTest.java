/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Moment/LST/sampling sweep across the distribution and arrival-process
 * classes of jline.lang.processes. For every distribution:
 * - moments are finite and mutually consistent (Var = SCV * mean^2);
 * - the LST is a completely monotone transform: LST(0)=1, decreasing in s;
 * - the empirical mean of a large sample matches the analytic mean.
 */
public class ProcessesSweepTest {

    private static final int SAMPLE_SIZE = 20000;
    private static final long SAMPLE_SEED = 23000L;
    private static final double SAMPLE_RTOL = 0.10;
    private static final double MOMENT_RTOL = 1e-8;

    /** Distribution under test plus which checks are applicable to it. */
    public static class Case {
        final String name;
        final Distribution dist;
        final boolean hasLst;
        final boolean sampleable;

        Case(String name, Distribution dist, boolean hasLst, boolean sampleable) {
            this.name = name;
            this.dist = dist;
            this.hasLst = hasLst;
            this.sampleable = sampleable;
        }

        @Override
        public String toString() {
            return name;
        }
    }

    @BeforeAll
    public static void setUpClass() {
        Maths.setRandomNumbersMatlab(true);
    }

    public static List<Case> distributions() {
        List<Case> cases = new ArrayList<Case>();
        // Continuous distributions
        cases.add(new Case("Exp", new Exp(2.0), true, true));
        cases.add(new Case("Erlang", new Erlang(4.0, 2), true, true));
        cases.add(new Case("HyperExp", new HyperExp(0.3, 4.0, 1.5), true, true));
        cases.add(new Case("Cox2", new Cox2(3.0, 1.5, 0.4), true, true));
        List<Double> coxMu = Arrays.asList(3.0, 1.5);
        List<Double> coxPhi = Arrays.asList(0.4, 1.0);
        cases.add(new Case("Coxian", new Coxian(coxMu, coxPhi), true, true));
        cases.add(new Case("Det", new Det(0.5), true, true));
        cases.add(new Case("Uniform", new Uniform(0.2, 0.8), true, true));
        cases.add(new Case("Gamma", new Gamma(2.0, 0.25), true, true));
        // hasLst was false, which skipped lstIsValidTransform for Pareto entirely --
        // that is why evalLST returning A*(0)=0.98558 for this very case went
        // unnoticed. The transform is exact now (2.9e-12 worst case), so the LST
        // assertions are enabled. Weibull/Lognormal/Normal stay opted out: their
        // evalLST is still the inherited numeric fallback, not a checked transform.
        cases.add(new Case("Pareto", new Pareto(3.0, 0.334), true, true));
        cases.add(new Case("Weibull", new Weibull(2.0, 0.5), false, true));
        cases.add(new Case("Lognormal", new Lognormal(-1.0, 0.5), false, true));
        cases.add(new Case("Normal", new Normal(5.0, 1.0), false, true));

        // Phase-type family
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 1.0);
        Matrix T = new Matrix(2, 2);
        T.set(0, 0, -3.0);
        T.set(0, 1, 3.0);
        T.set(1, 1, -2.0);
        cases.add(new Case("PH", new PH(alpha.copy(), T.copy()), true, true));
        cases.add(new Case("APH", new APH(alpha.copy(), T.copy()), true, true));

        // Markov arrival processes
        Matrix d0 = new Matrix(2, 2);
        d0.set(0, 0, -4.0);
        d0.set(0, 1, 1.0);
        d0.set(1, 0, 0.5);
        d0.set(1, 1, -2.0);
        Matrix d1 = new Matrix(2, 2);
        d1.set(0, 0, 2.0);
        d1.set(0, 1, 1.0);
        d1.set(1, 0, 1.0);
        d1.set(1, 1, 0.5);
        cases.add(new Case("MAP", new MAP(d0.copy(), d1.copy()), false, true));
        cases.add(new Case("MMPP2", new MMPP2(4.0, 1.0, 0.5, 0.2), false, true));

        // Discrete distributions
        cases.add(new Case("Bernoulli", new Bernoulli(0.3), false, true));
        cases.add(new Case("Binomial", new Binomial(10, 0.3), false, true));
        cases.add(new Case("Geometric", new Geometric(0.4), false, true));
        cases.add(new Case("Poisson", new Poisson(3.0), false, true));
        cases.add(new Case("DiscreteUniform", new DiscreteUniform(1.0, 6.0), false, true));
        return cases;
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("distributions")
    public void momentsAreConsistent(Case c) {
        double mean = c.dist.getMean();
        double scv = c.dist.getSCV();
        double var = c.dist.getVar();
        assertTrue(Double.isFinite(mean), c.name + ": mean not finite");
        assertTrue(mean > 0 || "Normal".equals(c.name) || "Bernoulli".equals(c.name),
                c.name + ": mean must be positive, got " + mean);
        assertTrue(Double.isFinite(scv) && scv >= 0, c.name + ": invalid SCV " + scv);
        assertTrue(Double.isFinite(var) && var >= 0, c.name + ": invalid variance " + var);
        if (mean != 0) {
            assertEquals(scv * mean * mean, var,
                    MOMENT_RTOL * Math.max(1.0, Math.abs(var)) + 1e-12,
                    c.name + ": Var and SCV*mean^2 disagree");
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("distributions")
    public void lstIsValidTransform(Case c) {
        if (!c.hasLst) {
            return;
        }
        double lst0 = c.dist.evalLST(0.0);
        assertEquals(1.0, lst0, 1e-8, c.name + ": LST(0) must equal 1");
        double prev = lst0;
        for (double s : new double[]{0.1, 0.5, 1.0, 2.0, 5.0}) {
            double v = c.dist.evalLST(s);
            assertTrue(v > 0 && v <= 1.0 + 1e-12,
                    c.name + ": LST(" + s + ") outside (0,1]: " + v);
            assertTrue(v <= prev + 1e-12,
                    c.name + ": LST must be nonincreasing in s");
            prev = v;
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("distributions")
    public void sampleMeanMatchesAnalyticMean(Case c) {
        if (!c.sampleable) {
            return;
        }
        double[] samples = c.dist.sample(SAMPLE_SIZE, new Random(SAMPLE_SEED));
        assertNotNull(samples, c.name + ": sample returned null");
        assertEquals(SAMPLE_SIZE, samples.length, c.name + ": wrong sample count");
        double sum = 0;
        for (double x : samples) {
            assertTrue(Double.isFinite(x), c.name + ": non-finite sample " + x);
            sum += x;
        }
        double empirical = sum / SAMPLE_SIZE;
        double mean = c.dist.getMean();
        double tol = Math.max(SAMPLE_RTOL * Math.abs(mean), 5.0
                * Math.sqrt(c.dist.getVar() / SAMPLE_SIZE));
        assertEquals(mean, empirical, tol,
                c.name + ": empirical mean " + empirical + " vs analytic " + mean);
    }
}
