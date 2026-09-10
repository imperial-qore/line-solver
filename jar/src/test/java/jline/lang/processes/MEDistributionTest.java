/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ldes.SolverLDES;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static jline.lib.butools.ph.CheckMERepresentation.checkMERepresentation;
import static jline.TestTools.*;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for ME (Matrix Exponential) distribution
 */
public class MEDistributionTest {


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
        assertArrayEquals(new double[]{0.4, 0.6}, me.getAlpha().toArray1D(), LOOSE_FINE_TOL);
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
        assertEquals(1.0 / rate, mean, LOOSE_FINE_TOL);

        // SCV should be 1 (for exponential)
        double scv = me.getSCV();
        assertEquals(1.0, scv, LOOSE_FINE_TOL);
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
        assertEquals((double) k / rate, mean, LOOSE_FINE_TOL);

        // SCV should be 1/k (for Erlang)
        double scv = me.getSCV();
        assertEquals(1.0 / k, scv, LOOSE_FINE_TOL);
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
        assertEquals(expectedMean, mean, LOOSE_FINE_TOL);
    }

    @Test
    public void testMEMeanComputation() {
        // Create ME with known mean
        Matrix alpha = new Matrix(new double[]{1.0});
        Matrix A = new Matrix(new double[][]{{-2.0}});

        ME me = new ME(alpha, A);

        // Mean = -alpha * A^(-1) * e = -1 * (-1/2) * 1 = 0.5
        double mean = me.getMean();
        assertEquals(0.5, mean, LOOSE_FINE_TOL);
    }

    @Test
    public void testMEVarianceComputation() {
        // Exponential with rate 2.0
        ME me = ME.fromExp(2.0);

        // For Exp(lambda): variance = 1/lambda^2 = 1/4 = 0.25
        double variance = me.getVar();
        assertEquals(0.25, variance, LOOSE_FINE_TOL);
    }

    @Test
    public void testMESCVComputation() {
        // Erlang-2 with rate 1.0
        ME me = ME.fromErlang(2, 1.0);

        // For Erlang-k: SCV = 1/k = 1/2 = 0.5
        double scv = me.getSCV();
        assertEquals(0.5, scv, LOOSE_FINE_TOL);
    }

    @Test
    public void testMECDFEvaluation() {
        // Exponential with rate 1.0
        ME me = ME.fromExp(1.0);

        // For Exp(1): CDF(t) = 1 - exp(-t)
        double cdf0 = me.evalCDF(0.0);
        assertEquals(0.0, cdf0, LOOSE_FINE_TOL);

        double cdf1 = me.evalCDF(1.0);
        assertEquals(1.0 - Math.exp(-1.0), cdf1, LOOSE_FINE_TOL);

        double cdf2 = me.evalCDF(2.0);
        assertEquals(1.0 - Math.exp(-2.0), cdf2, LOOSE_FINE_TOL);
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
        assertEquals(expectedMean, empiricalMean, VERY_COARSE_TOL * expectedMean);
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
                assertEquals(A.get(i, j), D0.get(i, j), LOOSE_FINE_TOL);
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
        assertArrayEquals(alpha.toArray1D(), retrievedAlpha.toArray1D(), LOOSE_FINE_TOL);

        // Test getA()
        Matrix retrievedA = me.getA();
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                assertEquals(A.get(i, j), retrievedA.get(i, j), LOOSE_FINE_TOL);
            }
        }

        // Test getNumberOfPhases()
        assertEquals(2, me.getNumberOfPhases());
    }

    // ------------------------------------------------------------------
    // Sampling tests. ME sampling inverts the exact CDF; the CTMC walk of
    // Map_sample is only valid on a phase-type (alpha, A) pair.
    // ------------------------------------------------------------------

    /** Number of draws used by the Monte Carlo sampling tests. */
    private static final int SAMPLE_SIZE = 200000;

    /** Empirical mean of a sample. */
    private static double sampleMean(double[] x) {
        double s = 0.0;
        for (int i = 0; i < x.length; i++) {
            s += x[i];
        }
        return s / x.length;
    }

    /** Empirical k-th raw moment of a sample. */
    private static double sampleMoment(double[] x, int k) {
        double s = 0.0;
        for (int i = 0; i < x.length; i++) {
            s += Math.pow(x[i], k);
        }
        return s / x.length;
    }

    /** Empirical squared coefficient of variation of a sample. */
    private static double sampleSCV(double[] x) {
        double m1 = sampleMean(x);
        double m2 = sampleMoment(x, 2);
        return (m2 - m1 * m1) / (m1 * m1);
    }

    /**
     * A matrix-exponential distribution that is not phase-type.
     *
     * <p>Representation: {@code A} is the real Jordan form with eigenvalues
     * {-1, -2 + 2i, -2 - 2i}, so the dominant eigenvalue -1 is real and simple
     * and {@code checkMERepresentation} accepts the pair. The initial vector is
     * {@code alpha = [0.6105899..., -0.1554714..., 0.5448815...]}, which sums to
     * one but has a negative second entry, so it is not a phase-type starting
     * distribution and {@code A} is not a sub-generator (it has the negative
     * off-diagonal entry A(2,1) = -2).</p>
     *
     * <p>How we know the distribution itself is not phase-type (not merely that
     * this particular representation is not a PH representation): the parameters
     * were solved so that the density
     * {@code f(t) = a e^-t + e^-2t ((2P - 2Q) cos 2t + (2Q + 2P) sin 2t)},
     * with a = alpha(1) and P, Q the cosine/sine coefficients implied by alpha,
     * has a double root at t0 = 1.2, i.e. f(1.2) = f'(1.2) = 0 while f >= 0
     * everywhere else. By O'Cinneide's characterisation a distribution is
     * phase-type if and only if its density is strictly positive on (0, infinity)
     * and its dominant pole is real and simple. The density here vanishes at an
     * interior point, so no phase-type representation of any order exists.
     * {@link #testMENonPhaseTypeDensityHasInteriorZero()} checks that root
     * numerically.</p>
     */
    private static ME nonPhaseTypeME() {
        Matrix alpha = new Matrix(new double[]{
                0.61058991931158258, -0.15547146730086722, 0.54488154798928464});
        Matrix A = new Matrix(new double[][]{
                {-1.0, 0.0, 0.0},
                {0.0, -2.0, 2.0},
                {0.0, -2.0, -2.0}
        });
        return new ME(alpha, A);
    }

    @Test
    public void testMESamplingMomentsExp() {
        ME me = ME.fromExp(2.0);
        double[] x = me.sample(SAMPLE_SIZE, new Random(7919));

        assertEquals(SAMPLE_SIZE, x.length);
        assertEquals(me.getMean(), sampleMean(x), 0.02 * me.getMean());
        assertEquals(me.getSCV(), sampleSCV(x), 0.05 * me.getSCV());
    }

    @Test
    public void testMESamplingMomentsErlang() {
        ME me = ME.fromErlang(3, 2.0);
        double[] x = me.sample(SAMPLE_SIZE, new Random(7919));

        assertEquals(me.getMean(), sampleMean(x), 0.02 * me.getMean());
        assertEquals(me.getSCV(), sampleSCV(x), 0.05 * me.getSCV());
    }

    @Test
    public void testMESamplingMomentsHyperExp() {
        ME me = ME.fromHyperExp(new double[]{0.3, 0.7}, new double[]{1.0, 4.0});
        double[] x = me.sample(SAMPLE_SIZE, new Random(7919));

        assertEquals(me.getMean(), sampleMean(x), 0.02 * me.getMean());
        assertEquals(me.getSCV(), sampleSCV(x), 0.05 * me.getSCV());
    }

    @Test
    public void testMENonPhaseTypeRepresentationIsAccepted() {
        ME me = nonPhaseTypeME();

        // The representation is valid as an ME but is not a PH representation:
        // alpha has a negative entry and A has a negative off-diagonal entry.
        assertTrue(checkMERepresentation(me.getAlpha(), me.getA(), 1e-14));
        assertTrue(me.getAlpha().get(0, 1) < 0.0, "alpha must have a negative entry");
        assertTrue(me.getA().get(2, 1) < 0.0, "A must have a negative off-diagonal entry");
        assertEquals(1.0, me.getAlpha().elementSum(), 1e-12);
    }

    @Test
    public void testMENonPhaseTypeDensityHasInteriorZero() {
        ME me = nonPhaseTypeME();

        // f(t) = dF/dt evaluated by central differences on the exact CDF. The
        // density is non-negative everywhere and touches zero at t0 = 1.2, which
        // is what rules out any phase-type representation (see nonPhaseTypeME).
        // The argmin is taken over t in [0.2, 3] because beyond that the density
        // decays monotonically to zero and the global minimum is the tail rather
        // than the interior root of interest.
        double dt = 1e-4;
        double fMin = Double.POSITIVE_INFINITY;
        double tMin = -1.0;
        for (int i = 1; i <= 40000; i++) {
            double t = i * 1e-3;
            double f = (me.evalCDF(t + dt) - me.evalCDF(t - dt)) / (2 * dt);
            assertTrue(f > -1e-6, "density must be non-negative, f(" + t + ") = " + f);
            if (t >= 0.2 && t <= 3.0 && f < fMin) {
                fMin = f;
                tMin = t;
            }
        }
        assertEquals(1.2, tMin, 1e-9);
        assertTrue(fMin < 1e-5, "density minimum should be zero, got " + fMin);
        // The zero is a genuine interior dip, not a flat tail: the density is
        // two orders of magnitude larger on either side of t0.
        double fLeft = (me.evalCDF(0.6 + dt) - me.evalCDF(0.6 - dt)) / (2 * dt);
        double fRight = (me.evalCDF(2.0 + dt) - me.evalCDF(2.0 - dt)) / (2 * dt);
        assertTrue(fLeft > 1e-2, "f(0.6) = " + fLeft);
        assertTrue(fRight > 1e-2, "f(2.0) = " + fRight);
    }

    @Test
    public void testMENonPhaseTypeSamplingMoments() {
        ME me = nonPhaseTypeME();

        // Analytic moments m_k = k! * alpha * (-A)^-k * e.
        double m1 = 0.532854185661149;
        double m2 = 1.0460915848006274;
        assertEquals(m1, me.getMean(), LOOSE_FINE_TOL);

        double[] x = me.sample(SAMPLE_SIZE, new Random(20260719));
        for (int i = 0; i < x.length; i++) {
            assertTrue(x[i] >= 0.0, "Sample should be non-negative");
        }
        assertEquals(m1, sampleMean(x), 0.02 * m1);
        assertEquals(m2, sampleMoment(x, 2), 0.05 * m2);
    }

    @Test
    public void testMENonPhaseTypeCtmcWalkIsWrong() {
        ME me = nonPhaseTypeME();

        // Regression guard for the defect the inverse-CDF sampler fixes: the CTMC
        // walk of Map_sample presumes a probabilistic reading of (D0, D1), which
        // this pair does not admit (D1 = (-A*e)*alpha has negative entries because
        // alpha does). Running it on the same process must NOT reproduce the
        // analytic moments; if it did, the ME-specific sampler would be pointless.
        double m1 = 0.532854185661149;
        double m2 = 1.0460915848006274;
        double[] bad = jline.api.mam.Map_sample.map_sample(
                me.D(0), me.D(1), SAMPLE_SIZE, new Random(20260719));

        double badM1 = sampleMean(bad);
        double badM2 = sampleMoment(bad, 2);
        boolean meanWrong = Math.abs(badM1 - m1) > 0.05 * m1;
        boolean secondWrong = Math.abs(badM2 - m2) > 0.05 * m2;
        assertTrue(meanWrong || secondWrong,
                "CTMC walk unexpectedly reproduced the non-PH ME moments: mean "
                        + badM1 + " vs " + m1 + ", second moment " + badM2 + " vs " + m2);
    }

    @Test
    public void testMESamplingEmpiricalCDF() {
        ME me = nonPhaseTypeME();
        double[] x = me.sample(SAMPLE_SIZE, new Random(4241));
        java.util.Arrays.sort(x);

        // Kolmogorov-Smirnov statistic against F(t) = 1 - alpha*expm(A*t)*e. The
        // asymptotic 1% critical value is 1.63/sqrt(n) = 3.6e-3 at n = 200000;
        // the bound below is deliberately looser so that the test measures
        // sampler correctness rather than the realisation of one fixed seed.
        double ksMax = 0.0;
        for (int i = 0; i < x.length; i++) {
            double f = me.evalCDF(x[i]);
            double dPlus = (i + 1.0) / x.length - f;
            double dMinus = f - (double) i / x.length;
            ksMax = Math.max(ksMax, Math.max(dPlus, dMinus));
        }
        assertTrue(ksMax < 0.01, "KS deviation " + ksMax + " exceeds 0.01");

        // Spot-check individual quantiles of the empirical CDF.
        double[] probes = {0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
        for (int k = 0; k < probes.length; k++) {
            double t = probes[k];
            int count = 0;
            for (int i = 0; i < x.length; i++) {
                if (x[i] <= t) {
                    count++;
                }
            }
            assertEquals(me.evalCDF(t), (double) count / x.length, 0.01,
                    "empirical CDF at t = " + t);
        }
    }

    @Test
    public void testMEErlangMatchesErlangObject() {
        // Regression guard: an ME that happens to be phase-type must sample like
        // the equivalent Erlang object, which uses the inherited CTMC walk.
        ME me = ME.fromErlang(3, 2.0);
        Erlang erl = new Erlang(2.0, 3);

        assertEquals(erl.getMean(), me.getMean(), LOOSE_FINE_TOL);
        assertEquals(erl.getSCV(), me.getSCV(), LOOSE_FINE_TOL);

        double[] xMe = me.sample(SAMPLE_SIZE, new Random(1301));
        double[] xErl = erl.sample(SAMPLE_SIZE, new Random(1301));

        // Both estimates carry Monte Carlo error of order 1/sqrt(n) = 2.2e-3
        // relative, so compare them at a few percent of the analytic value.
        assertEquals(sampleMean(xErl), sampleMean(xMe), 0.02 * erl.getMean());
        assertEquals(sampleSCV(xErl), sampleSCV(xMe), 0.05 * erl.getSCV());
        assertEquals(erl.getMean(), sampleMean(xMe), 0.02 * erl.getMean());
        assertEquals(erl.getSCV(), sampleSCV(xMe), 0.05 * erl.getSCV());
    }

    // ------------------------------------------------------------------
    // Phase decomposition. mu_i = -A(i,i) and phi = D1*e ./ mu are exit rates
    // and exit probabilities only when (alpha, A) is a phase-type pair, which an
    // ME representation need not be.
    // ------------------------------------------------------------------

    @Test
    public void testMEComputesPhaseDecompositionLikeMatlab() {
        // ME inherits the Markovian formulas rather than refusing them, matching
        // MATLAB and native Python. An earlier revision threw here; the refusal
        // was dropped on 2026-07-20 because it also discarded the correct answer
        // in the degenerate case, which is what this test now pins.
        //
        // Degenerate case: an ME that IS a phase-type has independently known
        // values, so this is an identity and not a matter of interpretation.
        ME phMe = ME.fromErlang(3, 2.0);
        Erlang erl = new Erlang(2.0, 3);
        Matrix meMu = phMe.getMu();
        Matrix mePhi = phMe.getPhi();
        Matrix erlMu = erl.getMu();
        Matrix erlPhi = erl.getPhi();
        assertEquals(erlMu.getNumRows(), meMu.getNumRows());
        for (int i = 0; i < erlMu.getNumRows(); i++) {
            assertEquals(erlMu.get(i, 0), meMu.get(i, 0), LOOSE_FINE_TOL,
                    "mu must match the equivalent Erlang at phase " + i);
            assertEquals(erlPhi.get(i, 0), mePhi.get(i, 0), LOOSE_FINE_TOL,
                    "phi must match the equivalent Erlang at phase " + i);
        }

        // A genuine (non phase-type) ME still returns values; they simply carry
        // no probabilistic reading, and may leave [0,1]. The contract is that
        // the call succeeds and applies the same formulas, not that the result
        // is a probability.
        ME me = nonPhaseTypeME();
        Matrix mu = me.getMu();
        Matrix phi = me.getPhi();
        assertEquals(3, mu.getNumRows());
        assertEquals(3, phi.getNumRows());
        for (int i = 0; i < 3; i++) {
            assertEquals(-me.getProcess().get(0).get(i, i), mu.get(i, 0), LOOSE_FINE_TOL,
                    "mu_i must be -D0(i,i) at phase " + i);
        }
        assertNotNull(me.getInitProb());
    }

    @Test
    public void testMEStructRecordsProcessAndNaNPhases() {
        Network model = new Network("M/ME/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1", 0);
        source.setArrival(oclass, Exp.fitRate(0.5));
        queue.setService(oclass, nonPhaseTypeME());
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.addConnection(oclass, oclass, source, queue, 1.0);
        rm.addConnection(oclass, oclass, queue, sink, 1.0);
        model.link(rm);

        NetworkStruct sn = model.getStruct(false);

        // The (D0, D1) matrices reach the struct in full: they are the exact and
        // complete description of the process and are what LDES reads.
        MatrixCell proc = sn.proc.get(queue).get(oclass);
        assertEquals(3, proc.get(0).getNumRows());
        assertEquals(-1.0, proc.get(0).get(0, 0), LOOSE_FINE_TOL);
        assertEquals(-2.0, proc.get(0).get(2, 1), LOOSE_FINE_TOL);

        // The phase decomposition is recorded, not sentinelled: Station applies
        // the same Markovian formulas it applies to any other Markovian process,
        // matching MATLAB Station.m, which lists ME and RAP in the same arm as
        // Exp/Coxian/Erlang/PH/APH. An earlier revision wrote NaN here.
        Matrix mu = sn.mu.get(queue).get(oclass);
        Matrix phi = sn.phi.get(queue).get(oclass);
        assertEquals(3, mu.getNumRows());
        assertEquals(3, phi.getNumRows());
        for (int i = 0; i < 3; i++) {
            assertFalse(Double.isNaN(mu.get(i, 0)), "mu must be populated at phase " + i);
            assertFalse(Double.isNaN(phi.get(i, 0)), "phi must be populated at phase " + i);
            assertEquals(-proc.get(0).get(i, i), mu.get(i, 0), LOOSE_FINE_TOL,
                    "mu_i must be -D0(i,i) at phase " + i);
        }
    }

    /**
     * End-to-end regression: an ME that is really a phase-type distribution must
     * still solve under LDES and agree with the equivalent native object.
     */
    @Test
    public void testMEPhaseTypeSolvesUnderLDES() {
        double qLenErlang = solveMM1QueueLength(new Erlang(3.0, 3));
        double qLenMeErlang = solveMM1QueueLength(ME.fromErlang(3, 3.0));
        assertEquals(qLenErlang, qLenMeErlang, COARSE_TOL);

        double[] p = {0.3, 0.7};
        double[] rates = {1.0, 4.0};
        double qLenHyper = solveMM1QueueLength(new HyperExp(p, rates));
        double qLenMeHyper = solveMM1QueueLength(ME.fromHyperExp(p, rates));
        assertEquals(qLenHyper, qLenMeHyper, COARSE_TOL);
    }

    /** Solves an M/G/1 with the given service process and returns the queue length. */
    private static double solveMM1QueueLength(Distribution service) {
        Network model = new Network("M/G/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1", 0);
        source.setArrival(oclass, Exp.fitRate(0.5));
        queue.setService(oclass, service);
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.addConnection(oclass, oclass, source, queue, 1.0);
        rm.addConnection(oclass, oclass, queue, sink, 1.0);
        model.link(rm);

        SolverOptions options = new SolverOptions();
        options.seed = 23000;
        options.samples = 50000;
        options.verbose = VerboseLevel.SILENT;
        NetworkAvgTable table = new SolverLDES(model, options).getAvgTable();
        return table.getQLen().get(1);
    }
}
