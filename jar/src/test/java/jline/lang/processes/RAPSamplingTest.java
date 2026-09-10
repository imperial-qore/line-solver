/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.api.mam.Map_acf;
import jline.api.mam.Rap_sample;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static jline.TestTools.*;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Sampling tests for the Rational Arrival Process (RAP).
 *
 * <p>A RAP is sampled by conditional inversion: the conditional phase vector v is
 * propagated across events as {@code v <- v*expm(H0*x)*H1 / (v*expm(H0*x)*H1*e)}.
 * The vector update is what makes successive inter-event times correlated, so the
 * decisive test is that a sampled sequence reproduces the analytic lag-1
 * autocorrelation of the process, not merely its marginal moments.</p>
 */
public class RAPSamplingTest {

    /** Number of draws used by the Monte Carlo tests. */
    private static final int SAMPLE_SIZE = 200000;

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
        RandomManager.setMasterSeed(23000);
    }

    /**
     * A two-state MAP with a substantial positive lag-1 autocorrelation.
     *
     * <p>Analytic values: mean 0.55, SCV 2.3388, acf(1) 0.22898.</p>
     */
    private static Matrix correlatedD0() {
        return new Matrix(new double[][]{
                {-10.0, 0.0},
                {0.0, -1.0}
        });
    }

    private static Matrix correlatedD1() {
        return new Matrix(new double[][]{
                {9.0, 1.0},
                {0.1, 0.9}
        });
    }

    private static double sampleMean(double[] x) {
        double s = 0.0;
        for (int i = 0; i < x.length; i++) {
            s += x[i];
        }
        return s / x.length;
    }

    private static double sampleSCV(double[] x) {
        double m1 = sampleMean(x);
        double m2 = 0.0;
        for (int i = 0; i < x.length; i++) {
            m2 += x[i] * x[i];
        }
        m2 /= x.length;
        return (m2 - m1 * m1) / (m1 * m1);
    }

    /** Empirical lag-k autocorrelation of a sequence. */
    private static double sampleACF(double[] x, int lag) {
        double m = sampleMean(x);
        double num = 0.0;
        double den = 0.0;
        for (int i = 0; i < x.length; i++) {
            den += (x[i] - m) * (x[i] - m);
            if (i + lag < x.length) {
                num += (x[i] - m) * (x[i + lag] - m);
            }
        }
        return (num / (x.length - lag)) / (den / x.length);
    }

    @Test
    public void testRAPConstructionFromMAP() {
        MAP map = new MAP(correlatedD0(), correlatedD1());
        RAP rap = RAP.fromMAP(map);

        assertNotNull(rap);
        assertEquals(2, rap.getNumberOfPhases());
        assertEquals(0.55, rap.getMean(), LOOSE_FINE_TOL);
        assertEquals(2.338842975206611, rap.getSCV(), LOOSE_FINE_TOL);
    }

    /**
     * The decisive test: the sampled sequence must carry the autocorrelation of
     * the process. Getting this wrong (for instance by re-initialising the phase
     * vector from the embedded stationary vector at every draw) leaves the
     * marginal moments correct while destroying the correlation structure, which
     * is exactly the failure mode that silently corrupts M/RAP/1 simulations.
     */
    @Test
    public void testRAPSampleReproducesLag1Autocorrelation() {
        Matrix D0 = correlatedD0();
        Matrix D1 = correlatedD1();
        RAP rap = new RAP(D0, D1);

        double acf1 = Map_acf.map_acf(D0, D1, 1).get(0);
        assertEquals(0.2289752650176679, acf1, LOOSE_FINE_TOL);

        double[] x = rap.sample(SAMPLE_SIZE, new Random(90210));

        // Marginal: the conditional-vector update must not disturb the ME marginal.
        assertEquals(rap.getMean(), sampleMean(x), 0.02 * rap.getMean());
        assertEquals(rap.getSCV(), sampleSCV(x), 0.05 * rap.getSCV());

        // Correlation: the standard error of a lag-1 autocorrelation estimate at
        // n = 200000 is of order 1/sqrt(n) = 2.2e-3, so 0.02 absolute is a wide
        // but still discriminating band (an i.i.d. sequence would give ~0).
        assertEquals(acf1, sampleACF(x, 1), 0.02);
    }

    @Test
    public void testRAPSampleReproducesHigherLagAutocorrelations() {
        Matrix D0 = correlatedD0();
        Matrix D1 = correlatedD1();
        RAP rap = new RAP(D0, D1);
        double[] x = rap.sample(SAMPLE_SIZE, new Random(90210));

        int[] lags = {2, 3, 5};
        for (int k = 0; k < lags.length; k++) {
            double expected = Map_acf.map_acf(D0, D1, lags[k]).get(0);
            assertEquals(expected, sampleACF(x, lags[k]), 0.02,
                    "lag-" + lags[k] + " autocorrelation");
        }
    }

    /**
     * A renewal RAP (Poisson) has zero autocorrelation at every lag. This is the
     * negative control for the previous test: it shows the estimator is not
     * manufacturing correlation of its own.
     */
    @Test
    public void testRAPPoissonHasNoAutocorrelation() {
        RAP rap = RAP.fromPoisson(2.0);
        double[] x = rap.sample(SAMPLE_SIZE, new Random(31337));

        assertEquals(0.5, sampleMean(x), 0.02 * 0.5);
        assertEquals(1.0, sampleSCV(x), 0.05);
        assertEquals(0.0, sampleACF(x, 1), 0.02);
    }

    @Test
    public void testRAPErlangRenewalMoments() {
        RAP rap = RAP.fromErlang(3, 2.0);
        double[] x = rap.sample(SAMPLE_SIZE, new Random(31337));

        // Erlang(3, 2): mean 1.5, SCV 1/3, and no autocorrelation (renewal).
        assertEquals(1.5, sampleMean(x), 0.02 * 1.5);
        assertEquals(1.0 / 3.0, sampleSCV(x), 0.05 / 3.0);
        assertEquals(0.0, sampleACF(x, 1), 0.02);
    }

    /**
     * A genuine RAP whose H0 has a negative off-diagonal entry, so it admits no
     * MAP representation and the CTMC walk of Map_sample cannot be applied. The
     * conditional-inversion sampler must still reproduce the marginal moments.
     */
    @Test
    public void testRAPWithNegativeOffDiagonalSampling() {
        Matrix H0 = new Matrix(new double[][]{
                {-1.0, 0.0, 0.0},
                {0.0, -2.0, 2.0},
                {0.0, -2.0, -2.0}
        });
        // H1 = (-H0*e)*alpha with the non-phase-type alpha of MEDistributionTest,
        // so the RAP is the ME renewal process and its moments are known exactly.
        // new Matrix(double[]) builds a column vector, so transpose to get the
        // 1 x n row vector required by the outer product below.
        Matrix alpha = new Matrix(new double[]{
                0.61058991931158258, -0.15547146730086722, 0.54488154798928464}).transpose();
        Matrix exit = H0.mult(Matrix.ones(3, 1)).scale(-1.0);
        Matrix H1 = exit.mult(alpha);

        RAP rap = new RAP(H0, H1);
        assertTrue(H0.get(2, 1) < 0.0, "H0 must have a negative off-diagonal entry");

        double[] x = Rap_sample.rap_sample(H0, H1, SAMPLE_SIZE, new Random(5150));
        double m1 = 0.532854185661149;
        double m2 = 1.0460915848006274;
        double empM2 = 0.0;
        for (int i = 0; i < x.length; i++) {
            assertTrue(x[i] >= 0.0, "Sample should be non-negative");
            empM2 += x[i] * x[i];
        }
        empM2 /= x.length;

        assertEquals(m1, sampleMean(x), 0.02 * m1);
        assertEquals(m2, empM2, 0.05 * m2);
        // A renewal process: the conditional vector resets to alpha after every
        // event, so there must be no autocorrelation.
        assertEquals(0.0, sampleACF(x, 1), 0.02);
    }

    /**
     * mu_i = -H0(i,i) and phi = H1*e ./ mu are rates and probabilities only for a
     * MAP. A RAP must refuse the decomposition rather than return entries outside
     * [0,1], and the refusal is unconditional so that no caller depends on a
     * quantity that is only sometimes meaningful.
     */
    @Test
    public void testRAPComputesPhaseDecompositionLikeMAP() {
        // A RAP whose matrices are nonnegative IS a MAP, so the two objects must
        // return identical phase quantities on identical input. That identity is
        // the reason RAP no longer refuses the decomposition: an earlier
        // revision threw here, which cannot reproduce a value that is known
        // independently. Native Python had the mirror-image defect, where RAP
        // returned H1 row sums and a vector of ones and so contradicted its own
        // MAP; both now use the shared Markovian formulas.
        Matrix d0 = correlatedD0();
        Matrix d1 = correlatedD1();
        RAP rap = new RAP(d0, d1);
        MAP map = new MAP(d0, d1);

        Matrix rapMu = rap.getMu();
        Matrix mapMu = map.getMu();
        Matrix rapPhi = rap.getPhi();
        Matrix mapPhi = map.getPhi();
        Matrix rapAlpha = rap.getInitProb();
        Matrix mapAlpha = map.getInitProb();

        assertEquals(mapMu.getNumRows(), rapMu.getNumRows());
        for (int i = 0; i < mapMu.getNumRows(); i++) {
            assertEquals(mapMu.get(i, 0), rapMu.get(i, 0), 1e-12,
                    "mu must match the equivalent MAP at phase " + i);
            assertEquals(mapPhi.get(i, 0), rapPhi.get(i, 0), 1e-12,
                    "phi must match the equivalent MAP at phase " + i);
            // mu_i = -D0(i,i) is the analytic value, pinned so the identity
            // cannot pass by both objects sharing one wrong implementation.
            assertEquals(-d0.get(i, i), rapMu.get(i, 0), 1e-12,
                    "mu_i must be -D0(i,i) at phase " + i);
        }
        assertEquals(mapAlpha.getNumElements(), rapAlpha.getNumElements());
        for (int i = 0; i < mapAlpha.getNumElements(); i++) {
            assertEquals(mapAlpha.get(0, i), rapAlpha.get(0, i), 1e-12,
                    "alpha must match the equivalent MAP at phase " + i);
        }
    }
}
