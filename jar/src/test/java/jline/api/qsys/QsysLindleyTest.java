/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

/**
 * Tests for the conditional Lindley recursion algorithms.
 *
 * <p>The expected values are the ones the MATLAB twins produce, and were
 * themselves validated against 4e6 Monte Carlo replications of the corresponding
 * one-step experiment.
 */
public class QsysLindleyTest {
    private static final double TOL = 1e-9;

    @Test
    public void testMm1ConditionalMomentsMatchTheorem1() {
        QsysLindleyResult r = Qsys_mm1_lindley.qsys_mm1_lindley(0.8, 1.0,
                new double[] {2.0}, 3);
        assertEquals(1.890206, r.moments[0][0], 1e-6);
        assertEquals(5.274485, r.moments[0][1], 1e-6);
        assertEquals(18.220680, r.moments[0][2], 1e-6);
        assertEquals(1.890206, r.mean[0], 1e-6, "the mean agrees with corollary 2");
        assertEquals(r.moments[0][1] - r.moments[0][0] * r.moments[0][0], r.var[0], TOL,
                "the variance is the second central moment");
    }

    @Test
    public void testMm1ConditionalMomentsFromAnEmptyQueue() {
        QsysLindleyResult r = Qsys_mm1_lindley.qsys_mm1_lindley(0.8, 1.0,
                new double[] {0.0}, 3);
        assertEquals(0.444444, r.moments[0][0], 1e-6);
        assertEquals(0.888889, r.moments[0][1], 1e-6);
        assertEquals(2.666667, r.moments[0][2], 1e-6);
    }

    @Test
    public void testMm1ConditionalMeanAgreesWithTheMomentForm() {
        // corollary 2 and theorem 1 at m = 1 are algebraically the same expression
        double[] wn = {0.0, 0.1, 1.0, 3.7, 25.0};
        for (double lambda : new double[] {0.3, 0.8, 1.0, 1.9}) {
            for (double mu : new double[] {0.5, 1.0, 2.0}) {
                QsysLindleyResult r = Qsys_mm1_lindley.qsys_mm1_lindley(lambda, mu, wn);
                for (int i = 0; i < wn.length; i++) {
                    assertEquals(r.moments[i][0], r.mean[i], 1e-9,
                            "lambda=" + lambda + " mu=" + mu + " Wn=" + wn[i]);
                    assertTrue(r.var[i] >= -1e-12, "the variance is nonnegative");
                }
            }
        }
    }

    @Test
    public void testMm1IsFiniteAboveSaturation() {
        // the conditional moments are a one-step expectation, so they exist even
        // when the queue is unstable
        QsysLindleyResult r = Qsys_mm1_lindley.qsys_mm1_lindley(2.0, 1.0,
                new double[] {1.0}, 3);
        for (int m = 0; m < 3; m++) {
            assertTrue(Double.isFinite(r.moments[0][m]), "order " + (m + 1) + " is finite");
        }
        assertTrue(r.var[0] > 0.0);
    }

    @Test
    public void testMm1ConditionalMomentsAgainstDirectMonteCarlo() {
        Random rng = new Random(12345);
        int n = 400000;
        double lambda = 0.8;
        double mu = 1.0;
        double wn = 2.0;
        double s1 = 0.0;
        double s2 = 0.0;
        for (int i = 0; i < n; i++) {
            double a = -Math.log(1.0 - rng.nextDouble()) / lambda;
            double s = -Math.log(1.0 - rng.nextDouble()) / mu;
            double w = Math.max(wn + s - a, 0.0);
            s1 += w;
            s2 += w * w;
        }
        QsysLindleyResult r = Qsys_mm1_lindley.qsys_mm1_lindley(lambda, mu,
                new double[] {wn});
        assertEquals(r.moments[0][0], s1 / n, 0.02);
        assertEquals(r.moments[0][1], s2 / n, 0.10);
    }

    @Test
    public void testHh1ConditionalMomentsMatchTheorem3() {
        QsysLindleyResult r = Qsys_hh1_lindley.qsys_hh1_lindley(
                new double[] {0.5, 2.0}, new double[] {0.4, 0.6},
                new double[] {1.0, 4.0}, new double[] {0.7, 0.3},
                new double[] {1.0}, 2);
        assertEquals(1.048425, r.mean[0], 1e-6);
        assertEquals(2.141581, r.moments[0][1], 1e-6);
    }

    @Test
    public void testHh1DegenerateMixtureReducesToMm1() {
        double[] wn = {0.0, 1.0, 2.0, 7.5};
        QsysLindleyResult hh = Qsys_hh1_lindley.qsys_hh1_lindley(
                new double[] {0.8}, new double[] {1.0},
                new double[] {1.0}, new double[] {1.0}, wn, 3);
        QsysLindleyResult mm = Qsys_mm1_lindley.qsys_mm1_lindley(0.8, 1.0, wn, 3);
        for (int i = 0; i < wn.length; i++) {
            for (int m = 0; m < 3; m++) {
                assertEquals(mm.moments[i][m], hh.moments[i][m], 1e-12,
                        "Wn=" + wn[i] + " order " + (m + 1));
            }
        }
    }

    @Test
    public void testHh1VarianceExceedsTheMixtureOfPhaseVariances() {
        // the phase is itself random, so the between-phase spread of the means adds
        double[] lambda = {0.5, 2.0};
        double[] pa = {0.4, 0.6};
        double[] mu = {1.0, 4.0};
        double[] ps = {0.7, 0.3};
        double[] wn = {1.0};
        QsysLindleyResult mixed = Qsys_hh1_lindley.qsys_hh1_lindley(lambda, pa, mu, ps, wn);
        double avgPhaseVar = 0.0;
        for (int i = 0; i < lambda.length; i++) {
            for (int j = 0; j < mu.length; j++) {
                QsysLindleyResult one = Qsys_mm1_lindley.qsys_mm1_lindley(
                        lambda[i], mu[j], wn);
                avgPhaseVar += pa[i] * ps[j] * one.var[0];
            }
        }
        assertTrue(mixed.var[0] > avgPhaseVar,
                "the mixture variance exceeds the average phase variance");
    }

    @Test
    public void testTandemConditionalMean() {
        assertEquals(0.345679, Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(0.8, 1.0, 1.0, 0.0, 0.0).mean[0], 1e-6);
        assertEquals(2.906058, Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(0.8, 1.0, 1.0, 2.0, 3.0).mean[0], 1e-6);
        assertEquals(0.527153, Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(0.8, 1.0, 1.5, 5.0, 0.5).mean[0], 1e-6);
        assertEquals(3.486517, Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(0.5, 2.0, 1.2, 1.0, 4.0).mean[0], 1e-6);
    }

    @Test
    public void testTandemConditionalMeanEqualRateBranch() {
        // lambda == mu1 makes the interdeparture time Erlang(2,mu1); the branch must
        // agree with the limit of the unequal-rate expression
        double target = Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(1.0, 1.0, 1.0, 1.0, 1.0).mean[0];
        assertEquals(1.084585, target, 1e-6);
        double near = Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(1.0 + 1e-6, 1.0, 1.0, 1.0, 1.0).mean[0];
        assertEquals(target, near, 1e-6, "the two branches agree in the limit");
    }

    @Test
    public void testTandemConditionalMeanTendsToTheBusyServerLimit() {
        // as the upstream wait grows the server never idles and the interdeparture
        // time is its own service time
        double lambda = 0.8;
        double mu1 = 1.0;
        double mu2 = 1.3;
        double y = 2.0;
        QsysTandemLindleyResult far = Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(lambda, mu1, mu2, 500.0, y);
        // E[g(S1)] with S1 ~ Exp(mu1), evaluated by direct quadrature
        double acc = 0.0;
        int steps = 4000000;
        double upper = 60.0;
        double h = upper / steps;
        for (int i = 0; i < steps; i++) {
            double d = (i + 0.5) * h;
            double g = d <= y ? y - d + 1.0 / mu2 : Math.exp(-mu2 * (d - y)) / mu2;
            acc += mu1 * Math.exp(-mu1 * d) * g * h;
        }
        assertEquals(acc, far.mean[0], 1e-6, "the busy-server limit is E[g(S1)]");
    }

    @Test
    public void testTandemInterdepartureMean() {
        // an always-busy upstream server passes on its own service time
        QsysTandemLindleyResult far = Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(0.9, 1.0, 1.0, 200.0, 0.0);
        assertEquals(1.0, far.interdepMean[0], 1e-9);
        assertEquals(0.0, far.idleProb[0], 1e-9);
        // from an empty upstream queue it is E[max(A,S1)]
        QsysTandemLindleyResult zero = Qsys_mm1_tandem_lindley
                .qsys_mm1_tandem_lindley(0.8, 1.0, 1.0, 0.0, 0.0);
        assertEquals(1.0 / 0.8 + 1.0 - 1.0 / 1.8, zero.interdepMean[0], 1e-9);
    }

    @Test
    public void testTandemPathRecursionMatchesDirectEventSimulation() {
        Random rng = new Random(11);
        int n = 4000;
        int k = 4;
        double[] a = new double[n];
        double[][] s = new double[n][k];
        for (int i = 0; i < n; i++) {
            a[i] = -Math.log(1.0 - rng.nextDouble()) / 0.8;
            for (int j = 0; j < k; j++) {
                s[i][j] = -Math.log(1.0 - rng.nextDouble());
            }
        }
        QsysTandemPathResult r = Qsys_tandem_lindley.qsys_tandem_lindley(a, s);

        // direct event-driven tandem: arrival at station j is departure from j-1
        double[] epoch = new double[n];
        for (int i = 1; i < n; i++) {
            epoch[i] = epoch[i - 1] + a[i - 1];
        }
        double[][] wd = new double[n][k];
        double[][] dep = new double[n][k];
        for (int i = 0; i < n; i++) {
            double arrival = epoch[i];
            for (int j = 0; j < k; j++) {
                double prev = i > 0 ? dep[i - 1][j] : Double.NEGATIVE_INFINITY;
                wd[i][j] = Math.max(prev - arrival, 0.0);
                dep[i][j] = Math.max(arrival, prev) + s[i][j];
                arrival = dep[i][j];
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                assertEquals(wd[i][j], r.W[i][j], 1e-9,
                        "W at customer " + i + " station " + j);
                assertEquals(dep[i][j], r.departure[i][j], 1e-9,
                        "departure at customer " + i + " station " + j);
            }
        }
    }

    @Test
    public void testTandemPathRecursionFirstStationIsPlainLindley() {
        Random rng = new Random(3);
        int n = 500;
        double[] a = new double[n];
        double[][] s = new double[n][2];
        for (int i = 0; i < n; i++) {
            a[i] = -Math.log(1.0 - rng.nextDouble()) / 0.7;
            s[i][0] = -Math.log(1.0 - rng.nextDouble());
            s[i][1] = -Math.log(1.0 - rng.nextDouble());
        }
        QsysTandemPathResult r = Qsys_tandem_lindley.qsys_tandem_lindley(a, s);
        double w = 0.0;
        for (int i = 0; i < n - 1; i++) {
            w = Math.max(w + s[i][0] - a[i], 0.0);
            assertEquals(w, r.W[i + 1][0], 1e-12, "station 1 follows Lindley");
        }
    }

    @Test
    public void testTandemPathInterdepartureIdentity() {
        // the two equivalent forms of the interdeparture time must agree on the path
        Random rng = new Random(19);
        int n = 800;
        int k = 3;
        double[] a = new double[n];
        double[][] s = new double[n][k];
        for (int i = 0; i < n; i++) {
            a[i] = -Math.log(1.0 - rng.nextDouble()) / 0.75;
            for (int j = 0; j < k; j++) {
                s[i][j] = -Math.log(1.0 - rng.nextDouble());
            }
        }
        QsysTandemPathResult r = Qsys_tandem_lindley.qsys_tandem_lindley(a, s);
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < k - 1; j++) {
                double idle = Math.max(r.G[i][j] - r.W[i][j] - s[i][j], 0.0) + s[i + 1][j];
                double diff = r.G[i][j] + r.W[i + 1][j] - r.W[i][j] + s[i + 1][j] - s[i][j];
                assertEquals(r.G[i][j + 1], idle, 1e-12, "idle-time form");
                assertEquals(r.G[i][j + 1], diff, 1e-12, "difference form");
            }
        }
    }

    @Test
    public void testInputValidation() {
        assertThrows(IllegalArgumentException.class,
                () -> Qsys_mm1_lindley.qsys_mm1_lindley(0.8, 1.0, new double[] {-1.0}));
        assertThrows(IllegalArgumentException.class,
                () -> Qsys_hh1_lindley.qsys_hh1_lindley(new double[] {1.0, 2.0},
                        new double[] {0.5, 0.9}, new double[] {1.0}, new double[] {1.0},
                        new double[] {0.0}));
        assertThrows(IllegalArgumentException.class,
                () -> Qsys_mm1_tandem_lindley.qsys_mm1_tandem_lindley(0.8, 1.0, 1.0,
                        new double[] {1.0, 2.0}, new double[] {1.0}));
        assertThrows(IllegalArgumentException.class,
                () -> Qsys_tandem_lindley.qsys_tandem_lindley(new double[] {1.0, 2.0},
                        new double[][] {{1.0}}));
    }
}
