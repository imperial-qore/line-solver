/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

/**
 * Tests for the Ott-Yashkov sojourn time transform of the M/G/1-PS queue.
 *
 * <p>The references are exact where one exists: the conditional mean is x/(1-rho)
 * for any service law, the M/M/1-PS unconditional moments are the Coffman-Muntz-
 * Trotter values that {@code Qsys_mm1_ps} returns, and the k = 0 atom is
 * (1-rho)*exp(-lambda*x). The conditional CDF values were additionally checked
 * against 40000 replications of an exact event-driven processor-sharing
 * simulation, and the unconditional CDF against SolverLDES.
 */
public class QsysMg1PsTest {
    private static final double LAMBDA = 0.7;
    private static final double MU = 1.0;
    private static final double[] EXP_ALPHA = {1.0};
    private static final double[][] EXP_T = {{-1.0}};

    @Test
    public void testConditionalMeanIsInsensitive() {
        double[] x = {0.25, 1.0, 4.0};
        QsysMg1PsResult r = Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, EXP_ALPHA, EXP_T, x, null, null, 41);
        for (int i = 0; i < x.length; i++) {
            assertEquals(x[i] / (1 - r.rho), r.meanCond[i], 1e-12);
        }
        assertEquals(1.0 / (MU * (1 - LAMBDA / MU)), r.meanUncond, 1e-12);
    }

    @Test
    public void testMm1PsMomentsMatchCoffmanMuntzTrotter() {
        QsysMg1PsResult r = Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, EXP_ALPHA, EXP_T, null, null, null, 41);
        double rho = LAMBDA / MU;
        double m2exact = 4.0 / (MU * MU * (1 - rho) * (1 - rho) * (2 - rho));
        assertEquals(m2exact, r.m2Uncond, 1e-6 * m2exact);
        QsysMm1PsResult ref = Qsys_mm1_ps.qsys_mm1_ps(new double[]{LAMBDA}, new double[]{MU});
        assertEquals(ref.W[0], r.meanUncond, 1e-12);
        assertEquals(ref.W2[0], r.m2Uncond, 1e-6 * ref.W2[0]);
    }

    @Test
    public void testTransformAtTheOriginAndAtLargeArguments() {
        QsysMg1PsResult r = Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, EXP_ALPHA, EXP_T, null, null, null, 41);
        assertEquals(1.0, r.lstCond(Complex.ZERO, 2.0).getReal(), 1e-12,
                "the transform is one at the origin");
        double atom = (1 - r.rho) * Math.exp(-LAMBDA * 2.0);
        assertEquals(atom, r.lstExcess(new Complex(1e8, 0), 2.0).getReal(), 1e-7,
                "the excess transform tends to the atom at t = x");
    }

    @Test
    public void testConditionalDistributionAgainstSimulation() {
        double[] t = {1.2, 1.5, 1.8, 2.4, 3.0, 4.0, 6.0, 10.0, 16.0};
        double[] sim = {0.207350, 0.282000, 0.342050, 0.481000, 0.584625,
            0.715825, 0.866050, 0.969875, 0.996150};
        QsysMg1PsResult r = Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, EXP_ALPHA, EXP_T,
                new double[]{1.0}, null, t, 41);
        for (int j = 0; j < t.length; j++) {
            assertEquals(sim[j], r.cdfCond[0][j], 5e-3,
                    "conditional CDF at t = " + t[j] + " within the simulation error");
            if (j > 0) {
                assertTrue(r.cdfCond[0][j] >= r.cdfCond[0][j - 1], "the CDF is monotone");
            }
        }
        assertEquals((1 - r.rho) * Math.exp(-LAMBDA), r.atomCond[0], 1e-12);
    }

    @Test
    public void testNoMassBelowTheServiceRequirement() {
        QsysMg1PsResult r = Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, EXP_ALPHA, EXP_T,
                new double[]{2.0}, null, new double[]{0.5, 1.9, 2.0}, 41);
        assertEquals(0.0, r.cdfCond[0][0], 0.0, "V(x) >= x leaves no mass below x");
        assertEquals(0.0, r.cdfCond[0][1], 0.0);
        assertEquals(r.atomCond[0], r.cdfCond[0][2], 1e-12, "at t = x the CDF is the atom");
    }

    @Test
    public void testDensityIsUndefinedOnTheLattice() {
        QsysMg1PsResult r = Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, EXP_ALPHA, EXP_T,
                new double[]{1.0}, null, new double[]{1.5, 2.0, 3.0}, 41);
        assertTrue(Double.isFinite(r.pdfCond[0][0]), "a density exists off the lattice");
        assertTrue(Double.isNaN(r.pdfCond[0][1]), "V(x) has an atom at 2x, so no density");
        assertTrue(Double.isNaN(r.pdfCond[0][2]), "V(x) has an atom at 3x, so no density");
    }

    @Test
    public void testPhaseTypeAgreesWithTheTransformHandle() {
        final double mu = 4.0;
        double lam = 1.2;
        double[] x = {0.25, 1.0, 3.0};
        double[] s = {0.2, 2.0};
        QsysMg1PsResult ph = Qsys_mg1_ps.qsys_mg1_ps(lam, new double[]{1.0, 0.0},
                new double[][]{{-mu, mu}, {0.0, -mu}}, x, s, null, 41);
        Qsys_mg1_ps.ServiceLST bhat = new Qsys_mg1_ps.ServiceLST() {
            public Complex value(Complex tau) {
                return new Complex(mu, 0).divide(tau.add(mu)).pow(2.0);
            }
        };
        Qsys_mg1_ps.ServicePDF bpdf = new Qsys_mg1_ps.ServicePDF() {
            public double value(double y) {
                return mu * mu * y * Math.exp(-mu * y);
            }
        };
        QsysMg1PsResult hd = Qsys_mg1_ps.qsys_mg1_ps(lam, bhat, 2.0 / mu, x, s, 41, bpdf);
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < s.length; j++) {
                assertEquals(ph.lstCondVal[i][j], hd.lstCondVal[i][j], 1e-8,
                        "the numerical inversion in tau reproduces the residue sum");
            }
        }
        for (int j = 0; j < s.length; j++) {
            assertEquals(ph.lstUncondVal[j], hd.lstUncondVal[j], 1e-8);
        }
        assertEquals(ph.m2Uncond, hd.m2Uncond, 1e-4 * ph.m2Uncond);
    }

    @Test
    public void testRepeatedPoleFallsBackToTheCompanionMatrix() {
        // at s = -(sqrt(mu)-sqrt(lambda))^2 the two poles of P coincide, so the residue
        // formula is unusable and the companion-matrix branch must take over. There the
        // confluent partial fraction is exact: Ahat/(tau-r)^2 inverts to
        // (a1 + (a1*r + a0)*x)*exp(r*x)
        double rho = LAMBDA / MU;
        double sstar = -Math.pow(Math.sqrt(MU) - Math.sqrt(LAMBDA), 2);
        double[] x = {0.25, 1.0, 3.0, 10.0};
        QsysMg1PsResult r = Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, EXP_ALPHA, EXP_T, x,
                new double[]{sstar}, null, 41);
        double a1 = 1 - rho;
        double a0 = (1 - rho) * (MU - LAMBDA) + sstar * rho;
        double rr = -(MU - sstar - LAMBDA) / 2;
        for (int i = 0; i < x.length; i++) {
            double exact = (1 - rho) / ((a1 + (a1 * rr + a0) * x[i]) * Math.exp(rr * x[i]));
            assertEquals(exact, r.lstCondVal[i][0], 1e-12 * Math.abs(exact),
                    "the confluent branch is exact at the double pole, x = " + x[i]);
        }
        double lo = r.lstCond(new Complex(sstar + 1e-6, 0), 1.0).getReal();
        double hi = r.lstCond(new Complex(sstar - 1e-6, 0), 1.0).getReal();
        assertTrue(lo < r.lstCondVal[1][0] && r.lstCondVal[1][0] < hi,
                "the residue branch either side of the confluence brackets it");
    }

    @Test
    public void testUnstableAndMalformedInputsAreRejected() {
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qsys_mg1_ps.qsys_mg1_ps(1.5, EXP_ALPHA, EXP_T, null, null, null, 41);
            }
        }, "rho >= 1 is rejected");
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, new double[]{0.5}, EXP_T, null, null, null, 41);
            }
        }, "alpha must be a probability vector");
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qsys_mg1_ps.qsys_mg1_ps(LAMBDA, EXP_ALPHA, EXP_T, null, null, null, 40);
            }
        }, "nterms must be odd");
    }
}
