/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.pfqn.PfqnResptPsResult;
import jline.api.qsys.QsysMm1PsResult;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static jline.api.pfqn.Pfqn_respt_ps_moments.pfqn_respt_ps_moments;
import static jline.api.qsys.Qsys_mm1_ps.qsys_mm1_ps;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for the processor-sharing sojourn-time moments of Mitra and Morrison
 * (1983): {@code Qsys_mm1_ps} for the open station and
 * {@code Pfqn_respt_ps_moments} for the closed terminal-driven system.
 *
 * <p>The mathematics is anchored to results that do not come from that paper:
 * the single-class M/M/1-PS second moment of Coffman, Muntz and Trotter (1970),
 * the trivial one-job closed system whose sojourn time is the service time
 * itself, and the fact that a lone job cannot be delayed. The paper's own two
 * routes are then cross-checked against each other, the exact route being the
 * reference and the asymptotic one required to converge to it as the expansion
 * parameter grows. The literal values are the MATLAB and native-Python outputs,
 * so a divergence in any codebase shows up here.</p>
 */
public class PfqnResptPsTest {

    private static void assertRel(double expected, double actual, double relTol, String msg) {
        double denom = Math.max(Math.abs(expected), 1e-300);
        assertTrue(Math.abs(actual - expected) <= relTol * denom,
                msg + ": expected " + expected + " but got " + actual);
    }

    @Test
    @DisplayName("single-class open PS reproduces Coffman-Muntz-Trotter")
    public void testOpenSingleClassCmt() {
        double mu = 2.0;
        double[] rhos = {0.2, 0.5, 0.8};
        for (int i = 0; i < rhos.length; i++) {
            double rho = rhos[i];
            QsysMm1PsResult r = qsys_mm1_ps(new double[]{rho * mu}, new double[]{mu});
            assertRel(1.0 / (mu * (1 - rho)), r.W[0], 1e-12, "M/M/1-PS mean");
            double cmt = 4.0 / (mu * mu * (1 - rho) * (1 - rho) * (2 - rho));
            assertRel(cmt, r.W2[0], 1e-12, "M/M/1-PS second moment");
        }
    }

    @Test
    @DisplayName("multiclass open PS matches the MATLAB and Python values")
    public void testOpenMulticlassParity() {
        QsysMm1PsResult r = qsys_mm1_ps(new double[]{0.3, 0.4}, new double[]{1.0, 3.0});
        assertRel(1.764705882353, r.W[0], 1e-10, "W1");
        assertRel(0.588235294118, r.W[1], 1e-10, "W2");
        assertRel(7.750865051903, r.W2[0], 1e-10, "second moment class 1");
        assertRel(0.927201263144, r.W2[1], 1e-10, "second moment class 2");
        assertRel(1.0 - 0.3 - 0.4 / 3.0, r.alpha, 1e-12, "alpha");
    }

    @Test
    @DisplayName("an unstable open PS station is an error, not a negative moment")
    public void testOpenUnstableThrows() {
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                qsys_mm1_ps(new double[]{1.0, 1.0}, new double[]{1.0, 2.0});
            }
        });
    }

    @Test
    @DisplayName("a single closed job is never delayed: W = S, E[W^2] = 2 S^2")
    public void testClosedSingleJob() {
        double S = 0.5;
        PfqnResptPsResult r = pfqn_respt_ps_moments(new double[]{S}, new double[]{1},
                new double[]{5.0}, "exact");
        assertRel(S, r.W[0], 1e-12, "lone job mean");
        assertRel(2 * S * S, r.W2[0], 1e-12, "lone job second moment");
    }

    @Test
    @DisplayName("closed exact route matches the MATLAB and Python values")
    public void testClosedExactParity() {
        PfqnResptPsResult r = pfqn_respt_ps_moments(new double[]{1.0, 0.5},
                new double[]{4, 2}, new double[]{50.0, 100.0}, "exact");
        assertRel(1.072358918582, r.W[0], 1e-9, "W1");
        assertRel(0.544470487795, r.W[1], 1e-9, "W2");
        assertRel(2.372013893895, r.W2[0], 1e-9, "second moment class 1");
        assertRel(0.623810622242, r.W2[1], 1e-9, "second moment class 2");
        assertEquals("exact", r.method[0]);
        assertEquals(12.0, r.nstates[0], 0.0);
        assertEquals(10.0, r.nstates[1], 0.0);
    }

    @Test
    @DisplayName("closed asymptotic route matches the MATLAB and Python values")
    public void testClosedAsymptoticParity() {
        PfqnResptPsResult r = pfqn_respt_ps_moments(new double[]{1.0, 0.5},
                new double[]{4, 2}, new double[]{50.0, 100.0}, "asymptotic");
        assertRel(1.072160744545, r.W[0], 1e-9, "W1");
        assertRel(2.370214800752, r.W2[0], 1e-9, "second moment class 1");
        assertRel(0.623075482738, r.W2[1], 1e-9, "second moment class 2");
        assertRel(2.392420108971, r.c0[0], 1e-9, "c0");
        assertRel(-4.441061643908, r.c1[0], 1e-9, "c1");
        assertEquals(200.0, r.expansionParam, 1e-12);
    }

    @Test
    @DisplayName("the asymptotic route converges to the exact one as Nexp grows")
    public void testAsymptoticConverges() {
        double[] q = {1.0, 2.0};
        int[] K = {3, 2};
        double prevErr = Double.MAX_VALUE;
        for (int scale = 1; scale <= 8; scale *= 2) {
            double[] think = {50.0 * scale, 100.0 * scale};
            double[] S = {1 / q[0], 1 / q[1]};
            double[] N = {K[0] + 1, K[1]};
            PfqnResptPsResult ex = pfqn_respt_ps_moments(S, N, think, "exact");
            PfqnResptPsResult as = pfqn_respt_ps_moments(S, N, think, "asymptotic");
            double err = Math.abs(as.W2[0] / ex.W2[0] - 1);
            assertTrue(err < prevErr, "the expansion must improve with Nexp, but "
                    + err + " >= " + prevErr);
            prevErr = err;
        }
        assertTrue(prevErr < 1e-5, "two terms of the expansion should be tight at "
                + "Nexp = 1600, got " + prevErr);
    }

    @Test
    @DisplayName("the leading closed term is the open formula")
    public void testClosedLeadingTermIsOpen() {
        // lambda_j = K_j/Z_j held fixed while K_j grows: the closed system tends to
        // the open one, and c0 is already exactly the open second moment
        double[] q = {1.0, 3.0};
        double[] lam = {0.3, 0.4};
        QsysMm1PsResult open = qsys_mm1_ps(lam, q);
        int[] pops = {10, 100, 1000};
        for (int i = 0; i < pops.length; i++) {
            int K = pops[i];
            double[] S = {1 / q[0], 1 / q[1]};
            double[] N = {K + 1, K};
            double[] Z = {K / lam[0], K / lam[1]};
            PfqnResptPsResult r = pfqn_respt_ps_moments(S, N, Z, "asymptotic");
            assertRel(open.W2[0], r.c0[0], 1e-10, "leading term at K = " + K);
        }
    }

    @Test
    @DisplayName("auto takes the exact route when the state space is small")
    public void testAutoRoutes() {
        PfqnResptPsResult small = pfqn_respt_ps_moments(new double[]{1.0},
                new double[]{5}, new double[]{50.0});
        assertEquals("exact", small.method[0]);
        PfqnResptPsResult big = pfqn_respt_ps_moments(new double[]{1.0},
                new double[]{100000}, new double[]{5.0e6});
        assertEquals("asymptotic", big.method[0]);
    }

    @Test
    @DisplayName("an empty class reports no sojourn time rather than a zero")
    public void testEmptyClassIsNaN() {
        PfqnResptPsResult r = pfqn_respt_ps_moments(new double[]{1.0, 0.5},
                new double[]{3, 0}, new double[]{50.0, 100.0});
        assertTrue(Double.isNaN(r.W[1]), "an unpopulated class has no sojourn time");
        assertTrue(!Double.isNaN(r.W[0]));
        assertEquals("none", r.method[1]);
    }
}
