/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Analytic validation of the BMAP/PH/N/N retrial engine on the M/G/1 retrial
 * orbit, whose mean is known in closed form:
 *
 *   E[N_orbit] = rho^2/(1-rho) * ((1+cB^2)/2 + 1/nu)
 *
 * with constant per-customer retrial rate nu. The heavy-tailed cases (cB = 5,
 * 10) are the regression: a truncation level pinned to rho alone cuts the orbit
 * tail and underestimates the mean by a factor of several, so they only pass
 * with the residual-driven adaptive truncation.
 */
public class QsysBmapPhNnRetrialTest {

    private static final double LAMBDA = 0.8;
    private static final double MEAN_SERVICE = 1.0;
    private static final double NU = 1.0;
    private static final double RHO = LAMBDA * MEAN_SERVICE;

    /** Closed-form mean orbit length of the M/G/1 retrial queue. */
    private static double meanOrbit(double cB) {
        return RHO * RHO / (1.0 - RHO) * ((1.0 + cB * cB) / 2.0 + 1.0 / NU);
    }

    /** Poisson arrivals as a one-state BMAP {D0,D1}. */
    private static Matrix[] poisson(double lambda) {
        Matrix D0 = new Matrix(1, 1);
        D0.set(0, 0, -lambda);
        Matrix D1 = new Matrix(1, 1);
        D1.set(0, 0, lambda);
        return new Matrix[] { D0, D1 };
    }

    /**
     * Phase-type service with mean 1 and squared coefficient of variation cB^2:
     * an exponential for cB = 1, a balanced-mean two-phase hyperexponential
     * otherwise. Returns {beta, S}.
     */
    private static Matrix[] service(double cB) {
        if (cB == 1.0) {
            Matrix beta = new Matrix(1, 1);
            beta.set(0, 0, 1.0);
            Matrix S = new Matrix(1, 1);
            S.set(0, 0, -1.0 / MEAN_SERVICE);
            return new Matrix[] { beta, S };
        }
        double c2 = cB * cB;
        double s = Math.sqrt((c2 - 1.0) / (c2 + 1.0));
        double p1 = 0.5 * (1.0 + s);
        double p2 = 1.0 - p1;
        double mu1 = 2.0 * p1 / MEAN_SERVICE;
        double mu2 = 2.0 * p2 / MEAN_SERVICE;
        Matrix beta = new Matrix(1, 2);
        beta.set(0, 0, p1);
        beta.set(0, 1, p2);
        Matrix S = new Matrix(2, 2);
        S.set(0, 0, -mu1);
        S.set(1, 1, -mu2);
        return new Matrix[] { beta, S };
    }

    private static QsysRetrialResult solve(double cB) {
        Matrix[] D = poisson(LAMBDA);
        Matrix[] ph = service(cB);
        // N = 1 server, R = 0: an arrival finding the server busy joins the orbit.
        return Qsys_bmapphnn_retrial.qsys_bmapphnn_retrial(
                D, ph[0], ph[1], 1, NU, 0.0, 0.0, 0, -1, 1e-10, false);
    }

    @Test
    public void meanOrbitMatchesClosedFormExponentialService() {
        QsysRetrialResult res = solve(1.0);
        assertEquals(meanOrbit(1.0), res.L_orbit, 1e-3 * meanOrbit(1.0));
    }

    @Test
    public void meanOrbitMatchesClosedFormModerateVariability() {
        QsysRetrialResult res = solve(2.0);
        assertEquals(meanOrbit(2.0), res.L_orbit, 1e-3 * meanOrbit(2.0));
    }

    @Test
    public void meanOrbitMatchesClosedFormHighVariability() {
        QsysRetrialResult res = solve(5.0);
        assertEquals(meanOrbit(5.0), res.L_orbit, 1e-3 * meanOrbit(5.0));
    }

    @Test
    public void meanOrbitMatchesClosedFormVeryHighVariability() {
        QsysRetrialResult res = solve(10.0);
        assertEquals(meanOrbit(10.0), res.L_orbit, 1e-3 * meanOrbit(10.0));
    }

    /** Utilization of the single server is rho, and the orbit residual is met. */
    @Test
    public void utilizationAndResidual() {
        QsysRetrialResult res = solve(2.0);
        assertEquals(RHO, res.Utilization, 1e-6);
        assertEquals(LAMBDA, res.Throughput, 1e-6);
        assertTrue(res.truncError <= 1e-6, "truncation residual not met: " + res.truncError);
        assertTrue(res.truncLevel > 0);
    }

    /** A fixed truncation level is honoured verbatim. */
    @Test
    public void fixedTruncationLevelIsHonoured() {
        Matrix[] D = poisson(LAMBDA);
        Matrix[] ph = service(1.0);
        QsysRetrialResult res = Qsys_bmapphnn_retrial.qsys_bmapphnn_retrial(
                D, ph[0], ph[1], 1, NU, 0.0, 0.0, 0, 60, 1e-10, false);
        assertEquals(60, res.truncLevel);
        assertTrue(res.truncError > 0.0);
    }

    /** A non-Markovian (single-block) arrival representation is rejected. */
    @Test
    public void singleBlockArrivalIsRejected() {
        Matrix D0 = new Matrix(1, 1);
        D0.set(0, 0, -LAMBDA);
        final Matrix[] D = new Matrix[] { D0 };
        final Matrix[] ph = service(1.0);
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qsys_bmapphnn_retrial.qsys_bmapphnn_retrial(D, ph[0], ph[1], 1, NU, 0.0, 0.0, 0);
            }
        });
    }

    /** A NaN in the service subgenerator is rejected before the generator build. */
    @Test
    public void nanServiceIsRejected() {
        final Matrix[] D = poisson(LAMBDA);
        final Matrix beta = new Matrix(1, 1);
        beta.set(0, 0, 1.0);
        final Matrix S = new Matrix(1, 1);
        S.set(0, 0, Double.NaN);
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qsys_bmapphnn_retrial.qsys_bmapphnn_retrial(D, beta, S, 1, NU, 0.0, 0.0, 0);
            }
        });
    }

    /** Arrival blocks whose total does not have zero row sums are rejected. */
    @Test
    public void inconsistentBmapIsRejected() {
        Matrix D0 = new Matrix(1, 1);
        D0.set(0, 0, -LAMBDA);
        Matrix D1 = new Matrix(1, 1);
        D1.set(0, 0, 0.5 * LAMBDA);
        final Matrix[] D = new Matrix[] { D0, D1 };
        final Matrix[] ph = service(1.0);
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qsys_bmapphnn_retrial.qsys_bmapphnn_retrial(D, ph[0], ph[1], 1, NU, 0.0, 0.0, 0);
            }
        });
    }

    /** An out-of-range admission threshold is rejected. */
    @Test
    public void outOfRangeThresholdIsRejected() {
        final Matrix[] D = poisson(LAMBDA);
        final Matrix[] ph = service(1.0);
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qsys_bmapphnn_retrial.qsys_bmapphnn_retrial(D, ph[0], ph[1], 1, NU, 0.0, 0.0, 3);
            }
        });
    }

    /** An oversized per-level block is rejected rather than attempted. */
    @Test
    public void oversizedBlockIsRejected() {
        final Matrix[] D = poisson(LAMBDA);
        final Matrix[] ph = service(2.0);
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qsys_bmapphnn_retrial.qsys_bmapphnn_retrial(D, ph[0], ph[1], 1, NU, 0.0, 0.0, 0,
                        -1, 1e-10, false, 1e-6, 2e5, 2);
            }
        });
    }
}
