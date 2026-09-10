/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Transient sensitivity against the closed form of the two-state CTMC, which is
 * differentiable in the parameter by hand, and against a finite difference of
 * the transient solve itself.
 *
 * <p>For Q = [[-a, a], [b, -b]] and s = a + b,
 * pi_1(t) = b/s + (p0 - b/s) e^{-st}, so with theta = a,
 * dpi_1/da = -b/s^2 (1 - e^{-st}) - t (p0 - b/s) e^{-st} and dpi_2 = -dpi_1.
 */
public class CtmcTransientSensTest {

    private static Matrix gen(double a, double b) {
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -a);
        Q.set(0, 1, a);
        Q.set(1, 0, b);
        Q.set(1, 1, -b);
        return Q;
    }

    private static Matrix dGenDa() {
        Matrix dQ = new Matrix(2, 2);
        dQ.set(0, 0, -1.0);
        dQ.set(0, 1, 1.0);
        return dQ;
    }

    @Test
    public void twoStateMatchesTheClosedForm() {
        double a = 1.3;
        double b = 0.7;
        double p0 = 0.9;
        double s = a + b;
        Matrix pi0 = new Matrix(1, 2);
        pi0.set(0, 0, p0);
        pi0.set(0, 1, 1.0 - p0);
        Ctmc_transient_sens.Result r =
                Ctmc_transient_sens.ctmc_transient_sens(gen(a, b), dGenDa(), pi0, 0.0, 4.0);
        assertTrue(r.t.length > 2, "the integrator returned no interior grid");
        for (int k = 0; k < r.t.length; k++) {
            double t = r.t[k];
            double e = Math.exp(-s * t);
            double pi1 = b / s + (p0 - b / s) * e;
            double d1 = -b / (s * s) * (1.0 - e) - t * (p0 - b / s) * e;
            assertEquals(pi1, r.pi.get(k, 0), 1e-7, "pi1 at t=" + t);
            assertEquals(1.0 - pi1, r.pi.get(k, 1), 1e-7, "pi2 at t=" + t);
            assertEquals(d1, r.dpi.get(k, 0), 1e-6, "dpi1 at t=" + t);
            assertEquals(-d1, r.dpi.get(k, 1), 1e-6, "dpi2 at t=" + t);
        }
    }

    @Test
    public void sensitivityStartsAtZeroAndConservesMass() {
        Matrix pi0 = new Matrix(1, 2);
        pi0.set(0, 0, 0.25);
        pi0.set(0, 1, 0.75);
        Ctmc_transient_sens.Result r =
                Ctmc_transient_sens.ctmc_transient_sens(gen(2.0, 0.5), dGenDa(), pi0, 0.0, 3.0);
        assertEquals(0.0, r.dpi.get(0, 0), 1e-12);
        assertEquals(0.0, r.dpi.get(0, 1), 1e-12);
        for (int k = 0; k < r.t.length; k++) {
            // pi sums to one, so its derivative sums to zero
            assertEquals(1.0, r.pi.get(k, 0) + r.pi.get(k, 1), 1e-7);
            assertEquals(0.0, r.dpi.get(k, 0) + r.dpi.get(k, 1), 1e-7);
        }
    }

    @Test
    public void agreesWithACentralDifferenceOfTheTransientSolve() {
        double a = 0.9;
        double b = 1.4;
        double t1 = 2.0;
        double da = 1e-5;
        Matrix pi0 = new Matrix(1, 2);
        pi0.set(0, 0, 1.0);
        pi0.set(0, 1, 0.0);
        Ctmc_transient_sens.Result r =
                Ctmc_transient_sens.ctmc_transient_sens(gen(a, b), dGenDa(), pi0, 0.0, t1);
        int last = r.t.length - 1;
        double tEnd = r.t[last];
        double s;
        double fwd;
        double bwd;
        s = (a + da) + b;
        fwd = b / s + (1.0 - b / s) * Math.exp(-s * tEnd);
        s = (a - da) + b;
        bwd = b / s + (1.0 - b / s) * Math.exp(-s * tEnd);
        assertEquals((fwd - bwd) / (2 * da), r.dpi.get(last, 0), 1e-5);
    }
}
