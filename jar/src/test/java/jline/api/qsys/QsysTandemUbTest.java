/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

import org.junit.jupiter.api.Test;

import java.util.function.DoubleUnaryOperator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Ciucu-Mehri tail bounds for a two-station tandem.
 *
 * <p>The load-bearing check is the exact reduction: for M/M/1 -&gt; ./M/1 the five
 * inequalities of Lemma 4 hold as equalities, so the bound must return the exact
 * tails, (1 + theta x) e^{-theta x} for the sojourn time and the Kraemer form for
 * the waiting time. The remaining cases pin the digits produced by
 * matlab/src/api/qsys/qsys_tandem_ub_ciucu.m, which were themselves checked
 * against an exact CTMC reference for the Erlang(2)/M/1 tandem.
 */
public class QsysTandemUbTest {

    private static final double TOL = 1e-9;

    private static DoubleUnaryOperator expLst(double rate) {
        return s -> rate / (rate + s);
    }

    private static DoubleUnaryOperator expDlst(double rate) {
        return s -> rate / ((rate + s) * (rate + s));
    }

    private static DoubleUnaryOperator detLst(double d) {
        return s -> Math.exp(-s * d);
    }

    private static DoubleUnaryOperator detDlst(double d) {
        return s -> d * Math.exp(-s * d);
    }

    private static DoubleUnaryOperator erlang2Lst(double rate) {
        return s -> (rate / (rate + s)) * (rate / (rate + s));
    }

    private static DoubleUnaryOperator erlang2Dlst(double rate) {
        return s -> 2.0 * rate * rate / ((rate + s) * (rate + s) * (rate + s));
    }

    @Test
    public void testExactForMM1Tandem() {
        double mu = 1.0, lambda = 0.5, theta = 0.5;
        double[] x = {1.0, 5.0, 10.0};
        QsysTandemUbResult r = Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(
                x, expLst(lambda), new double[]{1.0}, new double[]{mu}, expDlst(lambda));

        assertEquals(theta, r.theta, 1e-12);
        assertEquals(1.0, r.A, TOL);
        assertEquals(0.0, r.B, 1e-12);
        assertEquals(0.0, r.D, 0.0);
        for (int k = 0; k < x.length; k++) {
            double S = (1.0 + theta * x[k]) * Math.exp(-theta * x[k]);
            double W = (1.0 - 2.0 * theta * theta / (mu * (mu + theta))
                    + x[k] * (mu - theta) * theta / (mu + theta)) * Math.exp(-theta * x[k]);
            assertEquals(S, r.S[k], TOL);
            assertEquals(W, r.W[k], TOL);
        }
    }

    @Test
    public void testDM1TandemMatchesMatlab() {
        double[] x = {5.0, 10.0, 20.0};
        double d = 4.0 / 3.0;                       // utilization 3/4, unit service rate
        QsysTandemUbResult r = Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(
                x, detLst(d), new double[]{1.0}, new double[]{1.0}, detDlst(d));

        assertEquals(0.454394983439251, r.theta, 1e-10);
        // Arrivals less variable than Poisson put the bound on the D > 0 branch.
        assertTrue(r.D > 0.0);
        assertEquals(0.493699069414895, r.S[0], TOL);
        assertEquals(0.09117765981667, r.S[1], TOL);
        assertEquals(0.00182565464670795, r.S[2], TOL);
        assertEquals(0.249922122549679, r.W[0], TOL);
    }

    @Test
    public void testErlangM1TandemMatchesMatlab() {
        double[] x = {10.0, 20.0};
        QsysTandemUbResult r = Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(
                x, erlang2Lst(1.0), new double[]{1.0}, new double[]{1.0}, erlang2Dlst(1.0));

        assertEquals(0.618033988749894, r.theta, 1e-10);
        assertEquals(0.0170463897138194, r.S[0], TOL);
        assertEquals(6.62788942098923e-05, r.S[1], TOL);
    }

    @Test
    public void testHyperexponentialServiceMatchesMatlab() {
        double p2 = 0.9, p1 = 0.1, m2 = 1.69;
        double m1 = p1 * m2 / (m2 - p2);            // CV(Y) = 2 with E[Y] = 1
        double meanY = p1 / m1 + p2 / m2;
        double rate = 2.0 / (meanY / 0.5);          // Erlang(2) arrivals at rho = 1/2
        double[] x = {10.0, 25.0, 50.0};
        QsysTandemUbResult r = Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(
                x, erlang2Lst(rate), new double[]{p1, p2}, new double[]{m1, m2},
                erlang2Dlst(rate));

        assertEquals(0.1500514416369, r.theta, 1e-10);
        assertEquals(0.524484692983612, r.S[0], TOL);
        assertEquals(0.0877479694886629, r.S[1], TOL);
        assertEquals(0.0033336271735292, r.S[2], TOL);
        for (int k = 0; k < x.length; k++) {
            assertTrue(Double.isNaN(r.W[k]));       // the W form is Exp-service only
        }
    }

    @Test
    public void testNumericalDerivativeMatchesAnalytic() {
        double[] x = {5.0, 10.0};
        double d = 4.0 / 3.0;
        QsysTandemUbResult ref = Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(
                x, detLst(d), new double[]{1.0}, new double[]{1.0}, detDlst(d));
        QsysTandemUbResult num = Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(
                x, detLst(d), new double[]{1.0}, new double[]{1.0});

        assertEquals(ref.alpha, num.alpha, 1e-10);
        for (int k = 0; k < x.length; k++) {
            assertEquals(ref.S[k], num.S[k], 1e-10);
        }
    }

    @Test
    public void testRejectsUnstableAndMalformedInput() {
        double[] x = {1.0};
        assertThrows(IllegalArgumentException.class, () ->
                Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(x, expLst(2.0),
                        new double[]{1.0}, new double[]{1.0}, expDlst(2.0)));
        assertThrows(IllegalArgumentException.class, () ->
                Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(new double[]{-1.0}, expLst(0.5),
                        new double[]{1.0}, new double[]{1.0}, expDlst(0.5)));
        assertThrows(IllegalArgumentException.class, () ->
                Qsys_tandem_ub_ciucu.qsys_tandem_ub_ciucu(x, expLst(0.5),
                        new double[]{0.3, 0.3}, new double[]{1.0, 2.0}));
    }
}
