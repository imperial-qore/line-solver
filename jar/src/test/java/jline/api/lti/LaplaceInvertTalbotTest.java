/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.lti;

import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Talbot inversion on the api.lti surface, against closed forms and against the
 * Euler path that was previously the only alternative for a complex transform.
 */
public class LaplaceInvertTalbotTest {

    /** L{2 exp(-2t)} = 2/(s+2). */
    private static final UnaryOperator<Complex> EXP2 = new UnaryOperator<Complex>() {
        @Override
        public Complex apply(Complex s) {
            return new Complex(2.0, 0.0).divide(s.add(new Complex(2.0, 0.0)));
        }
    };

    /** L{t exp(-t)} = 1/(s+1)^2, an Erlang-2 density up to scale. */
    private static final UnaryOperator<Complex> ERL2 = new UnaryOperator<Complex>() {
        @Override
        public Complex apply(Complex s) {
            Complex d = s.add(new Complex(1.0, 0.0));
            return d.multiply(d).reciprocal();
        }
    };

    @Test
    public void talbotMatchesExponentialDensity() {
        double[] ts = {0.1, 0.5, 1.0, 2.0, 5.0};
        for (int i = 0; i < ts.length; i++) {
            double exact = 2.0 * Math.exp(-2.0 * ts[i]);
            double got = Laplace_invert.laplace_invert_talbot(EXP2, ts[i], 32);
            assertEquals(exact, got, 1e-9 * Math.max(1.0, Math.abs(exact)));
        }
    }

    @Test
    public void talbotMatchesErlangDensity() {
        double[] ts = {0.25, 1.0, 3.0};
        for (int i = 0; i < ts.length; i++) {
            double exact = ts[i] * Math.exp(-ts[i]);
            double got = Laplace_invert.laplace_invert_talbot(ERL2, ts[i], 32);
            assertEquals(exact, got, 1e-9);
        }
    }

    @Test
    public void dispatchAcceptsTalbotAndAgreesWithEuler() {
        double[] ts = {0.1, 1.0, 4.0};
        for (int i = 0; i < ts.length; i++) {
            double tal = Laplace_invert.laplace_invert(EXP2, ts[i], "talbot", 0);
            double eul = Laplace_invert.laplace_invert(EXP2, ts[i], "euler", 0);
            assertEquals(eul, tal, 1e-8 * Math.max(1.0, Math.abs(eul)));
        }
    }

    @Test
    public void talbotDrivesThePdfAndCdfGrids() {
        double[] ts = {0.0, 0.5, 1.0, 2.0, 4.0};
        double[] pdf = Laplace_invert.laplace_invert_pdf(EXP2, ts, "talbot", 0);
        double[] cdf = Laplace_invert.laplace_invert_cdf(EXP2, ts, "talbot", 0);
        assertEquals(0.0, pdf[0], 0.0);
        assertEquals(0.0, cdf[0], 0.0);
        for (int i = 1; i < ts.length; i++) {
            assertEquals(2.0 * Math.exp(-2.0 * ts[i]), pdf[i], 1e-9);
            assertEquals(1.0 - Math.exp(-2.0 * ts[i]), cdf[i], 1e-9);
            assertTrue(cdf[i] >= cdf[i - 1]);
        }
    }

    @Test
    public void namedWrappersAgreeWithTheDispatch() {
        double t = 1.5;
        assertEquals(Laplace_invert.laplace_invert(EXP2, t, "euler", 0),
                Laplace_invert.laplace_invert_euler(EXP2, t, 0), 0.0);
        assertEquals(Laplace_invert.laplace_invert(EXP2, t, "gaver-stehfest", 0),
                Laplace_invert.laplace_invert_gaver_stehfest(EXP2, t, 0), 0.0);
        assertEquals(Laplace_invert.laplace_invert(EXP2, t, "cme", 0),
                Laplace_invert.laplace_invert_cme(EXP2, t, 0), 0.0);
    }
}
