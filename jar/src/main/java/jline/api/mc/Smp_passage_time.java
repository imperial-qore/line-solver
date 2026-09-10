/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import java.util.List;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import jline.api.lti.Laplace_invert;
import jline.util.matrix.Matrix;

/**
 * Cumulative distribution and density of the SEMI-MARKOV first passage time, by
 * inverting {@link Smp_passage_lst} through {@code jline.api.lti}.
 *
 * <p>There is no matrix-exponential route here: a semi-Markov chain has no
 * generator to exponentiate, which is exactly the case uniformization does not
 * reach and the transform does. That is the argument the paper makes for
 * preferring transform inversion.
 *
 * <p>The inverter defaults to "euler" RATHER THAN "weeks". Semi-Markov passage
 * densities are the case Sec. 4.2 singles out as slow-converging for a Laguerre
 * series, and the Weeks scaling search then refuses by name rather than
 * returning noise.
 */
public final class Smp_passage_time {

    private Smp_passage_time() {
    }

    public static PassageCurve smp_passage_time(Matrix P, List<UnaryOperator<Complex>> hlst,
                                                Matrix pi0, int[] target, double[] tset) {
        return smp_passage_time(P, hlst, pi0, target, tset, "euler");
    }

    public static PassageCurve smp_passage_time(Matrix P, List<UnaryOperator<Complex>> hlst,
                                                Matrix pi0, int[] target, double[] tset,
                                                String ltiMethod) {
        int n = P.getNumRows();
        int[] tgt = PassageSupport.uniqueTarget(target, n, "smp_passage_time");
        double a = 0.0;
        if (pi0 != null && pi0.length() == n) {
            for (int k : tgt) {
                a += pi0.get(k);
            }
        }
        final double atom = a;
        final Matrix Pf = P;
        final List<UnaryOperator<Complex>> hf = hlst;
        final Matrix p0 = pi0;
        final int[] tf = target;

        UnaryOperator<Complex> L = new UnaryOperator<Complex>() {
            @Override
            public Complex apply(Complex s) {
                return Smp_passage_lst.smp_passage_lst(Pf, hf, p0, tf, s);
            }
        };
        UnaryOperator<Complex> Ld = new UnaryOperator<Complex>() {
            @Override
            public Complex apply(Complex s) {
                return L.apply(s).subtract(atom);
            }
        };
        double[] F = Laplace_invert.laplace_invert_cdf(L, tset, ltiMethod, 0);
        double[] f = Laplace_invert.laplace_invert_pdf(Ld, tset, ltiMethod, 0);
        return new PassageCurve(tset, F, f, atom);
    }
}
