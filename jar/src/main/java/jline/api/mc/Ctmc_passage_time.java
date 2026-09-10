/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import jline.api.lti.Laplace_invert;
import jline.util.matrix.Matrix;

/**
 * Cumulative distribution and density of the first passage time into a target
 * state set: F(t) = 1 - alpha exp(St) 1 and f(t) = alpha exp(St) s0.
 *
 * <p>The method is "expm" (default) or "lt". The transform route exists for
 * chains whose non-target block is too large for a dense exp(St), NOT because it
 * needs fewer time points: every Abate-Whitt inverter places its nodes at
 * s = beta/t, so a grid of T points costs T*|beta| solves. On a small chain
 * "expm" is both faster and more accurate, which is why it is the default.
 *
 * <p>Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time
 * Distributions in Large Markov Chains", 2002.
 */
public final class Ctmc_passage_time {

    private Ctmc_passage_time() {
    }

    public static PassageCurve ctmc_passage_time(Matrix Q, Matrix pi0, int[] target,
                                                 double[] tset) {
        return ctmc_passage_time(Q, pi0, target, tset, "expm", "euler");
    }

    public static PassageCurve ctmc_passage_time(Matrix Q, Matrix pi0, int[] target, double[] tset,
                                                 String method, String ltiMethod) {
        CtmcPassagePh ph = Ctmc_passage_ph.ctmc_passage_ph(Q, pi0, target);
        int nA = ph.S.getNumRows();
        double[] F = new double[tset.length];
        double[] f = new double[tset.length];

        if ("expm".equalsIgnoreCase(method)) {
            double[] s0 = new double[nA];
            double[] a0 = new double[nA];
            for (int i = 0; i < nA; i++) {
                s0[i] = ph.s0.get(i, 0);
                a0[i] = ph.alpha.get(0, i);
            }

            boolean uniform = tset.length > 2;
            double dt = tset.length > 1 ? tset[1] - tset[0] : 0.0;
            if (uniform && !(dt > 0.0)) {
                uniform = false;
            }
            for (int i = 1; uniform && i + 1 < tset.length; i++) {
                if (Math.abs((tset[i + 1] - tset[i]) - dt) > 1e-12 * Math.max(1.0, Math.abs(dt))) {
                    uniform = false;
                }
            }

            if (uniform) {
                // One exponential, then propagate: recomputing expm(S*t) at
                // every grid point is the same answer at a cost linear in the
                // grid.
                Matrix E = ph.S.scale(dt).expm();
                Matrix E0 = ph.S.scale(tset[0]).expm();
                double[] v = rowTimes(a0, E0, nA);
                for (int i = 0; i < tset.length; i++) {
                    if (i > 0) {
                        v = rowTimes(v, E, nA);
                    }
                    double sF = 0.0;
                    double sf = 0.0;
                    for (int j = 0; j < nA; j++) {
                        sF += v[j];
                        sf += v[j] * s0[j];
                    }
                    F[i] = 1.0 - sF;
                    f[i] = sf;
                }
            } else {
                for (int i = 0; i < tset.length; i++) {
                    if (tset[i] < 0.0) {
                        continue;
                    }
                    Matrix E = ph.S.scale(tset[i]).expm();
                    double[] v = rowTimes(a0, E, nA);
                    double sF = 0.0;
                    double sf = 0.0;
                    for (int j = 0; j < nA; j++) {
                        sF += v[j];
                        sf += v[j] * s0[j];
                    }
                    F[i] = 1.0 - sF;
                    f[i] = sf;
                }
            }
        } else if ("lt".equalsIgnoreCase(method)) {
            final CtmcPassagePh phf = ph;
            UnaryOperator<Complex> L = new UnaryOperator<Complex>() {
                @Override
                public Complex apply(Complex s) {
                    return Ctmc_passage_lst.ctmc_passage_lst(phf, s);
                }
            };
            final double atom = ph.atom;
            UnaryOperator<Complex> Ld = new UnaryOperator<Complex>() {
                @Override
                public Complex apply(Complex s) {
                    return L.apply(s).subtract(atom);
                }
            };
            F = Laplace_invert.laplace_invert_cdf(L, tset, ltiMethod, 0);
            f = Laplace_invert.laplace_invert_pdf(Ld, tset, ltiMethod, 0);
        } else {
            throw new RuntimeException("ctmc_passage_time: unknown method '" + method
                    + "', expected expm or lt");
        }

        for (int i = 0; i < F.length; i++) {
            F[i] = Math.min(1.0, Math.max(0.0, F[i]));
            f[i] = Math.max(0.0, f[i]);
        }
        return new PassageCurve(tset, F, f, ph.atom);
    }

    private static double[] rowTimes(double[] v, Matrix E, int nA) {
        double[] out = new double[nA];
        for (int j = 0; j < nA; j++) {
            double acc = 0.0;
            for (int i = 0; i < nA; i++) {
                acc += v[i] * E.get(i, j);
            }
            out[j] = acc;
        }
        return out;
    }
}
