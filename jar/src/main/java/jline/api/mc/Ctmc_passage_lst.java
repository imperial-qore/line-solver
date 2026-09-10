/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

/**
 * Laplace-Stieltjes transform of the first passage time into a target state
 * set: L(s) = alpha (sI-S)^-1 s0 + atom.
 *
 * <p>Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time
 * Distributions in Large Markov Chains", 2002, Eqs. 1-2: one linear system per
 * value of s, of the size of the non-target block.
 *
 * <p>ONE SOLVE PER s, NOT PER (s,t) PAIR. The saving over a dense matrix
 * exponential is that the solves are sparse, so this route reaches chains a
 * dense expm cannot hold. It is NOT a saving in the number of time points:
 * every Abate-Whitt inverter places its nodes at s = beta/t, so a grid of T
 * points costs T*|beta| solves. On a small chain the exponential route of
 * {@code Ctmc_passage_time} is faster.
 */
public final class Ctmc_passage_lst {

    private Ctmc_passage_lst() {
    }

    /** L(s) at a single (possibly complex) point. */
    public static Complex ctmc_passage_lst(Matrix Q, Matrix pi0, int[] target, Complex s) {
        return ctmc_passage_lst(Ctmc_passage_ph.ctmc_passage_ph(Q, pi0, target), s);
    }

    /** L(s) from a phase-type form already built, so it is not rebuilt per s. */
    public static Complex ctmc_passage_lst(CtmcPassagePh ph, Complex s) {
        int nA = ph.S.getNumRows();
        Complex[][] A = new Complex[nA][nA];
        Complex[] b = new Complex[nA];
        for (int i = 0; i < nA; i++) {
            for (int j = 0; j < nA; j++) {
                double d = (i == j) ? 1.0 : 0.0;
                A[i][j] = s.multiply(d).subtract(new Complex(ph.S.get(i, j), 0.0));
            }
            b[i] = new Complex(ph.s0.get(i, 0), 0.0);
        }
        Complex[] x = PassageSupport.solveComplex(A, b, "ctmc_passage_lst");
        Complex acc = new Complex(ph.atom, 0.0);
        for (int i = 0; i < nA; i++) {
            acc = acc.add(x[i].multiply(ph.alpha.get(0, i)));
        }
        return acc;
    }

    /** L(s) on a vector of points. */
    public static Complex[] ctmc_passage_lst(Matrix Q, Matrix pi0, int[] target, Complex[] s) {
        CtmcPassagePh ph = Ctmc_passage_ph.ctmc_passage_ph(Q, pi0, target);
        Complex[] out = new Complex[s.length];
        for (int i = 0; i < s.length; i++) {
            out[i] = ctmc_passage_lst(ph, s[i]);
        }
        return out;
    }
}
