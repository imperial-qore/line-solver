/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

/**
 * Phase-type representation of the first passage time into a target STATE SET.
 *
 * <p>This is the primitive behind the whole {@code ctmc_passage_*} family.
 * Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
 * in Large Markov Chains", 2002, Eqs. 1-2, which write the same system as n
 * scalar equations with L_i = 1 on the target.
 *
 * <p>THIS IS NOT THE SPLIT {@code SolverCTMC.getCdfRespT} USES. That one is by
 * EVENT (the tagged job arriving at or departing from a station, through the
 * filtration); this one is by STATE SET. The two are complementary and must not
 * be merged.
 *
 * <p>Target indices are 0-based here and 1-based in the MATLAB reference.
 */
public final class Ctmc_passage_ph {

    private Ctmc_passage_ph() {
    }

    /**
     * @param Q      generator, rows summing to zero
     * @param pi0    initial distribution over the state space; {@code null}
     *               starts from the conditional stationary law on the
     *               complement of the target set, as the reference and the
     *               python twin do
     * @param target 0-based target state indices
     */
    public static CtmcPassagePh ctmc_passage_ph(Matrix Q, Matrix pi0, int[] target) {
        int n = Q.getNumRows();
        if (Q.getNumCols() != n) {
            throw new RuntimeException("ctmc_passage_ph: the generator must be square");
        }
        int[] tgt = PassageSupport.uniqueTarget(target, n, "ctmc_passage_ph");

        double scale = 1.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                scale = Math.max(scale, Math.abs(Q.get(i, j)));
            }
        }
        for (int i = 0; i < n; i++) {
            double rs = 0.0;
            for (int j = 0; j < n; j++) {
                rs += Q.get(i, j);
            }
            if (Math.abs(rs) > 1e-8 * scale) {
                throw new RuntimeException("ctmc_passage_ph: Q is not an infinitesimal generator: "
                        + "its rows do not sum to zero. Pass it through ctmc_makeinfgen first");
            }
        }

        int[] keep = PassageSupport.complement(tgt, n);
        int nA = keep.length;
        Matrix S = new Matrix(nA, nA);
        for (int a = 0; a < nA; a++) {
            for (int c = 0; c < nA; c++) {
                S.set(a, c, Q.get(keep[a], keep[c]));
            }
        }
        Matrix s0 = new Matrix(nA, 1);
        for (int a = 0; a < nA; a++) {
            double r = 0.0;
            for (int c = 0; c < nA; c++) {
                r += S.get(a, c);
            }
            s0.set(a, 0, -r);
        }

        if (pi0 == null) {
            // An empty initial law selects the conditional stationary one on the
            // complement of the target set, the python api twin's contract.
            Matrix p = Ctmc_solve.ctmc_solve(Q);
            double mass = 0.0;
            for (int a = 0; a < nA; a++) {
                mass += p.get(keep[a]);
            }
            if (mass <= 0.0) {
                throw new RuntimeException("ctmc_passage_ph: the stationary law puts no mass outside "
                        + "the target set, so there is no passage to time");
            }
            Matrix alphaStat = new Matrix(1, nA);
            for (int a = 0; a < nA; a++) {
                alphaStat.set(0, a, p.get(keep[a]) / mass);
            }
            return new CtmcPassagePh(alphaStat, S, s0, keep, 0.0);
        }
        if (pi0.length() != n) {
            throw new RuntimeException("ctmc_passage_ph: pi0 must be a distribution over the "
                    + "state space, one entry per state");
        }
        Matrix alpha = new Matrix(1, nA);
        for (int a = 0; a < nA; a++) {
            alpha.set(0, a, pi0.get(keep[a]));
        }
        double atom = 0.0;
        for (int k : tgt) {
            atom += pi0.get(k);
        }
        return new CtmcPassagePh(alpha, S, s0, keep, atom);
    }
}
