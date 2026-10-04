/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import java.util.List;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

/**
 * Laplace-Stieltjes transform of the first passage time into a target state set
 * for a SEMI-MARKOV chain.
 *
 * <p>Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time
 * Distributions in Large Markov Chains", 2002, Eqs. 4-5:
 *
 * <pre>
 *     L_i(s) = sum_{k not in B} r*_ik(s) L_k(s) + sum_{k in B} r*_ik(s),
 * </pre>
 *
 * so (I - R*_AA(s)) L_A(s) = R*_AB(s) 1, one linear system per value of s. The
 * sojourn transforms are given per state, so r*_ik(s) = P(i,k) h*_i(s) and the
 * complex numbers stay on the DIAGONAL of the system, which is the easier of
 * the two cases the paper distinguishes.
 *
 * <p>Distribution objects supply their own transform: {@code Markovian.evalLST}
 * gives the closed form pie (sI-D0)^-1 (-D0) e for the phase-type family, so
 * {@code s -> new Complex(dist.evalLST(s.getReal()), 0)} is NOT enough for a
 * complex argument -- pass a genuinely complex transform.
 */
public final class Smp_passage_lst {

    private Smp_passage_lst() {
    }

    /**
     * @param hlst one sojourn transform h*_i(s) per state, complex-argument
     */
    public static Complex smp_passage_lst(Matrix P, List<UnaryOperator<Complex>> hlst, Matrix pi0,
                                          int[] target, Complex s) {
        int n = P.getNumRows();
        int[] tgt = PassageSupport.uniqueTarget(target, n, "smp_passage_lst");
        if (hlst == null || hlst.size() != n) {
            throw new RuntimeException(
                    "smp_passage_lst: hlst must carry one transform per state");
        }
        int[] A = PassageSupport.complement(tgt, n);
        int nA = A.length;

        double atom = 0.0;
        boolean weighted = pi0 != null && pi0.length() == n;
        if (weighted) {
            for (int k : tgt) {
                atom += pi0.get(k);
            }
        }

        Complex[] h = new Complex[nA];
        for (int a = 0; a < nA; a++) {
            h[a] = hlst.get(A[a]).apply(s);
        }
        Complex[][] M = new Complex[nA][nA];
        Complex[] b = new Complex[nA];
        for (int a = 0; a < nA; a++) {
            for (int c = 0; c < nA; c++) {
                Complex diag = new Complex((a == c) ? 1.0 : 0.0, 0.0);
                M[a][c] = diag.subtract(h[a].multiply(P.get(A[a], A[c])));
            }
            double pb = 0.0;
            for (int k : tgt) {
                pb += P.get(A[a], k);
            }
            b[a] = h[a].multiply(pb);
        }
        Complex[] x = PassageSupport.solveComplex(M, b, "smp_passage_lst");
        Complex acc = new Complex(atom, 0.0);
        for (int a = 0; a < nA; a++) {
            double w = weighted ? pi0.get(A[a]) : 1.0 / n;
            acc = acc.add(x[a].multiply(w));
        }
        return acc;
    }
}
