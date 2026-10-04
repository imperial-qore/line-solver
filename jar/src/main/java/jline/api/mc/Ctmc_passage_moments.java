/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

/**
 * Moments of order 1..nmax of the first passage time into a target state set.
 *
 * <p>This is Eq. 3 of P. G. Harrison and W. J. Knottenbelt, "Passage Time
 * Distributions in Large Markov Chains", 2002,
 *
 * <pre>
 *     -q_ii M_i(n) = sum_{k not in B} q_ik M_k(n) + n M_i(n-1),
 * </pre>
 *
 * i.e. (-S) M(n) = n M(n-1) with M(0) = 1: nmax linear solves and no transform
 * inversion at all. The equivalent closed form n! alpha (-S)^-n 1 is NOT how it
 * is evaluated here -- forming the inverse of the sub-generator destroys the
 * sparsity the recursion preserves.
 */
public final class Ctmc_passage_moments {

    private Ctmc_passage_moments() {
    }

    public static PassageMomentsResult ctmc_passage_moments(Matrix Q, Matrix pi0, int[] target,
                                                            int nmax) {
        if (nmax < 1) {
            throw new RuntimeException("ctmc_passage_moments: nmax must be positive");
        }
        CtmcPassagePh ph = Ctmc_passage_ph.ctmc_passage_ph(Q, pi0, target);
        int n = Q.getNumRows();
        int[] tgt = PassageSupport.uniqueTarget(target, n, "ctmc_passage_moments");
        int nA = ph.keep.length;

        Matrix mall = new Matrix(n, nmax);
        double[] m = new double[nmax];
        if (nA == 0) {
            return new PassageMomentsResult(mall, m);
        }

        // A state that cannot reach the target has an infinite passage time; the
        // sub-generator is singular on that block, and a solve that ignored this
        // would return a finite number instead of saying so.
        boolean[] reach = PassageSupport.reachesTarget(Q, ph.keep, tgt);
        List<Integer> idx = new ArrayList<Integer>();
        for (int i = 0; i < nA; i++) {
            if (reach[i]) {
                idx.add(i);
            }
        }
        int nR = idx.size();

        Matrix A = new Matrix(nR, nR);
        for (int a = 0; a < nR; a++) {
            for (int c = 0; c < nR; c++) {
                A.set(a, c, -ph.S.get(idx.get(a), idx.get(c)));
            }
        }

        double[] x = new double[nR];
        for (int i = 0; i < nR; i++) {
            x[i] = 1.0;
        }
        for (int k = 1; k <= nmax; k++) {
            Matrix rhs = new Matrix(nR, 1);
            for (int i = 0; i < nR; i++) {
                rhs.set(i, 0, k * x[i]);
            }
            Matrix sol = PassageSupport.solveReal(A, rhs, "ctmc_passage_moments");
            for (int i = 0; i < nR; i++) {
                x[i] = sol.get(i, 0);
            }
            for (int i = 0; i < nA; i++) {
                mall.set(ph.keep[i], k - 1, reach[i] ? 0.0 : Double.POSITIVE_INFINITY);
            }
            for (int i = 0; i < nR; i++) {
                mall.set(ph.keep[idx.get(i)], k - 1, x[i]);
            }
        }

        boolean unreachableStart = false;
        for (int i = 0; i < nA; i++) {
            if (!reach[i] && ph.alpha.get(0, i) > 0.0) {
                unreachableStart = true;
            }
        }
        for (int k = 0; k < nmax; k++) {
            if (unreachableStart) {
                m[k] = Double.POSITIVE_INFINITY;
                continue;
            }
            double acc = 0.0;
            for (int i = 0; i < nA; i++) {
                if (reach[i]) {
                    acc += ph.alpha.get(0, i) * mall.get(ph.keep[i], k);
                }
            }
            m[k] = acc;
        }
        return new PassageMomentsResult(mall, m);
    }
}
