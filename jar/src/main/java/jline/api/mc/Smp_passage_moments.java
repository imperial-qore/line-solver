/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

/**
 * Moments of order 1..nmax of the first passage time into a target state set for
 * a SEMI-MARKOV chain with embedded transition matrix P.
 *
 * <p>Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time
 * Distributions in Large Markov Chains", 2002, Sec. 3.2. This class implements
 * Eq. 7 with the u_i(r) recurrence of Eq. 8,
 *
 * <pre>
 *     u_i(r) = -sum_{j=1..r} C(r,j) m_i(j) u_i(r-j),   u_i(0) = 1,
 * </pre>
 *
 * which are the derivatives at the origin of 1/h*_i(s). It is the branch for a
 * holding time that depends only on the CURRENT state, and is cheaper than the
 * full Eq. 6 kernel because it needs no per-pair moments.
 *
 * <p>Unlike the Markov case the n-th moment needs every moment from 1 to n, so
 * nmax cannot be raised for free.
 *
 * @see Ctmc_passage_moments the Markov case, which this reproduces exactly when
 *      the holding times are exponential
 */
public final class Smp_passage_moments {

    private Smp_passage_moments() {
    }

    /**
     * @param P     embedded transition matrix, rows summing to one
     * @param hmom  (nstates x nmax): hmom(i,r-1) is the r-th moment of the
     *              sojourn in state i
     * @param pi0   initial distribution, or {@code null} for per-state results
     *              only
     * @param target 0-based target state indices
     */
    public static PassageMomentsResult smp_passage_moments(Matrix P, Matrix hmom, Matrix pi0,
                                                           int[] target, int nmax) {
        int n = P.getNumRows();
        if (P.getNumCols() != n) {
            throw new RuntimeException(
                    "smp_passage_moments: the embedded transition matrix must be square");
        }
        if (nmax < 1) {
            throw new RuntimeException("smp_passage_moments: nmax must be positive");
        }
        for (int i = 0; i < n; i++) {
            double rs = 0.0;
            for (int j = 0; j < n; j++) {
                rs += P.get(i, j);
            }
            if (Math.abs(rs - 1.0) > 1e-8) {
                throw new RuntimeException(
                        "smp_passage_moments: the embedded transition matrix rows must sum to one");
            }
        }
        if (hmom.getNumRows() != n) {
            throw new RuntimeException("smp_passage_moments: hmom must carry one row per state");
        }
        if (hmom.getNumCols() < nmax) {
            throw new RuntimeException("smp_passage_moments: hmom must carry at least nmax "
                    + "holding-time moments per state");
        }
        int[] tgt = PassageSupport.uniqueTarget(target, n, "smp_passage_moments");
        int[] A = PassageSupport.complement(tgt, n);
        int nA = A.length;

        Matrix mall = new Matrix(n, nmax);
        double[] m = new double[nmax];
        if (nA == 0) {
            return new PassageMomentsResult(mall, m);
        }

        // Eq. 8, per state.
        Matrix u = new Matrix(nA, nmax);
        for (int r = 1; r <= nmax; r++) {
            for (int a = 0; a < nA; a++) {
                double acc = 0.0;
                for (int j = 1; j <= r; j++) {
                    double base = (r - j == 0) ? 1.0 : u.get(a, r - j - 1);
                    acc += PassageSupport.binom(r, j) * hmom.get(A[a], j - 1) * base;
                }
                u.set(a, r - 1, -acc);
            }
        }

        Matrix IPAA = new Matrix(nA, nA);
        for (int a = 0; a < nA; a++) {
            for (int c = 0; c < nA; c++) {
                IPAA.set(a, c, ((a == c) ? 1.0 : 0.0) - P.get(A[a], A[c]));
            }
        }

        Matrix M = new Matrix(nA, nmax);
        for (int q = 1; q <= nmax; q++) {
            Matrix b = new Matrix(nA, 1);
            for (int r = 1; r <= q; r++) {
                double c = PassageSupport.binom(q, r);
                for (int a = 0; a < nA; a++) {
                    double base = (r < q) ? M.get(a, q - r - 1) : 1.0;
                    b.set(a, 0, b.get(a, 0) - c * u.get(a, r - 1) * base);
                }
            }
            Matrix x = PassageSupport.solveReal(IPAA, b, "smp_passage_moments");
            for (int a = 0; a < nA; a++) {
                M.set(a, q - 1, x.get(a, 0));
            }
        }

        for (int a = 0; a < nA; a++) {
            for (int q = 0; q < nmax; q++) {
                mall.set(A[a], q, M.get(a, q));
            }
        }
        if (pi0 != null && pi0.length() == n) {
            for (int q = 0; q < nmax; q++) {
                double acc = 0.0;
                for (int i = 0; i < n; i++) {
                    acc += pi0.get(i) * mall.get(i, q);
                }
                m[q] = acc;
            }
        }
        return new PassageMomentsResult(mall, m);
    }
}
