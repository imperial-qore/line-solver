/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

/**
 * Shared machinery for the passage-time family: target validation, the
 * reachability closure, and a complex linear solve.
 *
 * <p>Package-private on purpose. The public surface is one class per function,
 * as elsewhere in {@code jline.api}.
 */
final class PassageSupport {

    private PassageSupport() {
    }

    /** Sorted, de-duplicated and bounds-checked target indices (0-based). */
    static int[] uniqueTarget(int[] target, int n, String fn) {
        if (target == null || target.length == 0) {
            throw new RuntimeException(fn + ": the target state set is empty: a first passage "
                    + "time into no state is undefined");
        }
        int[] t = target.clone();
        Arrays.sort(t);
        int m = 1;
        for (int i = 1; i < t.length; i++) {
            if (t[i] != t[m - 1]) {
                t[m++] = t[i];
            }
        }
        t = Arrays.copyOf(t, m);
        if (t[0] < 0 || t[m - 1] >= n) {
            throw new RuntimeException(fn + ": a target state index is outside the state space");
        }
        return t;
    }

    /** The complement of the target set, in increasing order. */
    static int[] complement(int[] target, int n) {
        boolean[] isTarget = new boolean[n];
        for (int k : target) {
            isTarget[k] = true;
        }
        List<Integer> keep = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            if (!isTarget[i]) {
                keep.add(i);
            }
        }
        int[] out = new int[keep.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = keep.get(i);
        }
        return out;
    }

    /** Backward reachability closure over the transition graph. */
    static boolean[] reachesTarget(Matrix Q, int[] keep, int[] target) {
        int n = Q.getNumRows();
        boolean[] seen = new boolean[n];
        List<Integer> frontier = new ArrayList<Integer>();
        for (int k : target) {
            seen[k] = true;
            frontier.add(k);
        }
        while (!frontier.isEmpty()) {
            List<Integer> next = new ArrayList<Integer>();
            for (int j : frontier) {
                for (int i = 0; i < n; i++) {
                    if (i != j && !seen[i] && Q.get(i, j) != 0.0) {
                        seen[i] = true;
                        next.add(i);
                    }
                }
            }
            frontier = next;
        }
        boolean[] out = new boolean[keep.length];
        for (int a = 0; a < keep.length; a++) {
            out[a] = seen[keep[a]];
        }
        return out;
    }

    /**
     * Gaussian elimination with partial pivoting for a COMPLEX system.
     *
     * <p>{@code Matrix.solve} is real-valued, so the transform routes carry
     * their own solve rather than splitting each system into real and imaginary
     * blocks, which would double its dimension for no gain.
     */
    static Complex[] solveComplex(Complex[][] A, Complex[] b, String fn) {
        int n = b.length;
        Complex[][] M = new Complex[n][];
        for (int i = 0; i < n; i++) {
            M[i] = A[i].clone();
        }
        Complex[] r = b.clone();
        for (int k = 0; k < n; k++) {
            int piv = k;
            double best = M[k][k].abs();
            for (int i = k + 1; i < n; i++) {
                double mag = M[i][k].abs();
                if (mag > best) {
                    best = mag;
                    piv = i;
                }
            }
            if (!(best > 0.0)) {
                throw new RuntimeException(fn + ": singular transform matrix");
            }
            if (piv != k) {
                Complex[] tmp = M[k];
                M[k] = M[piv];
                M[piv] = tmp;
                Complex tb = r[k];
                r[k] = r[piv];
                r[piv] = tb;
            }
            for (int i = k + 1; i < n; i++) {
                Complex f = M[i][k].divide(M[k][k]);
                if (f.getReal() == 0.0 && f.getImaginary() == 0.0) {
                    continue;
                }
                for (int j = k; j < n; j++) {
                    M[i][j] = M[i][j].subtract(f.multiply(M[k][j]));
                }
                r[i] = r[i].subtract(f.multiply(r[k]));
            }
        }
        Complex[] x = new Complex[n];
        for (int i = n - 1; i >= 0; i--) {
            Complex s = r[i];
            for (int j = i + 1; j < n; j++) {
                s = s.subtract(M[i][j].multiply(x[j]));
            }
            x[i] = s.divide(M[i][i]);
        }
        return x;
    }

    /** Solve A x = b for a real square system, returning a column Matrix. */
    static Matrix solveReal(Matrix A, Matrix b, String fn) {
        Matrix x = new Matrix(A.getNumRows(), 1);
        if (!Matrix.solve(A, b, x)) {
            throw new RuntimeException(fn + ": the linear system is singular");
        }
        return x;
    }

    /** Binomial coefficient C(n,k) as a double, exact for the sizes used here. */
    static double binom(int n, int k) {
        double c = 1.0;
        for (int i = 0; i < k; i++) {
            c = c * (n - i) / (i + 1);
        }
        return c;
    }
}
