/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import jline.util.matrix.Matrix;

public final class PmfFromMG {
    private PmfFromMG() {}

    /**
     * Returns the probability mass function of a matrix-geometric distribution.
     *
     * @param alpha The initial vector of the matrix-geometric distribution.
     * @param A The matrix parameter of the matrix-geometric distribution.
     * @param x Vector of non-negative integers at which to compute the PMF.
     * @return The probabilities that the matrix-geometrically distributed random variable
     *         takes the corresponding "x" values.
     */
    public static double[] pmfFromMG(Matrix alpha, Matrix A, int[] x) {
        int N = A.getNumRows();

        // a = 1 - sum(A, 2) = closing vector
        Matrix a = new Matrix(N, 1);
        for (int i = 0; i < N; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < N; j++) {
                rowSum += A.get(i, j);
            }
            a.set(i, 0, 1.0 - rowSum);
        }

        double[] pmf = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            if (x[i] == 0) {
                // PMF(0) = 1 - sum(alpha)
                pmf[i] = 1.0 - alpha.elementSum();
            } else {
                // PMF(x) = alpha * A^(x-1) * a
                Matrix APower = Matrix.pow(A, x[i] - 1);
                Matrix term = alpha.mult(APower).mult(a);
                pmf[i] = term.elementSum();
            }
        }

        return pmf;
    }

    public static double[] pmfFromMG(double[] alpha, Matrix A, int[] x) {
        return pmfFromMG(new Matrix(alpha), A, x);
    }

    public static double pmfFromMG(Matrix alpha, Matrix A, int x) {
        return pmfFromMG(alpha, A, new int[] {x})[0];
    }

    public static double pmfFromMG(double[] alpha, Matrix A, int x) {
        return pmfFromMG(new Matrix(alpha), A, new int[] {x})[0];
    }
}
