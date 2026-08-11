/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc;

import jline.util.matrix.Matrix;

public final class CTMCSolve {
    private CTMCSolve() {}

    /**
     * Computes the stationary solution of a continuous time Markov chain.
     *
     * @param Q The generator matrix of the Markov chain.
     * @param prec Numerical precision.
     * @return The stationary probability vector (row vector).
     */
    public static Matrix ctmcSolve(Matrix Q, double prec) {
        int n = Q.getNumRows();

        // Create a copy of Q and modify the last column
        Matrix A = Q.copy();
        for (int i = 0; i < n; i++) {
            A.set(i, n - 1, 1.0);
        }

        // Create right-hand side vector
        Matrix b = Matrix.zeros(1, n);
        b.set(0, n - 1, 1.0);

        // Solve the linear system b = pi * A, i.e., pi = b * inv(A)
        Matrix pi = b.mult(A.inv());

        // Normalize to ensure sum is 1
        double sum = pi.elementSum();
        if (sum > 0) {
            for (int i = 0; i < n; i++) {
                pi.set(0, i, pi.get(0, i) / sum);
            }
        }

        return pi;
    }

    public static Matrix ctmcSolve(Matrix Q) {
        return ctmcSolve(Q, 1e-14);
    }
}
