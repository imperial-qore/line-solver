/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc;

import jline.util.matrix.Matrix;

public final class CRPSolve {
    private CRPSolve() {}

    /**
     * Computes the stationary solution of a continuous time
     * rational process (CRP).
     *
     * @param Q The generator matrix of the rational process
     * @param prec Numerical precision.
     * @return The vector that satisfies pi*Q = 0, sum(pi) = 1
     */
    public static Matrix crpSolve(Matrix Q, double prec) {
        int n = Q.getNumRows();

        // Check rowsums
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += Q.get(i, j);
            }
            if (Math.abs(rowSum) > prec) {
                throw new IllegalArgumentException("CRPSolve: The matrix has a rowsum which isn't zero!");
            }
        }

        // M = Q with first column replaced by ones
        Matrix M = Q.copy();
        for (int i = 0; i < n; i++) {
            M.set(i, 0, 1.0);
        }

        // m = [1, 0, 0, ...]
        Matrix m = Matrix.zeros(1, n);
        m.set(0, 0, 1.0);

        // pi = m / M = m * inv(M)
        Matrix pi = m.mult(M.inv());

        return pi;
    }

    public static Matrix crpSolve(Matrix Q) {
        return crpSolve(Q, 1e-14);
    }

    /**
     * Computes the stationary solution of a discrete time rational process (DRP).
     *
     * @param P The matrix parameter of the rational process
     * @return The vector that satisfies pi*P = pi, sum(pi) = 1
     */
    public static Matrix drpSolve(Matrix P) {
        int n = P.getNumRows();
        // Q = P - I
        Matrix Q = P.sub(Matrix.eye(n));
        return crpSolve(Q);
    }
}
