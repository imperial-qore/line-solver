/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc;

import jline.util.matrix.Matrix;

public final class DTMCSolve {
    private DTMCSolve() {}

    /**
     * Computes the stationary solution of a discrete time
     * Markov chain.
     *
     * @param P The transition probability matrix of the Markov chain.
     * @param prec Numerical precision. The default value is 1e-14.
     * @return The stationary probability vector (row vector).
     */
    public static Matrix dtmcSolve(Matrix P, double prec) {
        int n = P.getNumRows();

        // Compute Q = P - I
        Matrix Q = P.sub(Matrix.eye(n));

        // Use CTMC solve on the embedded generator
        return CTMCSolve.ctmcSolve(Q, prec);
    }

    public static Matrix dtmcSolve(Matrix P) {
        return dtmcSolve(P, 1e-14);
    }
}
