/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam;

import jline.util.matrix.Matrix;

public final class QBDStationaryDistr {
    private QBDStationaryDistr() {}

    /**
     * Returns the stationary distribution of a QBD up to a given level K.
     *
     * @param pi0 The stationary probability vector of level zero
     * @param R The matrix parameter of the matrix geometrical distribution of the QBD
     * @param K The stationary distribution is returned up to this level
     * @return The stationary probability vector up to level K (length (K+1)*N)
     */
    public static Matrix qbdStationaryDistr(Matrix pi0, Matrix R, int K) {
        int m = R.getNumRows();
        Matrix pi = Matrix.zeros(1, (K + 1) * m);

        // Copy pi0 to first m elements
        for (int i = 0; i < m; i++) {
            pi.set(0, i, pi0.get(0, i));
        }

        Matrix pix = pi0.copy();
        for (int k = 1; k <= K; k++) {
            pix = pix.mult(R);
            for (int i = 0; i < m; i++) {
                pi.set(0, k * m + i, pix.get(0, i));
            }
        }

        return pi;
    }

    /**
     * Overload for double[] pi0.
     */
    public static Matrix qbdStationaryDistr(double[] pi0, Matrix R, int K) {
        return qbdStationaryDistr(new Matrix(pi0), R, K);
    }
}
