/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam;

import jline.util.matrix.Matrix;

public final class FluidStationaryDistr {
    private FluidStationaryDistr() {}

    /**
     * Returns the stationary distribution of a Markovian fluid model
     * at specific points.
     *
     * @param mass0 The stationary probability vector of zero level
     * @param ini The initial vector of the stationary density
     * @param K The matrix parameter of the stationary density
     * @param clo The closing matrix of the stationary density
     * @param x The points at which the distribution is evaluated
     * @return The stationary distribution at the specified points
     */
    public static Matrix fluidStationaryDistr(Matrix mass0, Matrix ini, Matrix K, Matrix clo, double[] x) {
        int N = mass0.getNumCols();
        int Np = K.getNumRows();
        Matrix result = Matrix.zeros(x.length, N);

        // closing = -K^{-1} * clo
        Matrix closing = K.neg().inv().mult(clo);
        Matrix I = Matrix.eye(Np);

        for (int i = 0; i < x.length; i++) {
            // y(x) = mass0 + ini * (I - e^(K*x)) * closing
            Matrix expKx = K.scale(x[i]).expm();
            Matrix IminusExp = I.sub(expKx);
            Matrix contrib = ini.mult(IminusExp).mult(closing);
            for (int j = 0; j < N; j++) {
                result.set(i, j, mass0.get(0, j) + contrib.get(0, j));
            }
        }

        return result;
    }
}
