/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam;

import java.util.Arrays;
import java.util.List;

import jline.lib.butools.mc.DTMCSolve;
import jline.util.matrix.Matrix;

public final class GM1StationaryDistr {
    private GM1StationaryDistr() {}

    /**
     * Returns the stationary distribution of the G/M/1 type Markov chain
     * up to a given level K.
     *
     * @param B List of matrix blocks of the G/M/1 type generator at the boundary.
     * @param R Matrix R of the G/M/1 type Markov chain.
     * @param K The stationary distribution is returned up to this level.
     * @return The stationary probability vector up to level K.
     */
    public static Matrix gm1StationaryDistr(List<Matrix> B, Matrix R, int K) {
        int m = R.getNumRows();
        Matrix I = Matrix.eye(m);

        // Check that spectral radius of R is below 1
        Matrix IminusRinv = I.sub(R).inv();

        int maxb = B.size();

        // Compute BR = sum_{i} R^i * B[i] using Horner's method in reverse
        Matrix BR = B.get(maxb - 1).copy();
        for (int i = maxb - 1; i >= 1; i--) {
            BR = R.mult(BR).add(1.0, B.get(i - 1));
        }

        // Solve pi0 as the stationary vector of BR
        Matrix pix = DTMCSolve.dtmcSolve(BR);

        // Normalize: pi0 * inv(I - R) * ones = 1
        Matrix ones = Matrix.ones(m, 1);
        double nr = pix.mult(IminusRinv).mult(ones).get(0, 0);
        if (nr > 0) {
            pix = pix.scale(1.0 / nr);
        }

        // Build the distribution up to level K
        Matrix pi = Matrix.zeros(1, (K + 1) * m);
        for (int j = 0; j < m; j++) {
            pi.set(0, j, pix.get(0, j));
        }

        double sumpi = pix.elementSum();
        int numit = 1;
        while (sumpi < 1.0 - 1e-10 && numit <= K) {
            pix = pix.mult(R);
            sumpi += pix.elementSum();
            for (int j = 0; j < m; j++) {
                pi.set(0, numit * m + j, pix.get(0, j));
            }
            numit++;
        }

        return pi;
    }

    /**
     * Overload accepting Matrix[].
     */
    public static Matrix gm1StationaryDistr(Matrix[] B, Matrix R, int K) {
        return gm1StationaryDistr(Arrays.asList(B), R, K);
    }
}
