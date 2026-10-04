/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import java.util.Random;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class SamplesFromDMAP {
    private SamplesFromDMAP() {}

    /**
     * Generates random samples from a discrete Markovian arrival process.
     *
     * @param D0 The D0 matrix of the discrete MAP
     * @param D1 The D1 matrix of the discrete MAP
     * @param K The number of samples to generate
     * @param initial Optional initial state (1-indexed; null for random)
     * @param prec Numerical precision for validation
     * @param random Random number generator
     * @return Vector of random samples (inter-arrival times)
     */
    public static int[] samplesFromDMAP(Matrix D0, Matrix D1, int K, Integer initial, double prec, Random random) {
        if (!CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, prec)) {
            throw new IllegalArgumentException("SamplesFromDMAP: Input isn't a valid DMAP representation!");
        }

        MatrixCell D = new MatrixCell(2);
        D.set(0, D0);
        D.set(1, D1);

        return (int[]) SamplesFromDMMAP.samplesFromDMMAP(D, K, initial, prec, random);
    }

    public static int[] samplesFromDMAP(Matrix D0, Matrix D1, int K, Integer initial, double prec) {
        return samplesFromDMAP(D0, D1, K, initial, prec, new Random());
    }

    public static int[] samplesFromDMAP(Matrix D0, Matrix D1, int K, Integer initial) {
        return samplesFromDMAP(D0, D1, K, initial, 1e-14, new Random());
    }

    public static int[] samplesFromDMAP(Matrix D0, Matrix D1, int K) {
        return samplesFromDMAP(D0, D1, K, null, 1e-14, new Random());
    }
}
