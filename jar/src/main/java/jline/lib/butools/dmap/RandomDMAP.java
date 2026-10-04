/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import java.util.Random;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class RandomDMAP {
    private RandomDMAP() {}

    /**
     * Returns a random discrete Markovian arrival process.
     *
     * @param order The size of the DMAP
     * @param mean The mean inter-arrival time of the DMAP (default: 10.0)
     * @param zeroEntries The number of zero entries in the D0 and D1 matrices (default: 0)
     * @param maxTrials Maximum number of trials to find a proper DMAP (default: 1000)
     * @param prec Numerical precision for checking irreducibility (default: 1e-7)
     * @param random Random number generator
     * @return Pair of (D0, D1) matrices of the DMAP
     */
    public static Pair<Matrix, Matrix> randomDMAP(int order, double mean, int zeroEntries, int maxTrials, double prec, Random random) {
        MatrixCell D = RandomDMMAP.randomDMMAP(order, 1, mean, zeroEntries, maxTrials, prec, random);
        return new Pair<Matrix, Matrix>(D.get(0), D.get(1));
    }

    public static Pair<Matrix, Matrix> randomDMAP(int order, double mean, int zeroEntries, int maxTrials, double prec) {
        return randomDMAP(order, mean, zeroEntries, maxTrials, prec, new Random());
    }

    public static Pair<Matrix, Matrix> randomDMAP(int order, double mean, int zeroEntries, int maxTrials) {
        return randomDMAP(order, mean, zeroEntries, maxTrials, 1e-7, new Random());
    }

    public static Pair<Matrix, Matrix> randomDMAP(int order, double mean, int zeroEntries) {
        return randomDMAP(order, mean, zeroEntries, 1000, 1e-7, new Random());
    }

    public static Pair<Matrix, Matrix> randomDMAP(int order, double mean) {
        return randomDMAP(order, mean, 0, 1000, 1e-7, new Random());
    }

    public static Pair<Matrix, Matrix> randomDMAP(int order) {
        return randomDMAP(order, 10.0, 0, 1000, 1e-7, new Random());
    }
}
