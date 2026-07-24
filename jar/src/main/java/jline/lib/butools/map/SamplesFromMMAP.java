/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import java.util.Random;

import jline.lib.butools.ph.PHRepresentation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class SamplesFromMMAP {
    private SamplesFromMMAP() {}

    public static Object samplesFromMMAP(MatrixCell D, int K) {
        return samplesFromMMAP(D, K, null, 1e-14, new Random());
    }

    public static Object samplesFromMMAP(MatrixCell D, int K, Integer initial) {
        return samplesFromMMAP(D, K, initial, 1e-14, new Random());
    }

    public static Object samplesFromMMAP(MatrixCell D, int K, Integer initial, double prec) {
        return samplesFromMMAP(D, K, initial, prec, new Random());
    }

    /**
     * Generates random samples from a continuous marked Markovian arrival process.
     */
    public static Object samplesFromMMAP(MatrixCell D, int K, Integer initial, double prec, Random random) {
        if (!CheckMMAPRepresentation.checkMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("SamplesFromMMAP: Input isn't a valid MMAP representation!");
        }

        int N = D.get(0).getNumRows();
        int numTypes = D.size();

        int state;
        if (initial != null) {
            state = initial - 1;
        } else {
            PHRepresentation margDist =
                    MarginalDistributionFromMMAP.marginalDistributionFromMMAP(D, prec);
            Matrix stst = margDist.getAlpha();

            double[] cummInitial = new double[N];
            cummInitial[0] = stst.get(0, 0);
            for (int i = 1; i < N; i++) {
                cummInitial[i] = cummInitial[i - 1] + stst.get(0, i);
            }

            double r = random.nextDouble();
            int s = 0;
            while (s < N - 1 && cummInitial[s] <= r) {
                s++;
            }
            state = s;
        }

        int totalCols = N * numTypes;
        Matrix nextpr = new Matrix(N, totalCols);

        for (int i = 0; i < N; i++) {
            double rate = -D.get(0).get(i, i);
            if (rate <= 0) continue;
            double invRate = 1.0 / rate;

            int colOffset = 0;
            for (int dIdx = 0; dIdx < numTypes; dIdx++) {
                for (int j = 0; j < N; j++) {
                    double value;
                    if (dIdx == 0 && i == j) {
                        value = 0.0;
                    } else {
                        value = D.get(dIdx).get(i, j) * invRate;
                    }
                    nextpr.set(i, colOffset + j, value);
                }
                colOffset += N;
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 1; j < totalCols; j++) {
                nextpr.set(i, j, nextpr.get(i, j) + nextpr.get(i, j - 1));
            }
        }

        if (numTypes > 2) {
            double[][] samples = new double[K][2];
            for (int n = 0; n < K; n++) {
                double time = 0.0;
                while (state < N) {
                    double rate = -D.get(0).get(state, state);
                    time += -Math.log(random.nextDouble()) / rate;
                    double r = random.nextDouble();
                    int nstate = 0;
                    while (nstate < totalCols - 1 && nextpr.get(state, nstate) <= r) {
                        nstate++;
                    }
                    state = nstate;
                }
                samples[n][0] = time;
                samples[n][1] = (double) (state / N);
                state = state % N;
            }
            return samples;
        } else {
            double[] samples = new double[K];
            for (int n = 0; n < K; n++) {
                double time = 0.0;
                while (state < N) {
                    double rate = -D.get(0).get(state, state);
                    time += -Math.log(random.nextDouble()) / rate;
                    double r = random.nextDouble();
                    int nstate = 0;
                    while (nstate < totalCols - 1 && nextpr.get(state, nstate) <= r) {
                        nstate++;
                    }
                    state = nstate;
                }
                samples[n] = time;
                state = state % N;
            }
            return samples;
        }
    }

    public static Object samplesFromMMAP(Matrix[] D, int K) {
        return samplesFromMMAP(D, K, null, 1e-14, new Random());
    }

    public static Object samplesFromMMAP(Matrix[] D, int K, Integer initial, double prec, Random random) {
        MatrixCell cell = new MatrixCell(D.length);
        for (int i = 0; i < D.length; i++) {
            cell.set(i, D[i]);
        }
        return samplesFromMMAP(cell, K, initial, prec, random);
    }
}
