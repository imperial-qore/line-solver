/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import java.util.Random;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class SamplesFromDMMAP {
    private SamplesFromDMMAP() {}

    public static Object samplesFromDMMAP(MatrixCell D, int K) {
        return samplesFromDMMAP(D, K, null, 1e-14, new Random());
    }

    public static Object samplesFromDMMAP(MatrixCell D, int K, Integer initial) {
        return samplesFromDMMAP(D, K, initial, 1e-14, new Random());
    }

    public static Object samplesFromDMMAP(MatrixCell D, int K, Integer initial, double prec) {
        return samplesFromDMMAP(D, K, initial, prec, new Random());
    }

    /**
     * Generates random samples from a discrete marked Markovian arrival process.
     */
    public static Object samplesFromDMMAP(MatrixCell D, int K, Integer initial, double prec, Random random) {
        if (!CheckDMMAPRepresentation.checkDMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("SamplesFromDMMAP: Input isn't a valid DMMAP representation!");
        }

        int N = D.get(0).getNumRows();
        int numTypes = D.size();

        int state;
        if (initial != null) {
            state = initial - 1;
        } else {
            MGRepresentation margDist =
                    MarginalDistributionFromDMMAP.marginalDistributionFromDMMAP(D, prec);
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

        double[] diag0 = new double[N];
        double[] sojourn = new double[N];
        double[] logp = new double[N];
        for (int i = 0; i < N; i++) {
            diag0[i] = D.get(0).get(i, i);
            sojourn[i] = 1.0 / (1.0 - diag0[i]);
            logp[i] = Math.log(diag0[i]);
        }

        int totalCols = N * numTypes;
        Matrix nextpr = new Matrix(N, totalCols);

        for (int i = 0; i < N; i++) {
            int colOffset = 0;
            for (int dIdx = 0; dIdx < numTypes; dIdx++) {
                for (int j = 0; j < N; j++) {
                    double value;
                    if (dIdx == 0 && i == j) {
                        value = 0.0;
                    } else {
                        value = sojourn[i] * D.get(dIdx).get(i, j);
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
            int[][] samples = new int[K][2];
            for (int n = 0; n < K; n++) {
                int time = 0;
                while (state < N) {
                    time += 1 + (int) Math.floor(Math.log(random.nextDouble()) / logp[state]);
                    double r = random.nextDouble();
                    int nstate = 0;
                    while (nstate < totalCols - 1 && nextpr.get(state, nstate) <= r) {
                        nstate++;
                    }
                    state = nstate;
                }
                samples[n][0] = time;
                samples[n][1] = state / N;
                state = state % N;
            }
            return samples;
        } else {
            int[] samples = new int[K];
            for (int n = 0; n < K; n++) {
                int time = 0;
                while (state < N) {
                    time += 1 + (int) Math.floor(Math.log(random.nextDouble()) / logp[state]);
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

    public static Object samplesFromDMMAP(Matrix[] D, int K) {
        return samplesFromDMMAP(D, K, null, 1e-14, new Random());
    }

    public static Object samplesFromDMMAP(Matrix[] D, int K, Integer initial, double prec, Random random) {
        MatrixCell cell = new MatrixCell(D.length);
        for (int i = 0; i < D.length; i++) {
            cell.set(i, D[i]);
        }
        return samplesFromDMMAP(cell, K, initial, prec, random);
    }
}
