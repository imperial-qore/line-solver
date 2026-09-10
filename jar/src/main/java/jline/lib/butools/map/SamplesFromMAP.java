/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import java.util.Random;

import jline.lib.butools.mc.CRPSolve;
import jline.util.matrix.Matrix;

public final class SamplesFromMAP {
    private SamplesFromMAP() {}

    /**
     * Generates random samples from a Markovian arrival process.
     */
    public static double[] samplesFromMAP(Matrix D0, Matrix D1, int K, Integer initial, Random random) {
        int N = D0.getNumRows();

        Matrix P = D0.neg().inv().mult(D1);
        Matrix pi = CRPSolve.drpSolve(P);

        double[] cummInitial = new double[N];
        double cumSum = 0.0;
        for (int i = 0; i < N; i++) {
            cumSum += pi.get(0, i);
            cummInitial[i] = cumSum;
        }

        double[] sojourn = new double[N];
        for (int i = 0; i < N; i++) {
            sojourn[i] = -1.0 / D0.get(i, i);
        }

        double[][] nextprD0 = new double[N][N];
        for (int i = 0; i < N; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    double prob = sojourn[i] * D0.get(i, j);
                    rowSum += prob;
                    nextprD0[i][j] = rowSum;
                } else {
                    nextprD0[i][j] = rowSum;
                }
            }
        }

        double[] arrivalProb = new double[N];
        for (int i = 0; i < N; i++) {
            double d1sum = 0.0;
            for (int j = 0; j < N; j++) {
                d1sum += D1.get(i, j);
            }
            arrivalProb[i] = sojourn[i] * d1sum;
        }

        double[][] nextprD1 = new double[N][N];
        for (int i = 0; i < N; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < N; j++) {
                double prob = D1.get(i, j);
                rowSum += prob;
                nextprD1[i][j] = rowSum;
            }
            if (rowSum > 0) {
                for (int j = 0; j < N; j++) {
                    nextprD1[i][j] /= rowSum;
                }
            }
        }

        double[] samples = new double[K];

        int state;
        if (initial != null) {
            state = initial.intValue();
        } else {
            double r = random.nextDouble();
            int s = 0;
            while (s < N - 1 && cummInitial[s] <= r) {
                s++;
            }
            state = s;
        }

        for (int n = 0; n < K; n++) {
            double time = 0.0;

            while (true) {
                time -= Math.log(random.nextDouble()) * sojourn[state];
                double r = random.nextDouble();

                if (r < arrivalProb[state]) {
                    double r2 = random.nextDouble();
                    int nstate = 0;
                    while (nstate < N - 1 && nextprD1[state][nstate] <= r2) {
                        nstate++;
                    }
                    state = nstate;
                    break;
                } else {
                    double r2 = random.nextDouble();
                    int nstate = 0;
                    while (nstate < N - 1 && nextprD0[state][nstate] <= r2) {
                        nstate++;
                    }
                    state = nstate;
                }
            }
            samples[n] = time;
        }

        return samples;
    }

    public static double[] samplesFromMAP(Matrix D0, Matrix D1, int K, Integer initial) {
        return samplesFromMAP(D0, D1, K, initial, new Random());
    }

    public static double[] samplesFromMAP(Matrix D0, Matrix D1, int K) {
        return samplesFromMAP(D0, D1, K, null, new Random());
    }
}
