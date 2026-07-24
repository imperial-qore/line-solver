/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import java.util.Random;

import jline.util.matrix.Matrix;

public final class SamplesFromDPH {
    private SamplesFromDPH() {}

    public static int[] samplesFromDPH(Matrix alpha, Matrix A, int K) {
        return samplesFromDPH(alpha, A, K, new Random());
    }

    /**
     * Generates random samples from a discrete phase-type distribution.
     */
    public static int[] samplesFromDPH(Matrix alpha, Matrix A, int K, Random random) {
        int N = alpha.length();

        // Compute cumulative initial distribution
        double[] cummInitial = new double[N];
        double sum = 0.0;
        for (int i = 0; i < N; i++) {
            sum += alpha.get(i);
            cummInitial[i] = sum;
        }

        // logp = log(diag(A))
        double[] logp = new double[N];
        for (int i = 0; i < N; i++) {
            logp[i] = (A.get(i, i) > 0) ? Math.log(A.get(i, i)) : Double.NEGATIVE_INFINITY;
        }

        // sojourn = 1 / (1 - diag(A))
        double[] sojourn = new double[N];
        for (int i = 0; i < N; i++) {
            sojourn[i] = 1.0 / (1.0 - A.get(i, i));
        }

        // nextpr = diag(sojourn) * A - diag(diag(nextpr))
        // Then add column for absorption: [nextpr, 1 - sum(nextpr, 2)]
        double[][] nextpr = new double[N][N + 1];
        for (int i = 0; i < N; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    nextpr[i][j] = sojourn[i] * A.get(i, j);
                    rowSum += nextpr[i][j];
                }
            }
            nextpr[i][N] = 1.0 - rowSum;
        }

        // Cumsum along rows
        double[][] nextprCum = new double[N][N + 1];
        for (int i = 0; i < N; i++) {
            double cumSum = 0.0;
            for (int j = 0; j <= N; j++) {
                cumSum += nextpr[i][j];
                nextprCum[i][j] = cumSum;
            }
        }

        // Generate samples
        int[] x = new int[K];
        for (int n = 0; n < K; n++) {
            int time = 0;

            // Draw initial state
            double r = random.nextDouble();
            int state = 0;
            while (state < N && cummInitial[state] <= r) {
                state++;
            }

            // Play state transitions until absorption
            while (state < N) {
                // Sojourn time in current state (geometric distribution)
                if (logp[state] > Double.NEGATIVE_INFINITY) {
                    int sojournTime = 1 + (int) Math.floor(Math.log(random.nextDouble()) / logp[state]);
                    time += sojournTime;
                } else {
                    time += 1;
                }

                // Next state transition
                double r2 = random.nextDouble();
                int nstate = 0;
                while (nstate <= N && nextprCum[state][nstate] <= r2) {
                    nstate++;
                }
                state = nstate;
            }

            x[n] = time;
        }

        return x;
    }

    /**
     * Overload for double[] alpha.
     */
    public static int[] samplesFromDPH(double[] alpha, Matrix A, int K, Random random) {
        return samplesFromDPH(new Matrix(alpha), A, K, random);
    }

    public static int[] samplesFromDPH(double[] alpha, Matrix A, int K) {
        return samplesFromDPH(new Matrix(alpha), A, K, new Random());
    }
}
