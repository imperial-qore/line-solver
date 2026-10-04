/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import java.util.Random;

import jline.util.matrix.Matrix;

public final class SamplesFromPH {
    private SamplesFromPH() {}

    /**
     * Generates random samples from a phase-type distribution.
     *
     * @param alpha  The initial probability vector of the phase-type distribution.
     * @param A      The transient generator matrix of the phase-type distribution.
     * @param K      The number of samples to generate.
     * @param random Random number generator.
     * @return The vector of random samples.
     */
    public static double[] samplesFromPH(Matrix alpha, Matrix A, int K, Random random) {
        int N = alpha.length();

        double[] cummInitial = new double[N];
        double cumSum = 0.0;
        for (int i = 0; i < N; i++) {
            cumSum += alpha.get(i);
            cummInitial[i] = cumSum;
        }

        double[] sojourn = new double[N];
        for (int i = 0; i < N; i++) {
            sojourn[i] = -1.0 / A.get(i, i);
        }

        double[][] nextpr = new double[N][N + 1];
        for (int i = 0; i < N; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    double prob = sojourn[i] * A.get(i, j);
                    rowSum += prob;
                    nextpr[i][j] = rowSum;
                } else {
                    nextpr[i][j] = rowSum;
                }
            }
            nextpr[i][N] = 1.0;
        }

        double[] samples = new double[K];
        for (int n = 0; n < K; n++) {
            double time = 0.0;

            double r1 = random.nextDouble();
            int state = 0;
            while (state < N - 1 && cummInitial[state] <= r1) {
                state++;
            }

            while (state < N) {
                time -= Math.log(random.nextDouble()) * sojourn[state];
                double r2 = random.nextDouble();
                int nstate = 0;
                while (nstate < N && nextpr[state][nstate] <= r2) {
                    nstate++;
                }
                state = nstate;
            }
            samples[n] = time;
        }

        return samples;
    }

    public static double[] samplesFromPH(Matrix alpha, Matrix A, int K) {
        return samplesFromPH(alpha, A, K, new Random());
    }

    public static double[] samplesFromPH(double[] alpha, Matrix A, int K, Random random) {
        return samplesFromPH(new Matrix(alpha), A, K, random);
    }

    public static double[] samplesFromPH(double[] alpha, Matrix A, int K) {
        return samplesFromPH(new Matrix(alpha), A, K, new Random());
    }
}
