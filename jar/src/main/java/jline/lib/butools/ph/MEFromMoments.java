/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.lib.butools.ReducedMomsFromMoms;
import jline.util.matrix.Matrix;

public final class MEFromMoments {
    private MEFromMoments() {}

    /**
     * Creates a matrix-exponential distribution that has the same moments as given.
     */
    public static MERepresentation meFromMoments(double[] moms) {
        Matrix K = appie(ReducedMomsFromMoms.ReducedMomsFromMoms(moms));
        int N = (int) Math.ceil(moms.length / 2.0);

        Matrix T = Matrix.zeros(N, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j <= i; j++) {
                T.set(i, j, 1.0);
            }
        }

        Matrix U = Matrix.zeros(N, N);
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                U.set(i, j, 1.0 / (N - i));
            }
        }

        Matrix alphaInit = Matrix.zeros(1, N);
        alphaInit.set(0, 0, 1.0);
        Matrix alpha = alphaInit.mult(T.inv()).mult(U);

        Matrix Uinv = U.inv();
        Matrix Tinv = T.inv();
        Matrix inner = Uinv.neg().mult(T).mult(K).mult(Tinv).mult(U);
        Matrix A = inner.inv();

        return new MERepresentation(alpha, A);
    }

    private static Matrix appie(double[] rmom) {
        int m = rmom.length;
        int actualM;
        double[] rm;

        if (m % 2 == 0) {
            rm = new double[m];
            rm[0] = 1.0;
            for (int i = 1; i < m; i++) rm[i] = rmom[i - 1];
            actualM = m / 2;
        } else {
            rm = new double[m + 1];
            rm[0] = 1.0;
            for (int i = 1; i < m + 1; i++) rm[i] = rmom[i - 1];
            actualM = (int) Math.ceil(m / 2.0);
        }

        int twoM = 2 * actualM;
        double[] f = new double[twoM];
        f[0] = 1.0;
        double[] y = new double[twoM];
        Matrix dd = Matrix.zeros(twoM, twoM);

        int n = 0;
        int k = 0;
        double q = 1.0;
        int[] d = new int[actualM];
        double[][] alphaArr = new double[actualM][actualM];
        double[] beta = new double[actualM];

        for (int i = 1; i < twoM; i++) {
            dd.set(i, i - 1, 1.0);
        }

        for (int i = 0; i < twoM; i++) {
            double ro = 0.0;
            int limit = Math.min(rm.length, f.length);
            for (int j = 0; j < limit; j++) {
                ro += q * rm[j] * f[j];
            }

            int nold = n;
            n = nold + 1;

            double[] yold = y.clone();

            if (n > 0 && ro != 0.0) {
                if (k > 0) {
                    int power = d[k - 1] + n - 1;
                    beta[k - 1] = ro / Math.pow(rm[1], (double) power);
                }
                k++;
                d[k - 1] = n;
                n = -n;
                q = q / ro;

                for (int j = 0; j < twoM; j++) {
                    y[j] = 0.0;
                    for (int l = 0; l < twoM; l++) {
                        y[j] += dd.get(j, l) * f[l];
                    }
                }
            } else if (n <= 0) {
                if (k > 0) {
                    int j = nold + d[k - 1] + 1;
                    if (j > 0 && j <= actualM) {
                        alphaArr[k - 1][j - 1] = ro / Math.pow(rm[1], (double) (j - 1));
                    }
                }
            }

            double[] newF = new double[twoM];
            for (int j = 0; j < twoM; j++) {
                for (int l = 0; l < twoM; l++) {
                    newF[j] += dd.get(j, l) * f[l];
                }
                newF[j] -= ro * yold[j];
            }
            for (int j = 0; j < twoM; j++) {
                f[j] = newF[j];
            }
        }

        int sumD = 0;
        for (int i = 0; i < actualM; i++) {
            sumD += d[i];
        }
        if (sumD != actualM) {
            throw new IllegalArgumentException("MEFromMoments: Insufficient matrix order!");
        }

        Matrix K = Matrix.zeros(actualM, actualM);
        K.set(0, 0, rm[1]);
        for (int i = 0; i < actualM - 1; i++) {
            K.set(i, i + 1, rm[1]);
        }

        int ind = d[0];
        for (int i = 1; i < actualM; i++) {
            if (ind < actualM) {
                int inc = d[i];
                ind += inc;
                if (ind <= actualM) {
                    K.set(ind - 1, ind - inc - d[i - 1], beta[i - 1]);
                    for (int j = 0; j < inc; j++) {
                        K.set(ind - 1, ind - j - 1, alphaArr[i][j]);
                    }
                }
            }
        }

        return K;
    }
}
