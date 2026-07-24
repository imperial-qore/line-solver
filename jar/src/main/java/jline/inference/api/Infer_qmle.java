/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import jline.util.matrix.Matrix;

public final class Infer_qmle {
    private Infer_qmle() {}

    /**
     * Quick Maximum Likelihood Estimation closed-form formula.
     *
     * D(i,j) = Q(i,j) / (N(j) - sum(Q(:,j))) * Z(j) / (1 + sum(Q(i,:)) - Q(i,j)/N(j))
     *
     * @param Q mean queue lengths matrix (M x R)
     * @param N population vector (1 x R or length-R array)
     * @param Z think time vector (1 x R or length-R array)
     * @return demand estimates matrix (M x R)
     */
    public static Matrix infer_qmle(Matrix Q, double[] N, double[] Z) {
        int M = Q.getNumRows();
        int R = Q.getNumCols();
        Matrix D = new Matrix(M, R);

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                // sum(Q(:,j))
                double sumQcolJ = 0.0;
                for (int ii = 0; ii < M; ii++) {
                    sumQcolJ += Q.get(ii, j);
                }

                // sum(Q(i,:))
                double sumQrowI = 0.0;
                for (int jj = 0; jj < R; jj++) {
                    sumQrowI += Q.get(i, jj);
                }

                double denom1 = N[j] - sumQcolJ;
                double denom2 = 1.0 + sumQrowI - Q.get(i, j) / N[j];

                if (Math.abs(denom1) > 1e-14 && Math.abs(denom2) > 1e-14) {
                    D.set(i, j, Q.get(i, j) / denom1 * Z[j] / denom2);
                }
            }
        }

        return D;
    }
}
