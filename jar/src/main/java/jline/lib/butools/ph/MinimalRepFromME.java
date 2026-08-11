/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import java.util.ArrayList;
import java.util.List;

import jline.util.Pair;

import jline.lib.butools.reptrans.MStaircase;
import jline.util.matrix.Matrix;

public final class MinimalRepFromME {
    private MinimalRepFromME() {}

    /**
     * Returns the minimal representation of the given ME distribution.
     */
    public static MERepresentation minimalRepFromME(Matrix alpha, Matrix A, String how, double prec) {
        int N = A.getNumRows();

        if ("cont".equals(how)) {
            Matrix H0 = A;
            Matrix negRowSums = new Matrix(N, 1);
            for (int i = 0; i < N; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < N; j++) {
                    rowSum += A.get(i, j);
                }
                negRowSums.set(i, 0, -rowSum);
            }
            Matrix H1 = negRowSums.mult(alpha);

            List<Matrix> matrices = new ArrayList<Matrix>();
            matrices.add(H0);
            matrices.add(H1);
            Matrix ones = Matrix.ones(N, 1);
            Pair<Matrix, Integer> result = MStaircase.mStaircase(matrices, ones, prec);
            Matrix B = result.getFirst();
            int n = result.getSecond();

            if (n < N) {
                Matrix Binv = B.inv();
                Matrix alphaB = alpha.mult(B);
                Matrix newAlpha = new Matrix(1, n);
                for (int j = 0; j < n; j++) {
                    newAlpha.set(0, j, alphaB.get(0, j));
                }
                Matrix BAB = Binv.mult(A).mult(B);
                Matrix newA = new Matrix(n, n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        newA.set(i, j, BAB.get(i, j));
                    }
                }
                return new MERepresentation(newAlpha, newA);
            }
            return new MERepresentation(alpha, A);
        } else if ("obs".equals(how)) {
            Matrix H0 = A;
            Matrix negRowSums = new Matrix(N, 1);
            for (int i = 0; i < N; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < N; j++) {
                    rowSum += A.get(i, j);
                }
                negRowSums.set(i, 0, -rowSum);
            }
            Matrix H1 = negRowSums.mult(alpha);

            List<Matrix> matrices = new ArrayList<Matrix>();
            matrices.add(H0.transpose());
            matrices.add(H1.transpose());
            Matrix alphaT = alpha.transpose();
            Pair<Matrix, Integer> result = MStaircase.mStaircase(matrices, alphaT, prec);
            Matrix B = result.getFirst();
            int n = result.getSecond();

            if (n < N) {
                Matrix Binv = B.inv();
                Matrix alphaB = alpha.mult(B);
                Matrix newAlpha = new Matrix(1, n);
                for (int j = 0; j < n; j++) {
                    newAlpha.set(0, j, alphaB.get(0, j));
                }
                Matrix BAB = Binv.mult(A).mult(B);
                Matrix newA = new Matrix(n, n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        newA.set(i, j, BAB.get(i, j));
                    }
                }
                return new MERepresentation(newAlpha, newA);
            }
            return new MERepresentation(alpha, A);
        } else if ("obscont".equals(how)) {
            MERepresentation rep1 = minimalRepFromME(alpha, A, "cont", prec);
            return minimalRepFromME(rep1.getAlpha(), rep1.getA(), "obs", prec);
        } else if ("moment".equals(how)) {
            int order = MEOrder.meOrder(alpha, A, "moment", prec);
            if (order < N) {
                double[] momsArr = MomentsFromME.momentsFromME(alpha, A, 2 * order - 1);
                return MEFromMoments.meFromMoments(momsArr);
            }
            return new MERepresentation(alpha, A);
        } else {
            throw new IllegalArgumentException("MinimalRepFromME: Unknown method '" + how + "'!");
        }
    }

    public static MERepresentation minimalRepFromME(Matrix alpha, Matrix A) {
        return minimalRepFromME(alpha, A, "obscont", 1e-12);
    }

    public static MERepresentation minimalRepFromME(Matrix alpha, Matrix A, String how) {
        return minimalRepFromME(alpha, A, how, 1e-12);
    }
}
