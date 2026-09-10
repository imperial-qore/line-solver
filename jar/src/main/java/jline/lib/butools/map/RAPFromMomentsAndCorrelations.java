/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.MomsFromReducedMoms;
import jline.lib.butools.ph.MEFromMoments;
import jline.lib.butools.ph.MERepresentation;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class RAPFromMomentsAndCorrelations {
    private RAPFromMomentsAndCorrelations() {}

    /**
     * Returns a rational arrival process that has the same moments
     * and lag autocorrelation coefficients as given.
     */
    public static Pair<Matrix, Matrix> rapFromMomentsAndCorrelations(double[] moms, double[] corr) {
        MERepresentation meResult = MEFromMoments.meFromMoments(moms);
        Matrix alpha = meResult.getAlpha();
        Matrix D0 = meResult.getA();

        int M = alpha.getNumCols();

        if (corr.length < 2 * M - 3) {
            throw new IllegalArgumentException(
                "RAPFromMomentsAndCorrelations: The number of correlations given is less than the required 2n-3!");
        }

        double corrFactor = (moms[1] / 2.0 - moms[0] * moms[0]) / (moms[1] - moms[0] * moms[0]);
        int numCorr = 2 * M - 3;
        double[] rcorrArray = new double[numCorr];
        for (int i = 0; i < numCorr; i++) {
            rcorrArray[i] = corr[i] / corrFactor;
        }

        Matrix rcorrMatrix = new Matrix(1, numCorr);
        for (int i = 0; i < numCorr; i++) {
            rcorrMatrix.set(0, i, rcorrArray[i]);
        }
        Matrix rcorrMoms = MomsFromReducedMoms.momsFromReducedMoms(rcorrMatrix);
        double[] rcorrMomsArray = new double[numCorr];
        for (int i = 0; i < numCorr; i++) {
            rcorrMomsArray[i] = rcorrMoms.get(0, i);
        }

        MERepresentation meResult2 = MEFromMoments.meFromMoments(rcorrMomsArray);
        Matrix X = meResult2.getA();

        int NN = X.getNumRows();

        if (NN + 1 != D0.getNumRows()) {
            throw new IllegalArgumentException(
                "RAPFromMomentsAndCorrelations: Correlation order is different from ME order!");
        }

        Matrix T1 = Matrix.zeros(NN, NN);
        for (int i = 0; i < NN; i++) {
            for (int j = 0; j <= i; j++) {
                T1.set(i, j, 1.0);
            }
        }

        Matrix U1 = Matrix.zeros(NN, NN);
        for (int i = 0; i < NN; i++) {
            for (int j = i; j < NN; j++) {
                U1.set(i, j, 1.0 / (NN - i));
            }
        }

        Matrix T2 = Matrix.zeros(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j <= i; j++) {
                T2.set(i, j, 1.0);
            }
        }

        Matrix U2 = Matrix.zeros(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = i; j < M; j++) {
                U2.set(i, j, 1.0 / (M - i));
            }
        }

        Matrix Y = T1.inv().neg().mult(U1).mult(X.inv()).mult(U1.inv()).mult(T1);

        Matrix II = Matrix.eye(NN + 1);
        for (int i = 0; i < NN; i++) {
            for (int j = 0; j < NN; j++) {
                II.set(i + 1, j + 1, Y.get(i, j));
            }
        }

        Matrix D1 = D0.neg().mult(U2.inv()).mult(T2).mult(II).mult(T2.inv()).mult(U2);

        return new Pair<Matrix, Matrix>(D0, D1);
    }
}
