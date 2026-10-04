/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.mc.CRPSolve;
import jline.lib.butools.ph.MomentsFromME;
import jline.util.matrix.Matrix;

public final class LagCorrelationsFromMAP {
    private LagCorrelationsFromMAP() {}

    /**
     * Returns the lag autocorrelations of a Markovian arrival process.
     *
     * @param D0 The D0 matrix of the Markovian arrival process
     * @param D1 The D1 matrix of the Markovian arrival process
     * @param L The number of lags to compute. The default value is 1
     * @return The lag autocorrelation function up to lag L
     */
    public static double[] lagCorrelationsFromMAP(Matrix D0, Matrix D1, int L) {
        return lagCorrelationsFromRAP(D0, D1, L);
    }

    public static double[] lagCorrelationsFromMAP(Matrix D0, Matrix D1) {
        return lagCorrelationsFromRAP(D0, D1, 1);
    }

    /**
     * Returns the lag autocorrelations of a rational arrival process.
     *
     * @param H0 The H0 matrix of the rational arrival process
     * @param H1 The H1 matrix of the rational arrival process
     * @param L The number of lags to compute. The default value is 1
     * @return The lag autocorrelation function up to lag L
     */
    public static double[] lagCorrelationsFromRAP(Matrix H0, Matrix H1, int L) {
        Matrix H0i = H0.neg().inv();
        Matrix P = H0i.mult(H1);
        Matrix pi = CRPSolve.drpSolve(P);
        double[] moms = MomentsFromME.momentsFromME(pi, H0, 2);

        pi = pi.mult(H0i).mult(P);

        double[] acf = new double[L];
        for (int i = 0; i < L; i++) {
            double sum = pi.mult(H0i).elementSum();
            acf[i] = (sum - moms[0] * moms[0]) / (moms[1] - moms[0] * moms[0]);
            pi = pi.mult(P);
        }

        return acf;
    }

    public static double[] lagCorrelationsFromRAP(Matrix H0, Matrix H1) {
        return lagCorrelationsFromRAP(H0, H1, 1);
    }
}
