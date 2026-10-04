/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MomentsFromMG;
import jline.lib.butools.mc.CRPSolve;
import jline.util.matrix.Matrix;

public final class LagCorrelationsFromDRAP {
    private LagCorrelationsFromDRAP() {}

    /**
     * Returns the lag autocorrelations of a discrete rational arrival process.
     *
     * @param H0 The H0 matrix of the discrete rational arrival process
     * @param H1 The H1 matrix of the discrete rational arrival process
     * @param L The number of lags to compute. The default value is 1.
     * @param prec Numerical precision to check if the input is valid
     * @return The lag autocorrelation function up to lag L
     */
    public static double[] lagCorrelationsFromDRAP(Matrix H0, Matrix H1, int L, double prec) {
        if (!CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, prec)) {
            throw new IllegalArgumentException("LagCorrelationsFromDRAP: Input isn't a valid DRAP representation!");
        }

        int n = H0.getNumRows();
        Matrix I = Matrix.eye(n);

        // H0i = inv(I - H0)
        Matrix H0i = I.sub(H0).inv();

        // P = H0i * H1
        Matrix P = H0i.mult(H1);

        // pi = DRPSolve(P)
        Matrix pi = CRPSolve.drpSolve(P);

        // moms = MomentsFromMG(pi, H0, 2)
        double[] moms = MomentsFromMG.momentsFromMG(pi, H0, 2);

        // pi = pi * H0i * P
        pi = pi.mult(H0i).mult(P);

        // Compute lag correlations
        double[] acf = new double[L];
        for (int i = 0; i < L; i++) {
            Matrix piH0i = pi.mult(H0i);
            double sumPiH0i = piH0i.elementSum();
            acf[i] = (sumPiH0i - moms[0] * moms[0]) / (moms[1] - moms[0] * moms[0]);

            // pi = pi * P
            pi = pi.mult(P);
        }

        return acf;
    }

    public static double[] lagCorrelationsFromDRAP(Matrix H0, Matrix H1, int L) {
        return lagCorrelationsFromDRAP(H0, H1, L, 1e-14);
    }

    public static double[] lagCorrelationsFromDRAP(Matrix H0, Matrix H1) {
        return lagCorrelationsFromDRAP(H0, H1, 1, 1e-14);
    }
}
