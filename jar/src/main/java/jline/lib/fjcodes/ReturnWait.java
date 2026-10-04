package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

/**
 * Compute waiting time distribution for Fork-Join queue.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class ReturnWait {
    private ReturnWait() {}

    /**
     * Compute Phase-Type representation of stationary waiting time.
     */
    public static ReturnWaitResult returnWait(double En1, Matrix pi0, Matrix T, Matrix phi, Matrix sum_Ajump) {
        // alfa = -pi0 / T
        Matrix alfa = pi0.scale(-1.0).rightMatrixDivide(T);

        int ds = phi.getNumRows();

        Matrix phiAlfa = phi.mult(alfa);
        Matrix alfaPhi = alfa.mult(phi);
        double alfaPhiSum = alfaPhi.get(0, 0);

        Matrix rhos = new Matrix(1, ds);
        for (int i = 0; i < ds; i++) {
            rhos.set(0, i, phiAlfa.get(i, 0) / alfaPhiSum);
        }

        double alfaSumAjump = alfa.mult(sum_Ajump).get(0, 0);
        double sumPi0 = elementSum(pi0);
        double En0 = alfaSumAjump / sumPi0;

        double prob_wait = (En0 - 1.0) / (En0 - 1.0 + En1);

        Matrix wait_alpha = rhos.scale(prob_wait);

        Matrix wait_Smat = new Matrix(ds, ds);
        for (int i = 0; i < ds; i++) {
            for (int j = 0; j < ds; j++) {
                double gij = alfa.get(0, j) * T.get(j, i) / alfa.get(0, i);
                wait_Smat.set(i, j, gij);
            }
        }

        return new ReturnWaitResult(wait_alpha, wait_Smat, prob_wait, alfa);
    }

    private static double elementSum(Matrix m) {
        double sum = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                sum += m.get(i, j);
            }
        }
        return sum;
    }
}
