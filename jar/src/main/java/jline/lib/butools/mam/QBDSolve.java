/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam;

import java.util.Map;

import jline.lib.butools.QBDFundamentalMatrices;
import jline.lib.butools.mc.DTMCSolve;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class QBDSolve {
    private QBDSolve() {}

    /**
     * Returns the parameters of the matrix-geometrically distributed stationary
     * distribution of a QBD.
     *
     * Using vector pi0 and matrix R provided by this function, the stationary
     * solution can be obtained by: pi_k = pi_0 * R^k
     *
     * @param B The matrix corresponding to backward transitions (N x N)
     * @param L The matrix corresponding to local transitions (N x N)
     * @param F The matrix corresponding to forward transitions (N x N)
     * @param L0 The matrix corresponding to local transitions at level zero (N x N)
     * @param prec The fundamental matrix R is computed up to this precision
     * @return Pair of (pi0, R) where pi0 is the stationary probability vector of level zero
     *         and R is the matrix parameter of the matrix geometrical distribution
     */
    public static Pair<Matrix, Matrix> qbdSolve(Matrix B, Matrix L, Matrix F, Matrix L0, double prec) {
        int m = L0.getNumRows();
        Matrix I = Matrix.eye(m);

        // Get fundamental matrix R
        Map<String, Matrix> fundMatrices =
                QBDFundamentalMatrices.QBDFundamentalMatrices(B, L, F, prec, null, "R", null);
        Matrix R = fundMatrices.get("R");

        // Convert to discrete time problem, if needed (continuous time has negative diagonal)
        Matrix Bdtmc = B;
        Matrix L0dtmc = L0;
        int negDiagCount = 0;
        for (int i = 0; i < m; i++) {
            if (L0.get(i, i) < 0) negDiagCount++;
        }
        if (negDiagCount > 0) {
            // Continuous time - uniformize
            double maxRate = 0.0;
            for (int i = 0; i < m; i++) {
                if (-L0.get(i, i) > maxRate) maxRate = -L0.get(i, i);
            }
            Bdtmc = B.scale(1.0 / maxRate);
            L0dtmc = L0.scale(1.0 / maxRate).add(I);
        }

        // pi0 = DTMCSolve(L0 + R*B) in discrete time
        Matrix transMatrix = L0dtmc.add(1.0, R.mult(Bdtmc));
        Matrix pi0 = DTMCSolve.dtmcSolve(transMatrix);

        // Normalize: pi0 * inv(I - R) * ones = 1
        Matrix IminusR = I.sub(R);
        Matrix IminusRinv = IminusR.inv();
        Matrix ones = Matrix.ones(m, 1);
        double nr = pi0.mult(IminusRinv).mult(ones).get(0, 0);
        if (nr > 0) {
            pi0 = pi0.scale(1.0 / nr);
        }

        return new Pair<Matrix, Matrix>(pi0, R);
    }

    public static Pair<Matrix, Matrix> qbdSolve(Matrix B, Matrix L, Matrix F, Matrix L0) {
        return qbdSolve(B, L, F, L0, 1e-14);
    }
}
