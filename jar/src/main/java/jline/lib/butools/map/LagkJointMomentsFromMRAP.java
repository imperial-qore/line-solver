/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.mc.CRPSolve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class LagkJointMomentsFromMRAP {
    private LagkJointMomentsFromMRAP() {}

    public static MatrixCell lagkJointMomentsFromMRAP(MatrixCell H) {
        return lagkJointMomentsFromMRAP(H, 0, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromMRAP(MatrixCell H, int K) {
        return lagkJointMomentsFromMRAP(H, K, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromMRAP(MatrixCell H, int K, int L) {
        return lagkJointMomentsFromMRAP(H, K, L, 1e-14);
    }

    /**
     * Returns the lag-L joint moments of a continuous marked rational arrival process.
     */
    public static MatrixCell lagkJointMomentsFromMRAP(MatrixCell H, int K, int L, double prec) {
        if (!CheckMRAPRepresentation.checkMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException("LagkJointMomentsFromMRAP: Input isn't a valid MRAP representation!");
        }

        int M = H.size() - 1;
        int N = H.get(0).getNumRows();
        int actualK = (K == 0) ? N - 1 : K;

        // sumH = sum of H[1] ... H[M]
        Matrix sumH = new Matrix(N, N);
        for (int i = 1; i <= M; i++) {
            sumH = sumH.add(H.get(i));
        }

        // KEY CHANGE from DMRAP: iH0 = (-H0)^{-1} instead of (I-H0)^{-1}
        Matrix iH0 = H.get(0).neg().inv();

        // pi = DRPSolve(iH0 * sumH) -- stationary distribution of embedded DTMC
        Matrix P = iH0.mult(sumH);
        Matrix pi = CRPSolve.drpSolve(P);

        // Compute H0 powers
        Matrix[] H0p = new Matrix[actualK + 1];
        Matrix Pw = Matrix.eye(N);
        H0p[0] = Pw;

        Pw = Pw.mult(iH0);
        H0p[1] = Pw;

        for (int i = 2; i <= actualK; i++) {
            Pw = Pw.scale((double) i).mult(iH0);
            H0p[i] = Pw;
        }

        // Pl = (iH0*sumH)^(L-1)
        Matrix Pl = Matrix.eye(N);
        Matrix transitionMatrix = iH0.mult(sumH);
        for (int i = 0; i < L - 1; i++) {
            Pl = Pl.mult(transitionMatrix);
        }

        // Compute joint moments for each type
        MatrixCell Nm = new MatrixCell(M);
        for (int m = 0; m < M; m++) {
            Matrix Nmm = new Matrix(actualK + 1, actualK + 1);
            for (int i = 0; i <= actualK; i++) {
                for (int j = 0; j <= actualK; j++) {
                    Matrix temp = pi.mult(H0p[i]).mult(iH0).mult(H.get(m + 1)).mult(Pl).mult(H0p[j]);
                    Nmm.set(i, j, temp.elementSum());
                }
            }
            Nm.set(m, Nmm);
        }

        return Nm;
    }

    /**
     * Overload for Matrix[].
     */
    public static MatrixCell lagkJointMomentsFromMRAP(Matrix[] H, int K, int L, double prec) {
        MatrixCell cell = new MatrixCell(H.length);
        for (int i = 0; i < H.length; i++) {
            cell.set(i, H[i]);
        }
        return lagkJointMomentsFromMRAP(cell, K, L, prec);
    }
}
