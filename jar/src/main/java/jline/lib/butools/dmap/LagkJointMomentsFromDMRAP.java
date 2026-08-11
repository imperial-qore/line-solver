/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.JMomsFromJFactorialMoms;
import jline.lib.butools.MomsFromFactorialMoms;
import jline.lib.butools.mc.CRPSolve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class LagkJointMomentsFromDMRAP {
    private LagkJointMomentsFromDMRAP() {}

    /**
     * Returns the lag-L joint moments of a discrete marked rational arrival process.
     */
    public static MatrixCell lagkJointMomentsFromDMRAP(MatrixCell H, int K, int L, double prec) {
        if (!CheckDMRAPRepresentation.checkDMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException("LagkJointMomentsFromDMRAP: Input isn't a valid DMRAP representation!");
        }

        int M = H.size() - 1;
        int N = H.get(0).getNumRows();
        int actualK = (K == 0) ? N - 1 : K;

        Matrix sumH = new Matrix(N, N);
        for (int i = 1; i <= M; i++) {
            sumH = sumH.add(H.get(i));
        }

        Matrix I = Matrix.eye(N);
        Matrix iH0 = I.sub(H.get(0)).inv();

        Matrix pi = CRPSolve.drpSolve(iH0.mult(sumH));

        Matrix[] H0p = new Matrix[actualK + 1];
        Matrix Pw = Matrix.eye(N);
        H0p[0] = Pw;

        Pw = Pw.mult(iH0);
        if (actualK >= 1) {
            H0p[1] = Pw;
        }

        for (int i = 2; i <= actualK; i++) {
            Pw = Pw.scale((double) i).mult(iH0).mult(H.get(0));
            H0p[i] = Pw;
        }

        Matrix Pl = Matrix.eye(N);
        Matrix transitionMatrix = iH0.mult(sumH);
        for (int i = 0; i < L - 1; i++) {
            Pl = Pl.mult(transitionMatrix);
        }

        MatrixCell Nm = new MatrixCell(M);
        for (int m = 0; m < M; m++) {
            Matrix Nmm = new Matrix(actualK + 1, actualK + 1);
            for (int i = 0; i <= actualK; i++) {
                for (int j = 0; j <= actualK; j++) {
                    Matrix temp = pi.mult(H0p[i]).mult(iH0).mult(H.get(m + 1)).mult(Pl).mult(H0p[j]);
                    Nmm.set(i, j, temp.elementSum());
                }
            }

            Matrix row1Input = new Matrix(1, actualK);
            for (int j = 0; j < actualK; j++) {
                row1Input.set(0, j, Nmm.get(0, j + 1));
            }
            Matrix row1 = MomsFromFactorialMoms.MomsFromFactorialMoms(row1Input);

            Matrix col1Input = new Matrix(actualK, 1);
            for (int j = 0; j < actualK; j++) {
                col1Input.set(j, 0, Nmm.get(j + 1, 0));
            }
            Matrix col1 = MomsFromFactorialMoms.MomsFromFactorialMoms(col1Input);

            Matrix midInput = new Matrix(actualK, actualK);
            for (int r = 0; r < actualK; r++) {
                for (int c = 0; c < actualK; c++) {
                    midInput.set(r, c, Nmm.get(r + 1, c + 1));
                }
            }
            Matrix mid = JMomsFromJFactorialMoms.jMomsFromJFactorialMoms(midInput);

            Matrix result = new Matrix(actualK + 1, actualK + 1);
            result.set(0, 0, Nmm.get(0, 0));
            for (int j = 0; j < actualK; j++) {
                result.set(0, j + 1, row1.get(0, j));
                result.set(j + 1, 0, col1.get(j, 0));
            }
            for (int r = 0; r < actualK; r++) {
                for (int c = 0; c < actualK; c++) {
                    result.set(r + 1, c + 1, mid.get(r, c));
                }
            }

            Nm.set(m, result);
        }

        return Nm;
    }

    public static MatrixCell lagkJointMomentsFromDMRAP(MatrixCell H, int K, int L) {
        return lagkJointMomentsFromDMRAP(H, K, L, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMRAP(MatrixCell H, int K) {
        return lagkJointMomentsFromDMRAP(H, K, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMRAP(MatrixCell H) {
        return lagkJointMomentsFromDMRAP(H, 0, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMRAP(Matrix[] H, int K, int L, double prec) {
        MatrixCell cell = new MatrixCell(H.length);
        for (int i = 0; i < H.length; i++) {
            cell.set(i, H[i]);
        }
        return lagkJointMomentsFromDMRAP(cell, K, L, prec);
    }

    public static MatrixCell lagkJointMomentsFromDMRAP(Matrix[] H, int K, int L) {
        return lagkJointMomentsFromDMRAP(H, K, L, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMRAP(Matrix[] H, int K) {
        return lagkJointMomentsFromDMRAP(H, K, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMRAP(Matrix[] H) {
        return lagkJointMomentsFromDMRAP(H, 0, 1, 1e-14);
    }
}
