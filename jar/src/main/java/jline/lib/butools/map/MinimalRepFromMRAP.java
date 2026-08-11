/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import java.util.ArrayList;
import java.util.List;

import jline.lib.butools.ph.PHRepresentation;
import jline.lib.butools.reptrans.MStaircase;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MinimalRepFromMRAP {
    private MinimalRepFromMRAP() {}

    public static MatrixCell minimalRepFromMRAP(MatrixCell H) {
        return minimalRepFromMRAP(H, "obscont", 1e-12);
    }

    public static MatrixCell minimalRepFromMRAP(MatrixCell H, String how) {
        return minimalRepFromMRAP(H, how, 1e-12);
    }

    /**
     * Returns the minimal representation of a marked rational arrival process.
     */
    public static MatrixCell minimalRepFromMRAP(MatrixCell H, String how, double precision) {
        int N = H.get(0).getNumRows();
        int M = H.size();

        if ("cont".equals(how)) {
            // Controllability reduction
            List<Matrix> matrices = new ArrayList<Matrix>(M);
            for (int i = 0; i < M; i++) {
                matrices.add(H.get(i));
            }
            Matrix ones = Matrix.ones(N, 1);
            Pair<Matrix, Integer> result = MStaircase.mStaircase(matrices, ones, precision);
            Matrix B = result.getLeft();
            int n = result.getRight();
            if (n < N) {
                Matrix Bi = B.inv();
                MatrixCell ret = new MatrixCell(M);
                for (int i = 0; i < M; i++) {
                    Matrix transformed = Bi.mult(H.get(i)).mult(B);
                    ret.set(i, Matrix.getSubMatrix(transformed, 0, n, 0, n));
                }
                return ret;
            } else {
                MatrixCell ret = new MatrixCell(M);
                for (int i = 0; i < M; i++) {
                    ret.set(i, H.get(i).copy());
                }
                return ret;
            }
        } else if ("obs".equals(how)) {
            // Observability reduction
            // Get marginal distribution alpha from MRAP
            PHRepresentation alpha = MarginalDistributionFromMRAP.marginalDistributionFromMRAP(H, precision);

            List<Matrix> matrices = new ArrayList<Matrix>(M);
            for (int i = 0; i < M; i++) {
                matrices.add(H.get(i).transpose());
            }
            Pair<Matrix, Integer> result = MStaircase.mStaircase(matrices, alpha.alpha.transpose(), precision);
            Matrix B = result.getLeft();
            int n = result.getRight();
            if (n < N) {
                Matrix Bi = B.inv();
                MatrixCell ret = new MatrixCell(M);
                for (int i = 0; i < M; i++) {
                    Matrix transformed = Bi.mult(H.get(i)).mult(B);
                    ret.set(i, Matrix.getSubMatrix(transformed, 0, n, 0, n));
                }
                return ret;
            } else {
                MatrixCell ret = new MatrixCell(M);
                for (int i = 0; i < M; i++) {
                    ret.set(i, H.get(i).copy());
                }
                return ret;
            }
        } else if ("obscont".equals(how)) {
            // First controllability, then observability
            MatrixCell D = minimalRepFromMRAP(H, "cont", precision);
            return minimalRepFromMRAP(D, "obs", precision);
        } else {
            throw new IllegalArgumentException("MinimalRepFromMRAP: Unknown method '" + how + "'");
        }
    }

    /**
     * Overload for Matrix[].
     */
    public static MatrixCell minimalRepFromMRAP(Matrix[] H, String how, double precision) {
        MatrixCell cell = new MatrixCell(H.length);
        for (int i = 0; i < H.length; i++) {
            cell.set(i, H[i]);
        }
        return minimalRepFromMRAP(cell, how, precision);
    }
}
