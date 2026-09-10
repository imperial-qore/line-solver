/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * Andras Horvath, Gabor Horvath, Miklos Telek, "A traffic based decomposition
 * of two-class queueing networks with priority service," Computer Networks
 * 53:(8) pp. 1235-1248. (2009)
 */
package jline.lib.butools.map;

import jline.lib.butools.ph.MEFromMoments;
import jline.lib.butools.ph.MERepresentation;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MRAPFromMoments {
    private MRAPFromMoments() {}

    /**
     * Creates a continuous marked rational arrival process that
     * has the same marginal and lag-1 joint moments as given.
     */
    public static MatrixCell mrapFromMoments(double[] moms, MatrixCell Nm) {
        MERepresentation meResult = MEFromMoments.meFromMoments(moms);
        Matrix v = meResult.alpha;
        Matrix H0 = meResult.A;

        int N = H0.getNumRows();

        // KEY CHANGE from DMRAP: H0i = (-H0)^{-1} instead of (I-H0)^{-1}
        Matrix H0i = H0.neg().inv();
        Matrix Ge = new Matrix(N, N);
        Matrix G1 = new Matrix(N, N);

        Matrix H0ip = Matrix.eye(N);
        for (int i = 0; i < N; i++) {
            // Ge(i,:) = v * H0ip
            Matrix row = v.mult(H0ip);
            for (int j = 0; j < N; j++) {
                Ge.set(i, j, row.get(0, j));
            }

            // G1(:,i) = sum(H0ip, 2) -- row sums of H0ip
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int k = 0; k < N; k++) {
                    sum += H0ip.get(j, k);
                }
                G1.set(j, i, sum);
            }

            // H0ip = H0ip * (i+1) * H0i
            H0ip = H0ip.scale((double) (i + 1)).mult(H0i);
        }

        Matrix Gei = Ge.inv();
        Matrix G1i = G1.inv();

        int numTypes = Nm.size();
        MatrixCell H = new MatrixCell(numTypes + 1);
        H.set(0, H0);

        for (int i = 0; i < numTypes; i++) {
            // H{i+1} = (-H0) * Gei * Nm{i} * G1i
            // Continuous case: use Nm directly (no factorial moment conversion)
            H.set(i + 1, H0.neg().mult(Gei).mult(Nm.get(i)).mult(G1i));
        }

        return H;
    }

    /**
     * Overload for Matrix[].
     */
    public static MatrixCell mrapFromMoments(double[] moms, Matrix[] Nm) {
        MatrixCell cell = new MatrixCell(Nm.length);
        for (int i = 0; i < Nm.length; i++) {
            cell.set(i, Nm[i]);
        }
        return mrapFromMoments(moms, cell);
    }

    /**
     * Creates a rational arrival process that has the same
     * marginal and lag-1 joint moments as given.
     *
     * Convenience wrapper for single arrival type (K=1).
     */
    public static Pair<Matrix, Matrix> rapFromMoments(double[] moms, Matrix Nm) {
        MatrixCell result = mrapFromMoments(moms, new Matrix[]{Nm});
        return new Pair<Matrix, Matrix>(result.get(0), result.get(1));
    }
}
