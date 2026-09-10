/**
 * @file Distribution conversion for AoI analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Aoi_dist2ph {
    private Aoi_dist2ph() {}

    /**
     * Convert LINE process representation to PH format for AoI analysis.
     */
    public static Aoi_dist2phResult aoi_dist2ph(MatrixCell proc) {
        if (proc.size() < 2) {
            throw new IllegalArgumentException("proc must be a cell array {D0, D1}");
        }

        Matrix D0 = proc.get(0);
        Matrix D1 = proc.get(1);

        if (D0.hasNaN() || D1.hasNaN()) {
            throw new IllegalArgumentException("Process contains NaN - class may be disabled");
        }

        int n = D0.getNumRows();
        if (n != D0.getNumCols() || n != D1.getNumRows() || n != D1.getNumCols()) {
            throw new IllegalArgumentException("D0 and D1 must be square matrices of the same size");
        }

        Matrix T = D0.copy();

        Matrix Q = D0.add(D1);

        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += Q.get(i, j);
            }
            if (Math.abs(rowSum) > 1e-10) {
                Q.set(i, i, Q.get(i, i) - rowSum);
            }
        }

        Matrix A = new Matrix(n + 1, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A.set(j, i, Q.get(i, j));
            }
        }
        for (int j = 0; j < n; j++) {
            A.set(n, j, 1.0);
        }

        Matrix b = new Matrix(n + 1, 1);
        b.set(n, 0, 1.0);

        Matrix At = A.transpose();
        Matrix AtA = At.mult(A);
        Matrix Atb = At.mult(b);
        Matrix piCol = new Matrix(n, 1);
        Matrix.solveSafe(AtA, Atb, piCol);

        Matrix pi = new Matrix(1, n);
        double piSum = 0.0;
        for (int j = 0; j < n; j++) {
            double v = Math.max(piCol.get(j, 0), 0.0);
            pi.set(0, j, v);
            piSum += v;
        }
        if (piSum > 0) {
            for (int j = 0; j < n; j++) {
                pi.set(0, j, pi.get(0, j) / piSum);
            }
        }

        Matrix ones = Matrix.ones(n, 1);
        Matrix completionRates = D1.mult(ones);

        Matrix alpha = new Matrix(1, n);
        double alphaSum = 0.0;
        for (int j = 0; j < n; j++) {
            double v = pi.get(0, j) * completionRates.get(j, 0);
            alpha.set(0, j, v);
            alphaSum += v;
        }

        if (alphaSum > 0) {
            for (int j = 0; j < n; j++) {
                alpha.set(0, j, alpha.get(0, j) / alphaSum);
            }
        } else {
            for (int j = 0; j < n; j++) {
                alpha.set(0, j, pi.get(0, j));
            }
        }

        double finalSum = 0.0;
        for (int j = 0; j < n; j++) {
            finalSum += alpha.get(0, j);
        }
        if (finalSum > 0 && Math.abs(finalSum - 1.0) > 1e-14) {
            for (int j = 0; j < n; j++) {
                alpha.set(0, j, alpha.get(0, j) / finalSum);
            }
        }

        for (int i = 0; i < n; i++) {
            if (T.get(i, i) > 0) {
                throw new IllegalArgumentException("T(" + i + "," + i + ") = " + T.get(i, i) + " > 0, expected non-positive");
            }
        }

        return new Aoi_dist2phResult(alpha, T);
    }
}
