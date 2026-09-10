/**
 * Stationary vector at arrival epochs for discrete-time MAPs.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Dmap_pie {
    private Dmap_pie() {}

    /**
     * Computes the stationary vector at arrival epochs of a discrete-time MAP.
     *
     * The embedded chain is P = (I - D0)^{-1} * D1, and the stationary vector
     * pi satisfies pi * P = pi, pi * e = 1.
     */
    public static Matrix dmap_pie(Matrix D0, Matrix D1) {
        int N = D0.getNumRows();
        Matrix I = Matrix.eye(N);
        Matrix P = I.add(-1.0, D0).inv().mult(D1);
        // Solve (P^T - I) pi^T = 0 with sum(pi) = 1 by replacing the FIRST
        // ROW (equation 0) with the normalization sum(pi)=1. The row index is
        // the equation, so this must set row 0 to all ones, not column 0.
        Matrix M = P.transpose().add(-1.0, I);
        for (int j = 0; j < N; j++) {
            M.set(0, j, 1.0);
        }
        Matrix rhs = new Matrix(N, 1, N);
        rhs.set(0, 0, 1.0);
        return Matrix.robustLeftDivide(M, rhs).transpose();
    }

    public static Matrix dmap_pie(MatrixCell DMAP) {
        return dmap_pie(DMAP.get(0), DMAP.get(1));
    }
}
