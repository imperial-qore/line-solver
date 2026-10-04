/**
 * @file Fixed-point iteration for Logistic expansion method
 *
 * Implements the fixed-point iteration algorithm used in the Logistic expansion method
 * for computing normalizing constants. Iteratively solves for the saddle point of the
 * integrand to enable accurate asymptotic approximation.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_le_fpi {
    private Pfqn_le_fpi() {}

    /**
     * Fixed-point iteration used in the logistic expansion method.
     *
     * @param L demands at all stations
     * @param N number of jobs for each class
     * @return fixed point
     */
    public static Ret.pfqnLeFpi pfqn_le_fpi(Matrix L, Matrix N) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix u = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            u.set(i, 1.0 / M);
        }
        Matrix u_1 = new Matrix(M, 1);
        u_1.fill(GlobalConstants.Inf);
        Matrix d = new Matrix(0, M);
        Matrix u_abs_diff = new Matrix(M, 1);

        for (int i = 0; i < M; i++) {
            u_abs_diff.set(i, FastMath.abs(u.get(i) - u_1.get(i)));
        }

        while (u_abs_diff.elementSum() > 1e-10) {
            u_1 = u.copy();
            for (int i = 0; i < M; i++) {
                u.set(i, 1 / (N.elementSum() + M));
                for (int r = 0; r < R; r++) {
                    Matrix L_col_r = new Matrix(L.getNumRows(), 1);
                    Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
                    u.set(i, u.get(i)
                            + N.get(r) / (N.elementSum() + M) * L.get(i, r) * u_1.get(i)
                            / (u_1.transpose().mult(L_col_r).get(0)));
                }
            }

            for (int i = 0; i < M; i++) {
                u_abs_diff.set(i, FastMath.abs(u.get(i) - u_1.get(i)));
            }

            d = Matrix.concatRows(d, u_abs_diff.transpose(), null);
        }

        return new Ret.pfqnLeFpi(u, d);
    }
}
