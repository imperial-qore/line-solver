/**
 * @file Fixed-point iteration for Logistic expansion method with think times
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_le_fpiZ {
    private Pfqn_le_fpiZ() {}

    /**
     * Fixed-point iteration used in the logistic expansion method in models with delays.
     *
     * @param L demands at all stations
     * @param N number of jobs for each class
     * @param Z think time for each class
     * @return fixed point
     */
    public static Ret.pfqnLeFpiZ pfqn_le_fpiZ(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double eta = N.elementSum() + M;
        Matrix u = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            u.set(i, 1.0 / M);
        }
        // Note: eq. (35) in the SIGMETRICS 2017 paper has a spurious +1 in the v
        // equation; the correct stationary point is v = eta - sum_r xi_r*Z_r.
        double v = eta;
        Matrix u_1 = new Matrix(M, 1);
        u_1.fill(GlobalConstants.Inf);
        double v_1 = GlobalConstants.Inf;
        Matrix d = new Matrix(0, M);
        Matrix u_abs_diff = new Matrix(M, 1);

        for (int i = 0; i < M; i++) {
            u_abs_diff.set(i, FastMath.abs(u.get(i) - u_1.get(i)));
        }

        while (u_abs_diff.elementSum() > 1e-10) {
            u_1 = u.copy();
            v_1 = v;
            for (int i = 0; i < M; i++) {
                u.set(i, 1 / eta);
                for (int r = 0; r < R; r++) {
                    Matrix L_col_r = new Matrix(L.getNumRows(), 1);
                    Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
                    u.set(i, u.get(i) + N.get(r) / eta * (Z.get(r) + v * L.get(i, r)) * u_1.get(i)
                            / (Z.get(r) + v * u_1.transpose().mult(L_col_r).get(0)));
                }
            }

            Matrix xi = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                Matrix L_col_r = new Matrix(L.getNumRows(), 1);
                Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
                xi.set(r, N.get(r) / (Z.get(r) + v * u_1.transpose().mult(L_col_r).get(0)));
            }
            v = eta;
            for (int r = 0; r < R; r++) {
                v -= xi.get(r) * Z.get(r);
            }

            for (int i = 0; i < M; i++) {
                u_abs_diff.set(i, FastMath.abs(u.get(i) - u_1.get(i)));
            }

            d = Matrix.concatRows(d, u_abs_diff.transpose().elementIncrease(FastMath.abs(v - v_1)), null);
        }

        return new Ret.pfqnLeFpiZ(u, v, d);
    }
}
