/**
 * @file Hessian matrix computation for Logistic expansion method
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.util.matrix.Matrix;

public final class Pfqn_le_hessian {
    private Pfqn_le_hessian() {}

    /**
     * Auxiliary function to compute the Hessian used in the logistic expansion method.
     *
     * @param L  - demands at all stations
     * @param N  - number of jobs for each class
     * @param u0 - term appearing in the integrand
     * @return the Hessian matrix
     */
    public static Matrix pfqn_le_hessian(Matrix L, Matrix N, Matrix u0) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double Ntot = N.elementSum();
        Matrix hu = new Matrix(M - 1, M - 1);
        hu.fill(0.0);

        for (int i = 0; i < M - 1; i++) {
            for (int j = 0; j < M - 1; j++) {
                if (i != j) {
                    hu.set(i, j, -(Ntot + M) * u0.get(i) * u0.get(j));
                    for (int r = 0; r < R; r++) {
                        Matrix L_col_r = new Matrix(L.getNumRows(), 1);
                        Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
                        double denom = u0.mult(L_col_r).get(0) * u0.mult(L_col_r).get(0);
                        hu.set(i, j, hu.get(i, j) + N.get(r) * L.get(i, r) * L.get(j, r) * u0.get(i) * u0.get(j) / denom);
                    }
                } else {
                    hu.set(i, j, (Ntot + M) * u0.get(i) * Matrix.allbut(u0, i).elementSum());
                    for (int r = 0; r < R; r++) {
                        Matrix L_col_r = new Matrix(L.getNumRows(), 1);
                        Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
                        Matrix tmp_L = Matrix.allbut(L_col_r, i).transpose();
                        double denom = u0.mult(L_col_r).get(0) * u0.mult(L_col_r).get(0);
                        hu.set(i, j, hu.get(i, j) - (N.get(r) * L.get(i, r) * u0.get(i) * Matrix.allbut(u0, i).mult(tmp_L).get(0)) / denom);
                    }
                }
            }
        }
        return hu;
    }
}
