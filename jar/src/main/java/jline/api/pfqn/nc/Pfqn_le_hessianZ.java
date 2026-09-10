/**
 * @file Hessian matrix computation for Logistic expansion method with think times
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.util.matrix.Matrix;

final class Pfqn_le_hessianZ {
    private Pfqn_le_hessianZ() {}

    /**
     * Auxiliary function to compute the Hessian used in the logistic expansion method in models with delays.
     */
    static Matrix pfqn_le_hessianZ(Matrix L, Matrix N, Matrix Z, Matrix u, double v) {
        int K = L.getNumRows();
        int R = L.getNumCols();
        double Ntot = N.elementSum();
        Matrix A = new Matrix(K, K);
        A.fill(0.0);
        Matrix csi = new Matrix(1, R);
        // csi2N is csi(r)^2/N(r) rewritten as N(r)/c(r)^2. Identical where both are
        // defined, but 0 rather than 0/0 for an empty class, which oner() makes
        // routine in the mean-value pipeline of Pfqn_nc.
        Matrix csi2N = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double c = Z.get(r) + v * u.mult(Matrix.extractColumn(L, r, null)).get(0);
            csi.set(r, N.get(r) / c);
            csi2N.set(r, N.get(r) / (c * c));
        }
        Matrix Lhat = new Matrix(K, R);
        Lhat.fill(0.0);
        for (int k = 0; k < K; k++) {
            for (int r = 0; r < R; r++) {
                Lhat.set(k, r, Z.get(r) + v * L.get(k, r));
            }
        }
        double eta = Ntot + K;
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                if (i != j) {
                    A.set(i, j, -eta * u.get(i) * u.get(j));
                    for (int r = 0; r < R; r++) {
                        A.set(i, j, A.get(i, j) + csi2N.get(r) * Lhat.get(i, r) * Lhat.get(j, r) * u.get(i) * u.get(j));
                    }
                }
            }
        }
        for (int i = 0; i < K; i++) {
            A.set(i, i, -Matrix.allbut(Matrix.extractRows(A, i, i + 1, null), i).elementSum());
        }
        Matrix tmp_A = new Matrix(K, K);
        tmp_A.fill(0.0);
        Matrix.extract(A, 0, K - 1, 0, K - 1, tmp_A, 0, 0);
        A = tmp_A;
        A.set(K - 1, K - 1, 1.0);

        for (int r = 0; r < R; r++) {
            Matrix L_col_r = new Matrix(L.getNumRows(), 1);
            Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
            A.set(K - 1, K - 1, A.get(K - 1, K - 1)
                    - csi2N.get(r) * Z.get(r) * u.mult(L_col_r).get(0));
        }

        A.set(K - 1, K - 1, v * A.get(K - 1, K - 1));

        for (int i = 0; i < K - 1; i++) {
            for (int r = 0; r < R; r++) {
                Matrix L_col_r = new Matrix(L.getNumRows(), 1);
                Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
                A.set(i, K - 1, A.get(i, K - 1)
                        + v * u.get(i) * ((csi2N.get(r) * Lhat.get(i, r) * (u.mult(L_col_r).get(0))) - csi.get(r) * L.get(i, r)));
            }
            A.set(K - 1, i, A.get(i, K - 1));
        }

        return A;
    }
}
