/**
 * @file Markovian Arrival Process feasibility validation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_checkfeasible {
    private Map_checkfeasible() {}

    /**
     * Check the feasibility of a MAP with detailed validation.
     */
    public static boolean map_checkfeasible(MatrixCell MAP, double TOL) {
        int n = MAP.get(0).length();
        Matrix D0 = MAP.get(0);
        Matrix D1 = MAP.get(1);
        if (D0.hasNaN()) {
            return false;
        }
        boolean D0_has_inf = false;
        boolean D1_has_inf = false;

        for (int i = 0; i < D0.getNumRows(); i++) {
            for (int j = 0; j < D0.getNumCols(); j++) {
                if (Utils.isInf(D0.get(i, j))) {
                    D0_has_inf = true;
                    break;
                }
            }
        }

        for (int i = 0; i < D1.getNumRows(); i++) {
            for (int j = 0; j < D1.getNumCols(); j++) {
                if (Utils.isInf(D1.get(i, j))) {
                    D1_has_inf = true;
                    break;
                }
            }
        }

        if (D0.hasNaN() || D1.hasNaN() || D0_has_inf || D1_has_inf) {
            return false;
        }

        Matrix neg_D0 = D0.copy();
        neg_D0.scaleEq(-1.0);

        Matrix P = neg_D0.inv().mult(D1);
        Matrix Q = D0.add(1.0, D1);

        for (int i = 0; i < D0.getNumRows(); i++) {
            for (int j = 0; j < D0.getNumCols(); j++) {
                if (Math.abs(D0.get(i, j)) < TOL) {
                    D0.set(i, j, 0);
                }
            }
        }

        for (int i = 0; i < D1.getNumRows(); i++) {
            for (int j = 0; j < D1.getNumCols(); j++) {
                if (Math.abs(D1.get(i, j)) < TOL) {
                    D1.set(i, j, 0);
                }
            }
        }

        for (int i = 0; i < P.getNumCols(); i++) {
            for (int j = 0; j < P.getNumCols(); j++) {
                if (Math.abs(P.get(i, j)) < TOL) {
                    P.set(i, j, 0);
                }
            }
        }

        for (int i = 0; i < Q.getNumCols(); i++) {
            for (int j = 0; j < Q.getNumCols(); j++) {
                if (Math.abs(Q.get(i, j)) < TOL) {
                    Q.set(i, j, 0);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && D0.get(i, j) < 0) {
                    return false;
                }
                if (i == j && D0.get(i, j) > 0) {
                    return false;
                }
                if (D1.get(i, j) < 0) {
                    return false;
                }
                if (i != j && Q.get(i, j) < 0) {
                    return false;
                }
                if (i == j && Q.get(i, j) > 0) {
                    return false;
                }
                if (P.get(i, j) < 0) {
                    return false;
                }
            }

            if (Math.abs(P.sumRows(i)) < 1 - n * TOL) {
                return false;
            }

            if (Math.abs(P.sumRows(i)) > 1 + n * TOL) {
                return false;
            }
            if (Math.abs(Q.sumRows(i)) < 0 - n * TOL) {
                return false;
            }
            if (Math.abs(Q.sumRows(i)) > 0 + n * TOL) {
                return false;
            }
        }

        if (n < Map_largemap.map_largemap()) {
            double[][] P_data = P.toArray2D();
            RealMatrix P_matrix = MatrixUtils.createRealMatrix(P_data);
            EigenDecomposition P_eigenDecomposition = new EigenDecomposition(P_matrix);
            double[] P_eigenvalues = P_eigenDecomposition.getRealEigenvalues();
            int P_sum = 0;
            for (double e : P_eigenvalues) {
                if (e > 1 - TOL) {
                    P_sum++;
                }
                if (P_sum > 1) {
                    return false;
                }
            }

            double[][] Q_data = Q.toArray2D();
            RealMatrix Q_matrix = MatrixUtils.createRealMatrix(Q_data);
            EigenDecomposition Q_eigenDecomposition = new EigenDecomposition(Q_matrix);
            double[] Q_eigenvalues = Q_eigenDecomposition.getRealEigenvalues();
            int Q_sum = 0;
            for (double e : Q_eigenvalues) {
                if (e > 0 - TOL) {
                    Q_sum++;
                }
                if (Q_sum > 1) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Check the feasibility of a MAP given separate D0 and D1 matrices.
     */
    public static boolean map_checkfeasible(Matrix D0, Matrix D1, double TOL) {
        MatrixCell MAP = new MatrixCell(2);
        MAP.set(0, D0);
        MAP.set(1, D1);
        return map_checkfeasible(MAP, TOL);
    }

    public static boolean map_checkfeasible(Matrix D0, Matrix D1) {
        return map_checkfeasible(D0, D1, 1e-14);
    }
}
