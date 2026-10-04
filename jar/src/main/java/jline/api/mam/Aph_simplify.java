/**
 * @file Absorbing Phase-type distribution simplification and combination
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Aph_simplify {
    private Aph_simplify() {}

    /**
     * Simplifies and combines two APH distributions using different structural patterns.
     */
    public static Pair<Matrix, Matrix> aph_simplify(Matrix a1, Matrix T1, Matrix a2, Matrix T2,
                                                     double p1, double p2, int pattern) {
        if (pattern == 1) {
            return sequencePattern(a1, T1, a2, T2);
        } else if (pattern == 2) {
            return parallelPattern(a1, T1, a2, T2);
        } else if (pattern == 3) {
            return branchPattern(a1, T1, a2, T2, p1, p2);
        } else {
            throw new IllegalArgumentException("Pattern must be 1 (sequence), 2 (parallel), or 3 (branch)");
        }
    }

    /**
     * Sequence structure: first APH followed by second APH.
     */
    private static Pair<Matrix, Matrix> sequencePattern(Matrix a1, Matrix T1, Matrix a2, Matrix T2) {
        int order1 = a1.getNumCols();
        int order2 = a2.getNumCols();
        Matrix e1 = Matrix.ones(order1, 1);

        Matrix a1_e1 = a1.mult(e1);
        Matrix one_minus_a1_e1 = Matrix.ones(1, 1).sub(1.0, a1_e1);
        Matrix alpha_right = one_minus_a1_e1.mult(a2);
        Matrix alpha = Matrix.concatColumns(a1, alpha_right, null);

        Matrix minusT1_e1 = T1.scale(-1.0).mult(e1);
        Matrix T_top_right = minusT1_e1.mult(a2);
        Matrix T_bottom_left = Matrix.zeros(order2, order1);

        Matrix T_top = Matrix.concatColumns(T1, T_top_right, null);
        Matrix T_bottom = Matrix.concatColumns(T_bottom_left, T2, null);
        Matrix T = Matrix.concatRows(T_top, T_bottom, null);

        return new Pair<Matrix, Matrix>(alpha, T);
    }

    /**
     * Parallel structure: both APHs run in parallel.
     */
    private static Pair<Matrix, Matrix> parallelPattern(Matrix a1, Matrix T1, Matrix a2, Matrix T2) {
        int order1 = a1.getNumCols();
        int order2 = a2.getNumCols();
        Matrix e1 = Matrix.ones(order1, 1);
        Matrix e2 = Matrix.ones(order2, 1);

        Matrix kron_a1_a2 = kroneckerProduct(a1, a2);
        Matrix one_minus_a2_e2 = Matrix.ones(1, 1).sub(1.0, a2.mult(e2));
        Matrix alpha_middle = one_minus_a2_e2.mult(a1);
        Matrix one_minus_a1_e1 = Matrix.ones(1, 1).sub(1.0, a1.mult(e1));
        Matrix alpha_right = one_minus_a1_e1.mult(a2);

        Matrix alpha_temp = Matrix.concatColumns(kron_a1_a2, alpha_middle, null);
        Matrix alpha = Matrix.concatColumns(alpha_temp, alpha_right, null);

        Matrix eye_order1 = Matrix.eye(order1);
        Matrix eye_order2 = Matrix.eye(order2);
        Matrix kron_T1_eye2 = kroneckerProduct(T1, eye_order2);
        Matrix kron_eye1_T2 = kroneckerProduct(eye_order1, T2);
        Matrix kron_eye1_minusT2_e2 = kroneckerProduct(eye_order1, T2.scale(-1.0).mult(e2));
        Matrix kron_minusT1_e1_eye2 = kroneckerProduct(T1.scale(-1.0).mult(e1), eye_order2);

        Matrix Tr1_left = kron_T1_eye2.add(1.0, kron_eye1_T2);
        Matrix Tr1_temp = Matrix.concatColumns(Tr1_left, kron_eye1_minusT2_e2, null);
        Matrix Tr1 = Matrix.concatColumns(Tr1_temp, kron_minusT1_e1_eye2, null);

        Matrix Tr2_temp = Matrix.concatColumns(Matrix.zeros(order1, order1 * order2), T1, null);
        Matrix Tr2 = Matrix.concatColumns(Tr2_temp, Matrix.zeros(order1, order2), null);

        Matrix Tr3_temp = Matrix.concatColumns(Matrix.zeros(order2, order1 * order2), Matrix.zeros(order2, order1), null);
        Matrix Tr3 = Matrix.concatColumns(Tr3_temp, T2, null);

        Matrix T_temp = Matrix.concatRows(Tr1, Tr2, null);
        Matrix T = Matrix.concatRows(T_temp, Tr3, null);

        return new Pair<Matrix, Matrix>(alpha, T);
    }

    /**
     * Branch structure: select between two APHs with given probabilities.
     */
    private static Pair<Matrix, Matrix> branchPattern(Matrix a1, Matrix T1, Matrix a2, Matrix T2, double p1, double p2) {
        int order1 = a1.getNumCols();
        int order2 = a2.getNumCols();

        Matrix alpha_left = a1.scale(p1);
        Matrix alpha_right = a2.scale(p2);
        Matrix alpha = Matrix.concatColumns(alpha_left, alpha_right, null);

        Matrix T_top = Matrix.concatColumns(T1, Matrix.zeros(order1, order2), null);
        Matrix T_bottom = Matrix.concatColumns(Matrix.zeros(order2, order1), T2, null);
        Matrix T = Matrix.concatRows(T_top, T_bottom, null);

        return new Pair<Matrix, Matrix>(alpha, T);
    }

    /**
     * Computes the Kronecker product of two matrices.
     */
    private static Matrix kroneckerProduct(Matrix A, Matrix B) {
        Matrix result = new Matrix(A.getNumRows() * B.getNumRows(), A.getNumCols() * B.getNumCols());

        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                double aij = A.get(i, j);
                for (int k = 0; k < B.getNumRows(); k++) {
                    for (int l = 0; l < B.getNumCols(); l++) {
                        result.set(i * B.getNumRows() + k, j * B.getNumCols() + l, aij * B.get(k, l));
                    }
                }
            }
        }

        return result;
    }
}
