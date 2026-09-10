/**
 * @file Markovian Arrival Process probabilistic mixture models
 *
 * Creates probabilistic mixtures of MAP processes with specified mixture probabilities.
 * Used for modeling heterogeneous arrival patterns and traffic characterization.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_mixture {
    private Map_mixture() {}

    /**
     * Creates a probabilistic mixture of Markovian Arrival Processes (MAPs).
     */
    public static MatrixCell map_mixture(double[] alpha, MatrixCell[] MAPs) {
        if (alpha.length != MAPs.length) {
            throw new IllegalArgumentException("Alpha array and MAPs array must have the same size");
        }

        // Verify that alpha sums to 1.0 (within tolerance)
        double alphaSum = 0.0;
        for (double a : alpha) alphaSum += a;
        if (Math.abs(alphaSum - 1.0) > 1e-10) {
            throw new IllegalArgumentException("Alpha probabilities must sum to 1.0, got " + alphaSum);
        }

        // Pre-compute pie vectors for all MAPs
        Matrix[] pies = new Matrix[MAPs.length];
        for (int i = 0; i < MAPs.length; i++) {
            pies[i] = Map_pie.map_pie(MAPs[i]);
        }

        // Build D0 as block diagonal matrix
        Matrix D0 = null;
        for (int i = 0; i < MAPs.length; i++) {
            if (i == 0) {
                D0 = MAPs[i].get(0).copy();
            } else {
                D0 = blockDiag(D0, MAPs[i].get(0));
            }
        }

        // Build D1 matrix by vertically concatenating D1i blocks
        Matrix D1 = null;
        for (int i = 0; i < MAPs.length; i++) {
            int mapSize = MAPs[i].get(0).getNumRows();
            Matrix ones = Matrix.ones(mapSize, 1);

            // Start with alpha[0] * MAPs[i]{2} * ones * pie[0]
            Matrix D1i = MAPs[i].get(1).mult(ones.scale(alpha[0])).mult(pies[0]);

            // Horizontally concatenate with alpha[j] * MAPs[i]{2} * ones * pie[j] for j=1..n-1
            for (int j = 1; j < MAPs.length; j++) {
                Matrix D1i_j = MAPs[i].get(1).mult(ones.scale(alpha[j])).mult(pies[j]);
                D1i = Matrix.concatColumns(D1i, D1i_j, null);
            }

            // Vertically concatenate to build D1
            if (i == 0) {
                D1 = D1i;
            } else {
                D1 = Matrix.concatRows(D1, D1i, null);
            }
        }

        return Map_normalize.map_normalize(D0, D1);
    }

    /**
     * Creates a block diagonal matrix from two matrices.
     */
    private static Matrix blockDiag(Matrix A, Matrix B) {
        Matrix result = new Matrix(A.getNumRows() + B.getNumRows(), A.getNumCols() + B.getNumCols());

        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                result.set(i, j, A.get(i, j));
            }
        }

        for (int i = 0; i < B.getNumRows(); i++) {
            for (int j = 0; j < B.getNumCols(); j++) {
                result.set(A.getNumRows() + i, A.getNumCols() + j, B.get(i, j));
            }
        }

        return result;
    }
}
