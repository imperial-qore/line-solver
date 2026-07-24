/**
 * @file Marked Markovian Arrival Process normalization and sanitization
 *
 * Normalizes MMAP matrices to ensure feasibility and mathematical validity.
 * Essential for maintaining proper stochastic properties in multiclass models.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_normalize {
    private Mmap_normalize() {}

    /**
     * Normalizes a Markovian Arrival Process with marked arrivals (MMAP) to ensure feasibility.
     *
     * This method adjusts the MMAP by setting negative off-diagonal values in the D0 matrix to zero and ensuring that
     * all elements in the marking matrices are non-negative. It also recalculates the D1 matrix to reflect the sum of
     * the marking matrices. The diagonal elements of the D0 matrix are adjusted to ensure that each row sums to zero,
     * maintaining the properties of a valid generator matrix.
     *
     * @param MMAP the MatrixCell representing the MMAP to be normalized
     * @return the normalized MatrixCell, or null if the input MMAP is empty
     */
    public static MatrixCell mmap_normalize(MatrixCell MMAP) {
        if (MMAP.isEmpty()) {
            return null;
        }
        int K = MMAP.get(0).getNumRows();
        int C = MMAP.size() - 2;

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                if (i != j) {
                    MMAP.get(0).set(i, j, FastMath.max(MMAP.get(0).get(i, j), 0.0));
                }
            }
        }

        MMAP.set(1, new Matrix(MMAP.get(0).getNumRows(), MMAP.get(0).getNumCols(),
                MMAP.get(0).getNumRows() * MMAP.get(0).getNumCols()));

        for (int c = 0; c < C; c++) {
            MMAP.get(2 + c).removeNegative();
            // MATLAB: "if isnan(MMAP{2+c})" holds only when EVERY element is NaN
            // (an if-condition on a matrix is true iff all entries are true), and is
            // false for an empty matrix. Testing only element (0,0) diverges in both
            // directions: it zeroes a partially-NaN mark that MATLAB keeps, and keeps
            // a mark whose first entry is finite but whose others are NaN, letting the
            // NaN reach D1 and the D0 diagonal.
            if (isAllNaN(MMAP.get(2 + c))) {
                MMAP.set(2 + c, new Matrix(MMAP.get(2 + c).getNumRows(), MMAP.get(2 + c).getNumCols(),
                        MMAP.get(2 + c).getNumRows() * MMAP.get(2 + c).getNumCols()));
            }
            MMAP.set(1, MMAP.get(1).add(1.0, MMAP.get(2 + c)));
        }

        for (int k = 0; k < K; k++) {
            MMAP.get(0).set(k, k, 0);
            MMAP.get(0).set(k, k, -MMAP.get(0).sumRows(k) - MMAP.get(1).sumRows(k));
        }

        return MMAP;
    }

    /**
     * Mirrors the MATLAB "if isnan(X)" test on a matrix: true iff X is non-empty and
     * every entry is NaN.
     *
     * @param X the matrix to test
     * @return true if X is non-empty and all of its entries are NaN
     */
    private static boolean isAllNaN(Matrix X) {
        int rows = X.getNumRows();
        int cols = X.getNumCols();
        if (rows == 0 || cols == 0) {
            return false;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (!Double.isNaN(X.get(i, j))) {
                    return false;
                }
            }
        }
        return true;
    }
}
