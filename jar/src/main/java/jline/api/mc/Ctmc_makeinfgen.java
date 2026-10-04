/**
 * @file CTMC infinitesimal generator construction and validation
 *
 * Constructs and validates infinitesimal generator matrices for continuous-time
 * Markov chains. Ensures row sums equal zero and non-positive diagonal elements,
 * fundamental requirements for valid CTMC generator matrices.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

public final class Ctmc_makeinfgen {
    private Ctmc_makeinfgen() {}

    /**
     * Converts a matrix into a valid infinitesimal generator for a CTMC.
     * An infinitesimal generator has row sums equal to zero and non-positive diagonal elements.
     *
     * @param Q candidate infinitesimal generator matrix
     * @return valid infinitesimal generator matrix
     */
    public static Matrix ctmc_makeinfgen(Matrix Q) {
        // Extract diagonal elements and create off-diagonal matrix
        double[] diagonalValues = new double[Q.length()];
        for (int i = 0; i < Q.length(); i++) {
            diagonalValues[i] = Q.get(i, i);
        }
        Matrix offDiagonal = Q.sub(Matrix.diagMatrix(null, diagonalValues, 0, diagonalValues.length));

        // Create new diagonal from negative row sums to ensure row sums = 0
        // Use toArray1D() to get ALL row sums (including zeros) to create correct-sized diagonal
        Matrix rowSums = offDiagonal.sumRows();
        double[] rowSumsArray = rowSums.toArray1D();
        Matrix newDiagonal = Matrix.diagMatrix(null, rowSumsArray, 0, rowSumsArray.length);

        // Combine off-diagonal elements with corrected diagonal
        Matrix result = offDiagonal.sub(newDiagonal);
        result.removeZeros(0.0);
        return result;
    }
}
