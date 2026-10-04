/**
 * @file Markovian Arrival Process renewal process construction
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_renewal {
    private Map_renewal() {}

    /**
     * Creates a renewal MAP by removing all correlations from the input MAP.
     *
     * @param MAPIN The input MAP stored in a MatrixCell, containing the (D0, D1) matrices
     * @return A MatrixCell representing the renewal MAP with no correlations
     */
    public static MatrixCell map_renewal(MatrixCell MAPIN) {
        return map_renewal(MAPIN.get(0), MAPIN.get(1));
    }

    /**
     * Creates a renewal MAP by removing all correlations from the input MAP.
     *
     * @param D0 The hidden transition matrix of the input MAP
     * @param D1 The visible transition matrix of the input MAP
     * @return A MatrixCell representing the renewal MAP with no correlations
     */
    public static MatrixCell map_renewal(Matrix D0, Matrix D1) {
        // Compute the pie vector (steady-state probabilities of the embedded DTMC)
        Matrix pie = Map_pie.map_pie(D0, D1);

        // Create the output MAP
        MatrixCell MAPOUT = new MatrixCell();

        // D0 remains the same
        MAPOUT.set(0, D0.copy());

        // D1 = D1 * ones(size(D1,1), 1) * pie
        Matrix ones = Matrix.ones(D1.getNumRows(), 1);
        Matrix D1_ones = D1.mult(ones);
        Matrix D1_new = D1_ones.mult(pie);

        MAPOUT.set(1, D1_new);

        return MAPOUT;
    }
}
