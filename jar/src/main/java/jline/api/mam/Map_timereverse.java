/**
 * @file Markovian Arrival Process time reversal operations
 *
 * Computes time-reversed MAP by adjusting transition rates based on stationary distributions.
 * Used for analyzing reversibility properties and theoretical MAP characterization.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_timereverse {
    private Map_timereverse() {}

    /**
     * Computes the time-reversed MAP of a given MAP.
     *
     * <p>The time-reversed MAP has transition rates adjusted based on the stationary distribution of the original MAP.
     *
     * @param map the original MAP stored in a MatrixCell containing the D0 and D1 matrices
     * @return a MatrixCell representing the time-reversed MAP, containing the reversed D0 and D1 matrices
     */
    public static MatrixCell map_timereverse(MatrixCell map) {
        Matrix piq = Map_piq.map_piq(map.get(0), map.get(1));
        Matrix D = Matrix.diag(piq.toArray1D());
        Matrix iD = D.inv();
        Matrix MAPr0 = iD.mult(map.get(0).transpose()).mult(D);
        Matrix MAPr1 = iD.mult(map.get(1).transpose()).mult(D);
        return new MatrixCell(MAPr0, MAPr1);
    }
}
