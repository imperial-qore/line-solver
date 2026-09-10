/**
 * @file MAP to Phase-Type distribution conversion
 *
 * Converts a Markovian Arrival Process (MAP) to a Phase-Type (PH) distribution
 * and its associated PH-renewal process representation.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Triple;

public final class Map2ph {
    private Map2ph() {}

    /**
     * Converts a MAP to a Phase-Type (PH) distribution.
     *
     * Given a MAP {D0, D1}, extracts the PH distribution (alpha, T) where:
     * - T = D0 (the subgenerator)
     * - alpha = map_pie(MAP) (the embedded steady-state probability vector)
     *
     * Also returns the PH-renewal process PHR = {D0, D1'} where
     * D1' = D1 * ones * alpha, which is a MAP whose inter-arrival times
     * are i.i.d. with the extracted PH distribution.
     *
     * @param MAP a MatrixCell containing {D0, D1}
     * @return Triple of (alpha, T, PHR) where alpha is the initial probability vector,
     *         T is the subgenerator matrix, and PHR is the PH-renewal process as a MatrixCell
     */
    public static Triple<Matrix, Matrix, MatrixCell> map2ph(MatrixCell MAP) {
        Matrix T = MAP.get(0);
        Matrix alpha = Map_pie.map_pie(MAP);
        int n = MAP.get(1).getNumRows();
        Matrix PHR_D1 = MAP.get(1).mult(Matrix.ones(n, 1)).mult(alpha);
        MatrixCell PHR = new MatrixCell();
        PHR.set(0, T);
        PHR.set(1, PHR_D1);
        return new Triple<Matrix, Matrix, MatrixCell>(alpha, T, PHR);
    }
}
