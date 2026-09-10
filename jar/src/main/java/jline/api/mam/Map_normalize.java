/**
 * @file Markovian Arrival Process matrix normalization and sanitization
 *
 * Sanitizes MAP matrices by ensuring non-negativity constraints and proper diagonal adjustments.
 * Essential for maintaining mathematical validity and numerical stability in MAP algorithms.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_normalize {
    private Map_normalize() {}

    /**
     * Sanitizes the (D0, D1) matrices of a Markovian Arrival Process (MAP).
     *
     * @param D0 The hidden transition matrix of the MAP.
     * @param D1 The visible transition matrix of the MAP.
     * @return A MatrixCell containing the sanitized D0 and D1 matrices.
     */
    public static MatrixCell map_normalize(Matrix D0, Matrix D1) {
        MatrixCell D = new MatrixCell();
        D.set(0, D0.copy());
        D.set(1, D1 != null ? D1.copy() : null);
        int nPhases = D0.length();
        for (int i = 0; i < nPhases; i++) {
            for (int j = 0; j < nPhases; j++) {
                if (D.get(0).get(i, j) < 0) {
                    D.get(0).set(i, j, 0.0);
                }
                if (D.get(1).get(i, j) < 0) {
                    D.get(1).set(i, j, 0.0);
                }
            }
        }

        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < nPhases; j++) {
                if (j != i) {
                    rowSum = rowSum + D.get(0).get(i, j);
                }
                rowSum = rowSum + D.get(1).get(i, j);
                D.get(0).set(i, i, -rowSum);
            }
        }
        return D;
    }

    /**
     * Sanitizes the (D0, D1) matrices of a MAP stored in a MatrixCell.
     */
    public static MatrixCell map_normalize(MatrixCell MAP) {
        return map_normalize(MAP.get(0), MAP.get(1));
    }
}
