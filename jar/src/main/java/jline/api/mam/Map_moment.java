/**
 * @file Markovian Arrival Process raw moment computation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_moment {
    private Map_moment() {}

    /**
     * Computes the raw moments of the inter-arrival times of a MAP.
     *
     * @param D0    the hidden transition matrix of the MAP
     * @param D1    the visible transition matrix of the MAP
     * @param order the moment order
     * @return the raw moment of the inter-arrival times of the specified order
     */
    public static double map_moment(Matrix D0, Matrix D1, int order) {
        boolean isZeroMatrix = true;
        for (int i = 0; i < D0.getNumRows(); i++) {
            for (int j = 0; j < D0.getNumCols(); j++) {
                if (Math.abs(D0.get(i, j)) > 1e-14) {
                    isZeroMatrix = false;
                    break;
                }
            }
            if (!isZeroMatrix) break;
        }

        if (isZeroMatrix || Math.abs(D0.det()) < 1e-12) {
            return 0.0;
        }

        Matrix pie = Map_pie.map_pie(D0, D1);
        Matrix iD0 = D0.copy();
        iD0.scaleEq(-1.0);
        iD0 = iD0.inv();
        Matrix iD0k = new Matrix(iD0);
        for (int i = 2; i <= order; i++) {
            iD0k = iD0k.mult(iD0);
            iD0k.scaleEq((double) i);
        }
        Matrix e = Matrix.ones(D0.getNumRows(), 1);
        pie = pie.mult(iD0k);
        pie = pie.mult(e);
        return pie.toDouble();
    }

    /**
     * Computes the raw moments using a MatrixCell.
     */
    public static double map_moment(MatrixCell MAP, int order) {
        return map_moment(MAP.get(0), MAP.get(1), order);
    }
}
