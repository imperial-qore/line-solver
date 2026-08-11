/**
 * @file Markovian Arrival Process complementary CDF derivative analysis
 *
 * Computes derivatives of MAP complementary cumulative distribution functions at zero.
 * Used for advanced moment analysis and joint queue analysis in MAP/MAP/1 queueing systems.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_ccdf_derivative {
    private Map_ccdf_derivative() {}

    /**
     * Compute derivative at 0 of a MAP complementary CDF.
     *
     * <p>Based on: A. Horvath et al. A Joint Moments Based Analysis of Networks of MAP/MAP/1 Queues
     *
     * @param MAP MatrixCell containing D0 and D1 matrices of the MAP
     * @param i Derivative order
     * @return Scalar value representing the derivative
     */
    public static double map_ccdf_derivative(MatrixCell MAP, int i) {
        Matrix D0 = MAP.get(0);
        int n = D0.getNumRows();

        Matrix pie = Map_pie.map_pie(MAP);
        Matrix D0PowerI = Matrix.pow(D0, i);
        Matrix ones = Matrix.ones(n, 1);

        return pie.mult(D0PowerI).mult(ones).get(0, 0);
    }

    /**
     * Compute derivative at 0 of a MAP complementary CDF.
     *
     * @param D0 Hidden transition matrix of the MAP
     * @param D1 Visible transition matrix of the MAP
     * @param i Derivative order
     * @return Scalar value representing the derivative
     */
    public static double map_ccdf_derivative(Matrix D0, Matrix D1, int i) {
        MatrixCell MAP = new MatrixCell(new Matrix[]{D0, D1});
        return map_ccdf_derivative(MAP, i);
    }
}
