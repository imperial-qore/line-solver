/**
 * @file Markovian Arrival Process joint PDF derivative analysis
 *
 * Computes partial derivatives of MAP joint probability density functions at origin.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_jointpdf_derivative {
    private Map_jointpdf_derivative() {}

    /**
     * Compute partial derivative at 0 of a MAP's joint PDF.
     *
     * @param MAP MatrixCell containing D0 and D1 matrices of the MAP
     * @param iset Index set for the derivative computation
     * @return Scalar value representing the derivative
     */
    public static double map_jointpdf_derivative(MatrixCell MAP, int[] iset) {
        Matrix D0 = MAP.get(0);
        Matrix D1 = MAP.get(1);
        int n = D0.getNumRows();

        Matrix gamma = Map_pie.map_pie(MAP);

        for (int j : iset) {
            Matrix D0PowerJ = Matrix.pow(D0, j);
            gamma = gamma.mult(D0PowerJ).mult(D1);
        }

        Matrix ones = Matrix.ones(n, 1);
        return gamma.mult(ones).get(0, 0);
    }

    /**
     * Compute partial derivative at 0 of a MAP's joint PDF.
     *
     * @param D0 Hidden transition matrix of the MAP
     * @param D1 Visible transition matrix of the MAP
     * @param iset Index set for the derivative computation
     * @return Scalar value representing the derivative
     */
    public static double map_jointpdf_derivative(Matrix D0, Matrix D1, int[] iset) {
        MatrixCell MAP = new MatrixCell(new Matrix[]{D0, D1});
        return map_jointpdf_derivative(MAP, iset);
    }
}
