/**
 * @file Markovian Arrival Process autocorrelation function coefficients for counting processes
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_acfc {
    private Map_acfc() {}

    /**
     * Computes the autocorrelation function coefficients (ACFC) for a MAP counting process.
     *
     * @param D0   the hidden transition matrix of the MAP
     * @param D1   the visible transition matrix of the MAP
     * @param lags an array of integers representing the lags at which to compute the ACFC
     * @param u    the length of the timeslot (timescale)
     * @return an array of doubles containing the ACFC values for the specified lags
     */
    public static double[] map_acfc(Matrix D0, Matrix D1, int[] lags, double u) {
        int n = D0.getNumCols();
        Matrix Q = Map_infgen.map_infgen(D0, D1);
        Matrix I = Matrix.eye(n);
        Matrix piq = Map_piq.map_piq(D0, D1);
        Matrix PRE = piq.mult(D1).mult(I.sub(Maths.matrixExp(Q.scale(u))));
        Matrix e = Matrix.ones(n, 1);
        Matrix inv2_epiqQ = (e.mult(piq).sub(Q)).square().inv();
        Matrix POST = (I.sub(Maths.matrixExp(Q.scale(u)))).mult(inv2_epiqQ).mult(D1).mult(e);
        double vart = Map_varcount.map_varcount(D0, D1, u);
        double[] acfCoeffs = new double[lags.length];
        for (int i = 0; i < lags.length; i++) {
            acfCoeffs[i] = PRE.mult(Maths.matrixExp(Q.scale((lags[i] - 1) * u))).mult(POST).scale(1.0 / vart).value();
        }
        return acfCoeffs;
    }

    /**
     * Computes the autocorrelation function coefficients (ACFC) for a MAP counting process using a MatrixCell.
     */
    public static double[] map_acfc(MatrixCell MAP, int[] lags, double u) {
        return map_acfc(MAP.get(0), MAP.get(1), lags, u);
    }
}
