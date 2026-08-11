/**
 * @file Squared L2 distance between autocorrelation functions of two MAPs
 *
 * Computes the squared L2 distance between the autocorrelation functions of
 * two Markovian Arrival Processes using geometric sums and discrete Lyapunov equations.
 *
 * Reference:
 *   G. Horvath, "Measuring the distance between MAPs and some
 *   applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
 *   https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;

public final class Map_dist_acf {
    private Map_dist_acf() {}

    /**
     * Computes the geometric sum needed for the autocorrelation distance.
     */
    public static double map_geo_mul_sum(Matrix A0, Matrix A1, Matrix B0, Matrix B1, Matrix alA, Matrix alB) {
        Matrix negA0 = A0.copy(); negA0.scaleEq(-1.0);
        Matrix negB0 = B0.copy(); negB0.scaleEq(-1.0);
        Matrix A0i = negA0.inv();
        Matrix B0i = negB0.inv();
        int NA = A0.getNumRows();
        int NB = B0.getNumRows();

        Matrix PAh = A0i.mult(A1).add(-1.0, Matrix.ones(NA, 1).mult(alA));
        Matrix PBh = B0i.mult(B1).add(-1.0, Matrix.ones(NB, 1).mult(alB));

        Matrix M = Matrix.eye(NA * NB).add(-1.0, PBh.transpose().kron(PAh));
        if (Math.abs(M.det()) < 1e-10) {
            return Double.MAX_VALUE;
        } else {
            Matrix X = Matrix.dlyap(PAh, PBh, A0i.sumRows().mult(alB.mult(B0i)));
            // The result is the total inner product = sum of every entry of the
            // 1xNB row vector alA*A0i*X*B0i. sumCols().toDouble() returned only
            // the first column (dropping the rest); use elementSum instead.
            return alA.mult(A0i).mult(X).mult(B0i).elementSum();
        }
    }

    /**
     * Computes the squared L2 distance between autocorrelation functions of two MAPs.
     *
     * Stationary vectors at arrivals are computed internally.
     */
    public static double map_dist_acf(Matrix A0, Matrix A1, Matrix B0, Matrix B1) {
        Matrix alA = Map_pie.map_pie(A0, A1);
        Matrix alB = Map_pie.map_pie(B0, B1);
        return map_dist_acf(A0, A1, B0, B1, alA, alB);
    }

    /**
     * Computes the squared L2 distance between autocorrelation functions of two MAPs.
     */
    public static double map_dist_acf(Matrix A0, Matrix A1, Matrix B0, Matrix B1, Matrix alA, Matrix alB) {
        double momA1 = Map_moment.map_moment(A0, A1, 1);
        double momA2 = Map_moment.map_moment(A0, A1, 2);
        double momB1 = Map_moment.map_moment(B0, B1, 1);
        double momB2 = Map_moment.map_moment(B0, B1, 2);
        double varA = momA2 - momA1 * momA1;
        double varB = momB2 - momB1 * momB1;

        return (map_geo_mul_sum(A0, A1, A0, A1, alA, alA) - momA2 * momA2 / 4) / (varA * varA)
                - 2 * (map_geo_mul_sum(A0, A1, B0, B1, alA, alB) - momA2 * momB2 / 4) / (varA * varB)
                + (map_geo_mul_sum(B0, B1, B0, B1, alB, alB) - momB2 * momB2 / 4) / (varB * varB);
    }
}
