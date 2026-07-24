/**
 * @file Squared L2 distance between lag-L joint densities of two MAPs
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;

public final class Map_dist {
    private Map_dist() {}

    /**
     * Computes the squared L2 distance between lag-L joint densities of two MAPs.
     * The stationary vectors at arrivals are computed internally.
     *
     * @param A0 hidden transition matrix of the first MAP
     * @param A1 visible transition matrix of the first MAP
     * @param B0 hidden transition matrix of the second MAP
     * @param B1 visible transition matrix of the second MAP
     * @param L  lag parameter
     * @return the squared L2 distance
     */
    public static double map_dist(Matrix A0, Matrix A1, Matrix B0, Matrix B1, int L) {
        Matrix alA = Map_pie.map_pie(A0, A1);
        Matrix alB = Map_pie.map_pie(B0, B1);
        return map_dist(A0, A1, B0, B1, L, alA, alB);
    }

    /**
     * Computes the squared L2 distance between lag-L joint densities of two MAPs.
     *
     * @param A0  hidden transition matrix of the first MAP
     * @param A1  visible transition matrix of the first MAP
     * @param B0  hidden transition matrix of the second MAP
     * @param B1  visible transition matrix of the second MAP
     * @param L   lag parameter
     * @param alA stationary vector at arrivals of the first MAP
     * @param alB stationary vector at arrivals of the second MAP
     * @return the squared L2 distance
     */
    public static double map_dist(Matrix A0, Matrix A1, Matrix B0, Matrix B1, int L, Matrix alA, Matrix alB) {
        return Map_exp_mul_int.map_exp_mul_int(A0, A1, A0, A1, L + 1, alA, alA)
                - 2 * Map_exp_mul_int.map_exp_mul_int(A0, A1, B0, B1, L + 1, alA, alB)
                + Map_exp_mul_int.map_exp_mul_int(B0, B1, B0, B1, L + 1, alB, alB);
    }
}
